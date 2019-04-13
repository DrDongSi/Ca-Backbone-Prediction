"""REST API allowing providing the CA-backbone prediction service

The REST API allows the user to invoke the backbone prediction remotely through
HTTP requests. In particular, this is used by the Chimera client plug-in.

In order to run the prediction pipeline the user has to first create a new job,
create a new protein for each protein which should be predicted, upload the
density map and ground truth for each protein, optionally upload the thresholds,
and then start the prediction. Once the prediction is finished an email with a
link to download the results will be sent to the email specified for the job.
"""

from flask import Flask, request, abort, jsonify, send_file
from api.db import get_db, JobNotFoundError, ProteinNotFoundError, JobStatus
from prediction.prediction import run_predictions
from api.email_service import send_email
from threading import Thread

__author__ = 'Jonas Pfab'


app = Flask(__name__)
db = get_db()


@app.route('/jobs', methods=['POST'])
def create_job():
    """Creates empty job

    In order to run predictions we first have to create a job. This method will
    create the folder which will contain all data related to the job. Data about
    all proteins that should be predicted within the job are put within this
    folder.

    The request body must include a JSON object containing the email to which
    the link with the results are sent to once the prediction is finished.
    """
    data = request.get_json()
    if 'email' not in data:
        abort(400, 'Email is missing')

    email = data['email']
    job_id = db.create_job(email)

    return jsonify({'job_id': job_id})


@app.route('/jobs/<job_id>/proteins', methods=['POST'])
def create_protein(job_id):
    """Creates empty protein

    Creates an empty folder in the input folder of the specific job in which the
    density map and ground truth of the protein will be stored.

    The request body must include a JSON object containing the protein ID of the
    protein created.

    Parameters
    ----------
    job_id: str
        ID of job to which protein belongs to
    """
    data = request.get_json()
    if 'protein_id' not in data:
        abort(400, 'Protein ID is missing')

    try:
        protein_id = data['protein_id']
        db.create_protein(job_id, protein_id)
    except JobNotFoundError as e:
        abort(404, str(e))

    return 'Protein created'


@app.route('/jobs/<job_id>/proteins/<protein_id>/density_map', methods=['POST'])
def upload_density_map(job_id, protein_id):
    """Uploads density map for protein

    Density map must be included as a file in the request body. The file is
    then stored in the input folder of the corresponding protein and job.

    Parameters
    ----------
    job_id: str
        ID of the corresponding job

    protein_id: str
        ID of the corresponding protein
    """
    if 'density_map' not in request.files:
        abort(400, 'Density map is missing')

    try:
        density_map = request.files['density_map']
        print(type(density_map))
        db.store_density_map(job_id, protein_id, density_map)
    except (JobNotFoundError, ProteinNotFoundError) as e:
        abort(404, str(e))

    return 'Density map uploaded'


@app.route('/jobs/<job_id>/proteins/<protein_id>/ground_truth', methods=['POST'])
def upload_ground_truth(job_id, protein_id):
    """Uploads ground truth for protein

    The request body must contain the ground truth file. It is then stored in
    the input folder of the corresponding protein and job

    Parameters
    ----------
    job_id: str
        ID of the corresponding job

    protein_id: str
        ID of the corresponding protein
    """
    if 'ground_truth' not in request.files:
        abort(400, 'Ground truth is missing')

    try:
        ground_truth = request.files['ground_truth']
        db.store_ground_truth(job_id, protein_id, ground_truth)
    except (JobNotFoundError, ProteinNotFoundError) as e:
        abort(404, str(e))

    return 'Ground truth uploaded'


@app.route('/jobs/<job_id>/thresholds', methods=['POST'])
def upload_thresholds(job_id):
    """Uploads thresholds file for job

    Uploads thresholds file which contains the thresholds for all proteins
    (The file is not protein specific, but can contain the thresholds for
    multiple proteins).

    The request body must contain a JSON object describing the thresholds for
    each protein.

    Parameters
    ----------
    job_id: str
        ID of the corresponding job

    Note
    ----------
    This step is optional. If no thresholds are given they will be automatically
    determined during the prediction.
    """
    data = request.get_json()
    if 'thresholds' not in data:
        abort(400, 'Thresholds are missing')

    try:
        thresholds = data['thresholds']
        db.store_thresholds(job_id, thresholds)
    except JobNotFoundError as e:
        abort(404, str(e))

    return 'Thresholds uploaded'


@app.route('/jobs/<job_id>/run', methods=['POST'])
def run_job(job_id):
    """Starts prediction in background thread

    This method will start the prediction of all proteins contained in the job
    in a background thread. Once the prediction is finished it sends an email
    to the address linked to the job with a link from which the results can be
    downloaded

    Parameters
    ----------
    job_id: str
        ID of the job that should be ran
    """
    results_url = request.host_url + 'jobs/' + job_id + '/results'

    if db.get_job_status(job_id) != JobStatus.NOT_STARTED:
        abort(400, 'Job has already been started')

    try:
        input_path = db.get_job_input(job_id) + '/'
        output_path = db.get_job_output(job_id) + '/'
        thresholds_file = db.get_thresholds(job_id)
    except JobNotFoundError as e:
        abort(404, str(e))
    except FileNotFoundError:
        thresholds_file = None

    def background_job():
        db.set_job_status(job_id, JobStatus.STARTED)
        run_predictions(input_path, output_path, thresholds_file, 0, False)
        db.set_job_status(job_id, JobStatus.FINISHED)

        db.create_results(job_id)
        send_email(db.get_job_email(job_id), results_url)

    thread = Thread(target=background_job)
    thread.start()

    return 'Job started'


@app.route('/jobs/<job_id>/results', methods=['GET'])
def get_results(job_id):
    """Returns results.zip file for job

    This should only be called once the job has finished.

    Parameters
    ----------
    job_id: str
        ID of job for which results are returned
    """
    try:
        return send_file(db.get_results(job_id))
    except (JobNotFoundError, ProteinNotFoundError, FileNotFoundError) as e:
        abort(404, str(e))
