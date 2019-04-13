"""This module implements the database which maintains all data related to
prediction jobs

For simplicity, the database is implemented using only the File I/O of the
system.
"""

from datetime import datetime
import random
import string
import os
from os.path import join
import json
import zipfile
from enum import Enum

__author__ = 'Jonas Pfab'


def get_db():
    """Creates database instance

    This methods requires the DB_BASE_PATH environment variable to be set
    which defines the base path at which the database stores its data.
    """
    return DB(base_path=os.environ['DB_BASE_PATH'])


class JobNotFoundError(Exception):
    """Job with given job ID does not exist"""
    def __init__(self, job_id):
        super(Exception, self).__init__('Job ' + job_id + ' does not exist')


class ProteinNotFoundError(Exception):
    """Protein with given protein ID does not exist"""
    def __init__(self, job_id, protein_id):
        super(Exception, self).__init__('Protein ' + protein_id + ' does not exist in job ' + job_id)


class JobStatus(Enum):
    """Different statuses of a job"""
    NOT_STARTED = 0  # Job has not yet started
    STARTED = 1  # Job has started but not finished yet
    FINISHED = 2  # Job has finished


class DB:
    """Encapsulates database operations

    All data elements are currently stored through file I/O. In the future it
    might be better to use an actual database.
    """

    def __init__(self, base_path):
        """Creates DB instance

        Parameters
        ----------
        base_path: str
            Base path at which database maintains its data
        """
        self.base_path = base_path

    def create_job(self, email):
        """Creates a new job

        A new folder is created where all data for this job is maintained. An
        info.json file is added where the email, creation data, and status of
        the job are tracked.

        Parameters
        ----------
        email: str
            Email associated with this job. This will later on be used to send a
            notification when the job has finished.
        """
        job_id = ''.join(random.choices(string.ascii_uppercase + string.digits, k=8))
        os.makedirs(join(self.base_path, job_id, 'output'))
        info = {
            'email': email,
            'creation_date': str(datetime.today()),
            'status': JobStatus.NOT_STARTED.value
        }

        with open(join(self.base_path, job_id, 'info.json'), 'w') as f:
            json.dump(info, f)

        return job_id

    def get_job(self, job_id):
        """Returns base path of job

        Parameters
        ----------
        job_id: str
            ID of job that is retrieved

        Raises
        ----------
        JobNotFoundError
            If job with given ID does not exist
        """
        if not os.path.isdir(join(self.base_path, job_id)):
            raise JobNotFoundError(job_id)

        return join(self.base_path, job_id)

    def get_job_input(self, job_id):
        """Returns path to input folder of job

        Parameters
        ----------
        job_id: str
            ID of job containing input path

        Raises
        ----------
        JobNotFoundError
            If job with given ID does not exist
        """
        return join(self.get_job(job_id), 'input')

    def get_job_output(self, job_id):
        """Returns path to output folder of job

        Parameters
        ----------
        job_id: str
            ID of job containing output path

        Raises
        ----------
        JobNotFoundError
            If job with given ID does not exist
        """
        return join(self.get_job(job_id), 'output')

    def get_job_email(self, job_id):
        """Returns email associated with job

        Parameters
        ----------
        job_id: str
            ID of job for which email is fetched

        Raises
        ----------
        JobNotFoundError
            If job with given ID does not exist
        """
        with open(join(self.get_job(job_id), 'info.json')) as f:
            info = json.load(f)

            return info['email']

    def get_job_status(self, job_id):
        """Returns current status of job

        The current status of the job is maintained in the info.json file.
        Possible values are defined in the JobStatus class.

        Parameters
        ----------
        job_id: str
            ID of job for which status is fetched

        Raises
        ----------
        JobNotFoundError
            If job with given ID does not exist
        """
        with open(join(self.get_job(job_id), 'info.json')) as f:
            info = json.load(f)

            return JobStatus(info['status'])

    def set_job_status(self, job_id, status):
        """Updates current status of job

        The current status of the job is maintained in the info.json file.
        Possible values are defined in the JobStatus class.

        Parameters
        ----------
        job_id: str
            ID of job for which status is updated
        status: JobStatus
            New status of the job

        Raises
        ----------
        JobNotFoundError
            If job with given ID does not exist
        """
        with open(join(self.get_job(job_id), 'info.json')) as f:
            info = json.load(f)
            info['status'] = status.value

        with open(join(self.get_job(job_id), 'info.json'), 'w') as f:
            json.dump(info, f)

    def create_protein(self, job_id, protein_id):
        """Creates folder for new protein

        This method will create a new folder in the input path of the specified
        job in which the density map and ground truth will be stored.

        Parameters
        ----------
        job_id: str
            ID of job for which protein is created

        protein_id: str
            ID of protein which is created

        Raises
        ----------
        JobNotFoundError
            If job with given ID does not exist
        """
        os.makedirs(join(self.get_job_input(job_id), protein_id))

    def get_protein(self, job_id, protein_id):
        """Returns input path of protein

        Parameters
        ----------
        job_id: str
            ID of job containing protein

        protein_id: str
            ID of protein whose path is fetched

        Raises
        ----------
        JobNotFoundError
            If job with given ID does not exist

        ProteinNotFoundError
            If protein with given ID does not exist for job
        """
        if not os.path.isdir(join(self.get_job_input(job_id), protein_id)):
            raise ProteinNotFoundError(job_id, protein_id)

        return join(self.get_job_input(job_id), protein_id)

    def store_density_map(self, job_id, protein_id, density_map):
        """Stores density map of protein

        Parameters
        ----------
        job_id: str
            ID of job containing protein

        protein_id: str
            ID of protein to which density map belongs to

        density_map: FileStorage
            File containing density map of protein

        Raises
        ----------
        JobNotFoundError
            If job with given ID does not exist

        ProteinNotFoundError
            If protein with given ID does not exist for job
        """
        density_map.save(join(self.get_protein(job_id, protein_id), 'density_map.mrc'))

    def store_ground_truth(self, job_id, protein_id, ground_truth):
        """Stores ground truth of protein

        Parameters
        ----------
        job_id: str
            ID of job containing protein

        protein_id: str
            ID of protein to which ground truth belongs to

        ground_truth: FileStorage
            File containing ground truth of protein

        Raises
        ----------
        JobNotFoundError
            If job with given ID does not exist

        ProteinNotFoundError
            If protein with given ID does not exist for job
        """
        ground_truth.save(join(self.get_protein(job_id, protein_id), 'ground_truth.pdb'))

    def store_thresholds(self, job_id, thresholds):
        """Stores thresholds used for job

        This is not specific for each protein. One thresholds file should
        contain thresholds for all proteins of the job.

        Parameters
        ----------
        job_id: str
            ID of job to which thresholds belong to

        thresholds: JSON
            JSON object containing thresholds for each protein

        Raises
        ----------
        JobNotFoundError
            If job with given ID does not exist
        """
        with open(join(self.get_job_input(job_id), 'thresholds.json'), 'w') as f:
            json.dump(thresholds, f)

    def get_thresholds(self, job_id):
        """Returns path to thresholds file

        Parameters
        ----------
        job_id: str
            ID of job for which thresholds are fetched

        Raises
        ----------
        JobNotFoundError
            If job with given ID does not exist

        FileNotFoundError
            If thresholds file does not exist for job
        """
        if not os.path.isfile(join(self.get_job_input(job_id), 'thresholds.json')):
            raise FileNotFoundError('Thresholds file not found')

        return join(self.get_job_input(job_id), 'thresholds.json')

    def create_results(self, job_id):
        """Creates results.zip file

        Compresses the output folder to a results.zip file which is stored in
        the root folder of the job.

        This method does not check the contents of the output folder. It simply
        compresses whatever is contained in it.

        Parameters
        ----------
        job_id: str
            ID of job for which results.zip file is created

        Raises
        ----------
        JobNotFoundError
            If job with given ID does not exist
        """
        zipf = zipfile.ZipFile(join(self.get_job(job_id), 'results.zip'), 'w', zipfile.ZIP_DEFLATED)

        for root, dirs, files in os.walk(self.get_job_output(job_id)):
            for file in files:
                zipf.write(join(root, file),
                           join(root[len(self.get_job_output(job_id)):], file))

        zipf.close()

    def get_results(self, job_id):
        """Returns path to results.zip file

        Parameters
        ----------
        job_id: str
            ID of job for which results are fetched

        Raises
        ----------
        JobNotFoundError
            If job with given ID does not exist

        FileNotFoundError
            If results.zip file does not exist for job
        """
        if not os.path.isfile(join(self.get_job(job_id), 'results.zip')):
            raise FileNotFoundError('Results file not found')

        return join(self.get_job(job_id), 'results.zip')
