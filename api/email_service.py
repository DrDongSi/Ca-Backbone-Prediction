"""Module responsible for sending email notifications when prediction jobs have
finished"""

import os
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

__author__ = 'Jonas Pfab'


def send_email(to, download_link):
    """Sends email with download link to given address

    This method will send an email containing the download link for the
    prediction results to the given email. In order for this method to work
    several environment variables have to be set:
        - EMAIL_ADDRESS: Email address from emails will be sent
        - EMAIL_SERVER: SMTP server address
        - EMAIL_SERVER_PORT: Port of SMTP server
        - EMAIL_PASSWORD: Password for email account

    Parameters
    ----------
    to: str
        Email address to which email is sent to

    download_link: str
        Link from which prediction results can be downloaded
    """
    msg = MIMEMultipart()
    msg['From'] = os.environ['EMAIL_ADDRESS']
    msg['To'] = to
    msg['Subject'] = 'CA Backbone Prediction Finished'
    body = 'The result of your CA backbone prediction job can be downloaded from ' + download_link
    msg.attach(MIMEText(body, 'plain'))

    server = get_server()
    text = msg.as_string()
    server.sendmail(os.environ['EMAIL_ADDRESS'], to, text)
    server.quit()


def get_server():
    """Creates and connects to SMTP server"""
    server = smtplib.SMTP(os.environ['EMAIL_SERVER'], int(os.environ['EMAIL_SERVER_PORT']))
    server.connect(os.environ['EMAIL_SERVER'], int(os.environ['EMAIL_SERVER_PORT']))
    server.ehlo()
    server.starttls()
    server.ehlo()
    server.login(os.environ['EMAIL_ADDRESS'], os.environ['EMAIL_PASSWORD'])

    return server
