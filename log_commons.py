import sys
from os.path import join
from os import uname

from getpass import getuser
from socket import gethostname
from requests import get

import logging
formatter = logging.Formatter('%(asctime)s - %(levelname)8s: %(message)s', "%H:%M:%S")

file_formatter = logging.Formatter('%(asctime)s %(levelname)8s - %(module)20s %(funcName)30s : %(message)s',
                                   "%H:%M:%S")


class LogClass:
    """
    Main class to log information to stdout and ASCII logfile

    To use:
    log = LogClass(log_dir, logfile).get_logger()

    :param log_dir: Relative path for exported logfile directory
    :param logfile: Filename for exported log file
    """

    def __init__(self, log_dir: str, logfile: str):
        self.LOG_FILENAME = join(log_dir, logfile)

    def get_logger(self):
        file_log_level = logging.DEBUG  # This is for file logging
        log = logging.getLogger("main_logger")
        if not log.handlers:
            log.setLevel(file_log_level)

            sh = logging.StreamHandler(sys.stdout)
            sh.setLevel(logging.INFO)  # Only at INFO level
            sh.setFormatter(formatter)
            log.addHandler(sh)

            fh = logging.FileHandler(self.LOG_FILENAME)
            fh.setLevel(file_log_level)
            fh.setFormatter(file_formatter)
            log.addHandler(fh)

            log.handler_set = True
            log.propagate = False
        return log


def log_stdout():
    log_level = logging.INFO
    log = logging.getLogger("stdout_logger")
    if not log.handlers:
        log.setLevel(log_level)
        sh = logging.StreamHandler(sys.stdout)
        sh.setFormatter(formatter)
        log.addHandler(sh)

        log.handler_set = True
        log.propagate = False
    return log


def get_user_hostname() -> dict:
    """
    Retrieve user, hostname, IP, and OS configuration

    :return sys_info: dictionary with 'user' 'hostname' and 'ip' dict
    """

    sys_info = dict()

    sys_info['user'] = getuser()
    sys_info['hostname'] = gethostname()
    sys_info['ip'] = get('https://api.ipify.org').text

    os_name = uname()
    sys_info['os'] = f"{os_name[0]} {os_name[2]} {os_name[3]}"

    return sys_info
