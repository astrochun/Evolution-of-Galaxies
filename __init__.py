from astropy.io import ascii as asc

from .log_commons import log_stdout


def table_to_dict(tbl_name, log=None):
    """
    Purpose:
        This function takes an astropy ascii table and converts it to a dictionary.

    Parameters:
        tbl_name --> a string containing the full file path of the table.

    Returns:
        A dictionary that contains key-value pairs of data column name to data array.
    """

    if log is None:
        log = log_stdout()

    log.info(f"starting ...")

    log.info(f"Reading: {tbl_name}")
    tbl = asc.read(tbl_name)

    log.info("finished.")
    return tbl.to_pandas().to_dict(orient='list')
