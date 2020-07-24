from astropy.io import ascii as asc


def table_to_dict(tbl_name):
    '''
    Purpose: 
        This function takes an astropy ascii table and converts it to a dictionary.
        
    Parameters:
        tbl_name --> a string containing the full file path of the table.
        
    Returns:
        A dictionary that contains key-value pairs of data column name to data array.  
    '''
    
    tbl = asc.read(tbl_name)
    cols = tbl.colnames
    
    return {cols[ii]:tbl[cols[ii]].data for ii in range(len(cols))}
