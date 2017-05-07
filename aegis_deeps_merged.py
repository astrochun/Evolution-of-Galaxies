from astropy.io import fits, ascii
import numpy as np
import pdb
import matplotlib.pyplot as plt


aegis = ascii.read('../magfiles/aegis5.1_deep2_crossmatch.mag')
deeps = ascii.read('../magfiles/stacked_deep_fields.mag')

aid = np.array(aegis['id'])
did = np.array(deeps['id'])
temp = np.zeros(3889)
new_aegis = np.append(aid, temp)

new_id = np.array([])
for ag, dps in zip(new_aegis, did):
    dp = float(dps)
    pdb.set_trace()
    if ag == dp:
        np.append(new_id, ag)
    else:
        np.append(new_id, dp)

deeps['id'] = new_id

ascii.write(deeps, 'aegis_deeps_merged.py', overwrite = True)        
  
