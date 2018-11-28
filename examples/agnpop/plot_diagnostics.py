# coding: utf-8
import numpy as np
from astropy.io import fits, ascii
from matplotlib import pyplot as plt

data = ascii.read('results.txt')
fig=plt.figure()
ax=fig.subplots(2,1)

sigmas = []
slp = []
spl = []

for key in data.keys():
    if key.startswith('Sigma'):
        for s in data[key]:
            sigmas.append(s)
        if 'logp_ebl' in key:
            for s in data[key]:
                slp.append(s)
        elif 'pwl_ebl' in key:
            for s in data[key]:
                spl.append(s)

sigmas = np.array(sigmas)
slp = np.array(slp)
spl = np.array(spl)

sigmas = sigmas[~np.isnan(sigmas)]
slp = slp[~np.isnan(np.array(slp))]
spl = spl[~np.isnan(spl)]

#plt.title('All sites, all livetimes, all energy thresholds')

ax[0].set_xlabel('Sigma')
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].hist(sigmas,bins=1000, label='All')
ax[0].set_title('All sites, all livetimes, all energy thresholds')


ax[1].set_xlabel('Sigma')
ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].hist(spl, bins=500, label='Power-law+EBL')
ax[1].hist(slp, bins=500, label='Log-parabola+EBL')
ax[1].legend()
