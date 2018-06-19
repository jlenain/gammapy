# coding: utf-8
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt

data = fits.open('results.fits')
data=data[1].data
fig=plt.figure()
ax=fig.subplots(2,2)
ax[0][0].set_xlabel('Sigma')
ax[0][0].set_xscale('log')
ax[0][0].set_yscale('log')
ax[0][0].hist(data['Sigma'],bins=1000, label='All')
ax[0][0].set_title('All sites, all livetimes, all energy thresholds')
mask=np.where(data['AboveEthFlag'] == True)
d2=data[mask]
ax[0][0].hist(d2['Sigma'],bins=1000, label='Above Eth only')
ax[0][0].legend()

m1=(d2['Model_Name']=='logp_ebl')
m2=(d2['Model_Name']=='pwl_ebl')
dpl=d2[m1]
dlp=d2[m2]
ax[0][1].set_xlabel('Sigma (power-law+EBL)')
ax[0][1].set_ylabel('Sigma (log-parabola+EBL)')
ax[0][1].plot(dpl['Sigma'],dlp['Sigma'],'bo',label='Above Eth only')
ax[0][1].plot([np.min(dpl['Sigma']),np.max(dpl['Sigma'])],[np.min(dpl['Sigma']),np.max(dpl['Sigma'])],'g--')
ax[0][1].set_title('All sites, all livetimes, all energy thresholds')

d3 = fits.open('results_test_all_sp.fits')
d3 = d3[1].data
ax[1][0].set_xlabel('Sigma (gammapy 0.7.dev5114)')
ax[1][0].set_ylabel('Sigma (gammapy 0.8.dev5926)')
ax[1][0].plot(d3['Sigma'], data['Sigma'], 'bo')
ax[1][0].plot([np.min(d3['Sigma']),np.max(d3['Sigma'])],[np.min(d3['Sigma']),np.max(d3['Sigma'])],'g--')

ax[1][1].set_xlabel('Excess (gammapy 0.7.dev5114)')
ax[1][1].set_ylabel('Excess (gammapy 0.8.dev5926)')
ax[1][1].plot(d3['Excess'], data['Excess'], 'bo')
ax[1][1].plot([np.min(d3['Excess']),np.max(d3['Excess'])],[np.min(d3['Excess']),np.max(d3['Excess'])],'g--')

