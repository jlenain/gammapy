#!/bin/env python
# coding: utf-8

import numpy as np
import gammapy
from matplotlib import pyplot as plt
from astropy.io import fits

fig=plt.figure()
ax=fig.add_subplot(111)
ax.set_xlabel('Sigma (old gammapy)')
ax.set_ylabel('Sigma (new gammapy)')
ax.set_title('Comparison of sigmas in AGNPop between\nSantiago\'s gammapy version and gamampy {}'.format(gammapy.__version__))

old = fits.open('results_test_all_sp.fits')
old = old[1].data
ax.plot([np.min(old['Sigma']),np.max(old['Sigma'])], [np.min(old['Sigma']),np.max(old['Sigma'])], '--')

new = fits.open('results.fits')
new = new[1].data
ax.plot(old['Sigma'], new['Sigma'], 'bo')

