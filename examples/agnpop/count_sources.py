# coding: utf-8
from astropy.io import fits
import numpy as np

d = fits.open('results.fits')
d = d[1].data

livetime = 20.
Eth = 0.501187205314636
sigmin = 5.0

for site in ['North', 'South']:
    dN=d[(d['Irf_Site']==site) & (np.isclose(d['Eth'],Eth)) & (np.isclose(d['Livetime'],livetime)) & (d['Sigma']>=sigmin)]
    print('Number of sources detected above {} sigmas in the site {} above Eth={} in {}h: {}'.format(sigmin, site, Eth, livetime, len(np.unique(dN['Source_Name']))))
