#!/bin/env python

import numpy as np
import itertools
from matplotlib import pyplot as plt
from astropy.io import fits
from gammapy.stats import significance_on_off

siteTab = ['North', 'South']
modelTab = ['pwl_ebl', 'pwl_ebl_cutoff', 'logp_ebl']
eminTab = [0.03, 0.05, 0.1, 0.3, 0.5, 1.0]
livetimeTab = [5., 20.]

data = fits.open('results.fits')
dataorg = data[1].data

for s, m, e, t in itertools.product(siteTab, modelTab, eminTab, livetimeTab):
    #s, m, e, t='North', 'pwl_ebl', 0.05, 5.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #ax.set_xlabel('Excess')
    ax.set_xlabel('Sigma')
    ax.set_yscale('log')
    ax.set_title('%s: %s, E>%.2f TeV, t=%d h' % (s, m, e, t))
    mask = np.array(dataorg['Irf_Site']==s) & np.array(dataorg['Model_Name']==m) & np.array(dataorg['Emin']==e) & np.array(dataorg['Livetime']==t)
    data = dataorg[mask]

    excess = data['Excess']
    bkg = data['Bkg']
    sigmas = data['Sigma']
    alpha=0.2
    non = excess + alpha * bkg 

    ax.hist(sigmas, bins=100) # label='%s: %s, E>%d GeV, t=%d h' % (s, m, e, t))
    #ax.plot(excess, sigmas, 'o', label='%s: %s, E>%d GeV, t=%d h' % (s, m, e, t))
    # xt = np.logspace(np.log10(np.min(excess)), np.log10(np.max(excess)), num=100, base=10.0)
    # ax.plot(excess, significance_on_off(non,bkg,alpha), 'rx')
    outfig = './fig/sigma_%s_%s_%.2fTeV_%dh.png' % (s, m, e, t)
    fig.savefig(outfig)

# plt.legend()
# plt.show()
