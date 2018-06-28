# coding: utf-8
from astropy.io import ascii
import numpy as np
from itertools import product
from matplotlib import pyplot as plt

d = ascii.read('results.txt')

models = ['pwl_ebl_cutoff', 'logp_ebl']
sites = ['North', 'South']
livetime = ['5h', '20h']
emin = ['0.032TeV', '0.050TeV', '0.079TeV', '0.316TeV', '0.501TeV', '0.794TeV']

for m,s,l,e in product(models, sites, livetime, emin):
    mask = []
    for key in d.keys():
        if m in key and s in key and l in key and e in key and key.startswith('Detect'):
            mask.append(d[key])
    mask = np.array(mask)
    d1 = d[~np.isnan(mask)[0]]  # in a given config, all *simulated* sources, detected or not
    m2 = np.where(mask[0] == 1.)  # ~np.isnan(np.array(mask))
    d2 = d[m2]  # keep only detected sources
    print('[{} {} {} {}] Number of detected sources: {}'.format(m, s, l, e, len(d2)))

    if e == emin[0]:
        fig = plt.figure()
        ax = fig.subplots(1,1)
        ax.set_xlabel('Redshift')
        ax.set_ylabel('Source counts')
        ax.hist(d1['Redshift'], color='b', label='All simulated sources')
        ax.hist(d2['Redshift'], color='r', label='Detected sources')
        plt.legend()
        filename = './fig/histz_{}_{}_Eth{}_{}.png'.format(s, l, e, m)
        plt.savefig(filename)



# for m,s,l in product(models, sites, livetime):
#     z = []
#     for key in d.keys():
#         if m in key and s in key and l in key and emin[0] in key and key.startswith('Redshift'):
#             z.append(d['Redshift'])
