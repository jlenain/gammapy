# coding: utf-8
from astropy.io import ascii
import numpy as np
from itertools import product
from matplotlib import pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u

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

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='aitoff')
        c1 = SkyCoord(ra=d1['RAJ2000']*u.degree, dec=d1['DEJ2000']*u.degree, frame='icrs')
        l1 = c1.galactic.l.wrap_at(180*u.degree)
        b1 = c1.galactic.b.wrap_at(180*u.degree)
        ax.plot(-l1.radian,  # East to the left
                b1.radian,
                'o',
                color='b',
                label='All simulated sources')
        c2 = SkyCoord(ra=d2['RAJ2000']*u.degree, dec=d2['DEJ2000']*u.degree, frame='icrs')
        l2 = c2.galactic.l.wrap_at(180*u.degree)
        b2 = c2.galactic.b.wrap_at(180*u.degree)
        ax.plot(-l2.radian,  # East to the left
                b2.radian,
                'x',
                color='r',
                label='Detected sources')
        plt.legend(loc='upper right', bbox_to_anchor=(1., 1.3))
        ax.grid(True)
        filename = './fig/skyplot_{}_{}_Eth{}_{}.png'.format(s, l, e, m)
        plt.savefig(filename)



# for m,s,l in product(models, sites, livetime):
#     z = []
#     for key in d.keys():
#         if m in key and s in key and l in key and emin[0] in key and key.startswith('Redshift'):
#             z.append(d['Redshift'])
