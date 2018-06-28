"""
CTA AGN Population simulations


From: https://forge.in2p3.fr/projects/cta_science_phys/wiki/CTA-AGN-Population#Analysis-cross-checks

Each analysis group is requested to provide the pre-trials significance (or TS) for each of these cases:
- Observation times 5h and 20h
- For each CTA site (if observable)
- Test 6 different analysis energy thresholds: 30, 50, 100, 300, 500, 1000 GeV. Note only the thresholds above the expression Daniel calculated (see 'GetThreshold.C' at the end of this wiki) will be treated as trials, but all significance values should be provided pre-trials.
- 3 different extrapolations: Power-law + EBL, Power-law + Exponential cutoff at 1/(1+z) TeV + EBL and log-parabola + EBL. Note the spectral parameters provided describe the intrinsic spectra. For the first two cases, the 'PL_*' columns should be used. For the third model, 'LP_*' columns should be used.
"""

from astropy import units as u
from astropy.io import ascii
#from astropy.units import Quantity as Q
from astropy.table import Table, Column

import sys, os
import numpy as np
#import matplotlib.pyplot as plt
import collections
from tqdm import tqdm


from gammapy.scripts import CTAPerf
from gammapy.utils.modeling import Parameter, ParameterList
from gammapy.spectrum.models import AbsorbedSpectralModel, Absorption, SpectralModel, PowerLaw, LogParabola
from gammapy.scripts.cta_utils import CTAObservationSimulation, Target, ObservationParameters

################################
######## Configuration #########
################################

outputtype = 'both'  # fits or ascii or both
verbose = False

# Generate a ~numpy.random.RandomState object
myrand = np.random.RandomState()

# EBL
ebl_model = 'dominguez'

# Cut-off
cut_off_energy = 1. * u.TeV

# Observation parameters
alpha = 0.2 * u.Unit('')
livetimes = [5. * u.h, 20. * u.h]
sites = ['North','South']
sites_lat = dict()
sites_lat['North']= 29. * u.deg
sites_lat['South']= -25. * u.deg


#EminTab = [30., 50., 100., 300., 500., 1000.] * u.GeV
EminTab = [0.03, 0.05, 0.1, 0.3, 0.5, 1.] * u.TeV
emax = 100 * u.TeV

ModelsTab = ['pwl_ebl','pwl_ebl_cutoff','logp_ebl']

################################
################################
################################

class AbsorbedSpectralModelExpoCutOff(SpectralModel):
    """
    Class to handle any spectra with a cut-off
    """

    def __init__(self, spectral_model, cut_off):
        self.spectral_model = spectral_model
        self.cut_off = cut_off

        param_list = []
        for param in spectral_model.parameters.parameters:
            param_list.append(param)

        # Add parameter to the list
        par = Parameter('cut_off', cut_off,
                        parmin=10 * u.MeV, parmax=emax,
                        frozen=True)
        param_list.append(par)

        self.parameters = ParameterList(param_list)


    def evaluate(self, energy, **kwargs):
        """Evaluate the model at a given energy."""
        flux = self.spectral_model(energy=energy)
        absorption = np.exp(- energy / (self.cut_off.to(energy.unit)))
        return flux * absorption


def recursively_default_dict():
    return collections.defaultdict(recursively_default_dict)

ETHRESH0 = 0.02  # energy threshold in TeV at 0 deg zenith angle
def GetThreshold(zen):
    """
    zen: zenith angle in degrees
    """
    z2 = zen*zen
    z4 = z2*z2
    res = 1.2e-4*z2 + 3.5e-8*z4
    res = ETHRESH0*np.power(10.,res)
    return res



###############################################################
########################## START ##############################
###############################################################

catalog = ascii.read("./gll_psch_v12_cross_check.dat")
#print(catalog)

absorption = Absorption.read_builtin(ebl_model)

ctaperf = recursively_default_dict()
for site in ['North','South']:
    for zen in ['20','40']:
        for irftime in ['05h', '50h']:
            #THIS IS IRF, not livetime
            if site == 'North':
                irfname='./irf/irf_'+site+'/CTA-Performance-'+site+'-'+zen+'deg-average-'+irftime+'_20170627.fits.gz'
            else:
                irfname='./irf/irf_'+site+'/CTA-Performance-'+site+'-'+zen+'deg-'+irftime+'_20170627.fits.gz'
            ctaperf[site][zen][irftime] = CTAPerf.read(irfname)

# output = dict()
# output = {'Source_Name' : [],
#           'Redshift': [],
#           'Irf_Site' : [],
#           'Irf_Zen' : [],
#           'Irf_Time' : [],
#           'Livetime' : [],
#           'Emin' : [],
#           'Emin_Real' : [],
#           'AboveEthFlag' : [],
#           'Model_Name' : [],
#           'Bkg' : [],
#           'Excess' : [],
#           'Sigma' : []
# }


for idx, isrc in enumerate(tqdm(catalog)):
    # if isrc['Source_Name'] != 'J2158.8-3013':
    #     continue
    # if idx == 3:
    #     break
    if verbose:
        print(isrc)

    models = dict()
    
    pwlmodel = PowerLaw(amplitude=isrc['PL_Flux_Density']*u.Unit('cm^-2 s^-1 GeV^-1'),
                        reference=isrc['Pivot_Energy']*u.Unit('GeV'),
                        index=isrc['PL_Index']*u.Unit(''))
    
    lpmodel = LogParabola(amplitude=isrc['LP_Flux_Density']*u.Unit('cm^-2 s^-1 GeV^-1'),
                          reference=isrc['Pivot_Energy']*u.Unit('GeV'),
                          alpha=isrc['LP_Index']*u.Unit(''),
                          beta=isrc['LP_Beta']*u.Unit(''))

    redshift = isrc['Redshift'] * u.Unit('')

    # Absorbed spectral model
    """
    abs_pwlmodel = AbsorbedSpectralModel(spectral_model=pwlmodel,
                                         absorption=absorption,
                                         parameter=redshift)

    abs_pwlmodel_cutoff = AbsorbedSpectralModelExpoCutOff(abs_pwlmodel,
                                                          cut_off_energy / (1 + redshift))

    abs_lpmodel = AbsorbedSpectralModel(spectral_model=lpmodel,
                                      absorption=absorption,
                                      parameter=redshift)
    """

    
    models['pwl_ebl'] = AbsorbedSpectralModel(spectral_model=pwlmodel,
                                              absorption=absorption,
                                              parameter=redshift)

    models['pwl_ebl_cutoff'] = AbsorbedSpectralModelExpoCutOff(spectral_model=models['pwl_ebl'],
                                                               cut_off=cut_off_energy / (1 + redshift))

    models['logp_ebl'] = AbsorbedSpectralModel(spectral_model=lpmodel,
                                               absorption=absorption,
                                               parameter=redshift)

    #pwlmodel.plot(energy_range=[0.01, 10] * u.TeV,label='pwl')
    #abs_pwlmodel.plot(energy_range=[0.01, 10] * u.TeV,label='pwl+ebl')
    #abs_pwlmodel_cutoff.plot(energy_range=[0.01, 10] * u.TeV,label='pwl+ebl+cutoff')
    #abs_lpmodel.plot(energy_range=[0.01, 10] * u.TeV,label='lp+ebl')

    #Loops over time (5h, 20h), site (N,S), thresholds (var min, fixed max)

    dec = isrc['DEJ2000']*u.deg

    #Loop over both sites
    for mysite in sites:

        #Determine allowed Irf_zen for each site
        zen_at_culmination_plus = np.abs(dec-sites_lat[mysite])
        # zen_at_culmination_plus += 5 * u.deg

        irfzen = None
        if zen_at_culmination_plus <= 25 * u.deg:
            irfzen = '20'
        elif 25 * u.deg < zen_at_culmination_plus <= 45 * u.deg:
            irfzen = '40'
        else:
            if verbose:
                print('Source {} never above 45 deg. in site {}, skipping it...'.format(isrc['Source_Name'], mysite))
            continue
                    
        #Determine E threshold from zenith angle
        emin_from_zen = GetThreshold(zen_at_culmination_plus.value)*u.TeV

        for mylivetime in livetimes:
            if mylivetime == 5.*u.h:
                irftime = '05h'
            elif mylivetime == 20.*u.h:
                irftime = '50h'
            else:
                irftime = None
            for emin in EminTab:
                AboveEthFlag = True
                if emin < emin_from_zen:
                    AboveEthFlag = False
        
                if verbose:
                    print(mysite,dec,sites_lat[mysite],zen_at_culmination_plus,irfzen,mylivetime,emin)
        
                obs_param = ObservationParameters(alpha=alpha, livetime=mylivetime,
                                                  emin=emin, emax=emax)


                reco_energy = ctaperf[mysite][irfzen][irftime].bkg.energy
                idx_min = np.abs(reco_energy.lo - emin).argmin()
                emin_real = reco_energy.lo[idx_min]


                for mymodel in ModelsTab:
                    sigmadictkey = 'Sigma_{}_{:d}h_Eth{:.3f}TeV_{}'.format(mysite, int(mylivetime.value), emin_real.value, mymodel)
                    excessdictkey = 'Excess_{}_{:d}h_Eth{:.3f}TeV_{}'.format(mysite, int(mylivetime.value), emin_real.value, mymodel)
                    bkgdictkey = 'Bkg_{}_{:d}h_Eth{:.3f}TeV_{}'.format(mysite, int(mylivetime.value), emin_real.value, mymodel)
                    detectdictkey = 'Detected_{}_{:d}h_Eth{:.3f}TeV_{}'.format(mysite, int(mylivetime.value), emin_real.value, mymodel)
                    if sigmadictkey not in catalog.keys():
                         catalog[sigmadictkey] = np.full(len(catalog), np.nan)
                    if excessdictkey not in catalog.keys():
                         catalog[excessdictkey] = np.full(len(catalog), np.nan)
                    if bkgdictkey not in catalog.keys():
                         catalog[bkgdictkey] = np.full(len(catalog), np.nan)
                    if detectdictkey not in catalog.keys():
                         catalog[detectdictkey] = np.full(len(catalog), np.nan)

                    target = Target(name=isrc['Source_Name'],model=models[mymodel])

                    simu = CTAObservationSimulation.simulate_obs(perf=ctaperf[mysite][irfzen][irftime],
                                                                 target=target,
                                                                 obs_param=obs_param,
                                                                 random_state=myrand)  #'random-seed')

                    sigmas = simu.total_stats_safe_range.sigma
                    excess = simu.total_stats_safe_range.excess
                    bkg = simu.total_stats_safe_range.background

                    detected = sigmas>=5. and excess > 10 and excess > 0.05*bkg*alpha
                    
                    # output['Source_Name'].append(isrc['Source_Name'])
                    # output['Redshift'].append(isrc['Redshift'])
                    # output['Irf_Site'].append(mysite)
                    # output['Irf_Zen'].append(float(irfzen))
                    # output['Irf_Time'].append(float(irftime.split('h')[0]))
                    # output['Livetime'].append(mylivetime.value) # pq je ne peux pas mettre une grandeur avec unite?
                    # output['Emin'].append(emin.value)
                    # output['Emin_Real'].append(emin_real.value)
                    # output['Model_Name'].append(mymodel)
                    # output['AboveEthFlag'].append(AboveEthFlag)
                    # output['Bkg'].append(bkg)
                    # output['Excess'].append(excess)
                    # output['Sigma'].append(sigmas)
                    catalog[sigmadictkey][idx] = '{:.2f}'.format(sigmas)
                    catalog[excessdictkey][idx] = '{:.2f}'.format(excess)
                    catalog[bkgdictkey][idx] = '{:.2f}'.format(bkg)
                    catalog[detectdictkey][idx] = detected

# results = Table()
# results['Source_Name'] = Column(output['Source_Name'], description='Source name')
# results['Redshift'] = Column(output['Redshift'], format='{:.3f}', unit='', description='Redshift')
# results['Irf_Site'] = Column(output['Irf_Site'], description='Site of the considered CTA array')
# results['Irf_Zen'] = Column(output['Irf_Zen'], format='{:.1f}', unit='degree', description='Zenith used for the IRF')
# results['Irf_Time'] = Column(output['Irf_Time'], format='{:.1f}', unit='hour', description='Livetime used for the IRF')
# results['Livetime'] = Column(output['Livetime'], format='{:.1f}', unit='hour', description='Observation time')
# results['Eth'] = Column(output['Emin_Real'], format='{:.3f}', unit='TeV', description='Energy threshold')
# results['Model_Name'] = Column(output['Model_Name'], unit='', description='Name of the spectral model')
# results['Sigma'] = Column(output['Sigma'], format='{:.2f}', unit='', description='Significance')

# results['Emin'] = Column(output['Emin'], unit='', description='')
# results['Emin_Real'] = Column(output['Emin_Real'], unit='', description='')
# results['AboveEthFlag'] = Column(output['AboveEthFlag'], unit='', description='')
## results['SigmaErr'] = Column(src_sigma_err, unit='', description='')
# results['Excess'] = Column(output['Excess'], unit='', description='')
## results['ExcessErr'] = Column(src_excess_err, unit='', description='')
# results['Bkg'] = Column(output['Bkg'], unit='', description='')
## results['BkgErr'] = Column(src_bkg_err, unit='', description='')

# Save results
outdir = './'                  
if not os.path.exists(outdir):
    os.makedirs(outdir)

if outputtype == 'ascii':
    filename = 'results.txt'
    catalog.write(outdir + filename, format='ascii', overwrite=True)
elif outputtype == 'fits':
    filename = 'results.fits'
    catalog.write(outdir + filename, format='fits', overwrite=True)
elif outputtype == 'both':
    filename = 'results.txt'
    catalog.write(outdir + filename, format='ascii', overwrite=True)
    filename = 'results.fits'
    catalog.write(outdir + filename, format='fits', overwrite=True)
else:
    print('Error: output format \'{}\' not defined, no output written'.format(outputtype))
    sys.exit(1)
    
print('output file: {}'.format(outdir + filename))

# pwlmodel.plot(energy_range=[0.1, 100] * u.TeV)
# plt.legend()
# plt.show()
