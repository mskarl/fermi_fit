import numpy as np
import os
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.integrate import quad
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
MeV_to_erg = 1.60218e-6


def get_luminosity(basepath, emin=100, emax=100000, z=0, verbose=False, z_correction=True):
    bowtie_path = os.path.join(basepath, 'bowtie.npy')
    llh_path = os.path.join(basepath, 'llh.npy')
    L_D = cosmo.luminosity_distance(z).to_value(u.cm)
    if z_correction:
        z_surface_area = 4*np.pi*L_D**2
    else:
        z_surface_area = 1
        z=0
    if (not os.path.exists(bowtie_path)) or (not os.path.exists(llh_path)):
        print('Wrong path or missing file')
        return
    bowtie = np.load(os.path.join(basepath, 'bowtie.npy'),
                     allow_pickle=True, encoding='bytes')[()]
    llh = np.load(os.path.join(basepath, 'llh.npy'),
                  allow_pickle=True, encoding='bytes')[()]
    source = llh['config']['selection']['target']
    if llh['sources'][source]['ts'] < 4:#13.25
        try:
            gamma = llh['sources'][source]['spectral_pars']['Index']['value']
            up_lim = llh['sources'][source]['eflux100_ul95']/(1+z)**(2-gamma) * MeV_to_erg * z_surface_area
        except:
            up_lim =0
        return up_lim, 0, 0
    e2 = 10 ** (2 * bowtie['log_energies'])
    energies = np.array(10 ** bowtie['log_energies'])
    interp = InterpolatedUnivariateSpline(energies, bowtie['dnde'], k=1,ext=0)
    bf=quad(lambda x: interp(x/(1+z))*x/(1+z)**2, emin, emax)[0] * MeV_to_erg * z_surface_area
    interp = InterpolatedUnivariateSpline(energies, bowtie['dnde_lo'], k=1, ext=0)
    low = quad(lambda x: interp(x/(1+z))*x/(1+z)**2, emin, emax)[0] * MeV_to_erg * z_surface_area
    interp = InterpolatedUnivariateSpline(energies, bowtie['dnde_hi'], k=1, ext=0)
    up = quad(lambda x: interp(x/(1+z))*x/(1+z)**2, emin, emax)[0] * MeV_to_erg * z_surface_area
    if verbose:
        try:
            gamma = llh['sources'][source]['spectral_pars']['Index']['value']
            bf2 = llh['sources'][source]['eflux100']/(1+z)**(2-gamma) * MeV_to_erg * z_surface_area
            print('Redshift {}'.format(z))
            print('Luminosity Distance {} cm'.format(L_D))
            print('Gamma : {}'.format(gamma))
            print('Powerlaw K-Correction: {}'.format(np.log10(bf2)))
            print('Integral {}'.format(np.log10(bf)))
        except:
            print('No powerlaw')
            pass
    print('----')
    return bf, low, up


if __name__ == '__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("--basepath", "-b", help="basepath of the Fermi files",
                        type=str)
    parser.add_argument("--emin", help="minimum energy in MeV",
                        type=float, default=100.)
    parser.add_argument("--emax", help="maximum energy in MeV",
                        type=float, default=100000.)
    parser.add_argument("--redshift", "-z", help='redshift of the source',
                        type=float)
    args=parser.parse_args()

    bf, low, up = get_luminosity(args.basepath, args.emin, args.emax, args.redshift)
    print('Best-Fit {}, Lower Limit {}, Upper Limit {}'.format(bf, low, up))

