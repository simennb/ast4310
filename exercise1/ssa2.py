from pylab import *
import astropy as ap
from astropy import units as u

def partfunc_E(temp):
    '''
    Calculates the partition function of Schadeenium for a given temperature
    '''
    chiion = [7, 16, 31, 51]  # Schaade ionization energies
    k = 8.61734e-5 # Boltzmann constant in eV/deg
    U = zeros(4)
    for r in range(4):
        for s in range(0, chiion[r]):
            U[r] += exp(-s/(k*temp))
            
    return U

#U = partfunc_E(5000)
#print U

###########################

def boltz_E(temp,r,s):
    '''
    Computes Boltzmann distribution
    '''
    U_r = partfunc_E(temp)
    k = 8.61734e-5 # Boltzmann constant in eV/deg
    relnrs = 1./U_r[r-1]*exp(-(s-1)/(k*temp))

    return relnrs

#for s in range(1,11):
 #   print boltz_E(5000,1,s)


def saha_E(temp,r):
    #need constants
    m_e = 9.10939e-28 # mass of electron in grams
    h = 6.62607e-27   # Planck constant in erg s
    k = 8.61734e-5    # Boltzmann constant in eV/deg
    #need N_e
    
    N_r = zeros(5)
    N_r[0] = 1

    pass