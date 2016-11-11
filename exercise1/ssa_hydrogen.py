from pylab import *
import numpy as np
import astropy as ap
from astropy import units
rc('font',**{'family':'serif'}) 

def sahabolt_H(temp,elpress,level):
    '''
    Computes Saha-Boltzmann population n_(r,s)/N for level r,s of H
    '''
    m_e = ap.constants.m_e # mass of electron in kg
    h = ap.constants.h   # Planck constant in J s
    k = ap.constants.k_B # Boltzmann constant in J/(K)


    # kT and electron density
    hergs = h.to('erg s')
    kevT = k.to('eV/K')*temp
    kergT = k.to('erg/K')*temp
    elmass = m_e.to('g')
    eldens = elpress/kergT

    # energy levels and weights for hydrogen
    nrlevels = 100                # reasonable partition function cut-off value
    g = zeros((2,nrlevels))       # declaration weights (too many for protons)
    chiexc = zeros((2,nrlevels))  # declaration excit. energies

    for s in range(nrlevels):
        g[0,s] = 2*(s+1)**2                  # statistical weights
        chiexc[0,s] = 13.598*(1-1./(s+1)**2) # excitation weights
    
    g[1,0] = 1.                    # statistical weights free proton 
    chiexc[1,0] = 0.

    # partition functions
    u = zeros(2)
    for s in range(nrlevels):
        u[0] += g[0,s]*exp(-chiexc[0,s]/kevT.value)
    u[1] = g[1,0]

    # Saha
    sahaconst = (2*pi*elmass.value*kergT.value / (hergs.value*hergs.value))**(1.5)*2./eldens.value

    nstage = zeros(2)
    nstage[0] = 1.
    nstage[1] = nstage[0]*sahaconst*u[1]/u[0] * exp(-13.598/kevT.value)
    ntotal = sum(nstage)            # sum both stages = total hydrogen density

    # Boltzmann
    nlevel = nstage[0]*g[0,level-1]/u[0]*np.exp(-chiexc[0,level-1]/kevT.value)
    nlevelrel = nlevel/ntotal       # fraction of total hydrogen density


    # Print for checking
    '''
    for s in range(6):
        print "%2d %5d %7.4f %9.7e" % (s+1, g[0,s], chiexc[0,s], 
                                       g[0,s]*np.exp(-chiexc[0,s]/kevT.value))

    print 
    for s in range(0,nrlevels,10):
        print "%2d %5d %6.4f %8.7e" %(s+1, g[0,s], chiexc[0,s], 
                                      g[0,s]*np.exp(-chiexc[0,s]/kevT.value))
    '''

    return nlevelrel

