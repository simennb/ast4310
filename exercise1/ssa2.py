from pylab import *
import numpy as np
import astropy as ap
from astropy import units
rc('font',**{'family':'serif'}) 

def partfunc_E(temp):
    '''
    Calculates the partition function of Schadeenium for a given temperature
    '''
    chiion = [7, 16, 31, 51]  # Schaade ionization energies
    k = ap.constants.k_B # Boltzmann constant in J/(K)
    U = zeros(4)
    for r in range(4):
        for s in range(0, chiion[r]):
            U[r] += exp(-s*units.eV/(k.to('eV/K')*temp))
            
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
    relnrs = 1./U_r[r-1]*exp(-(s-1)/(k*temp.value))

    return relnrs

#for s in range(1,11):
 #   print boltz_E(5000,1,s)

############################

def saha_E(temp,elpress,ionstage):
    # Some constants
    m_e = ap.constants.m_e # mass of electron in kg
    h = ap.constants.h   # Planck constant in J s
    k = ap.constants.k_B # Boltzmann constant in J/(K)
    
    # kT and electron density
    kevT = k.to('eV/K')*temp
    kergT = k.to('erg/K')*temp
    eldens = elpress/kergT

    # Partition function
    chiion = [7, 16, 31, 51]*units.eV  # Schaade ionization energies in eV
    
    u = partfunc_E(temp)
    u = np.append(u,2.0)
    
    # Saha constant
    sahaconst = (2*pi*m_e.to('g')*kergT/(h.to('erg*s')**2))**1.5*2./eldens

    # N(r)
    nstage = zeros(5)
    nstage[0] = 1

    for r in range(4):
        nstage[r+1] = (nstage[r]*sahaconst*u[r+1]/u[r]*exp(-chiion[r]/kevT)).value
    ntotal = sum(nstage)
    nstagerel = nstage/ntotal

    return nstagerel[ionstage-1]

#for r in range(1,6):
 #   print saha_E(20000*units.K,1e1,r)


##########################################
def sahabolt_E(temp,elpress,ion,level):
    '''
    Computes Saha-Boltzmann population n_(r,s)/N for level r,s of E
    '''
    return saha_E(temp,elpress,ion)*boltz_E(temp,ion,level)


if __name__=='__main__':
    task = '2.5'


    if task=='2.5':
        temp = np.arange(0,30001,1000)*units.K
        pop = zeros((5,31))

        '''
        for T in np.arange(1,31):
            for r in np.arange(1,5):
                pop[r,T] = sahabolt_E(temp[T],131.,r,1)

        # ground-state plot
        figure()
        for i in range(1,5):
            plot(temp,pop[i,:])
        title(r'Ground-state population for $P_e = 131$ dyne cm${}^{-2}$',size=15)
        xlabel(r'Temperature [K]',size=15)
        ylabel(r'Population $n_{r,1}/N$',size=15)
        yscale('log')
        grid('on')
        ylim([1e-3,1.1])
        legend(['ground stage','first ion stage','second ion stage','third ion stage'],loc='best')
        savefig('schadeenium_pop.png')
        show()
        '''

        # adding more lines
        figure()
        grid('on')
        colors = ['b','g','r','k']
        lines = ['-','--','--']
        s = [1,2,4]
        for j in range(len(s)):
            for T in np.arange(1,31):
                for r in np.arange(1,5):
                    pop[r,T] = sahabolt_E(temp[T],131.,r,s[j])

            # plot
            for i in range(1,5):
                plot(temp,pop[i,:],lines[j]+colors[i-1])
        title(r'Level population for $P_e = 131$ dyne cm${}^{-2}$',size=15)
        xlabel(r'Temperature [K]',size=15)
        ylabel(r'Population $n_{r,1}/N$',size=15)
        yscale('log')
        ylim([1e-3,1.1])
#        legend(['ground stage','first ion stage','second ion stage','third ion stage'],loc='best')
        savefig('schadeenium_pop_higher.png')
        show()
        
