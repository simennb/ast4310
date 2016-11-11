from pylab import *
import numpy as np
import astropy as ap
from astropy import units
rc('font',**{'family':'serif'}) 

from ssa2 import *
from ssa_hydrogen import *

if __name__=='__main__':
    temp = np.arange(2000,12001,100)
    dNCadT = np.zeros(temp.shape)
    dNHdT = np.zeros(temp.shape)
    dT = 1.
    chiion = np.array([6,11,50,67])
    
    for i in range(101):
        NCa = sahabolt_E(temp[i],1e2,2,1,chiion)
        NCa2 = sahabolt_E(temp[i]-dT,1e2,2,1,chiion)
        dNCadT[i] = (NCa - NCa2)/(dT*NCa)
        NH = sahabolt_H(temp[i],1e2,2)
        NH2 = sahabolt_H(temp[i]-dT,1e2,2)
        dNHdT[i] = (NH-NH2)/(dT*NH)

    # Plotting
    figure()
    plot(temp,np.absolute(dNHdT), label=r'H')
    plot(temp,np.absolute(dNCadT), label=r'Ca$^+$K')
    yscale('log')
    ylim(1e-9,1)
    grid('on')
    xlabel(r'temperature $T/K$', size=15)
    ylabel(r"$\left| \left( \Delta n(r,s) / \Delta T \right) / n(r,s) \right|$", size=18)
#    show()

    NCa = np.zeros(temp.shape)
    NH = np.zeros(temp.shape)
    for i in range(101):
        NCa[i] = sahabolt_E(temp[i],1e2,2,1,chiion)
        NH[i] = sahabolt_H(temp[i],1e2,2)
    
    plot(temp,NH/np.amax(NH),ls='--',label='rel. pop. H')
    plot(temp,NCa/np.amax(NCa),ls='--',label=r'rel. pop. Ca$^+$')
    legend(loc=4, fontsize=12)
    savefig('temp_sens_Ca_Halpha.png')
    show()
                            
