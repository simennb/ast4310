from pylab import *
import numpy as np
import astropy as ap
from astropy import units
rc('font',**{'family':'serif'}) 

from ssa2 import *
from ssa_hydrogen import *

if __name__=='__main__':
    chiion = np.array([6,11,50,67])
    # use our functions for Schadeenium as Calsium
    temp = np.arange(1000,20001,100)

    CaH = zeros(temp.shape)
    Caabund = 2.0e-6
    for i in range(0,191):
        NCa = sahabolt_E(temp[i],1e2,2,1,chiion) # is equal to sahabolt_Ca
        
        NH = sahabolt_H(temp[i],1e2,2)
        CaH[i] = NCa*Caabund/NH

    print 'Ca/H ratio at 5000 K = ', CaH[np.argwhere(temp==5000)][0][0]

    # Plotting
    plot(temp,CaH, label=r'strength ratio Ca$^+$K / H$\alpha$')
    yscale('log')
    grid('on')
    xlabel(r'temperature $T / K$',size=14)
    ylabel(r'Ca II K / H$\alpha$',size=14)
    legend(fontsize=14)
    savefig('ssa_Ca_Halpha.png')
    show()

