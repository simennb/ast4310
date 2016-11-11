from pylab import *
import numpy as np
import astropy as ap
from astropy import units
rc('font',**{'family':'serif'}) 

from ssa2 import *
from ssa_hydrogen import *

for T in arange(2e3,2e4+1,2e3):
    print T, sahabolt_H(T,1e2,1)

temp = arange(1e3,2e4+1,1e2)
nH = zeros(temp.shape)
for i in range(191):
    nH[i] = sahabolt_H(temp[i],1e2,1)

plot(temp,nH)
xlabel(r'temperature $T/K$',size=16)
ylabel(r'neutral hydrogen fraction',size=16)
grid('on')
axhline(0.5,linestyle='--',color='black')
savefig('ssa_2_10.png')
show()
