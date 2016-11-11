from pylab import *
import numpy as np
from ssa_3_3 import *

def profile(a, tau0, u):
    Ts = 5700.
    T1 = 4200.
    wav = 5000.0e-8
    usize = u.size
    intensity = np.zeros(usize)
    
    for i in range(usize):
        tau = tau0 * voigt(a, u[i])
        intensity[i] = planck(Ts,wav)*exp(-tau) + planck(T1,wav)*(1.-exp(-tau))

    return intensity

if __name__=='__main__':
    
    # Checking the profile
    u = np.arange(-200,200.4,0.4)
    a = 0.1
    tau0 = 1.0e2
    intensity = profile(a,tau0,u)

    plot(u,intensity)
    xlabel('u',size=16)
    ylabel('intensity',size=16)
    grid('on')
    savefig('ssa_3_4_1.png')
    show()

    # relative
    reldepth = (intensity[0]-intensity)/intensity[0]
    plot(u,reldepth)
    xlabel('u',size=16)
    ylabel('relative intensity',size=16)
    grid('on')
    savefig('ssa_3_4_2.png')

    eqw = sum(reldepth)*0.4
    print eqw
    show()
    
