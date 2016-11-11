from pylab import *
from scipy import special
import numpy as np
from ssa_3_1 import *

#########################
# Voigt
def voigt(gamma,x):
    z = (x+1j*gamma)
    V = special.wofz(z).real
    return V

if __name__=='__main__':
    u = np.arange(-10,10.01,0.1)
    a = np.array([0.001,0.01,0.1,1])
    vau = zeros((a.shape[0],u.shape[0]))

    for i in range(4):
        vau[i,:] = voigt(a[i],u[:])
        plot(u[:], vau[i,:], label='a = ' + np.str(a[i]))

    xlim(-10,10)
    legend(fontsize=14)
    grid('on')
    xlabel('u', size=16)
    ylabel('voigt profile', size=16)
    savefig('ssa_voigt.png')

    # ylog
    ylabel('logarithmic voigt profile', size=16)
    yscale('log')
    savefig('ssa_voigt_ylog.png')
    show()

    ##########################
    # Schuster-Schwarzchild line profile

    Ts = 5700.       # solar surface temperature
    T1 = 4200.       # solar T-min temperature = 'reversing layer'
    a = 0.1          # damping parameter
    wav = 5000.0e-8  # wavelength in cm
    tau0 = 1.        # reversing layer thickness at line center
    u = arange(-10,10.0,0.1)
    intensity = zeros(u.shape)

    for i in range(200):
        tau = tau0 * voigt(a,u[i])
        intensity[i] = planck(Ts,wav)*exp(-tau) + planck(T1,wav)*(1.-exp(-tau))

    plot(u,intensity)
    grid('on')
    xlabel('u', size=15)
    ylabel('intensity', size=15)
    savefig('ssa_3_3_1.png')
    show()


    logtau0 = arange(-2,2.1,0.5)
    for itau in range(9):
        for i in range(200):
            tau = 10.**(logtau0[itau])* voigt(a, u[i])
            intensity[i] = planck(Ts,wav)*exp(-tau) + planck(T1,wav)*(1.-exp(-tau))
        plot(u,intensity, label=r'$\log{(\tau_0)} = $'+np.str(logtau0[itau]))

    legend(loc=3, fontsize=13)
    xlabel('u', size=15)
    ylabel('intensity', size=15)
    grid('on')
    savefig('ssa_3_3_2.png')
    show()

    # different wavelengths
    for iwav in range(1,4):
        wav = (iwav**2+1.)*1.0e-5  # wav = 2000, 5000, 10000 angstrom
        
        for itau in range(8):
            for i in range(200):
                tau = 10.**(logtau0[itau]) * voigt(a,u[i])
                intensity[i] = planck(Ts,wav)*exp(-tau) + planck(T1,wav)*(1.-exp(-tau))
            intensity = intensity/intensity[0]
            plot(u,intensity[:], linewidth=1.)
        grid('on')
        xlabel('u',size=15)
        ylabel('intensity',size=15)
        title(r'$\lambda_0 =$ %d Angstrom'%(wav*1e8),size=16)
        savefig('ssa_3_3_wave_%d.png'%iwav)
        show()
