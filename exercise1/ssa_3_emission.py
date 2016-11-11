from ssa_3_3 import *

def profile(a, tau0, u):
    Ts = 4200.
    T1 = 5700.
    wav = 5000.0e-8
    usize = u.size
    intensity = np.zeros(usize)
    
    for i in range(usize):
        tau = tau0 * voigt(a, u[i])
        intensity[i] = planck(Ts,wav)*exp(-tau) + planck(T1,wav)*(1.-exp(-tau))

    return intensity

if __name__=='__main__':
    ##########################
    # Schuster-Schwarzchild line profile

    Ts = 4200.       # solar surface temperature
    T1 = 5700.       # solar T-min temperature = 'reversing layer'
    a = 0.1          # damping parameter
    wav = 5000.0e-8  # wavelength in cm
    tau0 = 1.        # reversing layer thickness at line center
    u = arange(-10,10.0,0.1)
    intensity = zeros(u.shape)

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
    savefig('ssa_3_5_1.png')
    show()

    # curve of growth
    u = np.arange(-200,200.4,0.4)
    tau0 = np.logspace(-2,4,61)
    a = 0.1
    eqw = np.zeros(tau0.size)

    for i in range(61):
        intensity = profile(a,tau0[i],u)
        reldepth = (intensity[0]-intensity) / intensity[0]
        eqw[i] = sum(reldepth)*0.4

    plot(tau0,abs(eqw))
    xlabel(r'$\tau_0$',size=16)
    ylabel(r'equivalent width $W_{\lambda}$',size=15)
    xscale('log')
    yscale('log')
    grid('on')
    savefig('ssa_3_5_2.png')
    show()

