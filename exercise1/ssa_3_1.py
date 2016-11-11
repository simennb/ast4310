from pylab import *

def planck(temp, wav):
    h = 6.62607e-27  # planck in erg s
    k = 1.38065e-16  # boltzmann in erg/K
    c = 2.99792e10   # speed of light cm/s

    return 2*h*c**2/wav**5*1./(exp(h*c/(wav*k*temp))-1)

if __name__=='__main__':

    wav = arange(1000,20801,200)
    b = zeros(wav.shape)

    figure()
    xlabel(r'wavelength $\lambda / \AA$', size=15)
    ylabel(r'Planck function', size=15)
    grid('on')

    for T in range(8000,5000-1,-200):
        b[:] = planck(T, wav[:]*1e-8)
        plot(wav,b,'-',label='T=%d'%T)

    legend(loc='best')
    savefig('ssa_3_1.png')

    # log y-axis
    yscale('log')
    legend().set_visible(False)
    savefig('ssa_3_1_logy.png')

    # log-log axes
    xscale('log')
    savefig('ssa_3_1_loglog.png')

    show()
