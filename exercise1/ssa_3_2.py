from pylab import *

B = 2.

tau = arange(0.01,10.01,0.01)
intensity = zeros(tau.shape)
for I0 in range(4,-1,-1):
    intensity[:] = I0 * exp(-tau[:]) + B*(1-exp(-tau[:]))
    plot(tau, intensity, label='intensity I0 = ' + str(I0))

xlabel(r'optical depth $\tau$', size=15)
ylabel('intensity', size=15)
legend(fontsize=14)
grid('on')
savefig('ssa_3_2.png')

# loglog
xscale('log')
yscale('log')
legend().set_visible(False)
savefig('ssa_3_2_log.png')
show()
