from ssa_3_4 import *

u = np.arange(-200,200.4,0.4)
tau0 = np.logspace(-2,4,61)
a = 0.1
eqw = np.zeros(tau0.size)

for i in range(61):
    intensity = profile(a,tau0[i],u)
    reldepth = (intensity[0]-intensity) / intensity[0]
    eqw[i] = sum(reldepth)*0.4

plot(tau0,eqw)
xlabel(r'$\tau_0$',size=16)
ylabel(r'equivalent width $W_{\lambda}$',size=15)
xscale('log')
yscale('log')
grid('on')
savefig('ssa_3_5.png')
show()
