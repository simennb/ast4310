from ssb import *

path = '../figures/task2/'
fs = 16

ssb = SSB()
wav,F, F_c, I, I_c = ssb.read_solspect()

###############################################################################
# Plot the four spectral distributions together in one figure lambda [0,2]
wave = wav.value
figure()
plot(wav[wave<2.1],F[wave<2.1],lw=2,label=r'$F_{\lambda}$')
plot(wav[wave<2.1],F_c[wave<2.1],lw=2,label=r"$F_{\lambda}'$")
plot(wav[wave<2.1],I[wave<2.1],lw=2,label=r'$I_{\lambda}$')
plot(wav[wave<2.1],I_c[wave<2.1],lw=2,label=r"$I_{\lambda}'$")

xlabel('wavelength [micron]',size=fs)
ylabel(r'spectral distributions',size=fs)
title(r'Spectral distributions per $\Delta\lambda$',size=fs)
legend(loc='best',fontsize=15)
grid('on')
savefig(path+'solspect_lambda_spec_dist.png')
print 'max(I_c) = %.3e, at lambda = %.3f\n'%(max(I_c.value),wave[argmax(I_c.value)])

################################################################################
# Convert spectral distributions into values per frequency bandwidth
nu = (ssb.c/wav).to('Hz')
F_nu  = F/nu*wav
F_cnu = F_c/nu*wav
I_nu  = I/nu*wav
I_cnu = I_c/nu*wav

figure()
plot(wav[wave<2.1],F_nu[wave<2.1],lw=2,label=r'$F_{\nu}$')
plot(wav[wave<2.1],F_cnu[wave<2.1],lw=2,label=r"$F_{\nu}'$")
plot(wav[wave<2.1],I_nu[wave<2.1],lw=2,label=r'$I_{\nu}$')
plot(wav[wave<2.1],I_cnu[wave<2.1],lw=2,label=r"$I_{\nu}'$")

xlabel('frequency [Hz]',size=fs)
ylabel(r'spectral distributions',size=fs)
title(r'Spectral distributions per $\Delta\nu$',size=fs)
legend(loc='best',fontsize=15)
grid('on')
savefig(path+'solspect_lambda_spec_dist_nu.png')
print 'max(I_c) = %.3e, at nu = %.3f\n'%(max(I_cnu.value),wave[argmax(I_cnu.value)])


################################################################################
# Trying to fit a Planck function  to the solar continuum intensity
temp = linspace(4000,8000,200)*u.K
min_i = 0
delta = sum(abs(I_c[wave<2.1])**2) # starting value for finding optimal temperature
planck_best = 0
for i in range(len(temp)):
    planck = ssb.planck(temp[i],wav[wave<2.1])
    #planck *= 1e-10   # scaling in units of 1e10 erg

    delta_current = sum(abs(planck-I_c[wave<2.1])**2)
    if delta_current < delta:
        min_i = i
        delta = delta_current
        planck_best = planck
print 'Best fit temperature, T[%d] = %.3f\n' %(min_i, temp[min_i].value)

figure()
plot(wav[wave<2.1],I_c[wave<2.1],lw=2,label=r"$I_{\lambda}'$")
plot(wav[wave<2.1],planck_best.to(I_c.unit),lw=2,label=r'Planck, T = %.1f K'%temp[min_i].value)

xlabel('wavelength [micron]',size=fs)
ylabel('spectral distributions',size=fs)
title(r'$B_{\lambda}$ fitted to the solar continuum intensity',size=fs)
legend(loc='best',fontsize=15)
grid('on')
savefig(path+'solspect_h_planck.png')

################################################################################
# Inverting the Planck function
T_b = ssb.brightness_temp(I_c,wav)
figure()
plot(wav[wave<2.1],T_b[wave<2.1],lw=2,label='brightness temperature')
xlabel('wavelength [micron]',size=fs)
ylabel('brightness temperature [K]',size=fs)
title('Brightness temperature versus wavelength',size=fs)
grid('on')
savefig(path+'solspect_h_Tb.png')

show()
