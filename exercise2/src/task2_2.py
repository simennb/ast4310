from ssb import *
# np, ap, u and pylab *

# some plotting parameters
path = '../figures/task2/'
fs = 16  # font size

ssb = SSB()

h,tau5,colm,temp,vturb,nhyd,nprot,nel,ptot,pgasptot,dens = ssb.read_falc()
nhel = 0.1*nhyd 
wav,F, F_c, I, I_c = ssb.read_solspect()
wave = wav.value

################################################################################
'''
temp_b = ssb.brightness_temp(I_c,wav)
kappa = ssb.exthmin(wav, temp_b, nel[h==0])

plot(wav[wave<2.1],kappa[wave<2.1],lw=2)
xlabel('wavelength [micron]',size=fs)
ylabel(r'H$^-$ extinction [dyn cm$^{-2}$]',size=fs)
title('Extinction profile versus wavelength',size=fs)
grid('on')
#ylim([1.0e-25,1.8e-25])
savefig(path+'exthmin_lam_kappa.png')
'''

################################################################################
# height dependent extinction
'''
kappa = ssb.exthmin(wav[wave==0.5],temp,nel)
# Make it not per neutral hydrogen
n_neutral = nhyd-nprot
kappa *= n_neutral

figure()
plot(h,kappa,lw=2)
xlabel('height [km]',size=fs)
ylabel('extinction ',size=fs)
title(r'Height dependent extinction, $\lambda=0.5\mu m$',size=fs)
yscale('log')
grid('on')
xlim([-100,2100])
savefig(path+'exthmin_2.png')
'''

################################################################################
# Thompson scattering
'''
kappa2 = kappa.value + (ssb.sigma_T*nel).value

figure()
plot(h,kappa,lw=2,label=r'H$^-$')
plot(h,ssb.sigma_T*nel,lw=2,ls='--',label='Thompson')
plot(h,kappa2,lw=2,label='total')
xlabel('height [km]',size=fs)
ylabel('extinction ',size=fs)
title(r'Height dependent extinction, $\lambda=0.5\mu m$',size=fs)
legend(loc='best',fontsize=16)
yscale('log')
xlim([-100,2100])
ylim([1e-18,1e-5])
grid('on')
savefig(path+'exthmin_3.png')
'''

################################################################################
# Integration
'''
tau = ssb.optical_depth(h, temp, nel, n_neutral)

figure()
plot(h,tau5,lw=2,label=r'$\tau_{500}$ from FALC')
plot(h,tau,lw=2,label=r'$\tau_{500}$ from H$^-$ and Thompson extinction')
xlabel('height [km]',size=fs)
ylabel('opacity',size=fs)
title('FALC',size=fs)
legend(loc='best',fontsize=14)
yscale('log')
grid('on')
savefig(path+'exthmin_h_tau.png')
'''

################################################################################
# Emergent intensity and height of formation
'''
wl = 0.5*u.micron
intt, contfunc, hmean = ssb.emergent_intensity(h,tau5,temp,nhyd,nprot,nel,wl)

print ('computed continuum intensity wl =%g : %g erg s-1 cm-2 ster-1 cm-1'%(wl.value, intt.value))
w = np.where(wav.value == wl.value)
print ('observed continuum intensity wav=%g : %g erg s-1 cm-2 ster-1 cm-1'%(wav[w].value, I_c[w].value*1e4))
'''

################################################################################
# Plot one wl
'''
figure()
contfunc /= max(contfunc.value)
plot(h,contfunc,lw=2)
plot([hmean.value], [contfunc[argmin(abs((h-hmean).value))].value], marker='o',ms=9,color='blue')
xlabel('height [km]',size=fs)
ylabel('contribution function',size=fs)
title('FALC')
grid('on')
xlim([-100,500])
annotate('mean height of formation = %.1f km'%hmean.value, (10,0.83),size=15)
savefig(path+'emergent_intensity_1.png')
'''

################################################################################
# Plot multiple wl
'''
cl = ['blue','red','green','cyan']
figure()

for wl,i in zip([0.5, 1.0, 1.6, 5.0],range(4)):
    intt, contfunc, hmean = ssb.emergent_intensity(h,tau5,temp,nhyd,nprot,nel,wl*u.micron)
    contfunc /= max(contfunc.value)
    plot(h, contfunc/max(contfunc.value),lw=2,color=cl[i],label=r'$\lambda$=%.1f $\mu$m, $h_{mean}$=%.1f km'%(wl,hmean.value))
    plot([hmean.value], [contfunc[argmin(abs((h-hmean).value))].value], marker='o',ms=9,color=cl[i])

xlabel('height [km]',size=fs)
ylabel('contribution function',size=fs)
title('FALC',size=fs)
legend(loc='best',fontsize=14)
xlim([-100,500])
grid('on')
savefig(path+'emergent_intensity_2.png')
'''


################################################################################
# Disk-center intensity
'''
em_int = np.zeros(len(wave))
for i in range(len(wave)):
    intt, contfunc, hmean = ssb.emergent_intensity(h,tau5,temp,nhyd,nprot,nel,wav[i])
    em_int[i] = intt.value

figure()
plot(wav,em_int/1e14,lw=2,label='computed from FALC')
plot(wav,I_c/1e10,lw=2,label='observed (Allen 1978)')
xlabel(r'wavelength [$\mu$m]',size=fs)
ylabel(r'intensity [$10^{14}$ erg s$^{-1}$ cm$^{-2}$ ster$^{-1}$ cm$^{-1}$]',size=fs)
title('Observed and computed continuum intensity',size=fs)
legend(loc='best',fontsize=15)
grid('on')
savefig(path+'disc_center.png')
'''

################################################################################
# Limb darkening
'''
figure()
mu = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
int_mu = np.zeros(len(wave))
for i in np.arange(len(mu)):
    for j in np.arange(len(wave)):
        intt, contfunc, hmean = ssb.emergent_intensity(h,tau5,temp,nhyd,nprot,nel,wav[j],mu[i])
        int_mu[j] = intt.value
    
    plot(wav,int_mu,lw=2,label=r'$\mu$=%.1f'%mu[i])
    
xlabel(r'wavelength [$\mu$m]',size=fs)
ylabel(r'intensity [erg s$^{-1}$ cm$^{-2}$ ster$^{-1}$ cm$^{-1}$]',size=fs)
title('Limb darkening in FALC',size=fs)
legend(loc='best',fontsize=15)
grid('on')
savefig(path+'limb_darkening_1.png')
'''

################################################################################
# Plot 2
'''
figure()
wav = array([0.2, 0.5, 1.6, 2.5, 5.0])*u.micron
sin_theta = h/ap.constants.R_sun
theta = arcsin(sin_theta)
mu = cos(theta)

int_mu = np.zeros(len(mu))
for i in range(len(wav.value)):
    intt, contfunc, hmean = ssb.emergent_intensity(h,tau5,temp,nhyd,nprot,nel,wav[i],mu)
    int_mu = intt.value
    
    plot(wav,int_mu,lw=2,label=r'$\lambda$=%.1f $\mu$m'%wav[i].value)
    
xlabel(r'$\sin(\theta) = r/R_{\odot}$',size=fs)
ylabel(r'intensity [erg s$^{-1}$ cm$^{-2}$ ster$^{-1}$ cm$^{-1}$]',size=fs)
title('Limb darkening in FALC',size=fs)
legend(loc='best',fontsize=15)
grid('on')
savefig(path+'limb_darkening_2.png')
'''

################################################################################
# Flux integration
'''
fluxspec = ssb.flux_integration(h,tau5,temp,nhyd,nprot,nel, wav)

figure()
plot(wav,fluxspec*1e-14,lw=2,label='computed from FALC')
plot(wav,F_c*1e-10,lw=2,label='observed (Allen 1978')
xlabel('wavelength [micron]',size=15)
ylabel(r'astrophysical flux [10$^{14}$ erg s$^{-1}$ cm$^{-2}$ ster$^{-1}$ cm$^{-1}$]',size=15)
title('Observed and computed continuum flux',size=fs)
grid('on')
legend(loc='best',fontsize=15)
savefig(path+'flux.png')
'''

show()
