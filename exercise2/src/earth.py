from ssb import *
# np, ap, u and pylab *

# some plotting parameters
path = '../figures/earth/' # folder where figures are stored
fs = 16  # font size

ssb = SSB()

h,logP,temp,logrho,logN = ssb.read_earth()


################################################################################
# Plot temperature against height
'''
plot(h,temp,lw=2)
xlabel('height [km]',size=fs)
ylabel('temperature [K]',size=fs)
title('Earth atmosphere',size=fs)
grid('on')
savefig(path+'earth_h_temp.png')
'''

################################################################################
# Particle density against height
'''
figure()
plot(h,10**(logN),lw=2)
xlabel('height [km]',size=fs)
ylabel(r'particle density [cm$^{-3}$]',size=fs)
title('Earth atmosphere',size=fs)
grid('on')
yscale('log')
savefig(path+'earth_h_N.png')
'''

################################################################################
# Gas density against height
'''
figure()
plot(h,10**(logrho),lw=2)
xlabel('height [km]',size=fs)
ylabel('gas density [g cm$^{-3}$]',size=fs)
title('Earth atmosphere',size=fs)
grid('on')
yscale('log')
savefig(path+'earth_h_rho.png')
'''

###############################################################################
# Pressure against height
'''
figure()
plot(h,10**logP,lw=2)
xlabel('height [km]',size=fs)
ylabel(r'pressure [dyn/cm$^2$]',size=fs)
title('Earth atmosphere',size=fs)
grid('on')
yscale('log')
savefig(path+'earth_h_P.png')
'''

################################################################################
# Pressure and density stratifications together in normalized units
'''
P = 10**(logP)
N = 10**(logN)

P /= max(P)
N /= max(N)

figure()
plot(h,P,lw=2,label='pressure')
plot(h,N,lw=2,label='density')
xlabel('height [km]',size=fs)
ylabel('pressure and density in normalized units',size=fs)
title('Earth atmosphere',size=fs)
grid('on')
savefig(path+'earth_h_PN.png')
'''

################################################################################
# Mean molecular weight against height
'''
mu = 10**(logrho)/(10**(logN)*ssb.m_H)
figure()
plot(h,mu,lw=2)
xlabel('height [km]',size=fs)
ylabel('mean molecular weight',size=fs)
title('Earth atmosphere',size=fs)
grid('on')
savefig(path+'earth_h_mu.png')
'''

################################################################################
# Estimate the density scale height of the lower terrestrial atmosphere
'''
rho = 10**logrho
rho1 = rho[0]/e
index = argmin(abs(rho-rho1))
height = h[index]-h[0]
print 'Density scale height = %.2f km'%height.value

# Mount Everest
h_ME = 8.848*u.km
rho_func = lambda h: rho[0]*exp(-h/height)
everest = argmin(abs(h-h_ME))
#print 'Everest to ground ratio = %.3f' %(rho[everest]/rho[0])
print 'Everest to ground ratio = %.3f' %(rho_func(h_ME)/rho[0]) # using function
print 'It is %.3f times harder to breath at Mount Everest than at the surface level\n'%(rho[0]/rho_func(h_ME)-1)
'''

################################################################################
# Compare the terrestrial parameter values to the solar ones, at the base
# of each atmosphere. Ratio of the particle densities at h=0 in both.
'''
h_sun,tau5,colm_sun,temp_sun,vturb,nhyd,nprot,nel,ptot,pgasptot,dens = ssb.read_falc()
nhel = 0.1*nhyd
n_sun = nhyd+nel+nhel
print 'Ratio N_earth/N_sun = %.3f'%(10**logN[0]/n_sun[h_sun==0]).value
'''

################################################################################
# Estimate atmospheric column mass at the Earth's surface and compare to the Sun
'''
colm = (10**logP*u.dyn*u.cm**(-2)/ssb.g_E).decompose().cgs
print '\nAtmospheric column mass at Earth surface = %.3f g cm^-2'%colm[h==0].value
print 'Column mass at base of stellar photosphere= %.3f g cm^-2'%colm_sun[h_sun==0].value
'''

################################################################################
# Irradiance
'''
temp_eff = 5770*u.K
n_photon   = (20*temp_eff**3/(2*pi))#.value
D = ap.constants.au.cgs
R = ap.constants.R_sun.cgs
N_p_E = pi*R**2/D**2*n_photon

print '\nSunshine proton density = %.3e'%N_p_E.value
print 'Particle density in air = %.3e'%10**logN[h==0]
print 'Local thermal photon prod = %.3e'%(20*temp[h==0]**3).value
'''

show()


