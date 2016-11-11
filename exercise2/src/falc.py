from ssb import *
# np, ap, u and pylab *

# some plotting parameters
path = '../figures/falc/' # folder where figures are stored
fs = 16  # font size

ssb = SSB()

h,tau5,colm,temp,vturb,nhyd,nprot,nel,ptot,pgasptot,dens = ssb.read_falc()

nhel = 0.1*nhyd 

'''
################################################################################
# Temperature stratification
figure()
plot(h, temp,lw=2)
grid('on')
xlabel('height [km]',size=fs)
ylabel('temperature [K]',size=fs)
title('Temperature stratification of the FALC model',size=fs)
ylim([0,20000])
savefig(path+'falc_h_temp.png')
'''

################################################################################
# Plot total pressure against column mass, linearly and logarithmically
'''
figure()
plot(colm,ptot,lw=2)
xlabel(r'column mass [g cm$^{-2}$]',size=fs)
ylabel(r'total pressure [dyn cm$^{-2}$]',size=fs)
grid('on')
savefig(path+'falc_colm_ptot.png')
yscale('log')
savefig(path+'falc_colm_ptot_ylog.png')
'''

################################################################################
# Ratio of hydrogen mass density versus height
'''
figure()
plot(h,nhyd*ssb.m_H/dens,label='Hydrogen',lw=2)
plot(h,nhel*ssb.m_He/dens,label='Helium',lw=2)
#plot(h,(nhel*ssb.m_He+nhyd*ssb.m_H)/dens)
xlabel('height [km]',size=fs)
ylabel('ratio of mass density',size=fs)
title('Ratio of mass density vs total mass density',size=fs)
legend(loc='best',fontsize=14)
grid('on')
savefig(path+'falc_h_dens.png')

# Fraction of metal
print 'Fraction of metal in the model mix = %.3e'% mean(1-(nhyd*ssb.m_H+nhel*ssb.m_He)/dens)
'''

################################################################################
# Plot column mass against height
'''
figure()
plot(h,colm,lw=2)
xlabel('height [km]',size=fs)
ylabel('column mass [g cm$^{-2}$]',size=fs)
title('Column mass against height',size=fs)
grid('on')
savefig(path+'falc_h_colm_lin.png')
yscale('log')
savefig(path+'falc_h_colm_ylog.png')
'''

################################################################################
# Density scale height
'''
height = 2017*u.km
rho0 = dens[h==0]
rhoN = dens[h==height]
H_p = height/(log(rho0/rhoN))
print 'Density scale height = %.2f km\n'%H_p.value

# Plot gas density against height
figure()
plot(h,dens,lw=2,label='FALC')
plot(h,rho0*exp(-h/H_p),'--',lw=2,label=r'$H_{\rho}$ = %.1f km'%H_p.value)
xlabel('height [km]',size=fs)
ylabel('gas density [g cm$^{-3}$]',size=fs)
title('Gas density against height',size=fs)
grid('on')
legend(loc='best',fontsize=14)
savefig(path+'falc_h_dens2.png')
yscale('log')
savefig(path+'falc_h_dens3.png')
'''

################################################################################
# Gas pressure plotted against height
'''
for i in range(3): # because multiple plots.........
    pgas = pgasptot*ptot
    p_str = 'n_H + n_e'
    if i==0:
        prod = (nhyd + nel)*ssb.k_B*temp
        yl = r'gas pressure [dyn cm$^{-2}$]'
    elif i==1:
        prod = (nhyd + nel)*ssb.k_B*temp
        prod /= ptot
        pgas /= ptot
        yl = r'relative gas pressure'
    elif i==2:
        prod = (nhyd + nel + nhel)*ssb.k_B*temp
        prod /= ptot
        pgas /= ptot
        p_str = 'n_H + n_e + n_{He}'
        yl = r'relative gas pressure'
    figure()
    plot(h,pgas,lw=2,label='Gas pressure')
    plot(h,prod,lw=2,label=r'$(%s)kT$'%p_str)
    xlabel(r'height [km]',size=fs)
    ylabel(yl,size=fs)
    yscale('log')
    xscale('log')
    if i!=0:
        ylim([0.5,1.5])
    title('Plot of gas pressure against height',size=fs)
    legend(loc='best',fontsize=14)
    grid('on')
    savefig(path+'falc_h_pgas_%d.png'%i)
'''

################################################################################
# Plot of total hydrogen density against height
'''
figure()
plot(h,nhyd,lw=2,label='hydrogen')
plot(h,nel,lw=2,label='all electrons')
plot(h,nprot,lw=2,label='protons')
plot(h,(nel-nprot),lw=2,label='electrons not from H')

xlabel('height [km]',size=fs)
ylabel(r'density [cm$^{-3}$]',size=fs)
legend(loc='best',fontsize=14)
grid('on')
savefig(path+'falc_h_n.png')
yscale('log')
savefig(path+'falc_h_n_log.png')
'''

################################################################################
# Plot of the hydrogen ionization fraction, log
'''
figure()
plot(h,nprot/nhyd,lw=2)
xlabel('height [km]',size=fs)
ylabel('fraction of ionized hydrogen',size=fs)
title('Plot of fraction of ionized hydrogen',size=fs)
yscale('log')
grid('on')
savefig(path+'falc_h_Hion.png')
'''

################################################################################
# do final part of 1.2, it makes little sense
# at deepest part of the atmosphere

print 'At h = %.3f :'%(h[-1]).value
print ' n_hydrogen = %.3e'%(nhyd[-1]).value
print ' n_photon   = %.3e'%(20*temp[-1]**3).value

# at highest part
print '\nAt h = %.3f :'%(h[0]).value
print ' n_hydrogen = %.3e'%(nhyd[0]).value
temp_eff = 5770*u.K
print ' n_photon   = %.3e'%(20*temp_eff**3/(2*pi)).value


show()
