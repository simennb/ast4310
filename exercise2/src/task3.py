from ssb import *

# some plotting parameters
path = '../figures/task3/'
fs = 16  # font size

ssb = SSB()

h,tau5,colm,temp,vturb,nhyd,nprot,nel,ptot,pgasptot,dens = ssb.read_falc()
nhel = 0.1*nhyd 
wav,F, F_c, I, I_c = ssb.read_solspect()
wave = wav.value

sigma, wav_vac, lines = ssb.read_NaID()

################################################################################
# Plot solar NaID lines against vacuum wavelength at various dispersions
wav_vac = wav_vac.to('AA')
'''
figure()
plot(wav_vac,lines)#,lw=2)
xlabel(r'vacuum wavelength [$\AA$]',size=fs)
#ylabel(r'wavenumber $\sigma$ [cm$^{-1}$]',size=fs)
title('Solar NaID lines against vacuum wavelength',size=fs)
grid('on')
xlim([min(wav_vac.value),max(wav_vac.value)])
ylim([0.0,1.07])
savefig(path+'NaID_wav_vac_sigma.png')
'''

# Find the minimas for the NaID lines
ind = argmin(abs(5895-wav_vac.value))
min1 = wav_vac[argmin(lines[0:ind])]
min2 = wav_vac[ind-1+argmin(lines[ind-1:])]
print min2, min1
print ssb.vactoair(min2.value), ssb.vactoair(min1.value)

# Plot NaID lines against air wavelength
wav_air = ssb.vactoair(wav_vac.value)
wav_air *= u.AA

'''
figure()
plot(wav_air,lines)#,lw=2)
xlabel(r'air wavelength [$\AA$]',size=fs)
#ylabel(r'wavenumber $\sigma$ [cm$^{-1}$]',size=fs)
title('Solar NaID lines against air wavelength',size=fs)
grid('on')
xlim([min(wav_air.value),max(wav_air.value)])
ylim([0.0,1.07])
savefig(path+'NaID_wav_air_sigma.png')

show()
'''

################################################################################
# Check Saha-Boltzmann functions    WORKS
'''
saha = zeros(len(temp))
boltz = zeros(len(temp))
figure()
for s in range(1,4):
    for i in range(len(temp)):
        boltz[i] = ssb.boltz_Na(temp[i].value, s)
    plot(h,boltz,lw=2)
title('Boltzmann distribution of Na I in FALC',size=fs)
xlabel('height [km]',size=fs)
ylabel('population fraction n_1,s/N_1',size=fs)
xlim([-100,2000])
grid('on')
savefig(path+'boltzmann.png')
figure()
for r in range(1,3):
    for i in range(len(temp)):
        saha[i]  = ssb.saha_Na(temp[i].value, nel[i].value,r)
    plot(h,saha,lw=2)
title('Saha distribution of Na in FALC',size=fs)
xlabel('height [km]',size=fs)
ylabel('ionization fraction N_r/N_total',size=fs)
xlim([-100,2000])
yscale('log')
ylim([1e-4,1e1])
grid('on')
savefig(path+'saha.png')
'''

################################################################################
# NaI D1 line D:

lambda0 = ssb.vactoair(min1.value)#5895.94
wav = linspace(lambda0-2,lambda0+2,250)
I = zeros(len(wav))

h = h.value
temp = temp.value
eldens = nel.value
nhyd = nhyd.value
vmicro = (vturb).value
pgas = (pgasptot*ptot).value

for i in range(len(wav)):
    I[i] = ssb.emergent_intensity_Na(h,tau5,temp,nhyd,nprot.value,eldens,wav[i]*1e-4,vmicro,pgas)

figure()
plot(wav,I/max(I),lw=2)
plot(wav_air,lines/max(lines),'--',lw=2)
xlim([wav[0],wav[-1]])
xlabel(r'wavelength [$\AA$]',size=fs)
ylabel('normalized intensity',size=fs)
title('Na I D1 in LTE in FALC',size=fs)
grid('on')
savefig(path+'Na_D1.png')

show()
