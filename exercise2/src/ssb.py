from pylab import *
from matplotlib import rc
rc('font',**{'family':'serif'})
import numpy as np
import astropy as ap
from astropy import units as u
from scipy import special

# Useful constants and other things
class SSB:
    def __init__(self):
        self.m_H = ap.constants.m_p.cgs   # mass of hydrogen/proton [g]
        self.m_He = 3.97*self.m_H         # mass of helium [g]
        self.m_e = ap.constants.m_e.cgs   # mass of electron [g]
        self.h = ap.constants.h.cgs       # planck's constant [erg s]
        self.k_B = ap.constants.k_B.cgs   # Boltzmann's constant [erg/K]
        self.k_eV = self.k_B.to('eV/K')   # Boltzmann's constant [eV/K]
        self.c = ap.constants.c.cgs       # speed of light [cm/s]
        self.g_E = ap.constants.g0.cgs    # gravitation acc. on earth [cm/s^2]
        self.sigma_T = 6.648e-25*u.cm**2  # Thompson e scattering cross section 
        self.hc = (self.h*self.c).to('eV AA') # in order to get photon energy

        # Initializing some things for Na I lines
        self.g_rs = [2.,2.,4.] # ground state, 3pS1/2, 3pS3/2 level
        self.chi_ion = [5.139, 47.29, 71.64] # ionization energies in eV
        self.chi_level = zeros(3) # chi[0] = 0 cause shared ground state
        self.chi_level[1] = self.hc.value/5895.94   # energy of NaI D1 [eV]
        self.chi_level[2] = self.hc.value/5889.97   # energy of NaI D2 [eV]
        self.wav0 = 5895.94

    ############################################################################
    ############################################################################
    # Functions for reading from file

    def read_falc(self):
        h,tau5,colm,temp,vturb,nhyd,nprot,nel,ptot,pgasptot,dens = np.loadtxt('../data/falc.dat',unpack=True)
        # give each an astropy unit
        h *= u.km
        colm  *= u.g*u.cm**(-2)
        temp  *= u.K
        vturb *= u.km/u.s
        nhyd  *= u.cm**(-3)
        nprot *= u.cm**(-3)
        nel   *= u.cm**(-3)
        ptot  *= u.dyn*u.cm**(-2)
        dens  *= u.g*u.cm**(-3)
        
        return h, tau5, colm, temp, vturb, nhyd, nprot, nel, ptot, pgasptot, dens

    def read_earth(self):
        h,logP,temp,logrho,logN = np.loadtxt('../data/earth.dat',unpack=True)
        # giving each an astropy unit
        h *= u.km
        #logP *= u.dyn*u.cm**(-2)
        temp *= u.K
        #logrho *= u.g*u.cm**(-3)
        #logN *= u.cm**(-3)

        return h, logP, temp, logrho, logN

    def read_solspect(self):
        wav, F, F_c, I, I_c = np.loadtxt('../data/solspect.dat',unpack=True)
        wav *= u.micron
        F   *= u.erg/(u.cm**2*u.s)
        F_c *= u.erg/(u.cm**2*u.s) #1 # 10**10 erg cm**-2 s**-1
        I   *= u.erg/(u.cm**2*u.s*u.micron) # neglecting steradian
        I_c *= u.erg/(u.cm**2*u.s*u.micron)

        # for some reason not allowed to do in same operation as units
        # Scaling back from 1e10 erg to erg
        F   *= 1e10
        F_c *= 1e10
        I   *= 1e10
        I_c *= 1e10

        return wav, F, F_c, I, I_c
    
    def read_NaID(self):
        sigma,lines= np.loadtxt('../data/int_nad.dat',usecols=(0,2),unpack=True)

        sigma /= u.cm
        wav_vac = 1./sigma
        
        return sigma, wav_vac, lines

    ############################################################################
    ############################################################################
    # Functions for calculating various things

    def planck(self, temp, wav):
        # making sure wavelength is in cgs
        wav = wav.to('cm')
        h = self.h
        k_B = self.k_B
        c = self.c

        return 2*h*c**2/wav**5*1./(exp(h*c/(wav*k_B*temp))-1)
    
    def brightness_temp(self, I_c, wav):
        # making sure wavelength is in cgs
        wav = wav.to('cm')
        h = self.h
        k_B = self.k_B
        c = self.c
        
        return h*c/(wav*k_B)/log(2*h*c**2/(I_c*wav**5)+1)

    def exthmin(self, wav, temp, eldens):
        '''
        H-minus extinction, from Gray 1992
        input:
         wav = wavelength [Angstrom] (float or float array)
         temp = temperature [K]
         eldens = electron density [electrons cm-3]
        output:
         H-minus bf+ff extinction [cm^2 per neutral hydrogen atom]
         assuming LTE ionization H/H-min
        '''
        wave = (wav.to('AA')).value # for simplicity
        # other parameters
        theta = 5040./temp.value
        elpress = eldens*self.k_B*temp

        # evaluate H-min bound-free per H-min ion = Gray (8.11)
        # his alpha = my sigma in NGSB/AFYC (per particle without stimulated)
        sigmabf = (1.99654 -1.18267E-5*wave +2.64243E-6*wave**2
                   -4.40524E-10*wave**3 +3.23992E-14*wave**4
                   -1.39568E-18*wave**5 +2.78701E-23*wave**6)
        sigmabf *= 1e-18  # cm^2 per H-min ion
        if size(wave) > 1:
            sigmabf[wave > 16444] = 0 # H-min ionization limit at lambda=1.6444 micron
        elif (size(wave) == 1):
            if wave > 16444:
                sigmabf = 0

        # Astropy unit, adding because simplicity
        sigmabf *= u.cm**2

        # convert into bound-free per neutral H atom assuming Saha = Gray p135
        # units: cm2 per neutral H atom in whatever level (whole stage)
        graysaha = 4.158e-10*elpress*theta**2.5*10.**(0.754*theta) # Gray (8.12)

        kappabf = sigmabf*graysaha     # per neutral H atom
        kappabf = kappabf*(1.-np.exp(-self.h*self.c/(wav*1e-8*self.k_B*temp))) # correct stimulated
        kappabf = kappabf.value/u.cm # shhh...

        # check Gray's Saha-Boltzmann with AFYC (edition 1999) p168
        #logratio=-0.1761-np.log10(elpress.value)+np.log10(2.)+2.5*np.log10(temp.value)-theta*0.754
        #print 'Hmin/H ratio=',1/(10.**logratio) # OK, same as Gray factor SB
        # evaluate H-min free-free including stimulated emission = Gray p136
        lwav=np.log10(wave)
        f0 = -2.2763 -1.6850*lwav +0.76661*lwav**2 -0.0533464*lwav**3
        f1 = 15.2827 -9.2846*lwav +1.99381*lwav**2 -0.142631*lwav**3
        f2 = (-197.789 +190.266*lwav -67.9775*lwav**2 +10.6913*lwav**3
               -0.625151*lwav**4)
        ltheta=np.log10(theta)
        kappaff = 1E-26*elpress*10**(f0+f1*ltheta+f2*ltheta**2) # Gray (8.13)
        kappaff = kappaff.value/u.cm

        return kappabf+kappaff

    def vactoair(self,wav_vac):
        return 0.99972683*wav_vac + 0.0107-196.25/wav_vac

    ############################################################################
    ############################################################################
    # Integrations
    
    def optical_depth(self, h, temp, eldens,n_neutral, wav=500*u.nm):
        '''
        calculates the optical depth integral
        '''
        tau = np.zeros(len(h), dtype=float) # initializing tau array
        ext = self.exthmin(wav, temp, eldens)
        # add thompson and stuff
        ext = ext*n_neutral.value + (self.sigma_T*eldens)

        for i in range(1,len(tau)):
            tau[i] = tau[i-1] + 0.5*(ext[i]+ext[i-1])*(h[i-1]-h[i])
        # index zero is not accounted for, so tau[0] = 0 because we have already initialized

        return tau

    def emergent_intensity(self,h,tau5,temp,nhyd,nprot,nel, wl=0.5*u.micron,mu=1.0):
        ext = np.zeros(len(tau5))/u.cm
        tau = np.zeros(len(tau5))
        integrand = np.zeros(len(tau5))*u.erg/(u.cm**3*u.s)
        contfunc = np.zeros(len(tau5))*u.erg/(u.cm**4*u.s)
        intt = 0.0*integrand.unit
        hint = 0.0*h.unit*integrand.unit
        
        for i in range(1, len(tau5)):
            ext[i] = (self.exthmin(wl, temp[i], nel[i])*(nhyd[i]-nprot[i]).value
                      + self.sigma_T*nel[i])
            tau[i] = tau[i-1] + 0.5 * (ext[i] + ext[i-1]) * (h[i-1]-h[i])#*1E5
            integrand[i] = (self.planck(temp[i],wl)*np.exp(-tau[i]/mu))
            intt += 0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])/mu
            hint += h[i]*0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])#/mu
            contfunc[i] = integrand[i]*ext[i]#.value
        # note : exthmin has wavelength in [Angstrom], planck in [cm]
        hmean = hint / intt

        return intt, contfunc, hmean

    def flux_integration(self,h,tau5,temp,nhyd,nprot,nel, wav):
        #remove units
        h = h.value; temp = temp.value; nhyd = nhyd.value
        nprot = nprot.value; nel = nel.value; wav = wav.value

        # SSB 2.7 page 17: flux integration
        # ===== three-point Gaussian integration intensity -> flux
        # abscissae + weights n=3 Abramowitz & Stegun page 916
        xgauss=[-0.7745966692,0.0000000000,0.7745966692]
        wgauss=[ 0.5555555555,0.8888888888,0.5555555555]
        fluxspec = np.zeros(len(wav),dtype=float)
        intmu = np.zeros((3,len(wav)), dtype=float)
        for imu in range(3):
            mu=0.5+xgauss[imu]/2. # rescale xrange [-1,+1] to [0,1]
            wg=wgauss[imu]/2. # weights add up to 2 on [-1,+1]
            for iw in range(0,len(wav)):
                wl=wav[iw]
                ext = np.zeros(len(tau5))
                tau = np.zeros(len(tau5))
                integrand = np.zeros(len(tau5))
                intt = 0.0
                for i in range(1, len(tau5)):
                    ext[i] = (self.exthmin(wl*u.micron, temp[i]*u.K, nel[i]/u.cm**(-3)).value*(nhyd[i]-nprot[i])
                              + self.sigma_T.value*nel[i])
                    tau[i] = (tau[i-1] + 0.5 * (ext[i] + ext[i-1]) *
                              (h[i-1]-h[i])*1E5)

                    integrand[i] = self.planck(temp[i]*u.K,wl*u.micron).value*np.exp(-tau[i]/mu)
                    intt += 0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])/mu

                intmu[imu,iw]=intt
                fluxspec[iw]=fluxspec[iw] + wg*intmu[imu,iw]*mu
        fluxspec *= 2 # no np.pi, Allen 1978 has flux F, not {\cal F}

        return fluxspec

    ############################################################################
    ############################################################################
    # NaID model
    
    def partfunc_Na(self,temp):
        theta = 5040./temp
        c0 =  0.30955
        c1 = -0.17778
        c2 =  1.10594
        c3 = -2.42847
        c4 =  1.70721
        logtheta = log10(theta)
        logU = c0+c1*logtheta+c2*logtheta**2+c3*logtheta**3+c4*logtheta**4

        U = zeros(3)
        U[0] = 10**logU
        U[1] = 1.  # Allen 1976
        U[2] = 6.  # Allen 1976
        
        return U
    
    def saha_Na(self, temp, eldens,rstage):
        h = self.h.value
        k_erg = self.k_B.value
        k_eV  = self.k_eV.value
        m_e = self.m_e.value

        U = self.partfunc_Na(temp)
        U = np.append(U,2)
        sahaconst = 2./eldens*(2*pi*m_e*k_erg*temp/h**2)**1.5

        nstage = zeros(len(self.chi_ion)+1)
        nstage[0] = 1
        for r in range(len(self.chi_ion)):
            nstage[r+1] = nstage[r]*sahaconst*U[r+1]/U[r]*exp(-self.chi_ion[r]/(k_eV*temp))
            
        ntotal = sum(nstage)
        nstagerel = nstage/ntotal
        
        return nstagerel[rstage-1]
    
    def boltz_Na(self, temp, s):
        k_eV = self.k_eV.value
        U = self.partfunc_Na(temp)

        relnrs =  self.g_rs[s-1]/U[0]*exp(-self.chi_level[s-1]/(k_eV*temp))
        return relnrs

    def sahaboltz_Na(self, temp, eldens, r, s):
        return self.saha_Na(temp,eldens,r)*self.boltz_Na(temp,s)

    def dopplerwidth(self, wav, temp, vmicro, atmass):
        c = self.c.value
        k_B = self.k_B.value


        return wav/c*sqrt(2*k_B*temp/atmass + vmicro**2*1e10)

    def voigt(self, gamma, x):
        z = (x+1j*gamma)
        V = special.wofz(z).real

        return V

    def gammavdw(self, temp, pgas, s):
        # Van der Waals broadening for Na D1 and Na D2
        # s=2 : Na D1
        # s=3 : Na D2
        # using classical recipe by Unsold
        # following recipe in SSB
        rsq_u = self.rsq(s)
        rsq_1 = self.rsq(1) # lower level D1 and D2 lines is ground state s=1
        loggvdw = (6.33 + 0.4*log10(rsq_u - rsq_1)
                   + log10(pgas) - 0.7 * log10(temp))

        return 10**loggvdw

    def rsq(self, s):
        # compute mean square radius of level s of Na D1 and Na D2 transitions
        # -> needed for van der Waals broadening in SSB
        # s=1 : ground state, angular momentum l=0
        # s=2 : Na D1 upper level l=1
        # s=3 : Na D2 upper level l=1
        E_ionization = self.chi_ion[0]
        E_n = self.chi_level
        Z = 1.                # ionization stage, neutral Na: NaI
        Rydberg = 13.6        # Rydberg constant [eV]
        l = [0.,1.,1.]        # angular quantum number
        nstar_sq = Rydberg * Z**2 / (E_ionization - E_n[s-1])
        rsq=nstar_sq / 2. / Z**2 * (5*nstar_sq + 1 - 3*l[s-1]*(l[s-1] + 1))
        
        return rsq

    def line_broadening(self, wav, temp, pgas, vmicro, atmass, s,E=1.0):
        c = self.c.value
        dopplerwidth = self.dopplerwidth(wav, temp, vmicro, atmass)
        gammavdw = self.gammavdw(temp, pgas, s)
        a_voigt = wav**2 / (4*pi*c) * E*gammavdw / dopplerwidth
        v_voigt = (wav-self.wav0*1e-8)/dopplerwidth
        voigt_NaD = self.voigt(a_voigt,abs(v_voigt))
        voigt_NaD /= dopplerwidth
        return dopplerwidth, voigt_NaD

    def NaD1_ext(self, wav, temp, eldens, nhyd, vmicro, pgas):
        # Wavelength given in [cm]
        # not using astropy in order to increase speed :3 unsuccessful
        # Functions called are modified to simplify and only work for the
        # Na D1 line, as that is what we are interested in
        m_e = self.m_e.value
        c = self.c.value
        hc = (self.h*c).value
        k_B = self.k_B.value

        m_Na = 22.99*1.6605e-24  # natrium mass [g]
        A_Na = 1.8e-6   # sodium abundance (Allen 1976)
        f_D1 = 0.318    # oscillator strength for NaI D1
        f_D2 = 0.631    # oscillator strength for NaI D2
        e = 4.803204e-10  # electron charge in cgs [statcoloumb]
        # Natrium I D1 line
        r = 1
        s = 2
        bl = 1.0
        bu = 1.0

        # Function calls
        dopplerwidth, H = self.line_broadening(wav,temp,pgas,vmicro,m_Na,s,3)

        alpha1 = sqrt(pi)*e**2/(m_e*c)*wav**2/c*bl
        alpha2 = self.sahaboltz_Na(temp, eldens, r, 1)  # s=1 cause gs?!?!?!?
        alpha3 = nhyd*A_Na*f_D1
        alpha4 = H
        alpha5 = 1. - bu/bl*exp(-hc/(wav*k_B*temp))
        alpha = alpha1*alpha2*alpha3*alpha4*alpha5
        
        return alpha

    def emergent_intensity_Na(self,h,tau5,temp,nhyd,nprot,nel,wl,vmicro,pgas):
        # New almost identical function since i removed units in my Na
        # calculations as astropy seems to really slow things down

        # SSB 2.4 page 16 emergent intensity, contribution function and mean height of formation in FALC
        sigma_Thomson = 6.648E-25 # Thomson cross-section [cm^2]

        ext = np.zeros(len(tau5))
        tau = np.zeros(len(tau5))
        integrand = np.zeros(len(tau5))
        contfunc = np.zeros(len(tau5))
        intt = 0.0
        hint = 0.0
        for i in range(1, len(tau5)):
            ext[i] = (self.exthmin(wl*1e4*u.AA, temp[i]*u.K, nel[i]*u.cm**(-3)).value*(nhyd[i]-nprot[i])
                      + sigma_Thomson*nel[i] 
                      + self.NaD1_ext(wl*1e-4,temp[i-1],nel[i-1],nhyd[i-1],vmicro[i-1],pgas[i-1]))
            tau[i] = tau[i-1] + 0.5 * (ext[i] + ext[i-1]) * (h[i-1]-h[i])*1E5
            integrand[i] = self.planck(temp[i]*u.K,wl*u.micron).value*np.exp(-tau[i])
            intt += 0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])
            hint += h[i]*0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])
            contfunc[i] = integrand[i]*ext[i]
        # note : exthmin has wavelength in [Angstrom], planck in [cm]
        hmean = hint / intt

        return intt
