; file: exthmin.pro = H-minus bound-free + free-free extinction per H atom
; last: June 4 1999
; note: taken from Gray 1992 p135 ff
;       plots Gray p140-141 are in 10^-26 cm^2/etc

function exthmin,wav,temp,eldens
  ; in:  wav = wavelength [Angstrom] (float or fltarr)
  ;      temp = temperature [K]
  ;      eldens = electron density [electrons cm-3]
  ; out: H-minus bf+ff extinction [cm^2 per neutral hydrogen atom]
  ;      assuming LTE ionization H/H-min

; physics constants in cgs (all cm)
k=1.380658D-16   ; Boltzmann constant (erg/K; double precision)
h=6.626076D-27   ; Planck constant (erg s)
c=2.997929D10    ; velocity of light (cm/s)

; other parameters
wav=float(wav)
theta=5040./temp
elpress=eldens*k*temp

; evaluate H-min bound-free per H-min ion = Gray (8.11)
; his alpha = my sigma in NGSB/AFYC (per particle without stimulated)
sigmabf = 1.99654 -1.18267E-5*wav +2.64243E-6*wav^2 -4.40524E-10*wav^3 $
          +3.23992E-14*wav^4 -1.39568E-18*wav^5 +2.78701E-23*wav^6
sigmabf=sigmabf*1E-18  ; cm^2 per H-min ion
if total((wav ge 16444)) then sigmabf(where(wav ge 16444))=0.  
    ; H-min ionization limit 
    ; Thijs Krijger trick to permit array input

; convert into bound-free per neutral H atom assuming Saha = Gray p135
; units: cm2 per neutral H atom in whatever level (whole stage)
graysaha=4.158E-10*elpress*theta^2.5*10.^(0.754*theta)    ; Gray (8.12)
kappabf=sigmabf*graysaha                         ; per neutral H atom
kappabf=kappabf*(1.-exp(-h*c/(wav*1E-8*k*temp))) ; correct stimulated

; check Gray's Saha-Boltzmann with AFYC (edition 1999) p168
; logratio=-0.1761-alog10(elpress)+alog10(2.)+2.5*alog10(temp)-theta*0.754
; print,'Hmin/H ratio=',1/(10.^logratio) ; OK, same as Gray factor SB

; evaluate H-min free-free including stimulated emission = Gray p136
lwav=alog10(wav)
f0 =  -2.2763 -1.6850*lwav +0.76661*lwav^2 -0.0533464*lwav^3
f1 =  15.2827 -9.2846*lwav +1.99381*lwav^2 -0.142631*lwav^3
f2 = -197.789 +190.266*lwav -67.9775*lwav^2 +10.6913*lwav^3 -0.625151*lwav^4
ltheta=alog10(theta)
kappaff = 1E-26*elpress*10^(f0+f1*ltheta+f2*ltheta^2)   ; Gray (8.13)

return,kappabf+kappaff
end

