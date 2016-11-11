; file: planck.pro = Planck function, now wavelengths in cm
; last: Nov 19 1999

function planck,temp,wav
; computes Planck function in erg / cm2 sec [delta lambda=1 micron] ster
; input: temp = temperature (K)

; physics constants in cgs (all cm)
k=1.380658D-16        ; Boltzmann constant (erg K; double precision)
h=6.626076D-27        ; Planck constant (erg s)
c=2.997929D10          ; velocity of light (cm/s)

blambda = 2*h*c^2/(wav^5*(exp(h*c/(wav*k*temp))-1))

return,blambda
end
