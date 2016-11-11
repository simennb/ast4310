; file: brighttemp.pro = brightness temperature = Planck inversion
; last: May 13 1999

function brighttemp,wav,intensity
; in: wav = wavelength in micron
;     intensity = ergs / cm2 s ster [micron bandwith]
; out: brightness temperature

; physics constants in cgs (all cm)
k=1.380658D-16        ; Boltzmann constant (erg K; double precision)
h=6.626076D-27        ; Planck constant (erg s)
c=2.997929D10         ; velocity of light (cm/s)

wavcm=wav*1E-4        ; change wav into cm
intcm=intensity*1E4   ; change into per cm bandwidth
tempbright=h*c/(wavcm*k)/alog(2*h*c^2/(intcm*wavcm^5)+1)

return,tempbright
end
