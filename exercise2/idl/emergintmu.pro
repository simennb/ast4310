
; emergent intensity at wavelength wl (micron) and angle mu
ext=fltarr(nh)
tau=fltarr(nh)
integrand=fltarr(nh)
contfunc=fltarr(nh)
int=0.
hint=0.
for ih=1,nh-1 do begin 
  ext[ih]=exthmin(wl*1E4,temp[ih],nel[ih])*(nhyd[ih]-nprot[ih])$
          +0.664E-24*nel[ih]
  tau[ih]=tau[ih-1]+0.5*(ext[ih]+ext[ih-1])*(h[ih-1]-h[ih])*1E5
  integrand[ih]=planck(temp[ih],wl)*exp(-tau[ih]/mu)
  int=int+0.5*(integrand[ih]+integrand[ih-1])*(tau[ih]-tau[ih-1])/mu
  hint=hint+h[ih]*0.5*(integrand[ih]+integrand[ih-1])*(tau[ih]-tau[ih-1])/mu
  contfunc[ih]=integrand[ih]*ext[ih]/mu
endfor
hmean=hint/int



