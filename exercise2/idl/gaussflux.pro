; file: gaussflux.pro = last part of file ssb2.pro
; last: Aug  2 2010

; ===== three-point Gaussian integration intensity -> flux
; abscissae + weights n=3 Abramowitz & Stegun page 916
xgauss=[-0.7745966692,0.0000000000,0.7745966692]
wgauss=[ 0.5555555555,0.8888888888,0.5555555555]
fluxspec=fltarr(nwav)
intmu=fltarr(3,nwav)
for imu=0,2 do begin
  mu=0.5+xgauss[imu]/2.   ; rescale xrange [-1,+1] to [0,1]
  wg=wgauss[imu]/2.       ; weights add up to 2 on [-1,+1]
  for iw=0,nwav-1 do begin
    wl=wav[iw]
    @emergintmu.pro       ; old trapezoidal integration I(0,mu)
    intmu[imu,iw]=int
    fluxspec[iw]=fluxspec[iw]+wg*intmu[imu,iw]*mu
  endfor
endfor
fluxspec=2*fluxspec    ; no !pi, AQ has flux F, not {\cal F}
plot,wav,fluxspec,$
  xrange=[0,2],yrange=[0,5E10],$
  xtitle='wavelength [micron]',ytitle='solar flux'
oplot,wav,Fcont,linestyle=2
xyouts,0.5,4E10,'computed'
xyouts,0.35,1E10,'observed'
xyouts,1.2,3E10,'GAUSSIAN'

