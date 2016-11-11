; file: readfalc.pro = IDL main to read & plot falc.dat = FALC model
; last: Aug  2 2010 
; note: file FALC.DAT from email Han Uitenbroek Apr 29 1999 
;       = resampled Table 2 of Fontenla et al ApJ 406 319 1993
;       = 2 LaTeX tables in falc.tex

close,1
openr,1,'falc.dat'
nh=80                ; nr FALC height values
falc=fltarr(11,nh)   ; 11 FALC columns
readf,1,falc
  h=falc[0,*]
  tau5=falc[1,*]
  colm=falc[2,*]
  temp=falc[3,*]
  vturb=falc[4,*]
  nhyd=falc[5,*]
  nprot=falc[6,*]
  nel=falc[7,*]
  ptot=falc[8,*]
  pgasptot=falc[9,*]
  dens=falc[10,*]

; plot falc model
plot,h,temp,yrange=[3000,10000],$
  xtitle='height [km]',ytitle='temperature [K]'

; repeat plot on postscript file
set_plot,'ps'      ; repeat plot on PostScript file
device,filename='falc_temp_height.ps'
plot,h,temp,yrange=[3000,10000],$
  xtitle='height [km]',ytitle='temperature [K]'
device,/close
set_plot,'x'       ; return to screen; set_plot,'win' for Windows
end









