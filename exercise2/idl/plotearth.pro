; file: plotearth.pro = IDL main to plot temp-height Earth atmosphere
; last: Aug  2 2010
; note: from Allen but here are much better tables

; read earth.dat
close,2
openr,2,'earth.dat'
nhE=25
earth=fltarr(5,nhE)
readf,2,earth
hE=earth[0,*]
pgasE=10.^earth[1,*]
tempE=earth[2,*]
densE=10.^earth[3,*]
npartE=10.^earth[4,*]

;@psplot_params.idl   ; RJR x/ps switching params
;@psplot_begin.idl    ; RJR x/ps loop start

; plot Earth temperature stratification
filename='../figs/earth_temp_height.ps'
plot,hE,tempE,$
  xtitle='height [km]',ytitle='temperature [K]
xyouts,0.3,0.8,'EARTH ATMOSPHERE',/norm

;@psplot_end.idl       ; RJR x/ps loop end
end




