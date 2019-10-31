function doppler,pic,m,ber=ber,plot=plot

;**************************************************************************
;
; NAME:
;	DOPPLER
; PURPOSE:
;	Compute the dopplershift(s) of a lineprofile, resp. spectrum
;	opional interpolation of the linecore intensity
; CATEGORY:
;	data analysis
; CALLING SEQUENCE:
;	result = DOPPLER(data,[I_c],[ber=ber])
; INPUTS:
;	data     : 1 or 2-dim array with profile/spectrum
; OPTIONAL INPUT PARAMETER:
;	I_c      : returns the interpolated intensity in the linecore
; KEYWORDS:
;       ber      : number of points used for polynomial fit
;		   default: 5
; OUTPUTS:
; 	scalar or 1-dim floatarray with the abscissavalues of the 
;	linecore(s)
; RESTRICTIONS:
; 	parabolic fit might be bad for strongly asymmetric lines
; PROCEDURE:
; 	find minimum m of dataset
;	2nd order polyfit using the range [m-ber,m+ber] around the
;	minimum.
;	return minimum of fit as value for dopplershift, optional 
;	fit(minimum) as central intensity
; MODIFICATION HISTORY:
;	8-Aug-1991 PS
;
;*************************************************************************


if not keyword_set(ber) then ber=5

s=size(pic)
if (s(0) eq 1) then begin
  resx=s(1)
  resy=1
  p=fltarr(resx,2)
  p(*,0)=pic
endif else begin
  resx=s(1)
  resy=s(2)
  p=pic
endelse

dv=fltarr(resy)
m=fltarr(resy)

for i=0,resy-1 do begin
  zeile=p(*,i)
  mitte = where(zeile eq min(zeile))
  mitte=fix(mitte(0))
  xu=max([mitte-ber,0])
  xo=min([mitte+ber,resx])
  zeile=zeile(xu:xo)
  x=indgen(xo-xu+1)
  c=poly_fit(indgen(xo-xu+1),zeile,2,yfit=fit)
  if c(2) ne 0 then zent=-1*c(1)/(2*c(2))
  if zent ge 0 and zent le (xo-xu+1) then $
    dv(i)=xu+zent $
  else $
    dv(i)=xu+mitte
  m(i)=c(0)+c(1)*zent+c(2)*zent^2
end

if n_elements(dv) eq 1 then begin
  dv=dv(0)
  m=m(0)
endif

if keyword_set(plot) then begin

plot,x,p(xu:xo)
oplot,x,fit,psym=-4,lines=1

endif


return ,dv
end
