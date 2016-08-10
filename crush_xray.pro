pro crush_xray


tic
;testing commentsOB

for i=0,1 do begin

FMT='X,F,F,F,F,'

readcol,'XrayCOSMOS_flux.dat', F=FMT, ra, dec, flux1, flux2


;;;;;Cutting past the maglim
if i EQ 13 then begin
flux=flux1
endif else begin
flux=flux2
endelse


img=readfits("template.fits",head)
adxy, head, ra, dec, x0, y0

x0=round(x0)
y0=round(y0)
;;ymax(i,k)=max(y0)

ns= n_elements(flux)
;get the degrees per pixel count from header
cdelt=sxpar(head,'CDELT1')
cdeltsquare= abs(cdelt)^2
if i EQ 0 then lambda=9.91e-10

if i EQ 1 then lambda=2.06e-10

;lambda in meters
;begin loop for adding catalog data into template file by coordinates
for h=0,ns-1 do begin
;check to see if there are flux values less than -99  because these
;mean no data. if there are, add in 0 into the template coordinate
  if (flux[h] LE -99) then begin
   img[x0[h],y0[h]] += flux[h]*0l
endif else begin
   if x0[h] LE 6322 && x0[h] GE 0 && y0[h] LE 6195 && y0[h] GE 0 then begin
         img[x0[h],y0[h]] += flux[h]*1e-23*1e6*1d6*10d-23*3d8/lambda/(cdeltsquare/3282.8)
     endif
  ;img[x0[h],y0[h]] += flux[h]*1d6*1d9*(cdeltsquare)/3282.8
  endelse
endfor


;convert the sigma of 2arcseconds into pixels using cdelt from header
sigma=3/(3600*abs(cdelt))
;apply gaussian smoothing with a sigma of 2 arcseconds
  img=GAUSS_SMOOTH(img, sigma)
;write out new fits image with name redhisftimage+the bottomrange of
;the redshiftrange+.fits
  writefits,  'xray_'+strtrim(i,1)+'.fits', img, head;
;;write_csv, 'ymax.dat', ymax
endfor

end
