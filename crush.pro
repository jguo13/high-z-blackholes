pro crush
tic
;because the optical images are slightly smaller than the east field
;of the cosmos field, we have to crop them so we create empty arrays
;to store their max value, above which we cut
ymaxbe=Make_Array(10,10,10, value=0)
ymaxbw=Make_Array(10,10,10, value=0)
ymax=Make_Array(10,10,10, value=0)

for i=0,8 do begin
for j=0,1 do begin
      for k=0,4 do begin

FMT='F,F,X,F,F'

readcol, strtrim(i,1) + '_Mag_' + strtrim(j,1) + '_4.dat', F=FMT, raold, decold, fluxold , mag



;start the magnitude limit inputs

if i EQ 0 then begin
   if j EQ 0 then begin
      m0=27.0-6
   endif else begin
      m0=26.5-6
   endelse
endif

if i EQ 1 then begin
   if j EQ 0 then begin
      m0=26.9-6
   endif else begin
      m0=26.2-6
   endelse
endif

if i EQ 2 then begin
   if j EQ 0 then begin
      m0=26.4-6
   endif else begin
      m0=25.9-6
   endelse
endif

if i EQ 3 then begin
   if j EQ 0 then begin
      m0=25.8-6
   endif else begin
      m0=25.3-6
   endelse
endif

if i EQ 4 then begin
   if j EQ 0 then begin
      m0=24.9-6
   endif else begin
      m0=24.4-6
   endelse
endif

if i EQ 5 then begin
  if j EQ 0 then begin
      m0=25.4-6
   endif else begin
      m0=24.9-6
   endelse
endif

if i EQ 6 then begin
   if j EQ 0 then begin
      m0=25.0-6
   endif else begin
      m0=24.6-6
   endelse
endif

if i EQ 7 then begin
   if j EQ 0 then begin
      m0=25.2-6
   endif else begin
      m0=24.7-6
   endelse
endif

if i EQ 8 then begin
   if j EQ 0 then begin
      m0=25.5-6
   endif else begin
      m0=25.5-6
   endelse
endif


;start the increments
if k EQ 0 then begin
mlim=m0
endif else if k EQ 1 then begin
mlim=m0+2
endif else if k EQ 2 then begin
mlim=m0+4
endif else if k EQ 3 then begin
mlim=m0+6
endif else if k EQ 4 then begin
mlim=max(mag)
endif





magindex= WHERE(mag GE mlim, /NULL)
magarr=mag(magindex)

ra=raold(magindex)
dec=decold(magindex)
flux=fluxold(magindex)



img=readfits("template.fits",head)
adxy, head, ra, dec, x0, y0

x0=round(x0)
y0=round(y0)
ymax(i,j,k)=max(y0)

ns= n_elements(flux)
;get the degrees per pixel count from header
cdelt=sxpar(head,'CDELT1')
cdeltsquare= abs(cdelt)^2
;these are the bandwidths in meters
if i EQ 0 then lambda=6288.7e-10
if i EQ 1 then lambda=7683.9e-10
if i EQ 2 then lambda=9105.7e-10
if i EQ 3 then lambda=10214.2e-10
if i EQ 4 then lambda=9791.4e-10
if i EQ 5 then lambda=12534.6e-10
if i EQ 6 then lambda=16453.4e-10
if i EQ 7 then lambda=21539.9e-10
if i EQ 8 && j EQ 0 then lambda=35634.3e-10
if i EQ 8 && j EQ 1 then lambda=45110.1e-10
;begin loop for adding catalog data into template file by coordinates
for h=0,ns-1 do begin
;convert the value from microJansky to nw/m^2/sr
         img[x0[h],y0[h]] += flux[h]*10d-23*3d8/lambda/(cdeltsquare/3282.8)

endfor


;convert the sigma of arcseconds(determined by FHW from Irac) into
;pixels using cdelt from header 
  sigma=(1.6/2.3555)/(3600*abs(cdelt))
;apply gaussian smoothing 
  img=GAUSS_SMOOTH(img, sigma)
;write out new fits image with name redhisftimage+the bottomrange of
;the redshiftrange+.fits
  writefits,  strtrim(i,1) + '_Mag_' + strtrim(j,1) + '_' + strtrim(k,1) + '.fits', img, head
write_csv, 'ymax.dat', ymax
endfor
endfor
endfor

;align to east subfield
for i=0,8 do begin
for j=0,1 do begin
for k=0,4 do begin

 img=readfits( strtrim(i,1) + '_Mag_' + strtrim(j,1) + '_' + strtrim(k,1) + '.fits', oldhead)
;read in east reference image
   im=readfits("ch1_epoch1.skymap.00.rot.e.fits", refhead)

;align the old flux image with new subfield template east image
   hastrom, img, oldhead, neweastimg, neweasthead, refhead, missing=0

;write new east image with new header
   writefits, 'alignedeastimage' + strtrim(i,1) + '_Mag_' + strtrim(j,1) + '_' + strtrim(k,1) + '.fits', neweastimg, neweasthead


ya=ymax(i,j,k)
xyxy, oldhead, refhead, 0, ya, xb, ybe
ymaxbe(i,j,k)=ybe


;read in west reference image
   im=readfits("ch1_epoch1.skymap.00.rot.w.fits", refhead)

;align the old flux image with new subfield template west image
   hastrom, img, oldhead, newwestimg, newwesthead, refhead, missing=0

;write new west image with new headr
   writefits, 'alignedwestimage' +  strtrim(i,1) + '_Mag_' + strtrim(j,1) + '_' + strtrim(k,1) + '.fits', newwestimg, newwesthead

xyxy, oldhead, refhead, 0, ya, xb, ybw
ymaxbw(i,j,k)=ybw
endfor
endfor
endfor
write_csv, 'ymaxbw.dat', ymaxbw
write_csv, 'ymaxbe.dat', ymaxbe

;crop images
for i=0,8 do begin


for j=0,1 do begin
for k=0,4 do begin

im=readfits('alignedeastimage' +  strtrim(i,1) + '_Mag_' + strtrim(j,1) + '_' + strtrim(k,1) + '.fits', head)

ymaxe=ymaxbe(i,j,k)
if ymaxe LE 3499 then begin
hextract, im, head, newim, newhead, 0, 1699, 0, ymaxe
endif else begin

hextract, im, head, newim, newhead, 0, 1699, 0, 3499
endelse

writefits, 'cropalignedeastimage' +  strtrim(i,1) + '_Mag_' + strtrim(j,1) + '_' + strtrim(k,1) + '.fits', newim, newhead



im=readfits('alignedwestimage' +  strtrim(i,1) + '_Mag_' + strtrim(j,1) + '_' + strtrim(k,1) + '.fits', head)

ymaxw=ymaxbw(i,j,k)

if ymaxw LT 3499 then begin
hextract, im, head, newim, newhead, 0, 1849, 0, ymaxw
endif else begin

hextract, im, head, newim, newhead, 0, 1849, 0, 3499
endelse
writefits, 'cropalignedwestimage' +  strtrim(i,1) + '_Mag_' + strtrim(j,1) + '_' + strtrim(k,1) + '.fits', newim, newhead


endfor
endfor
endfor




;align mask and then apply mask

for i=0,8 do begin

for j=0,1 do begin
for k=0,4 do begin


im=readfits('cropalignedeastimage' + strtrim(i,1) + '_Mag_' + strtrim(j,1) + '_' + strtrim(k,1) + '.fits', refhead)
img=readfits('ch1_epoch1.mask.00.fits', oldhead)

hastrom, img, oldhead, eastmaskimg, eastmaskhead, refhead, missing=0

writefits,  'eastmask'+ strtrim(i,1) + '_Mag_' + strtrim(j,1) + '_' + strtrim(k,1)+'.fits', eastmaskimg, eastmaskhead



im=readfits('cropalignedwestimage' + strtrim(i,1) + '_Mag_' + strtrim(j,1) + '_' + strtrim(k,1) + '.fits', refhead)
img=readfits('ch1_epoch1.mask.00.fits', oldhead)

hastrom, img, oldhead, westmaskimg, westmaskhead, refhead, missing=0

writefits, 'westmask'+ strtrim(i,1) + '_Mag_' + strtrim(j,1) + '_' + strtrim(k,1)+'.fits', westmaskimg, westmaskhead






;read in image
img=readfits('cropalignedwestimage' + strtrim(i,1) + '_Mag_' + strtrim(j,1) + '_' + strtrim(k,1) + '.fits', refhead)
;read in mask
mask=readfits('westmask'+ strtrim(i,1) + '_Mag_' + strtrim(j,1) + '_' + strtrim(k,1)+'.fits', maskhead)

image=mask*img

;write out new masked fit

writefits, 'cropalignedmaskedwest'+ strtrim(i,1) + '_Mag_' + strtrim(j,1) + '_' + strtrim(k,1) + '.fits' , image, refhead



img=readfits('cropalignedeastimage' + strtrim(i,1) + '_Mag_' + strtrim(j,1) + '_' + strtrim(k,1) + '.fits', refhead)
;read in mask
mask=readfits('eastmask'+ strtrim(i,1) + '_Mag_' + strtrim(j,1) + '_' + strtrim(k,1)+'.fits', maskhead)

image=mask*img

;write out new masked fit

writefits, 'cropalignedmaskedeast'+ strtrim(i,1) + '_Mag_' + strtrim(j,1) + '_' + strtrim(k,1) + '.fits' , image, refhead

endfor
endfor
endfor


toc

end
