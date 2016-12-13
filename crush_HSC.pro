pro crush_HSC
tic
;creates the HSC images
for i=0,3 do begin
     for k=0,4 do begin

FMT='X,X,X,X,X,X,X,X,X,F,F,F,F'

readcol,'HSC_' + strtrim(i,1) + '_n.dat', F=FMT, raold, decold, mag, fluxold

;start the magnitude limit inputs
if i EQ 0 then m0=26.5-6
if i EQ 1 then m0=26.7-6
if i EQ 2 then m0=26.2-6
if i EQ 3 then m0=24.6-6


;start the increments
if k EQ 0 then mlim=m0
if k EQ 1 then mlim=m0+2
if k EQ 2 then mlim=m0+4
if k EQ 3 then mlim=m0+6
if k EQ 4 then begin
mlim=0
endif

maglim=m0+6


;;;;;;



magindex= WHERE(mag GE mlim, /NULL)
magarr=mag(magindex)

ra=raold(magindex)
dec=decold(magindex)
flux=fluxold(magindex)

;;;;;Cutting past the maglim
magindex= WHERE(mag LE maglim, /NULL);
mag=mag(magindex)

ra=ra(magindex)
dec=dec(magindex)
flux=flux(magindex)


;read in COSMOS field template, which is a cosmos field image but all
;pixels have a value of zero
img=readfits("template.fits",head)
adxy, head, ra, dec, x0, y0

x0=round(x0)
y0=round(y0)

ns= n_elements(flux)
;get the degrees per pixel count from header
cdelt=sxpar(head,'CDELT1')
cdeltsquare= abs(cdelt)^2

;input in the wavelength of the bands in meters
if i EQ 0 then lambda=4750e-10

if i EQ 1 then lambda=1050e-10
if i EQ 2 then lambda=7700e-10
if i EQ 3 then lambda=9105e-10

for h=0,ns-1 do begin
;check to see if there are flux values less than -99  because these
;mean no data. if there are, add in 0 into the template coordinate
  if (flux[h] LE -99) then begin
   img[x0[h],y0[h]] += flux[h]*0l
endif else begin
   if x0[h] LE 6322 && x0[h] GE 0 && y0[h] LE 6195 && y0[h] GE 0 then begin
         img[x0[h],y0[h]] += flux[h]*1d6*10d-23*3d8/lambda/(cdeltsquare/3282.8)
     endif
  endelse
endfor


;convert the sigma (based on FWM of IRAC)into pixels using cdelt from header
sigma=(1.6/2.3555)/(3600*abs(cdelt))
;apply gaussian smoothing with sigma
  img=GAUSS_SMOOTH(img, sigma)
;write out new fits image with name redhisftimage+the bottomrange of
;the redshiftrange+.fits
  writefits,  'HSC_' + strtrim(i,1) + '_' + strtrim(k,1)+ '_n.fits', img, head;
;;write_csv, 'ymax.dat', ymax
endfor
endfor

;align the image to east and wes subfield

for i=0,3 do begin

for k=0,4 do begin

img=readfits( 'HSC_' + strtrim(i,1) + '_' + strtrim(k,1)+ '_n.fits', oldhead)
;read in east reference image
   im=readfits("ch1_epoch1.skymap.00.rot.e.fits", refhead)

;align the old flux image with new subfield template east image
   hastrom, img, oldhead, neweastimg, neweasthead, refhead, missing=0

;write new east image with new header
   writefits, 'alignedeastimageHSC_'+strtrim(i,1)+ '_'+ strtrim(k,1)+'_n.fits',neweastimg, neweasthead


img=readfits( 'HSC_' + strtrim(i,1) + '_' + strtrim(k,1)+ '_n.fits', oldhead)
;read in east reference image
   im=readfits("ch1_epoch1.skymap.00.rot.w.fits", refhead)

;align the old flux image with new subfield template east image
   hastrom, img, oldhead, newwestimg, newwesthead, refhead, missing=0

;write new east image with new header
   writefits, 'alignedwestimageHSC_'+strtrim(i,1)+ '_'+ strtrim(k,1)+'_n.fits',newwestimg, newwesthead


endfor
endfor
;;we don't need to crop the HSC data becaseu it's actually bigger than
;;the subfield, so hastrom will take care of everything



; align the mask and apply the mask
for i=0,3 do begin

for k=0,4 do begin


im=readfits( 'alignedeastimageHSC_'+strtrim(i,1)+ '_'+ strtrim(k,1)+'_n.fits', refhead)
;;;;;img=readfits('ch1_epoch1.mask.00.fits', oldhead)
mask=readfits('y_east_mask.fits', oldhead)
hastrom, mask, maskhead, refhead, missing=0l

image=mask*im

writefits, 'alignedmaskedeastHSC_'+ strtrim(i,1) + '_'+strtrim(k,1)+'_n.fits', image, refhead




im=readfits('alignedwestimageHSC_'+strtrim(i,1)+ '_'+ strtrim(k,1)+'_n.fits' , refhead)
img=readfits('ch1_epoch1.mask.00.fits', oldhead)

hastrom, img, oldhead, westmaskimg, westmaskhead, refhead, missing=0

writefits, 'westmaskHSC_' + strtrim(i,1) + '_'+ strtrim(k,1)+'_n.fits', westmaskimg, westmaskhead



;read in image
img=readfits('alignedwestimageHSC_'+strtrim(i,1)+ '_'+ strtrim(k,1)+'_n.fits', refhead)
;read in mask
mask=readfits('westmaskHSC_' + strtrim(i,1) + '_'+ strtrim(k,1)+'_n.fits', maskhead)

image=mask*img

;write out new masked fit

writefits, 'alignedmaskedwestHSC_'+ strtrim(i,1) + '_'+strtrim(k,1)+'_n.fits', image, refhead




endfor

endfor


toc

end

