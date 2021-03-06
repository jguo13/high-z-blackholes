pro walle_HSC
;device,retain=2,decompose=0
;colours
;s_e=size(msk_ch1)
;nx_tot=s_e(1)
;ny_tot=s_e(2)

;	set the field's (x,y) size

q_p_1gridarr=dblarr(29,10,10,10)
power1arr=dblarr(29,10,10,10)
sig_p1arr=dblarr(29,10,10,10)

del1arr=dblarr(29,10,10,10)
sig_del1arr=dblarr(29,10,10,10)




for a=0,3 do begin
for c=0,4 do begin

flux1=readfits('alignedmaskedeastHSC_'+ strtrim(a,1) +'_'+strtrim(c,1)+ '_n.fits',hea2)
cmask=readfits('eastmaskHSC_' + strtrim(a,1) +'_'+strtrim(c,1)+ '_n.fits')
t_ch1=cmask*0+1
flux1=flux1/t_ch1

s_e=size(flux1)
nx_tot=s_e(1)

ny_tot=s_e(2)
sec2rad=!pi/180/3600.
pixel=1.1988; this is the cdelt version converted to arcsec					;	in arcsec
;pixel=1
tpixel=pixel*sec2rad				; in radian
;area=1
area=float(nx_tot)*float(ny_tot)*(pixel/3600.)^2/3282.8
;area=float(nx_tot)*float(ny_tot)*(pixel/3600.)^2*3282.8
;mm=where(cmask LT 1.0 and pippo ge 3.  , count)
;if count gt -1 then cmask(mm)=0
;writefits, 'mask_rep_corr2.fits',cmask, hea5
flux1=flux1*cmask


nx21=nx_tot/2+1
ny21=ny_tot/2+1

n_k=1000
;n_k=400
nx_center=nx21-1  
ny_center=ny21-1  
k_w=fltarr(nx_tot,ny_tot)
exp=fltarr(nx_tot,ny_tot)
rx=fltarr(nx_tot,ny_tot)
media=fltarr(nx_tot,ny_tot)
f1=fltarr(nx_tot,ny_tot)
;k_w=shift(dist(nx_tot,ny_tot),-nx21,-ny21)
k_x=shift(dist(nx_tot,1),-nx21)
k_y=shift(dist(ny_tot,1),-ny21)
for i=0l,nx_tot-1 do k_w(i,*)=sqrt((k_x(i)/nx_tot)^2+(k_y(*)/ny_tot)^2)/pixel
k_0=1./(sqrt(2.)*nx_tot*pixel)	;	0.
k_f=1.01*max(k_w)

tet_1grid=[1.2+10.^(.1*(findgen(30)))]
;tet_1grid=[20l, 400l,1000l]

nq_1grid=n_elements(tet_1grid)
q_p_minmax=2.*!pi/tet_1grid
k_1_minmax=q_p_minmax/2./!pi








q_p_1grid=fltarr(nq_1grid-1)



power1=dblarr(nq_1grid-1)



pairs=dblarr(nq_1grid-1)
cmask(*,0)=0.
cmask(*,1)=0.
mm=where(cmask LT 1.0 , count)
if count gt -1 then cmask(mm)=0
for ik=0,nq_1grid-2 do q_p_1grid(ik)=0.5*(q_p_minmax(ik)+q_p_minmax(ik+1))
f1=flux1

hn0=where(cmask ne 0.,nn)
f_clip_1=float(nn)/nx_tot/ny_tot
weight_1=fltarr(nx_tot,ny_tot)
weight_1(hn0)=t_ch1(hn0)/mean(t_ch1(hn0))
f1(hn0)=f1(hn0)*weight_1(hn0)
f1(hn0)=f1(hn0)-mean(f1(hn0))
amp1=shift(fft(f1,-1,double=1),-nx21,-ny21)


for iq=0,nq_1grid-2 do begin
hp=where(k_w ge k_1_minmax(iq+1) and k_w lt k_1_minmax(iq))
	if(n_elements(hp) gt 1) then begin
power1(iq)=0.5*mean(abs(amp1(hp))^2)*area/f_clip_1
pairs(iq)=n_elements(hp)


	endif
     endfor


sig_p1=(power1/sqrt(0.5*pairs))

;q_p_1gridsqua=(q_p_1grid/sec2rad)^2

del1=(q_p_1grid/sec2rad)^2*power1/2./!pi

sig_del1=(q_p_1grid/sec2rad)^2*sig_p1/2./!pi

print, n_elements(sig_p1)
print, n_elements(del1)
print, n_elements(sig_del1)
print, n_elements(power1)
print, n_elements(q_p_1grid)


sig_p1arr(*,c,a)=sig_p1
del1arr(*,c,a)=sqrt(del1)
sig_del1arr(*,c,a)=sqrt(sig_del1)
power1arr(*,c,a)=power1
q_p_1gridarr(*,c,a)=q_p_1grid

endfor
endfor
sq='sqrt'
radical=TEXTOIDL('q^2P_{2,X}(q)/2\pi')
write_csv, 'power1.dat', power1arr(*,1)
write_csv, 'power2.dat', power1arr(*,3)
write_csv, 'power1arr.dat', power1arr
write_csv, 'q_p_1grid.dat', q_p_1grid
write_csv, 'del1arr.dat', q_p_1grid




set_plot,'PS'

DEVICE, FILE='HSCflucteastmulti.ps', /COLOR, BITS=8,  XSIZE=27, ysize=27
multiplot, /reset
!y.tickname=''
multiplot,[2,2], /square, mytitle=sq+radical, mxtitle='2!7p!5/q  arcsec', /doyaxis
for i=0,3 do begin

;plot_oo,2.*!pi/q_p_1gridarr(*,0),del1arr(*,0),psym=sym(1),xtitle='2!7p!5/q  arcsec',ytitle=TEXTOIDL('q^2P_{2,X}(q)/2\pi'),symsize=1
;errplot,2.*!pi/q_p_1gridarr(*,0),del1arr(*,0)-sig_del1arr(*,0),del1arr(*,0)+sig_del1arr(*,0)
plot_oo,2.*!pi/q_p_1gridarr(*,0,i),del1arr(*,0,i),psym=sym(10),xran=[1.,2000],/xs,yran=[10d-3,10^3],/ys,symsize=1
errplot,2.*!pi/q_p_1gridarr(*,0,i),del1arr(*,0,i)-sig_del1arr(*,0,i),del1arr(*,0,i)+sig_del1arr(*,0,i)

oplot,2.*!pi/q_p_1gridarr(*,1,i),del1arr(*,1,i),psym=sym(2),symsize=1
errplot,2.*!pi/q_p_1gridarr(*,1,i),del1arr(*,1,i)-sig_del1arr(*,1,i),del1arr(*,1,i)+sig_del1arr(*,1,i)

oplot,2.*!pi/q_p_1gridarr(*,2,i),del1arr(*,2,i),psym=sym(6),symsize=1
errplot,2.*!pi/q_p_1gridarr(*,2,i),del1arr(*,2,i)-sig_del1arr(*,2,i),del1arr(*,2,i)+sig_del1arr(*,2,i)

oplot,2.*!pi/q_p_1gridarr(*,3,i),del1arr(*,3,i),psym=sym(4),symsize=1
errplot,2.*!pi/q_p_1gridarr(*,3,i),del1arr(*,3,i)-sig_del1arr(*,3,i),del1arr(*,3,i)+sig_del1arr(*,3,i)

;plot,2.*!pi/q_p_1gridarr(*,4,i),del1arr(*,4,i),psym=sym(5),symsize=1
;errplot,2.*!pi/q_p_1gridarr(*,4,i),del1arr(*,4,i)-sig_del1arr(*,4,i),del1arr(*,4,i)+sig_del1arr(*,4,i)
if i EQ 0 then legend, ['G'], /right
if i EQ 1 then legend, ['Y'],/right
if i EQ 2 then legend, ['I'],/right
if i EQ 3 then legend, ['Z'], /right

multiplot
endfor

device,/close
multiplot,/reset



for a=0,3 do begin
;for b=0,1 do begin
for c=0,4 do begin

flux1=readfits('alignedmaskedwestHSC_'+ strtrim(a,1) +'_'+strtrim(c,1)+ '_n.fits',hea2)
cmask=readfits('westmaskHSC_' + strtrim(a,1) +'_'+strtrim(c,1)+ '_n.fits')
t_ch1=cmask*0+1
flux1=flux1/t_ch1

s_e=size(flux1)
nx_tot=s_e(1)
ny_tot=s_e(2)
sec2rad=!pi/180/3600.
pixel=1.1988; this is the cdelt version converted to arcsec					;	in arcsec
;pixel=1
tpixel=pixel*sec2rad				; in radian
;area=1
area=float(nx_tot)*float(ny_tot)*(pixel/3600.)^2/3282.8
;area=float(nx_tot)*float(ny_tot)*(pixel/3600.)^2*3282.8
;mm=where(cmask LT 1.0 and pippo ge 3.  , count)
;if count gt -1 then cmask(mm)=0
;writefits, 'mask_rep_corr2.fits',cmask, hea5
flux1=flux1*cmask


nx21=nx_tot/2+1
ny21=ny_tot/2+1

n_k=1000
;n_k=400
nx_center=nx21-1  
ny_center=ny21-1  
k_w=fltarr(nx_tot,ny_tot)
exp=fltarr(nx_tot,ny_tot)
rx=fltarr(nx_tot,ny_tot)
media=fltarr(nx_tot,ny_tot)
f1=fltarr(nx_tot,ny_tot)
;k_w=shift(dist(nx_tot,ny_tot),-nx21,-ny21)
k_x=shift(dist(nx_tot,1),-nx21)
k_y=shift(dist(ny_tot,1),-ny21)
for i=0l,nx_tot-1 do k_w(i,*)=sqrt((k_x(i)/nx_tot)^2+(k_y(*)/ny_tot)^2)/pixel
k_0=1./(sqrt(2.)*nx_tot*pixel)	;	0.
k_f=1.01*max(k_w)

tet_1grid=[1.2+10.^(.1*(findgen(30)))]
;tet_1grid=[20l, 400l,1000l]

nq_1grid=n_elements(tet_1grid)
q_p_minmax=2.*!pi/tet_1grid
k_1_minmax=q_p_minmax/2./!pi








q_p_1grid=fltarr(nq_1grid-1)



power1=dblarr(nq_1grid-1)



pairs=dblarr(nq_1grid-1)
cmask(*,0)=0.
cmask(*,1)=0.
mm=where(cmask LT 1.0 , count)
if count gt -1 then cmask(mm)=0
for ik=0,nq_1grid-2 do q_p_1grid(ik)=0.5*(q_p_minmax(ik)+q_p_minmax(ik+1))
f1=flux1

hn0=where(cmask ne 0.,nn)
f_clip_1=float(nn)/nx_tot/ny_tot
weight_1=fltarr(nx_tot,ny_tot)
weight_1(hn0)=t_ch1(hn0)/mean(t_ch1(hn0))
f1(hn0)=f1(hn0)*weight_1(hn0)
;f1(hn0)=f1(hn0)-mean(f1(hn0))
amp1=shift(fft(f1,-1,double=1),-nx21,-ny21)


for iq=0,nq_1grid-2 do begin
hp=where(k_w ge k_1_minmax(iq+1) and k_w lt k_1_minmax(iq))
	if(n_elements(hp) gt 1) then begin
power1(iq)=0.5*mean(abs(amp1(hp))^2)*area/f_clip_1
pairs(iq)=n_elements(hp)


	endif
     endfor


sig_p1=(power1/sqrt(0.5*pairs))

;q_p_1gridsqua=(q_p_1grid/sec2rad)^2

del1=(q_p_1grid/sec2rad)^2*power1/2./!pi

sig_del1=(q_p_1grid/sec2rad)^2*sig_p1/2./!pi

print, n_elements(sig_p1)
print, n_elements(del1)
print, n_elements(sig_del1)
print, n_elements(power1)
print, n_elements(q_p_1grid)


sig_p1arr(*,c,a)=sig_p1
del1arr(*,c,a)=sqrt(del1)
sig_del1arr(*,c,a)=sqrt(sig_del1)
power1arr(*,c,a)=power1
q_p_1gridarr(*,c,a)=q_p_1grid

endfor
endfor
sq='sqrt'
radical=TEXTOIDL('q^2P_{2,X}(q)/2\pi')
write_csv, 'power1.dat', power1arr(*,1)
write_csv, 'power2.dat', power1arr(*,3)
write_csv, 'power1arr.dat', power1arr
write_csv, 'q_p_1grid.dat', q_p_1grid
write_csv, 'del1arr.dat', q_p_1grid




set_plot,'PS'
DEVICE, FILE='HSCfluctwestmulti.ps', /COLOR, BITS=8,  XSIZE=27, ysize=27
!y.tickname=''
multiplot,[2,2], /square, mytitle=sq+radical,/doyaxis,mxtitle='2!7p!5/q  arcsec'
for i=0,3 do begin

;plot_oo,2.*!pi/q_p_1gridarr(*,0),del1arr(*,0),psym=sym(1),xtitle='2!7p!5/q  arcsec',ytitle=TEXTOIDL('q^2P_{2,X}(q)/2\pi'),symsize=1
;errplot,2.*!pi/q_p_1gridarr(*,0),del1arr(*,0)-sig_del1arr(*,0),del1arr(*,0)+sig_del1arr(*,0)
plot_oo,2.*!pi/q_p_1gridarr(*,0,i),del1arr(*,0,i),psym=sym(10),xran=[1.,2000],/xs,yran=[10d-3,10^3],/ys,symsize=1
errplot,2.*!pi/q_p_1gridarr(*,0,i),del1arr(*,0,i)-sig_del1arr(*,0,i),del1arr(*,0,i)+sig_del1arr(*,0,i)

oplot,2.*!pi/q_p_1gridarr(*,1,i),del1arr(*,1,i),psym=sym(2),symsize=1
errplot,2.*!pi/q_p_1gridarr(*,1,i),del1arr(*,1,i)-sig_del1arr(*,1,i),del1arr(*,1,i)+sig_del1arr(*,1,i)

oplot,2.*!pi/q_p_1gridarr(*,2,i),del1arr(*,2,i),psym=sym(6),symsize=1
errplot,2.*!pi/q_p_1gridarr(*,2,i),del1arr(*,2,i)-sig_del1arr(*,2,i),del1arr(*,2,i)+sig_del1arr(*,2,i)

oplot,2.*!pi/q_p_1gridarr(*,3,i),del1arr(*,3,i),psym=sym(4),symsize=1
errplot,2.*!pi/q_p_1gridarr(*,3,i),del1arr(*,3,i)-sig_del1arr(*,3,i),del1arr(*,3,i)+sig_del1arr(*,3,i)

;plot,2.*!pi/q_p_1gridarr(*,4,i),del1arr(*,4,i),psym=sym(5),symsize=1
;errplot,2.*!pi/q_p_1gridarr(*,4,i),del1arr(*,4,i)-sig_del1arr(*,4,i),del1arr(*,4,i)+sig_del1arr(*,4,i)
if i EQ 0 then legend, ['G'], /right
if i EQ 1 then legend, ['Y'],/right
if i EQ 2 then legend, ['I'],/right
if i EQ 3 then legend, ['Z'], /right

multiplot
endfor
device,/close
multiplot,/reset



end
