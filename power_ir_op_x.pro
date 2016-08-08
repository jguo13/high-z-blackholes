pro roz

;create arrays to store the magnitude interval data
q_p_1gridarr=dblarr(29,10,10)
sig_pxarr=dblarr(29,10,10)
sig_pirarr=dblarr(29,10,10)
sig_poparr=dblarr(29,10,10)
delirarr=dblarr(29,10,10)

delabarr=dblarr(29,10,10)

deloparr=dblarr(29,10,10)
delxarr=dblarr(29,10,10)
sig_delxarr=dblarr(29,10,10)
sig_delirarr=dblarr(29,10,10)
sig_deloparr=dblarr(29,10,10)
coherencearr=dblarr(29,10,10)


for a=0,3 do begin
for c=0,4 do begin

;fluxir= yanxia's ir map
fluxir=readfits('ir.fits', irhead)
;mask is the ir mask
mask=readfits('eastmaskHSC_'+strtrim(a,1)+'_'+strtrim(c,1)+'_n.fits', maskhead)
;flux from optical galaxy
fluxop=readfits('alignedmaskedeastHSC_'+strtrim(a,1)+'_'+strtrim(c,1)+'_n.fits',ophead)
;a-b map
delta_f=readfits('ch1_epoch1.skymap.00.rot.e.a_b.fits', abhead)


;;make sure all images are same size
hastrom, fluxop, ophead, newop, newophead, irhead, missing=0

fluxop=newop

hastrom, delta_f, abhead, newab, newabhead, irhead, missing=0

delta_f=newab

hastrom, mask, maskhead, newmask, newmaskhead, irhead, missing=0

mask=newmask


; convert the flux and ab maps from MJy/Sr to nw/m^2/sr
fluxir=fluxir*1e-3*1e9*1e-26*3d8/(3.6e-6)
delta_f=delta_f*1e-3*1e9*1e-26*3d8/(3.6e-6)



fluxir=fluxir*mask
fluxop=fluxop*mask
delta_f=delta_f*mask

t_ch1=mask*0+1
;fluxir=fluxir/t_ch1
s_e=size(mask)
nx_tot=s_e(1)
ny_tot=s_e(2)
sec2rad=!pi/180/3600.
pixel=1.1988;arcsec/pixel
pix2rad=pixel*sec2rad	

			; in radian
;area=1
area=float(nx_tot)*float(ny_tot)*(pixel/3600.)^2/3282.8
;area=float(nx_tot)*float(ny_tot)*(pixel/3600.)^2*3282.8
;mm=where(cmask LT 1.0 and pippo ge 3.  , count)
;if count gt -1 then cmask(mm)=0
;writefits, 'mask_rep_corr2.fits',cmask, hea5


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



;;;;;DONT NEED AXIS MASKING?
;mask(*,0)=0.
;mask(*,1)=0.
mm=where(mask LT 1.0 , count)
if count gt -1 then mask(mm)=0
for ik=0,nq_1grid-2 do q_p_1grid(ik)=0.5*(q_p_minmax(ik)+q_p_minmax(ik+1))


hn0=where(mask ne 0.,nn)
f_clip_1=float(nn)/nx_tot/ny_tot
weight_1=fltarr(nx_tot,ny_tot)
weight_1(hn0)=t_ch1(hn0)/mean(t_ch1(hn0))



fluxir(hn0)=fluxir(hn0)*weight_1(hn0)
fluxir(hn0)=fluxir(hn0)-mean(fluxir(hn0))
fluxop(hn0)=fluxop(hn0)*weight_1(hn0)
fluxop(hn0)=fluxop(hn0)-mean(fluxop(hn0))
delta_f(hn0)=delta_f(hn0)*weight_1(hn0)
delta_f(hn0)=delta_f(hn0)-mean(delta_f(hn0))

;FFT of ir map, optical map, a-b map, and cross images between ir and optical
ampir=shift(fft(fluxir,-1,double=1),-nx21,-ny21)
ampop=shift(fft(fluxop,-1,double=1),-nx21,-ny21)
ampab=shift(fft(delta_f,-1,double=1),-nx21,-ny21)
ampx=(real_part(ampir)*real_part(ampop)+imaginary(ampir)*imaginary(ampop))

;define variables

delir=dblarr(nq_1grid-1)
delop=dblarr(nq_1grid-1)
delx=dblarr(nq_1grid-1)
coherence=dblarr(nq_1grid-1)

powerx=dblarr(nq_1grid-1)
powerir=dblarr(nq_1grid-1)
powerop=dblarr(nq_1grid-1)
powerab_auto=dblarr(nq_1grid-1)

pairs=dblarr(nq_1grid-1)

; compute the power
for iq=0,nq_1grid-2 do begin
hp=where(k_w ge k_1_minmax(iq+1) and k_w lt k_1_minmax(iq))
	if(n_elements(hp) gt 1) then begin
powerx(iq)=mean(abs(ampx(hp)))*area/f_clip_1
powerir(iq)=mean(abs(ampir(hp))^2)*area/f_clip_1
powerop(iq)=mean(abs(ampop(hp))^2)*area/f_clip_1
powerab_auto(iq)=mean(abs(ampab(hp))^2)*area/f_clip_1
pairs(iq)=n_elements(hp)
	endif
     endfor

delir=(q_p_1grid/sec2rad)^2*powerir/2./!pi
delab=(q_p_1grid/sec2rad)^2*powerab_auto/2./!pi

powerir=powerir-powerab_auto

sig_px=(powerx/sqrt(0.5*pairs))
sig_pir=(powerir/sqrt(0.5*pairs))
sig_pop=(powerop/sqrt(0.5*pairs))


delx=(q_p_1grid/sec2rad)^2*powerx/2./!pi
delop=(q_p_1grid/sec2rad)^2*powerop/2./!pi

sig_delx=(q_p_1grid/sec2rad)^2*sig_px/2./!pi
sig_delir=(q_p_1grid/sec2rad)^2*sig_pir/2./!pi
sig_delop=(q_p_1grid/sec2rad)^2*sig_pop/2./!pi

coherence=(abs(powerx))^2/abs(powerir*powerop)


delabarr(*,c,a)=sqrt(delab)


delxarr(*,c,a)=sqrt(delx)
delirarr(*,c,a)=sqrt(delir)
deloparr(*,c,a)=sqrt(delop)

sig_pxarr(*,c,a)=sig_px
sig_pirarr(*,c,a)=sig_pir
sig_poparr(*,c,a)=sig_pop

sig_delirarr(*,c,a)=sqrt(sig_delir)
sig_deloparr(*,c,a)=sqrt(sig_delop)
sig_delxarr(*,c,a)=sqrt(sig_delx)

coherencearr(*,c,a)=coherence
q_p_1gridarr(*,c,a)=q_p_1grid

write_csv, 'delx.csv', delx
write_csv, 'ampx.dat', ampx
write_csv, 'powerir.dat',powerir
write_csv, 'powerx.dat', powerx
write_csv, 'powerop.dat',powerop

endfor
endfor
write_csv, 'coherence.dat', coherencearr

sq='sqrt'
radical=TEXTOIDL('q^2P_{2,X}(q)/2\pi')

set_plot,'PS'
cgerase
DEVICE, FILE='HSC_east_crosscorr.ps', /COLOR, BITS=8,  XSIZE=27, ysize=27
multiplot, /reset
!y.tickname=''
multiplot,[2,2], /square, mytitle=sq+radical, mxtitle='2!7p!5/q  arcsec', /doyaxis
for i=0,3 do begin

;plot_oo,2.*!pi/q_p_1gridarr(*,0),del1arr(*,0),psym=sym(1),xtitle='2!7p!5/q  arcsec',ytitle=TEXTOIDL('q^2P_{2,X}(q)/2\pi'),symsize=1
;errplot,2.*!pi/q_p_1gridarr(*,0),del1arr(*,0)-sig_del1arr(*,0),del1arr(*,0)+sig_del1arr(*,0)

plot_oo,2.*!pi/q_p_1gridarr(*,0,i),delxarr(*,0,i),psym=sym(6),xtitle='2!7p!5/q  arcsec',ytitle=TEXTOIDL('q^2P_{2,X}(q)/2\pi'),symsize=1,yran=[10d-8,10d-4],/ys;,xran=[10,1000],/xs;,yran=[10d-3,10^2],/ys,symsize=.35
errplot,2.*!pi/q_p_1gridarr(*,0,i),delxarr(*,0,i)-sig_delxarr(*,0,i),delxarr(*,0,i)+sig_delxarr(*,0,i)

oplot,2.*!pi/q_p_1gridarr(*,1,i),delxarr(*,1,i),psym=sym(2),symsize=1
errplot,2.*!pi/q_p_1gridarr(*,1,i),delxarr(*,1,i)-sig_delxarr(*,1,i),delxarr(*,1,i)+sig_delxarr(*,1,i)

oplot,2.*!pi/q_p_1gridarr(*,2,i),delxarr(*,2,i),psym=sym(6),symsize=1
errplot,2.*!pi/q_p_1gridarr(*,2,i),delxarr(*,2,i)-sig_delxarr(*,2,i),delxarr(*,2,i)+sig_delxarr(*,2,i)

oplot,2.*!pi/q_p_1gridarr(*,3,i),delxarr(*,3,i),psym=sym(4),symsize=1
errplot,2.*!pi/q_p_1gridarr(*,3,i),delxarr(*,3,i)-sig_delxarr(*,3,i),delxarr(*,3,i)+sig_delxarr(*,3,i)

;oplot,2.*!pi/q_p_1gridarr(*,4,i),del1arr(*,4,i),psym=sym(5),symsize=1
;errplot,2.*!pi/q_p_1gridarr(*,4,i),del1arr(*,4,i)-sig_del1arr(*,4,i),del1arr(*,4,i)+sig_del1arr(*,4,i)



if i EQ 0 then legend, ['G'], /right
if i EQ 1 then legend, ['Y'],/right
if i EQ 2 then legend, ['I'],/right
if i EQ 3 then legend, ['Z'], /right

multiplot
endfor
device,/close
multiplot,/reset


set_plot,'PS'
cgerase
DEVICE, FILE='IR_autopower.ps', /COLOR, BITS=8,  XSIZE=27, ysize=27
multiplot, /reset
!y.tickname=''
multiplot,[2,2], /square, mytitle=sq+radical, mxtitle='2!7p!5/q  arcsec', /doyaxis

for i=0,3 do begin

;plot_oo,2.*!pi/q_p_1gridarr(*,0),del1arr(*,0),psym=sym(1),xtitle='2!7p!5/q  arcsec',ytitle=TEXTOIDL('q^2P_{2,X}(q)/2\pi'),symsize=1
;errplot,2.*!pi/q_p_1gridarr(*,0),del1arr(*,0)-sig_del1arr(*,0),del1arr(*,0)+sig_del1arr(*,0)

plot_oo,2.*!pi/q_p_1gridarr(*,0,i),delirarr(*,0,i),psym=sym(6),xtitle='2!7p!5/q  arcsec',ytitle=TEXTOIDL('q^2P_{2,X}(q)/2\pi'),symsize=1,yran=[1e-11,1e-9],/ys;,xran=[10,1000],/xs;,yran=[10d-3,10^2],/ys,symsize=.35
errplot,2.*!pi/q_p_1gridarr(*,0,i),delirarr(*,0,i)-sig_delirarr(*,0,i),delirarr(*,0,i)+sig_delirarr(*,0,i)

oplot,2.*!pi/q_p_1gridarr(*,1,i),delabarr(*,1,i),psym=sym(2),symsize=1
errplot,2.*!pi/q_p_1gridarr(*,1,i),delirarr(*,1,i)-sig_delirarr(*,1,i),delirarr(*,1,i)+sig_delirarr(*,1,i)

oplot,2.*!pi/q_p_1gridarr(*,2,i),delirarr(*,2,i),psym=sym(6),symsize=1
errplot,2.*!pi/q_p_1gridarr(*,2,i),delirarr(*,2,i)-sig_delirarr(*,2,i),delirarr(*,2,i)+sig_delirarr(*,2,i)

oplot,2.*!pi/q_p_1gridarr(*,3,i),delirarr(*,3,i),psym=sym(4),symsize=1
errplot,2.*!pi/q_p_1gridarr(*,3,i),delirarr(*,3,i)-sig_delirarr(*,3,i),delirarr(*,3,i)+sig_delirarr(*,3,i)

;oplot,2.*!pi/q_p_1gridarr(*,4,i),del1arr(*,4,i),psym=sym(5),symsize=1
;errplot,2.*!pi/q_p_1gridarr(*,4,i),del1arr(*,4,i)-sig_del1arr(*,4,i),del1arr(*,4,i)+sig_del1arr(*,4,i)



if i EQ 0 then legend, ['G'], /right
if i EQ 1 then legend, ['Y'],/right
if i EQ 2 then legend, ['I'],/right
if i EQ 3 then legend, ['Z'], /right

multiplot
endfor
device,/close
multiplot,/reset



set_plot,'PS'
cgerase
DEVICE, FILE='HSC_east_coherence.ps', /COLOR, BITS=8,  XSIZE=27, ysize=27
multiplot, /reset
!y.tickname=''
multiplot,[2,2], /square, mytitle=sq+radical, mxtitle='2!7p!5/q  arcsec', /doyaxis

for i=0,3 do begin

plot_oo,2.*!pi/q_p_1gridarr(*,0,i),coherencearr(*,0,i),psym=sym(10),xtitle='2!7p!5/q  arcsec',ytitle=TEXTOIDL('q^2P_{2,X}(q)/2\pi'),symsize=1,yran=[10e-2,1],/ys;,xran=[10,1000],/xs;,symsize=.35
;;;errplot,2.*!pi/q_p_1gridarr(*,0,i),del1arr(*,0,i)-sig_del1arr(*,0,i),del1arr(*,0,i)+sig_del1arr(*,0,i)

oplot,2.*!pi/q_p_1gridarr(*,1,i),coherencearr(*,1,i),psym=sym(2),symsize=1
;;;errplot,2.*!pi/q_p_1gridarr(*,1,i),del1arr(*,1,i)-sig_del1arr(*,1,i),del1arr(*,1,i)+sig_del1arr(*,1,i)


oplot,2.*!pi/q_p_1gridarr(*,2,i),coherencearr(*,2,i),psym=sym(6),symsize=1
;;;errplot,2.*!pi/q_p_1gridarr(*,2,i),del1arr(*,2,i)-sig_del1arr(*,2,i),del1arr(*,2,i)+sig_del1arr(*,2,i)


oplot,2.*!pi/q_p_1gridarr(*,3,i),coherencearr(*,3,i),psym=sym(4),symsize=1
;;;errplot,2.*!pi/q_p_1gridarr(*,3,i),del1arr(*,3,i)-sig_del1arr(*,3,i),del1arr(*,3,i)+sig_del1arr(*,3,i)


;oplot,2.*!pi/q_p_1gridarr(*,4,i),coherence(*,4,i),psym=sym(5),symsize=1
;errplot,2.*!pi/q_p_1gridarr(*,4,i),del1arr(*,4,i)-sig_del1arr(*,4,i),del1arr(*,4,i)+sig_del1arr(*,4,i)


if i EQ 0 then legend, ['G'], /right
if i EQ 1 then legend, ['Y'],/right
if i EQ 2 then legend, ['I'],/right
if i EQ 3 then legend, ['Z'], /right

multiplot
endfor
device,/close
multiplot,/reset



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;WEST IMAGE;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

for a=0,3 do begin
for c=0,4 do begin


fluxir=readfits('ch1_epoch1.rot.w.fluc.cut.fits', irhead)
mask=readfits('westmaskHSC_'+strtrim(a,1)+'_'+strtrim(c,1)+'_n.fits', maskhead)
fluxop=readfits('alignedmaskedwestHSC_'+strtrim(a,1)+'_'+strtrim(c,1)+'_n.fits',ophead)
delta_f=readfits('ch1_epoch1.skymap.00.rot.w.a_b.fits', abhead)


;;make sure all images are same size
hastrom, fluxop, ophead, newop, newophead, irhead, missing=0

fluxop=newop

hastrom, delta_f, abhead, newab, newabhead, irhead, missing=0

delta_f=newab

hastrom, mask, maskhead, newmask, newmaskhead, irhead, missing=0

mask=newmask



fluxir=fluxir*1e-3*1e9*1e-26*3d8/(3.6e-6)
delta_f=delta_f*1e-3*1e9*1e-26*3d8/(3.6e-6)



fluxir=fluxir*mask
fluxop=fluxop*mask
delta_f=delta_f*mask

t_ch1=mask*0+1
;fluxir=fluxir/t_ch1
s_e=size(mask)
nx_tot=s_e(1)
ny_tot=s_e(2)
sec2rad=!pi/180/3600.
pixel=1.1988; this is the cdelt version converted to arcsec				;	in arcsec
;pixel=1
pix2rad=pixel*sec2rad	

			; in radian
;area=1
area=float(nx_tot)*float(ny_tot)*(pixel/3600.)^2/3282.8
;area=float(nx_tot)*float(ny_tot)*(pixel/3600.)^2*3282.8
;mm=where(cmask LT 1.0 and pippo ge 3.  , count)
;if count gt -1 then cmask(mm)=0
;writefits, 'mask_rep_corr2.fits',cmask, hea5


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



;;;;;DONT NEED AXIS MASKING?
mask(*,0)=0.
mask(*,1)=0.
mm=where(mask LT 1.0 , count)
if count gt -1 then mask(mm)=0
for ik=0,nq_1grid-2 do q_p_1grid(ik)=0.5*(q_p_minmax(ik)+q_p_minmax(ik+1))


hn0=where(mask ne 0.,nn)
f_clip_1=float(nn)/nx_tot/ny_tot
weight_1=fltarr(nx_tot,ny_tot)
weight_1(hn0)=t_ch1(hn0)/mean(t_ch1(hn0))

;fluxir(hn0)=fluxir(hn0)*weight_1(hn0)
;fluxir(hn0)=fluxir(hn0)-mean(fluxir(hn0))
;fluxop(hn0)=fluxop(hn0)*weight_1(hn0)
;fluxop(hn0)=fluxop(hn0)-mean(fluxop(hn0))
;delta_f(hn0)=delta_f(hn0)*weight_1(hn0)
;delta_f(hn0)=delta_f(hn0)-mean(delta_f(hn0))

ampir=shift(fft(fluxir,-1,double=1),-nx21,-ny21)
ampop=shift(fft(fluxop,-1,double=1),-nx21,-ny21)
ampab=shift(fft(delta_f,-1,double=1),-nx21,-ny21)
ampx=(real_part(ampir)*real_part(ampop)+imaginary(ampir)*imaginary(ampop))

;define variables

delir=dblarr(nq_1grid-1)
delop=dblarr(nq_1grid-1)
delx=dblarr(nq_1grid-1)
coherence=dblarr(nq_1grid-1)

powerx=dblarr(nq_1grid-1)
powerir=dblarr(nq_1grid-1)
powerop=dblarr(nq_1grid-1)
powerab_auto=dblarr(nq_1grid-1)

pairs=dblarr(nq_1grid-1)


for iq=0,nq_1grid-2 do begin
hp=where(k_w ge k_1_minmax(iq+1) and k_w lt k_1_minmax(iq))
	if(n_elements(hp) gt 1) then begin
powerx(iq)=mean(abs(ampx(hp)))*area/f_clip_1
powerir(iq)=mean(abs(ampir(hp))^2)*area/f_clip_1
powerop(iq)=mean(abs(ampop(hp))^2)*area/f_clip_1
powerab_auto(iq)=mean(abs(ampab(hp))^2)*area/f_clip_1
pairs(iq)=n_elements(hp)
	endif
     endfor


powerir=powerir-powerab_auto

sig_px=(powerx/sqrt(0.5*pairs))
sig_pir=(powerir/sqrt(0.5*pairs))
sig_pop=(powerop/sqrt(0.5*pairs))


delx=(q_p_1grid/sec2rad)^2*powerx/2./!pi
delir=(q_p_1grid/sec2rad)^2*powerir/2./!pi
delop=(q_p_1grid/sec2rad)^2*powerop/2./!pi

sig_delx=(q_p_1grid/sec2rad)^2*sig_px/2./!pi
sig_delir=(q_p_1grid/sec2rad)^2*sig_pir/2./!pi
sig_delop=(q_p_1grid/sec2rad)^2*sig_pop/2./!pi

coherence=(abs(powerx))^2/abs(powerir*powerop)




delxarr(*,c,a)=sqrt(delx)
delirarr(*,c,a)=sqrt(delir)
deloparr(*,c,a)=sqrt(delop)

sig_pxarr(*,c,a)=sig_px
sig_pirarr(*,c,a)=sig_pir
sig_poparr(*,c,a)=sig_pop

sig_delirarr(*,c,a)=sqrt(sig_delir)
sig_deloparr(*,c,a)=sqrt(sig_delop)
sig_delxarr(*,c,a)=sqrt(sig_delx)

coherencearr(*,c,a)=coherence
q_p_1gridarr(*,c,a)=q_p_1grid

write_csv, 'delx.csv', delx
write_csv, 'ampx.dat', ampx
write_csv, 'powerir.dat',powerir
write_csv, 'powerx.dat', powerx
write_csv, 'powerop.dat',powerop

endfor
endfor
write_csv, 'coherence.dat', coherencearr

sq='sqrt'
radical=TEXTOIDL('q^2P_{2,X}(q)/2\pi')

set_plot,'PS'
cgerase
DEVICE, FILE='HSC_west_crosscorr.ps', /COLOR, BITS=8,  XSIZE=27, ysize=27
multiplot, /reset
!y.tickname=''
multiplot,[2,2], /square, mytitle=sq+radical, mxtitle='2!7p!5/q  arcsec', /doyaxis

for i=0,3 do begin

;plot_oo,2.*!pi/q_p_1gridarr(*,0),del1arr(*,0),psym=sym(1),xtitle='2!7p!5/q  arcsec',ytitle=TEXTOIDL('q^2P_{2,X}(q)/2\pi'),symsize=1
;errplot,2.*!pi/q_p_1gridarr(*,0),del1arr(*,0)-sig_del1arr(*,0),del1arr(*,0)+sig_del1arr(*,0)

plot_oo,2.*!pi/q_p_1gridarr(*,0,i),delxarr(*,0,i),psym=sym(6),xtitle='2!7p!5/q  arcsec',ytitle=TEXTOIDL('q^2P_{2,X}(q)/2\pi'),symsize=1,yran=[10d-8,10d-4],/ys;,xran=[10,1000],/xs;,yran=[10d-3,10^2],/ys,symsize=.35
errplot,2.*!pi/q_p_1gridarr(*,0,i),delxarr(*,0,i)-sig_delxarr(*,0,i),delxarr(*,0,i)+sig_delxarr(*,0,i)

oplot,2.*!pi/q_p_1gridarr(*,1,i),delxarr(*,1,i),psym=sym(2),symsize=1
errplot,2.*!pi/q_p_1gridarr(*,1,i),delxarr(*,1,i)-sig_delxarr(*,1,i),delxarr(*,1,i)+sig_delxarr(*,1,i)

oplot,2.*!pi/q_p_1gridarr(*,2,i),delxarr(*,2,i),psym=sym(6),symsize=1
errplot,2.*!pi/q_p_1gridarr(*,2,i),delxarr(*,2,i)-sig_delxarr(*,2,i),delxarr(*,2,i)+sig_delxarr(*,2,i)

oplot,2.*!pi/q_p_1gridarr(*,3,i),delxarr(*,3,i),psym=sym(4),symsize=1
errplot,2.*!pi/q_p_1gridarr(*,3,i),delxarr(*,3,i)-sig_delxarr(*,3,i),delxarr(*,3,i)+sig_delxarr(*,3,i)

;oplot,2.*!pi/q_p_1gridarr(*,4,i),del1arr(*,4,i),psym=sym(5),symsize=1
;errplot,2.*!pi/q_p_1gridarr(*,4,i),del1arr(*,4,i)-sig_del1arr(*,4,i),del1arr(*,4,i)+sig_del1arr(*,4,i)



if i EQ 0 then legend, ['G'], /right
if i EQ 1 then legend, ['Y'],/right
if i EQ 2 then legend, ['I'],/right
if i EQ 3 then legend, ['Z'], /right

multiplot
endfor
device,/close
multiplot,/reset






set_plot,'PS'
cgerase
DEVICE, FILE='HSC_west_coherence.ps', /COLOR, BITS=8,  XSIZE=27, ysize=27
multiplot, /reset
!y.tickname=''
multiplot,[2,2], /square, mytitle=sq+radical, mxtitle='2!7p!5/q  arcsec', /doyaxis

for i=0,3 do begin

plot_oo,2.*!pi/q_p_1gridarr(*,0,i),coherencearr(*,0,i),psym=sym(10),xtitle='2!7p!5/q  arcsec',ytitle=TEXTOIDL('q^2P_{2,X}(q)/2\pi'),symsize=1,yran=[10e-2,1],/ys;,xran=[10,1000],/xs;,symsize=.35
;;;errplot,2.*!pi/q_p_1gridarr(*,0,i),del1arr(*,0,i)-sig_del1arr(*,0,i),del1arr(*,0,i)+sig_del1arr(*,0,i)

oplot,2.*!pi/q_p_1gridarr(*,1,i),coherencearr(*,1,i),psym=sym(2),symsize=1
;;;errplot,2.*!pi/q_p_1gridarr(*,1,i),del1arr(*,1,i)-sig_del1arr(*,1,i),del1arr(*,1,i)+sig_del1arr(*,1,i)


oplot,2.*!pi/q_p_1gridarr(*,2,i),coherencearr(*,2,i),psym=sym(6),symsize=1
;;;errplot,2.*!pi/q_p_1gridarr(*,2,i),del1arr(*,2,i)-sig_del1arr(*,2,i),del1arr(*,2,i)+sig_del1arr(*,2,i)


oplot,2.*!pi/q_p_1gridarr(*,3,i),coherencearr(*,3,i),psym=sym(4),symsize=1
;;;errplot,2.*!pi/q_p_1gridarr(*,3,i),del1arr(*,3,i)-sig_del1arr(*,3,i),del1arr(*,3,i)+sig_del1arr(*,3,i)


;oplot,2.*!pi/q_p_1gridarr(*,4,i),coherence(*,4,i),psym=sym(5),symsize=1
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
