pro hehe
;;loop through ch1 and ch2 of the IRAC maps
;;0=ch1
;;1=ch2
for t=0,1 do begin
if t EQ 0 then ch='1'
if t EQ 1 then ch='2'
;;loop through east and west field
;;0=east
;;1=west
for u=0,1 do begin
if u EQ 0 then field='e'
if u EQ 1 then field='w'
;create arrays to store the magnitude interval data
;the arrays are sized by: (#of bins, #of mag intervals, #of bands)
q_p_1gridarr=dblarr(12,4,16);stores the wavenumber that will be converted to arcsec and used as x-axis
sig_pxarr=dblarr(12,4,16);stores error of cross power
sig_pirarr=dblarr(12,4,16);stores error of ir (a+b)
sig_poparr=dblarr(12,4,16); stores error of optical + xray
delirarr=dblarr(12,4,16); stores y values to be graphed of autopower of ir
poweroparr=dblarr(12,4,16); stores autopower of optical +xray
powerxarr=dblarr(12,4,16); stores cross power of optical or xray with ir
delabarr=dblarr(12,4,16); stores y values to be graphed of autopower f a-b of ir

deloparr=dblarr(12,4,16); stores y values to be graphed of autopower of optical+xray
delxarr=dblarr(12,4,16); stores y values of cross power to be graphed
sig_delxarr=dblarr(12,4,16);stores error for x power for cross spectrum
sig_delirarr=dblarr(12,4,16); stores error for ir auto power for ir autopower spectrum
sig_deloparr=dblarr(12,4,16); stores error for optical power for optical autopower spectrum
coherencearr=dblarr(12,4,16); stores error for coherence for coherence power



;;;; I named these files from three different
;;;; catalogs/observations(COSMOS catalog, HSC observation from
;;;; profesor Hasinger, and xray observation from nico) in a kinda counterintuitive way so here's
;;;; the key to which band you are working on for the ath iteration:
;;;;a=0=g band (HSC)= 4750 A
;;;;a=1=y band (HSC)= 10150 A 
;;;;a=2=I band (HSC)= 7700 A
;;;;a=3=z band (HSC)= 8875 A
;;;;a=4=r band (COSMOS)= 6288 A
;;;;a=5=I band (COSMOS)= 7683 A 
;;;;a=6=z band (cosmos)= 9105 A
;;;;a=7=y band (cosmos)= 10214 A 
;;;;a=8=Y band (Cosmos)= 9791 A
;;;;a=9=j band (COSMOS)= 12534 A 
;;;;a=10=h band (cosmos)= 16453 A
;;;;a=11=Ks band (cosmos)= 21539 A 
;;;;a=12=Ch2 IRAC band (COSMOS)= 35634 A 
;;;;a=13=0.5-3 Kev xray band (Xray)= 9.91 A
;;;;a=14=2-10 Kev xray band (Xray)= 2.06 A
;;;;a=15=Ch2 IRAC band (COSMOS)= 45110.1 A   

;;;and the c loop loops over different magnitude intervals, with the
;;;oth c loop being the brightest and the 3th c loop being the
;;;faintest (ie at the magnitude limit) however it should be noted
;;;that only the COSMOS and HSC data have magnitude limits, for the
;;;xray data the loop just creates 3dim of the same array in the 13th
;;;and 14th dim of the powerop array 
for a=0,15 do begin
for c=0,3 do begin


;fluxir= yanxia's ir map, this is also the a+b map, later we
;will halve this
fluxir=readfits('ch'+ch+'_'+field+'.final.fits', irhead)
mask=readfits('ch'+ch+'_'+field+'.rot.mask.cut.sex.fits', maskhead)
hastrom, mask, maskhead, irhead, missing=0l
fluxir=fluxir*mask
;read in COSMOS data
   if a GE 4 && a LE 12 then begin
b=a-4
fluxop=readfits(strtrim(b,1)+'_Mag_1_'+strtrim(c,1)+'.fits', ophead)
hastrom,fluxop, ophead, irhead, missing=0l
fluxop=fluxop*mask

endif
;read in COSMOS CH2 data
   if a EQ 15 then begin


fluxop=readfits('8_Mag_0_'+strtrim(c,1)+'.fits', ophead)
hastrom,fluxop, ophead, irhead, missing=0l
fluxop=fluxop*mask


endif
;read in XRAY data
   if a GE 13 && a NE 15 then begin
b=a-13
;xray data needs to be fitted to the east field
fluxop=readfits('xray_'+strtrim(b,1)+'.fits', ophead)
hastrom,fluxop, ophead, irhead, missing=0l
;apply mask
fluxop=fluxop*mask

endif
;read in HSC data
      if a LE 3 then begin
;fluxfrom optical galaxy
fluxop=readfits('HSC_'+strtrim(a,1)+'_'+strtrim(c,1)+'_n.fits',ophead)
hastrom,fluxop, ophead, irhead, missing=0l
;apply mask                                                                                                                                           
fluxop=fluxop*mask
endif


;a-b ir map
delta_f=readfits('ch'+ch+'_'+field+'.a_b.final.fits', abhead)
hastrom,delta_f, abhead, irhead, missing=0l
;apply mask
delta_f=delta_f*mask
;
                                                                                                                                                      


;convert the ir maps from MJy/sr to nw/m^2/sr
;fluxir=fluxir*1e-17*(3e14/3.6)*1e6
;delta_f=delta_f*1e-17*(3e14/3.6)*1e6





t_ch1=mask*0+1
s_e=size(fluxir)
nx_tot=s_e(1)
ny_tot=s_e(2)
;arcsec to radian conversion
sec2rad=!pi/180/3600.
pixel=1.1988;arcsec/pixel
pix2rad=pixel*sec2rad

;calculate the area in radian
area=float(nx_tot)*float(ny_tot)*(pix2rad)^2


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
k_0=1./(sqrt(2.)*nx_tot*pixel)  ;0.
k_f=1.01*max(k_w)
;this creates the binning, where it is bins 0.3 arcsec size with 13
;bins total
tet_1grid=[1.2+10.^(.3*(findgen(13)))]

nq_1grid=n_elements(tet_1grid)
q_p_minmax=2.*!pi/tet_1grid
k_1_minmax=q_p_minmax/2./!pi



q_p_1grid=fltarr(nq_1grid-1)


mm=where(mask LT 1.0 , count)
if count gt -1 then mask(mm)=0
;this creates the x axis, which is eqaul to the middle part of
;the bin it's in 
for ik=0,nq_1grid-2 do q_p_1grid(ik)=0.5*(q_p_minmax(ik)+q_p_minmax(ik+1))

; this calculates % of map unmasked, which will be used to normalize
; data later
hn0=where(mask ne 0.,nn)
f_clip_1=float(nn)/nx_tot/ny_tot
weight_1=fltarr(nx_tot,ny_tot)
weight_1(hn0)=t_ch1(hn0)/mean(t_ch1(hn0))



fluxir(hn0)=fluxir(hn0)-mean(fluxir(hn0))
fluxop(hn0)=fluxop(hn0)-mean(fluxop(hn0))
delta_f(hn0)=delta_f(hn0)-mean(delta_f(hn0))

;FFT of ir map, optical map, a-b map, and cross images between ir and
;optical
ampir=shift(fft(fluxir,-1,double=1),-nx21,-ny21)
ampop=shift(fft(fluxop,-1,double=1),-nx21,-ny21)
ampab=shift(fft(delta_f,-1,double=1),-nx21,-ny21)
ampx=(real_part(ampir)*real_part(ampop)+imaginary(ampir)*imaginary(ampop))


;define variables
delab=dblarr(nq_1grid-1)
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
;powerx(iq)=mean(abs(ampx(hp)))*area/f_clip_1
powerir(iq)=mean(abs(ampir(hp))^2)*area/f_clip_1
powerop(iq)=mean(abs(ampop(hp))^2)*area/f_clip_1
powerab_auto(iq)=mean(abs(ampab(hp))^2)*area/f_clip_1
powerx(iq)=mean(ampx(hp))*area/f_clip_1

;powerir(iq)=mean(ampir(hp)^2)*area/f_clip_1
;powerop(iq)=mean(ampop(hp)^2)*area/f_clip_1
;powerab_auto(iq)=mean(ampab(hp)^2)*area/f_clip_1
pairs(iq)=n_elements(hp)
endif
endfor
powerx=0.25*powerx
powerir=0.5*powerir

delab=(q_p_1grid/sec2rad)^2*powerab_auto/2./!pi
delir=(q_p_1grid/sec2rad)^2*powerir/2./!pi

powerir=powerir-powerab_auto





sig_px=(powerx/sqrt(0.5*pairs))
sig_pir=(powerir/sqrt(0.5*pairs))
sig_pop=(powerop/sqrt(0.5*pairs))


delx=(q_p_1grid/sec2rad)^2*powerx/2./!pi
delop=(q_p_1grid/sec2rad)^2*powerop/2./!pi

sig_delx=(q_p_1grid/sec2rad)^2*sig_px/2./!pi
sig_delir=(q_p_1grid/sec2rad)^2*sig_pir/2./!pi
sig_delop=(q_p_1grid/sec2rad)^2*sig_pop/2./!pi

poweroparr(*,c,a)=powerop
powerxarr(*,c,a)=powerx
;coherence=(abs(powerx))^2/abs(powerir*powerop)
coherence=((powerx))^2/(powerir*powerop)


delxarr(*,c,a)=delx
deloparr(*,c,a)=delop

sig_pxarr(*,c,a)=sig_px
sig_pirarr(*,c,a)=sig_pir
sig_poparr(*,c,a)=sig_pop

sig_delirarr(*,c,a)=sqrt(sig_delir)
sig_deloparr(*,c,a)=sqrt(sig_delop)
sig_delxarr(*,c,a)=sqrt(sig_delx)

coherencearr(*,c,a)=coherence
q_p_1gridarr(*,c,a)=q_p_1grid


endfor
endfor
write_csv, 'delx.csv', delxarr
write_csv, 'delir.csv', delir
write_csv, 'delop.csv', deloparr
write_csv, 'sigx.csv', sig_delxarr
write_csv, 'sigir.csv', sig_delirarr
write_csv, 'sigop.csv', sig_deloparr
write_csv, 'coherence.dat', coherencearr
mus=String(108B)
sqrt=String(83B)
sq='!9'+ sqrt + '!X'
mu='!4'+ mus + '!X'
radical=TEXTOIDL('q^2P_{2,X}(q)/2\pi')
set_plot,'PS'

DEVICE, FILE='ch'+ch+'_'+field+'.fluctmulti.ps', /COLOR, BITS=8,  XSIZE=22, ysize=22, xoffset=-1, yoffset=1,/BOLD
multiplot, /reset
!y.tickname=''
!P.Thick=3
multiplot,[4,3], /square, mytitle=sq+radical, mxtitle='2!7p!5/q  arcsec', /doyaxis,mtitsize=1.4, mytitsize=1.4, mxtitsize=1.4
for i=0,15 do begin
   if i EQ 5 or i EQ 6 or i EQ 7 or i EQ 8 or i EQ 13 or i EQ 14 then continue

;plot_oo,2.*!pi/q_p_1gridarr(*,0),del1arr(*,0),psym=sym(1),xtitle='2!7p!5/q  arcsec',ytitle=TEXTOIDL('q^2P_{2,X}(q)/2\pi'),symsize=1
;errplot,2.*!pi/q_p_1gridarr(*,0),del1arr(*,0)-sig_del1arr(*,0),del1arr(*,0)+sig_del1arr(*,0)
plot_oo,2.*!pi/q_p_1gridarr(*,0,i),deloparr(*,0,i),psym=sym(10),symsize=1,xran=[1.,2000],/xs,yran=[10d-4,10^3],/ys
errplot,2.*!pi/q_p_1gridarr(*,0,i),deloparr(*,0,i)-sig_deloparr(*,0,i),deloparr(*,0,i)+sig_deloparr(*,0,i)

oplot,2.*!pi/q_p_1gridarr(*,1,i),deloparr(*,1,i),psym=sym(2),symsize=1
errplot,2.*!pi/q_p_1gridarr(*,1,i),deloparr(*,1,i)-sig_deloparr(*,1,i),deloparr(*,1,i)+sig_deloparr(*,1,i)

oplot,2.*!pi/q_p_1gridarr(*,2,i),deloparr(*,2,i),psym=sym(6),symsize=1
errplot,2.*!pi/q_p_1gridarr(*,2,i),deloparr(*,2,i)-sig_deloparr(*,2,i),deloparr(*,2,i)+sig_deloparr(*,2,i)

oplot,2.*!pi/q_p_1gridarr(*,3,i),deloparr(*,3,i),psym=sym(4),symsize=1
errplot,2.*!pi/q_p_1gridarr(*,3,i),deloparr(*,3,i)-sig_deloparr(*,3,i),deloparr(*,3,i)+sig_deloparr(*,3,i)

;plot,2.*!pi/q_p_1gridarr(*,4,i),del1arr(*,4,i),psym=sym(5),symsize=1
;errplot,2.*!pi/q_p_1gridarr(*,4,i),del1arr(*,4,i)-sig_del1arr(*,4,i),del1arr(*,4,i)+sig_del1arr(*,4,i)
if i EQ 0 then legend, ['G'], /right
if i EQ 1 then legend, ['Y'],/right ; this is HSC data
if i EQ 2 then legend, ['I'],/right
if i EQ 3 then legend, ['Z'], /right


if i EQ 4 then legend, ['R'], /right
if i EQ 5 then legend, ['I'],/right
if i EQ 6 then legend, ['Z'],/right
if i EQ 7 then legend, ['Y'], /right ;this is the other data from cosmos
if i EQ 8 then legend, ['Y'],/right ; this is the HSC data
if i EQ 9 then legend, ['J'],/right
if i EQ 10 then legend, ['H'],/right
if i EQ 11 then legend, ['Ks'],/right
if i EQ 12 then legend, ['4.5'+mu+'m'],/right
if i EQ 15 then legend, ['3.6'+mu+'m'],/right
multiplot
endfor

device,/close
multiplot,/reset




set_plot,'PS'
cgerase
DEVICE, FILE='ch'+ch+'_'+field+'_crosscorr.ps', /COLOR, BITS=8, XSIZE=22, ysize=22, xoffset=-1, yoffset=0,/BOLD

multiplot, /reset
!y.tickname=''
!P.thick=3
multiplot,[4,4], /square, mytitle=sq+radical, mxtitle='2!7p!5/q  arcsec', /doyaxis,mtitsize=1.4, mytitsize=1.4,mxtitsize=1.4
for i=0,15 do begin
  ; if i EQ 5 or i EQ 6 or i EQ 7 or i EQ 8 or i EQ 13 or i EQ 14 then continue
   if i EQ 5 or i EQ 6 or i EQ 7 or i EQ 8 or i then continue
;plot_oo,2.*!pi/q_p_1gridarr(*,0),del1arr(*,0),psym=sym(1),xtitle='2!7p!5/q  arcsec',ytitle=TEXTOIDL('q^2P_{2,X}(q)/2\pi'),symsize=1
;errplot,2.*!pi/q_p_1gridarr(*,0),del1arr(*,0)-sig_del1arr(*,0),del1arr(*,0)+sig_del1arr(*,0)

plot_oo,2.*!pi/q_p_1gridarr(*,0,i),delxarr(*,0,i),psym=sym(6),symsize=1;,yran=[1e-3,10],/ys;,xran=[10,1000],/xs;,yran=[10d-3,10^2],/ys,symsize=.35
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


if i EQ 4 then legend, ['R'], /right
if i EQ 5 then legend, ['I'],/right
if i EQ 6 then legend, ['Z'],/right
if i EQ 7 then legend, ['Y'], /right
if i EQ 8 then legend, ['YHS'],/right
if i EQ 9 then legend, ['J'],/right
if i EQ 10 then legend, ['H'],/right
if i EQ 11 then legend, ['Ks'],/right
if i EQ 12 then legend, ['4.5'+mu+'m'],/right
if i EQ 13 then legend, ['2-10kev'],/right
if i EQ 14 then legend, ['0.5-2kev'],/right
if i EQ 15 then legend, ['3.6'+mu+'m'],/right
multiplot
endfor
device,/close
multiplot,/reset


set_plot,'PS'
cgerase
DEVICE, FILE='ch'+ch+'_'+field+'IR_autopower.ps', /COLOR, BITS=8,  XSIZE=25, ysize=25, xoffset=-5, yoffset=-1


!y.tickname=''

;plot_oo,2.*!pi/q_p_1gridarr(*,0),del1arr(*,0),psym=sym(1),xtitle='2!7p!5/q  arcsec',ytitle=TEXTOIDL('q^2P_{2,X}(q)/2\pi'),symsize=1
;errplot,2.*!pi/q_p_1gridarr(*,0),del1arr(*,0)-sig_del1arr(*,0),del1arr(*,0)+sig_del1arr(*,0)

plot_oo,2.*!pi/q_p_1grid,delir,psym=sym(6),xtitle='2!7p!5/q  arcsec',ytitle=TEXTOIDL('q^2P_{2,X}(q)/2\pi'),symsize=1,yran=[1e-5,10],/ys, xran=[1,10000], /xs
errplot,2.*!pi/q_p_1grid,delir-sig_delir,delir+sig_delir

oplot,2.*!pi/q_p_1gridarr,delab,psym=sym(2),symsize=1
errplot,2.*!pi/q_p_1gridarr,delir-sig_delir,delir+sig_delir




device,/close


coherence_err_pos=dblarr(12,4,16)
for i=0,15 do begin
   for k=0,3 do begin
coherence_err_pos(*,k,i)=sqrt(((2*powerxarr(*,k,i))/(powerir(*)*poweroparr(*,k,i)))^2*(sig_pxarr(*,k,i))^2+(powerxarr(*,k,i)^2/(-1*powerir(*)^2*poweroparr(*,k,i)))^2*(sig_pir(*))^2+(powerxarr(*,k,i)^2/(-1*powerir(*)*poweroparr(*,k,i)^2))^2*(sig_poparr(*,k,i))^2)
endfor
endfor




set_plot,'PS'
cgerase
DEVICE, FILE='ch'+ch+'_'+field+'_coherence.ps', /COLOR, BITS=8,  XSIZE=25, ysize=25, xoffset=-2, yoffset=-1,/BOLD
multiplot, /reset
!y.tickname=''
!x.tickname=''
!p.THICK=3
print, coherencearr
multiplot,[3,4], /square, mytitle='Coherence', mxtitle='2!7p!5/q  arcsec',mtitsize=1.4, mytitsize=1.4,mxtitsize=1.4
print, 'THIS IS IMPORTANT',2.*!pi/q_p_1gridarr(*,0,1)
for i=0,15 do begin
   if i EQ 5 or i EQ 6 or i EQ 7 or i EQ 8 then continue

plot_oo,2.*!pi/q_p_1gridarr(*,0,i),coherencearr(*,0,i),psym=sym(10),symsize=1,yran=[10e-4,1],/ys;,xran=[10,1000],/xs;,symsize=.35
errplot,2.*!pi/q_p_1gridarr(*,0,i),coherencearr(*,0,i)-coherence_err_pos(*,0,i),coherencearr(*,0,i)+coherence_err_pos(*,0,i)

oplot,2.*!pi/q_p_1gridarr(*,1,i),coherencearr(*,1,i),psym=sym(2),symsize=1
errplot,2.*!pi/q_p_1gridarr(*,1,i),coherencearr(*,1,i)-coherence_err_pos(*,1,i),coherencearr(*,1,i)+coherence_err_pos(*,1,i)

oplot,2.*!pi/q_p_1gridarr(*,2,i),coherencearr(*,2,i),psym=sym(6),symsize=1
errplot,2.*!pi/q_p_1gridarr(*,2,i),coherencearr(*,2,i)-coherence_err_pos(*,2,i),coherencearr(*,2,i)+coherence_err_pos(*,2,i)

oplot,2.*!pi/q_p_1gridarr(*,3,i),coherencearr(*,3,i),psym=sym(4),symsize=1
errplot,2.*!pi/q_p_1gridarr(*,3,i),coherencearr(*,3,i)-coherence_err_pos(*,3,i),coherencearr(*,3,i)+coherence_err_pos(*,3,i)



if i EQ 0 then legend, ['G'], /left
if i EQ 1 then legend, ['Y'],/left
if i EQ 2 then legend, ['I'],/left
if i EQ 3 then legend, ['Z'], /left
if i EQ 4 then legend, ['R'], /left
if i EQ 5 then legend, ['I'],/left
if i EQ 6 then legend, ['Z'],/left
if i EQ 7 then legend, ['Y'], /left
if i EQ 8 then legend, ['YHS'],/left
if i EQ 9 then legend, ['J'],/left
if i EQ 10 then legend, ['H'],/left
if i EQ 11 then legend, ['Ks'],/left
if i EQ 12 then legend, ['4.5'+mu+'m'],/left
if i EQ 13 then legend, ['2-10kev'],/left
if i EQ 14 then legend, ['0.5-2kev'],/left
if i EQ 15 then legend, ['3.6'+mu+'m'],/left
multiplot
endfor
device,/close
multiplot,/reset
write_csv,'coherrg',coherence_err_pos(*,0,0)
write_csv,'coherry',coherence_err_pos(*,0,1)


;;create a plot of cohernece v wavelength @500 arcsec
wavelength=[2.06,9.91,4750,6288,7683,7700,8875,9105,9791,10150,10214,12534,16453,21539,35634,45110]
coherence_500=dblarr(16,4)
coherence_500(0,*)=coherencearr(9,*,0)
coherence_500(1,*)=coherencearr(9,*,4)
coherence_500(2,*)=coherencearr(9,*,5)
coherence_500(3,*)=coherencearr(9,*,2)

coherence_500(4,*)=coherencearr(9,*,3)
coherence_500(5,*)=coherencearr(9,*,6)
coherence_500(6,*)=coherencearr(9,*,8)
coherence_500(7,*)=coherencearr(9,*,1)
coherence_500(8,*)=coherencearr(9,*,7)
coherence_500(9,*)=coherencearr(9,*,9)
coherence_500(10,*)=coherencearr(9,*,10)
coherence_500(11,*)=coherencearr(9,*,11)
coherence_500(12,*)=coherencearr(9,*,12)
coherence_500(13,*)=coherencearr(9,*,14)
coherence_500(14,*)=coherencearr(9,*,13)
coherence_500(15,*)=coherencearr(9,*,15)

coherence_err_500=dblarr(16,4)

coherence_err_500(0,*)=coherence_err_pos(9,*,0)
coherence_err_500(1,*)=coherence_err_pos(9,*,4)
coherence_err_500(2,*)=coherence_err_pos(9,*,5)
coherence_err_500(3,*)=coherence_err_pos(9,*,2)

coherence_err_500(4,*)=coherence_err_pos(9,*,3)
coherence_err_500(5,*)=coherence_err_pos(9,*,6)
coherence_err_500(6,*)=coherence_err_pos(9,*,8)
coherence_err_500(7,*)=coherence_err_pos(9,*,1)
coherence_err_500(8,*)=coherence_err_pos(9,*,7)
coherence_err_500(9,*)=coherence_err_pos(9,*,9)
coherence_err_500(10,*)=coherence_err_pos(9,*,10)
coherence_err_500(11,*)=coherence_err_pos(9,*,11)
coherence_err_500(12,*)=coherence_err_pos(9,*,12)
coherence_err_500(13,*)=coherence_err_pos(9,*,14)
coherence_err_500(14,*)=coherence_err_pos(9,*,13)
coherence_err_500(15,*)=coherence_err_pos(9,*,15)

set_plot,'PS'
cgerase
DEVICE, FILE='ch'+ch+'_'+field+'_coherence_v_wavelength_500.ps', /COLOR, BITS=8, XSIZE=22, ysize=22, xoffset=-1, yoffset=0,/BOLD

multiplot, /reset
!y.tickname=''
!P.thick=3
multiplot,[2,2], /square, mytitle='Coherence',mxtitle='Angstrom',mtitsize=1.4, mytitsize=1.4, mxtitsize=1.4;, /doyaxis
for i=0,3 do begin

plot_oo,wavelength,coherence_500(*,i),psym=sym(5),symsize=1,yran=[10e-3,1],/ys;,xran=[10,1000],/xs;,symsize=.35
errplot,wavelength,coherence_500(*,i)-coherence_err_500(*,i),coherence_500(*,i)+coherence_err_500(*,i)

if i EQ 0 then legend, ['m0'], /left
if i EQ 1 then legend, ['m0+2'],/left
if i EQ 2 then legend, ['m0+4'],/left
if i EQ 3 then legend, ['m0+6'], /left
;do inset plot for the xray sources                                                                      
plot_oo,wavelength,coherence_500(*,0),psym=sym(5),symsize=0.5,xran=[1,15],/xs,yran=[10e-3,1],/ys,POSITION=[.43,.77,.53,.87]
errplot,wavelength,coherence_500(*,0)-coherence_err_500(*,0),coherence_500(*,0)+coherence_err_500(*,0)

plot_oo,wavelength,coherence_500(*,1),psym=sym(5),symsize=0.5,xran=[1,15],/xs,yran=[10e-3,1],/ys,POSITION=[.85,.77,.95,.87]
errplot,wavelength,coherence_500(*,1)-coherence_err_500(*,1),coherence_500(*,1)+coherence_err_500(*,1)


plot_oo,wavelength,coherence_500(*,2),psym=sym(5),symsize=0.5,xran=[1,15],/xs,yran=[10e-3,1],/ys,POSITION=[.43,.38,.53,.48]
errplot,wavelength,coherence_500(*,2)-coherence_err_500(*,2),coherence_500(*,2)+coherence_err_500(*,2)

plot_oo,wavelength,coherence_500(*,3),psym=sym(5),symsize=0.5,xran=[1,15],/xs,yran=[10e-3,1],/ys,POSITION=[.85,.38,.95,.48]
errplot,wavelength,coherence_500(*,3)-coherence_err_500(*,3),coherence_500(*,3)+coherence_err_500(*,3)
multiplot
endfor
device,/close
multiplot,/reset

;;create a plot of wavelength v coherence @1000 arcsec
wavelength=[2.06,9.91,4750,6288,7683,7700,8875,9105,9791,10150,10214,12534,16453,21539,35634,45110]
coherence_1000=dblarr(16,4)
coherence_1000(0,*)=coherencearr(10,*,0)
coherence_1000(1,*)=coherencearr(10,*,4)
coherence_1000(2,*)=coherencearr(10,*,5)
coherence_1000(3,*)=coherencearr(10,*,2)

coherence_1000(4,*)=coherencearr(10,*,3)
coherence_1000(5,*)=coherencearr(10,*,6)
coherence_1000(6,*)=coherencearr(10,*,8)
coherence_1000(7,*)=coherencearr(10,*,1)
coherence_1000(8,*)=coherencearr(10,*,7)
coherence_1000(9,*)=coherencearr(10,*,9)
coherence_1000(10,*)=coherencearr(10,*,10)
coherence_1000(11,*)=coherencearr(10,*,11)
coherence_1000(12,*)=coherencearr(10,*,12)
coherence_1000(13,*)=coherencearr(9,*,14)
coherence_1000(14,*)=coherencearr(9,*,13)
coherence_1000(15,*)=coherencearr(9,*,15)



coherence_err_1000=dblarr(16,4)

coherence_err_1000(0,*)=coherence_err_pos(10,*,0)
coherence_err_1000(1,*)=coherence_err_pos(10,*,4)
coherence_err_1000(2,*)=coherence_err_pos(10,*,5)
coherence_err_1000(3,*)=coherence_err_pos(10,*,2)

coherence_err_1000(4,*)=coherence_err_pos(10,*,3)
coherence_err_1000(5,*)=coherence_err_pos(10,*,6)
coherence_err_1000(6,*)=coherence_err_pos(10,*,8)
coherence_err_1000(7,*)=coherence_err_pos(10,*,1)
coherence_err_1000(8,*)=coherence_err_pos(10,*,7)
coherence_err_1000(9,*)=coherence_err_pos(10,*,9)
coherence_err_1000(10,*)=coherence_err_pos(10,*,10)
coherence_err_1000(11,*)=coherence_err_pos(10,*,11)
coherence_err_1000(12,*)=coherence_err_pos(10,*,12)

coherence_err_1000(13,*)=coherence_err_pos(9,*,14)
coherence_err_1000(14,*)=coherence_err_pos(9,*,13)
coherence_err_1000(15,*)=coherence_err_pos(9,*,15)







set_plot,'PS'
cgerase
DEVICE, FILE='ch'+ch+'_'+field+'_coherence_v_wavelength_1000.ps', /COLOR, BITS=8,  XSIZE=22, ysize=22, xoffset=-1, yoffset=0, /BOLD

multiplot, /reset
!y.tickname=''
!P.thick=3
multiplot,[2,2], /square, mytitle='Coherence', mxtitle='Angstrom',mtitsize=1.4, mytitsize=1.4,mxtitsize=1.4
for i=0,3 do begin

plot_oo,wavelength,coherence_1000(*,i),psym=sym(5),symsize=1,yran=[10e-3,1],/ys,xran=[1800,100000],/xs
errplot,wavelength,coherence_1000(*,i)-coherence_err_1000(*,i),coherence_1000(*,i)+coherence_err_1000(*,i)

if i EQ 0 then legend, ['m0'], /left
if i EQ 1 then legend, ['m0+2'],/left
if i EQ 2 then legend, ['m0+4'],/left
if i EQ 3 then legend, ['m0+6'], /left
;do inset plot for the xray sources
plot_oo,wavelength,coherence_1000(*,0),psym=sym(5),symsize=0.5,xran=[1,15],/xs,yran=[10e-3,1],/ys,POSITION=[.43,.77,.53,.87]
errplot,wavelength,coherence_1000(*,0)-coherence_err_1000(*,0),coherence_1000(*,0)+coherence_err_1000(*,0)

plot_oo,wavelength,coherence_1000(*,1),psym=sym(5),symsize=0.5,xran=[1,15],/xs,yran=[10e-3,1],/ys,POSITION=[.85,.77,.95,.87]
errplot,wavelength,coherence_1000(*,1)-coherence_err_1000(*,1),coherence_1000(*,1)+coherence_err_1000(*,1)


plot_oo,wavelength,coherence_1000(*,2),psym=sym(5),symsize=0.5,xran=[1,15],/xs,yran=[10e-3,1],/ys,POSITION=[.43,.38,.53,.48]
errplot,wavelength,coherence_1000(*,2)-coherence_err_1000(*,2),coherence_1000(*,2)+coherence_err_1000(*,2)

plot_oo,wavelength,coherence_1000(*,3),psym=sym(5),symsize=0.5,xran=[1,15],/xs,yran=[10e-3,1],/ys,POSITION=[.85,.38,.95,.48]
errplot,wavelength,coherence_1000(*,3)-coherence_err_1000(*,3),coherence_1000(*,3)+coherence_err_1000(*,3)


multiplot
endfor
device,/close
multiplot,/reset

endfor
endfor

end
