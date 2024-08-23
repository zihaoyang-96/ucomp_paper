;pro prewavetracking_ucomp
set_plot,'x'

Device, RETAIN=2 


date='20220805';'20120409';'20111230';
cd,'/System/Volumes/Data/Data/UCoMP/verified/'+date+'/level2/Dynamics/wave/'
; cd,'/Volumes/Elements/'+date+'/level2/Dynamics/wave/'

wavelen='1074'

restore,'./'+wavelen+'_cube.sav' ; for level 2 data

dt=33.5   ;mean cadence of data
sigma=10.   ;tension of spline
s=size(fit_arr)  &  nx=s(1)  &  ny=s(2)  &  ntime=s(4)

time_sec=anytim(time_obs) ;in the unit of sec
old_time=time_sec-time_sec[0]
time=old_time(0)+findgen(ntime)*dt;change the cadence into uniform 30s

; goto,findangle

for i=0,nx-1 do begin
for j=0,ny-1 do begin
  fit_arr(i,j,1,*)=spline(old_time,reform(fit_arr(i,j,1,*)),time,sigma,/double)
endfor
endfor
velo=reform(fit_arr[*,*,1,*])
fit_arr=0

save,velo,time,file='./velocity_interp.sav'



filtervelo:
restore,'./velocity_interp.sav'
t1=0 & t2=t1+ntime-1
ntap=30     ;number of points to taper each end of time series
freq=findgen(ntime)/(float(ntime)*dt)
for i=0,ntime/2-1 do freq(ntime-i-1)=freq(i+1)
freq0=.0035  ;Hz, ~5 min
fwidth=.0015 ;Hz
filter=exp( -(freq-freq0)^2/fwidth^2 )
filter(0)=0.
for i=0,nx-1 do begin
;print,i
for j=0,ny-1 do begin
  y=reform(velo(i,j,*))
  y=y-median(y)
  ; taper,y,ntap   ;taper points at start and end of timeseries
  trans=filter*fft(y)
  velo(i,j,0:ntime-1)=fft(trans,/inverse)
endfor
endfor
time=time(t1:t2)
save,velo,time,file='./velocity_filtered.sav'



velofft:
; ; compute spectra of velocity data cube for UCoMP waves analysis
restore,'./velocity_interp.sav'
nspec=ntime/2
spec=complexarr(nx,ny,nspec,/nozero)
h=do_apod(ntime,cpg)
for ix=0,nx-1 do for iy=0,ny-1 do begin
    d=reform(velo[ix,iy,*]-mean(velo[ix,iy,*])) 
    ; d=reform(velo[ix,iy,*])
    sp=fft(d*h)
    ; sp=fft(d)
    spec(ix,iy,*)=sp(0:nspec-1)
endfor
save,spec,file='./velocity_unfiltered_spec.sav' ;velocity_filtered_spec.sav


; goto,findangle


velofftsmooth:
; ; compute spectra of velocity data cube for CoMP waves analysis
restore,'./velocity_interp.sav'
for i=0, ntime-1 do begin
    velo[*,*,i]=smooth(velo[*,*,i],3,/nan)
end
nspec=ntime/2
spec=complexarr(nx,ny,nspec,/nozero)
h=do_apod(ntime,cpg)
for ix=0,nx-1 do for iy=0,ny-1 do begin
    d=reform(velo[ix,iy,*]-mean(velo[ix,iy,*]))  ;-smooth(velo(ix,iy,*),40,/edge_truncate)
    ; d=reform(velo[ix,iy,*])
    sp=fft(d*h)
    ; sp=fft(d)
    spec(ix,iy,*)=sp(0:nspec-1)
endfor
save,spec,velo,file='./velocity_unfiltered_spec_smoothed.sav' ;velocity_filtered_spec.sav

; endprogram:
; end


; find wave propagation angle
findangle:
xyrange=[0,1279,0,1023] ;xmin,xmax,ymin,ymax; spatial range of calculation
; xyrange=[150,350,400,700]
dt=33.5
restore,'./velocity_unfiltered_spec_smoothed.sav'
find_angle_ucomp_new_err,spec,dt,xyrange,wave_angle,angle_error,debug='no',date=date
save,wave_angle,angle_error,file='./wave_angle_3.5mHz_smoothed_new_2024.sav'



endprogram:
end

