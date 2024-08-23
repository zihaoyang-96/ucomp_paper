;pro space_time_run_ucomp_new,delt=delt


;  Procedure to compute space-time diagram of velocity time series over a whole image.
;  The track is computed from the angle map. This routine calls compute_speed to
;  compute the phase speed.


set_plot,'x'

debug='no'
date='20221026';'20120409';'20111230';
; cd,'/System/Volumes/Data/Data/UCoMP/verified/'+date+'/level2/Dynamics/wave/'
cd,'/Volumes/Elements/'+date+'/level2/Dynamics/wave/'



ans=' '
torad=!pi/180.

restore,'./velocity_filtered.sav'



ntime=n_elements(time)

restore,'./wave_angle_3.5mHz_smoothed_new.sav'
angle=wave_angle


s=size(velo)
print,s
nx=s(1)          ;x dimension of velocity data cube
ny=s(2)          ;y dimension
ntime=s(3)     ;time dimension

if ((ntime mod 2) eq 0) then nt = ntime else nt = ntime-1    ;number of time points to use (make even)
max_npt=31
npt=max_npt    ;number of steps to map along track (make odd)

dx=1.0     ;step size along track


; fff=file_search('/System/Volumes/Data/Data/UCoMP/verified/'+date+'/level2/Dynamics/wave/*ucomp*.fts')
fff=file_search('/Volumes/Elements/'+date+'/level2/Dynamics/wave/*.ucomp.*.dynamics*.fts')
f=fff[1]
fits_read,f,data,header,exten_no=0
index = fitshead2struct(header)
xscale=index.cdelt1*0.725*dx

xtrack=fltarr(npt)   ;x and y position along track
ytrack=fltarr(npt)
ang_track=fltarr(npt)   ;angle along track


delt=33.5

f=findgen(nt)/(float(nt)*delt)               ;compute temporal frequency scale
for i=0,nt/2 do f(nt-i-1)=f(i+1)
freq=rebin(f,nt,npt)

w=0.001      ;create filter to remove low frequency signal
filter=1.-exp(-freq^2/w^2)

angle=median(angle, 3) ;median smoothing of wave angle to reduce noise (from R. Morton's code)


nx=1280 & ny=1024


r_inner = index.radius1+4;225
r_upper=600.
mask = intarr(nx, ny)
promask = bytarr(nx, ny)
rpos=fltarr(nx,ny)
for xx=0,nx-1 do begin
  for yy=0,ny-1 do begin
    rpos[xx,yy] = sqrt((xx-640.)^2 + (yy-512.)^2)
    if (rpos[xx,yy] ge r_inner) and (rpos[xx,yy] le r_upper) and angle[xx,yy] ne 0 then mask[xx,yy] = 1
  endfor
endfor

rad=rpos

nmask=intarr(nx,ny,nt)
for tt=0,nt-1 do begin
  nmask[*,*,tt]=mask
endfor

angle=angle*mask
bad_velo=where(finite(velo,/nan))
velo[bad_velo]=0
velo=velo*nmask

; stop

; display first velocity image for interactive input
if debug eq 'yes' then  begin
    window,0,xs=nx,ys=ny,retain=2,xpos=0,ypos=700
    window,1,xs=nx,ys=ny,retain=2,xpos=nx+10,ypos=700
    window,2,xs=nt*2,ys=npt*2,title='Space-Time Diagram'
    window,3,xs=nt*2,ys=npt*2,xpos=0,ypos=2*npt+45,title='Prograde Filtered'
    window,4,xs=nt*2,ys=npt*2,xpos=0,ypos=4*npt+90,title='Retrograde Filtered'
    window,8,xs=600,ys=800
    device,decomposed=0

    wset,0    ;display images
    img1=velo(*,*,0)
    tvscl,img1<9>(-9)
    wset,1
    img2=bytscl(angle,-90,90)
    tvscl,img2  
end


pro_speed=fltarr(nx,ny)
pro_speed_err=fltarr(nx,ny)

;  x and y limits for map

xstart=0 
xend=nx-1
ystart=0
yend=ny-1


pixcounter = 0.
old_perc_proc = 0
tot_pix = (xend-xstart)*1L*(yend-ystart)


pro_speed=fltarr(nx,ny)
pro_speed_err=fltarr(nx,ny)
pro_speed_err_u=fltarr(nx,ny)
pro_speed_err_l=fltarr(nx,ny)
pro_power=fltarr(nx,ny)
ret_power=fltarr(nx,ny)

old_pro_speed=fltarr(nx,ny)
old_pro_speed_err=fltarr(nx,ny)


for ix=xstart,xend do for iy=ystart,yend do if angle[ix,iy] ne 0. then begin

  img1=velo(*,*,0)   ;reinitialize velocity image

  npt=31              ;redefine npt and x/ytrack due to possible change of npt when clost to the boundary.(note: it is not always good to shorten the track, npt=31 is a good option for most times)
  xtrack=fltarr(npt)   ;x and y position along track
  ytrack=fltarr(npt)   ;............................

  imid=fix(npt/2)

  xtrack(imid)=float(ix) & ytrack(imid)=float(iy)
  ang_track(imid)=angle(ix,iy)


;  use angle at each point to define track, starting at cursor position as central point
                                                                                     
if ix le nx/2-1 then leftright=1 else leftright=-1

ntrack=0
if sqrt((float(ix)-nx/2.-0.5)^2.+(float(iy)-ny/2.-0.5)^2.)-imid lt r_inner then short_track=1 else short_track=-1  ;condition for position close to the boundary

;  first, move out from cursor position

  for i=imid+1,npt-1 do begin
    xtrack(i)=0. > xtrack(i-1)-leftright*dx*cos( ang_track(i-1)*torad) < float(nx-1)
    ytrack(i)=0. > ytrack(i-1)-leftright*dx*sin( ang_track(i-1)*torad) < float(ny-1)

    intx=0 > fix(xtrack(i)+0.5) < (nx-1)
    inty=0 > fix(ytrack(i)+0.5) < (ny-1)

    pixcheck=sqrt((intx-nx/2.+0.5)^2 + (inty-ny/2.+0.5)^2)
		     IF pixcheck lt r_inner THEN BEGIN
	       	        ;this will not allow to track further down in the inward direction!
		        xtrack(i)=0.
	 	        ytrack(i)=0.

		        BREAK
	                
		     ENDIF
	   
		     IF pixcheck GT r_upper THEN BEGIN
	 	        xtrack(i)=0.
		        ytrack(i)=0.

		        BREAK
		     ENDIF
    ntrack=ntrack+1
    if debug eq 'yes' then begin
    img1(intx,inty)=min(img1) ;put track into displayed images
    endif
    ang_track(i)=angle(intx,inty)  ;  nearest neighbor interp for angle (works better than bilinear at boundaries)
;  if angle is zero (not defined), use previous value

    if ang_track(i) eq 0. then ang_track(i)=ang_track(i-1)

;  if big difference in adjacent angles, resolve 180 degree ambiguity by choosing angle closest to ang_track(i-1)

    if abs(ang_track(i)-ang_track(i-1)) gt 90. then begin
      if ang_track(i)-ang_track(i-1) gt 0. then ang_track(i)=ang_track(i)-180. else ang_track(i)=ang_track(i)+180.
    endif

  endfor

;  next, move in from cursor position

  for i=imid-1,0,-1 do begin
    ; factor=(ang_track(i+1)+180.)*torad
    xtrack(i)=0. > xtrack(i+1)-leftright*dx*cos( (ang_track(i+1)+180.)*torad) < float(nx-1)
    ytrack(i)=0. > ytrack(i+1)-leftright*dx*sin( (ang_track(i+1)+180.)*torad) < float(ny-1)

    intx=0 > fix(xtrack(i)+0.5) < (nx-1)
    inty=0 > fix(ytrack(i)+0.5) < (ny-1)

	  pixcheck=sqrt((intx-nx/2.+0.5)^2 + (inty-ny/2.+0.5)^2)
		    IF pixcheck lt r_inner THEN BEGIN
		       ;this will not allow to track further down in the inward direction!
		      xtrack(i)=0.
		      ytrack(i)=0.
			  BREAK ; stops tracks when hits lower boundary
	        ENDIF
	   
		    IF pixcheck GT r_upper THEN BEGIN
		        xtrack(i)=0.
		        ytrack(i)=0.
			    BREAK ; stops tracks when hits upper boundary
    		ENDIF
    ; print,'r=',sqrt((intx-nx/2.-0.5)^2.+(inty-ny/2.-0.5)^2.)
    ntrack=ntrack+1

    if debug eq 'yes' then begin
    img1(intx,inty)=min(img1) ;put track into displayed images
    endif


    ang_track(i)=angle(intx,inty)  ;  nearest neighbor interp for angle (works better than bilinear at boundaries)

;  if angle is zero (not defined), use previous value

    if ang_track(i) eq 0. then ang_track(i)=ang_track(i+1)

;  if big difference in adjacent angles, resolve 180 degree ambiguity by choosing angle closest to ang_track(i-1)

    if abs(ang_track(i)-ang_track(i+1)) gt 90. then begin
      if ang_track(i)-ang_track(i+1) gt 0. then ang_track(i)=ang_track(i)-180. else ang_track(i)=ang_track(i)+180.
    endif

  endfor

  if ntrack lt 4 then continue
  if ntrack mod 2 eq 0 then ntrack=ntrack+1   ;make it odd

  if ntrack lt max_npt then begin
        in=where(xtrack ne 0)
        newxtrack=xtrack(in)
        xtrack=newxtrack
        newytrack=ytrack(in)
        ytrack=newytrack
        npt=n_elements(in)
  endif else begin
  npt=max_npt
  endelse
  vmap=fltarr(nt,npt)

;  interpolate velo onto track

  for i=0,nt-1 do begin
    ; vel=velo(*,*,i) ;modified by Xianyu Liu, skip this step, do 3D interpolation instead. This makes the calculation over 100x faster 
    v=interpolate(velo,xtrack,ytrack,replicate(i,npt))
    vmap(i,*)=v
  endfor

  if debug eq 'yes' then begin
    wset,0       ;redisplay velocity image with track
  tvscl,img1
  endif

vmap=vmap-mean(vmap)
vmap=vmap-transpose(rebin(mean(vmap,dim=1),npt,nt),[1,0])


;filter velocity map

  for i=0,npt-1 do vmap(*,i)=vmap(*,i)-mean(vmap(*,i)) ;remove temporal mean


  trans=fft(vmap,-1)       ;compute fourier transform
  pro_trans=trans             ;select prograde waves (assume nt even, npt odd)
  pro_trans(1:nt/2-1,0:npt/2)=0.
  pro_trans(nt/2+1:nt-1,npt/2+1:npt-1)=0.

  filtrans=pro_trans


;filtering the low-speed components from the k-omega diagram
  filtvel=0.2

  for xx=nt/2,nt-1 do begin
    for yy=0,npt/2 do begin
      if yy-(xscale/delt)*(1./filtvel)*nt+(xscale/delt)*(1./filtvel)*xx ge 0 then filtrans[xx,yy]=0
    end
  end

  for xx=0,nt/2 do begin
    for yy=npt/2+1,npt-1 do begin
      if xx+(delt/xscale)*filtvel*yy-(delt/xscale)*filtvel*npt le 0 then filtrans[xx,yy]=0
    end
  end

  pro_vel=real_part(fft(filtrans,1))

  r=compute_speed_ucomp_linearfit_new(pro_vel,xscale,delt)
  

  pro_speed(ix,iy)=abs(r(0))
  pro_speed_err(ix,iy)=abs(r(1))




  !p.multi=0
endif



save,file='./'+date+'_phase_speed_new.sav',pro_speed,pro_speed_err

end
