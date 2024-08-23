
pro find_angle_ucomp_new,spec,dt,xyrange,wave_angle,angle_error,angle_err,debug=debug,date=date



ans=' '
torad=!pi/180.
; 'yes';  ;debug mode, yes or no

s=size(spec)
if debug eq 'yes' then print,s
nx=s(1)
ny=s(2)
nspec=s(3)

;  create filter
;  gaussian filter peaked at freq0 with width of fwidth
freq=findgen(nspec)/(float(nspec*2)*dt) ;note that nspec=ntime/2, see prewavetracking.pro
nfreq=n_elements(freq)
freq0=.0035
fwidth=.0015
filter=exp( -(freq-freq0)^2/fwidth^2 )
filter(0)=0.    ;set dc to zero
filter(where (freq lt .001))=0.     ;set low frequencies to zero
filter=filter/total(filter)

if debug eq 'yes' then begin
  window
  plot,freq,filter
;   read,'enter return',ans
endif

dx=20
nbox=2*dx+1   ;box size in pixels

nsmooth= 11;21 ;51  ;smoothing width for cross spectra (nominally 11)
; nsmooth=3
limit= 0.5;0.6; 0.9;    ;set coherence limit for acceptance


f=file_search('/System/Volumes/Data/Data/UCoMP/verified/'+date+'/level2/Dynamics/wave/*.ucomp.1074.dynamics.fts')
; f=file_search('/Volumes/Elements/'+date+'/level2/Dynamics/wave/*.ucomp.*.dynamics*.fts')

fits_read,f[0],data,header,exten_no=0
index = fitshead2struct(header)

nx=1280 & ny=1024

rinner = index.oradius1+3.;225
mask = intarr(nx, ny)
promask = bytarr(nx, ny)
for xx=0,nx-1 do begin
  for yy=0,ny-1 do begin
    rpos = sqrt((xx-640.)^2 + (yy-512.)^2)
    if (rpos ge rinner) then mask[xx,yy] = 1
  endfor
endfor
; stop
xstart=xyrange[0]
xend=xyrange[1]
ystart=xyrange[2]
yend=xyrange[3]

; open windows
if debug eq 'yes' then begin
  f=8
window,0,xs=nx,ys=ny,retain=2,xpos=400,ypos=0
window,1,xs=nbox*f,ys=nbox*f,xpos=0,ypos=0
window,2,xs=nbox*f,ys=nbox*f,xpos=0,ypos=f*nbox+50
window,3,xs=nbox*f,ys=nbox*f,xpos=0,ypos=2*f*nbox+100
window,4,xs=nx,ys=ny,retain=2,xpos=400+nx+20,ypos=0
window,5,xs=nx,ys=ny,retain=2,xpos=400,ypos=ny+50
window,6,xs=512,ys=512,xpos=400,ypos=800
device,decomposed=0
endif

if debug eq 'no' then begin
  pixcounter = 0L
  old_perc_proc=0
  tot_pix=long(n_elements(where (mask eq 1)))
endif

conjspec=conj(spec)

lar_g=fltarr(nx,ny,nspec)
FOR iy=0,ny-1 DO FOR ix=0,nx-1 DO $
    IF (mask[ix,iy] EQ 1) THEN lar_g[ix,iy,*]=smooth(real_part(reform(spec[ix,iy,*]*conjspec[ix,iy,*])),nsmooth)

xbox=rebin(findgen(nbox)-dx,nbox,nbox)  ;coordinates of points in box
ybox=transpose(xbox)

wave_angle=fltarr(nx,ny)
angle_error=fltarr(nx,ny)
angle_err=fltarr(nx,ny)
coh_measure = fltarr(nx,ny,2)


;calculates rms value of spectrum
;Used to provide secondary condition on for loops
;and exclude pixels with low signal
rms_spec=total(abs(spec),3)

for ix=xstart,xend do for iy=ystart,yend do if (mask(ix,iy) eq 1) and (rms_spec[ix,iy] gt 1e-7) then begin

  if debug eq 'no' then begin
    perc_proc=round((float(pixcounter)/float(tot_pix))*100.)
    pixcounter=pixcounter+1
    IF perc_proc GT old_perc_proc THEN BEGIN
            IF strmid(getenv('TERM'),0,5) EQ 'xterm' THEN BEGIN
     		progress_comp, perc_proc, 100, msg='Computing wave propagation angles'
    	    ENDIF ELSE message, /cont, strcompress(string(perc_proc,format='(I3)'),/rem)+'% processed.'
    	    old_perc_proc = perc_proc
     ENDIF
  endif
  
;  compute coherence for contiguous pixels with coherence greater than limit

  coh=fltarr(nbox,nbox)           ;initialize coherence array to zero
  coh_mask=fltarr(nbox,nbox)      ;initialize mask of good points in coherence

  spec1=reform(spec[ix,iy,*])
  
  g1=reform(lar_g[ix,iy,*])

  coh(nbox/2,nbox/2)=1.
  coh_mask(nbox/2,nbox/2)=1.

  count=1

  n_coh_points=1
  
  while count gt 0 and n_coh_points lt 150 do begin


        current_points=where(coh_mask eq 1.)                                       ; currently found points.
        new_to_current_points=[current_points-1,current_points+1,$
          current_points+nbox,current_points-nbox]               ; use the pixels adjacent to current points as test points.
        new_mask=fltarr(nbox,nbox)
        new_mask(new_to_current_points)=1.
        new_mask=new_mask*(1.-coh_mask)
        new=where(new_mask gt 0.,new_count)
        count=0

        for point=0,n_elements(new)-1 do begin
          point_x=new(point) mod nbox        ; the x coordinate
          point_y=fix(new(point)/nbox)       ; the y coordinate

          ii=ix-nbox/2+point_x                ; the x position of this point in the whole map.
          jj=iy-nbox/2+point_y                ; the y position of this point in the whole map.

          if ii lt 0 or ii gt nx-1 or jj lt 0 or jj gt ny-1 then coh[point_x,point_y]=0. else begin

            if coh(point_x,point_y) eq 0. then begin
              cspec=smooth(spec1*reform(conjspec[ii,jj,*]),nsmooth)   ;cross-spectrum: FFT of reference pixel * conjugation of FFt of other pixel
              g2=lar_g[ii,jj,*]
              coh[point_x,point_y]=mask[ii,jj]*total( filter*( abs(cspec)/sqrt(g1*g2) ) ,/nan) ;see McInotosh et al. 2008. Sol. Phys.
            endif
            if coh[point_x,point_y] gt limit then coh_mask[point_x,point_y]=1.
            if coh[point_x,point_y] gt limit then count=count+1

          endelse

        endfor

        coh_points=where(coh_mask eq 1.,n_coh_points)

  endwhile

      weight=(coh[coh_points]-limit)^2
      theta=min_dist_fit_new(xbox(coh_points),ybox(coh_points),error = error,error_slope=error_slope,weight=weight) ;analytical solution to the minimum distance problem
      wave_angle[ix,iy] =theta
      angle_error[ix,iy]=error
      angle_err[ix,iy]=error_slope
  endif   




print,'done'

end

