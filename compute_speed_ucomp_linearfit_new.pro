
FUNCTION compute_speed_ucomp_linearfit_new, vel, delx, delt, debug=debug, ret=ret, force_zero=force_zero


s=size(vel)
nt=s(1)
npt=s(2)
IF npt mod 2 EQ 0 THEN npt=npt-1 ;deals with even length tracks
icent=npt/2



nlag=5 ;number of points in cross correlation (make odd)
lag=findgen(nlag)-fix(nlag/2.)
dist_list0=findgen(npt)-fix(npt/2.)


n = 25 < npt  ;number of points along track for initial guess (make odd)

dist_list0=dist_list0[(npt/2-n/2):(npt/2+n/2)]

lag_list=fltarr(n)
ccor_list=fltarr(n)

FOR i=0,n-1 DO BEGIN
	ccor=c_correlate(vel(*,icent),vel(*,i),lag)
	; ccor=c_correlate(interpol(vel(tind1:tind2,i),(tind2-tind1+1)*10),interpol(vel(tind1:tind2,i+1),(tind2-tind1+1)*10),lag)
	m=max(ccor,imax)
	imax=1 > imax < (nlag-2)
	lag1=parabola(lag(imax-1:imax+1),1.-ccor(imax-1:imax+1))
    lag_list[i]=lag1
    ccor_list[i]=interpolate(ccor,lag1-lag[0])
ENDFOR


  flag=where(abs(lag_list) lt 5.0*median(abs(lag_list)))
  dist_list=dist_list0[flag]
  lag_list=lag_list[flag]
  w_list=(abs(ccor_list[flag])-0.5)^2



s=SPEARMAN(dist_list,lag_list)
number=0

while abs(s) lt 0.95 and n_elements(dist_list) gt 6 and number lt 100 do begin
    s=FAST_SPEARMAN(dist_list,lag_list,s_list,w=w_list)
    s_max=max(abs(s_list),point1)
    dist_list=dist_list[where(abs(s_list) ne s_max)]
    lag_list=lag_list[where(abs(s_list) ne s_max)]
    w_list=w_list[where(abs(s_list) ne s_max)]
	  number+=1
end

param1=(total(w_list*lag_list*dist_list)-total(w_list*lag_list)*total(w_list*dist_list)/total(w_list))/$
    (total(w_list*dist_list^2.)-total(w_list*dist_list)^2./total(w_list))*delt/delx

phase_speed=1./param1

err_param1=sqrt((total(w_list))/(total(w_list*(dist_list^2.))*total(w_list)-(total(w_list*dist_list))^2.))

speed_error=err_param1/param1^2.

return,[phase_speed,speed_error]

end

FUNCTION SPEARMAN,X,Y
  ;IF N_ELEMENTS(X) EQ 1 OR N_ELEMENTS(Y) EQ 1 THEN MESSAGE,'SPEARMAN: X OR Y SHOULD BE ARRAY.'
  ;IF N_ELEMENTS(X) NE N_ELEMENTS(Y) THEN MESSAGE,'SPEARMAN: X AND Y HAVE DIFFERNET SIZES!'
  
  N=N_ELEMENTS(X)
  
  TX=TOTAL(X)
  TY=TOTAL(Y)
  
  RETURN,(TOTAL(X*Y)-TX*TY/N)/SQRT((TOTAL(X*X)-TX*TX/N)*(TOTAL(Y*Y)-TY*TY/N))
  
  ;RETURN,total((x-total(x)/N)*(y-total(y)/N))/sqrt(total((x-total(x)/N)^2)*total((y-total(y)/N)^2))

END

FUNCTION FAST_SPEARMAN,X,Y,S_LIST,W=W

  TW=TOTAL(W)

  TX=TOTAL(X*W)
  TY=TOTAL(Y*W)
  
  M1=X*Y*W
  M2=X*X*W
  M3=Y*Y*W
  
  A1=TOTAL(M1)
  A2=TOTAL(M2)
  A3=TOTAL(M3)
  B1=TX*TY
  B2=TX*TX
  B3=TY*TY
  
  A1_LIST=A1-M1
  A2_LIST=A2-M2
  A3_LIST=A3-M3
  
  B1_LIST=B1-X*TY*W-Y*TX*W+M1*W
  B2_LIST=B2-2*X*TX*W+M2*W
  B3_LIST=B3-2*Y*TY*W+M3*W
  
  S_LIST=(A1_LIST-B1_LIST/(TW-W))/SQRT((A2_LIST-B2_LIST/(TW-W))*(A3_LIST-B3_LIST/(TW-W)))
  
  S=(A1-B1/TW)/SQRT((A2-B2/TW)*(A3-B3/TW))
  
  return,S

END

