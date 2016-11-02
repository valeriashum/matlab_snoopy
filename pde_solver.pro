;________________________
;_______Boundary_________
;_______Condition________
;________________________
PRO boundary, array, nx

  Compile_Opt DEFINT32

  array[0:5]=array[nx:nx+5]
  array[nx+6:nx+11]=array[6:11]

END 
;________________________
;_______derivatives______
;________________________
;________________________

PRO derivatives,t_derive, f_derive, B, R, nx,dx,x,o_derive,cos_flag=cos_flag  

  Compile_Opt DEFINT32

  if n_elements(cos_flag) eq 0 then cos_flag=1 

  f_derive=fltarr(nx+12)
  t_derive=fltarr(nx+12)
  o_derive=fltarr(nx+12)

;_ Fourth Derivative 

  f_derive[2:(nx+9)]=(b[0:(nx+7)]-4.*b[1:(nx+8)] +6.*b[2:(nx+9)] $
   -4.*b[3:(nx+10)] +b[4:(nx+11)])/dx^4
  boundary, f_derive, nx

;_ Second Derivative 
  if cos_flag eq 1 then begin
  t_derive[1:(nx+10)]=(cos(x[0:(nx+9)])*b[0:(nx+9)]-2.*cos(x[1:(nx+10)]) $
     *b[1:(nx+10)] +cos(x[2:(nx+11)])*b[2:(nx+11)])/dx^2


  endif else begin
  t_derive[1:(nx+10)]=(b[0:(nx+9)]-2.*b[1:(nx+10)] +b[2:(nx+11)])/dx^2

  endelse
  boundary, t_derive, nx
;_ First Derivative

  o_derive[1:(nx+10)]=(b[2:(nx+11)] - b[0:(nx+9)])/(2*dx)
  boundary, o_derive, nx

END


;________________________
;__________Initial_______
;_________Condition______
;________________________
PRO initial_condition, nx,dx, B, x

  x=dindgen(nx+12)*dx-11.*dx/2.
  b=fltarr(nx+12)
  b=cos(x)*0.01d0
boundary, b, nx

END 

;________________________
;_________Runge__________
;_________Kutta__________
;________________________

PRO renge_fed,length,nx,jump,T,R,b_out,time
 
  Compile_Opt DEFINT32

  dx=(length*2.*!pi)/double(nx) 			;grid step
  time=findgen(T+1)					;save time as an array 
 
  b_1=fltarr(nx+12)					;half step B 
  BB = fltarr(nx+12)					;store next step B
  B = fltarr(nx+12)	
  x = fltarr(nx+12)
  E = fltarr(nx+12)					;non-linear B^3
  E_1=fltarr(nx+12)					;half step B^3 
  b_out=fltarr(nx+12,T+1)				;store evolution of B as an array0
  meh = fltarr(nx+12)					;array to store cos*B^2 for V(t)
  slope = dblarr(T+1)					;growth/decay rates for B
  kk = fltarr(3)					;store choice of dt
  
  i=0							;loop variable
  ttt=0							;set initial current time 


    initial_condition, nx,dx, B, x

 ;renge_fed,64,1000,10,100,0.05,b_out,time

  ;_ FOR LOOP
  
  FOR n=0,T*jump DO BEGIN
    ;_Time step NOTE: CFL Condition
    kk = [dx^4*0.1, dx^2/R*0.1,dx^2*0.1/(9.*MAX(B^2))]
    dt = MIN(kk)
    ttt = ttt + dt 

    IF (n mod jump) eq 0 THEN BEGIN
	PRINT, ttt
       b_out[*,n/jump]=B
  	boundary, b_out, nx
	device, true=24, retain=2, decomposed=0 
       PLOT,x[6:nx+5], B[6:nx+5],XTITLE='x', ytitle='B(x,t='+ STRING(ttt)+')',title='R=1000', THICK=3, background=255, color=0, CHARSIZE=3,xstyle=1
       wait, 0.05
       i=i+1
       time[n/jump]=ttt
    ENDIF
 
	
   	derivatives,t_derive, f_derive, B, R, nx,dx,x,o_derive,cos_flag=0  	;call derivatives for B at t
	b_1=b-dt/2.*(R*t_derive + f_derive)	
    boundary, b_1, nx

				;calc half step B_1 
	derivatives,t_derive, f_derive, b_1, R, nx,dx,x,cos_flag=0		;call derivatives for B at t+1
	BB=B-dt*(R*t_derive + f_derive)		
    B = BB
    boundary, b, nx

  ENDFOR
END








                            

