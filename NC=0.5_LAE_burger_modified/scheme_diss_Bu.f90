 module variable
 implicit none
 integer,parameter :: n=200
 double precision :: dt,time,Tmax=22.0d0,dist,x0,x1,Amp2=0.0000d0,Amp1=5.0d0
 double precision :: Nc,C=1.00d0,kh,dx,Leng,kappa,eps=1.0d-08,fact,pi,kappa1
 double precision,dimension(1:n+1) ::disp,ddisp,predisp,exdisp,fdisp,d4disp
 double precision,dimension(n+1) :: disp_old,fdisp_old,ddisp_old,d4disp_old,res
 integer :: i,j,cnt,ii
 end module variable
 !------------------------------------------------------------------------------------------------------
 program noise
 use variable
 implicit none
 double precision,dimension(n+1) :: k1,k2,k3,k4
 character*11:: filestr
 real:: diss

 Nc=0.50d0
 Leng=2.0d0
 pi=ATAN(1.0d0)*4.0d0
 diss=0.0
 kh = 1.0d0

 dx=Leng/real(n)
 !kappa=kh/dx
 
 kappa = 0.250d0*pi

 
 dt=Nc*dx/C

 do i=1,n+1
    
    dist=(i-1)*dx - 1.0d0
    if (dist .lt. 0) then 
	disp(i) = -sin(pi*dist) - 0.50d0*(dist**3)
    else
	disp(i) = -sin(pi*dist) - 0.50d0*(dist**3) + 1.0d0 
    end if

 enddo
 

 cnt=0
 time=0.0d0

 do
 !   call ddisp_dx

    if (mod(cnt,nint(800.0))==0) then
     filestr(1:1)='e'
     filestr(2:2)=ACHAR(48 + mod((cnt/10000),10))
     filestr(3:3)=ACHAR(48 + mod((cnt/1000),10))
     filestr(4:4)=ACHAR(48 + mod((cnt/100),10))
     filestr(5:5)=ACHAR(48 + mod((cnt/10),10))
     filestr(6:6)=ACHAR(48 + mod((cnt/1),10))
     filestr(7:10)='.dat'
     open(1,file=filestr,status='replace')
     write(1,*)'VARIABLES= DISTANCE,DISPLACEMENT'
     do i=1,n+1
       dist=real(i-1)*dx - 1.0d0
       write(1,*)dist,',',disp(i)
     enddo
    endif

   do i=1,n+1
     predisp(i)=disp(i)
   enddo

   call ccdderivatives!ddisp_dx

   do i=1,n+1
     IF(ABS(d4disp(i)).gt.1.0d0) THEN 
        k1(i)=(fdisp(i))
     ELSE
        k1(i)=(fdisp(i))
     ENDIF
   enddo

   do i=1,n+1
     disp(i)=predisp(i)+(dt)*k1(i)
   enddo ! FIRST STAGE COMPLETE

   call ccdderivatives!ddisp_dx

   do i=1,n+1
     IF(ABS(d4disp(i)).gt.1.0d0) THEN
        k2(i)=(fdisp(i))
     ELSE
        k2(i)=fdisp(i)
     ENDIF
   enddo

   do i=1,n+1
     disp(i)=(3.0d0/4.0d0)*predisp(i) + (1.0d0/4.0d0)*disp(i) +(dt/4.0d0)*k2(i)
   enddo ! SECOND STAGE COMPLETE

   call ccdderivatives!ddisp_dx

   do i=1,n+1
     IF(ABS(d4disp(i)).gt.1.0d0) THEN
        k3(i)=(fdisp(i))
     ELSE
        k3(i)= fdisp(i)
     ENDIF
   enddo

   do i=1,n+1
     disp(i)=(1.0d0/3.0d0)*predisp(i)+(2.0d0/3.0d0)*disp(i)+ (2.0d0/3.0d0)*dt*k3(i)
   enddo ! THIRD STAGE COMPLETE

   call ccdderivatives!ddisp_dx

   do i=1,n+1
     IF(ABS(d4disp(i)).gt.1.0d0) THEN
        k4(i)= (fdisp(i))
     ELSE
        k4(i)= fdisp(i)
     ENDIF
   enddo

   do i=1,n+1
     !disp(i)=predisp(i)+(dt/6.0d0)*(k1(i)+2.0d0*k2(i)+2.0d0*k3(i)+k4(i))
   enddo ! FOURTH STAGE COMPLETE


   cnt=cnt+1
   time=real(cnt)*dt
   write(*,*) time , dt
   if (time>=Tmax) exit
 enddo

  end program noise
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ccdderivatives
use VARIABLE
implicit none

double precision,dimension(0:2) :: beta , d , alpha , w
double precision,dimension(0:2,0:2)::a
double precision,dimension(n+1):: fu 
double precision,dimension(-2:n+4)::f_plus, f_minus , f_plus_half , f_minus_half
double precision  :: epsilon , maxv , temp
integer :: l , k , r , z 

f_plus_half = 0.0d0
f_minus_half = 0.0d0
fdisp = 0.0d0
!fu = 0.50d0*disp*disp
f_plus = 0.0d0
f_minus = 0.0d0
epsilon = 0.0000010d0

a(0,0) = 2.0d0/6.0d0 ;  a(0,1) = -7.0d0/6.0d0 ;  a(0,2) = 11.0d0/6.0d0
a(1,0) = -1.0d0/6.0d0 ; a(1,1) = 5.0d0/6.0d0 ;   a(1,2) = 2.0d0/6.0d0
a(2,0) = 2.0d0/6.0d0 ;  a(2,1) = 5.0d0/6.0d0 ;   a(2,2) = -1.0d0/6.0d0


r = 3
beta = 0.0d0
d = 0.0d0
alpha = 0.0d0
w = 0.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!maxv = maxval(disp)


do i=1,n+1
        fu(i) = 1.0d0*disp(i)
	f_plus(i) = 0.50d0 * ( fu(i) + 1.0d0*disp(i) )
        f_minus(i) = 0.50d0 * ( fu(i) - 1.0d0*disp(i) )
enddo

f_plus(0) = f_plus(n) ; f_minus(0) = f_minus(n)
f_plus(-1) = f_plus(n-1) ; f_minus(-1) = f_minus(n-1)
f_plus(-2) = f_plus(n-2) ; f_minus(-2) = f_minus(n-2)
f_plus(n+2) = f_plus(2) ; f_minus(n+2) = f_minus(2)
f_plus(n+3) = f_plus(3) ; f_minus(n+3) = f_minus(3)
f_plus(n+4) = f_plus(4) ; f_minus(n+4) = f_minus(4)


do i=0,n+1
	beta(0) = (13.0d0/12.0d0)*(f_plus(i-2) - 2*f_plus(i-1) + f_plus(i))**2 
	beta(0) = beta(0) + (1/4.0d0)*(f_plus(i-2)-4.0d0*f_plus(i-1)+3.0d0*f_plus(i))**2
	beta(1) = (13.0d0/12.0d0)*(f_plus(i-1) - 2*f_plus(i) + f_plus(i+1))**2 
	beta(1) = beta(1) + (1/4.0d0)*(f_plus(i-1) - f_plus(i+1))**2
	beta(2) = (13.0d0/12.0d0)*(f_plus(i) - 2*f_plus(i+1) + f_plus(i+2))**2 
	beta(2) = beta(2) + (1/4.0d0)*(3.0d0*f_plus(i) - 4.0d0*f_plus(i+1) + f_plus(i+2))**2

	d(0) = 	1.0d0/10.0d0 ; d(1) = 3.0d0/5.0d0 ; d(2) = 3.0d0/10.0d0

	alpha(0) = d(0)/(epsilon + beta(0))**2 ; alpha(1) = d(1)/(epsilon + beta(1))**2 ; alpha(2) = d(2)/(epsilon + beta(2))**2

	w(0) = alpha(0)/(alpha(0)+alpha(1)+alpha(2))
	w(1) = alpha(1)/(alpha(0)+alpha(1)+alpha(2))  
	w(2) = alpha(2)/(alpha(0)+alpha(1)+alpha(2))

	f_plus_half(i) = w(0)*(a(0,0)*f_plus(i-2) + a(0,1)*f_plus(i-1) + a(0,2)*f_plus(i)) &
			 + w(1)*(a(1,0)*f_plus(i-1) + a(1,1)*f_plus(i) + a(1,2)*f_plus(i+1)) &
			 + w(2)*(a(2,0)*f_plus(i) + a(2,1)*f_plus(i+1) + a(2,2)*f_plus(i+2))
	

enddo

!!! Working on f_minus_half!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!!!!!!!!!!!

a(0,0) = -1.0d0/6.0d0 ; a(0,1) = 5.0d0/6.0d0; a(0,2) = 2.0d0/6.0d0
a(1,0) = 2.0d0/6.0d0 ; a(1,1) = 5.0d0/6.0d0 ; a(1,2) = -1.0d0/6.0d0
a(2,0) = 11.0d0/6.0d0 ; a(2,1) = -7.0d0/6.0d0  ; a(2,2) = 2.0d0/6.0d0




do i=0,n+1
	beta(0) = ((13.0d0/12.0d0)*(f_minus(i+1) - 2*f_minus(i+2) + f_minus(i+3)))**2 
	beta(0) = beta(0) + (1/4.0d0)*(3*f_minus(i+1)-4*f_minus(i+2)+f_minus(i+3))**2
	beta(1) = (13.0d0/12.0d0)*(f_minus(i) - 2*f_minus(i+1) + f_minus(i+2))**2 + (1/4.0d0)*(f_minus(i) - f_minus(i+2))**2
	beta(2) = (13.0d0/12.0d0)*(f_minus(i-1) - 2*f_minus(i) + f_minus(i+1))**2 
	beta(2) = beta(2) + (1/4.0d0)*(f_minus(i-1) - 4*f_minus(i) + 3*f_minus(i+1))**2

	d(0) = 	1.0d0/10.0d0 ; d(1) = 3.0d0/5.0d0 ; d(2) = 3.0d0/10.0d0

	alpha(0) = d(0)/(epsilon + beta(0))**2 ; alpha(1) = d(1)/(epsilon + beta(1))**2 ; alpha(2) = d(2)/(epsilon + beta(2))**2

	w(0) = alpha(0)/(alpha(0)+alpha(1)+alpha(2))
	w(1) = alpha(1)/(alpha(0)+alpha(1)+alpha(2))  
	w(2) = alpha(2)/(alpha(0)+alpha(1)+alpha(2))

	f_minus_half(i) =  w(2)*(a(0,0)*f_minus(i-1) + a(0,1)*f_minus(i) + a(0,2)*f_minus(i+1)) &
			 + w(1)*(a(1,0)*f_minus(i) + a(1,1)*f_minus(i+1) + a(1,2)*f_minus(i+2)) &
			 + w(0)*(a(2,0)*f_minus(i+1) + a(2,1)*f_minus(i+2) + a(2,2)*f_minus(i+3))
	

enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do i=0,n+1
        f_plus(i) =  f_plus_half(i) + f_minus_half(i)                                            !!!!!!!!  f _ j+1/2
end do

do i=1,n+1
	f_minus(i) = f_plus(i-1) 								!!!!!!!!  f _ j-1/2
        fdisp(i) =   -(1.0d0/dx)*(f_plus(i)-f_minus(i))
end do
        return

       end subroutine ccdderivatives

