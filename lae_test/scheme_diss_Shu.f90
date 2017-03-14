 module variable
 implicit none
 integer,parameter :: n=200
 double precision :: dt,time,Tmax=9.0d0,dist,x0,x1,Amp2=0.0000d0,Amp1=5.0d0
 double precision :: Nc,C=1.00d0,kh,dx,Len,kappa,eps=1.0d-08,fact,pi,kappa1
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
 Len=2.0d0
 pi=ATAN(1.0d0)*4.0d0
 diss=0.0
 kh = 1.0d0

 dx=Len/real(n)
 !kappa=kh/dx
 
 kappa = 0.250d0*pi

 
 dt=Nc*dx/C

 do i=1,n+1
    
    dist=real(i-1)*dx - 1.0d0
    disp(i) = -sin(pi*dist)
    

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
     write(1,*)'VARIABLES= DISTANCE,DISPLACEMENT,Forth_Der'
     do i=1,n+1
       dist=real(i-1)*dx
       write(1,*)dist,disp(i)
    endif

   do i=1,n+1
     predisp(i)=disp(i)
   enddo

   call ccdderivatives!ddisp_dx

   do i=1,n+1
     IF(ABS(d4disp(i)).gt.1.0d0) THEN 
        k1(i)=-C*(fdisp(i))
     ELSE
        k1(i)=-C*(fdisp(i))
     ENDIF
   enddo

   do i=1,n+1
     disp(i)=predisp(i)+(dt)*k1(i)
   enddo ! FIRST STAGE COMPLETE

   call ccdderivatives!ddisp_dx

   do i=1,n+1
     IF(ABS(d4disp(i)).gt.1.0d0) THEN
        k2(i)=-C*(fdisp(i))
     ELSE
        k2(i)=-C*fdisp(i)
     ENDIF
   enddo

   do i=1,n+1
     disp(i)=(3.0d0/4.0d0)*predisp(i) + (1.0d0/4.0d0)*disp(i) +(dt/4.0d0)*k2(i)
   enddo ! SECOND STAGE COMPLETE

   call ccdderivatives!ddisp_dx

   do i=1,n+1
     IF(ABS(d4disp(i)).gt.1.0d0) THEN
        k3(i)=-C*(fdisp(i))
     ELSE
        k3(i)=-C*fdisp(i)
     ENDIF
   enddo

   do i=1,n+1
     disp(i)=(1.0d0/3.0d0)*predisp(i)+(2.0d0/3.0d0)*disp(i)+ (2.0d0/3.0d0)*dt*k3(i)
   enddo ! THIRD STAGE COMPLETE

   call ccdderivatives!ddisp_dx

   do i=1,n+1
     IF(ABS(d4disp(i)).gt.1.0d0) THEN
        k4(i)=-C*(fdisp(i))
     ELSE
        k4(i)=-C*fdisp(i)
     ENDIF
   enddo

   do i=1,n+1
     disp(i)=predisp(i)+(dt/6.0d0)*(k1(i)+2.0d0*k2(i)+2.0d0*k3(i)+k4(i))
   enddo ! FOURTH STAGE COMPLETE

   do i=1,n+1
     exdisp(i)=Amp1*dexp(-0.1d0*(dist-C*time-x0)**2)*dcos(kappa*(dist-C*time-x0))+&
     Amp2*dexp(-(dist-C*time-x1)**2/eps)/dsqrt(pi*eps)
   enddo

   cnt=cnt+1
   time=real(cnt)*dt
   write(*,*) time
   if (time>=Tmax) exit
 enddo

  end program noise
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine ccdderivatives
use VARIABLE
implicit none


double precision,dimension(0:2) :: beta , d , alpha , w
double precision,dimension(0:2,0:2)::a
double precision,dimension(n+1)::f_plus,f_minus
double precision  :: epsilon 
integer :: l , k , r , z

f_plus = 0.0d0
f_minus = 0.0d0
fdisp = 0.0d0
epsilon = 0.0000010d0
a(0,0) = 2.0d0/6.0d0 ; a(0,1) = -7.0d0/6.0d0 ; a(0,2) = 11.0d0/6.0d0
a(1,0) = -1.0d0/6.0d0 ; a(1,1) = 5.0d0/6.0d0 ; a(1,2) = 2.0d0/6.0d0
a(2,0) = 2.0d0/6.0d0 ; a(2,1) = 5.0d0/6.0d0 ; a(2,2) = -1.0d0/6.0d0
r = 3
beta = 0.0d0
d = 0.0d0
alpha = 0.0d0
w = 0.0d0


do i=3,n-1
	beta(0) = (13.0d0/12.0d0)*(disp(i) - 2*disp(i+1) + disp(i+2))**2 + (1/4.0d0)*(3*disp(i)-4*disp(i+1)+disp(i+2))**2
	beta(1) = (13.0d0/12.0d0)*(disp(i-1) - 2*disp(i) + disp(i+1))**2 + (1/4.0d0)*(disp(i-1) - disp(i+1))**2
	beta(2) = (13.0d0/12.0d0)*(disp(i-2) - 2*disp(i-1) + disp(i))**2 + (1/4.0d0)*(disp(i-2) - 4*disp(i-1) + 3*disp(i))**2

	d(0) = 	1.0d0/10.0d0 ; d(1) = 3.0d0/5.0d0 ; d(2) = 3.0d0/10.0d0

	alpha(0) = d(0)/(epsilon + beta(0))**2 ; alpha(1) = d(1)/(epsilon + beta(1))**2 ; alpha(2) = d(2)/(epsilon + beta(2))**2

	w(0) = alpha(0)/(alpha(0)+alpha(1)+alpha(2))
	w(1) = alpha(1)/(alpha(0)+alpha(1)+alpha(2))  
	w(2) = alpha(2)/(alpha(0)+alpha(1)+alpha(2))

	do k=0,r-1
		do l=0,r-1
			f_plus(i) = f_plus(i) + w(k)*a(k,l)*disp(i-k+l)
		enddo
	enddo
enddo


do i=2,2
	beta(0) = (13.0d0/12.0d0)*(disp(i) - 2*disp(i+1) + disp(i+2))**2 + (1/4.0d0)*(3*disp(i)-4*disp(i+1)+disp(i+2))**2
	beta(1) = (13.0d0/12.0d0)*(disp(i-1) - 2*disp(i) + disp(i+1))**2 + (1/4.0d0)*(disp(i-1) - disp(i+1))**2
	beta(2) = (13.0d0/12.0d0)*(disp(n) - 2*disp(i-1) + disp(i))**2 + (1/4.0d0)*(disp(n) - 4*disp(i-1) + 3*disp(i))**2

	d(0) = 	1.0d0/10.0d0 ; d(1) = 3.0d0/5.0d0 ; d(2) = 3.0d0/10.0d0

	alpha(0) = d(0)/(epsilon + beta(0))**2 ; alpha(1) = d(1)/(epsilon + beta(1))**2 ; alpha(2) = d(2)/(epsilon + beta(2))**2

	w(0) = alpha(0)/(alpha(0)+alpha(1)+alpha(2))
	w(1) = alpha(1)/(alpha(0)+alpha(1)+alpha(2))  
	w(2) = alpha(2)/(alpha(0)+alpha(1)+alpha(2))



	do k=0,r-1
		do l=0,r-1
			z = i-k+l
			if(z .eq. 0) then
				z = n
			end if
			f_plus(i) = f_plus(i) + w(k)*a(k,l)*disp(z)
		enddo
	enddo
enddo

do i=1,1
	beta(0) = (13.0d0/12.0d0)*(disp(i) - 2*disp(i+1) + disp(i+2))**2 + (1/4.0d0)*(3*disp(i)-4*disp(i+1)+disp(i+2))**2
	beta(1) = (13.0d0/12.0d0)*(disp(n) - 2*disp(i) + disp(i+1))**2 + (1/4.0d0)*(disp(n) - disp(i+1))**2
	beta(2) = (13.0d0/12.0d0)*(disp(n-1) - 2*disp(n) + disp(i))**2 + (1/4.0d0)*(disp(n-1) - 4*disp(n) + 3*disp(i))**2

	d(0) = 	1.0d0/10.0d0 ; d(1) = 3.0d0/5.0d0 ; d(2) = 3.0d0/10.0d0

	alpha(0) = d(0)/(epsilon + beta(0))**2 ; alpha(1) = d(1)/(epsilon + beta(1))**2 ; alpha(2) = d(2)/(epsilon + beta(2))**2

	w(0) = alpha(0)/(alpha(0)+alpha(1)+alpha(2))
	w(1) = alpha(1)/(alpha(0)+alpha(1)+alpha(2))  
	w(2) = alpha(2)/(alpha(0)+alpha(1)+alpha(2))


	do k=0,r-1
		do l=0,r-1
			z = i-k+l	
			if(z .eq. 0) then
				z= n
			end if
			if(z .eq. -1) then
				z= n-1
			end if
			f_plus(i) = f_plus(i) + w(k)*a(k,l)*disp(z)
		enddo
	enddo
enddo

do i=n+1,n+1
	beta(0) = (13.0d0/12.0d0)*(disp(i) - 2*disp(2) + disp(3))**2 + (1/4.0d0)*(3*disp(i)-4*disp(2)+disp(3))**2
	beta(1) = (13.0d0/12.0d0)*(disp(i-1) - 2*disp(i) + disp(2))**2 + (1/4.0d0)*(disp(i-1) - disp(2))**2
	beta(2) = (13.0d0/12.0d0)*(disp(i-2) - 2*disp(i-1) + disp(i))**2 + (1/4.0d0)*(disp(i-2) - 4*disp(i-1) + 3*disp(i))**2

	d(0) = 	1.0d0/10.0d0 ; d(1) = 3.0d0/5.0d0 ; d(2) = 3.0d0/10.0d0

	alpha(0) = d(0)/(epsilon + beta(0))**2 ; alpha(1) = d(1)/(epsilon + beta(1))**2 ; alpha(2) = d(2)/(epsilon + beta(2))**2

	w(0) = alpha(0)/(alpha(0)+alpha(1)+alpha(2))
	w(1) = alpha(1)/(alpha(0)+alpha(1)+alpha(2))  
	w(2) = alpha(2)/(alpha(0)+alpha(1)+alpha(2))


	do k=0,r-1
		do l=0,r-1
			z=i-k+l
			if(z .eq. n+2) then
				z = 2
			end if

			if(z .eq. n+3) then
				z = 3
			end if

			f_plus(i) = f_plus(i) + w(k)*a(k,l)*disp(z)
		enddo
	enddo
enddo
	

do i=n,n
	beta(0) = (13.0d0/12.0d0)*(disp(i) - 2*disp(i+1) + disp(2))**2 + (1/4.0d0)*(3*disp(i)-4*disp(i+1)+disp(2))**2
	beta(1) = (13.0d0/12.0d0)*(disp(i-1) - 2*disp(i) + disp(i+1))**2 + (1/4.0d0)*(disp(i-1) - disp(i+1))**2
	beta(2) = (13.0d0/12.0d0)*(disp(i-2) - 2*disp(i-1) + disp(i))**2 + (1/4.0d0)*(disp(i-2) - 4*disp(i-1) + 3*disp(i))**2

	d(0) = 	1.0d0/10.0d0 ; d(1) = 3.0d0/5.0d0 ; d(2) = 3.0d0/10.0d0

	alpha(0) = d(0)/(epsilon + beta(0))**2 ; alpha(1) = d(1)/(epsilon + beta(1))**2 ; alpha(2) = d(2)/(epsilon + beta(2))**2

	w(0) = alpha(0)/(alpha(0)+alpha(1)+alpha(2))
	w(1) = alpha(1)/(alpha(0)+alpha(1)+alpha(2))  
	w(2) = alpha(2)/(alpha(0)+alpha(1)+alpha(2))

	

	do k=0,r-1
		do l=0,r-1
			z = i-k+l
			if(z .eq. n+2) then
				z = 2
			end if
			f_plus(i) = f_plus(i) + w(k)*a(k,l)*disp(z)
		enddo
	enddo
end do 


do i=1,n+1

	if(i .eq. 1) then
		f_minus(i) = f_plus(n)                                                         

	else
	        f_minus(i) = f_plus(i-1)                                                         !!!!!!!  f _ j-1/2
    	
	end if
end do

do i=1,n+1
	
	fdisp(i) = (1.0d0/dx)*(f_plus(i) - f_minus(i))

end do

        return

       end subroutine ccdderivatives

