implicit none
integer :: t,tmax,q,i,e,j,jmax
real (kind = 8) :: d,s,c1,c0,dt,z1,z2,z11,z22,u,w,xin,x,x1,sigr,sigv,crv,time,fac,dx,gam
real (kind = 8) , dimension (:) , allocatable :: rsqd
integer , dimension ( 8 ) :: value
integer , dimension ( : ) , allocatable :: seed
tmax = 10000
gam = 100.0d0						
dt = 0.01d0						
jmax = 10000						
allocate(rsqd(tmax+1))
call date_and_time ( values = value)
call random_seed ( size =q )	
allocate ( seed ( 1:q ) )
seed ( : ) = value ( 8 )
call random_seed ( put = seed )
open(50,file='out.txt')
	c0 = exp(-gam*dt)
      c1 = (1.0d0-c0)/(dt*gam)	
      sigr = dt*(2.0d0-(3.0d0-4.0d0*exp(-gam*dt)+exp(-2.0d0*gam*dt))/(dt*gam))/gam
      sigr = sqrt(sigr)
      sigv = (1.0d0-exp(-2.0d0*dt*gam))
      sigv = sqrt(sigv)
      crv = (1.0d0-exp(-gam*dt))*(1.0d0-exp(-gam*dt))/(gam*sigr*sigv)  
xin = 0.0d0  
x = 0.0d0    
x1 = 0.0d0   
do j=1,tmax
	rsqd(j) = 0.0d0
end do
do j=1,jmax
	print*, j						
	call random_number(w)					
	call random_number(u)					
	x1 = sqrt(-2.0*log(w))*cos(2*3.14*u)			
  	x = x1			
	do t=1,tmax
		call random_number(d)
		call random_number(s)
		z1 = sqrt(-2.0*log(d))*cos(2*3.1415d0*s)	
		z11 =  z1*sigr		
		x = x + z11				
		dx = x-x1
		rsqd(t) = rsqd(t) + dx*dx		
	end do							
end do			
time = 0.0d0					
do i =1,tmax
	time = time + dt
        fac = 2.0d0*(gam*time - 1.0d0 + exp(-time*gam))/(gam*gam)	
	write(50,*) time,fac,rsqd(i)/jmax    
end do
end	
