use omp_lib
implicit none
integer :: m,n,t,tmax,tequl,cnt,k,i,iseed,nx,ny
real (kind = 8) :: gam,gasdev,ran2,vp
real (kind = 8), dimension (:) , allocatable :: xx,xy,vx,vy,ax,ay,the
real (kind =8) :: box,Lx,Ly,dt,temp,rho,d
external :: gasdev,ran2
nx=32 ; ny = 32 ; dt = 0.01D0  ; gam = 100.0d0
m = nx*ny
allocate(xx(m),xy(m),the(m))
allocate(vx(m),vy(m))
allocate(ax(m),ay(m))
ax = 0.D0 ; ay = 0.D0 ; vx = 0.D0 
xx = 0.D0 ; xy = 0.D0 ; iseed = -887  
vy = 0.D0 ; the = 0.D0
temp = 1.4D0 ; vp = 90.0d0  ; d = temp/gam	      	
open (400,file='finvel.txt') 
open (500,file='accel.txt') 
tmax = 200 ; tequl = 1000 
rho = 0.80d0
call initial(xx,xy,vx,vy,the,nx,ny,m,Lx,Ly,temp,dt,iseed,rho)
do t = 1,tmax
	print*, t
	call force(ax,ay,xx,xy,m,Lx,Ly)
	call integrate(ax,ay,xx,xy,vx,vy,the,m,dt,Lx,Ly,iseed,gam,d,vp)
end do
do i=1,m
	write(400,*) xx(i),xy(i)
	write(500,*) ax(i),ay(i)
end do	
print*, 'Lx',Lx,'Ly',Ly,'rho',rho
end

subroutine initial(xx,xy,vx,vy,the,nx,ny,m,Lx,Ly,temp,dt,iseed,rho)	
implicit none
integer ,intent(in) :: nx,ny,iseed,m
integer :: i,j,k,p,g
real (kind = 8),intent(inout) :: xx(m),xy(m),vx(m),vy(m),the(m),temp
real (kind = 8)::gasdev,ax,ay,ran2
real (kind = 8), intent (in) :: dt,rho
real (kind = 8), intent (out) :: Lx,Ly
REAL (kind = 8), PARAMETER :: Pi = 3.14159265
external :: gasdev,ran2
ax = sqrt(2.0d0/(sqrt(3.0d0)*rho))
ay = ax*sqrt(3.0d0)/2
Lx = nx*ax						! length of box in x direction
Ly = ny*ay						! length of box in y direction
print*, 'ax',ax,'ay',ay
g = 1
do i =0,m-1						! position starts from zero and there are m particles
	if (g.le.m) then
		xx(g) = ax*i				! setting the intial layer positions (y=0)
		xy(g) = 0.0d0
		p = int((i)/nx)
		k = mod(p,2)
		if (i.gt.(nx-1)) then			! once the initial layer has finished, ie when i reaches the value nx
			xx(g) = xx(g-nx) + ax*0.5D0	! equilateral triangle,third point lies in the perpendicular bisector of the base
			if (k.eq. 0) then		
				xx(g) = xx(g) - ax*1.0d0
			end if
			xy(g) = xy(g-nx) + ay		! in each layer the y value remains the same and is incremented by ay
		end if
	end if	
	g = g + 1
end do
! assigning random orientations to each particle
do i =1,m
	the(i) = 2.0d0*Pi*ran2(iseed)			! the angle is a radom number between 0 and 2 pi
end do
do i=1,m
	vx(i) = cos(the(i))
	vy(i) = sin(the(i))
end do
end subroutine initial

subroutine force(ax,ay,xx,xy,m,Lx,Ly)
implicit none
integer , intent (in) :: m
integer :: i,j
real (kind = 8), intent (inout) :: ax(m),ay(m)
real (kind = 8), intent (in) :: xx(m),xy(m),Lx,Ly
real (kind = 8) :: rsqd,ff,r2i,r6i,rcut,rcut2,virij,rx,ry
rx = 0.D0 ; ry = 0.D0 
ax = 0.D0 ; ay = 0.D0 
rcut = 1.12246205			! cut off distance (rc = 4*sigma)
rcut2 = rcut*rcut
call omp_set_num_threads(10)
!$omp parallel do private(rx,ry,rsqd,r2i,r6i,virij,ff,j) schedule(dynamic) 
! calculation of pair-wise distances and forces on each particle
do i=1,m-1
	do j=i+1,m
		rsqd = 0.D0
		rx = xx(i) - xx(j)
		rx = rx-Lx*nint(rx/Lx)	! distance btw i and nearest periodic image of j
		ry = xy(i) - xy(j)
		ry = ry-Ly*nint(ry/Ly)						
		rsqd = rsqd + rx*rx +ry*ry 
	if (rsqd.lt.rcut2) then;
		r2i = 1/rsqd
		r6i = r2i**3
		virij = 48*(r6i*r6i-0.5D0*r6i)
                ff = virij*r2i		
		! force on the particle
		!$omp atomic
		ax(i) = ax(i) + rx*ff
		!$omp atomic				
		ay(i) = ay(i) + ry*ff	
		!$omp atomic									
		ax(j) = ax(j) - rx*ff
		!$omp atomic
		ay(j) = ay(j) - ry*ff					
	end if			
	end do
end do
!$omp end parallel do
end subroutine force

subroutine integrate(ax,ay,xx,xy,vx,vy,the,m,dt,Lx,Ly,iseed,gam,d,vp)
implicit none
integer , intent (in) :: m,iseed
integer :: i,j	
real (kind = 8), intent (in) :: ax(m),ay(m),dt,Lx,Ly,d,gam,vp  ! d and dr translational and rotational diffusion coefficients
real (kind = 8), intent (inout) :: vx(m),vy(m),xx(m),xy(m),the(m)
real (kind = 8) :: g1,g2,g3,gasdev,dr
external :: gasdev
dr = 3.0d0*d						! in the low Reynolds number regime
do i=1,m	
	g1 = gasdev(iseed)				! Gaussian random number with 0 mean and unit variance
	xx(i) = xx(i) + dt*(vp+ ax(i))/gam + sqrt(2*d)*g1	
        g2 = gasdev(iseed)				
	xy(i) = xy(i) + dt*(vp+ ay(i))/gam + sqrt(2*d)*g2        
        g3 = gasdev(iseed)
	the(i) = the(i) + sqrt(2*dr)*g3			! orientation of the particle is modified due to the the fluctuations 
	xx(i) = xx(i) - Lx*anint((xx(i)/Lx)-(0.5D0))				
	xy(i) = xy(i) - Ly*anint((xy(i)/Ly)-(0.5D0))
	vx(i) = vp*cos(the(i))				
	vy(i) = vp*sin(the(i))			
end do
end subroutine integrate

