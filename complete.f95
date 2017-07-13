use omp_lib
implicit none
integer :: m,n,t,tmax,tequl,cnt,k,i
real (kind = 8), dimension (:) , allocatable :: xx,xy,xz,vx,vy,vz,ax,ay,az
real (kind =8) :: box,l,dt,temp,rho,en,vir,ekin,press,etot,top,vsumx,vsumy,vsum,ptot
n = 10 ; m = n**3 ; l = 10.350D0 ; dt = 0.01D0 !; box = l-(l/n)
allocate(xx(m),xy(m),xz(m))
allocate(vx(m),vy(m),vz(m))
allocate(ax(m),ay(m),az(m))
ax = 0.D0 ; ay = 0.D0 ; az = 0.D0
xx = 0.D0 ; vx = 0.D0 ; xz = 0.D0 
xy = 0.D0 ; en = 0.D0 ; vir = 0.D0
vy = 0.D0 ; vz = 0.D0 
press = 0.D0 ; etot = 0.D0 ; top = 0.D0
vsumx = 0.D0 ; vsumy = 0.D0 ; vsum = 0.D0
ptot = 0.D0 ; cnt = 0
open (100,file='energy.txt')
open (200,file='temp.txt')
open (300,file='pressure.txt')
open (400,file='finalvel.txt') 
tmax = 20000 ; tequl = 10000 ; temp = 4.0D0 
call initial(xx,xy,xz,vx,vy,vz,m,n,l,temp,dt)
rho = m/l**3
print*, rho
do t =1,tmax
	print*, t
	call force(ax,ay,az,xx,xy,xz,m,l,en,vir)
	call integrate(ax,ay,az,xx,xy,xz,vx,vy,vz,m,dt,l,n,ekin,vsumx,vsumy,vsum)
	if (t.eq.tmax) then
		do i =1,m
			write(400,*) vx(i)
		end do	
	end if	
	etot = ekin + en
	write(100,*) t,vsumx,vsumy,vsum
	if (t .lt. tequl) then	
		call rescale(vx,vy,vz,m,temp,top)
		write(200,*) t,top
		
	end if	
	if (t .ge. tequl) then
		if (mod(t,100).eq. 0) then
			call pressure (press,temp,rho,xx,xy,xz,m,l,n)
			ptot = ptot + press
			cnt = cnt +1
		end if	
	end if
end do
ptot = ptot/cnt
write(300,*) rho,ptot
end

subroutine initial(xx,xy,xz,vx,vy,vz,m,n,l,temp,dt)	
implicit none
integer ,intent(in) :: m,n
integer :: q,i,j,k,p
real (kind = 8),intent(inout) :: xx(m),xy(m),xz(m),vx(m),vy(m),vz(m),temp
real (kind = 8):: r,s,t,vsqd,lambda,vsumx,vsumy,vsumz,vxsum,vysum,vzsum,b,box
real (kind = 8), intent (in) :: dt,l
integer , dimension ( 8 ) :: value
integer , dimension ( : ) , allocatable :: seed
! random seed assignment
call date_and_time ( values = value)
call random_seed ( size =q )	
allocate ( seed ( 1:q ) )
seed ( : ) = value ( 8 )
call random_seed ( put = seed )
open(50,file='inipos.txt')
open(70,file='inivel.txt')
xx = 0.D0 ; xy = 0.D0 ; xz = 0.D0 
box = 0.D0
p = 1
b = l/n
box = l-(l/dble(n))
! particles into simple cubic lattice from 0 to l
do i=0,n-1
	do j=0,n-1
		do k=0,n-1
			if (p<(	m+1)) then;
				xx(p) = i*b
				xy(p) = j*b
				xz(p) = k*b
			end if
			p = p+1
		end do
	end do
end do
! assigning random velocities to the particles
! set the center of mass velocity to be zero and rescale the velocities 
! to get the desired temperature 
vx = 0.D0 ; vy = 0.D0 ; vz = 0.D0
vsumx = 0.D0 ; vsumy = 0.D0 ; vsumz = 0.D0 ; vsqd = 0.D0
do i=1,m
	call random_number(r)
	vx(i) = r-0.5D0
	vsumx = vsumx + vx(i)	
	call random_number(s)
	vy(i) = s-0.5D0
	vsumy = vsumy + vy(i)	
	call random_number(t)
	vz(i) = t-0.5D0		
	vsumz = vsumx + vz(i)
	vsqd = vsqd + vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i)		
end do
vsumx = vsumx/m
vsumy = vsumy/m
vsumz = vsumz/m
lambda = sqrt(3*m*temp/vsqd)	
vxsum = 0.D0 ; vysum = 0.D0 ; vzsum = 0.D0
vsqd = 0.D0
do i=1,m
	vx(i) = (vx(i) - vsumx)*lambda
	vy(i) = (vy(i) - vsumy)*lambda
	vz(i) = (vz(i) - vsumz)*lambda		
	vxsum = vxsum + vx(i)
	vysum = vysum + vy(i)
	vzsum = vzsum + vz(i)
	vsqd = vsqd + vx(i)**2 + vy(i)**2 + vz(i)**2
end do		
vsqd = vsqd/dble(3*m)
vxsum = vxsum/m
vysum = vysum/m
vzsum = vzsum/m
temp = vsqd	
print*, 'temp',temp
do i=1,m
	write(50,*) i,xx(i),xy(i),xz(i)
	write(70,*) vx(i)
end do
end subroutine initial



subroutine integrate(ax,ay,az,xx,xy,xz,vx,vy,vz,m,dt,l,n,ekin,vsumx,vsumy,vsum)
implicit none
integer , intent (in) :: m,n
integer :: i,j	
real (kind = 8), intent (in) :: ax(m),ay(m),az(m),dt,l
real (kind = 8), intent (inout) :: vx(m),vy(m),vz(m),xx(m),xy(m),xz(m),ekin
real (kind = 8) :: vsqd,vxt,vyt,vzt,vsumz
real (kind = 8), intent (out) :: vsumx,vsumy,vsum
! Verlet leap-frog algorithm and updating positions and velocities
! applying PBC, inorder to confine the particles inside the box
! calculating KE per particle
ekin = 0.D0
vsqd = 0.D0
do i=1,m								
	vxt = vx(i)
	vyt = vy(i)
	vzt = vz(i)
	vx(i) = vx(i) + dt*ax(i)
	vy(i) = vy(i) + dt*ay(i)
	vz(i) = vz(i) + dt*az(i)
	xx(i) = xx(i) + dt*vx(i)
	xy(i) = xy(i) + dt*vy(i)
	xz(i) = xz(i) + dt*vz(i)				
	xx(i) = xx(i) - l*anint((xx(i)/l)-(0.5D0))
	xy(i) = xy(i) - l*anint((xy(i)/l)-(0.5D0))
	xz(i) = xz(i) - l*anint((xz(i)/l)-(0.5D0))	
	vsqd = vsqd + vx(i)**2 + vy(i)**2 + vz(i)**2 					
end do
ekin = vsqd/dble(2*m)			! kinetic energy per particle
! checking whether the center of mass velocity is zero or not
! coming out to be very very small
vsumx = 0.D0
vsumy = 0.D0
do i=1,m
	vsumx = vsumx + vx(i)
	vsumy = vsumy + vy(i)
	vsumz = vsumz + vz(i)
end do
vsumx = vsumx/m
vsumy = vsumy/m
vsumz = vsumz/m
vsum = 0.D0
do i=1,m
	vsum = vsum + vsumx**2 + vsumy**2 + vsumz**2
end do
vsum = sqrt(vsum)
vsum = vsum/m
end subroutine integrate


subroutine force(ax,ay,az,xx,xy,xz,m,l,en,vir)
implicit none
integer , intent (in) :: m
integer :: i,j
real (kind = 8), intent (inout) :: ax(m),ay(m),az(m),en,vir
real (kind = 8), intent (in) :: xx(m),xy(m),xz(m),l
real (kind = 8) :: rsqd,ff,r2i,r6i,rcut,rcut2,r2c,r6c,ecut,virij,rx,ry,rz,enij
rx = 0.D0 ; ry = 0.D0 ; rz = 0.D0
ax = 0.D0 ; ay = 0.D0 ; az = 0.D0
rcut = 4.0D0			! cut off distance (rc = 4*sigma)
rcut2 = rcut*rcut
r2c = 1/rcut2
r6c = r2c**3
ecut = 4*r6c*(r6c-1.0D0)	! shifting term
en = 0.D0
vir  = 0.D0
call omp_set_num_threads(10)
!$omp parallel do private(rx,ry,rz,rsqd,r2i,r6i,virij,ff,j,enij) reduction(+:vir,en) schedule(dynamic) 
! calculation of pair-wise distances and forces on each particle
do i=1,m-1
	do j=i+1,m
		rsqd = 0.D0
		rx = xx(i) - xx(j)
		rx = rx-l*nint(rx/l)	! distance btw i and nearest periodic image of j
		ry = xy(i) - xy(j)
		ry = ry-l*nint(ry/l)
		rz = xz(i) - xz(j)			
		rz = rz-l*nint(rz/l)						
		rsqd = rsqd + rx*rx +ry*ry + rz*rz	! why rsqd not a reduction variable (rsqd = 0.0 in the start for each thread)
	if (rsqd.le.rcut2) then;
		r2i = 1/rsqd
		r6i = r2i**3
                enij = 4*(r6i*r6i-r6i) - ECUT		! using the truncated and shifted potential
		virij = 48*(r6i*r6i-0.5D0*r6i)
                en = en + enij			
                vir = vir + virij		
                ff = virij*r2i		
		! force on the particle
		!$omp atomic
		ax(i) = ax(i) + rx*ff
		!$omp atomic				
		ay(i) = ay(i) + ry*ff
		!$omp atomic
		az(i) = az(i) + rz*ff	
		!$omp atomic									
		ax(j) = ax(j) - rx*ff
		!$omp atomic
		ay(j) = ay(j) - ry*ff
		!$omp atomic
		az(j) = az(j) - rz*ff						
	end if			
	end do
end do
!$omp end parallel do
en = en/m 		! potential energy per particle 
end subroutine force



subroutine rescale(vx,vy,vz,m,temp,top)
implicit none
integer , intent (in) :: m
real(kind=8), intent (inout) :: vx(m),vy(m),vz(m)
real (kind = 8) , intent (in) :: temp
real (kind = 8) , intent(out) :: top
real :: lambda,vsqd
integer :: i,j
vsqd = 0.D0
do i=1,m
	vsqd = vsqd + vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i)
end do	 
vsqd = vsqd/m
lambda = sqrt((3*temp)/vsqd)
do i=1,m
	vx(i) = vx(i)*lambda
	vy(i) = vy(i)*lambda
	vz(i) = vz(i)*lambda		
end do
vsqd = 0.D0
do i=1,m
	vsqd = vsqd + vx(i)**2 + vy(i)**2 + vz(i)**2
end do
vsqd = vsqd/dble(3*m)
top = vsqd
end subroutine rescale




subroutine pressure (press,temp,rho,xx,xy,xz,m,l,n)
implicit none
integer , intent (in) :: m,n
integer :: i,j
real (kind = 8), intent (in) :: xx(m),xy(m),xz(m),temp,l,rho
real (kind = 8) , intent (out) :: press
real (kind = 8) :: v,vir,rsqd,r2i,r6i,rcut,rcut2,rci,rc2i,rc6i,nn,ptail
real (kind = 8) :: rx,ry,rz
press = 0.D0
! tail correction to the pressure
v = l**3
vir = 0.D0
rcut = 4.D0
rci = 1/rcut
rcut2 = rcut*rcut
rc2i = 1/rcut2
rc6i = rc2i*rc2i*rc2i
nn = (rho**2)*3.14D0*32.D0*(1.D0/9.D0)
ptail = nn*rci*rc2i*(rc6i-dble(3/2))
! calculating pressure
do i=1,m-1
	do j=i+1,m
		rsqd = 0.D0
		rx = xx(i) - xx(j)
		rx = rx-l*anint(rx/l)	
		ry = xy(i) - xy(j)
		ry = ry-l*anint(ry/l)
		rz = xz(i) - xz(j)			
		rz = rz-l*anint(rz/l)						
		rsqd = rsqd + rx*rx +ry*ry + rz*rz		
		if (rsqd.le.rcut2) then;
			r2i = 1/rsqd
			r6i = r2i*r2i*r2i  
			vir = vir + 48*r6i*(r6i-0.5D0)
		end if	
	end do
end do
vir = vir/3
press = (rho*temp) + (vir/v) + ptail
end subroutine pressure
