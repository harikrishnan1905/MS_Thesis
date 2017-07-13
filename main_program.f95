! molecular dynamic simulation of active brownian particles subjected to WCA potential interactions and 
! confined between two boundaries in 2 dimensions
use omp_lib
implicit none
integer :: m,t,tmax,i,j,iseed,n,mtot,nwall,m_start,m_end
real (kind = 8) :: gasdev,ran2,vp,dm,a
real (kind = 8), dimension (:) , allocatable :: xx,xy,ax,ay,the
real (kind =8) :: L,dt,temp,rho,d,pe
external :: gasdev,ran2
!open (400,file='length100_310E5.txt') 
open(50,file='initial.txt')
!######################################################################################################
!######################################################################################################
!######################################################################################################
!###	n = 45 ==> m = 2025 ==> l = sqrt(2025/0.8) = 50.311529494	  a = 1.1180339887498949    ###
!###	n = 64 ==> m = 4096 ==> l = 71.554175279993274			  a = 1.1180339887498949    ###	 
!###	n = 90 ==> m = 8100 ==> l = sqrt(8100/0.8) = 100.62305898749054   a = 1.1180339887498949    ###
!###	n = 100 ==> m = 10000 &&& l = 112.0  rho = m/l**2 = 0.797193878   a = 1.1200000000000001    ###	rcut = 1.12246205
!###	n = 110 ==> m = 12100     l = 123.05d0 rho = 0.81366844100763458  a = 1.1186363636363637    ###
!###	n = 135 l = 151.10d0    rho = 0.81007533109283303 		  a = 1.1192592592592592    ###
!###    n = 124 l = 139.150d0   rho = 0.80691132148384603
!###    n = 115 l = 129.0d0     rho = 0.80854515954570039
!###    n = 80  l = 89.442719099991592 rho = 0.80d0
!###    n = 72  l = 80.498447189992433 rho = 0.80d0
!###    n = 180 l = 202.0d0            rho = 0.80286246446426823 
!###    n = 223 l = 250.0d0            rho = 0.80279999999999996
!###    n = 270 l = 301.86917696247161 rho = 0.80d0
!###    n = 160 l = 180.50000000000000 rho = 0.79563539260748462

!	n = 90  l = 100.62305898749054 rho = 0.80d0
! 	n = 112 l = 125.29806518554688 rho = 0.81326886165347345
!	n = 135 l = 151.10d0           rho = 0.81007533109283303
!	n = 157 l = 175.05000000000001 rho = 0.80d0
!	n = 180 l = 202.0d0            rho = 0.80286246446426823 
!######################################################################################################
!######################################################################################################
!************************************** PARAMETERS **************************************************##
!######################################################################################################

dt = 0.00001D0  					! time step for each updation
n = 64      						! number of particles in one dimension
m = n*n 						! total number of particles
dm = 1.12246205						! potential cut-off
Pe = 150.00d0						! velocity of individual particles
iseed = -887  						! seed for random number generation
temp = 1.0d0 						! effective temperature of the ambient medium initially	
vp = Pe  						! self propulsion velocity
d = 1.0!temp/gam	      				! translational diffusion coefficient
tmax = 1						! maximum number of time steps
rho = 0.80d0						! density
!l = 175.050d0 
l = sqrt(m/rho) 					! simulation box is square with side l
!l = 125.29806739988230
a = l/n							! a specifies the base of the triangle	
nwall = int(l/dm) + 1 					! number of particles at each boundary	
m_start = nwall + 1					! starting number of active particles
m_end = nwall + m					! ending number of active particles
mtot = m + 2*nwall					! total number of particles including the barrier ones
!rho = mtot/(l**2)
print*,l,rho,mtot,nwall
allocate(xx(mtot),xy(mtot),the(mtot))			! 'the(i)' specifies the direction of the i-th particle	
allocate(ax(mtot),ay(mtot))
the = 0.D0 ; ax = 0.D0 ; ay = 0.D0 
xx = 0.D0 ; xy = 0.D0 

!####################################################################################################
!************************************** MAIN ********************************************************
!####################################################################################################

call initial(xx,xy,the,m,mtot,n,l,dt,iseed,dm)						! initializing the system
do t = 1,tmax
	call force(ax,ay,xx,xy,mtot,L,m_start,m_end)					! calculating force btw particles
	call integrate(ax,ay,xx,xy,the,m,mtot,dt,L,iseed,d,vp,m_start,m_end)		! updating positions and directions
	if((t.gt. 10000).and.(mod(t,10000).eq. 0)) then	
		do i=1,mtot
			write(400,*) xx(i),xy(i)					! writing configurations into file
		end do
	end if	
end do	
end

!####################################################################################################
!**************************************** INITIALIZATION ********************************************
!####################################################################################################

subroutine initial(xx,xy,the,m,mtot,n,l,dt,iseed,d)	
implicit none
integer ,intent(in) :: iseed,m,n,mtot
integer :: i,j,k,p,g,mw,m_start,m_end,nwall
real (kind = 8),intent(inout) :: xx(mtot),xy(mtot),the(mtot)
real (kind = 8)::gasdev,a,ran2
real (kind = 8), intent (in) :: dt,l,d
REAL (kind = 8), PARAMETER :: Pi = 3.14159265
external :: gasdev,ran2					
nwall = (mtot-m)/2
a = l/n								! x direction separation distance btw adjacent active particle centers
! lower barrier
g = 1
do i = 0,nwall-1						! assigining positions for lower boundary particles
	xx(g) = i*d 						! separation btw centers of adjacent particles is set to be r_cut (d)
	xy(g) = 0.0d0						! so that other particles can not penetrate through
	the(g) = 0.0d0						! particle diameter is sigma = 1 
	g = g + 1
end do
! m particle in triangular lattice in lxl box 
g = nwall+1	
mw = m + nwall		
do i =0,m-1						
	if (g.le.mw) then
		xx(g) = a*i					! base of triangle is a
		xy(g) = 2.0d0					! first active layer 2.0 above the lower boundary
		p = int((i)/n)
		k = mod(p,2)
		if (i.gt.(n-1)) then			
			xx(g) = xx(g-n) + a*0.5D0		! triangular lattice property
			if (k.eq. 0) then		
				xx(g) = xx(g) - a*1.0d0		
			end if
			xy(g) = xy(g-n) + a			! height of triangular lattice
		end if
	end if	
	g = g + 1
end do
g = nwall + 1
do i=1,m
	the(g) = 2.0d0*Pi*ran2(iseed)				! orientation of particles
	g = g + 1						! random values btw 0 and 2*pi
end do
! upper barrier
g = m + nwall + 1
xx(m + nwall + 1) = 0.0d0					! upper boundary particles
xy(m + nwall + 1) = xy(m+nwall+1-n) + 2.0d0			! set 2.0 above the last active row 	
do i = 0, nwall-1						! inorder to avoid initial large repulsion
	xx(g) = xx(m + nwall + 1) + i*d
	xy(g) = xy(m+nwall+1) 
	the(g) = 0.0d0	
	g = g + 1 
end do
do i=1,mtot
	write(50,*) xx(i),xy(i)
end do
end subroutine initial

!####################################################################################################
!************************************** FORCE CALCULATION ******************************************* 
!####################################################################################################

subroutine force(ax,ay,xx,xy,mtot,L,m_start,m_end)
implicit none
integer , intent (in) :: mtot,m_start,m_end
integer :: i,j
real (kind = 8), intent (inout) :: ax(mtot),ay(mtot)
real (kind = 8), intent (in) :: xx(mtot),xy(mtot),L
real (kind = 8) :: rsqd,ff,r2i,r6i,rcut,rcut2,virij,rx,ry
rx = 0.D0 ; ry = 0.D0 
ax = 0.D0 ; ay = 0.D0 
rcut = 1.12246205 !2**(1.0d0/6.0d0)			! cut off distance for the WCA potential is 2^(1/6)
rcut2 = rcut*rcut !1.259921054
call omp_set_num_threads(10)
!$omp parallel do private(rx,ry,rsqd,r2i,r6i,virij,ff,j) schedule(dynamic) 
do i=1,mtot-1
	do j=i+1,mtot
		rsqd = 0.D0
		rx = xx(i) - xx(j)
		rx = rx - l*anint(rx/l)
		ry = xy(i) - xy(j)	
		rsqd = rsqd + rx*rx +ry*ry 
	if (rsqd.le.rcut2) then;
		r2i = 1/rsqd
		r6i = r2i**3
		virij = 48*(r6i*r6i-0.5D0*r6i) 		
	        ff = virij*r2i	
	        if (ff .gt. 10000) then	
			print*,'yup mate it sucks', i,j,rx,ff,r2i,virij   ! just to indicate blowing up
			EXIT
		end if
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
do i = 1,m_start-1
	ax(i) = 0.0d0						! boundary particles are stationary
	ay(i) = 0.0d0
end do
do i = m_end + 1 , mtot
	ax(i) = 0.0d0
	ay(i) = 0.0d0
end do
end subroutine force

!####################################################################################################
!************************************** INTEGRATION *************************************************
!####################################################################################################

subroutine integrate(ax,ay,xx,xy,the,m,mtot,dt,L,iseed,d,vp,m_start,m_end)
implicit none
integer , intent (in) :: m,iseed,mtot,m_start,m_end
integer :: i,j,ini,fin
real (kind = 8), intent (in) :: ax(mtot),ay(mtot),dt,L,d,vp 	 	 ! d and dr translational and rotational diffusion coefficients
real (kind = 8), intent (inout) :: xx(mtot),xy(mtot),the(mtot)
real (kind = 8) :: g1,g2,g3,gasdev,dr
external :: gasdev
dr = 3.0d0*d								 ! in the low Reynolds number regime
ini = m_start 								 ! assuming that active particles also follow the same rules 	
fin = m_end
do i=1,mtot
	if((i .ge. ini).and.(i .le. fin)) then	
		g1 = gasdev(iseed)				
		xx(i) = xx(i) + dt*(vp*cos(the(i))+ ax(i))*d + sqrt(2*d*dt)*g1	
		g2 = gasdev(iseed)				
		xy(i) = xy(i) + dt*(vp*sin(the(i))+ ay(i))*d + sqrt(2*d*dt)*g2        
		g3 = gasdev(iseed)
		the(i) = the(i) + sqrt(2*dr*dt)*g3			
		xx(i) = xx(i) - L*anint((xx(i)/L)-(0.5D0))		! Periodic Boundary Conditions in the x direction
	else								! no restriction in the x direction
		xx(i) = xx(i)						! stationary boundary particles
	end if	
end do
end subroutine integrate

!#########################################################################################################

FUNCTION ran2(idum)
INTEGER*4 idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
REAL*8 ran2,AM,EPS,RNMX
PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.0d0/IM1,IMM1=IM1-1,IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791)
PARAMETER (NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2d-7,RNMX=1.0d0-EPS)
INTEGER*4 idum2,j,k,iv(NTAB),iy
SAVE iv,iy,idum2
DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
return
END

!###############################################################################################################

DOUBLE PRECISION FUNCTION gasdev(iseed)
IMPLICIT NONE
DOUBLE PRECISION r, v1, v2, fac, gset, ran2
INTEGER iset, iseed 
SAVE gset, iset
DATA iset/0/
100 	 IF (iset.EQ.0) THEN
		 v1 = 2.d0*ran2(iseed) - 1.d0
		 v2 = 2.d0*ran2(iseed) - 1.d0
		 r = v1**2 + v2**2
         IF (r.GE.1) GOTO 100
		 fac = sqrt(-2.d0*log(r)/r)
		 gset = v1*fac
		 gasdev = v2*fac
		 iset = 1
     	 ELSE
		 gasdev = gset
		 iset = 0
	 END IF
Return
End

! rcut is the WCA potential cut-off. It is not same as the particle diameter. Diameter is 'sigma' and is set to be unity. 'a' is the initial spacing between particles in the triangular lattice and is greater than the individual diameter of the particle.
! smallest distance between particles after the formation of solid is coming to be smaller than the individual particle diameter (are the particles penetrating into each other???)
