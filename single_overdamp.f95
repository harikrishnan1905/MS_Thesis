implicit none
integer i,istep,nstep,iseed
double precision dt,gam,T,c0,c1,c2,sigr,sigv,crv,time,dx,dy,dz,x,y,z,vx,vy,vz,zetarx,zetary,zetarz,zetavx              
double precision zetavy,zetavz,g1,g2,Lx,Ly,Lz,gasdev,ran2,x0,y0,z0,rsq(1000),fac
external gasdev,ran2
dt = 0.01d0
gam = 100.0d0
T = 1.0d0                    
iseed = -887
nstep = 1000
rsq = 0.0d0
sigr = (dt*T*2.0d0)/gam     
sigr = sqrt(sigr)       
do i = 1,10000
      print*,i 
      x0 = ran2(iseed)
      x = x0
      do istep = 1,nstep
	g1 = gasdev(iseed)
	zetarx = sigr*g1		
	x = x + zetarx
	dx = x0 - x
	rsq(istep) = rsq(istep) + dx*dx 
      enddo
enddo
time = 0.0d0
do i = 1,nstep
time = time + dt
fac = 2.0d0*(gam*time - 1.0d0 + exp(-time*gam))/(gam*gam)
write(55,*)time,rsq(i)/10000.0d0,fac
enddo
end
