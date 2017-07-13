implicit none
integer :: i,j,n,k,m,p,no,n2o,n3o,n4o,nwall,num_pt,n5o
real (kind = 8) :: xx(2620800),xy(2620800),xx1(32760),xy1(32760)
real (kind = 8) , dimension(:) , allocatable :: xy2,h
integer , dimension(:) , allocatable :: grid_num,q
real (kind = 8) :: l,d,x_cut
open(12345,file='200_Pe_100.txt',status='old')
n = 80			! number of configurations
m = 2620800 		! total number of particles in all configurations
p = 32760		! number of particles in a single configuration
no = 100		! starting file number configurations
n2o = no + 400		! starting number for grided configurations
n3o = n2o + 400
n4o = n3o + 400
!nwall = 180		! number of boundary particles
l = 202.0d0 		! length of the box
x_cut = 1.12246205	! diameter of a particle
d = l/2	
nwall = int(l/x_cut) + 1
allocate(q(nwall),h(nwall))
do i=1,m
	read(12345,*) xx(i),xy(i)
end do		
call split(xx,xy,m,n)
do i = 1,n
	print*, i
	q = 0 
	open(no)
	rewind(no)
	n5o = 90416577
	write(55555,*) no,n2o,n3o,n4o
	num_pt = 0
	do j = 1,p
		read(no,*) xx1(j),xy1(j)	
	end do
	call grids(xx1,xy1,p,nwall,d,x_cut,n2o,num_pt)
	call array_allocation(num_pt,grid_num,xy2)
	rewind(n2o)
	do k =1,num_pt
		read(n2o,*) grid_num(k),xy2(k)
	end do
	call sort(grid_num,xy2,num_pt,n3o)
	rewind(n3o)
	do k =1,num_pt
		read(n3o,*) grid_num(k),xy2(k)
	end do		
	call height(grid_num,xy2,num_pt,n4o,q,nwall)	
	rewind(n4o)
	do k =1,nwall 
		read(n4o,*) h(k)
	end do
	call width(h,nwall,l,n5o)
	deallocate (grid_num,xy2)
	close(no)
	close(n2o)	
	no = no +1
	n2o = n2o +1
	n3o = n3o +1
	n4o = n4o +1
end do
contains
!#####################################################################
subroutine split(xx,xy,m,n)
implicit none
integer , intent (in) :: m,n
real (kind = 8) , intent(in) :: xx(m),xy(m)
integer :: i,k,j,p
p = 100
do i=0,n-1
	k = i*32760
	do j = k+1,k+32760
		write(p,*) xx(j),xy(j)
	end do
	p = p + 1	
end do
end subroutine split
!#####################################################################
subroutine grids(xx1,xy1,p,nwall,d,x_cut,n2o,num_pt)
implicit none
integer , intent (in) :: p,nwall,n2o
integer , intent (inout) :: num_pt
real (kind =8 ) ,intent (in) :: xx1(p),xy1(p),d,x_cut
integer :: k,i
num_pt = 0
do i =1,p
	if (xy1(i) .le. d) then 
		k = int(xx1(i)/x_cut) + 1
		write(n2o,*) k,xy1(i)
		num_pt = num_pt + 1
	end if			
end do
end subroutine grids
!#####################################################################
subroutine array_allocation(num_pt,grid_num,xy2)
implicit none
integer , intent(in) :: num_pt
integer , intent(inout) , dimension(:) ,allocatable :: grid_num
real (kind = 8) , intent(inout) , dimension(:) ,allocatable :: xy2
allocate(grid_num(num_pt),xy2(num_pt))
end subroutine array_allocation
!#####################################################################
subroutine sort(grid_num,xy2,num_pt,n3o)
implicit none
integer :: i,ii,jj
integer , intent (in) :: num_pt,n3o
integer , intent (inout) :: grid_num(num_pt)
real (kind = 8) , intent (inout) :: xy2(num_pt)
real (kind = 8) :: temp,tempy
logical :: swapped
DO jj = num_pt-1, 1, -1
	swapped = .FALSE.
	DO ii = 1, jj
		IF (grid_num(ii) > grid_num(ii+1)) THEN
			temp = xy2(ii)
			tempy = grid_num(ii)
			xy2(ii) = xy2(ii+1)
			grid_num(ii) = grid_num(ii+1)
			xy2(ii+1) = temp
			grid_num(ii+1) = tempy
			swapped = .TRUE.
		END IF
	END DO
	IF (.NOT. swapped) EXIT  
END DO
do i = 1,num_pt
	write(n3o,*) grid_num(i),xy2(i)
end do
end subroutine sort
!#####################################################################
subroutine height(p,xy,m,n4o,q,n)
implicit none
integer , intent (in) :: m,n4o,n
integer , intent (inout) :: p(m),q(n)
real (kind = 8) , intent (inout) :: xy(m)
real (kind = 8) :: temp,dist
logical :: me, swapped  
integer :: dum,start,en,i,j,ii,jj
do i =1,n
	q(i) = 0	
end do
do i=1,m
	dum = p(i)
	q(dum) = q(dum) + 1 		! number of particles in each column
end do
start = q(1)
en = 1
do j = 1,n
	do i =1,m
		if(p(i) .eq. j) then			! sorting each column
			! sorting the values in each column
		 	DO jj = start-1, en, -1
				swapped = .FALSE.
				DO ii = en, jj
			      		IF (xy(ii) > xy(ii+1)) THEN
						temp = xy(ii)
						xy(ii) = xy(ii+1)
						xy(ii+1) = temp
						swapped = .TRUE.
			      		END IF
			    	END DO
      			        IF (.NOT. swapped) EXIT  
		  	END DO
		end if
	end do
	en = start + 1	
	start = start + q(j+1)
end do	
start = 1
en = q(1)
do j =1,n
	me = .true.
	do i =1,m
		if (me) then
			if(p(i) .eq. j) then
				dist = xy(i+1)-xy(i)
				if(dist .le. 1.0d0) then
				
				else
					me = .false.		
				end if
			end if		
		end if	
		if (.not. me) exit
	end do 	
	write(n4o,*) xy(i)
	start = en + 1
	en = en + q(j+1)	
end do
end subroutine height
!####################################################################################################
subroutine width(h,m,l,n5o)
implicit none
integer :: i
integer , intent (in) :: m,n5o
real (kind=8) , intent (in) :: h(m)
real (kind=8) , intent (in) :: l
real (kind = 8) :: dum,std,mean_height,total,var
do i = 1,m
	write(4545454,*) h(i) 	
end do 
total = sum(h)
mean_height = total/l
std = 0.0d0
do i = 1,m
	dum = (h(i)-mean_height)**2	
	var = var + dum
end do
var = var/l
std = sqrt(var)
write(n5o,*) std
end subroutine width
end
