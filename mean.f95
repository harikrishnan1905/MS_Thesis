implicit none
integer :: i,p(60),k,j,cnt
real(kind=8) :: x(80),y(80),total,mean
open(100,file="200_Pe_100",status='old')
do i=1,80
	read(100,*) x(i),y(i)
end do
total = 0.0d0
 cnt = 0
do i = 1,80
	if(i.ge. 40) then
		total = total + y(i)
		cnt = cnt + 1
	end if
end do
mean = total/cnt
print*, mean
end
