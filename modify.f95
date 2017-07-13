integer :: i,k,j
real(kind=8) :: y(80)
open(100,file="fort.90416577",status='old')
open(200,file="200_Pe_100")
do i=1,80
	read(100,*) y(i)
end do
do i =1,80
	j = i-1
	k = 10000 + 10000*j
	write(200,*) k,y(i)
end do
end
