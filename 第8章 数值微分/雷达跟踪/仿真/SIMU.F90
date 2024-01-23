program main
!该文件为仿真程序，目的是产生数据
!假设目标运动学为  x(t)=sin(t)+t**2+3*t
!y方向运动学为 y(t)=exp(t)+sin(t)+5*t
!z方向运动学为 z(t)=cos(t)+t**3-t**2
implicit real*8(a-z)

integer i

real*8 t(-4:4),x(-4:4),y(-4:4),z(-4:4)


h=0.1d0

!建立文件用以存放离散数据
open(unit=file1,file='file1.txt')

do i=-4,4
   
  t(i)=6d0+i*h
  x(i)=dsin(t(i))+t(i)**2+3*t(i)
  y(i)=dexp(t(i))+dsin(t(i))+5*t(i)
  z(i)=dcos(t(i))+t(i)**3-t(i)**2
!记录下离散的函数数据，写入文件
  write(file1,'(f3.1,3f16.8)')t(i),x(i),y(i),z(i)
end do

end