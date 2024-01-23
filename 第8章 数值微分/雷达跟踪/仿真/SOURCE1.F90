module radardiff
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010.08.11
!-----------------------------------------------------
!  Description :   
!  
!  Post Script :
!      1.
!      2. 
!
!-----------------------------------------------------

contains
subroutine solve
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.08.04
!-----------------------------------------------------
!  Purpose     : 微分求速函数
! 
!  Post Script :
!       1.
!       2.
!       3.    
!-----------------------------------------------------

implicit real*8(a-z)

integer i
real*8 t(-4:4),x(-4:4),y(-4:4),z(-4:4)

open(unit=20,file='file1.txt')

open(unit=21,file='file2.txt')


!从测量资料里读数据，程序运行时，需要把测量资料放在当前目录下
do i=-4,4
 read(20,*)t(i),x(i),y(i),z(i)
end do 


!该段注释代码 是把读的数据重写一遍，检查看是否有错误
!open(unit=30,file='file2.txt')
!do i=-4,4
! write(30,'(f3.1,3f16.8)')t(i),x(i),y(i),z(i)
!end do 


!采用五点公式 微分求速
vx=x(-2)-8*x(-1)+8*x(1)-x(2)
vx=vx/12/0.1

vy=y(-2)-8*y(-1)+8*y(1)-y(2)
vy=vy/12/0.1

vz=z(-2)-8*z(-1)+8*z(1)-z(2)
vz=vz/12/0.1


ax=11*x(-1)-20*x(0)+6*x(1)+4*x(2)-x(3)
ax=ax/12/0.01

ay=11*y(-1)-20*y(0)+6*y(1)+4*y(2)-y(3)
ay=ay/12/0.01

az=11*z(-1)-20*z(0)+6*z(1)+4*z(2)-z(3)
az=az/12/0.01


!该段注释代码为验证程序-------------------------
!即严格按照运动学求其速度与加速度，为使之与数值法比较
!第一行为速度，第二行为加速度
write(21,'(3f16.8)')vx,vy,vz
write(21,'(3f18.8)')ax,ay,az

t1=6d0

write(21,*)'---'
write(21,*)'the real velocity'

write(21,*)(dcos(t1)+2*t1+3d0)
write(21,*)(dexp(t1)+dcos(t1)+5d0)
write(21,*)(-dsin(t1)+3*t1**2-2*t1)

write(21,*)'---'
write(21,*)'the real acceleration'

write(21,*)(-dsin(t1)+2d0)
write(21,*)(dexp(t1)-dsin(t1))
write(21,*)(-dcos(t1)+6*t1-2d0)
!-------------------------------------------------

end subroutine solve

end module radardiff

program main 
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose     : 主函数
!
!  Post Script :
!       1.
!       2.
!    
!-----------------------------------------------------

use radardiff
call solve

end