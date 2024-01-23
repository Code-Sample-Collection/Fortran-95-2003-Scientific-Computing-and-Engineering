module bezier
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
!  Description :   Bézier曲线模块
!    
!----------------------------------------------------- 
!  Contains    :
!      1.     方法函数  solve
!      2.
!-----------------------------------------------------

contains

subroutine solve(r1,r2,r3,r4,n,t,x,y)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :   Bézier曲线方法函数
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.   r1  r4 端点
!       2.   r2  r3 控制点
!       3.   n  要计算的节点个数
!       4    t  节点参数    t 在0到1之间
!  Output parameters  :
!       1.   x  坐标
!       2.   y  坐标
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.
!       2.
!----------------------------------------------------
implicit real*8(a-z)

integer::n,i

real*8::r1(2),r2(2),r3(2),r4(2)
real*8::c1(4),c2(4)
real*8::t(n),x(n),y(n)

bx=3d0*(r2(1)-r1(1))
cx=3d0*(r3(1)-r2(1))-bx

dx=r4(1)-r1(1)-bx-cx

by=3d0*(r2(2)-r2(1))
cy=3d0*(r3(2)-r2(2))-by
dy=r4(2)-r1(2)-by-cy


!已经解算出系数
do i=1,n
    x(i)=r1(1)+bx*t(i)+cx*(t(i))**2+dx*(t(i))**3
    y(i)=r1(2)+by*t(i)+cy*(t(i))**2+dy*(t(i))**3
end do

end subroutine solve

end module bezier


program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.05.13
!-----------------------------------------------------
!  Purpose   : Bézier曲线主函数
!    
!-----------------------------------------------------
!  Post Script :
!       1.
!       2.
!-----------------------------------------------------
use bezier

implicit real*8(a-z)

integer::n,i

real*8::r1(2),r2(2),r3(2),r4(2)
real*8::t(0:20),x(0:20),y(0:20)


open(unit=11,file='result.txt')


do i=0,20
  
  t(i)=1d0/20*i

end do

!端点
r1=(/1d0,1d0/)
r4=(/2d0,2d0/)

!控制点
r2=(/1d0,3d0/)
r3=(/3d0,3d0/)


call solve(r1,r2,r3,r4,21,t,x,y)

write(11,101)
write(11,102)((t(i),x(i),y(i)),i=0,20)


101 format(/,T16,'Bézier曲线',/)

102 format(3F12.6)

end program main