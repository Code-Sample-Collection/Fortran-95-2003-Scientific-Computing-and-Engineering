module richardson
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010.08.04
!-----------------------------------------------------
!  Description :   Richardson外推模块
!  
!  Post Script :
!      1.
!      2. 
!
!-----------------------------------------------------
!  Contains    :
!      1.
!      2.
!-----------------------------------------------------
!  Parameters  :
!      1.
!      2.
!----------------------------------------------------- 

contains

subroutine solve
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.08.04
!-----------------------------------------------------
!  Purpose     : Richardson外推
! 
!  Post Script :
!       1.
!       2.
!       3.    
!-----------------------------------------------------

implicit real*8(a-z)
!外推三级
real*8 h(3)

real*8 g(3,3)

h(1)=0.1d0
h(2)=h(1)/2
h(3)=h(2)/2

x=0.5d0

!求的底层差商
call fun(y1,x+h(1))
call fun(y2,x-h(1))
g(1,1)=(y1-y2)/2/h(1)
call fun(y1,x+h(2))
call fun(y2,x-h(2))
g(1,2)=(y1-y2)/2/h(2)
call fun(y1,x+h(3))
call fun(y2,x-h(3))
g(1,3)=(y1-y2)/2/h(3)


!进行外推
call gf(y,g(1,1),g(1,2),1)
g(2,1)=y

call gf(y,g(1,2),g(1,3),1)
g(2,2)=y

call gf(y,g(2,1),g(2,2),1)
g(3,1)=y

!记录整个外推结果
open(unit=20,file='result.txt')
write(20,101)
write(20,'(3f16.8)') g(1,1),g(1,2),g(1,3)
write(20,'(2f16.8)') g(2,1),g(2,2)
write(20,'(f16.8)')g(3,1)

101 format(/,T16,'Richardson外推结果',/)

!!----这段注释程序为验证程序-
! x=0.5d0
! y=x*2*dexp(-x)-x**2*dexp(-x)
! write(20,*)'----'
! write(20,*)y
!!----

end subroutine solve



subroutine fun(y,x)
implicit real*8(a-z)
!所求之函数
y=x**2*dexp(-x)
end subroutine fun



subroutine gf(y,a1,a2,n)
implicit real*8(a-z)
integer n
!外推公式
y=a2-a1*(1/2d0)**(2*n)
temp=1-(1/2d0)**(2*n)
y=y/temp
end subroutine gf

end module richardson

program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.08.04
!-----------------------------------------------------
!  Purpose     :  主函数
!
!  Post Script :
!       1.
!       2.
!    
!-----------------------------------------------------

use richardson

!调用外推函数
call solve

end program main

