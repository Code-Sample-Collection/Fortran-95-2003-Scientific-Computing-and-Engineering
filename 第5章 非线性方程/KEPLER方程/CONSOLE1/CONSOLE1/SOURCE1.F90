module kepler
!----------------------------------------module coment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-4
!-----------------------------------------------------
!  Purpose   :  计算kepler方程
!  
!-----------------------------------------------------
!  Parameters  :
!      1. MAX最大允许迭代次数
!      2. tol误差容限
!-----------------------------------------------------
!  Contains    :
!      1.  newton牛顿法计算
!      2.  picard不动点迭代法计算
!-----------------------------------------------------
implicit real*8(a-z)
integer::MAX=200
real*8::eps=1d-12

contains

subroutine newton(E,M,ec,i)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-4
!-----------------------------------------------------
!  Purpose   :  采用牛顿迭代法计算kepler方程
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1. M平近点角
!       2. ec轨道偏心率
!  Output parameters  :
!       1. E偏近点角
!       2. i迭代次数
!
!----------------------------------------------------
!  Post Script :
!       1.计算时，需要注意角度与弧度之间的换算
!       2.本来牛顿法需要提供方程函数与导函数，这里就直接在计算中给出计算格式
!----------------------------------------------------
implicit real*8(a-z)

integer i
parameter(pi=3.141592653589793d0)


M0=M

M0=M0*pi/180

E0=M0

!最大允许迭代50次
do i=1,MAX

!方程函数
  f=E0-ec*dsin(E0)-M0 
!导函数 
  df=1d0-ec*dcos(E0)
  
  dE=-f/df
  
  E1=E0+de
  
  
  E0=E1
  
   if (dE<eps) exit
   
end do

   E=E1*180/pi

end subroutine newton


subroutine picard(E,M,ec,i)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-4
!-----------------------------------------------------
!  Purpose   : 采用不动点迭代法计算kepler方程
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1. M平近点角
!       2. ec轨道偏心率
!  Output parameters  :
!       1. E偏近点角
!       2. i迭代次数
!  Common parameters  :
!
!----------------------------------------------------


implicit real*8(a-z)

integer i

parameter(pi=3.141592653589793d0)

M0=M*pi/180
!把角度化为弧度

E0=M0

!最大允许迭代次
do i=1,MAX

   E1=M0+ec*dsin(E0)
   tol=E1-E0
   E0=E1
   
   if (tol<eps) exit
   
end do

   E=E1*180/pi
!弧度还原为角度

end subroutine picard


end module kepler



program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :   计算kepler方程
!    
!-----------------------------------------------------
!  Output data files  :
!       1.  kepler.txt  计算结果
!       2.
!-----------------------------------------------------
!  Post Script :
!       1.    非线性方程的参数可以由外部读入
!       2.
!-----------------------------------------------------
use kepler
implicit real*8(a-z)

integer j,k

open(unit=101,file='kepler.txt')
!以文件保存结果

write(101,501)
501 format(30x,'牛顿迭代法计算结果',/,&
           10x,'M',T35,'E',T59,'iter')
!iter为迭代次数           
!此句为说明各量意义

do k=1,8

 M=k*10d0

 ec=0.01D0

 call newton(E,M,ec,j)
!调用牛顿迭代法迭代法函数
 write(101,*)M,E,j

end do



write(101,502)
502 format(///,30x,'不动点迭代计算结果',/,&
            10x,'M',T35,'E',T59,'iter')
!iter为迭代次数           
!此句为说明各量意义

do k=1,8

 M=k*10d0

 ec=0.01D0

 call picard(E,M,ec,j)
!调用牛顿迭代法迭代法函数

write(101,*)M,E,j

end do

end program main