module Gauss_Legendre

!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-5-30
!-----------------------------------------------------
!  Description :   高斯-勒让德计算数值积分模块
!    
!-----------------------------------------------------
!  Parameters  :
!      1.
!      2.
!----------------------------------------------------- 
!  Contains    :
!      1.   方法函数
!      2.   要计算的函数 
!-----------------------------------------------------

contains 

subroutine solve(func,s,a,b)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-5-30
!-----------------------------------------------------
!  Purpose   : Gauss-Legendre积分
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.    func  外部函数
!       2.    a,b   积分区间
!  Output parameters  :
!       1.    s 积分结果
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.  选用5节点的公式
!       2.  计算时先把一般区间变换到[-1,1]区间上利用公式
!----------------------------------------------------

implicit real*8(a-z)
external func

integer::i
real*8::node(5),w(5),t(5)

node=(/-0.9061798459,-0.5384693101,0,0.5384693101,0.9061798459/)
w=(/0.2369268851,0.4786286705,0.5688888889,&
  0.4786286705,0.2369268851/)
  
  t=(b+a)/2d0+(b-a)*node/2d0  

s=0

do i=1,5
  
  call func(fx,t(i))
  
  s=s+fx*w(i)

end do

s=s*(b-a)/2d0

end subroutine solve


subroutine fun1(fx,x)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-5-30
!-----------------------------------------------------
!  Purpose   :  待计算的函数
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.   x  自变量
!       2. 
!  Output parameters  :
!       1.  fx  因变量
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.
!       2.
!----------------------------------------------------
implicit real*8(a-z)
fx=dexp(x)*dcos(x)
end subroutine fun1


end module Gauss_Legendre


program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-5-30
!-----------------------------------------------------
!  Purpose   :  主函数
!    
!-----------------------------------------------------
!  In put data  files :
!       1.
!       2.
!  Output data files  :
!       1.  result.txt  计算结果存放文件
!       2.
!-----------------------------------------------------
!  Post Script :
!       1.
!       2.
!-----------------------------------------------------
use Gauss_Legendre

implicit real*8(a-z)

open(unit=11,file='result.txt')

pi=3.141592653589793238D0

call solve(fun1,s,0d0,pi)

write(11,101)s

101 format(/,T5,'Gauss-Legendre积分',//,&
            T3,'积分结果为：',F12.5  )

end program main