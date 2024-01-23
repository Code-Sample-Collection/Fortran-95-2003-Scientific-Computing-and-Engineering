module Gauss_Hermite

!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-5-30
!-----------------------------------------------------
!  Description :   Gauss-Hermite计算数值积分模块
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

subroutine solve(func,s)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-5-30
!-----------------------------------------------------
!  Purpose   : Gauss-Hermite积分
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
real*8::node(5),Ae(5)

node=(/-2.020182870456086d0,&
       -0.958572464613819d0,&
                        0d0,&
      0.958572464613819d0,&
      2.020182870456086d0  /)
Ae=(/1.1814886255360d0,&
     0.9865809967514d0,&
     0.9453087204829d0,&
     0.9865809967514d0,&
     1.1814886255360d0  /)
  
s=0d0

do i=1,5

  call func(fx,node(i))
  s=s+Ae(i)*fx
  
end do

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
fx=dexp(-x**2)
end subroutine fun1


end module Gauss_Hermite


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
use Gauss_Hermite

implicit real*8(a-z)

open(unit=11,file='result.txt')


call solve(fun1,s)

write(11,101)s

101 format(/,T5,'Gauss-Hermite积分',//,&
            T3,'积分结果为：',F12.8  )

end program main