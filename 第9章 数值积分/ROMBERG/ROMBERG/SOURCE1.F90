module Romberg
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-5-29
!-----------------------------------------------------
!  Description :   Romberg方法模块
!    
!-----------------------------------------------------
!  Contains    :
!      1.     方法函数
!      2.     测试函数
!-----------------------------------------------------
!  Post Script :
!      1.
!      2. 
!-----------------------------------------------------

contains 

subroutine solve(func,s,a,b,tol)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-5-29
!-----------------------------------------------------
!  Purpose   :  Romberg方积分方法函数
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.  func 外部函数
!       2. 
!       3.  a,b积分区间
!       4.  tol  积分误差容限
!
!  Output parameters  :
!       1.   s  积分结果
!    
!
!----------------------------------------------------
implicit real*8(a-z)
external func

integer::i,j,k,m

!定义较大矩阵存放数据
!此处亦可定义可变维数组
!因Romberg积分收敛较快，实际上用不了很大的矩阵
real*8::T(1:50,0:49)


call func(fa,a)
call func(fb,b)

T(1,0)=(b-a)/2d0*(fa+fb)


do i=1,40
  
     s=0
    do j=1,2**(i-1)
      
      x1=a+(2*j-1)*(b-a)/(2d0**i)
      call func(f1,x1)
      s=s+f1
    end do
!计算T(i,1) 
     T(1,i)=T(1,i-1)/2d0+(b-a)*s/2d0**i
 
!开始外推     
     do m=1,i
        do k=i,1,-1
          
         temp=4d0**m*T(m,k)-T(m,k-1)
         
         T(m+1,k-1)=temp/(4d0**m-1)
        end do
     end do

!计算最后两相邻元素绝对误差
       del=dabs(T(m,0)-T(m-1,0))
!如果误差小于给定的容限，则迭代停止       
       if (del<tol) exit
       
end do

!最后一个对角线元素作为计算结果
      s=T(m,0)

end subroutine solve



subroutine fun1(f,x)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-5-29
!-----------------------------------------------------
!  Purpose   :  需要计算的函数
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.     x  自变量
!       2. 
!  Output parameters  :
!       1.     f  因变量
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.
!       2.
!----------------------------------------------------
implicit real*8(a-z)

f=x/(4+x**2)
end subroutine fun1


end module romberg



program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-5-29
!-----------------------------------------------------
!  Purpose   :  Romberg方法计算数值积分主函数
!    
!-----------------------------------------------------
!  In put data  files :
!       1.
!       2.
!  Output data files  :
!       1.   result.txt计算结果
!       2.
!-----------------------------------------------------
!  Post Script :
!       1.  
!-----------------------------------------------------
use Romberg

implicit real*8(a-z)
integer::n

open(unit=11,file='result.txt')

write(11,101)

call solve(fun1,s,0d0,1.5d0,1d-7)

write(11,102)s


101 format(/,T5,'Romberg方法计算数值积分',/)
102 format(T5,T5,'积分结果：',F15.8)


end program main