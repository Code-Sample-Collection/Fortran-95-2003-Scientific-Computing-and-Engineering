module autoSimpson
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-5-29
!-----------------------------------------------------
!  Description :   自动变步长Simpson法模块
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

subroutine solve(func,s,a,b,tol,n)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-5-29
!-----------------------------------------------------
!  Purpose   :  自动变步长Simpson积分方法函数
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
!       2.   n  实际区间划分个数
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.  需要调用复合Simpson公式
!       2.
!----------------------------------------------------
implicit real*8(a-z)
external func
integer::n,i

!初始划分40个子区间
n=40

!最大允许划分20次
do i=1,20
call simp(func,s,a,b,n)

n=n*2
call simp(func,s1,a,b,n)

del=dabs(s-s1)

!满足精度后就停止循环
if(del<tol)  exit
end do

s=s1

end subroutine solve

subroutine simp(func,s,a,b,n)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-5-29
!-----------------------------------------------------
!  Purpose   :  复合Simpson公式计算数值积分
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.         func 为外部子程序
!       2.         
!       3.         a,b积分区间
!       4.         n 区间划分个数
!  Output parameters  :
!       1.         s 积分结果
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.
!       2.
!----------------------------------------------------
implicit real*8(a-z)
external func
integer::n,k

s=0d0
h=(b-a)/n/2d0

call func(f1,a)
call func(f2,b)

s=f1+f2

!k=0 情况
call func(f1,a+h)
s=s+4d0*f1

do k=1,n-1

t1=a+(2d0*k+1)*h

t2=a+2d0*k*h

call func(f3,t1)
call func(f4,t2)

s=s+f3*4d0+f4*2d0  

end do

s=s*h/3d0

end subroutine simp


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

f=1d0/(x**3-2*x-5)
end subroutine fun1


end module autoSimpson



program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-5-29
!-----------------------------------------------------
!  Purpose   :  自动变步长Simpson法计算数值积分主函数
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
use autoSimpson

implicit real*8(a-z)
integer::n

open(unit=11,file='result.txt')

write(11,101)

call solve(fun1,s,0d0,2d0,1d-7,n)

write(11,102)n,s




101 format(/,T5,'自动变步长Simpson法计算数值积分',/)
102 format(T5,'区间划分等份为：',I5,/,T5,'积分结果：',F15.8)


end program main