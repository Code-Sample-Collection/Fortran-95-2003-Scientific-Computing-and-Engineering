module Simpson
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-5-29
!-----------------------------------------------------
!  Description :   复合Simpson法模块
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

subroutine solve(func,s,a,b,n)
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

f=x**2+dsin(x)
end subroutine fun1


end module Simpson



program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-5-29
!-----------------------------------------------------
!  Purpose   :  复合Simpson法计算数值积分主函数
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
!       1.  可以通过公式预先估计在一定的精度下需要划分
!           多少节点
!-----------------------------------------------------
use Simpson

implicit real*8(a-z)


open(unit=11,file='result.txt')

write(11,101)

call solve(fun1,s,-2d0,2d0,40)

write(11,102)s
call solve(fun1,s,-2d0,2d0,80)

write(11,103)s

call solve(fun1,s,-2d0,2d0,200)

write(11,104)s


101 format(/,T5,'复合Simpson法计算数值积分',/)
102 format(T5,'n=40',F15.8)
103 format(T5,'n=80',F15.8)
104 format(T5,'n=200',F14.8)

end program main