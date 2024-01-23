module com_Gauss
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-5-31
!-----------------------------------------------------
!  Description :  复合高斯勒让德方法计算定积分模块
!    
!-----------------------------------------------------
!  Contains    :
!      1.   复合方法函数
!      2.   单区间5点高斯勒让德积分函数
!      3.   要计算的函数
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
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  复合5点高斯勒让德积分函数
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1. func 外部函数
!       2. a,b 积分区间
!       3. n 区间划分个数
!  Output parameters  :
!       1.  s 积分结果
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.   需要调用单区间高斯勒让德公式函数
!       2.
!----------------------------------------------------
implicit real*8(a-z)
external func
integer::n,i

hstep=(b-a)/n

s=0
do i=1,n

    c=a+(i-1)*hstep
    d=a+i*hstep

   call GL(func,s1,c,d)

   s=s+s1
end do

end subroutine solve


subroutine GL(func,s,a,b)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-5-30
!-----------------------------------------------------
!  Purpose   :5点 Gauss-Legendre积分
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

end subroutine GL


subroutine fun1(f,x)
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
!       1.  f 因变量
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.
!       2.
!----------------------------------------------------
implicit real*8(a-z)
f=1/x
end subroutine fun1


end module com_Gauss


program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-5-29
!-----------------------------------------------------
!  Purpose   :  复合高斯法主函数
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
use com_Gauss

implicit real*8(a-z)
integer::n

open(unit=11,file='result.txt')

write(11,101)

!划分4个区间
call solve(fun1,s,1d0,2d0,4)

write(11,102)s

101 format(/,T5,'复合高斯数值积分',/)
102 format(T5,'积分结果：',F15.8)


end program main
