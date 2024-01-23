module lagrange
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
!  Description :   拉格朗日插值法模块
!
!----------------------------------------------------- 
!  Contains    :
!      1.  solve 方法函数
!      2.   lag  单点插值
!-----------------------------------------------------
!  Post Script :
!      1.
!      2. 
!-----------------------------------------------------

contains

subroutine solve(n,x,y,m,t,ty)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  拉格朗日插值函数，可以直接针对向量插值
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.   n  节点个数
!       2.   x 节点自变量值
!       3.   y 节点因变量值
!       4.   m  要插值点的个数
!       5.   t 要插值点的值  为向量
!  Output parameters  :
!       1.   
!       2.   ty 插值点的结果  为向量  
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.
!       2.
!----------------------------------------------------
implicit real*8(a-z)
integer::n,m
integer:: i
real*8::x(n),y(n),t(m),ty(m)

do i=1,m

  call la(n,x,y,t(i),ty(i))

end do

end subroutine solve


subroutine  la(n,x,y,t,ty)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.05.11
!-----------------------------------------------------
!  Purpose   :  拉格朗日单点插值，为solve所调用
!    
!-----------------------------------------------------
!  Post Script :
!       1.
!       2.
!----------------------------------------------------
implicit real*8(a-z)
integer::n,i,k
!输入序列x,y  为N维向量

!要插值的点t  
!t对应的插值结果 ty
real*8::x(n),y(n),a(n)

!a(i)为各项


do i=1,n
    a(i)=y(i)
        
    do k=1,n
       
       if (k==i) cycle
       a(i)=a(i)*(t-x(k))/(x(i)-x(k))
       
    end do
end do

ty=0

do i=1,n
  ty=ty+a(i)
end do

end subroutine la


end module lagrange


program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   : 拉格朗日插值法主函数
!    
!-----------------------------------------------------
!  In put data  files :
!       1.
!       2.
!  Output data files  :
!       1.
!       2.
!-----------------------------------------------------
!  Post Script :
!       1.
!       2.
!-----------------------------------------------------
use lagrange

implicit real*8(a-z)
integer::N=4
real*8::x0(4),x(4),y(4),t0(3),t(3),ty(3),tyreal(3)
real*8,parameter::pi=3.141592653589793d0

open(unit=11,file='result.txt')

!插值节点及其函数值，化为弧度
x0=(/30,45,60,90/)
x=x0*pi/180
y=(/dsqrt(3d0)/2,dsqrt(2d0)/2,1d0/2,0d0/)

!待计算的点，化为弧度
t0=(/47,53,79/)
t=t0*pi/180

call solve(4,x,y,3,t,ty)


!调用系统逐元函数
tyreal=dcos(t)

write(11,101)
!输出标题

write(11,102)t0
!输出角度

write(11,103)ty
!输出插值结果

write(11,104)tyreal
!输出实际函数值

101 format(/,T22,'拉格朗日插值法',/)
102 format(T3,'角度：',T16,3F12.5)
103 format(T3,'插值结果：',T16,3F12.6)
104 format(T3,'实际结果:',T16,3F12.6)
end program main