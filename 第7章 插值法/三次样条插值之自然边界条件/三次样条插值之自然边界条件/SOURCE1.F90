
module  spline

!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010.05.12
!-----------------------------------------------------
!  Description :   三次样条插值之第二类边界条件模块
!    
!----------------------------------------------------- 
!  Contains    :
!      1.     solve函数 即方法函数
!      2.    
!      3.     chase 为用追赶法计算线性方程，因为插值中的
!             方程为三对角方程，用追赶法计算效率较高。
!-----------------------------------------------------
!  Post Script :
!      1.
!      2.     可以直接对向量插值
!-----------------------------------------------------


contains


subroutine solve(n,x,y,d2fa,d2fb,nt,t,ty)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.05.12
!-----------------------------------------------------
!  Purpose   :   三次样条之第二类边界条件
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.   n-----插值节点个数减1，如有九个节点 则N=8
!       2.   x ---节点自变量  为（0：N）维向量
!       3.   y----节点因变量  （0：N）维向量
!       4.   nt 要计算向量的维数
!       5.   t 要计算的向量  (1:nt) 维向量
!       6.   
!       7.    d2fa,d2fb  起点于终点处的二阶导数条件，
!              如果是自然边界条件，则二者为皆为0
!  Output parameters  :
!       1.   ty ---要计算的向量，（1：nt）维
!
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!          
!       1.   自然边界条件是第一类边界条件的特殊情况
!            本程序可以直接处理第一类边界条件
!       2.     
!            对于要插值的向量分量，可以不必按大小排列
!----------------------------------------------------

implicit real*8(a-z)

integer::n,nt
!n为插值节点数减1，即如果有9个节点，则n=8
!nt 要插值点的个数，及向量t,ty的维数

integer::i,j,k

real*8::x(0:n),y(0:n)

real*8::t(nt),ty(nt)

real*8::h(0:n-1)

real*8::f1(0:n-1),f2(1:n-1)

real*8::u(1:n-1),namda(1:n-1),d(1:n-1)

real*8::M(0:n),v(1:n-1)


real*8::A(1:n-1,1:n-1)


M(0)=d2fa
M(n)=d2fb


do i=0,n-1
  h(i)=x(i+1)-x(i)
  f1(i)=(y(i+1)-y(i))/h(i)
end do


!求得 u, namda, d
do i=1,n-1

 u(i)=h(i-1)/(h(i-1)+h(i))
 namda(i)=1-u(i)
 
 f2(i)=(f1(i-1)-f1(i))/(x(i-1)-x(i+1))
  
  d(i)=6d0*f2(i) 
  
end do




!设置A矩阵值
A=0
do i=1,n-1
a(i,i)=2d0
end do

do i=2,n-1
a(i,i-1)=u(i)
end do

do i=1,n-2
a(i,i+1)=namda(i)
end do


! 设置右向量值
d(1)=d(1)-u(1)*M(0)
d(n-1)=d(n-1)-namda(n-1)*M(n)


call chase(a,d,v,N-1)

do i=1,n-1

   M(i)=v(i)

end do


!--------以上以及求得系数
!已经完成插值多项式的建立

!------------------------------------------------
! 以下开始计算具体值
do k=1,nt

!------------
!  对要插值向量每个分量而言，先找到其在数据中的位置
do i=1,n-1

 if (t(k)<x(i+1)) exit 

end do


  ty(k)=M(i)*(x(i+1)-t(k))**3/6d0/h(i)+ M(i+1)*(t(k)-x(i))**3/6d0/h(i)+ &

     (y(i)-M(i)*h(i)**2/6d0)*(x(i+1)-t(k))/h(i)+ &
     (y(i+1)-M(i+1)*h(i)**2/6d0)*(t(k)-x(i))/h(i)

end do

!-------------------------------




!为了方便读者比对结果，这里输出中间值，实际应用时可以去掉输出部分
!-------------------------
write(11,101)
write(11,102)((i,h(i),m(i)),i=0,n-1)   !输出h


write(11,103)
write(11,104)((i,u(i),namda(i),d(i)),i=1,n-1)

101 format(/,T5,'          h             m')
102 format(I3,2F16.8)

103 format(/,T5,'         u             namda            d')
104 format(I3,3F16.8)
!----------------------------------------------------------------

end subroutine solve



 subroutine chase(A,f,x,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-9
!-----------------------------------------------------
!  Purpose   :  追赶法计算三对角方程组
!              Ax=f
!-----------------------------------------------------
!  Input  parameters  :
!       1.  A系数矩阵
!       2. f 右向量
!  Output parameters  :
!       1.  x方程的解
!       2.  N维数
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.   注意：该方法仅适用于三对角方程组
!       2.
!---------------------------------------------------
 
 implicit real*8(a-z)
 integer::N
 real*8::A(N,N),f(N),x(N)
 real*8::L(2:N),u(N),d(1:N-1)
 
 real*8::c(1:N-1),b(N),e(2:N)

 integer::i
 
 real*8::y(N)

!---------------把A矩阵复制给向量e,b,c 
do i=1,N
    b(i)=a(i,i)
end do 

do i=1,N-1
   c(i)=a(i,i+1)
end do

do i=2,N
  e(i)=a(i,i-1)
end do
!------------------------

do i=1,N-1
 d(i)=c(i)
end do 


u(1)=b(1)

do i=2,N
  L(i)=e(i)/u(i-1)
  u(i)=b(i)-L(i)*c(i-1)
end do

!------开始回带,求得y

y(1)=f(1)
do i=2,N
  y(i)=f(i)-L(i)*y(i-1)
end do

!-----开始回带，求得x

x(n)=y(n)/u(n)

do i=n-1,1,-1
x(i)=(y(i)-c(i)*x(i+1))/u(i)
end do

end subroutine chase



end module spline







program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  三次样条插值之自然边界条件
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
use spline


implicit  real*8(a-z)

integer::i

real*8::x(0:8),y(0:8)

real*8::t(8),ty(8)

open(unit=11,file='result.txt')

!插值节点，九个节点
x=(/1d0,2d0,5d0,6d0,7d0,8d0,10d0,13d0,17d0/)
!插值因变量
y=(/3d0,3.7d0,3.9d0,4.2d0,5.7d0,6.6d0,7.1d0,6.7d0,4.5d0/)


!要计算的点，可以不必按照大小排列
t=(/3.5d0,1.5d0,5.5d0,6.5d0,7.5d0,9.0d0,11.5d0,15d0/)


write(11,101)


!
!调用函数
call solve(8,x,y,0d0,0d0,8,t,ty)
! 作者提供的函数具有一般性，可以直接处理第一类边界条件，
!而如果是自然条件, 则 边界处的 二阶导数为0


write(11,102)

write(11,103)((i,t(i),ty(i)),i=1,8)



101 format(/,T5,'三次样条插值之自然边界条件')

102 format(/,T5,'        t               ty')

103 format(I3,2F16.8)


end program main