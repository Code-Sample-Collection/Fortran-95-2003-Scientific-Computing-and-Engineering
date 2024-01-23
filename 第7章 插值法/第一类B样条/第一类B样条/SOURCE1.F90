module spline
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
!  Description :   第一类标准B样条模块
!    
!-----------------------------------------------------
!  Parameters  :
!      1.
!      2.
!----------------------------------------------------- 
!  Contains    :
!      1.   方法函数
!      2.   3阶样条基
!      3.   追赶法解三对角方程
!-----------------------------------------------------
!  Post Script :
!      1.
!      2. 
!-----------------------------------------------------

contains


subroutine solve(n,x,y,da,db,m,t,ty)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  第一类标准B样条 方法函数
!               可以直接对多点插值
!-----------------------------------------------------
!  Input  parameters  :
!       1.n 节点个数 减1 ,比如九个节点N=8
!       2. x 节点向量   （0：n）维
!       3. y 节点处函数值 (0:n)维
!       4. da  起点处导数值
!       5. db  终点处导数值
!       6.  m  要计算点的个数
!       7.  t  要计算的点   （1：m）维向量        
!  Output parameters  :
!       1.  ty  计算结果  （1：m）维向量      
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.
!       2.
!----------------------------------------------------

implicit real*8(a-z)

integer::n,m,i,j,k

real*8::x(0:n),y(0:n)

real*8::t(m),ty(m)

real*8::A(0:n,0:n),c0(0:n),f(0:n),c(-1:n+1)

!设置A矩阵-------------------
!设初值
A=0
do i=0,n
   A(i,i)=4d0
end do

!主对角以下元素
do i=1,n-1
   A(i,i-1)=1d0
end do

A(n,n-1)=2d0


!主对角以上元素
do i=1,n-1

A(i,i+1)=1d0
end do

A(0,1)=2d0

!-------------------------


!设置f 向量

h=x(1)-x(0)  ! 步长间隔

do i=1,n-1

f(i)=6d0*y(i)
end do

f(0)=6d0*y(0)+2d0*h*da
f(n)=6d0*y(n)-2d0*h*db


!注意调用追赶法时，方程维数的值
call chase(A,f,c0,N+1)




do i=0,n
  c(i)=c0(i)
end do


c(-1)=c0(1)-2d0*h*da

c(n+1)=c0(n-1)+2d0*h*db

!至此已经给出完整的系数c


!以下开始计算值

do k=1,m

ty(k)=0

do j=-1,n+1
  
   ty(k)=ty(k)+c(j)*omiga3((t(k)-x(0))/h-j)

end do


end do

end subroutine solve



function omiga3(x)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  三阶B样条基
!    
!----------------------------------------------------
!  Post Script :
!       1.
!       2.
!----------------------------------------------------
implicit real*8(a-z)

if (dabs(x)>=2d0) then
 omiga3=0
else if (dabs(x)<=1d0) then
 omiga3=0.5d0*(dabs(x))**3-x**2+2d0/3d0
else if (dabs(x)>1d0.and.dabs(x)<2d0) then
 omiga3=-1d0/6d0*(dabs(x))**3+x**2-2d0*dabs(x)+4d0/3d0
end if

end function omiga3


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



program  main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.05.13
!-----------------------------------------------------
!  Purpose   :  第一类标准B样条主函数
!    
!-----------------------------------------------------
!  In put data  files :
!       1.
!       2.
!  Output data files  :
!       1.
!       2.    result.txt保存计算结果
!-----------------------------------------------------
!  Post Script :
!       1.
!       2.
!-----------------------------------------------------

use spline


implicit real*8(a-z)

real*8::x(0:10),y(0:10),t(4),ty(4),rty(4)


 open(unit=11,file='result.txt')

!节点
x=(/     1.0000d0,&
    1.2000d0,&
    1.4000d0,&
    1.6000d0,&
    1.8000d0,&
    2.0000d0,&
    2.2000d0,&
    2.4000d0,&
    2.6000d0,&
    2.8000d0,&
    3.0000d0/)

! 节点因变量    
y=(/    0.841470984807897d0,&
   0.932039085967226d0,&
   0.985449729988460d0,&
   0.999573603041505d0,&
   0.973847630878195d0,&
   0.909297426825682d0,&
   0.808496403819590d0,&
   0.675463180551151d0,&
   0.515501371821464d0,&
   0.334988150155905d0,&
   0.141120008059867d0/)

!边界导数值   这里实为dcos(1d0)、dcos(3d0)
   da=0.540302305868140d0
   db=-0.989992496600445d0  
   
 ! 要计算的点 ,不必按照大小排列
   t=(/2.1d0,1.7d0,1.9d0,2.5d0/)
  
! 调用方法函数
   call solve(10,x,y,da,db,4,t,ty)
   
  
  !真值
  rty=dsin(t)
  
  write(11,101)
  
 write(11,102)((t(i),ty(i),rty(i)),i=1,4)
   

101 format(/,T15,'第一类标准B样条',//,&
              '     节点       插值结果     真值',/)
102 format(3F12.6)

end program main