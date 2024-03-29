module spline
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
!  Description :   第三类标准B样条模块
!                  处理周期函数插值
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


subroutine solve(n,x,y,m,t,ty)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.05.13
!-----------------------------------------------------
!  Purpose   :  第三类标准B样条 方法函数
!               可以直接对多点插值
!               处理周期函数插值
!-----------------------------------------------------
!  Input  parameters  :
!       1.n 节点个数 减1 ,比如九个节点N=8
!       2. x 节点向量   （0：n）维
!       3. y 节点处函数值 (0:n)维
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

real*8::A(1:n,1:n),c0(1:n),f(1:n),c(-1:n+1)

!设置A矩阵-------注意与A矩阵下标
!设初值
A=0
do i=1,n
   A(i,i)=4d0
end do

!主对角以下元素
do i=2,n
   A(i,i-1)=1d0
end do




!主对角以上元素
do i=1,n-1

A(i,i+1)=1d0
end do



!-------------------------


!设置f 向量

do i=1,n-1

f(i)=6d0*y(i)
end do

f(n)=6d0*y(0)


!注意调用追赶法时，方程维数的值
!由参数意义决定
call chase(A,f,c0,N)




do i=1,n
  c(i)=c0(i)
end do


c(n+1)=c(1)
c(0)=c(n)
c(-1)=c(n-1)


!至此已经给出完整的系数c


h=x(1)-x(0)

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
!  Purpose   :  第三类标准B样条主函数
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

real*8::x(0:40),y(0:40),t(5),ty(5),rty(5)

real*8,parameter::pi=3.141592653589793d0

 open(unit=11,file='result.txt')

!节点
x=(/ 0d0,  &
     2d0,  &
     4d0,  &
     6d0,  &
     8d0,  &
    10d0,  &
    12d0,  &
    14d0,  &
    16d0,  &
    18d0,  &
    20d0,  &
    22d0,  &
    24d0,  &
    26d0,  &
    28d0,  &
    30d0,  &
    32d0,  &
    34d0,  &
    36d0,  &
    38d0,  &
    40d0,  &
    42d0,  &
    44d0,  &
    46d0,  &
    48d0,  &
    50d0,  &
    52d0,  &
    54d0,  &
    56d0,  &
    58d0,  &
    60d0,  &
    62d0,  &
    64d0,  &
    66d0,  &
    68d0,  &
    70d0,  &
    72d0,  &
    74d0,  &
    76d0,  &
    78d0,  &
    80d0  /)

! 节点因变量    
y=(/   1.000000000000000d0 ,&
   1.144122805635369d0 ,&
   1.260073510670101d0 ,&
   1.344997023927915d0 ,&
   1.396802246667421d0 ,&
   1.414213562373095d0 ,&
   1.396802246667421d0 ,&
   1.344997023927915d0 ,&
   1.260073510670101d0 ,&
   1.144122805635369d0 ,&
   1.000000000000000d0 ,&
   0.831253875554907d0 ,&
   0.642039521920206d0 ,&
   0.437016024448821d0 ,&
   0.221231742082474d0 ,&
   0.000000000000000d0 ,&
  -0.221231742082474d0 ,&
  -0.437016024448821d0 ,&
  -0.642039521920206d0 ,&
  -0.831253875554907d0 ,&
  -1.000000000000000d0 ,&
  -1.144122805635369d0 ,&
  -1.260073510670101d0 ,&
  -1.344997023927915d0 ,&
  -1.396802246667420d0 ,&
  -1.414213562373095d0 ,&
  -1.396802246667421d0 ,&
  -1.344997023927915d0 ,&
  -1.260073510670101d0 ,&
  -1.144122805635369d0 ,&
  -1.000000000000000d0 ,&
  -0.831253875554907d0 ,&
  -0.642039521920206d0 ,&
  -0.437016024448821d0 ,&
  -0.221231742082475d0 ,&
  -0.000000000000000d0 ,&
   0.221231742082474d0 ,&
   0.437016024448821d0 ,&
   0.642039521920206d0 ,&
   0.831253875554907d0 ,&
   1.000000000000000d0/)

   
 ! 要计算的点 ,不必按照大小排列
   t=(/35.2d0,27.6d0,19d0,65.7d0,55.7d0/)
  
! 调用方法函数
   call solve(40,x,y,5,t,ty)
   
  
  !真值
  rty=dsin(t*pi/40)+dcos(t*pi/40)
  
  write(11,101)
  
 write(11,102)((t(i),ty(i),rty(i)),i=1,5)
   

101 format(/,T15,'第三类标准B样条',//,&
              '     节点       插值结果     真值',/)
102 format(3F12.6)

end program main