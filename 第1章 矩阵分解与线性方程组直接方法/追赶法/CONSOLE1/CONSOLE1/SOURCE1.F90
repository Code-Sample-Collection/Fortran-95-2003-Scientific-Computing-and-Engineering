module chase
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-4-9
!-----------------------------------------------------
!  Description :  三对角线方程组方法模块
!    
!-----------------------------------------------------
!  Post Script :
!      1.
!      2. 
!-----------------------------------------------------

contains

 subroutine solve(A,f,x,N)
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

!---------------把A矩阵复制给向量 e,b,c 
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
 
end subroutine solve

end module chase


program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  追赶法计算三对角方程
!    
!-----------------------------------------------------
!  In put data  files :
!       1.   fin.txt输入数据
!       2.   fout.txt输出数据
!  Output data files  :
!       1.
!       2.
!-----------------------------------------------------

use chase

integer,parameter::N=4

real*8::A(N,N),f(N),x(N)

open(unit=11,file='fin.txt')
open(unit=12,file='fout.txt')

read(11,*)
read(11,*)((A(i,j),j=1,N),i=1,N)
read(11,*)f

call solve(A,f,x,N)

write(12,101)x
101 format(T5,'追赶法计算结果',/,T4,'x=',4(/F12.8))

end program main