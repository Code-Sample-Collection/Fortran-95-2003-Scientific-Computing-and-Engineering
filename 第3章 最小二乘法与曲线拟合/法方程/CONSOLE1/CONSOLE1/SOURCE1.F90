module  normal
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
!  Description :   超定方程的最小二乘计算模块
!                  
!-----------------------------------------------------
!  Post Script :
!      1.        模块主函数  solve
!      2. 
!-----------------------------------------------------

contains 

subroutine  solve(A,b,x,N,M)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  计算超定方程的最小二乘问题
!                  
!-----------------------------------------------------
!  Input  parameters  :
!       1.       A   (N,M)      N>M      
!       2.       b   (N)
!  Output parameters  :
!       1.       x  (M)
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.      方法是通过求解法方程
!       2.      而法方程的求解 基于Cholesky分解
!       3.    程序中 调用了2个基础函数
!            
!        transpose 矩阵转置，这个是非IMSL函数 ，Fortran自带的，可以直接使用 
!        Matmul    矩阵相乘，注意矩阵的维数
!----------------------------------------------------

implicit real*8(a-z)

integer::N,M,P
real*8::A(N,M),b(N),x(M)

real*8::AT(M,N)
real*8::ATA(M,M),ATb(M)


integer::i,j,k

! 转置  系统函数
AT=TRANSPOSE(A)
! 矩阵相乘 ，系统函数，注意维数匹配
ATA=MATMUL(AT,A)

ATb=MATMUL(AT,b)

!调用Cholesky分解方法计算法方程
call chol_eq(ATA,ATb,x,M)

end subroutine solve



subroutine chol_eq(A,b,x,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  用cholesky分解方法解方程
!    
!-----------------------------------------------------

implicit real*8(a-z)
integer::N
real*8::A(N,N),b(N),x(N)
real*8::L(N,N),y(N),LT(N,N)
!LT 为L的转置矩阵
integer::i,j

call  chol(A,L,N)

call  downtri(L,b,y,N)

do i=1,N
 do j=1,N
    LT(i,j)=L(j,i)
 end do
end do

call uptri(LT,y,x,N)  !这一步已经算出了x

end subroutine chol_eq

subroutine  chol(A,L,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  Cholesky分解子程序
!-----------------------------------------------------  
integer::N
real*8::A(N,N),L(N,N)
integer::i,j,k

L=0

L(1,1)=dsqrt(a(1,1))


L(2:,1)=a(2:,1)/L(1,1)


do j=2,N

   s=0
   do k=1,j-1
   s=s+L(j,k)**2
   end do
   
   L(j,j)=dsqrt(a(j,j)-s)
   
   !注意i范围
   do i=j+1,N
       
       s=0
       do k=1,j-1
         s=s+L(i,k)*L(j,k)
       end do
       
       L(i,j)=(a(i,j)-s)/L(j,j)      
     
   end do   

end do 

end subroutine chol


subroutine uptri(A,b,x,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-8
!-----------------------------------------------------
!  Purpose   :  上三角方程组的回带方法
!                 Ax=b
!-----------------------------------------------------

implicit real*8(a-z)

integer::i,j,k,N

real*8::A(N,N),b(N),x(N)

x(N)=b(N)/A(N,N)

!回带部分
do i=n-1,1,-1
   
    x(i)=b(i)
   do j=i+1,N
    x(i)=x(i)-a(i,j)*x(j)
   end do
    x(i)=x(i)/A(i,i)

end do

end subroutine uptri


subroutine downtri(A,b,x,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-9
!-----------------------------------------------------
!  Purpose   :  下三角方程组的回带方法
!                 Ax=b
!-----------------------------------------------------

implicit real*8(a-z)
integer::i,j,N
real*8::A(N,N),b(N),x(N)

x(1)=b(1)/a(1,1)

do k=2,N
   x(k)=b(k)
   do i=1,k-1
      x(k)=x(k)-a(k,i)*x(i)
   end do
   x(k)=x(k)/a(k,k)

end do

end subroutine downtri


end module normal




module driver
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
!  Description :   驱动程序模块
!    
!-----------------------------------------------------
!  Contains    :
!      1.    dri_main  读文件，并调用方法函数
!      2.
!-----------------------------------------------------
!  Post Script :
!      1.
!      2. 
!-----------------------------------------------------

contains


subroutine dri_main(N,M)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  驱动程序模块
!    
!----------------------------------------------------
!  Post Script :
!       1.      可以计算任意 N,M阶最小二乘  N>M  由用户输入
!       2.
!----------------------------------------------------

!引用方法模块
use normal

implicit real*8(a-z)

integer::N,M,i

real*8::A(N,M),b(N),x(M)

!读入系数A矩阵
read(11,*)((A(i,j),j=1,M),i=1,N)
!读入向量b
read(11,*)b

!调用方法函数
 call solve(A,b,x,N,M)
     
  write(12,101)
101 format(T5,'最小二乘解为：',/)

do i=1,m
 write(12,'(T5,F10.6)')x(i)
end do
     

end subroutine dri_main

end  module driver


program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-10
!-----------------------------------------------------
!  Purpose   :  生成计算最小二乘的软件，方法是采用解法方程
!                法方程的计算是基于cholesky分解
!-----------------------------------------------------
!  Post Script :
!       1.       软件适用于  超定方程，系数由文件读入
!       2.       由主函数启动驱动程序模块，驱动程序调用方法函数
!-----------------------------------------------------

use driver
integer::N,M

open(unit=11,file='fin.txt')
open(unit=12,file='fout.txt')


read(11,*)
!读入A矩阵

!读入方程维数系数
read(11,*)N,M

!调用驱动函数
call dri_main(N,M)

end program main 