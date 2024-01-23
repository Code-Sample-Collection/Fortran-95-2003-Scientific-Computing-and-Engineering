module sym_p
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
!  Description :   Cholesky分解计算对称正定方程模块
!    
!-----------------------------------------------------

contains 

subroutine solve(A,b,x,N)
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

end subroutine solve

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

end subroutine


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

end module sym_p


program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-8
!-----------------------------------------------------
!  Purpose   :  对称正定方程组的计算（Cholesky分解方法）
!              
!-----------------------------------------------------
! data files  :
!       1.  fin.out  输入文件 
!       2.  fout.txt 输出文件
!       2.
!-----------------------------------------------------

use sym_p

integer,parameter::N=3

real*8::A(N,N)
real*8::b(N),x(N)
open(unit=11,file='fin.txt')
open(unit=12,file='fout.txt')

read(11,*)
read(11,*)((A(i,j),j=1,N),i=1,N)
read(11,*)b

call solve(A,b,x,N)

write(12,21)
21 format(T5,'对称正定方程组的解为',/)

 write(12,22)x

22 format(T3,3F10.6)

end program main 