module LU

!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
!  Description :  LU 分解解方程
!    
!-----------------------------------------------------

contains  


subroutine solve(A,b,x,N)

implicit real*8(a-z)
integer::N
real*8::A(N,N),b(N),x(N)

real*8::L(N,N),U(N,N)

real*8::y(N)

 call doolittle(A,L,U,N)
  
 call  downtri(L,b,y,N)
 
 call uptri(U,y,x,N)

end subroutine solve


subroutine doolittle(A,L,U,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  LU分解之Doolittle函数
!              A=LU
!-----------------------------------------------------
!  Input  parameters  :
!       1.    A  方阵
!       2.    N  阶数
!  Output parameters  :
!       1.   L
!       2.   U
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.
!       2.
!----------------------------------------------------

implicit real*8(a-z)
integer::N,i,k,r

real*8::A(N,N),L(N,N),U(N,N)
!U的第一行


U(1,:)=A(1,:)

!L的第一列
L(:,1)=a(:,1)/U(1,1)

do k=2,N
   
    l(k,k)=1
   
   do j=k,n
       s=0
       do m=1,k-1
        s=s+l(k,m)*u(m,j)
       end do
       u(k,j)=a(k,j)-s
   end do
   
   
   do i=k+1,n
     s=0
     do m=1,k-1
      s=s+l(i,m)*u(m,k)
     end do
     l(i,k)=(a(i,k)-s)/u(k,k)
       
   end do
 
   
end do

end subroutine doolittle


subroutine uptri(A,b,x,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-8
!-----------------------------------------------------
!  Purpose   :  上三角方程组的回带方法
!                 Ax=b
!-----------------------------------------------------
!  Input  parameters  :
!       1.   A(N,N)系数矩阵
!       2.   b(N)右向量
!       3.   N方程维数
!  Output parameters  :
!       1.  x  方程的根
!       2.
!  Common parameters  :
!
!----------------------------------------------------

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
!  Input  parameters  :
!       1.   A(N,N)系数矩阵
!       2.   b(N)右向量
!       3.   N方程维数
!  Output parameters  :
!       1.  x  方程的根
!       2.
!  Common parameters  :
!
!----------------------------------------------------

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

end module LU



program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-8
!-----------------------------------------------------
!  Purpose   :   LU分解计算线性方程组
!              Ax=b
!-----------------------------------------------------
!  In put data  files :
!       1.    A,b
!       2.
!  Output data files  :
!       1.    x
!       2.
!-----------------------------------------------------
use LU

integer,parameter::N=4

real*8::A(n,n),L(N,N),U(N,N)
real*8::b(N),x(N)

open(unit=11,file='fin.txt')
open(unit=12,file='fout.txt')

read(11,*)
read(11,*)((A(i,j),j=1,N),i=1,N)
read(11,*)b

call solve(A,b,x,N)

write(12,101)x

101 format(T5,'LU分解计算线性方程组计算结果',//,4(/,F10.6))


end program main
