module doolittle

!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
!  Description :  LU 分解之doolittle模块
!    
!-----------------------------------------------------

contains  

subroutine solve(A,L,U,N)
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

end subroutine solve

end module doolittle



program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-8
!-----------------------------------------------------
!  Purpose   :   Doolittle 分解
!    
!-----------------------------------------------------
!  In put data  files :
!       1.    A,N
!       2.
!  Output data files  :
!       1.    L,U
!       2.
!-----------------------------------------------------
use doolittle

integer,parameter::N=3

real*8::A(n,n),L(N,N),U(N,N)

open(unit=11,file='fin.txt')
open(unit=12,file='fout.txt')

read(11,*)
read(11,*)((A(i,j),j=1,N),i=1,N)

call solve(A,L,U,N)

write(12,21)
21 format(T10,'Doolittle分解',/)


!输出L矩阵
write(12,*)'L='
do i=1,N
 write(12,22)L(i,:)
end do
22 format(3F10.6)

!输出U矩阵
write(12,*)'U='
do i=1,N
 write(12,22)U(i,:)
end do
23 format(3F10.6)

end program main 