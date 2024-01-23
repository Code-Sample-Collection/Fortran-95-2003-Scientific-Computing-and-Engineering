module crout

!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
!  Description :  LU 分解模块
!    
!-----------------------------------------------------

contains  

subroutine solve(A,L,U,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  LU之Crout分解
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
!L 第一列
L(:,1)=a(:,1)

!U 第一行
U(1,:)=a(1,:)/L(1,1)



do k=2,N

   do i=k,n
       s=0
       do r=1,k-1
        s=s+l(i,r)*u(r,k)
       end do
       l(i,k)=a(i,k)-s
   end do
   
   
   do j=k+1,n
     s=0
     do r=1,k-1
      s=s+l(k,r)*u(r,j)
     end do
     u(k,j)=(a(k,j)-s)/l(k,k)
       
   end do
   u(k,k)=1
   
end do

end subroutine solve

end module crout



program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-8
!-----------------------------------------------------
!  Purpose   :   Crot 分解
!    
!-----------------------------------------------------
!  In put data  files :
!       1.    A,N
!       2.
!  Output data files  :
!       1.    L,U
!       2.
!-----------------------------------------------------
use crout

integer,parameter::N=4

real*8::A(n,n),L(N,N),U(N,N)

open(unit=11,file='fin.txt')
open(unit=12,file='fout.txt')

read(11,*)
read(11,*)((A(i,j),j=1,N),i=1,N)

call solve(A,L,U,N)

write(12,21)
21 format(T10,'LU之Crout分解',/)


!输出L矩阵
write(12,*)'L='
do i=1,N
 write(12,22)L(i,:)
end do
22 format(4F10.6)

!输出U矩阵
write(12,*)'U='
do i=1,N
 write(12,22)U(i,:)
end do
23 format(4F10.6)

end program main 