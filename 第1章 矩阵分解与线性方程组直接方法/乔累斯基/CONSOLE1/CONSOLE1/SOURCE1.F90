module cholesky
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
!  Description :   Cholesky分解模块 
!    
!-----------------------------------------------------

contains 

subroutine  solve(A,L,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  Cholesky分解子程序
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.  A对称正定矩阵
!       2.  N矩阵阶数
!  Output parameters  :
!       1.  L 输出矩阵   A=L*L'
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.  Cholesky分解只适用于对称正定矩阵
!       2.
!----------------------------------------------------
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

end module cholesky


program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-8
!-----------------------------------------------------
!  Purpose   :
!              Cholesky
!-----------------------------------------------------
!  Output data files  :
!       1.   fout.txt 输出文件
!       2.
!-----------------------------------------------------
!  Post Script :
!       1.
!       2.  注意：Cholesky分解只适用于对称正定矩阵
!-----------------------------------------------------

use cholesky

integer,parameter::N=4

real*8::A(n,n),L(N,N)

open(unit=11,file='fin.txt')
open(unit=12,file='fout.txt')

read(11,*)
read(11,*)((A(i,j),j=1,N),i=1,N)

call solve(A,L,N)

write(12,21)
21 format(T5,'Cholesky分解之L矩阵',/)
do i=1,N
 write(12,22)L(i,:)
end do
22 format(4F10.6)

write(12,23)
23 format(/,1x,'其中  A=L*L’',/)

end program main 