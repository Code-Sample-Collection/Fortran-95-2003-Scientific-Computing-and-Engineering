module det
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
!  Description :  计算矩阵行列式模块
!    
!-----------------------------------------------------

contains

subroutine  solve(A,d,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  计算矩阵行列式函数
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.  A 矩阵
!       2.  N 矩阵维数
!  Output parameters  :
!       1.    d  矩阵行列式
!       2.
!  Common parameters  :
!
!----------------------------------------------------

implicit real*8(a-z)

integer::N,i
real*8::A(N,N),L(N,N),U(N,N)

call crout(A,L,U,N)

d=1d0

do i=1,N   
   d=d*L(i,i)
end do

end subroutine solve


subroutine crout(A,L,U,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  LU之Crout分解
!              A=LU
!-----------------------------------------------------

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

end subroutine crout

end module det



module driver
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
!  Description : 
!                驱动程序模块
!-----------------------------------------------------

contains

subroutine dri_main(N)

use det
integer::N
real*8::A(N,N),d

!读入矩阵
read(11,*)((A(i,j),j=1,N),i=1,N)

call solve(A,d,N)

write(12,101)

101 format(T5,'该矩阵的行列式为',/)

write(12,*)d


end subroutine dri_main

end module driver



program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :   计算矩阵行列式
!    
!-----------------------------------------------------
!  In put data  files :
!       1.   fin.txt   读入的数据
!       2.
!  Output data files  :
!       1.
!       2.   fout.txt 计算结果，矩阵的维数由输入卡片读取
!-----------------------------------------------------
use driver
integer::N

open(unit=11,file='fin.txt')
open(unit=12,file='fout.txt')


read(11,*)
!读入A矩阵

!读入方程维数系数
read(11,*)N

!调用驱动函数
call dri_main(N)

end program main 