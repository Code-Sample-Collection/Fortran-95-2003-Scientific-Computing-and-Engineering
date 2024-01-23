module tri_eq
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-4-8
!-----------------------------------------------------
!  Description : 用于解上、下三角形线性方程组的回带方法模块
!    
!-----------------------------------------------------
!  Contains    :
!      1.    
!      2.
!-----------------------------------------------------

contains


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


end module tri_eq

!##############################################################
module driver
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
!  Description :  驱动函数模块
!    
!-----------------------------------------------------
!  Parameters  :
!      1.
!      2. 
!  Contains    :
!      1.       dir_main  驱动函数入口
!      2.       dri_up    当读到关键字uptri 时启动该函数
!      3.       dri_down  当读到关键字 downtri时启动该函数
!-----------------------------------------------------
!  Post Script :
!      1.
!      2. 
!-----------------------------------------------------

contains 

subroutine dri_main()
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-9
!-----------------------------------------------------
!  Purpose   :  驱动程序入口函数
!    
!-----------------------------------------------------
!  Input  files  :
!       1.   fin.txt  准备数据
!       2.   
!  Output files  :
!       1.  fout.txt 输出结果文件
!       2. 
!
!----------------------------------------------------

implicit real*8(a-z)
integer::ioerr
character(20)::upordown

open(unit=11,file='fin.txt',iostat=ioerr)
open(unit=12,file='fout.txt')

do 
  read(11,*)upordown
  !  读输入文件
  !  当读到关键字uptri时启动上三角矩阵计算
  !  当读到关键字downtri时候启动下三角矩阵计算
  
  if ( upordown(1:5) == 'uptri' ) then
    call dri_up()
  else if (upordown(1:)=='downtri')then
    call dri_down()
  end if
  
  if (ioerr/=0) exit
  !读到文件结束时，退出读文件
end do

end subroutine dri_main


subroutine dri_up()
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-9
!-----------------------------------------------------
!  Purpose   : 启动上三角阵的计算
!    
!-----------------------------------------------------
use tri_eq
implicit real*8(a-z)

integer,parameter::N=4

integer::i,j
real*8::A(N,N),b(N),x(N)


read(11,*)((A(i,j),j=1,N),i=1,N)
!读入B向量
read(11,*) b

call uptri(A,b,x,N)

write(12,101)x
101 format(T5,'上三角形方程组的解',/,T4,'x=',4(/F12.8))

end subroutine dri_up


subroutine dri_down()
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-9
!-----------------------------------------------------
!  Purpose   :  启动下三角阵的计算
!    
!-----------------------------------------------------

use tri_eq

implicit real*8(a-z)

integer,parameter::N=4

integer::i,j
real*8::A(N,N),b(N),x(N)


read(11,*)((A(i,j),j=1,N),i=1,N)
!读入B向量
read(11,*) b

call downtri(A,b,x,N)

write(12,101)x
101 format(/,T5,'下三角形方程组的解',/,T4,'x=',4(/F12.8))

end subroutine dri_down

end module driver



!-----------------------------------------------------
program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-9
!-----------------------------------------------------
!  Purpose   :  计算上、下三角形方程组
!    
!-----------------------------------------------------
!  In put data  files :
!       1.  fin.txt  输入方程系数
!       2.
!  Output data files  :
!       1. fout.txt  计算结果
!       2.
!-----------------------------------------------------
!  Post Script :
!       1.    需要准备输入数据
!
!       2.    由主函数启动驱动程序进行计算   
!-----------------------------------------------------
!use tri_eq
use driver

!调用驱动函数
call dri_main()

end program main