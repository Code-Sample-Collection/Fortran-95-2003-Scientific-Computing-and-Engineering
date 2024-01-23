module iterprove
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-4-8
!-----------------------------------------------------
!  Description : 迭代改进模块
!    
!-----------------------------------------------------
!  Contains    :
!      1.   solve  方法函数
!      2.
!-----------------------------------------------------

contains


subroutine solve(A,b,x,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.08.10
!-----------------------------------------------------
!  Purpose     :  线性方程组的迭代改进
! 
!  Post Script :
!       1.    注意是先调用选主元消去法计算
!       2.    然后计算残差方程
!       3.    再迭代改进
!-----------------------------------------------------
!  Input  parameters  :
!       1.       A 系数矩阵
!       2.       b 右向量
!       3.       N 方程维数
!  Output parameters  :
!       1.       x 方程的解
!       2.
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------


implicit real*8(a-z)
integer::i,N


integer::itmax=5 
!默认迭代5次

real*8::A(n,n),b(n),x(N),x1(n),x2(n)

real*8::db(N),dx(N)

!先调用一次计算，得到初值
call gauss(A,b,x1,N)
   
!*****************************************   
!x1=1d0
! 该语句可以测试在X1非常差的情况下计算结果
! 如果要查看可以去掉这 x1=1d0的注释，
! 而把call gauss(A,b,x1,N)注释掉
!*****************************************

!输出初值
write(12,101)0,x1

do i=1,itmax

  !db为残差
  db=matmul(A,x1)-b

  call gauss(A,db,dx,N)
  
  !由此得到改进量dx
  
  
  !x2为改进后的新解
  x2=x1-dx
  
  !更新变量
  x1=x2
  
  !输出改进结果
  write(12,102)i,x1
  
end do

x=x1

101 format(/,T25,'迭代改进序列',//,I3,<n>(F10.5))
102 format(I3,<N>(F10.5))
end subroutine solve


subroutine gauss(A,b,x,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-8
!-----------------------------------------------------
!  Purpose   :  高斯列主元消去法
!                 Ax=b
!-----------------------------------------------------

implicit real*8(a-z)

integer::i,k,N
integer::id_max  !主元素标号

real*8::A(N,N),b(N),x(N)

real*8::Aup(N,N),bup(N)

!Ab为增广矩阵  [Ab]
real*8::Ab(N,N+1)

real*8::vtemp1(N+1),vtemp2(N+1)

Ab(1:N,1:N)=A

Ab(:,N+1)=b


do k=1,N-1

    elmax=dabs(Ab(k,k))
    id_max=k


	
	do i=k+1,n
      if (dabs(Ab(i,k))>elmax) then
         elmax=Ab(i,k)

         id_max=i
      end if          
    end do


    vtemp1=Ab(k,:)
    vtemp2=Ab(id_max,:)
   
    
    Ab(k,:)=vtemp2
    Ab(id_max,:)=vtemp1   
!

  
   do i=k+1,N
  
     temp=Ab(i,k)/Ab(k,k)
     
     Ab(i,:)=Ab(i,:)-temp*Ab(k,:)
   
   end do

end do

!-----------------------------
! 经过上一步，Ab已经化为如下形式的矩阵
!            | *  *  *  *  # |
!     [A b]= | 0  *  *  *  # |
!            | 0  0  *  *  # |
!            | 0  0  0  *  # |
!
Aup(:,:)=Ab(1:N,1:N)

bup(:)=Ab(:,N+1)

!调用用上三角方程组的回带方法
call uptri(Aup,bup,x,n)

end subroutine gauss



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

integer::i,j,N

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

end module iterprove



program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-8
!-----------------------------------------------------
!  Purpose   :  高斯列主元消去法
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
!-----------------------------------------------------
use iterprove

implicit real*8(a-z)

integer,parameter:: N=6

integer::i,j
real*8::A(N,N),b(N),x(N)



open(unit=11,file='fin.txt')
open(unit=12,file='fout.txt')

read(11,*)
!读入A矩阵

read(11,*)((A(i,j),j=1,N),i=1,N)
!读入B向量
read(11,*) b


!调用迭代模块
call solve(A,b,x,N)

end program main