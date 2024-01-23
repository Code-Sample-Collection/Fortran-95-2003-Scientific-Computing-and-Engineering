module inv_mat
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-4-8
!-----------------------------------------------------
!  Description : 计算逆矩阵
!    
!-----------------------------------------------------
!  Contains    :
!      1.   inv 计算逆矩阵
!-----------------------------------------------------

contains



subroutine inv(A,invA,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  计算逆矩阵
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.
!       2. 
!  Output parameters  :
!       1.
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.
!       2.
!----------------------------------------------------
implicit real*8(a-z)

integer::n
integer::i
real*8::A(n,n),invA(n,n),E(n,n)

E=0

!设置E为单位矩阵
do i=1,n
   E(i,i)=1
end do

call mateq(A,E,invA,N,N)


end subroutine inv


subroutine mateq(A,B,X,N,M)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-8
!-----------------------------------------------------
!  Purpose   :  高斯列主元消去法计算矩阵方程
!                 AX=B
!-----------------------------------------------------
!  Input  parameters  :
!       1.   A(N,N)系数矩阵
!       2.   B(N,M)右矩阵
!       3.   N方程维数， 
!       4.   M---B的列数
!  Output parameters  :
!       1.  X  方程的根（N,M）维
!       2.
!  Common parameters  :
!
!----------------------------------------------------

implicit real*8(a-z)
integer::N,M,i

real*8::A(N,N),B(N,M),X(N,M)

real*8::btemp(N),xtemp(N)

do i=1,M
    
    btemp=B(:,i)
    call elgauss(A,btemp,xtemp,N)
   
    X(:,i)=xtemp
end do

end subroutine mateq


subroutine elgauss(A,b,x,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-8
!-----------------------------------------------------
!  Purpose   :  高斯列主元消去法
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

integer::i,k,N
integer::id_max  !主元素标号

real*8::A(N,N),b(N),x(N)

real*8::Aup(N,N),bup(N)

!Ab为增广矩阵  [Ab]
real*8::Ab(N,N+1)

real*8::vtemp1(N+1),vtemp2(N+1)

Ab(1:N,1:N)=A

Ab(:,N+1)=b


!##########################################################
!  这段是 列主元消去法的核心部分
do k=1,N-1

    elmax=dabs(Ab(k,k))
    id_max=k
    
    !这段为查找主元素	
    !这段程序的主要目的不是为了赋值最大元素给elmax，而是为了找出最大元素对应的标号

	
	do i=k+1,n
      if (dabs(Ab(i,k))>elmax) then
         elmax=Ab(i,k)

         id_max=i
      end if          
    end do

    
 !至此，已经完成查找最大元素，查找完成以后与  第k行交换 
 !交换两行元素，其他不变
    vtemp1=Ab(k,:)
    vtemp2=Ab(id_max,:)
   
    
    Ab(k,:)=vtemp2
    Ab(id_max,:)=vtemp1   
!
!以上一大段是为交换两行元素，交换完成以后即按照消元法进行
!#########################################################
  
   do i=k+1,N
  
     temp=Ab(i,k)/Ab(k,k)
     
     Ab(i,:)=Ab(i,:)-temp*Ab(k,:)
   
   end do

end do

!-----------------------------
! 经过上一步，Ab已经化为上三角矩阵
!
Aup(:,:)=Ab(1:N,1:N)

bup(:)=Ab(:,N+1)

!调用用上三角方程组的回带方法
call uptri(Aup,bup,x,n)

end subroutine elgauss



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

end module inv_mat

module  simple_newton
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
!  Description :   要计算的方程相关模块
!    
!----------------------------------------------------- 
!  Contains    :
!      1.    函数文件
!      2.    偏导数文件
!      3.    solve函数
!-----------------------------------------------------
!  Post Script :
!      1.
!      2. 
!-----------------------------------------------------

contains

subroutine  solve()
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  简化牛顿法计算非线性方程组
!    
!-----------------------------------------------------
!  In put data  files :
!       1.
!       2.
!  Output data files  :
!       1.    report.txt 计算结果迭代序列
!       2.
!-----------------------------------------------------
!  Post Script :
!       1.    需要用到线性方程组的解法，这里用选主元消去法
!       2.    需要准备函数文件与偏导数文件
!-----------------------------------------------------

use inv_mat


implicit real*8(a-z)

integer::I,itmax=100
integer::N=2

real*8::x(2),f(2),dx(2)
real*8::df(2,2),invdf(2,2)

!itmax 最大允许迭代次数
!N 方程组维数
!df 偏导数矩阵
!invdf 偏导数矩阵的逆矩阵

open(unit=11,file='result.txt')
write(11,101)
101 format(/,T6,'简化牛顿法计算非线性方程组迭代序列',/)

x=(/2d0,2d0/)

tol=1d-8

!计算初值的偏导数矩阵
call jac(df,x)

!计算偏导数的逆矩阵
!仅进行一次矩阵求逆，不进入循环
call  inv(df,invdf,N)

do i=1,itmax

  call  func(f,x)
  
   dx=-matmul(invdf,f) 
 
  x=x+dx
  
    
  write(11,102)i,x
102 format(I5,2F16.10)
  
!判断计算精度，当满足误差容限时 退出循环。
    dx2=dsqrt(dx(1)**2+dx(2)**2)
	
	if (dx2<tol) exit 
!------

end do

end subroutine solve

subroutine func(f,x)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  方程函数
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.    x 自变量
!       2. 
!  Output parameters  :
!       1.    f 方程函数
!       2.
!  Common parameters  :
!
!----------------------------------------------------

implicit real*8(a-z)
real*8::x(2),f(2)

f(1)=6*x(1)**3+x(1)*x(2)-3*x(2)**3-4
f(2)=x(1)**2-18*x(1)*x(2)**2+16*x(2)**3+1

end subroutine func



subroutine  jac(df,x)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  偏导数矩阵
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.
!       2. 
!  Output parameters  :
!       1.
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.
!       2.
!----------------------------------------------------
implicit real*8(a-z)
real*8::x(2),df(2,2)

df(1,1)=18*x(1)**2+x(2)
df(2,1)=2*x(1)-18*x(2)**2

df(1,2)=x(1)-9*x(2)**2
df(2,2)=-36*x(1)*x(2)+48*x(2)**2
end subroutine jac


end module simple_newton



program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.08.04  
!-----------------------------------------------------
!  Purpose     :    简化牛顿法主函数
!
!  Post Script :
!       1.
!       2.
!    
!-----------------------------------------------------

use simple_newton

!调用方法函数
call solve

end program main
