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


module  homotopy
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
!-----------------------------------------------------
!  Post Script :
!      1.
!      2. 
!-----------------------------------------------------

use inv_mat

contains


subroutine solve(x0,N,step,xt)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :   homotopy方法函数
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

integer::k,n,itmax=50

integer::step  !划分序列

real*8::x1(n),x0(n),y(n),f0(n),f1(n),dx(n)

real*8::x2(n),xt(n),vect1(n)
!xt 计算终值

real*8::H0(n,n),H1(n,n),df(n,n)
real*8::m1(n,n),m2(n,n)

!itmax 最大允许迭代次数


!计算初始函数值
call func(f0,x0)

x1=x0

write(11,101)

do k=1,step-1
    
    
  
    call func(f1,x1)  
    
    call jac(df,x1)
    
    call inv(df,H1,n)
    
    vect1=f1-(1d0-k*1d0/step*1d0*f0)
    
    x2=x1-matmul(H1,vect1)    
    
    
    !把计算迭代序列写入文件
    write(11,102)k,x2  

	x1=x2
    
  
end do

101 format(/,T12,'数值延拓法获取初值序列',/)

102 format(I4,3F16.10)

 xt=x2
end subroutine solve


function vdot(a,b,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   : 计算向量a,b内积
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.   N 向量维数
!       2. 
!----------------------------------------------------
!  Post Script :
!       1.   vdot=a1(1)*(1)+a(2)*b(2)+...
!       2.
!----------------------------------------------------
integer::i,N
real*8::a(n),b(n),vdot

vdot=0

do i=1,N  
  vdot=vdot+a(i)*b(i)
end do

end function vdot



function vvmat(a,b,n)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  N维列向量a 乘以  N维横向量b
!               结果为   n*n维 矩阵
!----------------------------------------------------
!  Post Script :
!       1.
!       2.
!----------------------------------------------------
integer::n,i,j
real*8::a(n),b(n),vvmat(n,n)

do i=1,n
  do j=1,n
    vvmat(i,j)=a(i)*b(j)
  end do
end do

end function vvmat


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
real*8::x(2),f(2),pi=3.141592653589793

f(1)=x(1)**2-x(2)+1d0

f(2)=x(1)-dcos(pi/2d0*x(2))

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
real*8::x(2),df(2,2),pi=3.141592653589793

df(1,1)=2d0*x(1)
df(1,2)=-1d0


df(2,1)=1d0
df(2,2)=dsin(pi/2d0*x(2))*pi/2d0


end subroutine jac


end module homotopy





subroutine newton(x0,N,itmax)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  采用牛顿法 获得最终计算结果
!    
!----------------------------------------------------
!  Post Script :
!       1.   x0为初值，N为阶数，tol为误差容限
!       2.   计算结果保存在文件中
!----------------------------------------------------


use  inv_mat
use  homotopy

implicit real*8(a-z)
integer::i,n,itmax

real*8::x2(N),x1(N),x0(N),dx(N),f(N)

real*8::df(N,N),H(n,n)

x1=x0

write(12,101)

do i=1,itmax
   
   call func(f,x1)
   call jac(df,x1)
   call inv(df,H,N)

   dx=-matmul(H,f)

   x2=x1+dx

   write(12,102)i,x2

!更新x值
   x1=x2


end  do

101 format(/,T5,'由数值延拓法提供初值后牛顿法计算序列',/)
102 format(I4,2F16.9)

end subroutine newton






program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      : 2010.05.11 
!-----------------------------------------------------
!  Purpose   :   数值延拓法计算非线性方程组主函数
!    
!-----------------------------------------------------
!  In put data  files :
!       1.  
!       2.
!  Output data files  :
!       1.   result.txt 给出了计算结果
!       2.   Hmatrix.txt给出了计算中的H矩阵值
!-----------------------------------------------------
!  Post Script :
!       1.
!       2.
!-----------------------------------------------------
use  inv_mat
use  homotopy

implicit real*8(a-z)
integer::N=2,step=8  
!N 为方程组的维数，step为微分方程的步数
real*8::x0(2),x1(2),xt(2)

open(unit=11,file='result.txt')
open(unit=12,file='newton.txt')

x0=(/1d0,1d0/)


!直接调用牛顿法，则计算失败不收敛,可以保留此语句
!而注释下面两行语句，进行验证
!--------
!call newton(x0,n,1d-8)
!--------

!step为微分方程步数
! 通过数值延拓法 计算的中间结果放在文件result.txt中
! 最后终值 放在参数xt中 回代出来，后面作为牛顿法初值
call solve(x0,n,step,xt)


!把数值延拓法计算结果作为牛顿法初值，最大允许迭代5次
!计算结果放在文件newton.txt中
call newton(xt,n,5)

end program main


