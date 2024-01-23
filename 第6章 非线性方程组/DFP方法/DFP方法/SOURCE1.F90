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


module  DFP
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


subroutine solve(x0,N,tol)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :   DFP方法函数
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

integer::i,n,itmax=50

real*8::x1(n),x0(n),y(n),f0(n),f1(n),dx(n)

real*8::v1(n),v2(n)

real*8::H0(n,n),H1(n,n),df(n,n)
real*8::m1(n,n),m2(n,n)

!itmax 最大允许迭代次数


!tol 误差容限

call jac(df,x0)

!注意：设置H0初值为  偏导数矩阵的逆矩阵
!  也可以直接设置H0为单位矩阵，但是那样收敛较慢
!  设置H0初值为偏导数逆矩阵，但是在迭代过程中则不计算逆矩阵
call inv(df,H0,N)


write(11,101)
write(12,102)

do i=1,itmax
    
    !计算函数值
    call func(f0,x0)
    
    !更新方程的根
    x1=x0-matmul(H0,f0)
    
    
    dx=x1-x0
    
    call func(f1,x1)
    
    y=f1-f0
    
    m1=vvmat(dx,dx,n)
    t1=vdot(dx,y,N)
    
    !第一项矩阵
    m1=m1/t1
    
    
    !第二项分子 结果为矩阵
    m2=vvmat(y,y,n)
    m2=matmul(H0,M2)
    M2=matmul(m2,H0)
    
    !第二项分母
    v1=matmul(H0,y)
    t2=vdot(y,v1,N)
    
    m2=m2/t2
    
    H1=H0+m1-m2
    
    !把计算中的H矩阵写入文件
    write(12,103)i,H0
    
    x0=x1
    H0=H1
    
    !把计算迭代序列写入文件
    write(11,104)i,x0  
    
      
    
    !判断计算精度，当满足误差容限时 退出循环。
    dx2=dsqrt(dx(1)**2+dx(2)**2)
	
	if (dx2<tol) exit 
!------
  
end do

101 format(/,T18,'DFP方法计算序列',/)
102 format(/,T5,'DFP方法H矩阵序列为：',/)

103 format(T5,'iter=',I4,/,<N>(<N>F16.10/),/)
104 format(I4,3F16.10)


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
real*8::x(3),f(3),pi=3.141592653589793

f(1)=3*x(1)-dcos(x(2)*x(3))-1d0/2

f(2)=x(1)**2-81*(x(2)+0.1d0)**2+dsin(x(3))+1.06d0

f(3)=dexp(-x(1)*x(2))+20*x(3)+(10*pi-3d0)/3

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
real*8::x(3),df(3,3)

df(1,1)=3d0
df(1,2)=x(3)*dsin(x(2)*x(3))
df(1,3)=x(2)*dsin(x(2)*x(3))

df(2,1)=2*x(1)
df(2,2)=-162*(x(1)+0.1)
df(2,3)=dcos(x(3))

df(3,1)=-x(2)*dexp(-x(1)*x(2))
df(3,2)=-x(1)*dexp(-x(1)*x(2))
df(3,3)=20d0

end subroutine jac


end module DFP




program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :   DFP方法计算非线性方程组主函数
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
use  DFP

implicit real*8(a-z)
integer::N=3
real*8::x0(3)

open(unit=11,file='result.txt')
open(unit=12,file='Hmatrix.txt')

x0=(/0.1d0,0.1d0,-0.1d0/)


!1d-8  表示允许的误差容限
call solve(x0,n,1d-8)

end program main


