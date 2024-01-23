module lsqcurvefit
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010.07.09
!-----------------------------------------------------
!  Description :   任意阶多项式曲线拟合模块
!  
!  Post Script :
!      1.
!      2. 
!
!-----------------------------------------------------
!  Contains    :
!      1.    拟合函数，基函数
!      2.    最小二乘法函数
!-----------------------------------------------------
!  Parameters  :
!      1.
!      2.
!----------------------------------------------------- 
contains 
subroutine solve(x,y,N,c,m)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose     :  任意阶多项式拟合函数
! 
!  Post Script :
!       1.
!       2.
!       3.    
!-----------------------------------------------------
!  Input  parameters  :
!       1.  x,y --- 输入数据
!       2.  N --- x,y向量的维数
!       3.  m----希望用m阶多项式拟合
!  Output parameters  :
!       1.
!       2.
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------
implicit real*8(a-z)

integer::N,m
real*8::x(n),y(n),c(M)

real*8::bv(M)

integer::i

!A系数矩阵
real*8::A(N,M)


!如果多项式阶数加1高于实测数据个数则报错
!-------------
if (m>n) then

write(11,101)
stop
end if

101 format(/,'warning:the order of polynomial+1 > the ',/,&
              'number of date,that is forbidden')
!------------------------------



do i=1,n
    
    call basefunc(bv,x(i),m)
    
    A(i,:)=bv

end do

call leas_eq(A,y,c,N,m)


end subroutine solve



subroutine basefunc(bv,x,m)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :2010.07.09  
!-----------------------------------------------------
!  Purpose     :  任意阶多项式基函数
! 
!  Post Script :
!       1.
!       2.
!       3.    
!-----------------------------------------------------
!  Input  parameters  :
!       1.    x 标量
!       2.    m 多项式阶数
!  Output parameters  :
!       1.
!       2.
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------
implicit real*8(a-z)

integer::M

real*8::bv(m)

integer::i

bv(1)=1d0

do  i=2,m
    bv(i)=bv(i-1)*x
end do
end subroutine basefunc

subroutine  leas_eq(A,b,x,M,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  通过修正的Gram-Schmidt正交化求最小二问题
!               方法函数
!-----------------------------------------------------
!  Method    :
!               对超定方程 A进行QR分解后  方程变为
!                   QR x=b
!                => Rx=Q'b   R为上三角阵
!                => 回带，可以求得最小二乘意义下的解
!-----------------------------------------------------
!  Post Script :
!       1.       即求解超定方程组 Ax=b    其中A(M,N)  M>N    
!       2.
!----------------------------------------------------
implicit real*8(a-z)

integer::M,N
real*8::A(M,N),Q(M,N),R(N,N)
real*8::b(M)
real*8::QT(N,M)  !Q的转置矩阵
real*8::QTb(N)   !Q'b
real*8::x(N)
 
call gram_dec(A,Q,R,M,N)

QT=transpose(Q)
QTb=matmul(QT,b)  !  Rx=Q'b

call uptri(R,QTb,x,N) !回带

end subroutine leas_eq

subroutine gram_dec(A,Q,R,M,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :   采用修正的Gram-Schmidt分解求矩阵的QR分解
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.    A原始矩阵
!       2.    A(M,N)
!  Output parameters  :
!       1.    分解结果为  Q(M,N):注意Q不是方阵，Q列向量为标准正交基
!       2.                 R(N,N)：R是方阵
!       3.   
!----------------------------------------------------
!  Post Script :
!       1.  注意矩阵的维数，分解后Q列向量是正交的
!       2.  关于编程方法可以参看《矩阵分析与应用》张贤达编著
!       3.  详细的数学解释，可以参看麻省理工学院的
!           线性代数教材《Linear Algebra with Application》
!----------------------------------------------------
implicit real*8(a-z)

integer::M,N
integer::i,j,k

real*8::A(M,N),Q(M,N),R(N,N)

real*8::vec_temp(M)



R(1,1)=dsqrt(dot_product(a(:,1),a(:,1)))

Q(:,1)=a(:,1)/R(1,1)


do k=2,N

      do j=1,k-1
        R(j,k)=dot_product(Q(:,j),A(:,k))   
      end do
   
      vec_temp=A(:,k)
   
      do j=1,k-1
   
        vec_temp=vec_temp-Q(:,j)*R(j,k)
   
      end do


    R(k,k)=dsqrt(dot_product(vec_temp,vec_temp))

     
    Q(:,k)=vec_temp/R(k,k)

end do
 
end subroutine gram_dec

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

end module lsqcurvefit



program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.07.09
!-----------------------------------------------------
!  Purpose     : 多项式拟合主函数
!
!  Post Script :
!       1.
!       2.
!    
!-----------------------------------------------------
!  In put data  files :
!       1.
!       2.
!  Output data files  :
!       1.
!       2.
!
!-----------------------------------------------------
use  lsqcurvefit

implicit real*8(a-z)

real*8::x(7),y(7)

real*8::c1(3),c2(4),c3(9)

open(unit=11,file='result.txt')

x=(/-3d0,-2d0,-1d0,0d0,1d0,2d0,3d0/)
y=(/4d0,2d0,3d0,0d0,-1d0,-2d0,-5d0/)

call solve(x,y,7,c1,3)

call solve(x,y,7,c2,4)



write(11,101)c1,c2

101 format(/T4,'多项式拟合曲线拟合',//,&
              T3,'采用二阶多项式拟合，系数为：',3(/,F16.11),//,&
              T3,'采用三阶多项式拟合，系数为：',4(/,F16.11),/)
              

call solve(x,y,7,c3,9)              
              

end program main