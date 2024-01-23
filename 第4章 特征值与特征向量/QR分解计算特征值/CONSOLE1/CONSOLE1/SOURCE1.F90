module eig_qr

!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
!  Description :   QR分解计算全部特征值模块
!  
!----------------------------------------------------- 
!  Contains    :
!      1.        方法函数
!      2.        QR分解函数（采用G-S正交化方法）
!-----------------------------------------------------

contains 

subroutine solve(A,N,namda,tol)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose     :  QR分解计算全部特征值
! 
!  Post Script :
!       1.       QR分解采用修正的G-S分解方法，当然也可以
!       2.        使用Householder或者Givens变换方法
!       3.    
!-----------------------------------------------------
!  Input  parameters  :
!       1.   A 需要计算的矩阵
!       2.   N 矩阵的维数
!  Output parameters  :
!       1.   namda  特征值组成的向量N维
!       2.   tol    用户指定的误差容限
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------


implicit real*8(a-z)

integer::N
real*8::A(N,N),namda(N)



real*8::A1(N,N),Q(N,N),R(N,N)


integer::i,j,k


A1=A

!循环迭代,最大允许迭代200次
do i=1,200

   call gram_dec(A1,Q,R,N,N)
   
   A1=matmul(R,Q)
    
! 判断迭代停止标准
   do k=1,N
    ds=0d0
    ds=ds+A1(k,k)**2    
   end do

   do j=1,N
   namda(j)=A1(j,j)
   end do 

    if (ds<tol) exit
    
end do

end subroutine solve


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


end module eig_qr




program main

!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.07.06
!-----------------------------------------------------
!  Purpose     :  QR分解计算全部特征值
!
!  Post Script :
!       1.        QR分解采用的修正的 G-S正交化方法
!       2.
!    
!-----------------------------------------------------
!  In put data  files :
!       1.
!       2.
!  Output data files  :
!       1.     result.txt  结果输出文件
!       2.
!
!-----------------------------------------------------

use eig_qr

implicit real*8(a-z)


real*8::A(3,3),namda(3)


open(unit=11,file='result.txt')


A=reshape((/2d0,1d0,0d0,&
           1d0,3d0,1d0,&
           0d0,1d0,4d0/),(/3,3/))

!调用方法函数
call solve(A,3,namda,1d-7)

write(11,101)namda

101 format(/,T4,'QR分解计算全部特征值',//,&  
           T3,'特征值为：',3(/F16.10))

end program main


