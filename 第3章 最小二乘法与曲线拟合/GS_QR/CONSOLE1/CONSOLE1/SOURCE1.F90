module gram_sch
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010.08.10
!-----------------------------------------------------
!  Description :   修正的Gram-Schimdt正交化求QR分解
!    
!-----------------------------------------------------


contains 


subroutine solve(A,Q,R,M,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :   采用修正的 Gram-Schmidt分解求矩阵的QR分解
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.    A原始矩阵
!       2.    A(M,N)
!  Output parameters  :
!       1.    分解结果为   Q(M,N):注意 Q不是方阵，Q列向量为标准正交基
!       2.                 R(N,N)：R是方阵
!       3.   
!----------------------------------------------------
!  Post Script :
!       1.  注意矩阵的维数，分解后Q列向量是正交的
!       2.  关于编程方法可以参看《矩阵分析与应用》张贤达编著
!       3.  详细的数学解释，可以参看 麻省理工学院的
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
 
end subroutine solve


subroutine dri_main(M,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  驱动程序
!    
!-----------------------------------------------------

integer::M,N
real*8::A(m,n),Q(m,n),R(n,n)
real*8::b(M),x(N)

!读入矩阵
read(11,*)((A(i,j),j=1,n),i=1,m)




call solve(A,Q,R,m,n)

!-------------------这段程序用于输出QR分解
write(12,101)
101 format(T10,'修正的Gram-Schmidt方法QR分解',/,T3,'Q=',/)

write(12,102)((Q(i,j),j=1,N),i=1,M)

!变量输出格式只针对 IVF编译器，在CVF中不支持
102 format(<M>(<N>F10.5/))

write(12,103)
103 format(T3,'R=',/)

write(12,104)((R(i,j),j=1,n),i=1,n)

!变量输出格式只针对 IVF编译器，在CVF中不支持
104 format(<n>(<N>F10.5/))

!------------------------------

end subroutine dri_main

end module gram_sch



program main

!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-11
!-----------------------------------------------------
!  Purpose   :  1. 采用修正的Gram-Schimidt方法进行QR分解
!               2. 的主函数
!-----------------------------------------------------
!  In put data  files :
!       1.         fin.txt 输入文件
!       2.   
!  Output data files  :
!       1.         fou.txt 输出文件
!       2.
!-----------------------------------------------------
!  Post Script :
!       1.         由主函数引导驱动程序
!       2.       
!-----------------------------------------------------
use gram_sch
integer::M,N

open(unit=11,file='fin.txt')
open(unit=12,file='fout.txt')


read(11,*)
!读入A矩阵

!读入方程维数系数
read(11,*)M,N

!调用驱动函数
call dri_main(M,N)

end program main 