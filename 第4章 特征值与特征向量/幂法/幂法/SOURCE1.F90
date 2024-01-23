module power
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-5-31
!-----------------------------------------------------
!  Description :   幂法模块
!    
!-----------------------------------------------------
!  Parameters  :
!      1.
!      2.
!----------------------------------------------------- 
!  Contains    :
!      1.   方法函数
!      2.   取模最大分量函数
!-----------------------------------------------------
!  Post Script :
!      1.
!      2. 
!-----------------------------------------------------
contains

subroutine solve(A,N,namda,u,tol)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-5-31
!-----------------------------------------------------
!  Purpose   : 幂法计算主特征值及主特征向量
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.  N 矩阵维数 
!       2.  A 输入矩阵 N*N维
!       3.  tol 控制精度
!  Output parameters  :
!       1.  namda 主特征值
!       2.  u     主特征向量
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.
!       2.
!----------------------------------------------------

implicit real*8(a-z)
integer::n,i,k

real*8::A(n,n)
real*8::u(n),u0(n),v(n)

!设置迭代初值向量
do i=1,n
  u0(n)=1d0
end do

u=u0

!设置模最大分量初值为0，使之进入循环
m0=0

do k=1,500

v=matmul(A,u)
call max_rou(v,n,m1)
u=v/m1

!判断迭代停止标准
if (dabs(m1-m0)<tol) exit

!更新m值
m0=m1

end do

namda=m1
end subroutine solve


subroutine max_rou(r,n,ma)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  给出向量模最大的分量
!             *是给出分量本身，而非分量的绝对值
!-----------------------------------------------------
!  Input  parameters  :
!       1.    r 输入向量
!       2.    n 输入向量的维数
!  Output parameters  :
!       1.    ma 输出结果
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.    本函数取模最大分量本身，而非取其绝对值
!       2.    例如r=(1,0,-4,3),则结果取 ma= -4
!----------------------------------------------------
implicit real*8(a-z)
integer::n,i,k

!n 为向量维数
!i 为循环指标
!k 为第k个分量为模最大的分量，这里作为标记
real*8::r(n)

ma=dabs(r(1))
do i=2,n
  if (dabs(r(i))>ma) then
   ma=dabs(r(i))
   k=i   !用k 记录指标，但是ma 不取r(i)
  end if

end do
ma=r(k)
end subroutine max_rou

end module power


program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-5-31
!-----------------------------------------------------
!  Purpose   :  主函数
!    
!-----------------------------------------------------
!  In put data  files :
!       1.
!       2.
!  Output data files  :
!       1.  result.txt 计算结果文件
!       2.
!-----------------------------------------------------
!  Post Script :
!       1.
!       2.
!-----------------------------------------------------
use power

implicit real*8(a-z)


real*8::A(3,3),u(3)


open(unit=11,file='result.txt')


a=reshape((/-1d0,2d0,1d0,&
           2d0,-4d0,1d0,&
           1d0,1d0,-6d0/),(/3,3/))

call solve(A,3,namda,u,1d-7)


write(11,101)namda,u

101 format(/,T4,'幂法计算结果',//,&  
           T3,'主特征值为：',/,F12.7,//,&
           T3,'主特征向量为：',3(/F12.7))

end program main


