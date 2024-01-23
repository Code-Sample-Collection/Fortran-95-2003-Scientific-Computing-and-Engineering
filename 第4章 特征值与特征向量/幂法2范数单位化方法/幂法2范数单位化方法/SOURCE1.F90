module power
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-5-31
!-----------------------------------------------------
!  Description :   幂法2范数单位化方法模块
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
!  Purpose   : 幂法2范数单位化方法函数
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
real*8::u(n),u0(n),v(n),vt(n)

!设置迭代初值向量
do i=1,n
  u0(n)=1d0
end do

u=u0

!设置模最大分量初值为0，使之进入循环
m0=0

do k=1,500

v=matmul(A,u)
call norm(v,rou,n)
u=v/rou

vt=matmul(A,u)

call vvdot(u,vt,m1,n)


!判断迭代停止标准
if (dabs(m1-m0)<tol) exit

!更新m值
m0=m1

end do

namda=m1
end subroutine solve


subroutine norm(r,rou,n)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-6-1
!-----------------------------------------------------
!  Purpose   :  计算向量2范数
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.  r 输入向量
!       2.  n 向量维数
!  Output parameters  :
!       1.  rou 向量长度（2范数）
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.
!       2.
!----------------------------------------------------

implicit real*8(a-z)
integer::n,i

real*8::r(n)

rou=0d0

do i=1,n


rou=rou+r(i)**2

end do

rou=dsqrt(rou)

end subroutine norm


subroutine vvdot(r1,r2,dot,n)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-6-1
!-----------------------------------------------------
!  Purpose   :  向量内积
!               r1(1)*r2(1)+r1(2)*r2(2)+...
!-----------------------------------------------------
!  Input  parameters  :
!       1.  r1,r2 欲求之向量
!       2.  n 维数
!  Output parameters  :
!       1.  dot 向量内积
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.
!       2.
!----------------------------------------------------

implicit real*8(a-z)
integer::n,i
real*8::r1(n),r2(n)

dot=0d0
do i=1,n

dot=dot+r1(i)*r2(i)

end do

end subroutine vvdot


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


a=reshape((/2d0,10d0,3d0,&
           3d0,3d0,6d0,&
           2d0,4d0,1d0/),(/3,3/))

call solve(A,3,namda,u,1d-7)


write(11,101)namda,u

101 format(/,T4,'幂法2范数单位化计算结果',//,&  
           T3,'主特征值为：',/,F12.7,//,&
           T3,'主特征向量为：',3(/F12.7))

end program main


