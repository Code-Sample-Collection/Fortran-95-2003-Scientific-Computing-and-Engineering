module mody_Rayleigh
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-6-1
!-----------------------------------------------------
!  Description :   修正Rayleigh加速方法
!    
!-----------------------------------------------------
!  Parameters  :
!      1.
!      2.
!----------------------------------------------------- 
!  Contains    :
!      1.   方法函数
!      2.   内积函数
!      3.   范数函数
!-----------------------------------------------------

contains

subroutine solve(A,N,namda,u,tol)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-6-1
!-----------------------------------------------------
!  Purpose   : Rayleigh加速计算对称矩阵特征值函数
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

call norm(u0,rou,n)

!u 已经单位化
u=u0


m0=0

do k=1,500

v=matmul(A,u)


call vvdot(v,u,m1,n)

call norm(v,rou,n)

u=v/rou


!判断迭代停止标准
if (dabs(m1-m0)<tol) exit

!更新m值
m0=m1

end do

namda=m1

!特征向量归一化
call norm(u,rou,n)

u=u/rou
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


end module mody_Rayleigh


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
use mody_Rayleigh

implicit real*8(a-z)


real*8::A(3,3),u(3)


open(unit=11,file='result.txt')


a=reshape((/1.867310185256756d0,   0.676506748114136d0,   1.832446222432665d0,&
   0.676506748114136d0,   0.438313983411573d0,   0.927636560784582d0,&
   1.832446222432665d0,   0.927636560784582d0,   2.421070728983985d0/),(/3,3/))


call solve(A,3,namda,u,1d-7)


write(11,101)namda,u

101 format(/,T4,'修正Rayleigh加速方法',//,&  
           T3,'主特征值为：',/,F12.7,//,&
           T3,'主特征向量为：',3(/F12.7))

end program main


