module m_house
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2009.07.30
!-----------------------------------------------------
!  Description :   Houlseholder求QR分解之模块
!  
!  Post Script :
!      1.     
!      2. 
!
!----------------------------------------------------- 
contains
subroutine house_qr(A,Q,R,M,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2009.07.30
!-----------------------------------------------------
!  Purpose     : 采用Householder做QR分解
! 
!  Post Script :
!       1.     A=QR
!       2.     Q为正交矩阵
!       3.    
!-----------------------------------------------------
!  Input  parameters  :
!       1.   A 要求分解之矩阵
!       2.   M,N 矩阵维数
!  Output parameters  :
!       1.   Q  正交矩阵
!       2.   R  上三角矩阵
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------

implicit real*8(a-z)
integer::M,N
real*8::A(M,N),Q(M,M),R(M,N)

real*8::H0(M,M),H1(M,M),H2(M,M),qt(M,M)

real*8::A1(M,N),A2(M,N),u(M)

integer::i,k,j


A1=A

H1=0d0
do j=1,M
 H1(j,j)=1d0
end do

!k表示对所有的列
do k=1,n
  
 !设置H矩阵初值，这里设置为单位矩阵
  H0=0d0
  
  !设置为单位矩阵
  do i=1,M
  H0(i,i)=1d0
  end do
  
  s=0d0
  do i=k,M
  
  s=s+a1(i,k)*a1(i,k)
  
  end do

!算的向量的2范数  
  s=dsqrt(s)

  u=0d0

!---------------------------------
! 这段甚为重要，关系到数值稳定性问题
! 目的是使得u的2范数尽可能大
! 原则是，如果首元素大于零，则u的第一个元素是正+正
!         如果首元素小于零，则u的第一个元素是负-负
  if (a1(k,k)>=0) then
  u(k)=a1(k,k)+s
  else 
  u(k)=a1(k,k)-s
  end if
!-------------------------------  
  
  do i=k+1,m
  u(i)=a1(i,k)
  end do
  
  du=0
  do i=k,m
  !求的单位u 长度平方
  du=du+u(i)*u(i)
  end do
  
  
  !计算得到大的H矩阵
  do i=k,m
    do j=k,m
        
        H0(i,j)=-2d0*u(i)*u(j)/du
                  
        if (i==j) then
        H0(i,j)=1d0+H0(i,j)
        end if
        
    end do
  end do
  
!千万要注意矩阵相乘的次序
!先更新矩阵
A2=matmul(H0,A1)
A1=A2

!H1初值为单位矩阵，后逐步更新
H1=matmul(H1,H0)


!更新变换后的矩阵


!-------------------------------------------
!--这两行语句可以删除，这里为显示中间结果
!--方便比对
write(11,101)((H0(i,j),j=1,m),i=1,m)
write(11,102)((H1(i,j),j=1,m),i=1,m)
write(11,103)((A1(i,j),j=1,N),i=1,m)
!-------------------------------------------
end do

!中间结果格式控制
101 format(/,'H0=',/,<M>(<M>F10.5/))
102 format(/,'H=',/,<M>(<M>F10.5/))
103 format(/,'A=',/,<m>(<n>F10.5/))
!------------------------------------------
Q=H1

R=A1

end subroutine house_qr


end module m_house


subroutine drive(m,n)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose     : 驱动程序
! 
!  Post Script :
!       1.
!       2.
!       3.    
!-----------------------------------------------------

implicit real*8(a-z)

use m_house
implicit real*8(a-z)

integer::m,n
integer::i,j

real*8::A(m,n),Q(m,m),R(m,n)

read(10,*)((A(i,j),j=1,n),i=1,m)

write(11,100)
100 format('QR分解中间结果：')

call house_qr(A,Q,R,m,n)

write(11,101)
101 format(T10,'QR分解之结果：',/,T3,'Q=',/)

write(11,102)((Q(i,j),j=1,M),i=1,M)

!变量输出格式只针对IVF编译器，在CVF中不支持
102 format(<M>(<M>F10.5/))

write(11,103)
103 format(T3,'R=',/)

write(11,104)((R(i,j),j=1,n),i=1,m)

!变量输出格式只针对IVF编译器，在CVF中不支持
104 format(<m>(<N>F10.5/))

end subroutine drive

program main

!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2009.07.30
!-----------------------------------------------------
!  Purpose     : Househoulder镜像变幻求QR分解主函数
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
!       1.  result.txt  输出结果
!       2.
!
!-----------------------------------------------------
implicit real*8(a-z)
integer::m,n

open(unit=10,file='fin.txt')
open(unit=11,file='result.txt')

read(10,*)m,n
call drive(m,n)


end program main

