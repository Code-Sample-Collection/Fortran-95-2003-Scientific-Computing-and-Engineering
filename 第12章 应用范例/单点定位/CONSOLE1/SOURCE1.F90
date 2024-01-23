module navigation
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010.07.08
!-----------------------------------------------------
!  Description :   导航模块
!  
!  Post Script :
!      1.
!      2. 
!
!-----------------------------------------------------
!  Contains    :
!      1.      主要包含解算函数
!      2.      理论值计算函数
!      3.      偏导数矩阵计算
!      4.      正交分解的最小二乘相关函数 
!-----------------------------------------------------
!  Parameters  :
!      1.
!      2.
!----------------------------------------------------- 

contains

subroutine solve(rs,x,rou,x0,tol,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose     :  解算函数
! 
!  Post Script :
!       1.    需要调用理论值函数偏导数函数
!       2.    线性最小二乘函数
!       3.    
!-----------------------------------------------------
!  Input  parameters  :
!       1.  rs 卫星星历
!       2.  rou 测量数据
!       3.  x0 初值  包含用户位置和钟差
!       4.  tol 误差容限
!       5.  资料个数
!  Output parameters  :
!       1.  x 用户位置和钟差
!       2.
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------

implicit real*8(a-z)

integer::N

real*8::rou(N),rs(N,3),x(4),x0(4)

real*8::com(N),dx(N),oc(N),H(N,3)

integer::i

x=x0

do i=1,50
  

  call compu(x,rs,com,N)
  
  oc=rou-com
  
  call jacobi(x,rs,H,n)
  
  call  lineq(H,oc,dx,N,4)
   
  x=x+dx
  
  s=0d0
  s=dx(1)**2+dx(2)**2+dx(3)**2
  s=dsqrt(s)
  if (s<TOL) exit

end do

  

end subroutine solve



subroutine jacobi(x,rs,H,n)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.07.08
!-----------------------------------------------------
!  Purpose     :   偏导数矩阵
! 
!  Post Script :
!       1.
!       2.
!       3.    
!-----------------------------------------------------
!  Input  parameters  :
!       1.  x  用户位置与钟差 
!       2.  rs 卫星星历
!       3.  n  维数 即资料个数
!  Output parameters  :
!       1.   
!       2.  H 偏导数矩阵
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------

implicit real*8(a-z)

integer::n

real*8::x(4),rs(n,3),H(n,4)

integer::i

!光速
c=299792.458  

do i=1,n
   
   a=x(1)-rs(i,1)
   b=(x(1)-rs(i,1))**2+(x(2)-rs(i,2))**2+(x(3)-rs(i,3))**2
   b=dsqrt(b)
   
   H(i,1)=a/b
   
   a=x(2)-rs(i,2)
   H(i,2)=a/b
   
   a=x(3)-rs(i,3)
   H(i,3)=a/b
   
   H(i,4)=c
    
end do

end subroutine jacobi

subroutine compu(x,rs,rou,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose     :  理论值函数
! 
!  Post Script :
!       1.
!       2.
!       3.    
!-----------------------------------------------------
!  Input  parameters  :
!       1.  x,rs  用户位置钟差，   卫星星历
!       2.  N  维数
!  Output parameters  :
!       1.  rou  理论值
!       2.
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------

implicit real*8(a-z)

integer::N

real*8::x(4),rs(N,3),rou(N)

integer::i

!光速
c=299792.458  


do i=1,n
   rou(i)=(rs(i,1)-x(1))**2+(rs(i,2)-x(2))**2  &
           +(rs(i,3)-x(3))**2
   
   rou(i)=dsqrt(rou(i))+c*x(4)

end do
end subroutine compu



subroutine  lineq(A,b,x,M,N)
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

end subroutine lineq

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


end module navigation




subroutine drive

use navigation

implicit real*8(a-z)

real*8::rou(8),rs(8,3)

integer::i,j

real*8::x(4),x0(4)


open(unit=23,file='result.txt')



rou=(/23744.370148d0,&
24192.167976d0,&
24423.302089d0,&
24249.126973d0,&
26517.149564d0,&
26973.488734d0,&
26366.329636d0,&
27190.224074d0/)

rs=reshape((/-7134.529244d0,&
-22383.700040d0,& 
-5384.901317d0,&	
637.466571d0,&	  
-11568.199533d0,& 
-28908.916747d0 ,&
-1205.651181d0,&	
16456.527324d0,&	
16113.648836d0,&	 
18533.233168d0,&	
28971.622323d0,&	 
28016.053841d0,&	 
	-3328.511543d0,&	
	-577.061760d0	,& 
28296.890128d0,&	 
12347.282494d0,&	 
 23709.205570d0 ,&
5307.245613d0  ,&
 2079.796362d0  ,&
 9347.297933d0  ,&
26977.312423d0 ,&
 6051.375658d0  ,&
 -8397.025036d0 ,&
 21199.173063d0/),(/8,3/))

!设置初值
x0=(/ -2d3,5d3,3d3,1d-5/)
call solve(rs,x,rou,x0,1d-7,8)


write(23,101)x

101 format(/,T4,'导航定位结果',//,&
          T2, '用为位置为(km)：',/,&
          T2,  3(/F18.10),//,&
          T2,   '钟差为(s)：',/,&
          T2,   F18.12 )            


end subroutine drive


program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.07.08
!-----------------------------------------------------
!  Purpose     : 定位主函数
!
!  Post Script :
!       1.   调用驱动函数
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

call  drive

end program main