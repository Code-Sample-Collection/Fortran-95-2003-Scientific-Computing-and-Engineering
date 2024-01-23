module  spline

!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010.05.12
!-----------------------------------------------------
!  Description :   三次样条插值之周期条件模块
!    
!----------------------------------------------------- 
!  Contains    :
!      1.     solve函数 即方法函数
!      2.    
!      3.     
!-----------------------------------------------------
!  Post Script :
!      1.
!      2.     可以直接对向量插值
!-----------------------------------------------------

contains


subroutine solve(n,x,y,nt,t,ty)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.05.12
!-----------------------------------------------------
!  Purpose   :   三次样条之周期条件
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.   n-----插值节点个数减1，如有九个节点 则N=8
!       2.   x ---节点自变量  为（0：N）维向量
!       3.   y----节点因变量  （0：N）维向量
!       4.   nt 要计算向量的维数
!       5.   t 要计算的向量  (1:nt) 维向量
!  Output parameters  :
!       1.   ty ---要计算的向量，（1：nt）维
!
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!            
!            对于要插值的向量分量，可以不必按大小排列
!----------------------------------------------------

implicit real*8(a-z)

integer::n,nt
!n为插值节点数减1，即如果有9个节点，则n=8
!nt 要插值点的个数，及向量t,ty的维数

integer::i,j,k

real*8::x(0:n),y(0:n)

real*8::t(nt),ty(nt)

real*8::h(0:n-1)

real*8::f1(0:n-1),f2(1:n-1)

real*8::u(1:n),namda(1:n),d(1:n)

real*8::M(1:n)


real*8::A(1:n,1:n)



do i=0,n-1
  h(i)=x(i+1)-x(i)
  f1(i)=(y(i+1)-y(i))/h(i)
end do





!求得 u, namda, d
do i=1,n-1

 u(i)=h(i-1)/(h(i-1)+h(i))
 namda(i)=1-u(i)
 
 f2(i)=(f1(i-1)-f1(i))/(x(i-1)-x(i+1))
  
  d(i)=6d0*f2(i) 
  
end do


!对于周期条件而言，系数设置

u(n)=h(n-1)/(h(0)+h(n-1))

namda(n)=1-u(n)

!右向量
d(n)=6d0/(h(0)+h(n-1))*(f1(0)-f1(n-1))




!设置A矩阵值
A=0
do i=1,n
a(i,i)=2d0
end do

do i=2,n
a(i,i-1)=u(i)
end do

do i=1,n-1
a(i,i+1)=namda(i)
end do



!设置 A矩阵的首行元素
a(1,n)=u(1)


!设置A矩阵元素的末行元素
a(n,1)=namda(n)


!注意三种条件下，用Gauss法解线性方程时 N 大小设置
call gauss(a,d,M,N)



!--------以上以及求得系数
!已经完成插值多项式的建立

!------------------------------------------------
! 以下开始计算具体值
do k=1,nt

!------------
!  对要插值向量每个分量而言，先找到其在数据中的位置
do i=1,n-1

 if (t(k)<x(i+1)) exit 

end do


  ty(k)=M(i)*(x(i+1)-t(k))**3/6d0/h(i)+ M(i+1)*(t(k)-x(i))**3/6d0/h(i)+ &

     (y(i)-M(i)*h(i)**2/6d0)*(x(i+1)-t(k))/h(i)+ &
     (y(i+1)-M(i+1)*h(i)**2/6d0)*(t(k)-x(i))/h(i)

end do

!-------------------------------




!为了方便读者比对结果，这里输出中间值，需要这些数值时可以去掉注释
!!-------------------------
!write(11,101)
!write(11,102)((i,h(i),m(i)),i=0,n-1)   !输出h
!
!
!write(11,103)
!write(11,104)((i,u(i),namda(i),d(i)),i=1,n-1)
!
!101 format(/,T5,'          h             m')
!102 format(I3,2F16.8)
!
!103 format(/,T5,'         u             namda            d')
!104 format(I3,3F16.8)
!!----------------------------------------------------------------

end subroutine solve




subroutine gauss(A,b,x,N)
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

end subroutine gauss



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


end module spline




program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  周期边界条件下的三次样条插值
!    
!-----------------------------------------------------
!  In put data  files :
!       1.
!       2.
!  Output data files  :
!       1.
!       2.   reslut.txt文件保存了计算结果
!-----------------------------------------------------
!  Post Script :
!       1.
!       2.
!-----------------------------------------------------

use spline

implicit real*8(a-z)

integer::i,j

real*8::x(0:21),y(0:21),t(5),ty(5),simty(5)


open(unit=11,file='result.txt')


x=(/   0d0,&
    0.3000d0,&
    0.6000d0,&
    0.9000d0,&
    1.2000d0,&
    1.5000d0,&
    1.8000d0,&
    2.1000d0,&
    2.4000d0,&
    2.7000d0,&
    3.0000d0,&
    3.3000d0,&
    3.6000d0,&
    3.9000d0,&
    4.2000d0,&
    4.5000d0,&
    4.8000d0,&
    5.1000d0,&
    5.4000d0,&
    5.7000d0,&
    6.0000d0,&
    6.2832d0 /)

y=(/0.6000d0,&
    0.8687d0,&
    1.0598d0,&
    1.1563d0,&
    1.1495d0,&
    1.0399d0,&
    0.8375d0,&
    0.5603d0,&
    0.2330d0,&
   -0.1151d0,&
   -0.4529d0,&
   -0.7502d0,&
   -0.9806d0,&
   -1.1233d0,&
   -1.1657d0,&
   -1.1040d0,&
   -0.9437d0,&
   -0.6990d0,&
   -0.3919d0,&
   -0.0499d0,&
    0.2967d0,&
    0.6000d0  /)


!要进行计算的点  ------不需要按照大小排列
t=(/0.4d0,2.8d0,5.9d0,6.1d0,4.7d0/)

!调用方法函数
call solve(21,x,y,5,t,ty)                                                         


!仿真函数值
simty=dsin(t)+dcos(t)*0.6d0


write(11,101)

write(11,102)
write(11,103)((i,t(i),ty(i),simty(i)),i=1,5)

101  format(/,T9,'三次样条插值之周期条件',/)
102  format('  序列     插值点        插值结果        真值')
103  format(I3,3F16.8)

end program main