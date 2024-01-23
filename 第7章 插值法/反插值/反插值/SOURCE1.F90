module  newton
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
!  Description :   反插值问题
!    
!-----------------------------------------------------
!  Post Script :
!      1.   提供两个函数，一个是单点插值，一个是solve函数
!           注意两者的N 是不一样的。
!           在单点插值中N是节点的个数减1
!           为方便使用，在solve中N就直接指节点个数
!-----------------------------------------------------

contains

subroutine  solve(n,x,y,m,t,ty)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  方法主函数
!               可直接对多点插值
!-----------------------------------------------------
!  Input  parameters  :
!       1.  n  节点个数
!       2.  x 节点自变量值 N维向量
!       3.  y 节点因变量值  N维向量
!       4.  m  要计算点的个数
!       5.  t  要计算的点   M向量
!  Output parameters  :
!       1.  ty  计算结果，为M维向量
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.   注意这里的N  就直接指节点的个数
!       2.   而子函数  new中 N 是指节点个数减 1 
!----------------------------------------------------

implicit real*8(a-z)
integer::n,m

integer::i,j,k

real*8::x(n),y(n),t(m),ty(m)

do i=1,m

call new(n-1,x,y,t(i),ty(i))

end do

end subroutine solve




subroutine new(n,x,y,t,ty)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  单点插值程序
!               此函数目的为 solve所调用
!-----------------------------------------------------
!  Input  parameters  :
!       1.
!       2.   t 标量
!  Output parameters  :
!       1.   
!       2.    ty 标量
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.    N 是节点个数 减 1， 比如共4个节点  则N=3
!       2.
!----------------------------------------------------
implicit real*8(a-z)
integer::n

integer::i,j,k

real*8::x(0:n),y(0:n)


real*8::b(n+1)
real*8::Q(0:n,0:n)


do i=0,n

Q(i,0)=y(i)

end do


do i=1,n
  
  do j=1,i
   
   Q(i,j)=(Q(i,j-1)-Q(i-1,j-1))/(x(i)-x(i-j))
  
  end do

end do


b(n+1)=Q(n,n)


do k=n,1,-1
  
  b(k)=Q(k-1,k-1)+b(k+1)*(t-x(k-1))

end do

ty=b(1)

end subroutine new



end module newton


program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.05.12
!-----------------------------------------------------
!  Purpose   :  反插值问题主函数
!                solve函数可以直接对向量插值
!-----------------------------------------------------
!  In put data  files :
!       1.
!       2.
!  Output data files  :
!       1.
!       2.
!-----------------------------------------------------
!  Post Script :
!       1.
!       2.
!-----------------------------------------------------

call drive


end program main



subroutine  drive
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.05.12
!-----------------------------------------------------
!  Purpose   :  反插值问题驱动函数
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


use newton

implicit real*8(a-z)

integer::N=5,M=2

real*8::x(5),y(5),t(2),ty(2)

open(unit=11,file='result.txt')

!插值节点极其函数值
x=(/-1.3d0,-1.4d0,-1.5d0,-1.6d0,-1.7d0/)
y=(/3.033d0,1.776d0,0.375d0,-1.176d0,-2.883d0/)

!要计算的点
t=(/0d0,-1.6d0/)

call solve(n,y,x,m,t,ty)
write(11,101)
write(11,102)ty

101 format(/,T5,'反插值结果为:')
102 format(/,2F16.8)


end subroutine drive


