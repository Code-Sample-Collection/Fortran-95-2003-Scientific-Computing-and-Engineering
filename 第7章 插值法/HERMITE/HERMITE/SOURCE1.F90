module Herm
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
!  Description :    Hermite 插值模块
!                  提供两个函数
!                  solve 函数为方法函数
!                  single 函数为单点插值函数
!-----------------------------------------------------
!  Post Script :
!      1.
!      2. 
!-----------------------------------------------------

contains

subroutine solve(n,x,y,dy,m,t,ty)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :   Hermite 插值方法函数
!                可直接对向量插值
!-----------------------------------------------------
!  Input  parameters  :
!       1.   n 节点个数
!       2.   x  节点向量 N维
!       3.   y  节点函数值  N维
!       4.   dy  节点导数值 N维
!       5.   m  需要插值向量维数 
!       6    t  需要插值的点  M维
!  Output parameters  :
!       1.   ty  插值结果 M维  -----与t 对应
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.      需要注意的是这里的输入参数N  即表示节点个数
!       2.
!----------------------------------------------------
implicit real*8(a-z)
integer::n,m

integer::i,j,k

real*8::x(n),y(n),dy(n),t(m),ty(m)

do i=1,m

call single(n-1,x,y,dy,t(i),ty(i))

end do

end subroutine solve



subroutine  single(n,x,y,dy,t,ty)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  单点Hermite 插值
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
!       1.   需要注意 : 输入参数N表示 节点数减1
!       2.      如6个节点，则N=5
!             t,ty 为标量
!----------------------------------------------------

implicit real*8(a-z)

integer::N

integer::i,j

real*8::x(0:n),y(0:n),dy(0:n)
real*8::l(0:n),dl(0:n)


!---------这段程序求得l

do i=0,n
    l(i)=1
    
    do j=0,n
      
      if (j==i) cycle
    
     l(i)=l(i)*(t-x(j))/(x(i)-x(j))
    
    end do 
    
end do


!-------这段程序求得dl

do i=0,n
  
  dl(i)=0
  
  do j=0,n
    
    if (j==i) cycle

    dl(i)=dl(i)+1d0/(x(i)-x(j))

  end do 

end do



ty=0

do i=0,n

!& 表示续行
  ty=ty+y(i)*(1d0-2d0*(t-x(i))*dl(i))*l(i)**2  &
     +dy(i)*(t-x(i))*l(i)**2
end do


end subroutine single

end module Herm



program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.05.11
!-----------------------------------------------------
!  Purpose   :  Hermite 插值主函数
!               可直接对向量插值
!
!-----------------------------------------------------
!  Post Script :
!       1.   计算结果保存在文件result.txt中
!       2.
!-----------------------------------------------------
use Herm

implicit real*8(a-z)

integer::N=3

real*8::x(3),y(3),dy(3),t(2),ty(2)

open(unit=11,file='result.txt')

x=(/1.3d0,1.6d0,1.9d0/)
y=(/0.620086d0,0.4554022d0,0.2818186d0/)
dy=(/-0.5220232d0,-0.5698959d0,-0.5811571d0/)

t=(/1.5d0,1.7d0/)

call solve(3,x,y,dy,2,t,ty)

write(11,101)
write(11,102)t
write(11,103)ty

101 format(/,T10,'Hermite插值',/)
102 format('要计算的点：',T14,2F12.5)
103 format('插值结果：',T14,2F12.5)

end