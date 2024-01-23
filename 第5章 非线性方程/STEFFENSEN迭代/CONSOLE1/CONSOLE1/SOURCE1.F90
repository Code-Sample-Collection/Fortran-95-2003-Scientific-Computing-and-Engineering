module Steffensen
!----------------------------------------module coment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-3
!-----------------------------------------------------
!  Purpose   :  Steffensen加速计算方程的根
!    
!-----------------------------------------------------
!  Parameters  :
!      1.     MAX最大允许迭代次数
!      2.     tol误差容许限
!      3.     迭代初值
!  Contains    :
!      1.    solve Steffensen方法函数
!      2.    func  要计算的方程函数
!      3.    picard  不动点迭代法的方法函数
!-----------------------------------------------------

implicit real*8(a-z)
integer:: MAX=200
real*8::tol=1D-7
real*8::x0=0.5d0

contains  

subroutine solve(x,fx,iter)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  Steffensen加速方法处理函数
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.    
!       2. 
!  Output parameters  :
!       1.   x 方程的根
!       2.   fx 该点的函数值
!       3.   iter 实际迭代次数
!-----------------------------------------------------
implicit real*8(a-z)
integer i,iter

!这句为标题，计算可以注释掉
    write(102,*)'Steffensen迭代法中间结果'
!------------------

x1=x0
do i=1,MAX
    
    x2=func(x1)+x1
    x3=func(x2)+x2
    
    x1_bar=x1-(x2-x1)**2/(x3-2*x2+x1)
 
        
    dx=dabs(x1_bar-x1)
    
    
    !该段为输出中间结果，为了用于分析，实际计算时，可以注释掉
    write(102,*)i,x1_bar,func(x2)
    !---------
    
    
    if (dx<tol) exit
    
    x1=x1_bar
    
    
end do
    x=x1_bar
    fx=func(x2)
    iter=i
end subroutine solve


subroutine picard(x,fx,iter)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  不动点迭代计算方程的根
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.    
!       2. 
!  Output parameters  :
!       1.   x 方程的根
!       2.   fx 该点的函数值
!       3.   iter 实际迭代次数
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1. 此函数用于对比Steffensen方法，实际使用时可以删除该Subroutine
!       2.
!----------------------------------------------------
implicit real*8(a-z)
integer i,iter

!这句为标题，计算可以注释掉
    write(104,*)'picard迭代法中间结果'
!------------------
x1=x0
do i=1,MAX
    
    x2=func(x1)+x1
    
    dx=dabs(x2-x1)
    
    !该段为输出中间结果，为了用于分析，实际计算时，可以注释掉
    write(104,*)i,x2,func(x2)
    !---------
    
    
    if (dx<tol) exit
    x1=x2
    

    
end do
    x=x2
    fx=func(x2)
    iter=i
end subroutine picard


function func(x)
!---------------------------------function  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  方程函数
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.    x自变量
!       2. 
!  Output parameters  :
!       1.   func方程函数值
!-----------------------------------------------------

implicit real*8(a-z)
 func=dexp(-x)-x
end function func

end module  Steffensen


program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :   Steffensen迭代法主函数
!    
!-----------------------------------------------------
!  In put data  files :
!       1.   
!       2.
!  Output data files  :
!       1.   result.txt计算结果
!       2.   IM_result.txt计算的中间结果
!-----------------------------------------------------
!  Post Script :
!       1.
!       2.
!-----------------------------------------------------
use Steffensen
implicit real*8(a-z)
integer::iter_s,iter_p
open(unit=101,file='result_s.txt')
open(unit=102,file='Im_result_s.txt')

open(unit=103,file='result_p.txt')
open(unit=104,file='Im_result_p.txt')


call picard(x_p,fx_p,iter_p)
write(103,501)x_p,fx_p,iter_p
501 format(T20,'Picard不动点迭代法计算方程的根',//,&
            3x,'x=   ',F15.10,/,&
            3x,'F(x)=',F15.10,/,&
            3x,'iter=',I5)


call solve(x_s,fx_s,iter_s)
write(101,502)x_s,fx_s,iter_s
502 format(T20,'Steffensen加速计算方程的根',//,&
            3x,'x=   ',F15.10,/,&
            3x,'F(x)=',F15.10,/,&
            3x,'iter=',I5)

end program main
