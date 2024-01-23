module secant
!----------------------------------------module coment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-3
!-----------------------------------------------------
!  Purpose   :  割线法计算非线性方程
!    
!-----------------------------------------------------
!  Parameters  :
!      1.     MAX最大允许迭代次数
!      2.     tol误差容许限
!      3.     a,b区间初值
!  Contains    :
!      1.    solve 牛顿迭代方法函数
!      2.    func  要计算的方程函数
!-----------------------------------------------------

implicit real*8(a-z)
integer:: MAX=200
real*8::tol=1D-7
real*8::a=2d0,b=1d0

contains  


subroutine solve(x,fx,iter)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   : 割线法计算方程的根
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

implicit real*8(a-z)
integer i,iter

x0=a
x1=b

do i=1,MAX

    x2=x1-(x1-x0)*func(x1)/(func(x1)-func(x0))
    
    dx=dabs(x2-x0)
    if (dx<tol) exit
    
    x0=x1
    x1=x2
    
    !该段为输出中间结果，为了用于分析，实际计算时，可以注释掉
    write(102,*)i,x2,func(x2)
    !---------
    
    
end do
    x=x2
    fx=func(x2)
    iter=i-1
end subroutine solve


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
!----------------------------------------------------
implicit real*8(a-z)
   func=2*x**3-5*x-1
end function func

end module  secant


program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :   割线法计算方程的根主函数
!    
!-----------------------------------------------------
!  In put data  files :
!       1.   
!       2.
!  Output data files  :
!       1.   result.txt计算结果
!       2.   IM_result.txt计算的中间结果
!-----------------------------------------------------
use secant
implicit real*8(a-z)
integer::iter
open(unit=101,file='result.txt')
open(unit=102,file='Im_result.txt')

call solve(x,fx,iter)
write(101,501)x,fx,iter
501 format(T20,'割线法计算方程的根',//,&
            3x,'x=   ',F15.10,/,&
            3x,'F(x)=',F15.10,/,&
            3x,'iter=',I5)


end program main