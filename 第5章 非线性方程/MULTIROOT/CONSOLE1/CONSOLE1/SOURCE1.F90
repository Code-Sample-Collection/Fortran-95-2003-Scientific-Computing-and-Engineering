module multiroot
!----------------------------------------module coment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-3
!-----------------------------------------------------
!  Purpose   :  重根时改进算法模块
!    
!-----------------------------------------------------
!  Parameters  :
!      1.     MAX最大允许迭代次数
!      2.     tol误差容许限
!      3.     迭代初值
!  Contains    :
!      1.    solve 牛顿迭代方法函数
!      2.    func  要计算的方程函数
!      3.    dfunc 导函数
!-----------------------------------------------------

implicit real*8(a-z)
integer:: MAX=200
real*8::tol=1D-7
real*8::x0=1.5d0

contains  

subroutine solve(x,fx,iter)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  重根改进算法
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
!       1.   引用函数 func  dfunc  d2func
!       2.
!----------------------------------------------------
implicit real*8(a-z)
integer i,iter

x1=x0

!标题，实际计算时可以注释掉
write(102,501)
501 format(T20,'重根计算的中间结果')
! --------------

do i=1,MAX

    temp1=func(x1)*dfunc(x1)
    temp2=dfunc(x1)**2-func(x1)*d2func(x1)
    
    x2=x1-temp1/temp2
    
    dx=dabs(x2-x1)
    if (dx<tol) exit
    x1=x2
    
    !该段为输出中间结果，为了用于分析，实际计算时，可以注释掉
    write(102,*)i,x2,func(x2)
    !---------
    
    
end do
    x=x2
    fx=func(x2)
    iter=i-1
end subroutine solve


subroutine newton(x,fx,iter)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  牛顿法计算方程的根
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
!       1.   引用函数 func  dfunc
!      
!----------------------------------------------------
implicit real*8(a-z)
integer i,iter

!标题，实际计算时可以注释掉
write(104,502)
502 format(T20,'牛顿法计算的中间结果')
! --------------

x1=x0
do i=1,MAX

    x2=x1-func(x1)/dfunc(x1)
    
    dx=dabs(x2-x1)
    if (dx<tol) exit
    x1=x2
    
    !该段为输出中间结果，为了用于分析，实际计算时，可以注释掉
    write(104,*)i,x2,func(x2)
    !---------
    
    
end do
    x=x2
    fx=func(x2)
    iter=i-1
end subroutine newton


function func(x)
!---------------------------------subroutine  comment
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
   func=x**4-4*x**2+4d0
end function func

function dfunc(x)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  方程导函数
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.    x自变量
!       2. 
!  Output parameters  :
!       1.   dfunc方程导函数
!----------------------------------------------------
implicit real*8(a-z)
  dfunc=4*x**3-8*x
end function dfunc

function d2func(x)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  方程二阶导函数
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.    x自变量
!       2. 
!  Output parameters  :
!       1.   d2func方程二阶导数
!----------------------------------------------------
implicit real*8(a-z)
d2func=12*x**2-8d0
end function d2func


end module  multiroot


program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :   重根问题主函数
!    
!-----------------------------------------------------
!  In put data  files :
!       1.   
!       2.
!  Output data files  :
!       1.   result_m.txt计算结果
!       2.   IM_result_m.txt计算的中间结果
!       3.   result_n.txt牛顿法计算结果
!       4.   Im_result_n.txt 牛顿法计算的中间结果
!-----------------------------------------------------
use multiroot
implicit real*8(a-z)
integer::iter_m,iter_n
open(unit=101,file='result_m.txt')
open(unit=102,file='Im_result_m.txt')

open(unit=103,file='result_n.txt')
open(unit=104,file='Im_result_n.txt')

call solve(x_m,fx_m,iter_m)
write(101,501)x_m,fx_m,iter_m
501 format(T20,'重根改进计算结果',//,&
            3x,'x=   ',F15.10,/,&
            3x,'F(x)=',F15.10,/,&
            3x,'iter=',I5)


call newton(x_n,fx_n,iter_n)
write(103,502)x_n,fx_n,iter_n
502 format(T20,'牛顿法计算结果',//,&
            3x,'x=   ',F15.10,/,&
            3x,'F(x)=',F15.10,/,&
            3x,'iter=',I5)

end program main