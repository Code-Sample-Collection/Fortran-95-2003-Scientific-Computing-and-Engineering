module four_iter
!----------------------------------------module coment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-8
!-----------------------------------------------------
!  Purpose   :  多重迭代法迭代计算非线性方程
!    
!-----------------------------------------------------
!  Parameters  :
!      1.     MAX最大允许迭代次数
!      2.     tol误差容许限
!      3.     迭代初值
!  Contains    :
!      1.    solve 迭代方法函数
!      2.    func  要计算的方程函数
!      3.    dfunc 导函数
!-----------------------------------------------------

implicit real*8(a-z)
integer:: MAX=200
real*8::tol=1D-7
real*8::x0=2d0

contains  


subroutine solve(x,fx,iter)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  4阶收敛多重迭代法计算方程的根
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
!       1.
!       2.
!----------------------------------------------------
implicit real*8(a-z)
integer i,iter

x1=x0
do i=1,MAX
     
    !计算函数一次    
    fx=func(x1)
    !计算导函数一次
    dfx=dfunc(x1)
    
    
    y=x1-fx/dfx
    
    fy=func(y)
    
    x2=y-fy*(y-x1)/(2*fy-fx)
    
    dx=dabs(x2-x1)
    
    !该段为输出中间结果，为了用于分析，实际计算时，可以注释掉
    write(102,*)i,x2,func(x2)
    !---------
    
    
    if (dx<tol) exit
    x1=x2
    

    
end do
    x=x2
    fx=func(x2)
    iter=i
end subroutine solve


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
!       1.   func方程导函数值
!----------------------------------------------------
implicit real*8(a-z)
   func=dsin(x)-x/25d0
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
!       1.   dfunc方程函数值
!----------------------------------------------------
implicit real*8(a-z)
  dfunc=dcos(x)-1/25d0
end function dfunc


end module  four_iter


program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  4阶收敛多重迭代法迭代计算方程的根主函数
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
use four_iter
implicit real*8(a-z)
integer::iter
open(unit=101,file='result.txt')
open(unit=102,file='Im_result.txt')

call solve(x,fx,iter)
write(101,501)x,fx,iter
501 format(T20,'4阶收敛迭代法计算方程的根',//,&
            3x,'x=   ',F12.8,/,&
            3x,'F(x)=',F12.8,/,&
            3x,'iter=',I5)


end program main