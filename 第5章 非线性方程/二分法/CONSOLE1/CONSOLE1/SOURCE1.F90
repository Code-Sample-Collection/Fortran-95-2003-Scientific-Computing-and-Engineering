module bisect
!----------------------------------------module coment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-3
!-----------------------------------------------------
!  Purpose   :  不动点迭代法计算方程的根
!    
!-----------------------------------------------------
!  Parameters  :
!      1.    MAX最大允许迭代次数
!      2.    tol误差容限
!      3.    a,b  为区间 [a,b]
!  Contains    :
!      1.   方程函数     
!      2.   方法函数
!-----------------------------------------------------
!  Post Script :
!      1.   
!      2. 
!-----------------------------------------------------

implicit real*8(a-z)
integer ::MAX=200
real*8::tol=1D-7
real*8::a=2d0,b=4d0


contains

subroutine solve(x,fx,iter)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-3
!-----------------------------------------------------
!  Purpose   :  二分法函数
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.   
!       2. 
!  Output parameters  :
!       1.   x 方程的根
!       2.   fx 该点的函数值
!       3.   iter 实际迭代次数
!  Post Script : 
!       1.   
!       2.
!----------------------------------------------------
implicit real*8(a-z)
integer iter,i
do i=1,MAX
   c=(a+b)/2D0
   
   if (func(c)==0) exit
   
   
   if (func(a)*func(c)<0d0) then
        b=c
   else
        a=c
   end if
   
   dfx=dabs(func(a)-func(b))

   if (dfx<tol) exit
   
   
   !--------这段用于中间输出结果，这里用于分析，实际计算中可以注释掉
   write(102,*)i,c,func(c)
   !-----------     

end do
   
   x=c
!根
   fx=func(c)
!函数值
  iter=i-1
!实际迭代次数  

end subroutine 



function func(x)
!---------------------------------function  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-3
!-----------------------------------------------------
!  Purpose   :  要计算的方程函数 ，其中 f(x)=0
!    
!-----------------------------------------------------

 implicit real*8(a-z)
  func=(x-1d0)**3-3d0*x+2d0
end function func


end module bisect


program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010年4月3日
!-----------------------------------------------------
!  Purpose   :  二分法主函数
!    
!-----------------------------------------------------
!  In put data  files :
!       1.      
!       2.
!  Output data files  :
!       1.    result.txt   计算结果文件
!       2.   IM_result.txt  中间结果 
!-----------------------------------------------------
!  Post Script :
!       1.
!       2.
!-----------------------------------------------------


use bisect
implicit real*8(a-z)
integer iter

open(unit=101,file='result.txt')

!记录中间结果
open(unit=102,file='Im_result.txt')  

call solve(x,fx,iter)

write(101,501)x,fx,iter
501 format(T20,'二分法计算方程的根',//,&
            3x,'x=   ',F15.10,/,&
            3x,'F(x)=',F15.10,/,&
            3x,'iter=',I5)

end program main




