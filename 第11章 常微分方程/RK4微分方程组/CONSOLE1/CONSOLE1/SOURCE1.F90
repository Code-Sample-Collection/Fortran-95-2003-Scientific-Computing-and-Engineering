module m_rk_ods
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
!  Description :   RK4经典方法计算微分法方程组方法模块
!    
!-----------------------------------------------------
!  Parameters  :
!      1.     M ---积分区间划分的细度
!      2.
!----------------------------------------------------- 
!  Contains    :
!      1.      RK   RK方法
!      2.      fun  需要计算的函数
!-----------------------------------------------------

integer::M=100 !积分区间划分N步

contains 

subroutine  solve(func,t0,tt,y0,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :      4阶经典RK方法计算微分方程组方法函数
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.           func  传递的右函数名
!       2.           t0积分初时刻 
!       3.           tt积分结束时刻
!       4.           y0 初状态
!       5.           N 方程组的维数       
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

implicit real*8(a-z)
!方程组的维数
external func
integer::N    
real*8::y0(N),y(N)

real*8::k1(N),k2(N),k3(N),k4(N)

h=(tt-t0)/M

t=t0
y=y0

do i=1,M

     call func(k1,t,y,N)
     call func(k2,t+h/2,y+h/2*k1,N)
     
     call func(k3,t+h/2,y+h/2*k2,N)
     call func(k4,t+h,y+h*k3,N)
     
     y=y+(k1+2*k2+2*k3+k4)*h/6
     
     t=t0+i*h
     
     write(11,101)t,y
   
  end do 
  101 format(<M+1>F12.6)

end subroutine solve


subroutine fun1(f,t,y,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  需要计算的方程右函数
!               
!-----------------------------------------------------
!  Input  parameters  :
!       1.    N  方程阶数
!       2.    t,y  
!  Output parameters  :
!       1.    f 
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.  注意：在自治系统中，其实t是不起作用的
!       2.         如本例
!----------------------------------------------------

implicit real*8(a-z)
integer::N
real*8::f(N),y(N)

f(1)=y(2)*y(3)
f(2)=-y(1)*y(3)
f(3)=-0.51*y(1)*y(2)

end subroutine fun1



end module m_rk_ods




  program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-12
!-----------------------------------------------------
!  Purpose   :  采用RK4方法计算微分方程组
!    
!-----------------------------------------------------
!  In put data  files :
!       1.   fout.txt  存放计算结果
!       2.
!  Output data files  :
!       1.
!       2.
!-----------------------------------------------------
!  Post Script :
!       1.
!       2.
!-----------------------------------------------------

  use m_rk_ods
  
  implicit real*8(a-z)
  integer::N
  real*8::y(3),y0(3)

 
  open(unit=11,file='fout1.txt')
    
  write(11,101)
  101 format(/,T5,'采用RK4方法计算方程组',/)
  
  N=3
  
  t0=0
  tt=12
  
  y0=(/0,1,1/)
  
  call solve(fun1,t0,tt,y0,N)  
  
  end