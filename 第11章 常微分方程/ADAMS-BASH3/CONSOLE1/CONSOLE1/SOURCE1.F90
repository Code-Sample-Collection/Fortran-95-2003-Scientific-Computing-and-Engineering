module  m_adam3
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-4-13
!-----------------------------------------------------
!  Description :   3阶Adms方法计算微分法方程组模块
!    
!----------------------------------------------------- 
!  Contains    :
!      1.      adam3方法函数
!      2.      RK4 -修改的RK4阶方法提供起步
!      3.      func 需要计算的函数
!-----------------------------------------------------
!  Post Script :
!      1.
!      2. 
!-----------------------------------------------------
contains

subroutine solve(ta,tb,ya,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :3阶Adms方法计算微分法方程组方法函数
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.  ta 起步时间
!       2.  tb 介绍时间
!       3.  ya 初状态
!       4.  N  方程组维数
!
!----------------------------------------------------
!  Post Script :
!       1.
!       2.
!----------------------------------------------------


implicit real*8(a-z)

integer::N,M=100
integer::i
real*8::ya(N),y0(N),y1(N),y2(N),y3(N),y4(N)
real*8::f0(N),f1(N),f2(N),f3(N)

h=(tb-ta)/M


t0=ta
y0=ya

t1=ta+h
t2=t1+h
 
     
call RK4(t0,t1,y0,y1,N)
!算出y1

call RK4(t1,t2,y1,y2,N)
!算出y2



     call  func(f0,t0,y0,N)
     !算出f0
     
     call  func(f1,t1,y1,N)
     !算出f1
     
     call  func(f2,t2,y2,N)
     !算出f2
       


do i=3,M
  
     y3=y2+h/12d0*(23d0*f2-16d0*f1+5d0*f0) 

     t3=t2+h     
     
     !输出计算结果
  
     write(11,101)t3,y3
     
     !  以下为各量更新
     t2=t2+h
     
	 y2=y3
	  
     f0=f1
     f1=f2
 
	 call func(f2,t2,y2,N)
    
end do
  
101 format(T3,<N+1>f10.5)


end subroutine solve



subroutine  RK4(t0,tt,y0,yt,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-13
!-----------------------------------------------------
!  Purpose   :  修改后的4阶经典RK方法计算微分方程组方法函数
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.           
!       2.           t0积分初时刻 
!       3.           tt积分结束时刻
!       4.           y0 初状态
!       5.           方程组的维数       
!----------------------------------------------------

implicit real*8(a-z)
!方程组的维数

integer::N,M=20    
real*8::y0(N),y(N),yt(N)

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
    
   
  end do 
yt=y

end subroutine RK4




subroutine func(f,t,y,N)
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

f(1)=y(2)
f(2)=-y(1)
f(3)=-y(3)

end subroutine func



end module m_adam3



  program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  采用Adams3阶方法计算微分方程组的主函数
!    
!-----------------------------------------------------
!  In put data  files :
!       1.
!       2.
!  Output data files  :
!       1.
!       2.
!-----------------------------------------------------

  use m_adam3
  
  implicit real*8(a-z)
  integer::N
  real*8::y(3),y0(3),x(3)

 
  open(unit=11,file='fout1.txt')
    
  write(11,101)
  101 format(/,T5,'采用Adams3方法计算方程组',/)
  
  N=3
  
  t0=0
  tt=0.5
  
  y0=(/1,0,-1/)
  
  call solve(t0,tt,y0,N)  
  
  
  end program main