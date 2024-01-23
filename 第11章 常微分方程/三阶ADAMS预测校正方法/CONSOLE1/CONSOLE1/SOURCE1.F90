module  PECE3
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-4-13
!-----------------------------------------------------
!  Description :   三阶PECE方法计算微分法方程组模块
!    
!----------------------------------------------------- 
!  Contains    :
!      1.      三阶PECE方法函数
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
!  Purpose   :三阶Adams PECE方法计算微分法方程组方法函数
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
real*8::ya(N),y0(N),y1(N),y2(N),y3(N)
real*8::f0(N),f1(N),f2(N),f3(N)
real*8::y3_p(N),f3_e(N)

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


!----输出前几个结果,实际为RK法启动
 write(11,101)t0,y0
 101 format(T3,<N+1>f10.5)
 write(11,102)t1,y1
 102 format(T3,<N+1>f10.5)
  write(11,103)t2,y2
 103 format(T3,<N+1>f10.5)

!------------------ 

do i=3,M
    
    !P步:预测
     y3_p=y2+h/12*(23*f2-16*f1+5*f0) 
     
    
     t3=t2+h
     
     !E步：计算f
     call func(f3_e,t4,y3_p,N)
     
     
     !C步：计算下一点的值
     y3=y2+h/24d0*(9*f3_e+19*f2-5*f1+f0)
     
     !E:计算f3
    call func(f3,t3,y3,N)
   
     
     !输出计算结果
     write(11,105)t3,y3
     
     !  以下为各量更新
     t2=t2+h
     
	 y2=y3
	  
     f0=f1
     f1=f2
     f2=f3
    
end do
  
105 format(T3,<N+1>f10.5)


end subroutine solve



subroutine  RK4(t0,tt,y0,yt,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-13
!-----------------------------------------------------
!  Purpose   :  修改后的阶经典RK方法计算微分方程组方法函数
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

implicit real*8(a-z)
integer::N
real*8::f(N),y(N)

f(1)=y(2)*y(3)
f(2)=-y(1)*y(3)
f(3)=-0.4*y(1)*y(2)


end subroutine func


end module pece3



  program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  采用3阶Adams PECE方法计算微分方程组的主函数
!    
!-----------------------------------------------------
!  In put data  files :
!       1.
!       2.
!  Output data files  :
!       1.
!       2.
!-----------------------------------------------------

  use pece3
  
  implicit real*8(a-z)
  integer,parameter::N=3
  real*8::y(N),y0(N)

 
  open(unit=11,file='fout1.txt')
    
  write(11,101)
  101 format(/,T5,'3阶Adams 预测校正方法计算方程组',/)
 
  
  t0=1
  tt=6
  
  y0=(/1d0,1d0,1d0/)
  
  call solve(t0,tt,y0,N)  
  
end program main
