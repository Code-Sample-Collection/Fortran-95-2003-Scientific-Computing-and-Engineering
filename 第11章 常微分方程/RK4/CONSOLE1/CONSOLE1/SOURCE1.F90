
module m_rk4
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
!  Description :   RK-4经典方法计算微分方程
!    
!-----------------------------------------------------
!  Contains    :
!      1.    RK4
!      2.    RK4_2  传递函数版本  
!-----------------------------------------------------
!  Post Script :
!      1.
!      2. 
!-----------------------------------------------------

contains

  subroutine rk4(t0,tt,y0,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0          Vesion type 1
!  Coded by  :  syz 
!  Date      :  2010-4-11
!-----------------------------------------------------
!  Purpose   :  RK经典方法计算常微分方程
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.     t0 积分初时刻
!       2.     y0 初状态
!       3.     N  积分步数
!  Output parameters  :
!       1.
!
!----------------------------------------------------
!  Post Script :
!       1.    调用  func 子程序
!       2.
!----------------------------------------------------
  
  implicit real*8(a-z)
  
  integer::N,i
  
 !步长
  h=(tt-t0)/N
  
  t=t0
  y=y0
  
  do i=1,N
  
     call fun1(k1,t,y)
     call fun1(k2,t+h/2,y+h/2*k1)
     
     call fun1(k3,t+h/2,y+h/2*k2)
     call fun1(k4,t+h,y+h*k3)
  
     y=y+(k1+2*k2+2*k3+k4)*h/6
     
     t=t0+i*h
     
     write(11,101)i,t,y
   
  end do 
   
     101 format(t5,I4,2F15.8)
  end subroutine rk4
  
  
 subroutine rk4_2(func,t0,tt,y0,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0          Vesion type 2
!  Coded by  :  syz 
!  Date      :  2010-4-11
!-----------------------------------------------------
!  Purpose   :  RK经典方法计算常微分方程
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.     t0 积分初时刻
!       2.     y0 初状态
!       3.     N  积分步数
!  Output parameters  :
!       1.
!
!----------------------------------------------------
!  Post Script :
!       1.    调用  func 子程序
!       2.
!----------------------------------------------------
 
 implicit real*8(a-z)
 
  !声明func是外部函数，这句不能省略
  !把外部函数当初参数传递
 external func   
  integer::N,i
 
  h=(tt-t0)/N
  
  t=t0
  y=y0
 
   do i=1,N
  
     call func(k1,t,y)
     call func(k2,t+h/2,y+h/2*k1)
     
     call func(k3,t+h/2,y+h/2*k2)
     call func(k4,t+h,y+h*k3)
  
     y=y+(k1+2*k2+2*k3+k4)*h/6
     
     t=t0+i*h
     
     write(12,101)i,t,y
   
  end do 
  101 format(T5,I4,2F15.8)
 end subroutine rk4_2 
   
   subroutine fun1(f,t,y)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   : 待计算的函数  其中 f=(t,y)
!    
!-----------------------------------------------------
   implicit real*8(a-z)   
    f=-y+t**2+1
   end subroutine fun1   
   
   
   subroutine fun2(f,t,y)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   : 待计算的函数  其中 f=(t,y)
!    
!-----------------------------------------------------
   implicit real*8(a-z)
   
   f=2/t*y+t**2*exp(t)
   
   end subroutine fun2
   

end module m_rk4


  program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   : RK方法主程序
!    
!-----------------------------------------------------
!  Output data files  :
!       1.    fout1.txt
!       2.    fout2.txt
!-----------------------------------------------------
!  Post Script :
!       1.  提供两个版本
!       2.
!-----------------------------------------------------
  
  use m_rk4
  
  implicit real*8(a-z)
 
  integer::N=20
 
  open(unit=11,file='fout1.txt')
  open(unit=12,file='fout2.txt')
  
  t0=0
  t1=1
  y0=1
  
  write(11,101)
  101 format(/,T5,'采用版本1计算函数1结果为：',/)
  
  call rk4(t0,t1,y0,N)  
 
 
 write(12,102)
 102 format(/,T5,'采用版本2计算函数1结果为：',/)
  call rk4_2(fun1,t0,t1,y0,N)
  
  
  write(12,103)
  103 format(/,T5,'函数2计算结果为：',/)
  
 t0=1
 t1=2
 y0=0
  
  call rk4_2(fun2,t0,t1,y0,n)
  
  end program main
