module m_gill
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
!  Description :   Gill方法模块
!                     
!----------------------------------------------------- 
!  Contains    :
!      1.      gill---    Gill 方法函数
!      2.      func---    求解的函数
!-----------------------------------------------------


contains

subroutine solve(func,t0,tt,y0,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  Gill方法求解微分方程
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.         func为外部子程序
!       2. 
!  Output parameters  :
!----------------------------------------------------

implicit real*8(a-z)

!声明func为外部子程序
external func
integer::N,i

h=(tt-t0)/N

t=t0
y=y0

!先把根号2算出来 ，这样不用每递推一步都要进行开方，而每
! 递推一步要用到6次根号2
sq2=dsqrt(2d0)


do i=1,N

    call func(k1,t,y)
    call func(k2,t+h/2,y+h*k1/2)
    call func(k3,t+h/2,y+(sq2-1)/2d0*h+(1-sq2/2)*h*k2)
    call func(k4,t+h/2,y-sq2/2*h*k2+(1+sq2/2)*h*k3)

    y=y+(k1+(2-sq2)*k2+(2+sq2)*k3+k4)*h/6
    
    t=t+h
    
    write(11,101)i,t,y
    
end do 

101 format(5x,I4,2F15.8)

end subroutine solve


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
    f=0.5d0*(1+t)*y**2
   end subroutine fun1  



end module m_gill


  program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-12
!-----------------------------------------------------
!  Purpose   : Gill方法主程序
!    
!-----------------------------------------------------
!  Output data files  :
!       1.    fout.txt
!-----------------------------------------------------
!  Post Script :
!       1.  
!       2.
!-----------------------------------------------------
  
  use m_gill
  
  implicit real*8(a-z)
 
  integer::N=20
 
  open(unit=11,file='fout1.txt')
  
  t0=0
  t1=1
  y0=1
  
  write(11,101)
  101 format(/,T5,'采用Gill方法计算微分方程：',/)
  
  call solve(fun1,t0,t1,y0,N)  
  
  
  end program main