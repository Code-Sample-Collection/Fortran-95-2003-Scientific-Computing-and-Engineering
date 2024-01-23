module orbit
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
!  Description :   二体问题轨道外推模块
!    
!-----------------------------------------------------
!  Parameters  :
!      1.
!      2.
!----------------------------------------------------- 
!  Contains    :
!      1.       右函数
!      2.       经典的RK4积分器
!-----------------------------------------------------
!  Post Script :
!      1.
!      2. 
!-----------------------------------------------------
contains

subroutine  solve(t0,tt,y0)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :      4阶经典RK方法计算微分方程组方法函数
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.         
!       2.           t0积分初时刻 
!       3.           tt积分结束时刻
!       4.           y0 初状态    
!  Output parameters  :
!       1.
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.    计算结果保存在文件中，而没用用参数带出
!       2.
!----------------------------------------------------

implicit real*8(a-z)



real*8::y0(6),y(6)

real*8::k1(6),k2(6),k3(6),k4(6)

h=60d0  !积分步长60秒

t=t0
y=y0

do 

     call twobody(k1,y)
     call twobody(k2,y+h/2*k1)
     
     call twobody(k3,y+h/2*k2)
     call twobody(k4,y+h*k3)
     
     y=y+(k1+2*k2+2*k3+k4)*h/6
     
     t=t+h
     
     write(11,101)t,y/1000d0
    
     if (t>=tt) exit
    
 end do 
 
  101 format(F8.1,6F15.4)

end subroutine solve


subroutine twobody(f,x)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.05.17
!-----------------------------------------------------
!  Purpose   :  二体问题右函数
!               
!-----------------------------------------------------
!  Input  parameters  :
!       1.  x 状态向量
!       2.  t 时间
!  Output parameters  :
!       1.  f 右函数
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.  可以看到f 仅仅是x的函数  而与t无关
!       2.
!----------------------------------------------------
implicit real*8(a-z)

real*8::f(6),x(6)

real*8,parameter::gm=3.98600448D14 
!引力常数  单位：m**3/sec**2

r2=x(1)**2+x(2)**2+x(3)**2
r1=dsqrt(r2)
!r1 表示向量长度，r2表示平方

!速度向量
f(1)=x(4)
f(2)=x(5)
f(3)=x(6)

gm1=gm/(r1*r2)

!加速度向量
f(4)=-gm1*x(1)
f(5)=-gm1*x(2)
f(6)=-gm1*x(3)

end subroutine twobody


end module orbit


program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   : 二体问题轨道外推主函数
!             
!-----------------------------------------------------
!  In put data  files :
!       1.
!       2.
!  Output data files  :
!       1. report.txt 记录轨道数据文件
!       2.
!-----------------------------------------------------
!  Post Script :
!       1.   
!       2.  轨道根数初值为JASON卫星 1 Jul 2007 12:00:00.000状态
!       3.  输出节点为每间隔60秒
!-----------------------------------------------------
use orbit
implicit real*8(a-z)
real*8::x0(6)

open(unit=11,file='report.txt')

!单位km， km/s
x0=(/6678.137000d0,-0.000000d0, -0.000000d0,&
      -0.000000d0, 6.789530d0,  3.686414d0/)

write(11,101)
write(11,102)t0,x0

!化为标准单位(m,m/s)
x0=x0*1000d0

t0=0   
!单位s
tt=360*60
!单位s,即外推6小时共360分钟


!调用积分器
call solve(t0,tt,x0)

 101 format(/,T38,'航天器轨道外推',//,&
            T2,'(单位：秒，千米，千米/秒)')
 102 format(F8.1,6F15.4)
end program main
