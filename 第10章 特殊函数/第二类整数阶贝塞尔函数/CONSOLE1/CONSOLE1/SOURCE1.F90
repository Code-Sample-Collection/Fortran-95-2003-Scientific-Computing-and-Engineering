module bessy
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-7-25
!-----------------------------------------------------
!  Description :   第二类整数阶贝塞尔函数模块
!  
!-----------------------------------------------------
!  Contains    :
!      1.     零一阶第二类贝塞尔函数
!      2.     零一阶第一类贝塞尔函数
!-----------------------------------------------------

contains

subroutine solve(fx,x,n)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-7-25
!-----------------------------------------------------
!  Purpose     : 任意阶第二类贝塞尔函数
! 
!  Post Script :
!       1.      由0，1阶递推
!       2.      而0，1阶级需要调用0,1阶第一类函数
!       3.
!
!-----------------------------------------------------
!  Input  parameters  :
!       1.    x  自变量
!       2.    n  阶数    n=0,1,2,....
!  Output parameters  :
!       1.
!       2.
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------

implicit real*8(a-z)

integer::j,n


!如果N=0，则直接调用0阶函数
!调用结束后返回
if(n==0)then
call bessy0(fx,x)
return
end if

!如果N=1,则调用1阶函数
!调用结束后返回
if(n==1)then
call bessy1(fx,x)
return
end if

call bessy0(by0,x)
call bessy1(by1,x)

!采用递推公式
do j=1,n-1
  byn=j*by*2d0/x-by0
  by0=by1
  by1=byn
end do
fx=by1
end subroutine solve


subroutine bessy0(fx,x)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose     : 0阶第二类贝塞尔函数
! 
!  Post Script :
!       1.      需要调用0阶第一类贝塞尔函数
!       2.
!       3.
!
!-----------------------------------------------------
!  Input  parameters  :
!       1.  x 自变量
!       2. 
!  Output parameters  :
!       1.  fx 函数值
!       2.
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------


implicit real*8(a-z)

data p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,&
     -.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/&
	 -.1562499995d-1,.1430488765d-3,-.6911147651d-5,&
	 .7621095161d-6,-.934945152d-7/
data r1,r2,r3,r4,r5,r6/-2957821389.d0,7062834065.d0,&
     -512359803.6d0,10879881.29d0,-86327.92757d0,&
	 228.4622733d0/,s1,s2,s3,s4,s5,s6/40076544269.d0,&
	 745249964.8d0,7189466.438d0,47447.26470d0,&
	 226.1030244d0,1.d0/
if(x<8d0)then
  y=x*x
  
! 调用0阶第一类贝塞尔函数值
  call bessj0(tmp,x)
  fx=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/&
         (s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))+&
		 .636619772*tmp*dlog(x)
else
  z=8d0/x
  y=z*z
  xx=x-0.785398164d0
  fx=dsqrt(.636619772/x)*(dsin(xx)*(p1+y*(p2+y*&
         (p3+y*(p4+y*p5))))+z*dcos(xx)*(q1+y*(q2+y*&
		 (q3+y*(q4+y*q5)))))
end if
end subroutine bessy0

subroutine bessy1(fx,x)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-7-25
!-----------------------------------------------------
!  Purpose     :  1阶第二类贝塞尔函数
! 
!  Post Script :
!       1.   需要调用1阶第一类贝塞尔函数值
!-----------------------------------------------------
!  Input  parameters  :
!       1.  x 自变量
!       2. 
!  Output parameters  :
!       1. fx 函数值
!       2.
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------

implicit real*8(a-z)
SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,&
     s1,s2,s3,s4,s5,s6,s7
DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,&
     .2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/&
	 .04687499995d0,-.2002690873d-3,.8449199096d-5,&
	 -.88228987d-6,.105787412d-6/
DATA r1,r2,r3,r4,r5,r6/-.4900604943d13,.1275274390d13,&
     -.5153438139d11,.7349264551d9,-.4237922726d7,&
	 .8511937935d4/,s1,s2,s3,s4,s5,s6,s7/.2499580570d14,&
	 .4244419664d12,.3733650367d10,.2245904002d8,&
	 .1020426050d6,.3549632885d3,1.d0/
if(x<8.)then
  y=x**2
  call bessj1(tmp,x)
  fx=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/&
         (s1+y*(s2+y*(s3+y*(s4+y*(s5+y*(s6+y*s7))))))&
		 +.636619772*(tmp*dlog(x)-1./x)
else
  z=8./x
  y=z**2
  xx=x-2.356194491
  fx=dsqrt(.636619772/x)*(dsin(xx)*(p1+y*(p2+y*&
         (p3+y*(p4+y*p5))))+z*dcos(xx)*(q1+y*(q2+y*&
		 (q3+y*(q4+y*q5)))))
endif
end subroutine bessy1

subroutine bessj0(fx,x)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-7-23
!-----------------------------------------------------
!  Purpose     : 0 阶第一类贝塞尔函数
! 
!  Post Script :
!       1.
!       2.     0阶函数可以单独计算0阶结果
!       3.     这里为任意阶函数所调用
!
!-----------------------------------------------------
!  Input  parameters  :
!       1.   x 自变量
!       2. 
!  Output parameters  :
!       1.   fx 0阶计算结果
!       2.
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------
implicit real*8(a-z)

DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,&
     -.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5&
	 /-.1562499995d-1,.1430488765d-3,-.6911147651d-5,&
	 .7621095161d-6,-.934945152d-7/
DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,&
     651619640.7d0,-11214424.18d0,77392.33017d0,&
	 -184.9052456d0/,s1,s2,s3,s4,s5,s6/57568490411.d0,&
	 1029532985.d0,9494680.718d0,59272.64853d0,&
	 267.8532712d0,1.d0/
if(dabs(x)<8d0) then
  y=x**2
  fx=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/&
         (s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
else
  ax=dabs(x)
  z=8d0/ax
  y=z**2
  xx=ax-0.785398164d0
  fx=dsqrt(.636619772/ax)*(dcos(xx)*(p1+y*(p2+y*&
         (p3+y*(p4+y*p5))))-z*dsin(xx)*(q1+y*(q2+y*&
		 (q3+y*(q4+y*q5)))))
endif
end subroutine bessj0

subroutine bessj1(fx,x)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-7-23
!-----------------------------------------------------
!  Purpose     : 1阶第一类贝塞尔函数
! 
!  Post Script :
!       1.     可以单独计算1阶函数结果
!       2.     这里为任意阶函数所调用
!       3.
!
!-----------------------------------------------------
!  Input  parameters  :
!       1.  x自变量
!       2. 
!  Output parameters  :
!       1. fx 1阶计算结果
!       2.
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------

implicit real*8(a-z)

data r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,&
     242396853.1d0,-2972611.439d0,15704.48260d0,&
	 -30.16036606d0/,s1,s2,s3,s4,s5,s6/144725228442.d0,&
	 2300535178.d0,18583304.74d0,99447.43394d0,&
	 376.9991397d0,1.d0/
data p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,&
     .2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/&
	 .04687499995d0,-.2002690873d-3,.8449199096d-5,&
	 -.88228987d-6,.105787412d-6/
if(dabs(x)<8.) then
  y=x**2
  fx=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/&
         (s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
else
  ax=dabs(x)
  z=8d0/ax
  y=z**2
  xx=ax-2.356194491d0
  fx=dsqrt(0.636619772d0/ax)*(dcos(xx)*(p1+y*(p2+y*&
         (p3+y*(p4+y*p5))))-z*dsin(xx)*(q1+y*(q2+y*&
		 (q3+y*(q4+y*q5)))))*sign(1d0,x)
endif
end subroutine bessj1

end module bessy

program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-7-23
!-----------------------------------------------------
!  Purpose     :  第二类整数阶贝塞尔函数主函数
!
!  Post Script :
!       1.
!       2.
!    
!-----------------------------------------------------
!  In put data  files :
!       1.
!       2.
!  Output data files  :
!       1.
!       2.
!
!-----------------------------------------------------
use bessy
implicit real*8(a-z)

integer::i

open(unit=11,file='result.txt')

write(11,100)

do i=0,10
call solve(a,0.8d0,i)
call solve(b,5.6d0,i)
call solve(c,37.d0,i)

write(11,101)i,a
write(11,102)i,b
write(11,103)i,c

end do

100 format(/,T10,'第二类整数阶贝塞尔函数',/)

101 format(T2,'i=',I3,3x,'x=0.8',3x,'Bessj(x,i)=',E16.6)
102 format(T2,'i=',I3,3x,'x=5.6',3x,'Bessj(x,i)=',E16.6)
103 format(T2,'i=',I3,3x,'x=3 7',3x,'Bessj(x,i)=',E16.6)

end program main
