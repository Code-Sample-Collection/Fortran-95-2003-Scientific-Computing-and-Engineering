module bessj
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-7-23
!-----------------------------------------------------
!  Description :   第一类整数阶贝塞尔函数模块
!  
!  Post Script :
!      1.     入口函数solve 中包含0阶和1阶情况
!      2. 
!
!-----------------------------------------------------
!  Contains    :
!      1.  0 阶函数   
!      2.  1 阶函数     
!      3.  任意阶函数
!-----------------------------------------------------
!  Parameters  :
!      1.
!      2.
!----------------------------------------------------- 
contains

subroutine solve(fx,x,n)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose     : 第一类任意整数阶贝塞尔函数
! 
!  Post Script :
!       1.    计算方法实际是根据0阶和1阶而得出任意阶函数值
!       2.
!       3.
!
!-----------------------------------------------------
!  Input  parameters  :
!       1.     x 自变量
!       2.     n 阶数（n=0,1,2,3,.....）
!  Output parameters  :
!       1.     fx 贝塞尔函数值
!       2.
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------

implicit real*8(a-z)
integer:: n,index1


integer:: j,sum_j,m

! 如果输入是0阶，直接调用0阶函数
!调用结束后返回
  if (n==0) then
    call bessj0(bj0,x)
    fx=bj0
	return
  end if
!如果输入是1阶，直接调用1阶贝塞尔函数
!调用结束后返回
  if (n==1) then
   call bessj1(bj1,x)
   fx=bj1
   return
  end if

!以下计算N>=2的情况
index1=40
b_0=1d10
b_1=1d-10

ax=dabs(x)
if(ax==0d0) then
  fx=0d0
else if(ax>float(n)) then
  tox=2d0/ax
  call bessj0(bjm,ax)
  call bessj1(bj,ax)


 
!n>=2的情况
  do j=1,n-1
    bjp=j*tox*bj-bjm
    bjm=bj
    bj=bjp
  end do
  fx=bj
else
  tox=2d0/ax
  m=2*((n+int(sqrt(float(index1*n))))/2)
  fx=0d0
  sum_j=0d0
  sum1=0d0
  bjp=0d0
  bj=1d0
  do j=m,1,-1
    bjm=j*tox*bj-bjp
    bjp=bj
    bj=bjm
    if(dabs(bj)>b_0) then
      bj=bj*b_1
      bjp=bjp*b_1
      fx=fx*b_1
      sum1=sum1*b_1
    endif
    if(sum_j/=0) sum1=sum1+bj
    sum_j=1-sum_j
    if(j==n) fx=bjp
  end do
  sum1=2d0*sum1-bj
  fx=fx/sum1
endif

if(x<0d0.and.mod(n,2)==1) then
fx=-fx
end if
end subroutine solve


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

data p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,&
     -.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5&
	 /-.1562499995d-1,.1430488765d-3,-.6911147651d-5,&
	 .7621095161d-6,-.934945152d-7/
data r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,&
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

end module bessj


program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-7-23
!-----------------------------------------------------
!  Purpose     :  第一类整数阶贝塞尔函数主函数
!
!  Post Script :
!       1.
!       2.
!    
!-----------------------------------------------------

use bessj

implicit real*8(a-z)
integer::i

open(unit=11,file='result.txt')

write(11,100)

do i=0,10
   call solve(a,0.5d0,i)
   call solve(b,5d0,i)
   call solve(c,50d0,i)

write(11,101)i,a
write(11,102)i,b
write(11,103)i,c

end do

100 format(/,T10,'第一类整数阶贝塞尔函数',/)

101 format(T2,'i=',I3,3x,'x=0.5',3x,'Bessj(x,i)=',F11.8)
102 format(T2,'i=',I3,3x,'x=5.0',3x,'Bessj(x,i)=',F11.8)
103 format(T2,'i=',I3,3x,'x=5 0',3x,'Bessj(x,i)=',F11.8)

end program main