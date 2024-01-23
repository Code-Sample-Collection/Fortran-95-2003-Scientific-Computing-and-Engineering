module err_norm
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010.07.12
!-----------------------------------------------------
!  Description :   误差函数，互补形式，正态分布函数
!  
!  Post Script :
!      1.
!      2. 
!
!-----------------------------------------------------
!  Contains    :
!      1.     
!      2.     
!-----------------------------------------------------
!  Parameters  :
!      1.
!      2.
!----------------------------------------------------- 

contains

subroutine err(fx,x)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-7-15
!-----------------------------------------------------
!  Purpose     : 误差函数
! 
!  Post Script :
!       1.
!       2.
!       3.
!
!-----------------------------------------------------
!  Input  parameters  :
!       1.
!       2. 
!  Output parameters  :
!       1.
!       2.
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------


implicit real*8(a-z)

if (x<0d0) then
 call gp(tmp1,0.5d0,x*x)
 fx=-tmp1
else
 call gp(tmp1,0.5d0,x*x)
 fx=tmp1
end if
end subroutine err 

subroutine errc(fx,x)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010年7月15日
!-----------------------------------------------------
!  Purpose     : 误差函数的互补形式
! 
!  Post Script :
!       1.
!       2.
!       3.
!
!-----------------------------------------------------
!  Input  parameters  :
!       1.
!       2. 
!  Output parameters  :
!       1.
!       2.
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------


implicit real*8(a-z)

if(x<0d0) then
 call gp(tmp1,0.5d0,x*x)
 fx=1+tmp1
else
 call gq(tmp1,0.5d0,x*x)
 fx=tmp1
end if
end subroutine errc


subroutine  normal(fx,x)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-7-15
!-----------------------------------------------------
!  Purpose     : 标准正态分布
! 
!  Post Script :
!       1.
!       2.
!       3.
!
!-----------------------------------------------------
!  Input  parameters  :
!       1.
!       2. 
!  Output parameters  :
!       1.
!       2.
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------

implicit real*8(a-z)

t=x/dsqrt(2d0)

call err(fx,t)

fx=0.5d0+0.5d0*fx

end subroutine normal


subroutine gp(fx,a,x)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.07.12
!-----------------------------------------------------
!  Purpose     : 求P（afa,x）
! 
!  Post Script :
!       1.   当x<a+1时采用级数展开式计算
!       2.    
!       3.   当不满足x<a+1时，实际上是直接计算互补函数Q
!            然后，用1-Q
!            互补函数的计算是用连分式逼近 
!-----------------------------------------------------
!  Input  parameters  :
!       1.
!       2. 
!  Output parameters  :
!       1.
!       2.
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------

implicit real*8(a-z)

if(x<a+1d0)then
  call series(fx1,a,x,gln)
  fx=fx1
else
  call contfra(fx2,a,x,gln)
  fx=1d0-fx2
endif

end subroutine gp



subroutine gq(fx,a,x)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.07.12
!-----------------------------------------------------
!  Purpose     : 计算Q
! 
!  Post Script :
!       1.   当x<a+1时候，计算Q的互补形式P，然后用1-P
!       2.   计算P是用级数展开的方式
!       3.   不满足x<a+1时，采用连分式直接计算Q 
!-----------------------------------------------------
!  Input  parameters  :
!       1.
!       2. 
!  Output parameters  :
!       1.
!       2.
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------
implicit real*8(a-z)
if(x<a+1d0) then
  call series(fx1,a,x,gln)
  fx=1d0-fx1
else
  call contfra(fx2,a,x,gln)
  fx=fx2
endif
end subroutine gq

subroutine contfra(gammcf,a,x,gln)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.07.12
!-----------------------------------------------------
!  Purpose     : 用级数展开方法求P
! 
!  Post Script :
!       1.
!       2.
!       3.    
!-----------------------------------------------------
!  Input  parameters  :
!       1.
!       2. 
!  Output parameters  :
!       1.
!       2.
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------
implicit real*8(a-z)

integer:: i,Imax

Imax=100
eps=3d-7
fpmin=1d-30

call lngamma(gln,a)
b=x+1d0-a
c=1d0/fpmin
d=1d0/b
h=d
do i=1,Imax
  an=-i*(i-a)
  b=b+2d0
  d=an*d+b
  if(abs(d)<fpmin) d=fpmin
  c=b+an/c
  if(abs(c)<fpmin) c=fpmin
  d=1d0/d
  del=d*c
  h=h*del
  if(abs(del-1d0)<EPS) then
    gammcf=exp(-x+a*log(x)-gln)*h
    return
  endif
end do

end subroutine contfra

subroutine series(gamser,a,x,gln)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.07.12
!-----------------------------------------------------
!  Purpose     :  用连分式逼近Q
! 
!  Post Script :
!       1.
!       2.
!       3.    
!-----------------------------------------------------
!  Input  parameters  :
!       1.
!       2. 
!  Output parameters  :
!       1.
!       2.
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------
implicit real*8(a-z)
integer::Imax,n

Imax=100
eps=3d-7

call lngamma(gln,a)

if(x<=0.) then
  if(x<0d0) pause 'x < 0 in gser'
  gamser=0d0
  return
endif
ap=a
sum=1d0/a
del=sum
do n=1,Imax
  ap=ap+1d0
  del=del*x/ap
  sum=sum+del
  if(abs(del)<abs(sum)*EPS) then
    gamser=sum*exp(-x+a*log(x)-gln)
	return
  end if  
end do

end subroutine series 



subroutine lngamma(ln_gamma,x0)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.07.12  
!-----------------------------------------------------
!  Purpose     : 伽马函数的对数
! 
!----------------------------------------------------

implicit real*8(a-z) 

integer:: i
real*8::cof(6)

 cof=(/76.18009172947146d0,-86.50532032941677d0,&
     24.01409824083091d0,-1.231739572450155d0,&
	 0.1208650973866179d-2,-.5395239384953d-5/)

temp1=dsqrt(2*3.141592653589793d0)

x=x0
y=x
tmp=x+5.5d0
tmp=(x+0.5d0)*dlog(tmp)-tmp
ser=1.000000000190015d0
do i=1,6
  y=y+1.d0
  ser=ser+cof(i)/y
end do
ln_gamma=tmp+dlog(temp1*ser/x)
end subroutine lngamma

end module err_norm

program  main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose     : 误差函数，互补形式，及标准正态分布表
!                的制作
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

use err_norm
implicit real*8(a-z)

real*8::A(0:30,0:9),FX3(0:30,0:9)

integer::i,j

open(unit=11,file='result_err.txt')

open(unit=12,file='normal.txt')

call err(fx1,0.5d0)
call errc(fx2,0.6d0)

do i=0,30
  do j=0,9
    
   a(i,j)=i/10d0+j/100d0
    
   call normal(tmp1,a(i,j))

   fx3(i,j)=tmp1

  end do
end do

write(11,100)fx1,fx2


write(12,101)
write(12,102)((fx3(i,j),j=0,9),i=0,30)

100 format(/,t5,'误差函数及其互补形式计算结果',//,&
            T1,2F16.7)

101 format(/,T28,'正态分布表',//)
102 format(31(10F7.4,/))


end program main