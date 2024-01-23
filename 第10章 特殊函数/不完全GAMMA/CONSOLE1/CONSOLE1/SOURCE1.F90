module incompGamma
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010.07.12
!-----------------------------------------------------
!  Description :   不完全Gamma函数及其互补形式
!  
!  Post Script :
!      1.
!      2. 
!
!-----------------------------------------------------
!  Contains    :
!      1.     gp---求P
!      2.     gq---求Q
!-----------------------------------------------------
!  Parameters  :
!      1.
!      2.
!----------------------------------------------------- 

contains

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

end module incompGamma

program  main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose     : 不完全Gamma函数及其互补形式
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

use incompGamma
implicit real*8(a-z)

open(unit=11,file='result.txt')

call gp(a,1d0,1d0)
call gq(b,1d0,0.1d0)

write(11,101)a,b

101 format(/,T5,'不完全Gamma函数及其互补形式',//,&
            T3,'计算结果为：',/,2(/F16.8))
end program main