module beta_chi
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-7-25
!-----------------------------------------------------
!  Description :   贝塔函数及Chi_square函数
!  
!----------------------------------------------------- 

contains

subroutine beta(fx,z,w)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-7-25
!-----------------------------------------------------
!  Purpose     : 计算贝塔函数
! 
!  Post Script :
!       1.
!       2.  需要调用伽马函数对数
!       3.
!
!-----------------------------------------------------
!  Input  parameters  :
!       1.   z
!       2.   w
!  Output parameters  :
!       1.  fx 贝塔函数值
!       2.
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------

implicit real*8(a-z)

call lngamma(fz,z)
call lngamma(fw,w)
call lngamma(fzw,z+w)

fx=dexp(fz+fw-fzw)

end subroutine beta




subroutine chi_sq(fx,x,v)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose     : Chi_squre函数
! 
!----------------------------------------------------

implicit real*8(a-z)
integer::v
call gp(fx,v/2d0,x/2d0)

end subroutine chi_sq


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
!            然后，用-Q
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


subroutine contfra(gammcf,a,x,gln)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.07.12
!-----------------------------------------------------
!  Purpose     : 用级数展开方法求P
! 
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

end module beta_chi


program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-7-25
!-----------------------------------------------------
!  Purpose     : Beta函数与Chi_squre函数主函数
!
!-----------------------------------------------------
use beta_chi

implicit real*8(a-z)

integer::i


open(unit=11,file='beta_chi.txt')



call beta(fb1,1.5d0,2.5d0)
call beta(fb2,3.4d0,7.2d0)
call beta(fb3,4.6d0,2.4d0)

write(11,101)fb1,fb2,fb3

101 format(/,T5,'Beta及Chi-square函数',//,&
           T3,'Beta函数',/,&
           T2,'Beta(1.5,2.5)=',F10.6,/,&
           T2,'Beta(2.3,7.2)=',F10.6,/,&
           T2,'Beta(4.6,2.4)=',F10.6,//,&
           T3,'Chi_square函数')

do i=1,5
call chi_sq(a,0.5d0,i)
call chi_sq(b,5d0,i)

write(11,103)i,a
write(11,104)i,b
end do
103 format(T3,'P(0.5,',I2,')=',F12.6)
104 format(T3,'P(5.0,',I2,')=',F12.6)

end program main