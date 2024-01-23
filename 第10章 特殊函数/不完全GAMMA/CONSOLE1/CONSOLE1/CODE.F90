subroutine imcopGamma(gm,x,a,in_dex)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.07.12
!-----------------------------------------------------
!  Purpose     : 计算不完全Gamma函数或其互补形式
! 
!  Post Script :
!       1.     当In_dex=1时，返回P（afa,x）
!       2.     否则计算Q(afa,x)
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

integer::in_dex

if(in_dex==1)  then
	
	call gp(gm,a,x)
else then
	call gq(gm,a,x)
end  if

end subroutine impcopGamma

subroutine gp(gammp,a,x)
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
  call gser(gamser,a,x,gln)
  gammp=gamser
else
  call gcf(gammcf,a,x,gln)
  gammp=1d0-gammcf
endif

end subroutine gp



subroutine gq(gammq,a,x)
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
  call gser(gamser,a,x,gln)
  gammq=1d0-gamser
else
  call gcf(gammcf,a,x,gln)
  gammq=gammcf
endif
end subroutine gq

subroutine gcf(gammcf,a,x,gln)
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
b=x+1.-a
c=1./epmin
d=1./b
h=d
do i=1,Imax
  an=-i*(i-a)
  b=b+2.
  d=an*d+b
  if(abs(d)<fpmin) d=fpmin
  c=b+an/c
  if(abs(c)<fpmin) c=fpmin
  d=1./d
  del=d*c
  h=h*del
  if(abs(del-1.)<EPS) then
    gammcf=exp(-x+a*log(x)-gln)*h
    return
  endif
end do

end subroutine gcf

subroutine gser(gamser,a,x,gln)
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
  if(x<0.) pause 'x < 0 in gser'
  gamser=0.
  return
endif
ap=a
sum=1./a
del=sum
do n=1,Imax
  ap=ap+1.
  del=del*x/ap
  sum=sum+del
  if(abs(del)<abs(sum)*EPS) then
    gamser=sum*exp(-x+a*log(x)-gln)
	return
  end if  
end do

end subroutine gser 



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

