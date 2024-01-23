module gammlog
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
contains
subroutine solve(ln_gamma,x0)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.07.12  
!-----------------------------------------------------
!  Purpose     : 伽马函数的对数
! 
!  Post Script :
!       1.
!       2.
!       3.    
!-----------------------------------------------------
!  Input  parameters  :
!       1.    x0   自变量
!       2. 
!  Output parameters  :
!       1.   ln_gamma  伽马函数的对数
!       2.
!  Common parameters  :
!       1.
!       2.
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
end subroutine solve

end module gammalog


program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.07.12
!-----------------------------------------------------
!  Purpose     : 计算Gamma函数主程序
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
use gammalog
implicit real*8(a-z)


integer::i

open(unit=11,file='result.txt')

write(11,101)

do i=1,10
    
    x=i/10d0
        
    call solve(y,x)
   
   write(11,102)x,dexp(y)
   
end do

101 format(/,T12,'Gamma函数',//,&
            T9,'x',T18,'Gamma(x)',/)
102 format(T3,2F12.7)


end program main

