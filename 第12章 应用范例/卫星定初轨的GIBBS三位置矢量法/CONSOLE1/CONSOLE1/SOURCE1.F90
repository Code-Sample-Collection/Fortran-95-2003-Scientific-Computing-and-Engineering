module gibbs
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010.07.09
!-----------------------------------------------------
!  Description :   gibbs三位置矢量定初轨模块
!  
!  Post Script :
!      1.
!      2. 
!
!-----------------------------------------------------


contains

subroutine solve(v2,r1,r2,r3)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose     :  方法函数
! 
!  Post Script :
!       1.     需要用到矢量叉乘函数
!       2.
!       3.     注意定出的是v2    
!-----------------------------------------------------

implicit real*8(a-z)

real*8::r1(3),r2(3),r3(3),v2(3)

real*8::c12(3),c23(3),c31(3)

real*8::N(3),D(3),S(3)

real*8::Dr(3)

mu=398600

rou1=dsqrt(r1(1)**2+r1(2)**2+r1(3)**2)

rou2=dsqrt(r2(1)**2+r2(2)**2+r2(3)**2)

rou3=dsqrt(r3(1)**2+r3(2)**2+r3(3)**2)


call cross(c12,r1,r2)
call cross(c23,r2,r3)
call cross(c31,r3,r1)

N=rou1*c23+rou2*c31+rou3*c12

N_norm=dsqrt(n(1)**2+n(2)**2+n(3)**2)

D=c12+c23+c31

D_norm=dsqrt(d(1)**2+d(2)**2+d(3)**2)

S=r1*(rou2-rou3)+r2*(rou3-rou1)+r3*(rou1-rou2)


call cross(Dr,D,r2)

v2=dsqrt(mu/N_norm/D_norm)*(Dr/rou2+S)


end subroutine solve



subroutine cross(v,va,vb)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose     :   三维向量叉乘    V=Va * Vb
! 
!  Post Script :
!       1.
!       2.
!       3.    
!-----------------------------------------------------

implicit real*8(a-z)

real*8::v(3),va(3),vb(3)

v(1)=va(2)*vb(3)-va(3)*vb(2)
v(2)=va(3)*vb(1)-va(1)*vb(3)
v(3)=va(1)*vb(2)-va(2)*vb(1)

end subroutine cross


end module gibbs



program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose     :  主函数
!
!  Post Script :
!       1.    定出卫星v2 
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
use gibbs

implicit real*8(a-z)

real*8::r1(3),r2(3),r3(3),v2(3)


r1=(/-294.32d0,4265.1d0,5986.7d0/)

r2=(/-1365.d0,3637.6d0,6346.8d0/)

r3=(/-2940.3d0,2473.7d0,6555.8d0/)

call solve(v2,r1,r2,r3)



open(unit=11,file='result.txt')

write(11,101)v2

101 format(/,T5,'Gibbs 三位置矢量定初轨',///,&
                T5,'V2为(km/s)：',/,&
                3(/,T3,F12.8))
end program main