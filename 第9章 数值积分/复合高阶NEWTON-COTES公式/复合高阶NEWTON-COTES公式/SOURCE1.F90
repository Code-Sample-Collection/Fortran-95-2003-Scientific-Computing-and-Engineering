module newton_cotes
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-5-31
!-----------------------------------------------------
!  Description :  复合5点Newton-Cotes公式计算定积分模块
!    
!-----------------------------------------------------
!  Contains    :
!      1.   复合方法函数
!      2.   单区间5点Newton-Cotes积分函数
!      3.   要计算的函数
!-----------------------------------------------------
!  Post Script :
!      1.
!      2. 
!-----------------------------------------------------

contains

subroutine solve(func,s,a,b,n)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  复合5点Newton-Cotes积分函数
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1. func 外部函数
!       2. a,b 积分区间
!       3. n 区间划分个数
!  Output parameters  :
!       1.  s 积分结果
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.   需要调用单区间5点Newton-Cotes公式函数
!       2.
!----------------------------------------------------
implicit real*8(a-z)
external func
integer::n,i

hstep=(b-a)/n

s=0
do i=1,n

    c=a+(i-1)*hstep
    d=a+i*hstep

   call cotes(func,s1,c,d)

   s=s+s1
end do

end subroutine solve


subroutine cotes(func,s1,c,d)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  单区间5点Newton-Cotes公式
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.   func  外部函数
!       2.   c,d 积分区间
!  Output parameters  :
!       1.
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.
!       2.
!----------------------------------------------------
implicit real*8(a-z)
external func

h=(d-c)/5d0

call func(f1,c)
call func(f2,c+h)
call func(f3,c+2d0*h)
call func(f4,c+3d0*h)   !d-2*h
call func(f5,c+4d0*h)   !d-h
call func(f6,d)

s1=19d0*f1+75d0*f2+50d0*f3+50d0*f4+75d0*f5+19d0*f6
s1=s1*5d0*h/288d0

end subroutine cotes


subroutine fun1(f,x)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-5-29
!-----------------------------------------------------
!  Purpose   :  需要计算的函数
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.     x  自变量
!       2. 
!  Output parameters  :
!       1.     f  因变量
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.
!       2.
!----------------------------------------------------
implicit real*8(a-z)

f=x**2+dsin(x)
end subroutine fun1


end module newton_cotes


program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-5-29
!-----------------------------------------------------
!  Purpose   :  复合5点Newton-Cotes方法主函数
!    
!-----------------------------------------------------
!  In put data  files :
!       1.
!       2.
!  Output data files  :
!       1.   result.txt计算结果
!       2.
!-----------------------------------------------------
!  Post Script :
!       1.  
!-----------------------------------------------------
use newton_cotes

implicit real*8(a-z)
integer::n

open(unit=11,file='result.txt')

write(11,101)

!划分8个区间
call solve(fun1,s,-2d0,2d0,8)

write(11,102)s

101 format(/,T5,'复合高阶Newton-Cotes法计算数值积分',/)
102 format(T5,'积分结果：',F15.8)


end program main
