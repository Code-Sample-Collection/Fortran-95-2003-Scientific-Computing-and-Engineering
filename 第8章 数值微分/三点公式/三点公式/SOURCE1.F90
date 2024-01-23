module therepoint
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010.08.11
!-----------------------------------------------------
!  Description :   三点公式计算微分模块
!  
!-----------------------------------------------------

contains
subroutine solve()
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.08.04
!-----------------------------------------------------
!  Purpose     : 三点公式求微分
! 
!  Post Script :
!       1.
!       2.
!       3.    
!-----------------------------------------------------


implicit real*8(a-z)

integer::i

real*8::f(-2:2)


h=0.1d0

!建立文件用以存放离散数据
open(unit=11,file='file1.txt')

write(11,101)


do i=-2,2

  f(i)=dsin(3d0+i*h)+dcos(3d0+i*h)
!记录下离散的函数数据，写入文件1
  write(11,102)3d0+i*h,f(i)
end do


!以上为产生模拟数据，以下为三点公式计算微分

!公式1
df1=-3d0*f(0)+4d0*f(1)-f(2)
df1=df1/(2d0*h)

!公式2
df2=f(-2)-4d0*f(-1)+3d0*f(0)
df2=df2/(2d0*h)

!公式3
df3=(f(1)-f(-1))/(2d0*h)

!新建文件用以保存数值微分结果
open(unit=12,file='file2.txt')

write(12,103)
write(12,*)'公式1计算结果：',df1
write(12,*)'公式2计算结果：',df2
write(12,*)'公式3计算结果：',df3

write(12,104)dcos(3d0)-sin(3d0)


101 format(/,T14,'离散数据为：',/)
102 format(2(F16.10))
103 format(/,T14,'三点微分公式',/)
104 format('实际结果应该为：',F16.12)
end subroutine solve

end module therepoint

program main

!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.08.04 
!-----------------------------------------------------
!  Purpose     : 三点公式求微分
!
!  Post Script :
!       1.
!       2.
!    
!-----------------------------------------------------
use therepoint
call solve
end



