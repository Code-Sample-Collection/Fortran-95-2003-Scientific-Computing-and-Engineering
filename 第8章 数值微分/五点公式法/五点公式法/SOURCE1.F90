module fivediff
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010.08.11
!-----------------------------------------------------
!  Description :   
!  
!  Post Script :
!      1.
!      2. 
!
!-----------------------------------------------------

contains 
subroutine solve
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.08.04
!-----------------------------------------------------
!  Purpose     :  五点公式计算微分
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

integer i

real*8 f(-4:4)


h=0.1d0

!建立文件用以存放离散数据
open(unit=file1,file='file1.txt')

do i=-4,4

  f(i)=dsin(3d0+i*h)+dcos(3d0+i*h)
!记录下离散的函数数据，写入文件
  write(file1,*)3d0+i*h,f(i)
end do


!以上为产生模拟数据，以下为三点公式计算微分


!公式1
df1=-25*f(0)+48*f(1)-36*f(2)+16*f(3)-3*f(4)
df1=df1/12/h


!公式2
df2=f(-2)-8*f(-1)+8*f(1)-f(2)
df2=df2/12/h

!公式3
df3=-f(-3)+6*f(-2)-18*f(-1)+10*f(0)+3*f(1)
df3=df3/12/h

!公式4
df4=3*f(-4)-16*f(-3)+36*f(-2)-48*f(-1)+25*f(0)
df4=df4/12/h




!新建文件用以保存数值微分结果
open(unit=file2,file='file2.txt')

write(file2,*)df1
write(file2,*)df2
write(file2,*)df3
write(file2,*)df4


write(file2,*)'实际结果应该为：'
write(file2,*)dcos(3d0)-sin(3d0)

end subroutine solve

end module fivediff
program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose     : 五点公式计算微分
!
!  Post Script :
!       1.
!       2.
!    
!-----------------------------------------------------
use fivediff

call solve

end
