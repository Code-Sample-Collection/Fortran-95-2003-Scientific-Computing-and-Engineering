module middiff
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
subroutine  solve()
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose     : 方法函数
! 
!  Post Script :
!       1.
!       2.
!       3.    
!-----------------------------------------------------

implicit real*8(a-z)

integer i

s=dsqrt(2d0)
h=1d0

!建立文件

open(unit=file1,file='result.txt')


do i=1,50

   f2=datan(s+h)
   f1=datan(s-h)
   
   d=f2-f1
   
   r=d/(h*2d0)
!文件写入的内容为序号，步长与计算结果，准确值为1/3
   write(file1,*)i,h,r
   h=h/2
end do


end subroutine solve

end module middiff

program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.08.04
!-----------------------------------------------------
!  Purpose     : 主函数
!
!  Post Script :
!       1.
!       2.
!    
!-----------------------------------------------------
use middiff
call solve

end
