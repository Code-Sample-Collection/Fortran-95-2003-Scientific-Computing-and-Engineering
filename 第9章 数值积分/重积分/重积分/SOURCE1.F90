module dbquad

!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-5-30
!-----------------------------------------------------
!  Description :   重积分之高斯勒让德方法模块
!    
!-----------------------------------------------------
!  Parameters  :
!      1.
!      2.
!----------------------------------------------------- 
!  Contains    :
!      1.   方法函数
!      2.   要计算的函数 
!-----------------------------------------------------

contains 

subroutine solve(func,s,a,b,c,d)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-5-30
!-----------------------------------------------------
!  Purpose   : 重积分之高斯勒让德方法
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.    func  外部函数
!       2.    a,b   积分区间
!  Output parameters  :
!       1.    s 积分结果
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.  选用5节点的公式
!       2.  计算时先把一般区间变换到[-1,1]区间上利用公式
!----------------------------------------------------

implicit real*8(a-z)
external func

integer::i,j
real*8::node(5),w(5),t1(5),t2(5)

node=(/-0.9061798459,-0.5384693101,0,0.5384693101,0.9061798459/)
w=(/0.2369268851,0.4786286705,0.5688888889,&
  0.4786286705,0.2369268851/)
  
  t1=(b+a)/2d0+(b-a)*node/2d0  
  t2=(c+d)/2d0+(d-c)*node/2d0  

s=0

do i=1,5
   do j=1,5
  
  call func(fx,t1(i),t2(j))
  
  s=s+fx*w(i)*w(j)
  end do
end do

s=s*(b-a)/2d0*(d-c)/2d0

end subroutine solve


subroutine fun1(f,x,y)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-5-30
!-----------------------------------------------------
!  Purpose   :  待计算的二元函数
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.   x ，y 自变量
!       2. 
!  Output parameters  :
!       1.  f  因变量
!       2.
!  Common parameters  :
!
!----------------------------------------------------
!  Post Script :
!       1.
!       2.
!----------------------------------------------------
implicit real*8(a-z)
f=dlog(x+2d0*y)
end subroutine fun1


end module dbquad


program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-5-30
!-----------------------------------------------------
!  Purpose   :  主函数
!    
!-----------------------------------------------------
!  In put data  files :
!       1.
!       2.
!  Output data files  :
!       1.  result.txt  计算结果存放文件
!       2.
!-----------------------------------------------------
!  Post Script :
!       1.
!       2.
!-----------------------------------------------------
use dbquad

implicit real*8(a-z)

open(unit=11,file='result.txt')

pi=3.141592653589793238D0

call solve(fun1,s,1.4d0,2d0,1d0,1.5d0)

write(11,101)s

101 format(/,T5,'重积分之高斯方法',//,&
            T3,'积分结果为：',F12.5  )

end program main