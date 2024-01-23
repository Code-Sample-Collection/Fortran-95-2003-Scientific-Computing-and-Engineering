module RF
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-4-5
!-----------------------------------------------------
!  Description : RF迭代法模块
!    
!-----------------------------------------------------
!  Parameters  :
!      1.    IMAX--最大允许迭代次数  
!      2.    tol--误差容限
!      3.    omiga --RF因子
!  Contains    :
!      1.    solve RF迭代法方法函数
!      2.
!-----------------------------------------------------
!  Post Script :
!      1.
!      2. 
!-----------------------------------------------------

implicit real*8(a-z)
integer::IMAX=200
real*8::tol=1d-7
real*8::omiga=1.5

contains

subroutine solve(A,b,x,x0,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  RF迭代法函数
!               用于计算方程 AX=b
!-----------------------------------------------------
!  Input  parameters  :
!       1.  A,b 意义即  AX=b
!       2.  x0迭代初值
!       3.  N 方程的维数
!  Output parameters  :
!       1. x 方程的解
!       2.
!  Common parameters  :
!
!----------------------------------------------------

implicit real*8(a-z)
integer::N
integer::i,j,k

real*8::A(N,N),b(N),x(N),x0(N)

real*8::x1(N),x2(N)



!写入标题
  write(102,501)
  501 format(//,18x,'RF迭代法',//)

x1=x0


do k=1,IMAX
   
   do i=1,N
        s=0
        
        do j=1,N
        
        s=s+A(i,j)*x1(j)       
        
        end do 
       !------------------------

  !RF迭代更新       
       x2(i)=x1(i)-omiga*(s-b(i))
   
   end do
   

 !这段程序用于判断精度，满足精度时退出循环   
   dx2=0
   do i=1,N
    dx2=dx2+(x1(i)-x2(i))**2
   end do
   dx2=dsqrt(dx2)
   
   if (dx2<tol)  exit
!----------------------------------------      
   x1=x2
  
  !记录迭代中间值
     write(102,502)k,x1
     502 format(I3,3F12.8)
  !----
end do 

x=x2

end subroutine solve

end module RF



program  main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-5
!-----------------------------------------------------
!  Purpose   :  采用RF迭代法计算线性方程
!    
!-----------------------------------------------------
!  In put data  files :
!       1.   
!       2.
!  Output data files  :
!       1. Im_result.txt计算的中间数据
!       2.  result.txt计算结果
!-----------------------------------------------------


use RF
implicit real*8(a-z)

integer,parameter::N=3

real*8 ::A(N,N),b(N),x(N),x0(N)

  open(unit=101,file='result.txt')
  open(unit=102,file='Im_result.txt')


  x0=(/0d0,0d0,0d0/)

  b=(/1,1,1/)

  A=reshape((/0.8,-0.01,0.12,&
             -0.01,0.84,0.05,&
             0.12,0.05,0.88/),(/3,3/))
            
  
  call solve(A,b,x,x0,N)
  
  write(101,501)x                                                                                                            
  501 format(/,T10,'RF迭代法',//,&
                2x,'x(1)=',F15.8,/,&
                2x,'x(2)=',F15.8,/,&
                2x,'x(3)=',F15.8,/)


end program main
  