module mat_eq
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-4-8
!-----------------------------------------------------
!  Description : 矩阵方程
!    
!-----------------------------------------------------
!  Contains    :
!      1.   driver  驱动函数
!      2.   solve   方法函数 
!-----------------------------------------------------

contains


subroutine solve(A,B,X,N,M)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-8
!-----------------------------------------------------
!  Purpose   :  高斯列主元消去法计算矩阵方程组
!                 AX=B
!-----------------------------------------------------
!  Input  parameters  :
!       1.   A(N,N)系数矩阵
!       2.   B(N,M)右向量
!       3.   N方程维数
!       4.   M右矩阵的列数
!  Output parameters  :
!       1.  X 方程的解矩阵
!       2.
!  Common parameters  :
!
!----------------------------------------------------
implicit real*8(a-z)

integer::i,k,N,M
integer::id_max  !主元素标号

real*8::A(N,N),B(N,M),X(N,M)

real*8::Aup(N,N),Bup(N,M)

!Ab为增广矩阵  [AB]
real*8::AB(N,N+M)

real*8::vtemp1(N+M),vtemp2(N+M)
real*8::vtmp(N),xtmp(N)

AB(1:N,1:N)=A

AB(:,N+1:N+M)=B


!##########################################################
!  这段是 列主元消去法的核心部分
do k=1,N-1

    elmax=dabs(Ab(k,k))
    id_max=k
    
    !这段为查找主元素	
    !这段程序的主要目的不是为了赋值最大元素给elmax，而是为了找出最大元素对应的标号

	
	do i=k+1,n
      if (dabs(Ab(i,k))>elmax) then
         elmax=Ab(i,k)

         id_max=i
      end if          
    end do

    
 !至此，已经完成查找最大元素，查找完成以后与  第k行交换 
 !交换两行元素，其他不变
    vtemp1=Ab(k,:)
    vtemp2=Ab(id_max,:)
   
    
    Ab(k,:)=vtemp2
    Ab(id_max,:)=vtemp1   
!
!以上一大段是为交换两行元素，交换完成以后即按照消元法进行
!#########################################################
  
   do i=k+1,N
  
     temp=Ab(i,k)/Ab(k,k)
     
     Ab(i,:)=Ab(i,:)-temp*Ab(k,:)
   
   end do

end do

!-----------------------------
! 经过上一步，Ab已经化为如下形式的矩阵
!            | *  *  *  *  #  #|
!     [A b]= | 0  *  *  *  #  #|
!            | 0  0  *  *  #  #|
!            | 0  0  0  *  #  #|
!
Aup(:,:)=AB(1:N,1:N)




do i=1,m
!调用用上三角方程组的回带方法
vtmp=AB(:,N+i)
call uptri(Aup,vtmp,xtmp,n)
!把计算结果赋值给X
X(:,i)=xtmp
end do

end subroutine solve



subroutine uptri(A,b,x,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-8
!-----------------------------------------------------
!  Purpose   :  上三角方程组的回带方法
!                 Ax=b
!-----------------------------------------------------
!  Input  parameters  :
!       1.   A(N,N)系数矩阵
!       2.   b(N)右向量
!       3.   N方程维数
!  Output parameters  :
!       1.  x  方程的根
!       2.
!  Common parameters  :
!
!----------------------------------------------------

implicit real*8(a-z)

integer::i,j,N

real*8::A(N,N),b(N),x(N)

x(N)=b(N)/A(N,N)

!回带部分
do i=n-1,1,-1
   
    x(i)=b(i)
   do j=i+1,N
    x(i)=x(i)-a(i,j)*x(j)
   end do
    x(i)=x(i)/A(i,i)

end do

end subroutine uptri

subroutine driver(N,M)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-8
!-----------------------------------------------------
!  Purpose   : 驱动程序
!-----------------------------------------------------
!  Input  parameters  :
!       1.   N     描述 A(N,N)
!       2.   M     描述方程   X(N,M)，B(N,M)
!  Output parameters  :
!       1.  
!       2.
!  P.S  :
!      N,M  从文件中读取
!
!----------------------------------------------------

implicit real*8(a-z)

integer::N,M

integer::i,j
real*8::A(N,N),B(N,M),X(N,M)

!读入系数A矩阵
read(11,*)((A(i,j),j=1,N),i=1,N)
!读入B矩阵
read(11,*)((B(i,j),j=1,M),i=1,N)


call solve(A,B,X,N,M)

write(12,101)

write(12,102)((X(i,j),j=1,m),i=1,n)

!变量输出格式只针对IVF编译器，在CVF中不支持
102 format(<N>(<m>F16.10/))


!do i=1,m
!    
!	do j=1,n
!     write(12,102)j,i,x(j,i)
!	end do
!    
!end do

101 format(T4,'消去法计算矩阵方程，版本1.1',/)
!102 format(T4,'x(',I2,',',I2,')=',f10.5)

write(12,103)

103 format(/'P.S:本软件既可以计算线性方程组Ax=b，也可以计算矩阵方程AX=B')
end subroutine driver


end module mat_eq




program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-8
!-----------------------------------------------------
!  Purpose   :  消去法解矩阵方程
!    
!-----------------------------------------------------
!  In put data  files :
!       1.  fin.txt  输入方程系数
!       2.
!  Output data files  :
!       1. fout.txt  计算结果
!       2.
!-----------------------------------------------------
!  Post Script :
!       1.    需要准备输入数据
!       2.    由于驱动函数调用方法函数
!-----------------------------------------------------
use mat_eq

integer::N,M

open(unit=11,file='fin.txt')
open(unit=12,file='fout.txt')

read(11,*)
!读说明文字

!读入方程维数系数
read(11,*)N,M

!调用驱动函数
call driver(N,M)


end program main