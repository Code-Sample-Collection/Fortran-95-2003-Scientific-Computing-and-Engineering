module inv_mat
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-4-8
!-----------------------------------------------------
!  Description : 计算逆矩阵
!    
!-----------------------------------------------------
!  Contains    :
!      1.   driver  驱动函数
!      2.   solve   方法函数 
!-----------------------------------------------------

contains

subroutine drive(N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  驱动程序
!               用以读取数据，计算结果并把结果保存于文件中
!
!----------------------------------------------------
!  Post Script :
!       1.     fin   输入文件
!       2.     fout  输出文件
!----------------------------------------------------
implicit real*8(a-z)

integer::N,i,j
real*8::A(N,N),invA(N,N)

read(11,*)((A(i,j),j=1,N),i=1,N)

write(12,101)
101 format(/,T12,'逆矩阵为',/)

call solve(A,invA,N)


write(12,102)((invA(i,j),j=1,n),i=1,n)

!变量输出格式只针对 IVF编译器，在CVF中不支持
102 format(<N>(<N>F16.10/))

end subroutine drive




subroutine solve(A,invA,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  计算逆矩阵
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.
!       2. 
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

integer::n
integer::i
real*8::A(n,n),invA(n,n),E(n,n)

E=0

!设置E为单位矩阵
do i=1,n
   E(i,i)=1
end do

call mateq(A,E,invA,N,N)


end subroutine solve


subroutine mateq(A,B,X,N,M)
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
!       2.   B(N,M)右矩阵
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

!Ab为增广矩阵 [AB]
real*8::AB(N,N+M)

real*8::vtemp1(N+M),vtemp2(N+M)
real*8::vtmp(N),xtmp(N)

AB(1:N,1:N)=A

AB(:,N+1:N+M)=B


!##########################################################
!  这段是列主元消去法的核心部分
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

    
 !至此，已经完成查找最大元素，查找完成以后与 第k行交换
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

end subroutine mateq


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

end module inv_mat




program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-8
!-----------------------------------------------------
!  Purpose   :  逆矩阵的计算
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
use inv_mat

integer::N

!N表示矩阵的维数，由文件读入

open(unit=11,file='fin.txt')
open(unit=12,file='fout.txt')

read(11,*)
!读说明文字

!读入方程维数系数
read(11,*)N


call drive(N)


end program main