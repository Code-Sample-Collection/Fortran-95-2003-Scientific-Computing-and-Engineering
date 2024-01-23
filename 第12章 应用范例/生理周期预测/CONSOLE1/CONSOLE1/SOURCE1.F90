module predict
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
!  Description :   生理周期预测模块
!                   
!  Post Script :
!      1.      
!      2. 
!
!-----------------------------------------------------
!  Contains    :
!      1.   方法函数
!      2.   儒略日化公历
!      3.   公历化儒略日
!-----------------------------------------------------
!  Parameters  :
!      1.
!      2.
!----------------------------------------------------- 

contains

subroutine solve(bith,star_index,phy_cli,emo_cli,inte_cli)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose     :   生理(Physiological)、情绪(Emotional)
!                  智力(Intelligence)预测模块
!                  
!  Post Script :
!       1.      指数分为0-10分 
!       2.
!       3.    
!-----------------------------------------------------
!  Input  parameters  :
!       1.          bith（3）--输入生日 整型 年月日
!       2.    
!  Output parameters  :
!       1.    star_index（3）双精度  当前 生理、情绪、智力预测指数，
!                               最高10分，最低0分
!       2.   phy_cli(3) 整型 下一个生理高峰  日期  
!       3.   emo_cli(3)  整型 下一个情绪高峰日期  整型      
!       4.  inte_cli(3) 整型  下一个智力高峰日期 
!
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------

implicit real*8(a-z)

integer::bith(3),now(3),next_climax(3)

integer::phy_cli(3),emo_cli(3),inte_cli(3)


!span 为生日到当前日期的 总天数
integer::y,m,d,span1
real*8::star(3),span,jdnow,temp(3)

real*8::star_index(3)

!调用计算机函数，获取当前日期
call GETDAT (Y,M,D)

now(1)=y
now(2)=m
now(3)=d

!把生日及当前时间转化成儒略日
call gre2jule(bith,jdbith)
call gre2jule(now,jdnow)


twopi=6.2831853071795864

span=jdnow-jdbith

!让指标为0分到10分
star_index(1)=5d0*dsin(twopi*span/23d0)+5d0

star_index(2)=5d0*dsin(twopi*span/28d0)+5d0

star_index(3)=5*dsin(twopi*span/33d0)+5d0

!以上已经算出当日指数，下面计算下一个高峰日期





span1=int(span)

!直接算出下次高峰的儒略日
next_climax(1)=23-mod(span1+23*3/4,23)+int(jdnow)+1
next_climax(2)=28-mod(span1+28*3/4,28)+int(jdnow)+1
next_climax(3)=33-mod(span1+33*3/4,33)+int(jdnow)+1

!转为双精度
temp=DBLE(next_climax)


!把下一个高峰的儒略日化为公历
call  jule2gre(temp(1),phy_cli)
call  jule2gre(temp(2),emo_cli)
call  jule2gre(temp(3),inte_cli)

end subroutine solve

subroutine gre2jule(date,JD)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.07.07  
!-----------------------------------------------------
!  Purpose     : 公历（格里高利）化儒略日
! 
!  Post Script :
!       1.    
!       2.
!       3.    
!-----------------------------------------------------
!  Input  parameters  :
!       1.   date(3)  整型年月日
!       2.   
!  Output parameters  :
!       1.  JD  儒略日
!       2.
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------
integer date(3),Y,M,D  
real*8 j,jd
real*8 temp1,temp2,temp3

Y=date(1)
m=date(2)
d=date(3)

temp1=1461*(Y+4800+(M-14)/12)/4
temp2=367*(M-2-((m-14)/12)*12)/12
temp3=3*(Y+4900+(M-14)/12)/100/4
     
JD=D-32075+temp1+temp2-temp3-0.5d0
end subroutine gre2jule


subroutine jule2gre(JD,date)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010.07.07
!-----------------------------------------------------
!  Purpose     :  儒略日化成公历（格里高利）
! 
!  Post Script :
!       1.JD允许带小数  
!       2.
!       3. 编程方法可以参看《航天器轨道确定》，于志坚编
!-----------------------------------------------------
!  Input  parameters  :
!       1.  JD   浮点数，儒略日
!       2.  date(3)  整型 分别为年月日
!  Output parameters  :
!       1.
!       2.
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------

implicit real*8(a-z)
integer::date(3)

integer::J,N,L1,Y1,L2,M1,D,L3,M,Y

J=int(JD+0.5D0)
N=4*(J+68569)/146097
L1=J+68569-(N*146097+3)/4
Y1=4000*(L1+1)/1461001
L2=L1-1461*Y1/4+31

M1=80*L2/2447
D=L2-2447*M1/80
L3=INT(M1/11)
M=M1+2+12*L3
Y=100*(N-49)+Y1+L3

date(1)=Y
date(2)=M
date(3)=D

end subroutine jule2gre


end module predict



program main
!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose     :  生理、情绪、智力周期预测函数
!
!  Post Script :   
!       1.     今天的时间，从计算机读取
!       2.
!    
!-----------------------------------------------------
!  In put data  files :
!       1.     bithday.txt 生日 ----- 年月日
!       2.
!  Output data files  :
!       1.     predict.txt   当天生理、情绪、智力指数
!       2.      以及下次高峰的日期
!
!-----------------------------------------------------
use predict
implicit real*8(a-z)

integer::now(3),next_climax(3),bith(3) 

integer::a(3),b(3),c(3)

real*8::star_index(3)
integer::y,m,d

open(unit=11,file='birthday.txt')
open(unit=12,file='predict.txt')


read(11,*)
!读入生日
read(11,*)bith


!调用计算机函数，获取当前日期
call GETDAT (Y,M,D)

now(1)=y
now(2)=m
now(3)=d


call solve(bith,star_index,a,b,c)

write(12,101)now,star_index,a,b,c

101 format(/,T12,'生理周期预测',//,&
             T3,'今天是公元:',3(I6),//,&
             T3,'您今日生理指数为：',F6.1,/,&
             T3,'您今日情绪指数为：',F6.1,/,&
             T3,'您今日智力指数为：',F6.1,//,&
             T3,'您下次生理峰值日期为：',3(I6),/,&
             T3,'您下次情绪峰值日期为：',3(I6),/,&
             T3,'您下次智力峰值日期为：',3(I6),/ )

end program main