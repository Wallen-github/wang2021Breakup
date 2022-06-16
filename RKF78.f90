!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! @file: RkfModule.f90
! @brief: 此module里涵盖多种常用的积分器
! @author: Wang Hai-Shuo
! @time: 2020.3.3
    
    
! @sub RKF78
! @pram: YHC 右函数 
! @pram: H 初始步长
! @pram: T 初始积分时间
! @pram: x0 状态矢量   
! @pram: EE 误差限
! @pram: DD 变步长指示器，DD>5D-1时为固定步长积分
! @pram: N 状态矢量参数个数
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
SUBROUTINE RKF78(YHC,H,T,X,EE,DD,N)
	IMPLICIT REAL*8(A-H,O-Z)
	DOUBLE PRECISION :: X(N),X1(N),Y0(N),Y1(N),Y2(N),Y3(N),Y4(N),Y5(N),Y6(N),Y7(N),Y8(N),Y9(N),Y10(N),Y11(N),Y12(N)
    external :: YHC
    
	DO 16 I=1,N
16	X1(I)=X(I)
	T1=T
	CALL YHC(T,T1,X1,Y0)
	DO 1 I=1,N
1	X1(I)=X(I)+H*2D0/27D0*Y0(I)
	T1=T+H*2D0/27D0
	CALL YHC(T,T1,X1,Y1)
	DO 2 I=1,N
2	X1(I)=X(I)+H*(Y0(I)+3D0*Y1(I))/36D0
	T1=T+H*1D0/9D0
	CALL YHC(T,T1,X1,Y2)
	DO 3 I=1,N
3	X1(I)=X(I)+H*(Y0(I)+3D0*Y2(I))/24D0
	T1=T+H*1D0/6D0
	CALL YHC(T,T1,X1,Y3)
	DO 4 I=1,N
4	X1(I)=X(I)+H*(Y0(I)*20D0+(-Y2(I)+Y3(I))*75D0)/48D0
	T1=T+H*5D0/12D0
	CALL YHC(T,T1,X1,Y4)
	DO 5 I=1,N
5	X1(I)=X(I)+H*(Y0(I)+Y3(I)*5D0+Y4(I)*4D0)/20D0
	T1=T+H*1D0/2D0
	CALL YHC(T,T1,X1,Y5)
	DO 6 I=1,N
6	X1(I)=X(I)+H*(-Y0(I)*25D0+Y3(I)*125D0-Y4(I)*260D0+Y5(I)*250D0)/108D0
	T1=T+H*5D0/6D0
	CALL YHC(T,T1,X1,Y6)
	DO 7 I=1,N
7	X1(I)=X(I)+H*(Y0(I)*93D0+Y4(I)*244D0-Y5(I)*200D0+Y6(I)*13D0)/900D0
	T1=T+H*1D0/6D0
	CALL YHC(T,T1,X1,Y7)
	DO 8 I=1,N
8	X1(I)=X(I)+H*(Y0(I)*180D0-Y3(I)*795D0+Y4(I)*1408D0-Y5(I)*1070D0+Y6(I)*67D0+Y7(I)*270D0)/90D0
	T1=T+H*2D0/3D0
	CALL YHC(T,T1,X1,Y8)
	DO 9 I=1,N
9	X1(I)=X(I)+H*(-Y0(I)*455D0+Y3(I)*115D0-Y4(I)*3904D0+Y5(I)*3110D0-Y6(I)*171D0+Y7(I)*1530D0-Y8(I)*45D0)/540D0
	T1=T+H*1D0/3D0
	CALL YHC(T,T1,X1,Y9)
	DO 10 I=1,N
10	X1(I)=X(I)+H*(Y0(I)*2383D0-Y3(I)*8525D0+Y4(I)*17984D0-Y5(I)*15050D0+Y6(I)*2133D0+Y7(I)*2250D0+Y8(I)*1125D0+Y9(I)*1800D0)/4100D0
	T1=T+H
	CALL YHC(T,T1,X1,Y10)
	DO 11 I=1,N
11	X1(I)=X(I)+H*(Y0(I)*60D0-Y5(I)*600D0-Y6(I)*60D0+(Y8(I)-Y7(I)+2D0*Y9(I))*300D0)/4100D0
	T1=T
	CALL YHC(T,T1,X1,Y11)
	DO 12 I=1,N
12	X1(I)=X(I)+H*(-Y0(I)*1777D0-Y3(I)*8525D0+Y4(I)*17984D0-Y5(I)*14450D0+Y6(I)*2193D0+Y7(I)*2550D0+Y8(I)*825D0+Y9(I)*1200D0+Y11(I)*4100D0)/4100D0
	T1=T+H
	CALL YHC(T,T1,X1,Y12)
	DO 13 I=1,N
13  X(I)=X(I)+H*(Y5(I)*272D0+(Y6(I)+Y7(I))*216D0+(Y8(I)+Y9(I))*27D0+(Y11(I)+Y12(I))*41D0)/840D0
     
	EE=0D0
	DO 14 I=1,N
14	EE=EE+DABS((Y0(I)+Y10(I)-Y11(I)-Y12(I))*H*41D0/840D0)

	T=T+H
	IF(DD.GT.0.5D0) GO TO 15
	IF(EE.GT.DD) H=H/2D0
	IF(EE.LT.DD*1D-2) H=H*2D0
15  RETURN
    END