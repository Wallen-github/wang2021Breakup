!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! @file: ES_Method.f90
! @brief: 此module里包括多个计算周期轨道的必用方法
! @author: Wang Hai-Shuo
! @time: 2020.6.2
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    SUBROUTINE Monodromy(YHC,X,tf,M0,MM)
    EXTERNAL :: YHC
    REAL(8) :: M0(72),X(8),t0,tf,hh,err,EE,MM(8,8),Meye(8,8)
109 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18)
    N = 8
    ! 初始化Monodromy矩阵，单位矩阵
    Meye(:, :) = 0D0
    DO I = 1, N
        Meye(I, I) = 1D0
    ENDDO
    
    DO I = 1,N
        DO J = 1,N
            IF(I == J)THEN
                M0((I-1)*8+J) = 1.0D0
            ELSE
                M0((I-1)*8+J) = 0.0D0
            END IF
        END DO
    END DO
    
    M0(N*N+1:N*N+N) = X
    
    !积分
    
    err = 1D-12
    hh = 1D-5
    t0 = 0D0
    DO WHILE (t0<tf)
        Call RKF78(YHC,hh,t0,M0,EE,err,N*N+N)
        !write(14,109) M0(N*N+1:N*N+N),t0*Tunit
    ENDDO
    Call RKF78(YHC,tf-t0,t0,M0,EE,1D0,N*N+N) ! 积分积满，恰好积分到tf
    

    DO I=1,N*N
	    J=INT((I+N-1)/N)
	    MM(J,I-N*(J-1)) = M0(I)
	ENDDO
    MM = MM !- Meye
    END
    
    SUBROUTINE PeriodOrbitT(YHC,YHCM,TOF,X,T,eig,flag)
    use lapack95
    use global
    EXTERNAL :: YHC,YHCM
    integer, parameter :: N = 8
    integer, parameter :: NCoe = 4
    REAL(8) :: X(N),XT(N),dX(N),XB(N),dXB(N),MM(N,N),MMB(N,N),XX(N*N+N),T,TOF
    REAL(8) :: dF(NCoe),Coe(NCoe,NCoe),Anti_Coe(NCoe,NCoe),Sol(NCoe),JS(NCoe)
    REAL(8) :: eig
    real :: D
    INTEGER(4) :: flag,loopcount
102 format(1X, E25.18,1X,E25.18)    
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
104 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
105 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
106 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
107 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18)
108 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18)
109 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18)
    loopcount = 1
    TOF0 = 1D0
    flag = 1
    !open(14,file='G:\Master\科研工作\Binary system\SynBinary\orbit4New1.txt', status='REPLACE')
    DO WHILE (TOF0 > TOF)
        
        ! 将状态量拷贝进积分变量
        XB = X
        
        ! 计算Monodromy矩阵
        Call Monodromy(YHCM,X,T/2D0,XX,MM)
        
        ! 求解deltaX
        XT = XX(65:72)
        !CALL YHC(t0,t1,XT,dX)
        
        ! 赋值系数矩阵
        Coe(1,:) = (/MM(2,1),MM(2,6),MM(2,7),MM(2,8)/)
        Coe(2,:) = (/MM(3,1),MM(3,6),MM(3,7),MM(3,8)/)
        Coe(4,:) = (/MM(5,1),MM(5,6),MM(5,7),MM(5,8)/)
        dF(1) = X(2) - XT(2)
        dF(2) = X(3) - XT(3) + SIGN(1D0,XT(7))*PI
        dF(4) = 0D0 - XT(5)
        
        ! 计算Monodromy矩阵
        Call Monodromy(YHCM,X,2D0*PI*S0**(3D0/2D0),XX,MMB)
        
        ! 求解thetaB
        XB = XX(65:72)
        !CALL YHC(t0,t1,XB,dXB)
        
        ! 赋值系数矩阵
        Coe(3,:) = (/MMB(4,1),MMB(4,6),MMB(4,7),MMB(4,8)/)
        dF(3) = 0D0 - XB(4) + 2D0*PI
            
        ! 求解修正量
        CALL ANTIMATRIX(Coe,Anti_Coe,NCoe,flag)
        Sol = MATMUL(Anti_Coe,dF)
        
        ! 返回错误信息
        IF(flag == 0)THEN
            write(*,*) 'err: The Matrix can not inverse. locate: sub PeriodOrbitT'
            exit
        ENDIF
        loopcount  = loopcount + 1
        IF(loopcount > 50)THEN
            flag = 0
            write(*,*) 'err: Loop Count greater than ',loopcount,'. locate: sub PeriodOrbitT'
            exit
        ENDIF
        
        ! 修正初始状态量
        !write(*,104) X(1),X(6),X(7),X(8)
        X(1) = X(1) + Sol(1)
        X(6) = X(6) + Sol(2)
        X(7) = X(7) + Sol(3)
        X(8) = X(8) + Sol(4)
        TOF0 = DSQRT(Sol(1)**2+Sol(2)**2+Sol(3)**2+Sol(4)**2)
        !write(*,*) 'TOF0 = ',TOF0
        ! 返回错误信息
        IF(TOF0>1D1)THEN
            flag = 0
            write(*,*) 'err: the correction is not convergent. locate: sub PeriodOrbitT'
            exit
        ENDIF
        
    ENDDO
    eig = 0D0
    do I = 1,N
        eig = eig + MM(I,I)
    ENDDO
    eig = eig -2
    !close(14)
    end
    
    SUBROUTINE PeriodOrbitS(YHC,YHCM,TOF,X,T,eig,flag)
    use global
    EXTERNAL :: YHC,YHCM
    INTEGER(4),parameter :: N = 8, NCoe = 4
    REAL(8) :: X(N),XT(N),dX(N),XB(N),dXB(N),MM(N,N),MMB(N,N),XX(N*N+N),T,TOF,eig
    REAL(8) :: dF(NCoe),Coe(NCoe,NCoe),Anti_Coe(NCoe,NCoe),Sol(NCoe),JS(NCoe)

    INTEGER(4) :: flag,loopcount
102 format(1X, E25.18,1X,E25.18)    
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
104 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
105 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
106 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
107 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18)
108 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18)
109 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18)
    loopcount = 1
    TOF0 = 1D0
    flag = 1
    !open(14,file='G:\Master\科研工作\Binary system\SynBinary\orbit4New1.txt', status='REPLACE')
    DO WHILE (TOF0 > TOF)
        
        ! 计算Monodromy矩阵
        Call Monodromy(YHCM,X,T/2D0,XX,MM)
        
        ! 求解deltaX
        XT = XX(65:72)
        CALL YHC(T/2D0,t1,XT,dX)
        
        ! 赋值系数矩阵
        Coe(1,:) = (/MM(2,6),MM(2,7),MM(2,8),-dX(2)/)
        Coe(2,:) = (/MM(3,6),MM(3,7),MM(3,8),-dX(3)/)
        Coe(4,:) = (/MM(5,6),MM(5,7),MM(5,8),-dX(5)/)
        dF(1) = X(2) - XT(2)
        dF(2) = X(3) - XT(3) + SIGN(1D0,XT(7))*PI
        dF(4) = 0D0 - XT(5)
        
        ! 计算Monodromy矩阵
        Call Monodromy(YHCM,X,2D0*PI*S0**(3D0/2D0),XX,MMB)
        
        ! 求解thetaB
        XB = XX(65:72)
        CALL YHC(2D0*PI*S0**(3D0/2D0),t1,XB,dXB)
        
        ! 赋值系数矩阵
        Coe(3,:) = (/MMB(4,6),MMB(4,7),MMB(4,8),-dXB(4)/)
        dF(3) = 0D0 - XB(4) + 2D0*PI
            
        ! 求解修正量
        CALL ANTIMATRIX(Coe,Anti_Coe,NCoe,flag)
        Sol = MATMUL(Anti_Coe,dF)
        
        ! 返回错误信息
        IF(flag == 0)THEN
            write(*,*) 'err: The Matrix can not inverse. locate: sub PeriodOrbitS'
            exit
        ENDIF
        loopcount  = loopcount + 1
        IF(loopcount > 50)THEN
            flag = 0
            write(*,*) 'err: Loop Count greater than ',loopcount,'. locate: sub PeriodOrbitT'
            exit
        ENDIF
        
        ! 修正初始状态量
        !write(*,104) X(6),X(7),X(8),T
        X(6) = X(6) + Sol(1)
        X(7) = X(7) + Sol(2)
        X(8) = X(8) + Sol(3)
        T = T + Sol(4)
        TOF0 = DSQRT(Sol(1)**2+Sol(2)**2+Sol(3)**2+Sol(4)**2)
        !write(*,*) 'TOF0 = ',TOF0
        ! 返回错误信息
        IF(TOF0>1D1)THEN
            flag = 0
            write(*,*) 'err: the correction is not convergent. locate: sub PeriodOrbitS'
            exit
        ENDIF
        
    ENDDO
    eig = 0D0
    do I = 1,N
        eig = eig + MM(I,I)
    ENDDO
    eig = eig -2
    !close(14)
    end
    
    SUBROUTINE SOLVEEQUATIONSET(A,B,X,N,L,JS)
    ! @sub SOLVEEQUATIONSET 全选主元高斯消去法求解线性代数方程组 AX = B
    ! @pram: A  双精度实型二维数组，N*N，输入参数。存放方程组系数矩阵，返回时将被破坏
    ! @pram: B  双精度实型二维数组，N，输入参数。存放方程组右端向量，返回时将被破坏
    ! @pram: N  整型变量，输入参数。存放方程组的阶数。   
    ! @pram: X  双精度实型一维数组，长度为N，输出参数。返回方程组的解向量。
    ! @pram: L  整型变量，输出参数。若返回L=0，则说明方程组的系数矩阵奇异，求解失败，其余值说明返回正常
    ! @pram: JS 整型一维数组，N,本子程序的工作数组。
    DIMENSION :: A(N,N),X(N),B(N),JS(N)
    REAL(8) :: A,B,X,T,JS
    L = 1
    DO K = 1,N-1
        D = 0D0
        DO I = K,N
            DO J = K,N
                IF(ABS(A(I,J)) .GT. D)THEN
                    D = ABS(A(I,J))
                    JS(K) = J
                    IS = I
                ENDIF
            ENDDO
        ENDDO
        IF(D+1D0 == 1D0)THEN
            L = 0
        ELSE
            IF(JS(K) /= K)THEN
                DO I = 1,N
                    t = A(i,K)
                    A(I,K) = A(I,JS(K))
                    A(I,JS(K)) = t
                ENDDO
            ENDIF
            IF(IS /= K) THEN
                DO J = K,N
                    t = A(K,J)
                    A(K,J) = A(IS,J)
                    A(IS,J) = t
                ENDDO
                t = B(K)
                B(K) = B(IS)
                B(IS) = t
            ENDIF
        ENDIF
        IF(L == 0D0) THEN
            write(*,100)
            return
        ENDIF
        DO J = K+1,N
            A(K,J) = A(K,J)/A(K,K)
        ENDDO
        B(K) = B(K)/A(K,K)
        DO I = K+1,N
            DO J = K+1,N
                A(I,J) = A(I,J) - A(I,K)*A(K,J)
            ENDDO
            B(I) = B(I) - A(I,K)*B(K)
        ENDDO
    ENDDO
    IF(ABS(A(N,N))+1D0 == 1D0) THEN
        L = 0
        write(*,100)
        RETURN
    ENDIF
    X(N) = B(N)/A(N,N)
    DO I = N-1,1,-1
        t = 0D0
        DO J = I+1,N
            t = t + A(I,J)*X(j)
        ENDDO
        X(I) = B(I) - t
    ENDDO
100 FORMAT(1x,'err: Solve Equation Set Fail. locate: sub SOLVEEQUATIONSET')
    JS(N) = N
    DO K = N,1,-1
        IF(JS(K) /= K) THEN
            t = X(K)
            X(K) = X(JS(K))
            X(JS(K)) = t
        ENDIF
    ENDDO
    RETURN
    
    END
    
!-----------递归求行列式的值--------------------

    recursive function det(A,col,row) result(D)

    Implicit None

    integer row,col

    real(8) :: A(col,row),B(row-1,col-1)

    real :: D

    integer row_now,col_now,k,c,f

    row_now=row 

    col_now=col

    if (row_now>1) then

        D = 0.0;

        do k=1,row_now

            if(k==1)then

                B=A(2:col_now,2:row_now)

            elseif(k<row_now)then

                B(1:k-1,:)=A(1:k-1,2:row_now)

                B(k:col_now-1,:)=A(k+1:col_now,2:row_now)

            else

                B=A(1:col_now-1,2:row_now)

            endif

            c=col-1

            f=row-1

            D = D + (-1)**(1+k) * A(k,1) *&

                det(B,c,f)

        end do

    else

        D = A(1,1);

    end if

    end function det    
    
    SUBROUTINE FindInitialVlaue(X0,T0,max,min,step,TOF)
    use global
    EXTERNAL :: YHC4New,YHC4_Mon
    REAL(8) :: rAMid,TOF,eigvalue,max,min,step,epsilon
    REAL(8) :: T0,X0(8),TMid,XMid(8)
    INTEGER :: flag
    
    epsilon = min
    
    DO WHILE (ABS(epsilon-max)<TOF)
        
        rAMid = rA
        Tmid = T0
        Xmid = X0
        
201     epsilon = epsilon + step
        cB = 480D0!rB / (a_b * b_c ** 2D0) ** (1D0 / 3D0)
        bB = 480D0!cB * b_c
        aB = 5D2!bB * a_b
        cA = 96D1!rA / (a_bA * b_cA ** 2D0) ** (1D0 / 3D0)
        bA = 96D1!cA * b_cA
        aA = 1D3!bA * a_bA
        Length = R!(aA + aB)
        S0 = R / Length
        MA = 4D0 / 3D0 * PI * rA ** 3D0 * p
        MB = 4D0 / 3D0 * PI * rB ** 3D0 * p
        mu = MB / (MA + MB)
        m = mu * (1 - mu)
        alphaB = rB / Length
        alphaA = rA / length
        J2 = (aB ** 2D0 + bB ** 2D0 - 2D0 * cB ** 2D0) / (10D0 * rB ** 2D0)
        J2A = epsilon*(aA ** 2D0 + bA ** 2D0 - 2D0 * cA ** 2D0) / (10D0 * rA ** 2D0)
        J22 = (aB ** 2D0 - bB ** 2D0) / (20D0 * rB ** 2D0)
        J22A = epsilon*(aA ** 2D0 - bA ** 2D0) / (20D0 * rA ** 2D0)
        A1 = alphaB ** 2D0 * J2 + alphaA ** 2D0 * J2A
        A2 = 6D0 * alphaB ** 2D0 * J22
        A3 = 6D0 * alphaA ** 2D0 * J22A
        IzB = mu * (aB ** 2D0 + bB ** 2D0) / (5D0 * Length ** 2D0)
        IzA = (1-mu) * (aA ** 2D0 + bA ** 2D0) / (5D0 * Length ** 2D0)
        K = dsqrt((S0 ** 2D0 + 3D0 / 2D0 * (A1 + A2)) * (IzB + m * S0 ** 2D0) ** 2D0 / S0 ** 5D0)
        Tunit = dsqrt(Length ** 3D0 / (MA + MB) / G) / SecOfYear
        
        ! 修正周期轨道初值
        CALL PeriodOrbitT(YHC4New,YHC4_Mon,TOF,X0,T0, eigvalue, flag)
        IF (flag == 0) THEN
            write(*,*) 'warning: the step is too large. locate: sub FindInitialVlaue.'
            step = step /2D0
            rA = rAMid
            X0 = Xmid
            T0 = Tmid
            GO TO 201
        END IF
        write(*,*) 'epsilon = ',epsilon
        write(*,*) '======================================================='
        
    ENDDO
    
    end
    
    
    
    
    
    