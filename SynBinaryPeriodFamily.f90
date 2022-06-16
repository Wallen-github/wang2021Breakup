    PROGRAM main
    USE global
    IMPLICIT REAL(8) (A-H,O-Z)
    EXTERNAL :: YHC4New,YHC4_Mon
    REAL(8) :: XX(8),TT,TOF,step,Angular,max,min,Xmid(8),Tmid,thetaA1mid,eig
    REAL(8) :: EE,hh,err,X0(8),t0
    INTEGER(4),parameter :: Num = 200
    INTEGER(4) :: flag,flagindex(Num),sepflag = 10
    integer(4) :: time_begin, time_end
    
102 format(1X, E25.18,1X,E25.18)    
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
104 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
105 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
106 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
107 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18)
108 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18)
109 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18)
110 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
111 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
112 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
    
    call system_clock(time_begin)
    cB = 480D0!rB / (a_b * b_c ** 2D0) ** (1D0 / 3D0)
    bB = 480D0!cB * b_c
    aB = 5D2!bB * a_b
    cA = 96D1!rA / (a_bA * b_cA ** 2D0) ** (1D0 / 3D0)
    bA = 96D1!cA * b_cA
    aA = 1D3!bA * a_bA
    Length = R !(aA + aB)
    MA = 4D0 / 3D0 * PI * rA ** 3D0 * p
    MB = 4D0 / 3D0 * PI * rB ** 3D0 * p
    mu = MB / (MA + MB)
    m = mu * (1 - mu)
    alphaB = rB / Length
    alphaA = rA / length
    J2B = (aB ** 2D0 + bB ** 2D0 - 2D0 * cB ** 2D0) / (10D0 * rB ** 2D0)
    J2A = (aA ** 2D0 + bA ** 2D0 - 2D0 * cA ** 2D0) / (10D0 * rA ** 2D0)
    J22B = (aB ** 2D0 - bB ** 2D0) / (20D0 * rB ** 2D0)
    J22A = (aA ** 2D0 - bA ** 2D0) / (20D0 * rA ** 2D0)
    A1 = (alphaB ** 2D0 * J2B + alphaA ** 2D0 * J2A)/2D0
    A2 = 3D0 * alphaB ** 2D0 * J22B
    A3 = 3D0 * alphaA ** 2D0 * J22A
    IzB = mu * (aB ** 2D0 + bB ** 2D0) / (5D0 * Length ** 2D0)
    IzA = (1 - mu) * (aA ** 2D0 + bA ** 2D0) / (5D0 * Length ** 2D0)
    Izt = IzB+m*S0**2D0
    Tunit = dsqrt(Length ** 3D0 / (MA + MB) / G) / SecOfYear
    
    thetaA1_Const = 0.79D0
    S0 = R / Length
    theta = 0D0
    K = IzA*thetaA1_Const + dsqrt((S0 ** 2D0 + 3D0 * (A1 + A2 * dcos(2D0*theta))) * (IzB + m * S0 ** 2D0) ** 2D0 / S0 ** 5D0)
    
    write(*,*) '==========================Initial Value==========================='
    write(*,*) 'Tunit = ', Tunit, 'year'
    write(*,*) 'Tunit = ', Tunit*365D0, 'day'
    write(*,*) 'Tunit = ', Tunit*365D0*24D0, 'hour'
    write(*,*) 'Tunit = ', Tunit*365D0*24D0*6D1, 'min'
    
    ! 设置循环参数
    flag = 1
    TOF = 1D-10
    max = thetaA1_Const
    min = 0D0!(S0) ** (-3D0/2D0)
    step = (max-min)/Num
    
    ! 赋初值
    thetaA1_Const = max
    XX = (/S0, 0D0, 0D0, 0D0, 0D0, 0D0, thetaA1_Const-(S0) ** (-3D0/2D0), (S0) ** (-3D0/2D0)/)
    TT = 2D0*PI/ABS(XX(7))
    
    open(13,file='G:\Master\科研工作\Binary system\2020.7.9Paperwork\Pic\DynamicSubstitute\PeriodOrbitInfo2.txt', status='REPLACE')
    open(15,file='G:\Master\科研工作\Binary system\2020.7.9Paperwork\Pic\DynamicSubstitute\PeriodFamily2.txt', status='REPLACE')
    DO I = 1,Num
        
        
        ! 输出相关信息
        write(*,*) '===================================================================='
        write(*,*) 'Num,T,flag =',I,TT,flag
        write(*,*) 'XT = '
        write(*,104) XX
        write(*,*) '===================================================================='
        
        
        IF (flag == 1) Then
            
            Xmid = XX
            Tmid = TT
            thetaA1mid = thetaA1_Const
            
            ! 修正周期轨道初值
            CALL PeriodOrbitT(YHC4New,YHC4_Mon,TOF,XX,TT,eig,flag)
        
            ! 返回错误信息
            IF(flag == 0)THEN
                flag = 2
                XX = Xmid
                TT = Tmid
                thetaA1_Const = thetaA1mid
                write(*,*) 'err: Find Period Orbit with T Failed. locate: Main'
                GO TO 203
            ENDIF
        
            ! 计算角动量
            Angular = m * XX(1) ** 2D0 * (XX(6) + XX(8)) + IzA * thetaA1_Const + IzB * XX(8)
        
            ! 写入文件
            write(13,112) XX,TT,Angular, thetaA1_Const,eig
            
            ! 积分轨道
            X0 = XX
            t0 = 0D0
            hh = 1D-5
            err = 1D-12
            DO WHILE (t0<TT)
                Call RKF78(YHC4New,hh,t0,X0,EE,err,8)
                if (MOD(I,sepflag)==0) write(15,109) X0,t0
            ENDDO
            Call RKF78(YHC4New,TT-t0,t0,X0,EE,1D0,8) ! 积分积满，恰好积分到tf
            if (MOD(I,sepflag)==0) write(15,*) ''
            
            ! 更新状态量
            thetaA1_Const = max - step * I
            TT = 2D0*PI/DABS(thetaA1_Const-XX(8))
            XX = (/XX(1), 0D0, 0D0, 0D0, 0D0, XX(6), thetaA1_Const-XX(8), XX(8)/)
        
        ENDIF
        
        IF (flag == 2) Then
            
            Xmid = XX
            Tmid = TT
            thetaA1mid = thetaA1_Const
            
            ! 修正周期轨道初值
            CALL PeriodOrbitS(YHC4New,YHC4_Mon,TOF,XX,TT,eig,flag)
        
            ! 返回错误信息
            IF(flag == 0)THEN
                flag = 1
                XX = Xmid
                TT = Tmid
                thetaA1_Const = thetaA1mid
                write(*,*) 'err: Find Period Orbit with S Failed. locate: Main'
                GO TO 203
                
            ENDIF
        
            ! 计算角动量
            Angular = m * XX(1) ** 2D0 * (XX(6) + XX(8)) + IzA * thetaA1_Const + IzB * XX(8)
        
            ! 写入文件
            write(13,112) XX,TT,Angular,2D0*PI/TT+XX(8),eig
            
            ! 积分轨道
            X0 = XX
            t0 = 0D0
            hh = 1D-5
            err = 1D-12
            DO WHILE (t0<TT)
                Call RKF78(YHC4New,hh,t0,X0,EE,err,8)
                if (MOD(I,sepflag)==0) write(15,109) X0,t0
            ENDDO
            Call RKF78(YHC4New,TT-t0,t0,X0,EE,1D0,8) ! 积分积满，恰好积分到tf
            if (MOD(I,sepflag)==0) write(15,*) ''
            
            ! 更新状态量
            flag = 2
            XX(1) = XX(1) - step
        ENDIF
        
        ! 检测错误信息
203     flagindex(I) = flag
        IF (I>=3) THEN
            IF (flagindex(I) == flagindex(I-2) .AND. flagindex(I) /= flagindex(I-1)) THEN
                write(*,*) 'err: both methods are failed. locate: Main'
                exit
            ENDIF
        ENDIF
        
    ENDDO
    close(13)
    close(15)
    
    call system_clock(time_end)
    PRINT *, 'Time of operation was ', (time_end - time_begin)/1D3, ' seconds'
END