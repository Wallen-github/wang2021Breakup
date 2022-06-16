!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! @file: ES_Method.f90
! @brief: 此module里包括多个常用的方法
! @author: Wang Hai-Shuo
! @time: 2020.3.3
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE Propagate4(FILEID,YHC,t0,tf,tspan,err,hh,Num,x0)
    use global
    USE FindRoots
    external :: YHC
    REAL(8) :: t0,tf,err,step,Energy,Angular,w0(2),hh,tspan,x0(8),EE,the(10000000),max
    REAL(8) :: S,Theta,thetaA,thetaB,S1,Theta1,thetaA1,thetaB1
    integer :: FILEID
102 format(1X, E25.18,1X,E25.18)    
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
104 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
105 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
106 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
109 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18)
111 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18)
    do I=1,Num 
        J = 1
        do while (t0<tf)
            S = X0(1)
            Theta = X0(2)
            thetaA = X0(3)
            thetaB = X0(4)
            S1 = X0(5)
            Theta1 = X0(6)
            thetaA1 = X0(7)
            thetaB1 = X0(8)
            Energy = m / 2D0 * (S1 ** 2D0 + S ** 2D0 * Theta1 ** 2D0) + 1D0 / 2D0 * IzB * thetaB1 ** 2D0 + 1D0 / 2D0 * IzA * thetaA1 ** 2D0 - m / S - m / (2D0 * S ** 3D0) * (A1 + A2 * dcos(2D0 * Theta - 2D0 * thetaB) + A3 * dcos(2D0 * Theta - 2D0 * thetaA))
            Angular = m * S ** 2D0 * Theta1 + IzA * thetaA1 + IzB * thetaB1
            the(J) = theta-thetaB
            J = J+1
            Call RKF78(YHC,hh,t0,x0,EE,err,8)
            write(FILEID,111) x0,energy,Angular,t0*Tunit 
        enddo
        max = maxval(abs(the(1:J-1)))
        if (max>PI .OR. S<1D0) then
            write(*,*) 'breakup'
            GO TO 21
        end if
        write(*,*) 'I = ',I ,'/',Num, ';','T = ',t0*Tunit,'year',';','theta = ',max
        !write(FILEID,111) x0,energy,Angular,t0*Tunit 
        tf=t0+tspan
    enddo
21  write(*,*) 'DONE'

    end   
    
    SUBROUTINE Propagate4BT(FILEID,YHC,t0,tf,tspan,err,hh,Num,x0)
    use global
    USE FindRoots
    external :: YHC
    REAL(8) :: t0,tf,err,step,w0(2),hh,tspan,x0(8),EE,the(10000000),max,axis,ecc,axis_sum,ecc_sum,Energy,Angular,placeholder
    REAL(8) :: S,Theta,thetaA,thetaB,S1,Theta1,thetaA1,thetaB1
    integer :: FILEID
    REAL(8) :: thetaB_tid,thetaA_tid,ConstB,ThY
    
102 format(1X, E25.18,1X,E25.18)    
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
104 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
105 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
106 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
109 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18)
111 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18)
113 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X, E25.18,1X,E25.18)    
115 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X, E25.18,1X,E25.18)
    do I=1,Num 
        J = 1
        axis_sum = 0
        ecc_sum = 0
        do while (t0<tf)
            S = X0(1)
            Theta = X0(2)
            thetaA = X0(3)
            thetaB = X0(4)
            S1 = X0(5)
            Theta1 = X0(6)
            thetaA1 = X0(7)
            thetaB1 = X0(8)
            Energy = m / 2D0 * (S1 ** 2D0 + S ** 2D0 * Theta1 ** 2D0) + 1D0 / 2D0 * IzB * thetaB1 ** 2D0 + 1D0 / 2D0 * IzA * thetaA1 ** 2D0 - m / S - m / (2D0 * S ** 3D0) * (A1 + A2 * dcos(2D0 * Theta - 2D0 * thetaB) + A3 * dcos(2D0 * Theta - 2D0 * thetaA))
            Angular = m * S ** 2D0 * Theta1 + IzA * thetaA1 + IzB * thetaB1
            the(J) = Theta-thetaB

            !CALL Tide(S,Theta1,thetaA1,thetaB1,thetaA_tid,thetaB_tid)
            !CALL Byorp(ConstB,ThY)
            CALL mco_x2el (1D0,x0(1)*dcos(x0(2)),x0(1)*dsin(x0(2)),0D0,x0(4)*dcos(x0(2))-x0(1)*dsin(x0(2))*x0(5),x0(4)*dsin(x0(2))+x0(1)*dcos(x0(2))*x0(5),0D0,axis,ecc,placeholder,placeholder,placeholder,placeholder)
            !CALL GETOE4(x0,axis,ecc)
            axis_sum = axis_sum + axis
            ecc_sum = ecc_sum + ecc
            
            J = J+1
            Call RKF78(YHC,hh,t0,x0,EE,err,8)
            !write(FILEID,113) x0, IzB * thetaB_tid/(m*S**2D0), IzA * thetaA_tid/(m*S**2D0),t0*Tunit,axis,ecc,Energy,Angular
        enddo
        max = maxval(abs(the(1:J-1)))
        if (max>PI .OR. S<1D0) then
            write(*,*) 'breakup'
            GO TO 21
        end if
        write(*,*) 'I = ',I ,'/',Num, ';','T = ',t0*Tunit,'year',';','theta = ',max
        write(FILEID,113) x0,t0*Tunit,axis_sum/(J-1),ecc_sum/(J-1),Energy,Angular 
        tf=t0+tspan
    enddo
21  write(*,*) 'DONE'

    end

    SUBROUTINE Propagate3(FILEID,YHC,t0,tf,tspan,err,hh,Num,x0)
    ! @sub Propagate3 三个自由度运动方程的积分器
    ! @outfile FILEID： 输出数据文件FILEID，所有积分状态量，能量，角动量和以年为单位的积分时间x0,energy,Angular,t0*Tunit
    ! @pram: YHC    右函数 
    ! @pram: t0     积分初始时间
    ! @pram: tf     积分结束时间
    ! @pram: tspan  积分时间间隔（内循环，不输出数据），间隔tspan之后输出一次数据
    ! @pram: err    误差限
    ! @pram: hh     初始步长
    ! @pram: Num    积分时间间隔（外循环，输出数据），总计输出Num个数据
    ! @pram: x0     积分初始状态量和结束时输出状态量
    use global
    USE FindRoots
    external :: YHC
    REAL(8) :: t0,tf,err,step,Energy,Angular,w0(2),hh,tspan,x0(6),EE,the(10000000),max,placeholder
    REAL(8) :: S,Theta,thetaB,S1,Theta1,thetaB1
    integer :: FILEID
102 format(1X, E25.18,1X,E25.18)    
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
104 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
105 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
106 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
109 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18)
111 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18)
113 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X, E25.18,1X,E25.18)
    do I=1,Num
        J = 1
        do while (t0<tf)
            S = x0(1)
            Theta = x0(2)
            thetaB = x0(3)
            S1 = x0(4)
            Theta1 = x0(5)
            thetaB1 = x0(6)
            Energy = m / 2D0 * (S1 ** 2D0 + S ** 2D0 * Theta1 ** 2D0) + 1D0 / 2D0 * IzB * thetaB1 ** 2D0 - m / S - m / (2D0 * S ** 3D0) * (A1 + A2 * dcos(2D0 * Theta - 2D0 * thetaB))
            Angular = m * S ** 2D0 * Theta1 + IzB * thetaB1
            the(J) = Theta - thetaB
            J = J + 1
            Call RKF78(YHC,hh,t0,x0,EE,err,6)
            !write(FILEID,109) x0,energy,Angular,t0*Tunit 
        enddo
        max = maxval(abs(the(1:J-1)))
        if (max>PI .OR. S<1D0) then
            write(*,*) 'breakup'
            GO TO 21
        end if
        write(*,*) 'I = ',I ,'/',Num, ';','T = ',t0*Tunit,'year',';','theta = ',max
        write(FILEID,109) x0,energy,max,t0*Tunit 
        tf=t0+tspan
    enddo
21  write(*,*) 'DONE'

    end
    
    SUBROUTINE Propagate3BT(FILEID,YHC,t0,tf,tspan,err,hh,Num,x0)
    use global
    USE FindRoots
    external :: YHC
    REAL(8) :: t0,tf,err,step,Energy,Angular,w0(2),hh,tspan,x0(6),EE,the(10000000),max,axis,ecc,axis_sum,ecc_sum,placeholder
    REAL(8) :: S,Theta,thetaB,S1,Theta1,thetaB1,thetaA1
    integer :: FILEID
    REAL(8) :: thetaB_tid,thetaA_tid,ConstB,ThY
    
102 format(1X, E25.18,1X,E25.18)    
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
104 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
105 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
106 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
109 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18)
110 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
111 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18)
113 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X, E25.18,1X,E25.18)
    do I=1,Num 
        J = 1
        axis_sum = 0
        ecc_sum = 0
        do while (t0<tf)
            S = x0(1)
            Theta = x0(2)
            thetaB = x0(3)
            S1 = x0(4)
            Theta1 = x0(5)
            thetaB1 = x0(6)
            thetaA1 = thetaA1_Const
            Energy = m / 2D0 * (S1 ** 2D0 + S ** 2D0 * Theta1 ** 2D0) + 1D0 / 2D0 * IzB * thetaB1 ** 2D0 - m / S - m / (2D0 * S ** 3D0) * (A1 + A2 * dcos(2D0 * Theta - 2D0 * thetaB))
            Angular = m * S ** 2D0 * Theta1 + IzB * thetaB1
            the(J) = Theta - thetaB
    
            !CALL Tide(S,Theta1,thetaA1,thetaB1,thetaA_tid,thetaB_tid)
            !CALL Byorp(ConstB,ThY)
            CALL mco_x2el (1D0,x0(1)*dcos(x0(2)),x0(1)*dsin(x0(2)),0D0,x0(4)*dcos(x0(2))-x0(1)*dsin(x0(2))*x0(5),x0(4)*dsin(x0(2))+x0(1)*dcos(x0(2))*x0(5),0D0,axis,ecc,placeholder,placeholder,placeholder,placeholder)
            axis_sum = axis_sum + axis
            ecc_sum = ecc_sum + ecc
            J = J + 1
            
            Call RKF78(YHC,hh,t0,x0,EE,err,6)
            !write(FILEID,113) x0, IzB * thetaB_tid/(m*S**2D0),IzA * thetaA_tid/(m*S**2D0), t0*Tunit, axis, ecc ,Energy,Angular 
        enddo
        max = maxval(abs(the(1:J-1)))
        if (max>PI .OR. S<1D0) then
            write(*,*) 'breakup'
            GO TO 21
        end if
        write(*,*) 'I = ',I ,'/',Num, ';','T = ',t0*Tunit,'year',';','axis = ',axis_sum/(J-1)
        write(FILEID,111) x0, t0*Tunit, axis, ecc ,Energy,Angular 
        tf=t0+tspan
    enddo
21  write(*,*) 'DONE'

    end
    
    SUBROUTINE Propagate2(FILEID,YHC,t0,tf,tspan,err,hh,Num,x0)
    ! @sub Propagate2 两个自由度运动方程的积分器
    ! @outfile FILEID： 输出数据文件FILEID，所有积分状态量，能量，角动量和以年为单位的积分时间x0,energy,t0*Tunit,max
    ! @pram: YHC    右函数 
    ! @pram: t0     积分初始时间
    ! @pram: tf     积分结束时间
    ! @pram: tspan  积分时间间隔（内循环，不输出数据），间隔tspan之后输出一次数据
    ! @pram: err    误差限
    ! @pram: hh     初始步长
    ! @pram: Num    积分时间间隔（外循环，输出数据），总计输出Num个数据
    ! @pram: x0     积分初始状态量和结束时输出状态量
    use global
    USE FindRoots
    external :: YHC
    integer :: FILEID
    REAL(8) :: t0,tf,err,step,Energy,w0(2),hh,tspan,x0(4),EE,the(10000000),max
    REAL(8) :: S,theta,S1,theta1
102 format(1X, E25.18,1X,E25.18)    
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
104 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
105 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
106 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
107 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
    do I=1,Num 
        J = 1
        do while (t0<tf)
            S = x0(1)
            theta = x0(2)
            S1 = x0(3)
            theta1 = x0(4)
            Energy = 1D0 / 2D0 * K ** 2D0 / (m * S ** 2 + IzB) + 1D0 / 2D0 * IzB * m * S ** 2D0 * theta1 ** 2D0 / (IzB + m * S ** 2D0) + 1D0 / 2D0 * m * S1 ** 2D0 - m / S - m * (A1 + A2 * dcos(2 * theta)) / (2 * S ** 3)
            the(J) = theta
            J = J+1
            Call RKF78(YHC,hh,t0,x0,EE,err,4)
            write(FILEID,106) x0,Energy,t0*Tunit
        enddo
        max = maxval(abs(the(1:J-1)))
        if (max>PI) then
            write(*,*) 'breakup'
            GO TO 21
        end if
        write(*,*) 'I = ',I ,'/',Num, ';','T = ',t0*Tunit,'year',';'
        !write(FILEID,107) x0,energy,t0*Tunit,max
        tf=t0+tspan
    enddo
21  write(*,*) 'DONE'

    end

    SUBROUTINE Getx03(index,xi0,eta0,x0,w0)
    ! @sub Getx03   求解三自由度的初始状态量
    ! @pram: index  长短周期轨道指示器 
    ! @pram: xi0    初始横向振幅
    ! @pram: eta0   初始纵向振幅
    ! @pram: x0     输出状态量
    ! @pram: w0     输出长短周期轨道频率
    use global
    integer(4) :: index
    REAL(8) :: a42,a43,a31,a34,w0(2),x0(6),xi0,eta0,x1,x2,x3,x4,s1,s2
    
    a42 = -2D0 * A2 * (S0 ** 2D0 * m + IzB) / (IzB * S0 ** 5D0)
    a43 = -2D0 * K / (S0 * (S0 ** 2D0 * m + IzB))
    a31 = K ** 2D0 / (S0 ** 2D0 * m + IzB) ** 2D0 - 4D0 * S0 ** 2D0 * K ** 2D0 * m / (S0 ** 2D0 * m + IzB) ** 3D0 + 2D0 / S0 ** 3D0 + (6D0 * (A1 + A2)) / S0 ** 5D0
    a34 = 2D0 * S0 * K * IzB / (S0 ** 2D0 * m + IzB) ** 2D0
    s1 = (a31 + a42 + a34 * a43) / 2D0 + dsqrt((a31 + a42 + a34 * a43) ** 2D0 - 4D0 * a31 * a42) / 2D0
    s2 = (a31 + a42 + a34 * a43) / 2D0 - dsqrt((a31 + a42 + a34 * a43) ** 2D0 - 4D0 * a31 * a42) / 2D0
    w0(1) = dsqrt(-s1) ! real(aimag(sqrt(s1)),kind=8)
    w0(2) = dsqrt(-s2) ! real(aimag(sqrt(s2)),kind=8)
    alpha = -a43 * w0(index) / (a42 + w0(index) ** 2D0)
    x1 = S0 + xi0
    x2 = eta0
    x3 = eta0 * w0(index) / alpha
    x4 = -xi0 * alpha * w0(index)
    x0 = (/x1, x2, 0D0, x3, x4+(K - m * x1 ** 2D0 * x4) / (m * x1 ** 2D0 + IzB), (K - m * x1 ** 2D0 * x4) / (m * x1 ** 2D0 + IzB)/)
    return
    end

    SUBROUTINE Getx02(index,xi0,eta0,x0,w0)
    ! @sub Getx02   求解两自由度的初始状态量
    ! @pram: index  长短周期轨道指示器 
    ! @pram: xi0    初始横向振幅
    ! @pram: eta0   初始纵向振幅
    ! @pram: x0     输出状态量
    ! @pram: w0     输出长短周期轨道频率
    use global
    integer(4) :: index
    REAL(8) :: a42,a43,a31,a34,w0(2),x0(4),xi0,eta0,x1,x2,x3,x4,s1,s2
    
    a42 = -2D0 * A2 * (S0 ** 2D0 * m + IzB) / (IzB * S0 ** 5D0)
    a43 = -2D0 * K / (S0 * (S0 ** 2D0 * m + IzB))
    a31 = K ** 2D0 / (S0 ** 2D0 * m + IzB) ** 2D0 - 4D0 * S0 ** 2D0 * K ** 2D0 * m / (S0 ** 2D0 * m + IzB) ** 3D0 + 2D0 / S0 ** 3D0 + (6D0 * (A1 + A2)) / S0 ** 5D0
    a34 = 2D0 * S0 * K * IzB / (S0 ** 2D0 * m + IzB) ** 2D0
    s1 = (a31 + a42 + a34 * a43) / 2D0 + dsqrt((a31 + a42 + a34 * a43) ** 2D0 - 4D0 * a31 * a42) / 2D0
    s2 = (a31 + a42 + a34 * a43) / 2D0 - dsqrt((a31 + a42 + a34 * a43) ** 2D0 - 4D0 * a31 * a42) / 2D0
    w0(1) = dsqrt(-s1) ! real(aimag(sqrt(s1)),kind=8)
    w0(2) = dsqrt(-s2) ! real(aimag(sqrt(s2)),kind=8)
    alpha = -a43 * w0(index) / (a42 + w0(index) ** 2D0)
    x1 = S0 + xi0
    x2 = eta0
    x3 = eta0 * w0(index) / alpha
    x4 = -xi0 * alpha * w0(index)
    x0 = (/x1, x2, x3, x4/)
    return
    end
    
    SUBROUTINE GetAB(alpha,beta,phi1,phi2,x0)
    ! @sub GetAB    通过线性化参数求解两个自由度的状态量
    ! @pram: alpha  线性化参数，横向最大振幅
    ! @pram: beta   线性化参数，自由摆动角最大振幅
    ! @pram: phi1   线性化参数
    ! @pram: phi2   线性化参数
    ! @pram: x0     输出状态量
    use global
    REAL(8) :: a42,a43,a31,a34,w0(2),x0(4),xi0,eta0,x1,x2,x3,x4,s1,s2,alpha,beta,phi1,phi2
    
    a42 = -2D0 * A2 * (S0 ** 2D0 * m + IzB) / (IzB * S0 ** 5D0)
    a43 = -2D0 * K / (S0 * (S0 ** 2D0 * m + IzB))
    a31 = K ** 2D0 / (S0 ** 2D0 * m + IzB) ** 2D0 - 4D0 * S0 ** 2D0 * K ** 2D0 * m / (S0 ** 2D0 * m + IzB) ** 3D0 + 2D0 / S0 ** 3D0 + (6D0 * (A1 + A2)) / S0 ** 5D0
    a34 = 2D0 * S0 * K * IzB / (S0 ** 2D0 * m + IzB) ** 2D0
    s1 = (a31 + a42 + a34 * a43) / 2D0 + dsqrt((a31 + a42 + a34 * a43) ** 2D0 - 4D0 * a31 * a42) / 2D0
    s2 = (a31 + a42 + a34 * a43) / 2D0 - dsqrt((a31 + a42 + a34 * a43) ** 2D0 - 4D0 * a31 * a42) / 2D0
    w0(1) = dsqrt(-s1) ! real(aimag(sqrt(s1)),kind=8)
    w0(2) = dsqrt(-s2) ! real(aimag(sqrt(s2)),kind=8)
    alpha1 = -a43 * w0(1) / (a42 + w0(1) ** 2D0)
    alpha2 = -a43 * w0(2) / (a42 + w0(2) ** 2D0)
    x1 = S0 + alpha*dcos(phi2)+beta*dcos(phi1)/alpha1
    x2 = -alpha2*alpha*sin(phi2)-beta*dsin(phi1)
    x3 = -w0(2)*alpha*dsin(phi2)-beta*w0(1)*dsin(phi1)/alpha1
    x4 = -w0(2)*alpha2*alpha*dcos(phi2)-w0(1)*beta*dcos(phi1)
    x0 = (/x1, x2, x3, x4/)
    return
    end
    
    SUBROUTINE GetAB_re(x0,alpha,beta,phi1,phi2)
    ! @sub GetAB_re    通过两个自由度的状态量求解线性化模型参数
    ! @pram: alpha  线性化参数，横向最大振幅
    ! @pram: beta   线性化参数，自由摆动角最大振幅
    ! @pram: phi1   线性化参数
    ! @pram: phi2   线性化参数
    ! @pram: x0     输出状态量
    use global
    REAL(8) :: a42,a43,a31,a34,w0(2),x0(4),xi,eta,xi1,eta1,s1,s2,alpha,beta,phi1,phi2
    
    a42 = -2D0 * A2 * (S0 ** 2D0 * m + IzB) / (IzB * S0 ** 5D0)
    a43 = -2D0 * K / (S0 * (S0 ** 2D0 * m + IzB))
    a31 = K ** 2D0 / (S0 ** 2D0 * m + IzB) ** 2D0 - 4D0 * S0 ** 2D0 * K ** 2D0 * m / (S0 ** 2D0 * m + IzB) ** 3D0 + 2D0 / S0 ** 3D0 + (6D0 * (A1 + A2)) / S0 ** 5D0
    a34 = 2D0 * S0 * K * IzB / (S0 ** 2D0 * m + IzB) ** 2D0
    s1 = (a31 + a42 + a34 * a43) / 2D0 + dsqrt((a31 + a42 + a34 * a43) ** 2D0 - 4D0 * a31 * a42) / 2D0
    s2 = (a31 + a42 + a34 * a43) / 2D0 - dsqrt((a31 + a42 + a34 * a43) ** 2D0 - 4D0 * a31 * a42) / 2D0
    w0(1) = dsqrt(-s1) ! real(aimag(sqrt(s1)),kind=8)
    w0(2) = dsqrt(-s2) ! real(aimag(sqrt(s2)),kind=8)
    alpha1 = -a43 * w0(1) / (a42 + w0(1) ** 2D0)
    alpha2 = -a43 * w0(2) / (a42 + w0(2) ** 2D0)
    xi = x0(1)-S0
    eta = x0(2)
    xi1 = x0(3)
    eta1 = x0(4)
    phi1 = ATAN2(((w0(2)*eta-alpha2*xi1)/(w0(2)*alpha2*xi+eta1)),((w0(1)*alpha2-alpha1*w0(2))/(w0(2)*alpha2-w0(1)*alpha1)))
    phi2 = ATAN2(((w0(1)*eta-alpha1*xi1)/(alpha1*w0(1)*xi+eta1)),((alpha1*w0(2)-alpha2*w0(1))/(alpha1*w0(1)-alpha2*w0(2))))
    alpha = dsqrt(((alpha1*w0(1)*xi+eta1)/(alpha1*w0(1)-alpha2*w0(2)))**2D0+((alpha1*xi1-w0(1)*eta)/(alpha2*w0(1)-alpha2*w0(2)))**2D0)
    beta = dsqrt(((alpha1*(w0(2)*alpha2*xi+eta1))/(alpha2*w0(2)-alpha1*w0(1)))**2D0+((alpha1*(w0(2)*eta-alpha2*xi1))/(alpha1*w0(2)-alpha2*w0(1)))**2D0)
    end
    
    SUBROUTINE Pincare(FILEID1,FILEID2,YHC,t0,tf,err,hh,x0,xmid,xx)
    ! @sub Pincare    通过两个自由度的状态量求解线性化模型参数
    ! @outfile FILEID： 输出数据文件FILEID，所有满足要求的状态量x0
    ! @pram: YHC    右函数 
    ! @pram: t0     积分初始时间
    ! @pram: tf     积分结束时间
    ! @pram: err    误差限
    ! @pram: hh     初始步长
    ! @pram: xx(4)  中间变量
    ! @pram: xmid(4)中间变量
    ! @pram: w0(2)  输出长短周期轨道频率
    ! @pram: x0(4)  输出状态量
    use global
    external :: YHC
    REAL(8) :: t0,tf,err,step,w0(2),hh,tspan,x0(4),EE,xmid(4),hh1,xx(4),ymid
    REAL(8) :: S,theta,S1,theta1
    integer :: FILEID1,FILEID2
100 format(1X, E25.18,1X,E25.18,1X,E25.18)
101 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
102 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
103 format(1X, E25.18,1X,E25.18)
    hh1 = hh
    do while (t0<tf)
23      xmid = x0
        ymid = xmid(1)*dsin(xmid(2))
        tmid = t0
22      Call RKF78(YHC, hh1, t0, x0, EE, err, 4)
        !转换到直角坐标系
        xx(1) = x0(1)*dcos(x0(2))
        xx(2) = x0(1)*dsin(x0(2))
        xx(3) = x0(3)*dcos(x0(2))-x0(1)*dsin(x0(2))*x0(4)
        xx(4) = x0(3)*dsin(x0(2))+x0(1)*dcos(x0(2))*x0(4)
        write(FILEID1,101) x0
        if (ymid*xx(2)<0 .AND. xx(4)<0) then
            if (ABS(ymid)<err) then
                write(FILEID2,101) x0
                write(*,101) x0
                hh1 = hh
                GO TO 23
                !exit
            endif
            if (hh1 < 1D-15) then
                hh1 = hh
            endif
            hh1 = hh1/1D1
            x0 = xmid
            t0 = tmid
            GOTO 22
        end if
    end do
    end
    
    SUBROUTINE GETOE3(x0,axis,ecc)
    ! @sub GETOE    由轨道的状态量求解轨道根数
    ! @pram: x0(6)     输出状态量
    ! @pram: axis      半长轴
    ! @pram: ecc       偏心率
    use global
    REAL(8) :: x0(6)
    REAL(8) :: S,Theta,thetaB,S1,Theta1,thetaB1
    REAL(8) :: axis_1,axis,ecc
    
    S = x0(1)
    Theta = x0(2)
    thetaB = x0(3)
    S1 = x0(4)
    Theta1 = x0(5)
    thetaB1 = x0(6)
    
    axis_1 = 2D0/S-S1**2D0-S**2D0*Theta1**2D0
    axis = 1D0/axis_1
    ecc = dsqrt((1-S*axis_1)**2D0+(S*S1)**2D0*axis_1)

    end
    
    SUBROUTINE GETOE4(x0,axis,ecc)
    use global
    REAL(8) :: x0(8)
    REAL(8) :: S,Theta,thetaA,thetaB,S1,Theta1,thetaA1,thetaB1
    REAL(8) :: axis_1,axis,ecc
    
    S = X0(1)
    Theta = X0(2)
    thetaA = X0(3)
    thetaB = X0(4)
    S1 = X0(5)
    Theta1 = X0(6)
    thetaA1 = X0(7)
    thetaB1 = X0(8)
    
    axis_1 = 2D0/S-S1**2D0-S**2D0*Theta1**2D0
    axis = 1D0/axis_1
    ecc = dsqrt((1-S*axis_1)**2D0+(S*S1)**2D0*axis_1)

    end
    
    SUBROUTINE Contour(FILEID,alS_max,alS_min,beD_max,beD_min,Num,tf,x4)
    USE global
    external :: YHC3B,YHC2,YHC3,YHC3T
    REAL(8),allocatable :: Sta(:,:),xi(:),eta(:)
    REAL(8) :: h1,h2,alS_max,alS_min,beD_max,beD_min,x4(4),t0,tf,EE,hh,err,Energy,Angular
    integer(4) :: Num,FILEID
    allocate(xi(Num))
    allocate(eta(Num))
    allocate(Sta(Num,Num))
100 format(1X, E25.18,1X,E25.18,1X,E25.18)
101 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
102 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
103 format(1X, E25.18,1X,E25.18)    
    err=1D-10
    hh = 1D-2
    t0 = 0D0
    MM = 4

    h1 = (alS_max-alS_min)/Num
    h2 = (beD_max-beD_min)/Num
    
    do I = 1,Num
        xi(I) = alS_min + (I-1)*h1
        do J = 1,Num
            eta(J) = beD_min + (J-1)*h2
            CALL GetAB(xi(I)*S0,eta(J)*PI/1.8D2,0D0,0D0,x4)
            t0 = 0D0
            do while (t0<tf)
                Call RKF78(YHC2,hh,t0,x4,EE,err,MM)
                if (ABS(x4(2))>PI) then
                    Sta(I,J) = 0D0
                    GO TO 20
                end if
            end do
            Sta(I,J) = 1D0
20          write(FILEID,'(1x,E25.18,$)') Sta(I,J)
        end do
        WRITE(FILEID,*) ' '
        !write(16,103) xi(I)/S0,eta(I)*180/pi
    end do
    
    
    end SUBROUTINE
    
    SUBROUTINE ReadData(FILEID,filename,unitT,x)
    ! @sub ReadData    从轨道历表中读入轨道数据
    ! @pram: FILEID    文件编号
    ! @pram: FILENAME  文件路径
    ! @pram: unitT     第几行数据
    ! @pram: x0(9)     输出状态量与能量，角动量和时间，x0,energy,Angular,t0*Tunit
    implicit none
    integer, parameter :: IRECSZ = 20
    integer :: unitT,NRECL,fileid
    character(len=26)  :: hit(3)
    REAL(8) :: x(9)
    character(len=80) :: filename
    integer error
    
    NRECL = (unitT-1)*3+1
    open(unit=FILEID, file=filename, access="direct",form="unformatted", recl=IRECSZ, status="old")
    read(FILEID, rec=NRECL, IOSTAT=error) hit
    read(hit,"(1X,E25.18)")x(1:3)
    read(FILEID, rec=NRECL+1, IOSTAT=error) hit
    read(hit,"(1X,E25.18)")x(4:6)
    read(FILEID, rec=NRECL+2, IOSTAT=error) hit
    read(hit,"(1X,E25.18)")x(7:9)
    close(FILEID)

    end
    
    SUBROUTINE FindDot_version2(FILEID,YHC,t0,tf,err,DD,hh1,x4,xree,tree,xmid)
    USE global
    external :: YHC
    integer :: FILEID
    REAL(8) :: t0,tf,err,hh1,x4(4),EE,xmid(4),tmid,DD
    REAL(8) :: xree(4),tree,hree,alpha,beta,phi1,phi2
101 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18) 
    tree = t0
    xree = x4
    hree = hh1
    do while (tree<tf)
22      xmid = xree
        tmid = tree
        Call RKF78(YHC, hree, tree, xree, EE, DD, 4)
        if (MOD(xree(2),2D0*PI)*MOD(xmid(2),2D0*PI)<0 .AND. xree(4)>0) then
            if (ABS(xree(2))<1D-14) then
                write(FILEID,101) xree
                xree(2) = 0D0
                xree(3) = 0D0
                CALL GetAB_re(xree,alpha,beta,phi1,phi2)
                write(*,101) alpha/S0,beta*18D1/PI,phi1,phi2
                !write(*,101) alpha,beta,phi1,phi2
                hree = hh1
                GOTO 22
            end if
            xree = xmid
            tree = tmid
            hree = hree/1D1 
        end if
    end do

    end
    