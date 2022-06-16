!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! @file: YHCModule.f90
! @brief: 此module里涵盖多个右函数
! @author: Wang Hai-Shuo
! @time: 2020.3.3
! @pram: global 系统参数设置头文件
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    SUBROUTINE YHC2(t0,t1,X,Y)
    ! @sub YHC2 四维状态方程
    ! @pram: t0 初始时间
    ! @pram: t1 结束时间
    ! @pram: X  初始状态矢量   
    ! @pram: Y  输出状态矢量
    use global
    REAL(8) :: t0,t1,x(4),y(4),S,theta,S1,theta1

    S = X(1)
    theta = X(2)
    S1 = X(3)
    theta1 = X(4)
    
    Y(1) = S1
    Y(2) = theta1
    Y(3) = S * (IzB * theta1 + K) ** 2D0/(IzB + m * S ** 2D0) ** 2D0 - 1D0 / S ** 2D0 - 3D0 /2D0 / S ** 4D0 * (A1 + A2 * dcos(2D0*theta))
    Y(4) = -2D0 * S1 * (IzB * theta1 + K) / S / (IzB + m * S ** 2D0) - A2 * dsin(2D0 * theta) * (m * S ** 2D0 + IzB) / S ** 5D0 / IzB
    
    end
    

    
    SUBROUTINE YHC2B(t0,t1,X,Y)
    use global
    REAL(8) :: t0,t1,x(4),y(4),S,theta,S1,theta1
    REAL(8) :: ConstB,ThY
    
    CALL Byorp(ConstB,ThY)

    S = X(1)
    theta = X(2)
    S1 = X(3)
    theta1 = X(4)
    
    Y(1) = S1
    Y(2) = theta1
    Y(3) = S * (IzB * theta1 + K) ** 2D0/(IzB + m * S ** 2D0) ** 2D0 - 1D0 / S ** 2D0 - 3D0 /2D0 / S ** 4D0 * (A1 + A2 * dcos(2D0*theta))
    Y(4) = -2D0 * S1 * (IzB * theta1 + K) / S / (IzB + m * S ** 2D0) - A2 * dsin(2D0 * theta) * (m * S ** 2D0 + IzB) / S ** 5D0 / IzB + 1D5*ConstB/(m*S)
    
    end
    
    SUBROUTINE YHC3(t0,t1,X,Y)
    ! @sub YHC3 六维状态方程
    ! @pram: t0 初始时间
    ! @pram: t1 结束时间
    ! @pram: X  初始状态矢量   
    ! @pram: Y  输出状态矢量
    use global
    REAL(8) :: t0,t1,x(6),y(6),C,S,Theta,thetaB,S1,Theta1,thetaB1

    S = X(1)
    Theta = X(2)
    thetaB = X(3)
    S1 = X(4)
    Theta1 = X(5)
    thetaB1 = X(6)
    
    Y(1) = S1
    Y(2) = Theta1
    Y(3) = thetaB1
    Y(4) = S * (Theta1) ** 2D0 - 1D0 / S ** 2D0 - (3D0 / 2D0) * (A1 + A2 * dcos(2D0 * Theta - 2D0 * thetaB)) / S ** 4D0
    Y(5) = -A2 * dsin(2D0 * Theta - 2D0 * thetaB) / S ** 5D0 - 2D0 * S1 * Theta1 / S
    Y(6) = m * A2 * dsin(2D0 * Theta - 2D0 * thetaB) / (IzB * S ** 3D0)
    end
    
    SUBROUTINE YHC3B(t0,t1,X,Y)
    ! @sub YHC3B 六维状态方程+BYORP效应
    ! @pram: t0 初始时间
    ! @pram: t1 结束时间
    ! @pram: X  初始状态矢量   
    ! @pram: Y  输出状态矢量
    use global
    REAL(8) :: t0,t1,x(6),y(6),C,S,Theta,thetaB,S1,Theta1,thetaB1
    REAL(8) :: ConstB,ThY
    
    CALL Byorp(ConstB,ThY)
    S = X(1)
    Theta = X(2)
    thetaB = X(3)
    S1 = X(4)
    Theta1 = X(5)
    thetaB1 = X(6)
    
    Y(1) = S1
    Y(2) = Theta1
    Y(3) = thetaB1
    Y(4) = S * (Theta1) ** 2D0 - 1D0 / S ** 2D0 - (3D0 / 2D0) * (A1 + A2 * dcos(2D0 * Theta - 2D0 * thetaB)) / S ** 4D0
    Y(5) = -A2 * dsin(2D0 * Theta - 2D0 * thetaB) / S ** 5D0 - 2D0 * S1 * Theta1 / S - ConstB / (S * m)
    Y(6) = m * A2 * dsin(2D0 * Theta - 2D0 * thetaB) / (IzB * S ** 3D0)
    end
    
    SUBROUTINE YHC3T(t0,t1,X,Y)
    ! @sub YHC3T 六维状态方程+tidal效应
    ! @pram: t0 初始时间
    ! @pram: t1 结束时间
    ! @pram: X  初始状态矢量   
    ! @pram: Y  输出状态矢量
    use global
    REAL(8) :: t0,t1,x(6),y(6),C,S,Theta,thetaB,S1,Theta1,thetaB1,thetaA1
    REAL(8) :: thetaB_tid,thetaA_tid
    
    S = X(1)
    Theta = X(2)
    thetaB = X(3)
    S1 = X(4)
    Theta1 = X(5)
    thetaB1 = X(6)
    thetaA1 = thetaA1_Const
    
    CALL Tide(S,Theta1,thetaA1,thetaB1,thetaA_tid,thetaB_tid)
    
    Y(1) = S1
    Y(2) = Theta1
    Y(3) = thetaB1
    Y(4) = S * (Theta1) ** 2D0 - 1D0 / S ** 2D0 - (3D0 / 2D0) * (A1 + A2 * dcos(2D0 * Theta - 2D0 * thetaB)) / S ** 4D0
    Y(5) = -A2 * dsin(2D0 * Theta - 2D0 * thetaB) / S ** 5D0 - 2D0 * S1 * Theta1 / S  - IzB * thetaB_tid/(m*S**2D0) - IzA * thetaA_tid/(m*S**2D0)
    Y(6) = m * A2 * dsin(2D0 * Theta - 2D0 * thetaB) / (IzB * S ** 3D0) + thetaB_tid
    end
     
    SUBROUTINE YHC4(t0,t1,X,Y)
    ! @sub YHC4 八维状态方程
    ! @pram: t0 初始时间
    ! @pram: t1 结束时间
    ! @pram: X  初始状态矢量   
    ! @pram: Y  输出状态矢量
    use global
    REAL(8) :: t0,t1,x(8),y(8),C,S,Theta,thetaA,thetaB,S1,Theta1,thetaA1,thetaB1

    S = X(1)
    Theta = X(2)
    thetaA = X(3)
    thetaB = X(4)
    S1 = X(5)
    Theta1 = X(6)
    thetaA1 = X(7)
    thetaB1 = X(8)
    
    Y(1) = S1
    Y(2) = Theta1
    Y(3) = thetaA1
    Y(4) = thetaB1
    Y(5) = S * (Theta1) ** 2D0 - 1D0 / S ** 2D0 - (3D0 / 2D0) * (A1 + A2 * dcos(2D0 * Theta - 2D0 * thetaB) + A3 * dcos(2D0 * Theta - 2D0 * thetaA)) / S ** 4D0
    Y(6) = -A2 * dsin(2D0 * Theta - 2D0 * thetaB) / S ** 5D0 - A3 * dsin(2D0 * Theta - 2D0 * thetaA) / S ** 5D0 - 2D0 * S1 * Theta1 / S
    Y(7) = m * A3 * dsin(2D0 * Theta - 2D0 * thetaA) / (IzA * S ** 3D0)
    Y(8) = m * A2 * dsin(2D0 * Theta - 2D0 * thetaB) / (IzB * S ** 3D0)
    end
    
    SUBROUTINE YHC4T(t0,t1,X,Y)
    ! @sub YHC4T 八维状态方程+tidal效应
    ! @pram: t0 初始时间
    ! @pram: t1 结束时间
    ! @pram: X  初始状态矢量   
    ! @pram: Y  输出状态矢量
    use global
    REAL(8) :: t0,t1,x(8),y(8),C,S,Theta,thetaA,thetaB,S1,Theta1,thetaA1,thetaB1
    REAL(8) :: thetaB_tid,thetaA_tid

    S = X(1)
    Theta = X(2)
    thetaA = X(3)
    thetaB = X(4)
    S1 = X(5)
    Theta1 = X(6)
    thetaA1 = X(7)
    thetaB1 = X(8)
    
    CALL Tide(S,Theta1,thetaA1,thetaB1,thetaA_tid,thetaB_tid)
    
    Y(1) = S1
    Y(2) = Theta1
    Y(3) = thetaA1
    Y(4) = thetaB1
    Y(5) = S * (Theta1) ** 2D0 - 1D0 / S ** 2D0 - (3D0 / 2D0) * (A1 + A2 * dcos(2D0 * Theta - 2D0 * thetaB) + A3 * dcos(2D0 * Theta - 2D0 * thetaA)) / S ** 4D0
    Y(6) =  - 2D0 * S1 * Theta1 / S -A2 * dsin(2D0 * Theta - 2D0 * thetaB) / S ** 5D0 - A3 * dsin(2D0 * Theta - 2D0 * thetaA) / S ** 5D0 - IzB * thetaB_tid/(m*S**2D0) - IzA * thetaA_tid/(m*S**2D0)
    Y(7) = m * A3 * dsin(2D0 * Theta - 2D0 * thetaA) / (IzA * S ** 3D0) + thetaA_tid
    Y(8) = m * A2 * dsin(2D0 * Theta - 2D0 * thetaB) / (IzB * S ** 3D0) + thetaB_tid
    end
    
    SUBROUTINE YHC4B(t0,t1,X,Y)
    use global
    REAL(8) :: t0,t1,x(8),y(8),C,S,Theta,thetaA,thetaB,S1,Theta1,thetaA1,thetaB1
    REAL(8) :: ConstB,ThY
    
    CALL Byorp(ConstB,ThY)

    S = X(1)
    Theta = X(2)
    thetaA = X(3)
    thetaB = X(4)
    S1 = X(5)
    Theta1 = X(6)
    thetaA1 = X(7)
    thetaB1 = X(8)
    
    Y(1) = S1
    Y(2) = Theta1
    Y(3) = thetaA1
    Y(4) = thetaB1
    Y(5) = S * (Theta1) ** 2D0 - 1D0 / S ** 2D0 - (3D0 / 2D0) * (A1 + A2 * dcos(2D0 * Theta - 2D0 * thetaB) + A3 * dcos(2D0 * Theta - 2D0 * thetaA)) / S ** 4D0
    Y(6) = -A2 * dsin(2D0 * Theta - 2D0 * thetaB) / S ** 5D0 - A3 * dsin(2D0 * Theta - 2D0 * thetaA) / S ** 5D0 - 2D0 * S1 * Theta1 / S + ConstB / (S * m)
    Y(7) = m * A3 * dsin(2D0 * Theta - 2D0 * thetaA) / (IzA * S ** 3D0)+ThY*alphaA/(1D0-mu)
    Y(8) = m * A2 * dsin(2D0 * Theta - 2D0 * thetaB) / (IzB * S ** 3D0)+ThY*alphaB/mu
    end
    
    SUBROUTINE YHC4BT(t0,t1,X,Y)
    use global
    REAL(8) :: t0,t1,x(8),y(8)
    REAL(8) :: S,Theta,thetaA,thetaB,S1,Theta1,thetaA1,thetaB1
    REAL(8) :: thetaB_tid,thetaA_tid,ConstB,ThY

    S = X(1)
    Theta = X(2)
    thetaA = X(3)
    thetaB = X(4)
    S1 = X(5)
    Theta1 = X(6)
    thetaA1 = X(7)
    thetaB1 = X(8)
    
    CALL Tide(S,Theta1,thetaA1,thetaB1,thetaA_tid,thetaB_tid)
    CALL Byorp(ConstB,ThY)
    
    Y(1) = S1
    Y(2) = Theta1
    Y(3) = thetaA1
    Y(4) = thetaB1
    Y(5) = S * (Theta1) ** 2D0 - 1D0 / S ** 2D0 - (3D0 / 2D0) * (A1 + A2 * dcos(2D0 * Theta - 2D0 * thetaB) + A3 * dcos(2D0 * Theta - 2D0 * thetaA)) / S ** 4D0
    Y(6) =  - 2D0 * S1 * Theta1 / S -A2 * dsin(2D0 * Theta - 2D0 * thetaB) / S ** 5D0 - A3 * dsin(2D0 * Theta - 2D0 * thetaA) / S ** 5D0 - IzB * thetaB_tid/(m*S**2D0) - IzA * thetaA_tid/(m*S**2D0) + ConstB / (S * m)
    Y(7) = m * A3 * dsin(2D0 * Theta - 2D0 * thetaA) / (IzA * S ** 3D0) + thetaA_tid + ThY*alphaA/(1D0-mu)
    Y(8) = m * A2 * dsin(2D0 * Theta - 2D0 * thetaB) / (IzB * S ** 3D0) + thetaB_tid + ThY*alphaB/mu
    end
    
    SUBROUTINE YHC3BT(t0,t1,X,Y)
    use global
    REAL(8) :: t0,t1,x(6),y(6),C,S,Theta,thetaB,S1,Theta1,thetaB1,thetaA1
    REAL(8) :: thetaB_tid,thetaA_tid,ConstB,ThY
102 format(1X, E25.18,1X,E25.18) 
    
    S = X(1)
    Theta = X(2)
    thetaB = X(3)
    S1 = X(4)
    Theta1 = X(5)
    thetaB1 = X(6)
    thetaA1 = thetaA1_Const
    
    CALL Tide(S,Theta1,thetaA1,thetaB1,thetaA_tid,thetaB_tid)
    CALL Byorp(ConstB,ThY)
    
    Y(1) = S1
    Y(2) = Theta1
    Y(3) = thetaB1
    Y(4) = S * (Theta1) ** 2D0 - 1D0 / S ** 2D0 - (3D0 / 2D0) * (A1 + A2 * dcos(2D0 * Theta - 2D0 * thetaB)) / S ** 4D0
    Y(5) = -A2 * dsin(2D0 * Theta - 2D0 * thetaB) / S ** 5D0 - 2D0 * S1 * Theta1 / S - IzA * thetaA_tid/(m*S**2D0) - ConstB / (S * m) !- IzB * thetaB_tid/(m*S**2D0)
    Y(6) = m * A2 * dsin(2D0 * Theta - 2D0 * thetaB) / (IzB * S ** 3D0) !+ thetaB_tid
    !write(*,*) 'BYORP = ',C / (S * m),'Tidal = ',IzB * thetaB_tid/(m*S**2D0)
    !write(13,102) IzB * thetaB_tid/(m*S**2D0),C / (S * m)
    end
    
    SUBROUTINE YHC4New(t0,t1,X,Y)
    use global
    REAL(8) :: t0,t1,x(8),y(8),S,theta,phi,thetaB,S1,theta1,phi1,thetaB1

    S = X(1)
    theta = X(2)
    phi = X(3)
    thetaB = X(4)
    S1 = X(5)
    theta1 = X(6)
    phi1 = X(7)
    thetaB1 = X(8)
    
    Y(1) = S1
    Y(2) = theta1
    Y(3) = phi1
    Y(4) = thetaB1
    Y(5) = S * (theta1+thetaB1) ** 2D0 - 1D0 / S ** 2D0 - (3D0 / 2D0) * (A1 + A2 * dcos(2D0 * theta) + A3 * dcos(2D0 * (theta - phi))) / S ** 4D0
    Y(6) = -A2 * dsin(2D0 * theta) / S ** 5D0 - A3 * dsin(2D0 * (theta - phi)) / S ** 5D0 - 2D0 * S1 * (theta1 + thetaB1) / S - m * A2 * dsin(2D0 * theta) / (IzB * S ** 3D0)
    Y(7) = m * A3 * dsin(2D0 * (theta - phi)) / (IzA * S ** 3D0) - m * A2 * dsin(2D0 * theta) / (IzB * S ** 3D0)
    Y(8) = m * A2 * dsin(2D0 * theta) / (IzB * S ** 3D0)
    end