!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! @file: YHC2.f90
! @brief: 此module里涵盖多个右函数,相较于第一版，在球谐系数参数A1等做出调整，同时引入thetaA1
! @author: Wang Hai-Shuo
! @time: 2020.8.2
! @pram: global 系统参数设置头文件
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    SUBROUTINE YHC2(t0,t1,X,Y)
    use global
    REAL(8) :: t0,t1,x(4),y(4),S,theta,S1,theta1
    ! @sub YHC2 四维状态方程平均化系统
    ! @pram: t0 初始时间
    ! @pram: t1 结束时间
    ! @pram: X  初始状态矢量   
    ! @pram: Y  输出状态矢量
    S = X(1)
    theta = X(2)
    S1 = X(3)
    theta1 = X(4)
    
    Y(1) = S1
    Y(2) = theta1
    Y(3) = S * (IzB * theta1 + K) ** 2D0/(IzB + m * S ** 2D0) ** 2D0 - 1D0 / S ** 2D0 - 3D0 /2D0 / S ** 4D0 * (A1 + A2 * dcos(2D0*theta))
    Y(4) = -2D0 * S1 * (IzB * theta1 + K) / S / (IzB + m * S ** 2D0) - A2 * dsin(2D0 * theta) * (m * S ** 2D0 + IzB) / S ** 5D0 / IzB
    
    end
    
    SUBROUTINE YHC3(t0,t1,X,Y)
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
    
    SUBROUTINE YHC2A(t0,t1,X,Y)
    ! @sub YHC2A 四维状态方程平均化系统,这是一个错误的方程
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
    Y(3) = S * (IzB * theta1 + K - IzA*thetaA1_Const) ** 2D0/(IzB + m * S ** 2D0) ** 2D0 - 1D0 / S ** 2D0 - 3D0 / S ** 4D0 * (A1 + A2 * dcos(2D0*theta))
    Y(4) = -2D0 * S1 * (IzB * theta1 + K - IzA*thetaA1_Const) / S / (IzB + m * S ** 2D0) - 2D0 * A2 * dsin(2D0 * theta) * (m * S ** 2D0 + IzB) / S ** 5D0 / IzB
    
    end
    
    SUBROUTINE YHC3ABT(t0,t1,X,Y)
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
    Y(4) = S * (Theta1) ** 2D0 - 1D0 / S ** 2D0 - 3D0 * (A1 + A2 * dcos(2D0 * Theta - 2D0 * thetaB)) / S ** 4D0
    Y(5) = -2D0 * A2 * dsin(2D0 * Theta - 2D0 * thetaB) / S ** 5D0 - 2D0 * S1 * Theta1 / S - IzB * thetaB_tid/(m*S**2D0) !- ConstB / (S * m)
    Y(6) = 2D0 * m * A2 * dsin(2D0 * Theta - 2D0 * thetaB) / (IzB * S ** 3D0)
    end
    
    SUBROUTINE YHC3AB(t0,t1,X,Y)
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
    Y(4) = S * (Theta1) ** 2D0 - 1D0 / S ** 2D0 - 3D0 * (A1 + A2 * dcos(2D0 * Theta - 2D0 * thetaB)) / S ** 4D0
    Y(5) = -2D0 * A2 * dsin(2D0 * Theta - 2D0 * thetaB) / S ** 5D0 - 2D0 * S1 * Theta1 / S - ConstB / (S * m)
    Y(6) = 2D0 * m * A2 * dsin(2D0 * Theta - 2D0 * thetaB) / (IzB * S ** 3D0)
    end
    
    SUBROUTINE YHC3AT(t0,t1,X,Y)
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
    Y(4) = S * (Theta1) ** 2D0 - 1D0 / S ** 2D0 - 3D0 * (A1 + A2 * dcos(2D0 * Theta - 2D0 * thetaB)) / S ** 4D0
    Y(5) = -2D0 * A2 * dsin(2D0 * Theta - 2D0 * thetaB) / S ** 5D0 - 2D0 * S1 * Theta1 / S - IzA * thetaA_tid/(m*S**2D0) - IzB * thetaB_tid/(m*S**2D0)
    Y(6) = 2D0 * m * A2 * dsin(2D0 * Theta - 2D0 * thetaB) / (IzB * S ** 3D0)
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
    
    SUBROUTINE YHC4BT(t0,t1,X,Y)
    ! @sub YHC4 八维状态方程
    ! @pram: t0 初始时间
    ! @pram: t1 结束时间
    ! @pram: X  初始状态矢量   
    ! @pram: Y  输出状态矢量
    use global
    REAL(8) :: t0,t1,x(8),y(8),C,S,Theta,thetaA,thetaB,S1,Theta1,thetaA1,thetaB1
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
    Y(6) = -A2 * dsin(2D0 * Theta - 2D0 * thetaB) / S ** 5D0 - A3 * dsin(2D0 * Theta - 2D0 * thetaA) / S ** 5D0 - 2D0 * S1 * Theta1 / S
    Y(7) = m * A3 * dsin(2D0 * Theta - 2D0 * thetaA) / (IzA * S ** 3D0) - ThY * alphaA / (1 - mu)
    Y(8) = m * A2 * dsin(2D0 * Theta - 2D0 * thetaB) / (IzB * S ** 3D0)
    end
    
    SUBROUTINE YHC4New(t0,t1,X,Y)
    use global
    REAL(8) :: t0,t1,x(8),y(8),S,theta,phi,thetaB,S1,theta1,phi1,thetaB1,delta
    REAL(8) :: COST,SINT,COSD,SIND

    S = X(1)
    theta = X(2)
    phi = X(3)
    thetaB = X(4)
    delta = X(2) - X(3)
    S1 = X(5)
    theta1 = X(6)
    phi1 = X(7)
    thetaB1 = X(8)
    
    COST=DCOS(2D0 * X(2))
	SINT=DSIN(2D0 * X(2))
	COSD=DCOS(2D0 * delta)
	SIND=DSIN(2D0 * delta)
    
    Y(1) = S1
    Y(2) = theta1
    Y(3) = phi1
    Y(4) = thetaB1
    Y(5) = S * (theta1+thetaB1) ** 2D0 - 1D0 / S ** 2D0 - 3D0 * (A1 + A2 * COST + A3 * COSD) / S ** 4D0
    Y(6) =  - 2D0 * S1 * (theta1 + thetaB1) / S - 2D0 * A2 * SINT / S ** 5D0 - 2D0 * A3 * SIND / S ** 5D0 - 2D0 * m * A2 * SINT / (IzB * S ** 3D0)
    Y(7) = 2D0 * m * A3 * SIND / (IzA * S ** 3D0) - 2D0 * m * A2 * SINT / (IzB * S ** 3D0)
    Y(8) = 2D0 * m * A2 * SINT / (IzB * S ** 3D0)
    end
    
    SUBROUTINE YHC4_Mon(t0,t1,X,Y)
    ! sub YHC4_Mon Monodromy矩阵与思维运动方程的积分右函数
    ! pram: t0     积分初始时间
    ! pram: tf     积分结束时间
    ! pram: X      积分初始状态量和结束时输出状态量
    ! pram: Y      输出状态量 
    use global
    REAL(8) :: t0,t1,x(72),y(72),S,theta,phi,thetaB,S1,theta1,phi1,thetaB1,A(8,8),MM(8,8),AMM(8,8)

    S = X(65)
    theta = X(66)
    phi = X(67)
    thetaB = X(68)
    S1 = X(69)
    theta1 = X(70)
    phi1 = X(71)
    thetaB1 = X(72)
    
    DO I = 1,8
        DO J = 1,8
            MM(I,J) = X((I-1)*8+J)
        ENDDO
    ENDDO
    
    A(1,1:4) = 0D0
    A(1,5) = 1D0
    A(1,6:8) = 0D0
    A(2,1:5) = 0D0
    A(2,6) = 1D0
    A(2,7:8) = 0D0
    A(3,1:6) = 0D0
    A(3,7) = 1D0
    A(3,8) = 0D0
    A(4,1:7) = 0D0
    A(4,8) = 1D0
    A(5,1) = (theta1+thetaB1)**2D0 + 2D0/S**3D0 + 12D0/S**5D0 * (A1 + A2*dcos(2D0*theta) + A3*dcos(2D0*(theta-phi)))
    A(5,2) = 6D0/S**4D0*(A2*dsin(2*theta)+A3*dsin(2D0*(theta-phi)))
    A(5,3) = -6D0*A3*dsin(2D0*(theta-phi))
    A(5,4) = 0D0
    A(5,5) = 0D0
    A(5,6) = 2D0*S*(theta1+thetaB1)
    A(5,7) = 0D0
    A(5,8) = 2D0*S*(theta1+thetaB1)
    A(6,1) = 2D0*S1*(theta1+thetaB1)/S**2D0 + 10D0*(A2*dsin(2*theta)+A3*dsin(2D0*(theta-phi)))/S**6D0 + 6D0*m*A2*dsin(2D0*theta)/(IzB*S**4D0)
    A(6,2) = -4D0*(A2*dcos(2D0*theta)+A3*dcos(2D0*(theta-phi)))/S**5D0 - 4D0*m*A2*dcos(2D0*theta)/(IzB*S**3D0)
    A(6,3) = 4D0*A3*dcos(2D0*(theta-phi))/S**5D0
    A(6,4) = 0D0
    A(6,5) = -2D0*(theta1+thetaB1)/S
    A(6,6) = -2D0*S1/S
    A(6,7) = 0D0
    A(6,8) = -2D0*S1/S
    A(7,1) = 6D0*m*A2*dsin(2D0*theta)/(IzB*S**4D0) - 6D0*m*A2*dsin(2D0*(theta-phi))/(IzA*S**4D0)
    A(7,2) = 4D0*m*A3*dcos(2D0*(theta-phi))/(IzA*S**3D0) - 4D0*m*A2*dcos(2D0*theta)/(IzB*S**3D0)
    A(7,3) = -4D0*m*A3*dcos(2D0*(theta-phi))/(IzA*S**3D0)
    A(7,4:8) = 0D0
    A(8,1) = -6D0*m*A2*dsin(2D0*theta)/(IzB*S**4D0)
    A(8,2) = 4D0*m*A2*dcos(2D0*theta)/(IzB*S**3D0)
    A(8,3:8) = 0D0
    
    AMM = MATMUL(A,MM)
    
    DO I = 1,8
        DO J = 1,8
            Y((I-1)*8+J) = AMM(I,J)
        ENDDO
    ENDDO
    
    Y(65) = S1
    Y(66) = theta1
    Y(67) = phi1
    Y(68) = thetaB1
    Y(69) = S * (theta1+thetaB1) ** 2D0 - 1D0 / S ** 2D0 - 3D0 * (A1 + A2 * dcos(2D0 * theta) + A3 * dcos(2D0 * (theta - phi))) / S ** 4D0
    Y(70) = -2D0 * A2 * dsin(2D0 * theta) / S ** 5D0 - 2D0 * A3 * dsin(2D0 * (theta - phi)) / S ** 5D0 - 2D0 * S1 * (theta1 + thetaB1) / S - 2D0 * m * A2 * dsin(2D0 * theta) / (IzB * S ** 3D0)
    Y(71) = 2D0 * m * A3 * dsin(2D0 * (theta - phi)) / (IzA * S ** 3D0) - 2D0 * m * A2 * dsin(2D0 * theta) / (IzB * S ** 3D0)
    Y(72) = 2D0 * m * A2 * dsin(2D0 * theta) / (IzB * S ** 3D0)
    end