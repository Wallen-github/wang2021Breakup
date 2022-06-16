   
PROGRAM main
    USE global
    USE FindRoots
    IMPLICIT REAL(8) (A-H,O-Z)
    external :: YHC3ABT,YHC3AB
    REAL(8) :: thetaB_tid,thetaA_tid,ConstB,ThY
    REAL(8) :: S,Theta1,thetaA1,thetaB1
    REAL(8) :: Se,Semax,Semin,h
    integer(4) :: Num
    integer(4) :: time_begin, time_end
    
102 format(1X, E25.18,1X,E25.18)    
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
104 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
105 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
106 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
109 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18)
111 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18)
113 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X, E25.18,1X,E25.18)
    
    call system_clock(time_begin)
    cB = rB / (a_bB * b_cB ** 2D0) ** (1D0 / 3D0)
    bB = cB * b_cB
    aB = bB * a_bB
    cA = rA / (a_bA * b_cA ** 2D0) ** (1D0 / 3D0)
    bA = cA * b_cA
    aA = bA * a_bA
    Length = (aA + aB)
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
    A3 = 3D0 * alphaB ** 2D0 * J22A
    IzB = mu * (aB ** 2D0 + bB ** 2D0) / (5D0 * Length ** 2D0)
    IzA = (1 - mu) * (aA ** 2D0 + bA ** 2D0) / (5D0 * Length ** 2D0)
    Izt = IzB+m*S0**2D0
    Tunit = dsqrt(Length ** 3D0 / (MA + MB) / G) / SecOfYear
    
    thetaA1_Const = 0D0
    S0 = R / Length
    theta = 0D0
    K = IzA*thetaA1_Const + dsqrt((S0 ** 2D0 + 3D0 * (A1 + A2 * dcos(2D0*theta))) * (IzB + m * S0 ** 2D0) ** 2D0 / S0 ** 5D0)
    
    write(*,*) '==========================Initial Value==========================='
    write(*,*) 'Tunit = ', Tunit, 'year'
    write(*,*) 'Tunit = ', Tunit*365D0, 'day'
    write(*,*) 'Tunit = ', Tunit*365D0*24D0, 'hour'
    write(*,*) 'Tunit = ', Tunit*365D0*24D0*6D1, 'min'
    CALL Byorp(ConstB,ThY)
    write(*,*) 'ConstB = ',ConstB,'ThY = ',ThY
    
    
    write(*,*) '========================Propagate Loading==========================='
    open(12,file='G:\Master\科研工作\Binary system\2020.7.9Paperwork\Pic\Migration\BT.txt', status='REPLACE')
    Theta1 = 1D0
    thetaA1 = thetaA1_COnst
    thetaB1 = Theta1
    Semax = 5D0
    Semin = 1D0
    Num = 100
    h = (Semax-Semin)/Num
    DO I = 1,Num
        Se = Semin + (I-1)*h
        CALL Tide(Se,Theta1,thetaA1,thetaB1,thetaA_tid,thetaB_tid)
        write(12,105) Se, thetaA_tid, thetaB_tid, IzA * thetaA_tid/(m*Se**2D0) + IzB * thetaB_tid/(m*Se**2D0), ConstB / (Se * m)
    enddo
    
    close(12)
    write(*,*) '===============================End=================================='

    
    call system_clock(time_end)
    PRINT *, 'Time of operation was ', (time_end - time_begin)/1D3, ' seconds'
END 