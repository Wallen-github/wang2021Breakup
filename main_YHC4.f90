   
PROGRAM main
    USE global
    USE FindRoots
    IMPLICIT REAL(8) (A-H,O-Z)
    external :: YHC4BT
    REAL(8) :: t0,tf,err,step,w0(2),hh,tspan,EE,x0(8),alphaS0,beta
    integer(4) :: time_begin, time_end
    
102 format(1X, E25.18,1X,E25.18)
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
104 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
106 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
107 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
    
    call system_clock(time_begin)
    cB = rB / (a_bB * b_cB ** 2D0) ** (1D0 / 3D0)
    bB = cB * b_cB
    aB = bB * a_bB
    cA = rA / (a_bA * b_cA ** 2D0) ** (1D0 / 3D0)
    bA = cA * b_cA
    aA = bA * a_bA
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
    
    thetaA1_Const = 5.123D0
    S0 = R / Length
    theta = 0D0
    K = dsqrt((S0 ** 2D0 + 3D0 * (A1 + A2 * dcos(2D0*theta))) * (IzB + m * S0 ** 2D0) ** 2D0 / S0 ** 5D0)
    
    write(*,*) '==========================Initial Value==========================='
    write(*,*) 'Tunit = ', Tunit, 'year'
    write(*,*) 'Tunit = ', Tunit*365D0, 'day'
    write(*,*) 'Tunit = ', Tunit*365D0*24D0, 'hour'
    write(*,*) 'Tunit = ', Tunit*365D0*24D0*6D1, 'min'
    write(*,*) 'S0 = ',S0,'K = ',K
    
    alphaS0 = 1D-1
    beta = 5D0*PI/18D1
    write(*,*) '========================Test 1==========================='
    write(*,*) 'alpha,beta = '
    write(*,102) alphaS0*S0,beta
    CALL GetAB4(alphaS0*S0,beta,0D0,0D0,x0,w0)
    write(*,*) 'x0 ='
    write(*,104) x0
    
    !积分得到x0
    tspan = 5D4
    t0=0D0
    tf=t0+tspan
    err=1D-10
    hh = 1D-2
    Num = 2D4
    write(*,*) '========================Propagate Loading==========================='
    open(12,file='G:\Master\科研工作\Binary system\2020.7.9Paperwork\Pic\RotationOfA\ORBIT2.DAT', status='REPLACE')
    CALL Propagate4BT(12,YHC4BT,t0,tf,tspan,err,hh,Num,x0)
    close(12)
    write(*,*) '===============================End=================================='
    write(*,*) 'xt ='
    write(*,104) x0
    
    call system_clock(time_end)
    PRINT *, 'Time of operation was ', (time_end - time_begin)/1D3, ' seconds'
END 