PROGRAM main
    USE global
    USE FindRoots
    IMPLICIT REAL(8) (A-H,O-Z)
    external :: YHC2A
    character(len=80) :: filename
    real(8) :: x0(4),w0(2),alpha,beta,x6(6),phi1,phi2
    integer(4) :: time_begin, time_end,unitT,FILEID
    
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
    A3 = 3D0 * alphaA ** 2D0 * J22A
    IzB = mu * (aB ** 2D0 + bB ** 2D0) / (5D0 * Length ** 2D0)
    IzA = (1 - mu) * (aA ** 2D0 + bA ** 2D0) / (5D0 * Length ** 2D0)
    Izt = IzB+m*S0**2D0
    Tunit = dsqrt(Length ** 3D0 / (MA + MB) / G) / SecOfYear
    
    S0 = R / Length
    theta = 0D0
    K = dsqrt((S0 ** 2D0 + 3D0 * (A1 + A2 * dcos(2D0*theta))) * (IzB + m * S0 ** 2D0) ** 2D0 / S0 ** 5D0)
    
    write(*,*) '========================Test 1==========================='
    alpha = 12D0
    beta = 2D1*PI/18D1
    phi1 = 1D-2
    phi2 = 2D-2
    write(*,*) 'alpha,beta,phi1,phi2 = '
    write(*,104) alpha,beta,phi1,phi2
    CALL GetAB(alpha,beta,phi1,phi2,x0)
    write(*,*) 'x0 = '
    write(*,104) x0
    CALL GetAB_re(x0,alpha,beta,phi1,phi2)
    !CALL GetAB_re2(x0,alpha,beta)
    write(*,*) 'alpha,beta,phi1,phi2 = '
    write(*,104) alpha,beta,phi1,phi2
    CALL GetAB(alpha,beta,phi1,phi2,x0)
    write(*,*) 'x0 = '
    write(*,104) x0
    
    write(*,*) '========================Test 2==========================='
    x0 = (/1D1/3D0,0D0,0D0,-3D-1/)
    write(*,*) 'x0 = '
    write(*,104) x0
    CALL GetAB_re(x0,alpha,beta,phi1,phi2)
    write(*,*) 'alpha,beta,phi1,phi2 = '
    write(*,104) alpha,beta,phi1,phi2
    CALL GetAB(alpha,beta,phi1,phi2,x0)
    write(*,*) 'x0 = '
    write(*,104) x0
    CALL GetAB_re(x0,alpha,beta,phi1,phi2)
    write(*,*) 'alpha,beta,phi1,phi2 = '
    write(*,104) alpha,beta,phi1,phi2
    
    write(*,*) '========================Test 3==========================='
    write(*,*) 'alpha,beta = '
    write(*,102) alpha,beta
    CALL GetAB3(alpha,beta,phi1,phi2,x6,w0)
    write(*,*) 'x6 = '
    write(*,103) x6
    !将6维的x0转化为4维的x0
    x0(1) = x6(1)
    x0(2) = x6(2) - x6(3)
    x0(3) = x6(4)
    x0(4) = x6(5) - x6(6)
    write(*,*) 'x0 = '
    write(*,104) x0
    CALL GetAB_re2(x0,alpha,beta)
    write(*,*) 'alpha,beta = '
    write(*,102) alpha,beta
    
    
    call system_clock(time_end)
    PRINT *, 'Time of operation was ', (time_end - time_begin)/1D3, ' seconds'
    end