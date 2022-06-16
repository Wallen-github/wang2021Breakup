PROGRAM main
    USE global
    USE FindRoots
    IMPLICIT REAL(8) (A-H,O-Z)
    character(len=80) :: filename
    REAL(8) :: x0(6),x4(4),state(9),t0,tf,err,hh1,x1(6),w1(2),w0(2),xi0,eta0
    REAL(8) :: root0(1),roott(1)
    REAL(8) :: alpha,beta,phi1,phi2,xree(4),tree,xmid(4),xx(4)
    integer(4) :: time_begin, time_end,unitT,FILEID
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
    Length = (rA + aB)
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
    
    S0 = 1.5D0!R / Length
    theta = 0D0
    K = dsqrt((S0 ** 2D0 + 3D0 * (A1 + A2 * dcos(2D0*theta))) * (IzB + m * S0 ** 2D0) ** 2D0 / S0 ** 5D0)
    
    
    !循环contour图
    alS_max = 30D-2
    alS_min = 0D0
    beD_max = 8D1
    beD_min = 0D0
    Num = 500
    write(*,*) '=========================Contour Loading============================'
    
    tf = 3D2
    write(*,*) 'tf = ',tf
    write(*,*) '==============================='
    open(15,file='G:\Master\科研工作\Binary system\2020.7.9Paperwork\Pic\Test\Sta.dat',status='REPLACE')
    open(16,file='G:\Master\科研工作\Binary system\2020.7.9Paperwork\Pic\Test\AB.dat',status='REPLACE')
    CALL Contour_thetamax(15,alS_max,alS_min,beD_max,beD_min,Num,tf)
    close(15)
    close(16)
    write(*,*) '===============================End=================================='
    
    call system_clock(time_end)
    PRINT *, 'Time of operation was ', (time_end - time_begin)/1D3, ' seconds'
    end
    