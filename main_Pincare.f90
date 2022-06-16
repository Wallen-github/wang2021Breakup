PROGRAM main
    USE global
    USE FindRoots
    IMPLICIT REAL(8) (A-H,O-Z)
    external :: YHC2A,YHC3ABT
    character(len=200) :: filename
    REAL(8) :: t0,tf,err,DD,hh1
    REAL(8) :: x0(4),V,Energy,DSmax,DSmin
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
    
    ! Get initial values
    CALL GetAB(0D-2,0D0*PI/18D1,0D0,0D0,x0)
    thetaB1 = (K-m*x0(1)**2D0*x0(4))/(m*x0(1)**2D0+IzB)
    V = m/x0(1) + m/(2D0*x0(1)**3D0)*(A1+A2*dcos(2D0*theta))
    Energy = m/2D0*(x0(3) ** 2D0 + x0(1) ** 2D0 * (x0(4)+thetaB1) ** 2D0) + 1D0/2D0 * IzB * thetaB1 ** 2D0 - V
    write(*,*) 'x0 = '
    write(*,*) x0
    write(*,*) 'S0 = ',S0
    write(*,*) -dsqrt((2D0*Energy*m*x0(1)**2D0 + 2D0*V*m*x0(1)**2D0 + 2D0*Energy*IzB + 2D0*IzB*V - K**2D0) / (IzB*m*x0(1)**2D0))
    
    ! Set Parameters
    write(*,*)  '===========================Set Parameters=============================='
    open(15,file='G:\Master\科研工作\Binary system\2020.7.9Paperwork\Pic\Migration4\Dot.dat',status='REPLACE')
    open(16,file='G:\Master\科研工作\Binary system\2020.7.9Paperwork\Pic\Migration4\Pin.dat',status='REPLACE')
    open(12,file='G:\Master\科研工作\Binary system\2020.7.9Paperwork\Pic\Migration4\orbit1.dat', status='REPLACE')
    err = 1D-14
    DD = 1D-5
    hh1 = 1D-15
    t0=0D0
    tf = 1.5D4
    
    DSmax = 0.5D-1
    DSmin = 0D-1
    Num = 100
    !Energy = -8D-3
    write(*,*) 'Integrated time : ',tf
    write(*,*) 'Energy parameter : ',Energy
    
    DO I = 1,Num
        x0(1) = S0 + DSmin + (I)*DSmax/Num
        x0(2) = 0D0
        x0(3) = 0D0
        V = m/x0(1) + m/(2D0*x0(1)**3D0)*(A1+A2*dcos(2D0*theta))
        x0(4) = -dsqrt((2D0*Energy*m*x0(1)**2D0 + 2D0*V*m*x0(1)**2D0 + 2D0*Energy*IzB + 2D0*IzB*V - K**2D0) / (IzB*m*x0(1)**2D0))
        write(*,*) x0(1),x0(4)
        IF (x0(4)<0) THEN
            CALL FindDot(YHC2A,t0,tf,err,DD,hh1,x0)
        ELSE
            GO TO 801 
        ENDIF
        
801 ENDDO
    write(*,*)  '===========================End=============================='
    close(15)
    close(12)
    close(16)
    
    call system_clock(time_end)
    PRINT *, 'Time of operation was ', (time_end - time_begin)/1D3, ' seconds'
    end