PROGRAM main
    USE global
    USE FindRoots
    IMPLICIT REAL(8) (A-H,O-Z)
    external :: YHC2A,YHC3ABT
    character(len=200) :: filename
    REAL(8) :: x0(6),x4(4),state(12),t0,tf,err,hh1,x1(6),w1(2),w0(2),alpha(2),alphaS0,beta
    REAL(8) :: hh,tspan,EE
    integer(4) :: Num
    REAL(8) :: alS_max,alS_min,beD_max,beD_min,DD
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
    
    
    !给出初值
    alphaS0 = 0.1D0
    beta = 5D0
    CALL GetAB3(alphaS0*S0,beta*PI/18D1,0D0,0D0,x0,w0)
    write(*,*) 'x0 = '
    write(*,103) x0
    write(*,*) 'w0 = '
    write(*,102) w0
    
    
    !将6维的x0转化为4维的x0
    x4(1) = x0(1)
    x4(2) = x0(2) - x0(3)
    x4(3) = x0(4)
    x4(4) = x0(5) - x0(6)
    write(*,*) 'x4 = '
    write(*,104) x4
    write(*,*) 'S0 = '
    write(*,*) S0
    
    !求解dots
    t0=0D0
    tf = 1.5D3
    err = 1D-14
    DD = 1D-5
    hh1 = 1D-15
    write(*,*)  '===========================Find dots=============================='
    write(*,*) 'tf = ',tf
    open(15,file='G:\Master\科研工作\Binary system\2020.7.9Paperwork\Pic\Test\Dot.dat',status='REPLACE')
    open(16,file='G:\Master\科研工作\Binary system\2020.7.9Paperwork\Pic\Test\Pin.dat',status='REPLACE')
    open(12,file='G:\Master\科研工作\Binary system\2020.7.9Paperwork\Pic\Test\orbit1.dat', status='REPLACE')
    write(15,104) alphaS0,beta,0D0,0D0
    CALL FindDot(YHC2A,t0,tf,err,DD,hh1,x4)
    close(15)
    close(12)
    close(16)
    
    ! 积分轨道
    tspan = tf
    t0=0D0
    tf=t0+tspan
    err=1D-10
    hh = 1D-2
    Num = 1
    write(*,*) '========================Propagate Loading==========================='
    open(12,file='G:\Master\科研工作\Binary system\2020.7.9Paperwork\Pic\Test\orbit.dat', status='REPLACE')
    CALL Propagate2(12,YHC2A,t0,tf,tspan,err,hh,Num,x4)
    close(12)
    write(*,*) '===============================End=================================='
    write(*,*) 'xt ='
    write(*,104) x4

    
    !循环contour图
    alS_max = 30D-2
    alS_min = 0D0
    beD_max = 8D1
    beD_min = 0D0
    Num = 100
    write(*,*) '=========================Contour Loading============================'
    tf = 1.5D2
    write(*,*) 'tf = ',tf
    open(15,file='G:\Master\科研工作\Binary system\2020.7.9Paperwork\Pic\Test\Sta.dat',status='REPLACE')
    open(16,file='G:\Master\科研工作\Binary system\2020.7.9Paperwork\Pic\Test\AB.dat',status='REPLACE')
    CALL Contour(15,alS_max,alS_min,beD_max,beD_min,Num,tf)
    close(15)
    close(16)
    write(*,*) 'T = ',state(9),'year'
    write(*,*) '===============================End=================================='
    
    call system_clock(time_end)
    PRINT *, 'Time of operation was ', (time_end - time_begin)/1D3, ' seconds'
    end