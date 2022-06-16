PROGRAM main
    USE global
    USE FindRoots
    IMPLICIT REAL(8) (A-H,O-Z)
    external :: YHC2A
    character(len=80) :: filename
    REAL(8) :: x0(6),x4(4),state(9),t0,tf,err,hh1,x1(6),w1(2),w0(2),xi0,eta0
    REAL(8) :: root0(1),roott(1)
    REAL(8) :: alpha,beta,phi1,phi2,xree(4),tree,xmid(4),xx(4)
    integer(4) :: time_begin, time_end,unitT,FILEID
    
100 format(1X, E25.18,1X,E25.18,1X,E25.18)
101 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
102 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
103 format(1X, E25.18,1X,E25.18)
104 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
    
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
    xi0 = 0.0D0
    eta0 = 0.0D0
    CALL Getx03(2, xi0, eta0,x0,w0)
    write(*,*) 'x0 = '
    write(*,100) x0
    write(*,*) 'w0 = '
    write(*,101) w0
    
    filename = "G:\Master\科研工作\Binary system\2020.7.9Paperwork\Pic\TestMigration\JointOut1.dat"
    open(14,file='G:\Master\科研工作\Binary system\2020.7.9Paperwork\Pic\TestMigration\DotsV1.txt',status='REPLACE')
    open(15,file='G:\Master\科研工作\Binary system\2020.7.9Paperwork\Pic\TestMigration\DotsP1.txt',status='REPLACE')
    Tmin = 1
    Tmax = 180000
    Num = 1000
    h = (Tmax-Tmin)/Num
    DO I = 1,Num
    
        !读取状态量
        unitT = Tmin + h*I !外迁演化参数
        FILEID = 12
        CALL ReadData(FILEID,filename,unitT,state)
        write(*,*) 'state(x0,energy,Angular,t0*Tunit) = '
        write(*,*) state
    
        !将6维的x0转化为4维的x0
        x0 = state(1:6)
        x4(1) = x0(1)
        x4(2) = x0(2) - x0(3)
        x4(3) = x0(4)
        x4(4) = x0(5) - x0(6)
        write(*,*) 'x4 = '
        write(*,103) x4
    
        !求解新的系统参数
        err = 1D-10
        hh1 = 1D-2
        K = state(8)
        root0 = (/x4(1)/)
        CALL NewRaf(equ,1,root0,err,roott)
        S0 = roott(1)
        write(*,*) 'S0 = ',S0,'K = ',K
    
        ! 求解新系统的长短周期频率
        CALL Getx02(2, xi0, eta0, x1, w1)
    
        !求解dots
        t0=0D0
        tf = 1.5D2*(w0(2)/w1(2))
        err = 1D-10
        DD = 1D-5
        hh1 = 1D-2
        write(*,*)  '=========================Find dots============================'
        write(*,*) 'tf = ',tf
        CALL FindDot(YHC2A,t0,tf,err,DD,hh1,x4,xree)
        write(*,*) 'xree = '
        write(*,103) xree
    
    ENDDO
    close(14)
    close(15)
    
    call system_clock(time_end)
    PRINT *, 'Time of operation was ', (time_end - time_begin)/1D3, ' seconds'
    end