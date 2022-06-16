PROGRAM main
    USE global
    USE FindRoots
    IMPLICIT REAL(8) (A-H,O-Z)
    !external :: YHCB,YHC2,YHC3,YHCT 
    integer(4) :: time_begin, time_end
    real(8) :: w0(2),theta,abBmin,abBmax,h,alpha(2),abAmin,abAmax,S0min,S0max
    integer(4) :: Num

102 format(1X, E25.18,1X,E25.18)
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
104 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
106 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
107 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
    
    call system_clock(time_begin)
    
    !S0min = 1D0
    !S0max = 1D1
    !Num = 100
    !h = (S0max - S0min) / Num
    
    !rBmin = 0D0
    !rBmax = rA
    !Num = 100
    !h = (rBmax - rBmin) / Num
    
    abBmin = 1D0
    abBmax = 3D0
    Num = 100
    h = (abBmax - abBmin) / Num
    
    open(12,file='G:\Master\科研工作\Binary system\2020.7.9Paperwork\Pic\EP Of Averaged System\W0vsPar.txt',status='REPLACE')
    do I = 1,Num
        write(*,*) 'I = ',I
        a_bB = abBmin + (I - 1) * h
        rB = 8D2
        !rB = rBmin + (I-1) * h
        !a_bA = abAmin + (I - 1) * h
        !S0 = S0min + (I -1) * h
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
        J2A = 0D-1!(aA ** 2D0 + bA ** 2D0 - 2D0 * cA ** 2D0) / (10D0 * rA ** 2D0)
        J22B = (aB ** 2D0 - bB ** 2D0) / (20D0 * rB ** 2D0)
        J22A = (aA ** 2D0 - bA ** 2D0) / (20D0 * rA ** 2D0)
        A1 = (alphaB ** 2D0 * J2B + alphaA ** 2D0 * J2A)/2D0
        A2 = 3D0 * alphaB ** 2D0 * J22B
        A3 = 3D0 * alphaB ** 2D0 * J22A
        IzB = mu * (aB ** 2D0 + bB ** 2D0) / (5D0 * Length ** 2D0)
        IzA = (1 - mu) * (aA ** 2D0 + bA ** 2D0) / (5D0 * Length ** 2D0)
        Izt = IzB+m*S0**2D0
        Tunit = dsqrt(Length ** 3D0 / (MA + MB) / G) / SecOfYear
    
        S0 = R / Length
        theta = 0D0
        K = dsqrt((S0 ** 2D0 + 3D0 * (A1 + A2 * dcos(2D0*theta))) * (IzB + m * S0 ** 2D0) ** 2D0 / S0 ** 5D0)
        
        CALL Getw0(w0, alpha)
        write(*,*) 'w0 = '
        write(*,102) w0
        write(12,104) a_bB, w0(1), w0(2),mu
    enddo
    close(12)
    
    call system_clock(time_end)
    PRINT *, 'Time of operation was ', (time_end - time_begin)/1D3, ' seconds'
    end