!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! @file: ES_Method.f90
! @brief: 此module里包括多个常用的方法
! @author: Wang Hai-Shuo
! @time: 2020.3.3
    
    
    
    SUBROUTINE Tide(S,Theta1,thetaA1,thetaB1,thetaA_tid,thetaB_tid)
    use global
    REAL(8) :: S,Theta1,thetaA1,thetaB1
    REAL(8) :: rho,QloveB,deltaB,thetaB_tid,thetaA_tid,QloveA,deltaA
    
    rho = p*Length**3D0/(MA + MB)
    QloveB = 6D5*(rB/1D3)
    QloveA = 6D5*(rA/1D3)
    deltaB = dsqrt(6D0/(PI*rho))*3D0/2D0/QloveB*(1D0-mu)**2D0/S**6D0*(aB/length)**5D0/IzB
    deltaA = dsqrt(6D0/(PI*rho))*3D0/2D0/QloveA*mu**2D0/S**6D0*(aA/length)**5D0/IzA

    if (abs(Theta1-thetaB1)>deltaB) then
        thetaB_tid = SIGN (1D0, Theta1-thetaB1) * 3D0 * (1D0-mu)**2D0 / (2D0 * QloveB * S**6D0) * (aB/length) ** 5D0 / IzB
    else
        thetaB_tid = (Theta1-thetaB1)/deltaB * 3D0 * (1D0-mu)**2D0 / (2D0 * QloveB * S**6D0) * (aB/length) ** 5D0 / IzB
    end if
    if (abs(Theta1-thetaA1)>deltaA) then
        thetaA_tid = SIGN (1D0, Theta1-thetaA1) * 3D0 * mu**2D0 / (2D0 * QloveA * S**6D0) * (aA/length) ** 5D0 / IzA
    else
        thetaA_tid = (Theta1-thetaA1)/deltaA * 3D0 * mu**2D0 / (2D0 * QloveA * S**6D0) * (aA/length) ** 5D0 / IzA
    end if
    
    end
    
    SUBROUTINE Byorp(ConstB,ThY)
    use global
    REAL(8) :: fBY = 1D-2,CY = 0.25D-1,Bs = 2D0 / 3D0,Fs_unit = 1.0D17, es = 0D0, Fs,ThY,ConstB
    
    Fs = Fs_unit * Length ** 2D0 / (G * (MA + MB) ** 2D0)
    ConstB = Bs * PI * rB ** 2D0 * fBY * Fs / (au ** 2D0 * dsqrt(1 - es))
    ThY = Bs*Fs*CY/((au/Length)**2D0*dsqrt(1D0-es**2D0))
    
    end
    
    