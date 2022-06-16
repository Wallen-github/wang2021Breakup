!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! @file: ES_Method2.f90
! @brief: 此module里包括多个常用的方法
! @version 2: 重新编辑与梳理了相关的符合和变量表达
! @author: Wang Hai-Shuo
! @time: 2020.3.3
    
    
    
    SUBROUTINE Getw0(w0, alpha)
    ! @sub Getw0        求解平均化二维运动方程的长短周期频率（初始值）
    ! @output: w0       输出长短周期轨道频率
    ! @output: alpha    线性化模型中的系数参数
    ! @note:            其中代入的平衡点状态量（S0,0,0,0）
    use global
    REAL(8) :: a42,a43,a31,a34,w0(2),s1,s2,alpha(2)
    
    a31 = (K) ** 2D0 / (S0 ** 2D0 * m + IzB) ** 2D0 - 4D0 * S0 ** 2D0 * (K) ** 2D0 * m / (S0 ** 2D0 * m + IzB) ** 3D0 + 2D0 / S0 ** 3D0 + (12D0 * (A1 + A2)) / S0 ** 5D0
    a34 = 2D0 * S0 * (K) * IzB / (S0 ** 2D0 * m + IzB) ** 2D0
    a42 = -2D0 * A2 * (S0 ** 2D0 * m * 2D0 + IzB) / (IzB * S0 ** 5D0)
    a43 = -2D0 * (K) / (S0 * (S0 ** 2D0 * m + IzB))
    
    s1 = (a31 + a42 + a34 * a43) / 2D0 + dsqrt((a31 + a42 + a34 * a43) ** 2D0 - 4D0 * a31 * a42) / 2D0
    s2 = (a31 + a42 + a34 * a43) / 2D0 - dsqrt((a31 + a42 + a34 * a43) ** 2D0 - 4D0 * a31 * a42) / 2D0
    w0(1) = dsqrt(-s1) ! real(aimag(sqrt(s1)),kind=8)
    w0(2) = dsqrt(-s2) ! real(aimag(sqrt(s2)),kind=8)
    alpha(1) = -a43 * w0(1) / (a42 + w0(1) ** 2D0)
    alpha(2) = -a43 * w0(2) / (a42 + w0(2) ** 2D0)
    return
    end
    
    SUBROUTINE Getx02(index,xi0,eta0,x0,w0)
    ! @sub Getx02   求解两自由度的初始状态量
    ! @input: index  长短周期轨道指示器 
    ! @input: xi0    初始横向振幅
    ! @input: eta0   初始纵向振幅
    ! @output: x0     输出状态量
    ! @output: w0     输出长短周期轨道频率
    use global
    integer(4) :: index
    REAL(8) :: w0(2),x0(4),xi0,eta0,x1,x2,x3,x4,alpha(2)
    
    CALL Getw0(w0,alpha)
    x1 = S0 + xi0
    x2 = eta0
    x3 = eta0 * w0(index) / alpha(index)
    x4 = -xi0 * alpha(index) * w0(index)
    x0 = (/x1, x2, x3, x4/)
    return
    end
    
    SUBROUTINE Getx03(index,xi0,eta0,x0,w0)
    use global
    integer(4) :: index
    REAL(8) :: w0(2),x0(6),xi0,eta0,x1,x2,x3,x4,alpha(2)
    
    CALL Getw0(w0,alpha)
    x1 = S0 + xi0
    x2 = eta0
    x3 = eta0 * w0(index) / alpha(index)
    x4 = -xi0 * alpha(index) * w0(index)
    x0 = (/x1, x2, 0D0, x3, x4 + (K - m * x1 ** 2D0 * x4) / (m * x1 ** 2D0 + IzB), (K - m * x1 ** 2D0 * x4) / (m * x1 ** 2D0 + IzB)/)
    return
    end
    
    SUBROUTINE GetAB(al,be,phi1,phi2,x0)
    ! @sub GetAB    通过线性化参数求解两个自由度的状态量
    ! @pram: alpha  线性化参数，横向最大振幅
    ! @pram: beta   线性化参数，自由摆动角最大振幅
    ! @pram: phi1   线性化参数
    ! @pram: phi2   线性化参数
    ! @pram: x0     输出状态量
    use global
    REAL(8) :: w0(2),x0(4),xi,eta,xi1,eta1,al,be,phi1,phi2,alpha(2)
    
    CALL Getw0(w0,alpha)
    xi = al*dcos(phi2)+be*dcos(phi1)/alpha(1)
    eta = -alpha(2)*al*sin(phi2)-be*dsin(phi1)
    xi1 = -w0(2)*al*dsin(phi2)-be*w0(1)*dsin(phi1)/alpha(1)
    eta1 = -w0(2)*alpha(2)*al*dcos(phi2)-w0(1)*be*dcos(phi1)
    x0 = (/xi + S0, eta, xi1, eta1/)
    return
    end
    
    SUBROUTINE GetAB3(alpha0,beta0,phi1,phi2,x0,w0)
    ! @sub GetAB3    通过线性化参数求解三个自由度的状态量
    ! @pram: alpha  线性化参数，横向最大振幅
    ! @pram: beta   线性化参数，自由摆动角最大振幅
    ! @pram: phi1   线性化参数
    ! @pram: phi2   线性化参数
    ! @pram: x0     输出状态量
    use global
    REAL(8) :: w0(2),x0(6),xi0,eta0,alpha0,beta0,phi1,phi2,alpha(2),x1,x2,x3,x4
    
    CALL Getw0(w0,alpha)
    x1 = S0 + alpha0*dcos(phi2)+beta0*dcos(phi1)/alpha(1)
    x2 = -alpha(2)*alpha0*sin(phi2)-beta0*dsin(phi1)
    x3 = -w0(2)*alpha0*dsin(phi2)-beta0*w0(1)*dsin(phi1)/alpha(1)
    x4 = -w0(2)*alpha(2)*alpha0*dcos(phi2)-w0(1)*beta0*dcos(phi1)
    x0 = (/x1, x2, 0D0, x3, x4 + (K - m * x1 ** 2D0 * x4) / (m * x1 ** 2D0 + IzB), (K - m * x1 ** 2D0 * x4) / (m * x1 ** 2D0 + IzB)/)
    return
    end
    
    SUBROUTINE GetAB4(alpha0,beta0,phi1,phi2,x0,w0)
    ! @sub GetAB4    通过线性化参数求解四个自由度的状态量
    ! @Glob: thetaA1_Const     均匀自转的主星频率
    ! @pram: alpha  线性化参数，横向最大振幅
    ! @pram: beta   线性化参数，自由摆动角最大振幅
    ! @pram: phi1   线性化参数
    ! @pram: phi2   线性化参数
    ! @pram: x0     输出状态量
    use global
    REAL(8) :: w0(2),x0(8),xi0,eta0,alpha0,beta0,phi1,phi2,alpha(2),x1,x2,x3,x4
    
    CALL Getw0(w0,alpha)
    x1 = S0 + alpha0*dcos(phi2)+beta0*dcos(phi1)/alpha(1)
    x2 = -alpha(2)*alpha0*sin(phi2)-beta0*dsin(phi1)
    x3 = -w0(2)*alpha0*dsin(phi2)-beta0*w0(1)*dsin(phi1)/alpha(1)
    x4 = -w0(2)*alpha(2)*alpha0*dcos(phi2)-w0(1)*beta0*dcos(phi1)
    x0 = (/x1, x2, 0D0, 0D0, x3, x4 + (K - m * x1 ** 2D0 * x4) / (m * x1 ** 2D0 + IzB), thetaA1_Const ,(K - m * x1 ** 2D0 * x4) / (m * x1 ** 2D0 + IzB)/)
    return
    end
    
    SUBROUTINE Contour(FILEID,ecc_max,ecc_min,beD_max,beD_min,Num,tf)
    ! @sub Contour      计算稳定Contour图，由ecc，beta作为网格扫描参数，其中ecc=alpha/S0
    ! @input: FILEID    文件操作编号 
    ! @input: alS_max   alpha最大值
    ! @input: alS_min   alpha最小值
    ! @inout: beD_max   beta最大值
    ! @input: beD_min   beta最小值
    ! @input: Num       网格数目Num*Num
    ! @inout: tf        积分时间
    ! @output: sta      二维稳定图
    USE global
    external :: YHC2A
    REAL(8),allocatable :: Sta(:,:),xi(:),eta(:)
    REAL(8) :: h1,h2,ecc_max,ecc_min,beD_max,beD_min,x4(4),t0,tf,EE,hh,err,Energy,Angular
    integer(4) :: Num,FILEID
    allocate(xi(Num))
    allocate(eta(Num))
    allocate(Sta(Num,Num))
102 format(1X, E25.18,1X,E25.18)
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
104 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
106 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
107 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18) 
    err=1D-10
    hh = 1D-2
    t0 = 0D0
    MM = 4

    h1 = (ecc_max-ecc_min)/Num
    h2 = (beD_max-beD_min)/Num
    
    do I = 1,Num
        xi(I) = ecc_min + (I-1)*h1
        do J = 1,Num
            eta(J) = beD_min + (J-1)*h2
            CALL GetAB(xi(I)*S0,eta(J)*PI/1.8D2,0D0,0D0,x4)
            t0 = 0D0
            do while (t0<tf)
                Call RKF78(YHC2A,hh,t0,x4,EE,err,MM)
                if (ABS(x4(2))>PI) then
                    Sta(I,J) = 0D0
                    GO TO 20
                end if
            end do
            Sta(I,J) = 1D0
20          write(FILEID,'(1x,E25.18,$)') Sta(I,J)
        end do
        WRITE(FILEID,*) ' '
        write(16,102) xi(I),eta(I)
        write(*,'(f6.2,$)') 1D0*I/Num
    enddo
    write(*,*) 'Done'
    end
    
    SUBROUTINE Contour_thetamax(FILEID,ecc_max,ecc_min,beD_max,beD_min,Num,tf)
    ! @sub Contour      计算稳定Contour图，由ecc，beta作为网格扫描参数，其中ecc=alpha/S0
    ! @input: FILEID    文件操作编号 
    ! @input: alS_max   alpha最大值
    ! @input: alS_min   alpha最小值
    ! @inout: beD_max   beta最大值
    ! @input: beD_min   beta最小值
    ! @input: Num       网格数目Num*Num
    ! @inout: tf        积分时间
    ! @output: sta      二维稳定图
    USE global
    external :: YHC2A
    REAL(8),allocatable :: Sta(:,:),xi(:),eta(:)
    REAL(8) :: h1,h2,ecc_max,ecc_min,beD_max,beD_min,x4(4),t0,tf,EE,hh,err
    REAL(8) :: thetamax
    integer(4) :: Num,FILEID
    allocate(xi(Num))
    allocate(eta(Num))
    allocate(Sta(Num,Num))
102 format(1X, E25.18,1X,E25.18)
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
104 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
106 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
107 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18) 
    err=1D-10
    hh = 1D-2
    t0 = 0D0
    MM = 4

    h1 = (ecc_max-ecc_min)/Num
    h2 = (beD_max-beD_min)/Num
    
    do I = 1,Num
        xi(I) = ecc_min + (I-1)*h1
        do J = 1,Num
            eta(J) = beD_min + (J-1)*h2
            CALL GetAB(xi(I)*S0,eta(J)*PI/1.8D2,0D0,0D0,x4)
            t0 = 0D0
            thetamax = beD_min*PI/1.8D2
            do while (t0<tf)
                Call RKF78(YHC2A,hh,t0,x4,EE,err,MM)
                if (ABS(x4(2))>PI) then
                    Sta(I,J) = 9D1
                    GO TO 20
                end if
                if (thetamax<x4(2)) then
                    thetamax = x4(2)
                endif
            end do
            Sta(I,J) = thetamax*1.8D2/PI
20          write(FILEID,'(1x,E25.18,$)') Sta(I,J)
        end do
        WRITE(FILEID,*) ' '
        write(16,102) xi(I),eta(I)
        write(*,'(f6.2,$)') 1D0*I/Num
    enddo
    write(*,*) 'Done'
    end
    
    SUBROUTINE Contour_Energy(FILEID,ecc_max,ecc_min,beD_max,beD_min,Num,tf)
    ! @sub Contour      计算稳定Contour图，由ecc，beta作为网格扫描参数，其中ecc=alpha/S0
    ! @input: FILEID    文件操作编号 
    ! @input: alS_max   alpha最大值
    ! @input: alS_min   alpha最小值
    ! @inout: beD_max   beta最大值
    ! @input: beD_min   beta最小值
    ! @input: Num       网格数目Num*Num
    ! @inout: tf        积分时间
    ! @output: sta      二维稳定图
    USE global
    external :: YHC2A
    REAL(8),allocatable :: Sta(:,:),xi(:),eta(:)
    REAL(8) :: h1,h2,ecc_max,ecc_min,beD_max,beD_min,x4(4),t0,tf,EE,hh,err
    REAL(8) :: Energy,thetaB1
    integer(4) :: Num,FILEID
    allocate(xi(Num))
    allocate(eta(Num))
    allocate(Sta(Num,Num))
102 format(1X, E25.18,1X,E25.18)
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
104 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
106 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
107 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18) 
    err=1D-10
    hh = 1D-2
    t0 = 0D0
    MM = 4

    h1 = (ecc_max-ecc_min)/Num
    h2 = (beD_max-beD_min)/Num
    
    do I = 1,Num
        xi(I) = ecc_min + (I-1)*h1
        do J = 1,Num
            eta(J) = beD_min + (J-1)*h2
            CALL GetAB(xi(I)*S0,eta(J)*PI/1.8D2,0D0,0D0,x4)
            t0 = 0D0
            thetaB1 = (K - m*x4(1)**2D0*x4(4))/(m*x4(1)**2D0+IzB)
            Energy = m / 2D0 * (x4(3) ** 2D0 + x4(1) ** 2D0 * (x4(4)+thetaB1) ** 2D0) + 1D0 / 2D0 * IzB * thetaB1 ** 2D0 - m / x4(1) - m / (x4(1) ** 3D0) * (A1 + A2 * dcos(2D0 * x4(2)))
            do while (t0<tf)
                Call RKF78(YHC2A,hh,t0,x4,EE,err,MM)
                if (ABS(x4(2))>PI) then
                    Sta(I,J) = -26D-3
                    GO TO 20
                end if
            end do
            Sta(I,J) = Energy
20          write(FILEID,'(1x,E25.18,$)') Sta(I,J)
        end do
        WRITE(FILEID,*) ' '
        write(16,102) xi(I),eta(I)
        write(*,'(f6.2,$)') 1D0*I/Num
    enddo
    write(*,*) 'Done'
    end
    
    SUBROUTINE Contour_phi12(FILEID,ecc_max,ecc_min,beD_max,beD_min,phi1,phi2,Num,tf)
    ! @sub Contour      计算稳定Contour图，由ecc，beta作为网格扫描参数，其中ecc=alpha/S0
    ! @input: FILEID    文件操作编号 
    ! @input: alS_max   alpha最大值
    ! @input: alS_min   alpha最小值
    ! @inout: beD_max   beta最大值
    ! @input: beD_min   beta最小值
    ! @input: Num       网格数目Num*Num
    ! @inout: tf        积分时间
    ! @output: sta      二维稳定图
    USE global
    external :: YHC2A
    REAL(8),allocatable :: Sta(:,:),xi(:),eta(:)
    REAL(8) :: h1,h2,ecc_max,ecc_min,beD_max,beD_min,x4(4),t0,tf,EE,hh,err
    REAL(8) :: phi1,phi2
    integer(4) :: Num,FILEID
    allocate(xi(Num))
    allocate(eta(Num))
    allocate(Sta(Num,Num))
102 format(1X, E25.18,1X,E25.18)
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
104 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
106 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
107 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18) 
    err=1D-10
    hh = 1D-2
    t0 = 0D0
    MM = 4

    h1 = (ecc_max-ecc_min)/Num
    h2 = (beD_max-beD_min)/Num
    
    do I = 1,Num
        xi(I) = ecc_min + (I-1)*h1
        do J = 1,Num
            eta(J) = beD_min + (J-1)*h2
            CALL GetAB(xi(I)*S0,eta(J)*PI/1.8D2,phi1,phi2,x4)
            t0 = 0D0
            do while (t0<tf)
                Call RKF78(YHC2A,hh,t0,x4,EE,err,MM)
                if (ABS(x4(2))>PI) then
                    Sta(I,J) = 0D0
                    GO TO 20
                end if
            end do
            Sta(I,J) = 1D0
20          write(FILEID,'(1x,E25.18,$)') Sta(I,J)
        end do
        WRITE(FILEID,*) ' '
        write(16,102) xi(I),eta(I)
        write(*,'(f6.2,$)') 1D0*I/Num
    enddo
    write(*,*) 'Done'
    end
    
    SUBROUTINE Propagate3(FILEID,YHC,t0,tf,tspan,err,hh,Num,x0)
    ! @sub Propagate3 三个自由度运动方程的积分器，这一积分器与Propagate3Bt的不同之处在于，文件输出的格式，这里单次输出九个数据，三个一行，便于读取。
    ! @outfile FILEID： 输出数据文件FILEID，所有积分状态量，能量，角动量和以年为单位的积分时间x0,energy,Angular,t0*Tunit
    ! @pram: YHC    右函数 
    ! @pram: t0     积分初始时间
    ! @pram: tf     积分结束时间
    ! @pram: tspan  积分时间间隔（内循环，不输出数据），间隔tspan之后输出一次数据
    ! @pram: err    误差限
    ! @pram: hh     初始步长
    ! @pram: Num    积分时间间隔（外循环，输出数据），总计输出Num个数据
    ! @pram: x0     积分初始状态量和结束时输出状态量
    use global
    USE FindRoots
    external :: YHC
    REAL(8) :: t0,tf,err,step,Energy,Angular,w0(2),hh,tspan,x0(6),EE,the(10000000),max,placeholder
    REAL(8) :: S,Theta,thetaB,S1,Theta1,thetaB1
    integer :: FILEID
102 format(1X, E25.18,1X,E25.18)    
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
104 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
105 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
106 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
109 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18)
111 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18)
113 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X, E25.18,1X,E25.18)
    do I=1,Num
        J = 1
        do while (t0<tf)
            S = x0(1)
            Theta = x0(2)
            thetaB = x0(3)
            S1 = x0(4)
            Theta1 = x0(5)
            thetaB1 = x0(6)
            Energy = m / 2D0 * (S1 ** 2D0 + S ** 2D0 * Theta1 ** 2D0) + 1D0 / 2D0 * IzB * thetaB1 ** 2D0 + 1D0 / 2D0 * IzA * thetaA1 ** 2D0 - m / S - m / (S ** 3D0) * (A1 + A2 * dcos(2D0 * Theta - 2D0 * thetaB))
            Angular = m * S ** 2D0 * Theta1 + IzB * thetaB1 + IzA * thetaA1
            the(J) = Theta - thetaB
            J = J + 1
            Call RKF78(YHC,hh,t0,x0,EE,err,6)
            !write(FILEID,109) x0,energy,Angular,t0*Tunit 
            !WRITE(FILEID,103) x0
		    !WRITE(FILEID,103) energy,Angular,t0*Tunit
        enddo
        max = maxval(abs(the(1:J-1)))
        if (max>PI .OR. S<1D0) then
            write(*,*) 'breakup'
            GO TO 21
        end if
        write(*,*) 'I = ',I ,'/',Num, ';','T = ',t0*Tunit,'year',';','theta = ',max
        !write(FILEID,109) x0,energy,max,t0*Tunit 
		WRITE(FILEID,103) x0
		WRITE(FILEID,103) energy,Angular,t0*Tunit
        tf=t0+tspan
    enddo
21  write(*,*) 'DONE'

    end
    
    SUBROUTINE Propagate3BT(FILEID,YHC,t0,tf,tspan,err,hh,Num,x0)
    ! @sub Propagate3BT 三个自由度运动方程的积分器，这里的数据输出，单词输出一行是一个数据，无法进行格式化读取
    ! @outfile FILEID： 输出数据文件FILEID，所有积分状态量，能量，角动量和以年为单位的积分时间x0,energy,Angular,t0*Tunit
    ! @pram: YHC    右函数 
    ! @pram: t0     积分初始时间
    ! @pram: tf     积分结束时间
    ! @pram: tspan  积分时间间隔（内循环，不输出数据），间隔tspan之后输出一次数据
    ! @pram: err    误差限
    ! @pram: hh     初始步长
    ! @pram: Num    积分时间间隔（外循环，输出数据），总计输出Num个数据
    ! @pram: x0     积分初始状态量和结束时输出状态量
    use global
    USE FindRoots
    external :: YHC
    REAL(8) :: t0,tf,err,step,Energy,Angular,w0(2),hh,tspan,x0(6),EE,the(10000000),max,axis,ecc,axis_sum,ecc_sum,placeholder
    REAL(8) :: S,Theta,thetaB,S1,Theta1,thetaB1,thetaA1
    integer :: FILEID
    REAL(8) :: thetaB_tid,thetaA_tid,ConstB,ThY
    
102 format(1X, E25.18,1X,E25.18)    
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
104 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
105 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
106 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
109 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18)
110 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
111 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18)
113 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X, E25.18,1X,E25.18)
    do I=1,Num 
        J = 1
        axis_sum = 0
        ecc_sum = 0
        do while (t0<tf)
            S = x0(1)
            Theta = x0(2)
            thetaB = x0(3)
            S1 = x0(4)
            Theta1 = x0(5)
            thetaB1 = x0(6)
            thetaA1 = thetaA1_Const
            Energy = m / 2D0 * (S1 ** 2D0 + S ** 2D0 * Theta1 ** 2D0) + 1D0 / 2D0 * IzB * thetaB1 ** 2D0 - m / S - m / (S ** 3D0) * (A1 + A2 * dcos(2D0 * Theta - 2D0 * thetaB))
            Angular = m * S ** 2D0 * Theta1 + IzB * thetaB1 !+ IzA * thetaA1
            the(J) = Theta - thetaB
    
            !CALL Tide(S,Theta1,thetaA1,thetaB1,thetaA_tid,thetaB_tid)
            !CALL Byorp(ConstB,ThY)
            CALL mco_x2el (1D0,x0(1)*dcos(x0(2)),x0(1)*dsin(x0(2)),0D0,x0(4)*dcos(x0(2))-x0(1)*dsin(x0(2))*x0(5),x0(4)*dsin(x0(2))+x0(1)*dcos(x0(2))*x0(5),0D0,axis,ecc,placeholder,placeholder,placeholder,placeholder)
            axis_sum = axis_sum + axis
            ecc_sum = ecc_sum + ecc
            J = J + 1
            
            Call RKF78(YHC,hh,t0,x0,EE,err,6)
            !write(FILEID,111) x0, t0*Tunit, axis, ecc ,Energy,Angular
        enddo
        max = maxval(abs(the(1:J-1)))
        if (max>PI .OR. S<1D0) then
            write(*,*) 'breakup'
            GO TO 21
        end if
        write(*,*) 'I = ',I ,'/',Num, ';','T = ',t0*Tunit,'year',';','axis = ',axis
        !write(FILEID,111) x0, t0*Tunit, axis, ecc ,Energy,Angular 
        WRITE(FILEID,103) x0
		WRITE(FILEID,103) Energy,Angular,t0*Tunit
		WRITE(FILEID,103) axis,ecc,max
        tf=t0+tspan
    enddo
21  write(*,*) 'DONE'

    end
    
    SUBROUTINE Propagate4BT(FILEID,YHC,t0,tf,tspan,err,hh,Num,x0)
    use global
    USE FindRoots
    external :: YHC
    REAL(8) :: t0,tf,err,step,w0(2),hh,tspan,x0(8),EE,the(10000000),max,axis,ecc,axis_sum,ecc_sum,Energy,Angular,placeholder
    REAL(8) :: S,Theta,thetaA,thetaB,S1,Theta1,thetaA1,thetaB1
    integer :: FILEID
    REAL(8) :: thetaB_tid,thetaA_tid,ConstB,ThY
    
102 format(1X, E25.18,1X,E25.18)    
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
104 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
105 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
106 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
109 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18)
111 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18)
113 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X, E25.18,1X,E25.18)    
115 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X, E25.18,1X,E25.18)
    do I=1,Num 
        J = 1
        axis_sum = 0
        ecc_sum = 0
        do while (t0<tf)
            S = X0(1)
            Theta = X0(2)
            thetaA = X0(3)
            thetaB = X0(4)
            S1 = X0(5)
            Theta1 = X0(6)
            thetaA1 = X0(7)
            thetaB1 = X0(8)
            Energy = m / 2D0 * (S1 ** 2D0 + S ** 2D0 * Theta1 ** 2D0) + 1D0 / 2D0 * IzB * thetaB1 ** 2D0 + 1D0 / 2D0 * IzA * thetaA1 ** 2D0 - m / S - m / (2D0 * S ** 3D0) * (A1 + A2 * dcos(2D0 * Theta - 2D0 * thetaB) + A3 * dcos(2D0 * Theta - 2D0 * thetaA))
            Angular = m * S ** 2D0 * Theta1 + IzA * thetaA1 + IzB * thetaB1
            the(J) = Theta-thetaB

            !CALL Tide(S,Theta1,thetaA1,thetaB1,thetaA_tid,thetaB_tid)
            !CALL Byorp(ConstB,ThY)
            CALL mco_x2el (1D0,x0(1)*dcos(x0(2)),x0(1)*dsin(x0(2)),0D0,x0(5)*dcos(x0(2))-x0(1)*dsin(x0(2))*x0(6),x0(5)*dsin(x0(2))+x0(1)*dcos(x0(2))*x0(6),0D0,axis,ecc,placeholder,placeholder,placeholder,placeholder)
            axis_sum = axis_sum + axis
            ecc_sum = ecc_sum + ecc
            
            J = J+1
            Call RKF78(YHC,hh,t0,x0,EE,err,8)
            !write(FILEID,113) x0,t0*Tunit,axis,ecc,Energy,Angular
        enddo
        max = maxval(abs(the(1:J-1)))
        !if (max>PI .OR. S<(aA + aB)/R) then
        IF (S<(aA+aB)/R) THEN
            write(*,*) 'breakup'
            GO TO 21
        end if
        write(*,*) 'I = ',I ,'/',Num, ';','T = ',t0*Tunit,'year',';','omegaA = ',thetaA1
        write(FILEID,113) x0,t0*Tunit,axis_sum/(J-1),ecc_sum/(J-1),Energy,Angular 
        tf=t0+tspan
    enddo
21  write(*,*) 'DONE'

    end
    
    SUBROUTINE ReadData(FILEID,filename,unitT,x)
    ! @sub ReadData     读取数据文件，数据文件由Propagate3子函数给出，两个子函数必须配套使用
    ! @outfile FILEID： 输入数据文件FILEID
    ! @pram: filename   输入数据文件地址
    ! @pram: unitT      需要读取的第unitT行数据
    ! @pram: x          返回的状态量，九维
    implicit none
    integer, parameter :: IRECSZ = 20
    integer :: unitT,NRECL,fileid
    character(len=26)  :: hit(3)
    REAL(8) :: x(12)
    character(len=80) :: filename
    integer error
    
    NRECL = (unitT-1)*4+1
    open(unit=fileid, file=filename, access="direct",form="unformatted", recl=IRECSZ, status="old")
    read(fileid, rec=NRECL, IOSTAT=error) hit
    read(hit,"(1X,E25.18)")x(1:3)
    read(fileid, rec=NRECL+1, IOSTAT=error) hit
    read(hit,"(1X,E25.18)")x(4:6)
    read(fileid, rec=NRECL+2, IOSTAT=error) hit
    read(hit,"(1X,E25.18)")x(7:9)
    read(fileid, rec=NRECL+3, IOSTAT=error) hit
    read(hit,"(1X,E25.18)")x(10:12)
    close(fileid)

    end
    
    SUBROUTINE Propagate2(FILEID,YHC,t0,tf,tspan,err,hh,Num,x0)
    ! @sub Propagate2 两个个自由度运动方程的积分器，这一版本的该子函数将数据文件输出更改为了输出到单一文件，单次一行输出7个数据
    ! @outfile FILEID： 输出数据文件FILEID，所有积分状态量，能量，角动量和以年为单位的积分时间x0,energy,Angular,t0*Tunit
    ! @pram: YHC    右函数 
    ! @pram: t0     积分初始时间
    ! @pram: tf     积分结束时间
    ! @pram: tspan  积分时间间隔（内循环，不输出数据），间隔tspan之后输出一次数据
    ! @pram: err    误差限
    ! @pram: hh     初始步长
    ! @pram: Num    积分时间间隔（外循环，输出数据），总计输出Num个数据
    ! @pram: x0     积分初始状态量和结束时输出状态量
    use global
    USE FindRoots
    external :: YHC
    integer :: FILEID
    REAL(8) :: t0,tf,err,step,Energy,w0(2),hh,tspan,x0(4),EE,the(100000000),max
    REAL(8) :: S,theta,S1,theta1,thetaA1,thetaB1,thetaB
    REAL(8) :: axis,ecc,placeholder
102 format(1X, E25.18,1X,E25.18)    
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
104 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
105 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
106 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
107 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
109 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18)
110 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
111 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18)
113 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X, E25.18,1X,E25.18)
    do I=1,Num 
        J = 1
        do while (t0<tf)
            S = x0(1)
            theta = x0(2)
            S1 = x0(3)
            theta1 = x0(4)
            thetaA1 = thetaA1_Const
            thetaB1 = (K - m*S**2D0*theta1)/(m*S**2D0+IzB)
            thetaB = thetaB1*t0
            Energy = m / 2D0 * (S1 ** 2D0 + S ** 2D0 * (theta1+thetaB1) ** 2D0) + 1D0 / 2D0 * IzB * thetaB1 ** 2D0 - m / S - m / (S ** 3D0) * (A1 + A2 * dcos(2D0 * theta))
            the(J) = theta
            CALL mco_x2el (1D0,S*dcos(thetaB+theta),S*dsin(theta+thetaB),0D0,S1*dcos(theta+thetaB)-S*dsin(thetaB+theta)*(theta1+thetaB1),S1*dsin(theta+thetaB)+S*dcos(theta+thetaB)*(theta1+thetaB1),0D0,axis,ecc,placeholder,placeholder,placeholder,placeholder)
            !axis = 1D0/(2D0/S-(S1**2D0+S**2D0*(theta1+thetaB1)**2D0))
            !ecc=sqrt(1-S**4D0*(theta1+thetaB1)**2D0/(axis))
            J = J+1
            Call RKF78(YHC,hh,t0,x0,EE,err,4)
            if (S<1D0) then
                write(*,*) 'breakup'
                GO TO 21
            end if
            write(FILEID,107) x0,axis,ecc,t0*Tunit
        enddo
        max = maxval(abs(the(1:J-1)))
        if (max>PI.OR.S<1D0) then
            write(*,*) 'breakup'
            GO TO 21
        end if
        write(*,*) 'I = ',I ,'/',Num, ';','T = ',t0*Tunit,'year',';','S0 = ',S0
        !write(FILEID,107) x0,energy,t0*Tunit,max
        tf=t0+tspan
    enddo
21  write(*,*) 'DONE'

    end
    
    SUBROUTINE PincareX0(xi0,eta0,x0,Energy)
    ! @sub Getx02   求解两自由度的初始状态量
    ! @input: index  长短周期轨道指示器 
    ! @input: xi0    初始横向振幅
    ! @input: eta0   初始纵向振幅
    ! @output: x0     输出状态量
    ! @output: w0     输出长短周期轨道频率
    use global
    USE FindRoots
    integer(4) :: index
    REAL(8) :: x0(4),xi0,eta0,x1,x2,x3,x4,Energy
    REAL(8) :: min,coe
    
    x1 = S0 + xi0
    x2 = 0D0
    x3 = 0D0
    
    !求解新theta1
    min = 1D0/2D0*m*x3**2D0 - m/x1 - m*(A1+A2)/x1**3D0 &
        + (K**2D0)/(2D0*(m*x1**2D0 + IzB))
    coe = (IzB*m*x1**2D0)/(2D0*(m*x1**2D0 + IzB))
    x4 = -dsqrt((Energy-min)/coe)

    x0 = (/x1, x2, x3, x4/)
    return
    end
    
    SUBROUTINE Pincare(FILEID,YHC,tf,dx,Num,err,Energy)
    ! @sub Pincare      求解庞加莱截面
    ! @input: FILEID    数据输出文件标识符 
    ! @input: tf        最长积分时间
    ! @input: dx        振幅或者庞加莱截面半径
    ! @output: Num      绘制轨道的数目     
    ! @output: err      求解精度
    use global
    external :: YHC
    REAL(8) :: t0,tf,hh,dx,x0(4),w0(2),err,EE,xmid(4),ymid,xx(4),xi0,eta0
    REAL(8) :: S,theta,S1,theta1,Energy
    integer :: Num,FILEID
102 format(1X, E25.18,1X,E25.18)    
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
104 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
105 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
106 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
107 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
109 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18)
110 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
111 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18)
113 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X, E25.18,1X,E25.18)
    
    hh = dx/Num
    DO I =1,Num
        !给出初值
        xi0 = (I-1) * hh
        eta0 = 0D0
        t0 = 0D0
        CALL PincareX0(xi0,eta0,x0,Energy)
        !CALL Getx02(2, xi0, eta0,x0,w0)
        !write(*,*) 'x0'
        !write(*,104) x0
        DO while (t0<tf)
            
            !记录中间数据
            xmid = x0
            ymid = xmid(1)*dsin(xmid(2))
            tmid = t0
            
            !积分
22          Call RKF78(YHC, hh, t0, x0, EE, 1D0, 4)
            
            !转换到直角坐标系
            xx(1) = x0(1)*dcos(x0(2))
            xx(2) = x0(1)*dsin(x0(2))
            xx(3) = x0(3)*dcos(x0(2))-x0(1)*dsin(x0(2))*x0(4)
            xx(4) = x0(3)*dsin(x0(2))+x0(1)*dcos(x0(2))*x0(4)
            write(13,104) xx
            !write(*,104) xx
            
            !判定是否有解
            !if (hh<err) then
            !    hh = dx/Num
            !    GOTO 23
            !ENDIF 
            !if (ymid.GE.-1D-15 .AND. xx(2).LE.1D-15 .AND. xx(4)<0D0) then
            if (ymid*xx(2).LE.1D-15.AND. xx(4)<0D0) then
                if (ABS(ymid-xx(2))<err) then
                    write(FILEID,104) xx
                    write(*,104) xx
                    hh = dx/Num
                    GOTO 23
                endif
                hh = hh/2D0
                x0 = xmid
                t0 = tmid
                GOTO 22
            end if
            
             
23      ENDDO
        !write(13,*) ' '
    ENDDO
    end
    
    SUBROUTINE FindDot(YHC,t0,tf,err,DD,hh1,x4)
    ! @sub FindDot      寻找映射点
    ! @input: FILEID    数据输出文件标识符 
    ! @input: tf        最长积分时间
    ! @input: dx        振幅或者庞加莱截面半径
    ! @output: Num      绘制轨道的数目     
    ! @output: err      求解精度
    USE global
    external :: YHC
    REAL(8) :: t0,tf,err,hh1,x4(4),xmid(4),tmid,DD,xrecord(4)
    REAL(8) :: xree(4),tree,hree,alpha,beta,phi1,phi2
102 format(1X, E25.18,1X,E25.18)    
103 format(1X, E25.18,1X,E25.18,1X,E25.18)
104 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
105 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
    write(16,104) x4
    tree = t0
    xree = x4
    hree = hh1
    do while (tree<tf)
22      xmid = xree
        tmid = tree
        Call RKF78(YHC, hree, tree, xree, err, DD, 4)
        write(12,104) xree
        if (xree(1)*dsin(xree(2))*xmid(1)*dsin(xmid(2)) < 0D0 .AND. xree(4)<0) then
            if (ABS(xree(1)*dsin(xree(2)))<1D-14) then
                write(16,104) xree
                CALL GetAB_re(xree,alpha,beta,phi1,phi2)
                write(15,104) alpha/S0,beta*18D1/PI,phi1,phi2
                hree = hh1
                GOTO 22
            end if
            xree = xmid
            tree = tmid
            hree = hree/1D1
        end if
    end do
    end
    
    SUBROUTINE GetAB_re(x0,al,be,phi1,phi2)
    use global
    REAL(8) :: w0(2),alpha(2),x0(4),xi,eta,xi1,eta1,s1,s2,al,be,phi1,phi2
    
    CALL Getw0(w0,alpha)
    xi = x0(1) - S0
    eta = x0(2)
    xi1 = x0(3)
    eta1 = x0(4)
    phi1 = ATAN2(((w0(2)*eta-alpha(2)*xi1)/(w0(2)*alpha(2)*xi+eta1)),((w0(1)*alpha(2)-alpha(1)*w0(2))/(w0(2)*alpha(2)-w0(1)*alpha(1))))
    phi2 = ATAN2(((w0(1)*eta-alpha(1)*xi1)/(alpha(1)*w0(1)*xi+eta1)),((alpha(1)*w0(2)-alpha(2)*w0(1))/(alpha(1)*w0(1)-alpha(2)*w0(2))))
    al = dsqrt(((alpha(1)*w0(1)*xi+eta1)/(alpha(1)*w0(1)-alpha(2)*w0(2)))**2D0+((alpha(1)*xi1-w0(1)*eta)/(alpha(2)*w0(1)-alpha(1)*w0(2)))**2D0)
    be = dsqrt(((alpha(1)*(w0(2)*alpha(2)*xi+eta1))/(alpha(2)*w0(2)-alpha(1)*w0(1)))**2D0+((alpha(1)*(w0(2)*eta-alpha(2)*xi1))/(alpha(1)*w0(2)-alpha(2)*w0(1)))**2D0)
    end
    
    SUBROUTINE GetAB_re2(x0,al,be)
    use global
    REAL(8) :: w0(2),alpha(2),x0(4),xi,eta,xi1,eta1,s1,s2,al,be
    
    CALL Getw0(w0,alpha)
    xi = x0(1)-S0
    eta = 0D0!x0(2)
    xi1 = 0D0!x0(3)
    eta1 = x0(4)
    al = (alpha(1)*w0(1)*xi+eta1)/(alpha(1)*w0(1)-w0(2)*alpha(2))
    be = alpha(1)*(eta1+w0(2)*alpha(2)*xi)/(w0(2)*alpha(2)-w0(1)*alpha(1))
    end
    