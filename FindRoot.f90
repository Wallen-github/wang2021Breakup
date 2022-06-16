!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! @file: FindRoot.f90
! @brief: ��module�����������ⷽ�������ַ���ţ�ٵ����������������ⷽ��
! @author: Wang Hai-Shuo
! @time: 2020.3.3
! @pram: global ϵͳ��������ͷ�ļ�
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MODULE FindRoots
    IMPLICIT NONE
    CONTAINS
    
    SUBROUTINE equ(x,y)
    ! @sub equ ���ⷽ�̣���֪ϵͳ����K���ϵͳƽ���r0
    ! @pram: x  ƽ���r0
    ! @pram: y  ���������
        USE global
        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        REAL(iwp),INTENT(IN)::x
        REAL(iwp),INTENT(INOUT)::y
        y = (K / (IzB + m * x ** 2D0)) ** 2D0 * x ** 5D0 - x ** 2D0 - 3D0 * (A1 + A2)
    END SUBROUTINE
    
    SUBROUTINE equE(theta1,y)
    ! @sub equ ���ⷽ�̣���֪ϵͳ����K���ϵͳƽ���r0
    ! @pram: x  ƽ���r0
    ! @pram: y  ���������
        USE global
        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        REAL(iwp),INTENT(IN)::theta1
        REAL(iwp),INTENT(INOUT)::y
        y = (1D0/2D0)*m*(S1Const**2D0 + SConst**2D0 * (theta1-(m * r**2D0 * theta1 + IzA * thetaA1_Const-K) / (m * SConst ** 2D0 + IzB))**2D0) &
            + (1D0/2D0) * IzB * (m * SConst ** 2D0 * theta1+IzA*thetaA1_Const-K)**2D0/(m*SConst**2D0+IzB)**2D0 &
            + (1D0/2D0) * IzA * thetaA1_Const ** 2D0 &
            - m/SConst - m * (A1+A2)/SConst ** 3D0 &
            - Energy_Const
    END SUBROUTINE
        
    SUBROUTINE bisection(n,sec,tol,rt)
    ! @sub bisection ���ַ�
    ! @pram: n      ��ĸ���
    ! @pram: sec    �������
    ! @pram: tol    �����
    ! @pram: rt     ���������
        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        INTEGER,INTENT(IN)::n
        REAL(iwp),INTENT(IN)::sec(:,:)
        REAL(iwp),INTENT(IN)::tol
        REAL(iwp),INTENT(INOUT)::rt(:)
        REAL(iwp)::x1,x2,x3,y1,y2,y3
        INTEGER::i
        DO i=1,n
            x1=sec(1,i); x2=sec(2,i)
            CALL equ(x1,y1); call equ(x2,y2)
            IF(y1*y2<=0.0_iwp)THEN
                DO WHILE(x2-x1>tol)
                    x3=(x2+x1)/2.
                    CALL equ(x3,y3)
                    IF(y3*y2<=0.0_iwp)THEN
                        x1=x3; y1=y3
                    ELSE
                        x2=x3; y2=y3
                    END IF
                END DO
                rt(i)=(x1+x2)/2.
            ELSE
                WRITE(*,'(a,i1,a)')"ROOT MIGHT NOT IN THE ",i," SECTION"
                rt(i)=0
            END IF
        END DO
    END SUBROUTINE

    SUBROUTINE NewRaf(equ,n,x,tol,rt)
    ! @sub NewRaf   ţ�ٷ�
    ! @pram: equ    ����ⷽ��
    ! @pram: n      ��ĸ���
    ! @pram: x      ��ֵ
    ! @pram: tol    �����
    ! @pram: rt     ���������
        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        INTEGER,INTENT(IN)::n
        REAL(iwp),PARAMETER::dx=1.0e-3
        REAL(iwp),INTENT(IN)::x(:)
        REAL(iwp),INTENT(IN)::tol
        REAL(iwp),INTENT(INOUT)::rt(:)
        REAL(iwp)::x1,x2,y1,yd,dy1
        INTEGER::i
        DO i=1,n
            x1=x(i)
            x2=x1
            DO
                x1=x2
                CALL equ(x1,y1)
                CALL equ(x1+dx,yd)
                dy1=(yd-y1)/dx
                x2=x1-y1/dy1
                IF(ABS(x2-x1)<tol) EXIT
            END DO
            rt(i)=(x2+x1)/2.
        END DO
    END SUBROUTINE
    
    SUBROUTINE NewRafPar(equ,n,x,tol,rt,par)
    ! @sub NewRafPar   �����������ţ�ٷ�
    ! @pram: equ    ����ⷽ��
    ! @pram: n      ��ĸ���
    ! @pram: x      ��ֵ
    ! @pram: tol    �����
    ! @pram: rt     ���������
    ! @pram: par    �������
        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        INTEGER,INTENT(IN)::n
        REAL(iwp),PARAMETER::dx=1.0e-3
        REAL(iwp),INTENT(IN)::x(:),par(:)
        REAL(iwp),INTENT(IN)::tol
        REAL(iwp),INTENT(INOUT)::rt(:)
        REAL(iwp)::x1,x2,y1,yd,dy1
        INTEGER :: i,flag
        flag = 1
        DO i=1,n
            x1=x(i)
            x2=x1
            DO
                flag = flag+1
                x1=x2
                CALL equ(x1,y1,par)
                CALL equ(x1+dx,yd,par)
                dy1=(yd-y1)/dx
                x2=x1-y1/dy1
                IF(ABS(x2-x1)<tol) EXIT
                IF(FLAG>50) EXIT
            END DO
            rt(i)=(x2+x1)/2.
        END DO
    END SUBROUTINE
END MODULE