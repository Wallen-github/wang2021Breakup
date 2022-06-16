	SUBROUTINE ANTIMATRIX(A,Anti_A,N,L)
C-------------------------------------------------------------------
C	THE SUBROUTINE COMPUTE THE INVERSE OF THE MATRIX A
C	A------------THE INPUT  MATRIX 
C	Anti_A       the output matrix
C	N------------THE DIMENSION NUMBER OF THE MATRIX
C	L----------L=0, NO INVERSE;L=1,RETURN THE RIGHT INVERSE.
C--------------------------------------------------------------------
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(N,N),IS(N),JS(N)
	real*8 Anti_A(N,N)
	DO I=1,N
	  DO J=1,N
	    Anti_A(i,j) = A(i,j)
	  ENDDO
	ENDDO
	L=1
	DO 100 K=1,N
	  D=0D0
 	  DO 10 I=K,N
	  DO 10 J=K,N
	    IF(DABS(Anti_A(I,J)).GT.D)THEN
	      D=DABS(Anti_A(I,J))
	      IS(K)=I
	      JS(K)=J
	    ENDIF
10	  CONTINUE
        IF(D+1D0.EQ.1D0)THEN
	    L=0
	    WRITE(*,20)
	    RETURN
	  ENDIF
20	  FORMAT(1X,'ERR**NOT INVERSE!')
        DO 30 J=1,N
	  T=Anti_A(K,J)
	  Anti_A(K,J)=Anti_A(IS(K),J)
	  Anti_A(IS(K),J)=T
30	  CONTINUE
        DO 40 I=1,N
	  T=Anti_A(I,K)
	  Anti_A(I,K)=Anti_A(I,JS(K))
	  Anti_A(I,JS(K))=T
40	  CONTINUE
    	  Anti_A(K,K)=1D0/Anti_A(K,K)
	  DO 50 J=1,N
	    IF(J.NE.K)THEN
	      Anti_A(K,J)=Anti_A(K,J)*Anti_A(K,K)
	    ENDIF
50	  CONTINUE
        DO 70 I=1,N
	    IF(I.NE.K)THEN
	    DO 60 J=1,N
	      IF(J.NE.K)THEN
	        Anti_A(I,J)=Anti_A(I,J)-Anti_A(I,K)*Anti_A(K,J)
	      ENDIF
60	    CONTINUE
        ENDIF
70	  CONTINUE
        DO 80 I=1,N
	    IF(I.NE.K)THEN
	      Anti_A(I,K)=-Anti_A(I,K)*Anti_A(K,K)
	    ENDIF
80	  CONTINUE
100	CONTINUE
	DO 130 K=N,1,-1
	  DO 110 J=1,N
	    T=Anti_A(K,J)
	    Anti_A(K,J)=Anti_A(JS(K),J)
	    Anti_A(JS(K),J)=T
110	  CONTINUE
        DO 120 I=1,N
	    T=Anti_A(I,K)
	    Anti_A(I,K)=Anti_A(I,IS(K))
	    Anti_A(I,IS(K))=T
120	  CONTINUE
130	CONTINUE
	RETURN
	END

C---------------------------------------------------------------------------------------------
C	This subprogram calculates the product of two vectors.
C---------------------------------------------------------------------------------------------
	SUBROUTINE VECTOR(A,B,C,N)
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(N),B(N)
	C=0D0
	DO 102 I=1,N
	  C=C+A(I)*B(I)
102	CONTINUE
	RETURN
	END

C---------------------------------------------------------------------------------------------
C	This subprogram calculates the product of a matrix and a vector with matrix leftside.
C---------------------------------------------------------------------------------------------
	SUBROUTINE MATRIXVECTOR(A,B,C,M,N)
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(M,N),B(N),C(M)
	DO 101 I=1,M
	  C(I)=0D0
	  DO 102 J=1,N
		C(I)=C(I)+A(I,J)*B(J)
102	  CONTINUE
101	CONTINUE
	RETURN
	END

C---------------------------------------------------------------------------------------------
C	This subprogram calculates the product of a vector and a matrix with matrix rightside.
C---------------------------------------------------------------------------------------------
	SUBROUTINE VECTORMATRIX(A,B,C,M,N)
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(M),B(M,N),C(N)
	DO 101 I=1,N
	  C(I)=0D0
	  DO 102 J=1,M
		C(I)=C(I)+A(J)*B(J,I)
102	  CONTINUE
101	CONTINUE
	RETURN
	END

C---------------------------------------------------------------------------------------------
C	This subprogram calculates the product of two matrixes.
C---------------------------------------------------------------------------------------------
	SUBROUTINE MATRIXMATRIX(A,B,C,L,M,N)
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(L,M),B(M,N),C(L,N)
	DO 101 I=1,L
	DO 102 J=1,N
	  C(I,J)=0D0
	  DO 201 K=1,M
		C(I,J)=C(I,J)+A(I,K)*B(K,J)
201	  CONTINUE
102	CONTINUE
101	CONTINUE
	RETURN
	END


	SUBROUTINE ROTX(ANG,ROT)
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION ROT(3,3)
	ROT(1,1)=1D0
	ROT(1,2)=0D0
	ROT(1,3)=0D0
	ROT(2,1)=0D0
	ROT(2,2)=DCOS(ANG)
	ROT(2,3)=DSIN(ANG)
	ROT(3,1)=0D0
	ROT(3,2)=-DSIN(ANG)
	ROT(3,3)=DCOS(ANG)
	RETURN
	END


	SUBROUTINE ROTY(ANG,ROT)
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION ROT(3,3)
	ROT(1,1)=DCOS(ANG)
	ROT(1,2)=0D0
	ROT(1,3)=-DSIN(ANG)
	ROT(2,1)=0D0
	ROT(2,2)=1D0
	ROT(2,3)=0D0
	ROT(3,1)=DSIN(ANG)
	ROT(3,2)=0D0
	ROT(3,3)=DCOS(ANG)
	RETURN
	END


	SUBROUTINE ROTZ(ANG,ROT)
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION ROT(3,3)
	ROT(1,1)=DCOS(ANG)
	ROT(1,2)=DSIN(ANG)
	ROT(1,3)=0D0
	ROT(2,1)=-DSIN(ANG)
	ROT(2,2)=DCOS(ANG)
	ROT(2,3)=0D0
	ROT(3,1)=0D0
	ROT(3,2)=0D0
	ROT(3,3)=1D0
	RETURN
	END