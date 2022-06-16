    MODULE global
    implicit none
    !REAL(8) :: rB = 4D2,rA = 1.2D3,p = 2D3,R = 2.6D3, a_b = 1.32D0, b_c = 1.0D0 !KW4
    !REAL(8) :: rB = 163D0,rA = 780D0,p = 2100D0,R = 2100D0, a_bB = 1.5D0, b_cB = 1.2D0,a_bA = 1.2D0, b_cA = 1.2D0 !Didymos
    !REAL(8) :: rB = 8.58D2,rA = 3.3D3,p = 2D3,R = 5.94D3, a_b = 1.54D0, b_c = 1.0D0 !JQ58
    !REAL(8) :: rB = 2.46D3,rA = 8.2D3,p = 2D3,R = 1.968D4, a_b = 1.31D0, b_c = 1.0D0 !Mayall
    !REAL(8) :: rB = 1.26D3,rA = 3.6D3,p = 2D3,R = 7.92D3, a_b = 1.3D0, b_c = 1.0D0 !Kiuchi
    REAL(8) :: rB = 4.5D2,rA = 8D2,p = 2100D0,R = 5D3, a_bB = 1.2D0, b_cB = 1.2D0, a_bA = 1.2D0, b_cA = 1.2D0!µÚÒ»ÖÖ
    !EAL(8) :: rB = 486.5761596D0,rA = 973.1523D0,p = 2D3, R = 3D3, a_bB = 5D2/48D1 , b_cB = 1D0, a_bA = 1D3/96D1, b_cA = 1.0D0 !Case #1
    REAL(8) :: au = 1.496D11,G = 6.67D-11,SecOfYear = 31536D3
    REAL(8) :: cB,bB,aB,cA,bA,aA,Length,S0,MA,MB,mu,m,alphaB,alphaA,J2B,J2A,J22B,J22A, A1,A2,A3,IzB,IzA,K,Tunit,thetaA1_Const,Energy_Const,SConst,S1Const
    REAL(8) :: PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679D0
    end MODULE global
     