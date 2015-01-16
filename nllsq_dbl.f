C*************************************COMMON /DJA1/iteration, npts*********************************
C                                                                   ***
C      DOUBLE PRECISION VERSION OF NLLSQ                            ***
C                                                                   ***
C      Daniel Auerbach  2013-04-30                                  ***
C**********************************************************************
C      Modified the nllsq.f single precision version as follows     ***
C      IMPLICIT REAL(8) (A-H,O-Z) added to all procedures           ***
C                                                                   ***
C                                                                   ***
C**********************************************************************

CDJALABEL=NLLSQ.SGL  NON LINEAR LEAST SQUARES FROM BELL LABS
CDJALABEL MODIFIED AS FOLLOWS
C
C       LAST PARAMETER IS A VARIABLE PRINTED AS 10,A5 ON THE OUTPUT
C
C       FINAL PRINTOUT 1 => NO LIST OF VALUES AND RESIDUALS
C
C DJA 2013-0

      SUBROUTINE NLLSQ(Y,X,BB,RES,NARRAY,ARRAY,IBB,FILE)
C     NONLINEAR LEAST SQUARES FITTING ALGORITHM BY D W MARQUARD
C     ORIGINAL PROGRAM REWRITTEN BY W A BURNETTE  BTL JULY  1967
C     NARRAY CONTAINS PROGRAM PARAMETERS  ARRAY CONTAINS STATISTICAL
C     CONSTANTS   SET ARRAY EQUAL TO 0.FOR STANDARD SET OF CONSTANTS
C     MAXIMUM NUMBER OF PARAMETERS IS 20 THIS MAY BE CHANGED BY ALTERING
C     DIMENSION STATEMENTS AND MATRIX STORING STATEMENTS
C
C  DESCRIPTION OR 'ARRAY'
C
C   ELEMENT           PARAMETER                      DEFAULT
C
C   ARRAY(1)  ....    AL                             .1
C   ARRAY(2)  ....    DELTA - FOR DERIVATIVES        .00001
C   ARRAY(3)  ....    E     - CONVERGENCE CRITERION  .00005
C   ARRAY(4)  ....    FF                             4.0
C   ARRAY(5)  ....    GAMCR - CRITICAL ANGLE         45.
C   ARRAY(6)  ....    T                              2.
C   ARRAY(7)  ....    TAU                            .001
C   ARRAY(8)  ....    ZETA  - CRITERION FOR          1E-31
C                             SINGULAR MATRIX
C
C
C  DESCRIPTION OR 'NARRAY'
C
C   ELEMENT           PARAMETER
C
C   NARRAY(1) ....    N     - NUMBER OF DATA POINTS
C   NARRAY(2) ....    M     - NUMBER OF INDEPENDENT VARIABLES
C   NARRAY(3) ....    K     - ADJUSTABLE PARAMETERS
C   NARRAY(4) ....    IP    - NUMBER OF OMIITED PARAMETERS (CONSTANT)
C   NARRAY(5) ....    INTERMEDITATE PRINTOUT OPTION (0 TO 3)
C   NARRAY(6) ....    FINAL PRINTOUT OPTION         (-1TO 2)
C   NARRAY(7) ....    PRINT OUT UNIT NUMBER
C   NARRAY(8) ....    KITER - MAXIMUM NUMBER OF ITERATIONS
C
C
C**********************************************************************
C     Change impliit type of from REAL to REAL(8)                   ***
C**********************************************************************
      IMPLICITREAL(8) (A-H,O-Z)


      COMMON/BLK1/B(20),P(20),RE,N,M,K
      COMMON/BLK2/A(40,20),SA(20),K2,IK,NPMAX
      COMMON/BLK3/BS(20),DB(20),G(20),K3
      COMMON/BLK4/AL,DELTA,E,FF,GAMCR,T,TAU,ZETA,PHI,SE,PHICR
      COMMON/BLK5/IB(20),IP
      COMMON/DJA1/LJ
      DIMENSION Y(5000),X(5000,3,1000),RES(5000)
      DIMENSIONBB(20),IBB(20)
      DIMENSIONNARRAY(8),ARRAY(8)
      DIMENSION CONST(8),SCONST(8)
C
C  *****  FILE IS ARRAY USED FOR PRINTING HEADING
C  *****     DIM MUST MATCH FORMAT 747 AND MAIN PROG
C
C      DIMENSION FILE(5)
      character(20) FILE
C
C
      EQUIVALENCE(CONST(1),AL)
      DATA SCONST/0.1,1.0E-5,5.0E-5,4.0,45.0,2.0,0.001,1.0E-31/
      NPMAX=20
      N=NARRAY(1)
      M=NARRAY(2)
      K=NARRAY(3)
      K2=K

      K3=K
      IP=NARRAY(4)
      INTP=NARRAY(5)
      IFP=NARRAY(6)
      IK=NARRAY(7)
      IF(NARRAY(8).EQ.(-1)) NARRAY(8)=KITER
      KITER=NARRAY(8)
C     WHICH OF THE CONSTANTS HAVE BEEN DETERMINED BY USER
      DO 5 J=1,8
      IF(ARRAY(J).LE.0.)GO TO 4
      CONST(J)=ARRAY(J)
      GO TO 5
  4   CONST(J)=SCONST(J)
  5   CONTINUE
      IF(KITER.LE.0)KITER=30
      IF(IK.LE.0)IK=6
C label not used   26   IF(IFP.EQ.(-1))GO TO 100
      IF(IFP.EQ.(-1))GO TO 100
      WRITE(IK,2090)
      WRITE(IK,747) FILE
      WRITE(IK,2001)N,K,M,DELTA,E,FF,GAMCR,T,TAU,ZETA,AL
 100  DO120J=1,K
      BS(J)=BB(J)
      B(J)=BB(J)
  120 CONTINUE
      DO 121 J=1,IP
  121 IB(J)=IBB(J)
      CALL SUMSQ(PHI,Y,X,RES)
      LJ=0
 130  IF(LJ.GE.KITER)GOTO404
      LJ=LJ+1
      print *, lj
C     BEGIN LJTH ITERATION
      CALL XNEWAZ(Y,X,RES)
 1311 IF(AL.LT..1E-07) GO TO 131
      AL=AL/10.
 131  CALLXSCALX
      PHIOLD=PHI
C     STORE MATRIX
      DO132I=1,K
      II=I+NPMAX
      DO132J=1,K
 132  A(II,J)=A(I,J)
      CALLXSOLVX
 135  DO140J=1,K
  140 B(J)=BS(J)+DB(J)
C     COMPUTE GAMMA
C label not used   150  DD=0.
      DD=0.
      DG=0.
      GG=0.
      DO152J=1,K
      IF(SA(J).EQ.0.)GOTO152
      GG=GG+G(J)*G(J)/(SA(J)*SA(J))
      DD=DD+DB(J)*DB(J)*SA(J)*SA(J)
  152 DG=DG+DB(J)*G(J)
      XL=SQRT(DD)
      IF(DD*GG.GT.0.)GOTO160
C label not used    155  GAMMA=0.
      GAMMA=0.
      GOTO170
 160  CGAM=DG/SQRT(DD*GG)
      GAMMA=0
      IF(CGAM*CGAM.GE.1.)GO TO 170
      WS=SQRT(1.-CGAM*CGAM)
      GAMMA=57.2957795*ATAN2(WS,CGAM)
 170  CALL SUMSQ(PHI,Y,X,RES)
C label not used   171 IF(PHI.LE.PHIOLD)GOTO175
      IF(PHI.LE.PHIOLD)GOTO175
      IF(GAMMA-GAMCR)300,300,180
  175 DO 176 J=1,K
  176 BS(J)=B(J)
      IF (GAMMA.LT.90.) GO TO 190
C     GAMMA LAMBDA TEST
C label not used 178  IF(AL-1.)190,403,403
      IF(AL-1.)190,403,403
  180 AL=AL*10.
      CALLXSOLVX
      GOTO135
C     EPSILON TEST
  190 CALL XEPTST(L)
      GOTO(401,200),L
C     BEGIN INTERMEDIATE OUTPUT ROUTINE
 200  IF(INTP.EQ.0)GOTO130
      WRITE(IK,2000)
      WRITE(IK,2002)LJ,PHI,AL,(B(J),J=1,K)
      WRITE(IK,2003)GAMMA,XL,(DB(J),J=1,K)
      IF(INTP.EQ.1)GOTO130
      CALL XNEWAZ(Y,X,RES)
C     STORE MATRIX
      DO205I=1,K
      II=I+NPMAX
      DO205J=1,K
 205  A(II,J)=A(I,J)
      CALL  GJR (MS)
      GO TO (207,130),MS
 207  IF(INTP.EQ.2)GOTO210
      WRITE(IK,2004)
      CALLXPRT1X
 210  CALLXSCALX
      WRITE(IK,2006)
      CALLXPRT2X
C     GET MATRIX FROM STORAGE
      DO220I=1,K
      II=I+NPMAX
      DO220J=1,K
 220  A(I,J)=A(II,J)
      IF(LJ.GE.KITER)GO TO 404
      LJ=LJ+1
      GO TO 1311
 300  DO320J=1,K
      DB(J)=DB(J)/2.
  320 B(J)=BS(J)+DB(J)
C     GAMMA EPSILON TEST
      CALLXEPTST(L)
      GO TO (402,321),L
  321 CALL  SUMSQ(PHI,Y,X,RES)
      IF(PHIOLD.LT.PHI) GO TO 300
      DO330J=1,K
 330  BS(J)=B(J)
      GO TO 200
C     BEGIN FINAL  PRINTOUT ROUTINE
 401  IF(IFP.EQ.(-1))GO TO 600
      WRITE(IK,2090)
      WRITE(IK,747) FILE
      WRITE(IK,2010)
      GOTO405
  402 IF(IFP.EQ.(-1))GO TO 4025
      WRITE(IK,2090)
      WRITE(IK,747) FILE
      WRITE(IK,2011)
 4025 IF(IFP.NE.(-1).AND.PHIOLD.LT.PHI)WRITE(IK,2092)
      IF(PHIOLD.GE.PHI)GO TO 4029
      PHI=PHIOLD
      DO 4027 J=1,K
 4027 B(J)=BS(J)
 4029 IF(IFP.EQ.(-1)) GO TO 600
      GOTO405
 403  IF(IFP.EQ.(-1))GO TO 600
      WRITE(IK,2090)
      WRITE(IK,747) FILE
      WRITE(IK,2012)
      GOTO405
  404 WRITE(IK,2013)
      NARRAY(8)=-1
      IF(IFP.EQ.(-1))GO TO 600
      WRITE(IK,747) FILE
 405  DO406J=1,K
      BS(J)=B(J)
 406  BB(J)=B(J)
      WRITE(IK,2002)LJ,PHI,AL,(B(J),J=1,K)
      WRITE(IK,2003)GAMMA,XL,(DB(J),J=1,K)
      CALL XNEWAZ(Y,X,RES)
      IF(IFP.LE.1)GOTO430
      DO410I=1,K
      II=I+NPMAX
      DO410J=1,K
  410 A(II,J)=A(I,J)
      WRITE(IK,2022)
      CALL XPRT1X
      CALL XSCALX
      WRITE(IK,2023)
      CALL XPRT2X
C     GET MATRIX FROM STORAGE
      DO420I=1,K
      II=I+NPMAX
      DO420J=1,K
 420  A(I,J)=A(II,J)
 430  CALL  GJR (MS)
      GO TO (440,435),MS
 435  WRITE(IK,2060)
      GO TO 455
 440  IF(IFP.EQ.0)GOTO450
      WRITE(IK,2024)
      CALL XPRT1X
 450  CALL XSCALX
      WRITE(IK,2025)
      CALL XPRT2X
  455 IF(IFP.EQ.0) GO TO 461
      WRITE(IK,2030)
      DO460I=1,N
C      RESIDUAL ARRAY OPTION SATISFIED HERE
      J=4
      CALL MODEL(F,Y,X,RES,I,J)
  460  WRITE(IK,2031)I,X(I,:,1),Y(I),F,RE
C     ONE PARAMETER SUPPORT PLANE COMPUTATIONS
461      FNKW=N-K+IP
      IF(FNKW.LE.0.)GOTO589
      FKW=K-IP
      SE=SQRT(PHI/FNKW)
      WRITE(IK,2040)
      DO470I=1,K
C     CHECK FOR OMITTED PARAMETERS
      IF(IP.EQ.0)GOTO464
      DO462J=1,IP
      IF (I.EQ.IB(J)) GO TO 469
 462  CONTINUE
  464 STE=SA(I)*SE
      HJTD=SQRT(FF*FKW)*STE
      OPL=BS(I)-STE*T
      OPU=BS(I)+STE*T
      SPL=BS(I)-HJTD
      SPU=BS(I)+HJTD
      WRITE(IK,2041)I,STE,OPL,OPU,SPL,SPU
      GOTO470
 469  WRITE(IK,2042)I
 470  CONTINUE
      IF(IFP .EQ. 0) GOTO 590
      IF (IFP.EQ.1) GO TO 602
C     NONLINEAR CONFIDENCE REGION CALCULATIONS
      WS=FKW/FNKW
      PHICR=PHI*(1.+WS*FF)
      WRITE(IK,2049)PHICR
      CALL XCNFXX(Y,X,RES)
      IF(IFP.GE.0)WRITE(IK,2090)
      RETURN
 589  WRITE(IK,2060)
  590 IF(IFP.EQ.0) GO TO 602
 599  IF(IFP.GE.0)WRITE(IK,2090)
      RETURN
C     RETURNING PARAMETERS WITH NO OUTPUT
 600  DO601J=1,K
 601  BB(J)=B(J)
C      RESIDUAL ARRAY OPTION WITH NO OUTPUT
 602  J=4
      CALLMODEL(F,Y,X,RES,1,J)
      GOTO(599,599,599,604),J
 604  DO605I=2,N
 605  CALLMODEL(F,Y,X,RES,I,J)
      IF(IFP.GE.0)WRITE(IK,2090)
      RETURN
 2000 FORMAT(100H XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     1XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX/1H )
 2001 FORMAT(21H NO OF DATA POINTS IS,I4,22H   NO OF PARAMETERS IS,I3,3X
     1,30HNO OF INDEPENDENT VARIABLES IS,I3/7H DELTA=,E15.8,5H   E=,E15.
     28,6H   FF=,E15.8,9H   GAMCR=,E15.8/5H   T=,E15.8,7H   TAU=,E15.8,
     38H   ZETA=,E15.8,6H   AL=,E15.8 )
 2002 FORMAT(1X/19H NO OF ITERATIONS =I3/6H PHI =,E15.8,4X,8HLAMBDA =,
     1E15.8/11H PARAMETERS/(1X,7E17.8/))
 2003 FORMAT(8H GAMMA =,E15.8,4X,14HLENGTH OF DB =,E15.8/21H DB CORRECTI
     1ON VECTOR/(1X,7E17.8/))
 2004 FORMAT(1X/12H PTP INVERSE)
 2006 FORMAT(1X/25H CORRELATION COEFFICIENTS)
 2010 FORMAT(28H CONVERGENCE BY EPSILON TEST)
 2011 FORMAT(34H CONVERGENCE BY GAMMA EPSILON TEST)
 2012 FORMAT(33H CONVERGENCE BY GAMMA LAMBDA TEST)
 2013 FORMAT(10H FORCE OFF)
 2022 FORMAT(1X/11H PTP MATRIX)
 2023 FORMAT(1X/29H PTP CORRELATION COEFFICIENTS)
 2024 FORMAT(1X/12H PTP INVERSE)
 2025 FORMAT(1X/35H PARAMETER CORRELATION COEFFICIENTS)
 2030 FORMAT(1X/1H ,15X,'X',9X,11H   OBSERVED,11X,9HPREDICTED,10X,8HRESI
     1DUAL/1X)
 2031 FORMAT(1X,I4,3X,3(E15.8,4X),E15.8)
 2040 FORMAT(1X/1H ,12X,4H STD,18X,15H NE - PARAMETER,22X,13HSUPPORT PLA
     1NE/2X,4HPARA,7X,5HERROR,13X,5HLOWER,13X,5HUPPER,13X,5HLOWER,13X,5H
     2UPPER)
 2041 FORMAT(2X,I3,5E18.8)
 2042 FORMAT(2X,I3,5X,23HPARAMETER HELD CONSTANT)
 2049 FORMAT(1X/29H  NONLINEAR CONFIDENCE LIMITS/15H PHI CRITICAL =E15.8
     1/6H  PARA,6X,8H LOWER B,8X,10H LOWER PHI,10X,8H UPPER B,8X,10H UPP
     2ER PHI)
 2060 FORMAT(57H OUTPUT IS ABBREVIATED DUE TO MATHEMATICAL CONSIDERATION
     1S)
 2090 FORMAT(1H )
 2092 FORMAT(50H CORRECTION VECTOR FOR LAST ITERATION WAS NOT USED)
747     FORMAT(2X,(a))
      END


CNEWA         NEWA - CALCULATES PTP MATRIX, A, AND GRADIENT VECTOR, G.
      SUBROUTINE XNEWAZ(Y,X,RES)

C**********************************************************************
C     Change impliit type of from REAL to REAL(8)                   ***
C**********************************************************************
      IMPLICIT REAL(8) (A-H,O-Z)

      DIMENSION Y(5000),X(5000,3,1000),RES(5000)
      COMMON/BLK1/B(20),P(20),RE,N,M,K
      COMMON/BLK2/A(40,20),SA(20),K2,IK,NPMAX
      COMMON/BLK3/BS(20),DB(20),G(20),K3
      COMMON/BLK4/AL,DELTA,E,FF,GAMCR,T,TAU,ZETA,PHI,SE,PHICR
      COMMON/BLK5/IB(20),IP
      DO1J=1,K
      G(J)=0.
      P(J)=0.
      DO1I=1,K
 1    A(J,I)=0.
      DO 50 II=1,N
C     LOOK FOR PARTIALS
      !J=2 analytic derivatives
      J=2
      CALL MODEL(F,Y,X,RES,II,J)
      RD=RE
      DO30JJ=1,K
C     CHECK FOR OMITTED PARAMETERS
      IF(IP.GT.0)GOTO25
   10 GO TO(20,30,20),J
C     COMPUTE PARTIALS IF NECESSARY
 20   AB=B(JJ)
      B(JJ)=AB+DELTA*AB
      J=1
      CALL MODEL(FDEL,Y,X,RES,II,J)
      RE=RD
      P(JJ)=(FDEL-F)/(DELTA*AB)
      B(JJ)=AB
      GOTO30
   25 DO 26 I=1,IP
      IF(JJ.EQ.IB(I)) GO TO 29
 26   CONTINUE
      GOTO10
  29  P(JJ)=0.
C     USING PARTIALS AT ITH DATA POINT
   30 G(JJ)=G(JJ)+RE*P(JJ)
      DO 40 I=1,K
      DO 40 J=I,K
   40 A(I,J)=A(I,J)+P(I)*P(J)
 50   CONTINUE
      DO 55 I=1,K
      DO 55 J=I,K
   55 A(J,I)=A(I,J)
C         A(I,I)=1.0 FOR OMITTED PARAMETER I
      IF(IP.EQ.0)RETURN
      DO 60 I=1,IP
      DO 60 J=1,K
 60   IF(J.EQ.IB(I))A(J,J)=1.
      RETURN
      END


CSCALE      SUBROUTINE SCALE
      SUBROUTINE XSCALX
C     SCALES ACCORDING TO DIAGONAL ELEMENTS

C**********************************************************************
C     Change impliit type of from REAL to REAL(8)                   ***
C**********************************************************************
      IMPLICIT REAL(8) (A-H,O-Z)

      COMMON/BLK2/A(40,20),SA(20),K2,IK,NPMAX
      K=K2
      DO20I=1,K
      IF(A(I,I).GT.0.) GO TO 15
      SA(I)=0.
      GOTO20
 15   SA(I)=SQRT(A(I,I))
 20   CONTINUE
      DO50I=1,K
      DO40J=1,I
      WS=SA(I)*SA(J)
      IF(WS.GT.0.)GOTO30
      A(I,J)=0.
      GOTO40
 30   A(I,J)=A(I,J)/WS
 40   A(J,I)=A(I,J)
 50   A(I,I)=1.0
      RETURN
      END


CSOLVE      SOLVES (PTP)(DB)=(G)  WHERE PTP IS STORED IN A(I+20,J)
      SUBROUTINE XSOLVX
C     SOLVES  A SET OF LINEAR EQUATIONS IN DB DETERMINED BY MATRIX
C     A AND VECTOR G. USES SUBROUTINE GJR TO INVERT MATRIX

C**********************************************************************
C     Change impliit type of from REAL to REAL(8)                   ***
C**********************************************************************
      IMPLICIT REAL(8) (A-H,O-Z)

      COMMON/BLK2/A(40,20),SA(20),K2,IK,NPMAX
      COMMON/BLK3/BS(20), DB(20), G(20),K3
      COMMON/BLK4/AL,DELTA,E,FF,GAMCR,T,TAU,ZETA,PHI,SE,PHICR
      K=K2
      L=1
C     GET MATRIX FROM STORAGE
 1    DO10I=1,K
      II=I+20
      DO9J=1,K
 9    A(I,J)=A(II,J)
 10   A(I,I)=1.+AL
 20   CALLGJR(MS)
      GOTO(25,100),MS
 25   DO40I=1,K
      DB(I)=0.
      IF(SA(I).LE.0.)GOTO40
      DO30J=1,K
      IF(SA(J).LE.0.)GOTO30
      DB(I)=A(I,J)*G(J)/SA(J) +DB(I)
!      write(*,1000)'i,j,A,G,SA',i,j,A(i,j),G(J),SA(J)
!1000  format((a),2i3,3e10.3)
 30   CONTINUE
      DB(I)=DB(I)/SA(I)
      !debug write(*,1010) 'i,DB(i), SA(i)', i,db(i), sa(i)
      !debug write(7,1010) 'i,DB(i), SA(i)', i,db(i), sa(i)
1010  format((a),i5,3e10.3)
 40   CONTINUE
      RETURN
 100  AL=AL*10.
      L=L+1
      IF(L.GE.6)STOP
      GOTO1
      END

CETEST      SUBROUTINE ETEST
      SUBROUTINE XEPTST(ML)

C**********************************************************************
C     Change impliit type of from REAL to REAL(8)                   ***
C**********************************************************************
      IMPLICIT REAL(8) (A-H,O-Z)

      COMMON/BLK1/B(20),P(20),RE,N,M,K
      COMMON/BLK3/BS(20),DB(20),G(20),K3
      COMMON/BLK4/AL,DELTA,E,FF,GAMCR,T,TAU,ZETA,PHI,SE,PHICR
      EPS=E
      ML=1
      DO 20 I=1,K
      W=ABS(DB(I))/(TAU+ABS(B(I)))
      IF (W.GE.EPS) GO TO 30
   20 CONTINUE
      GO TO 40
   30 ML=2
   40 RETURN
      END


CGJR        GJR - INVERTS A MATRIX IN A(I,J), I=1,20,J=1,20
      SUBROUTINEGJR(MSING)
C     GAUSS-JORDAN-RUTISHAUSER MATRIX INVERSION WITH DOUBLE PIVOTING

C**********************************************************************
C     Change impliit type of from REAL to REAL(8)                   ***
C**********************************************************************
      IMPLICIT REAL(8) (A-H,O-Z)

      COMMON/BLK2/A(40,20),SA(20),K2,IK,NPMAX
      COMMON/BLK4/AL,DELTA,E,FF,GAMCR,T,TAU,ZETA,PHI,SE,PHICR
      DIMENSION P(20),Q(20),B(20),C(20)
      INTEGERP,Q
      EPS=ZETA
      N=K2
      MSING=1
      DO 10 K=1,N
C     DETERMINATION OF PIVOT ELEMENT
      PIVOT=0.
      DO20I=K,N
      DO20J=K,N
      IF(ABS(A(I,J))-ABS(PIVOT))20,20,30
 30   PIVOT=A(I,J)
      P(K)=I
      Q(K)=J
 20   CONTINUE
      IF(ABS(PIVOT)-EPS)40,40,50
C     EXCHANGE OF PIVOTAL ROW WITH KTH ROW
   50 IF (P(K).EQ.K) GO TO 80
      DO70J=1,N
      L=P(K)
      Z=A(L,J)
      A(L,J)=A(K,J)
 70   A(K,J)=Z
C     EXCHANGE OF COLUMN
 80   IF(Q(K).EQ.K)GOTO90
      DO100I=1,N
      L=Q(K)
      Z=A(I,L)
      A(I,L)=A(I,K)
  100 A(I,K)=Z
 90   CONTINUE
C     JORDAN STEP
      DO 110 J=1,N
      IF(J.EQ.K)GOTO120
      B(J)=-A(K,J)/PIVOT
      C(J)=A(J,K)
      GOTO140
 120  B(J)=1./PIVOT
      C(J)=1.
 140  A(K,J)=0.
 110  A(J,K)=0.
      DO 10 I=1,N
      DO 10 J=1,N
   10 A(I,J)=A(I,J)+C(I)*B(J)
C     REORDERING THE MATRIX
      DO155M=1,N
      K=N-M+1
      IF (P(K).EQ.K) GO TO 170
      DO180I=1,N
      L=P(K)
      Z=A(I,L)
      A(I,L)=A(I,K)
 180  A(I,K)=Z
 170  IF(Q(K).EQ.K)GOTO155
      DO150J=1,N
      L=Q(K)
      Z=A(L,J)
      A(L,J)=A(K,J)
 150  A(K,J)=Z
 155  CONTINUE
      RETURN
 40   WRITE(IK,45)P(K),Q(K),PIVOT
 45   FORMAT(20H SINGULAR MATRIX  I=,I3,4H  J=,I3,8H  PIVOT=,E16.8/)
      MSING=2
      RETURN
      END


CPRINT1     SUBROUTINE PRINT1
      SUBROUTINE XPRT1X
C     PRINTS A K BY K SINGLE PRECISION MATRIX

C**********************************************************************
C     Change impliit type of from REAL to REAL(8)                   ***
C**********************************************************************
      IMPLICIT REAL(8) (A-H,O-Z)

      COMMON/BLK2/A(40,20),SA(20),K2,IK,NPMAX
      K=K2
      L=1
 5    JJ=7*L
      LL=JJ-6
      IF(K.LT.LL)GOTO30
      IF(K.LT.JJ)GOTO20
      WRITE(IK,105)LL,JJ
      DO15I=1,K
 15   WRITE(IK,106)(A(I,J),J=LL,JJ)
      L=L+1
      GOTO5
 20   WRITE(IK,105)LL,K
      DO25I=1,K
 25   WRITE(IK,106)(A(I,J),J=LL,K)
 30   RETURN
 105  FORMAT(1X/8H COLUMNS,I4,9H  THROUGH,I4)
 106  FORMAT(1X,7E17.8)
      END


CPRINT2     SUBROUTINE PRINT2
      SUBROUTINE XPRT2X
C     PRINTS A K BY K CORRELATION COEFFICIENT MATRIX

C**********************************************************************
C     Change impliit type of from REAL to REAL(8)                   ***
C**********************************************************************
      IMPLICIT REAL(8) (A-H,O-Z)

      COMMON/BLK2/A(40,20),SA(20),K2,IK,NPMAX
      L=1
      K=K2
    5 JJ=13*L
      LL=JJ-12
      IF(K.LT.LL)GOTO30
      IF(K.LT.JJ)GOTO20
      WRITE(IK,105)LL,JJ
      DO15I=1,K
 15   WRITE(IK,107)(A(I,J),J=LL,JJ)
      L=L+1
      GOTO5
 20   WRITE(IK,105)LL,K
      DO25I=1,K
 25   WRITE(IK,107)(A(I,J),J=LL,K)
 30   RETURN
 105  FORMAT(1X/8H COLUMNS,I4,9H  THROUGH,I4)
  107 FORMAT(1X,13F9.4)
      END


CCONFRG               CONFRG - NON LINEAR CONFIDENCE REGION CALCULATIONS
      SUBROUTINE XCNFXX(Y,X,RES)

C**********************************************************************
C     Change impliit type of from REAL to REAL(8)                   ***
C**********************************************************************
      IMPLICIT REAL(8) (A-H,O-Z)

      COMMON/BLK1/B(20),P(20),RE,N,M,K
      COMMON/BLK2/A(40,20),SA(20),K2,IK,NPMAX
      COMMON/BLK3/BS(20),DB(20),G(20),K3
      COMMON/BLK4/AL,DELTA,E,FF,GAMCR,T,TAU,ZETA,PHI,SE,PHICR
      COMMON/BLK5/IB(20),IP
      DIMENSION Y(5000),X(5000,3,1000),RES(5000)
      LOGICAL NOLO
      DO580J=1,K
      NOLO=.FALSE.
C     CHECK FOR OMITTED PARAMETERS
      IF(IP.EQ.0) GO TO 509
      DO504I=1,IP
      IF(J.EQ.IB(I))GOTO506
 504  CONTINUE
      GOTO509
 506  WRITE(IK,2042)J
      GO TO 580
 509  DDS=-1.
 510  D=DDS
      DJ=SE*SA(J)
      B(J)=BS(J)+D*DJ
      CALLSUMSQ(PH,Y,X,RES)
      IF(PH.LT.PHICR)GOTO530
 520  D=D/2.
      IF(ABS(D).LE..001)GOTO570
      B(J)=BS(J)+D*DJ
      CALLSUMSQ(PPH,Y,X,RES)
      IF(PPH-PHICR)540,540,520
 530  D=D+DDS
      IF(ABS(D).GE.5.0) GO TO 570
      B(J)=BS(J)+D*DJ
      CALLSUMSQ(PPH,Y,X,RES)
      IF(PPH.LT.PHICR)GOTO530
 540  Q=1.-D
      XK1=PHI/D+PH/Q-PPH/(D*Q )
      XK2=-PHI*(1.+D)/D-PH*D/Q+PPH/(D*Q)
      XK3=PHI-PHICR
      BC=(-XK2+SQRT(XK2**2 -4.*XK1*XK3))/(2.*XK1)
      IF(DDS.GT.0.) GO TO 550
      B(J)=BS(J)-BC*DJ
      BL=B(J)
      CALLSUMSQ(PL,Y,X,RES)
 548  DDS=1.
      GOTO510
  550 B(J)=BS(J)+BC*DJ
      BU=B(J)
      CALLSUMSQ(PU,Y,X,RES)
      GOTO576
  570 IF(DDS.GT.0.) GO TO 571
      NOLO=.TRUE.
      GOTO548
 571  IF(NOLO)GOTO575
C     OMITTING UPPER LIMITS
      WRITE(IK,2055)J,BL,PL
      GOTO580
C     OMITTING BOTH
 575  WRITE(IK,2056)J
      GOTO580
 576  IF(NOLO)GOTO578
      WRITE(IK,2052)J,BL,PL,BU,PU
      GOTO580
C     OMITTING LOWER LIMITS
 578  WRITE(IK,2053)J,BU,PU
 580  B(J)=BS(J)
 2042 FORMAT(2X,I3,5X,23HPARAMETER HELD CONSTANT)
 2052 FORMAT(2X,I3,4E18.8)
 2055 FORMAT(2X,I3,2E18.8,11H  NOT FOUND)
 2053 FORMAT(2X,I3,11H  NOT FOUND,25X,2E18.8)
 2056 FORMAT(2X,I3,18X,11H  NOT FOUND)
      RETURN
      END


CMMGR       MICROFILM PLOTTING - RESIDUALS FROM NLLSQ
      SUBROUTINE MMGR(X,RES,N,NDIM,J,GAR)

C**********************************************************************
C     Change impliit type of from REAL to REAL(8)                   ***
C**********************************************************************
      IMPLICIT REAL(8) (A-H,O-Z)

      RETURN
      END


      SUBROUTINE SUMSQ(PHI,Y,X,RES)

C**********************************************************************
C     Change impliit type of from REAL to REAL(8)                   ***
C**********************************************************************
      IMPLICIT REAL(8) (A-H,O-Z)

      DIMENSION Y(5000),X(5000,3,1000),RES(5000)
      COMMON/BLK1/B(20),P(20),RE,N,M,K
      PHI=0.
      DO 10 I=1,N
      CALL MODEL(F,Y,X,RES,I,3)
   10 PHI=PHI+RE*RE
      RETURN
      END
