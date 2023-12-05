! provided by Xiaohua YANG, DMI, 2012-09-06
      SUBROUTINE RAYMOND (F,NLON,NLAT,NLEV,EPS)
      IMPLICIT NONE
      REAL*8 F(NLON*NLAT,NLEV)
      INTEGER NLON,NLAT,NLEV
      REAL*8 EPS
CF2PY INTENT(IN,OUT,COPY) F
      CALL IMPFILA(F,NLON,NLAT,NLEV,EPS,NLON,1)
      CALL IMPFILA(F,NLON,NLAT,NLEV,EPS,NLAT,NLON)
      END

      SUBROUTINE IMPFILA(F,NLON,NLAT,NLEV,EPS,NX,NSTEP)
      IMPLICIT NONE
      INTEGER K,NSTEP,NLON,NLAT,NLEV,NX
      REAL*8 F(NLON*NLAT,NLEV)
      REAL*8 EPS
      REAL*8 WRK(NX,13)
CF2PY INTENT(IN,OUT,COPY) F
c     DO I=1,NX*13
c      WRK(I,1)=0.
c     ENDDO
      WRK=0.
      CALL INVLOW_V(1,NX,F,F,F,EPS,.TRUE.,
     + wrk(1,1),
     + wrk(1,2),
     + wrk(1,3),
     + wrk(1,4),
     + wrk(1,5),
     + wrk(1,6),
     + wrk(1,7),
     + wrk(1,8),
     + wrk(1,9),
     + wrk(1,10),
     + wrk(1,11),
     + wrk(1,12),
     + wrk(1,13)
     + )
C
CDIR_PARALLEL_LOOP
C
      DO K=1,NLEV
      CALL FILSUB(NLON,NLAT,NSTEP,F(1,K),EPS,WRK,NX)
      ENDDO
      END

      SUBROUTINE FILSUB(NLON,NLAT,NSTEP,XY,EPS,WRK,NX)
C
C       SIXTH-ORDER LOW-PASS IMPLICIT TANGENT FILTER
C       (RAYMOND, MWR, 116, 2132-2141)
C
C*************************************************************
C***     THIS CODE IS COPIED FROM A LISTING PROVIDED       ***
C***     BY WILLIAM H RAYMOND. SOME NOTATIONAL CHANGES     ***
C***     HAVE BEEN MADE IN THE ROUTINE LOWPAS. THE         ***
C***     ROUTINE INVLOW HAS BEEN COPIED ALMOST VERBATIM.   ***
C*************************************************************
C
C       XY     UNFILTERED VALUES ON INPUT
C              FILTERED VALUES ON OUTPUT.
C       NLON AND NLAT   NUMBER OF VALUES.
C       EPS    FILTER PARAMETER
C              (DETERMINES CUTOFF)
C       NSTEP  STEP LENGTH (FOR MULTIDIMENSIONAL ARRAYS)
C
C---------------------------------------------------------------
C
      IMPLICIT NONE
      REAL*8 EPS
      INTEGER NLON,NLAT,NX,I,J,NSTEP,JSX,JSY,NV,NF,JSF,JSV,K
      REAL*8 RHS,XDASH,H,WRK,XY
      DIMENSION XY(NLON,NLAT)
      DIMENSION RHS(NLON*NLAT),XDASH(NLON*NLAT),H(NLON*NLAT)
      DIMENSION WRK(NX,*)
CF2PY INTENT(IN,OUT,COPY) XY
C
C---------------------------------------------------------------
C
C  NV  number of points in the vector
C  NF  number of points in the filter-direction
C  JSX the distance in words to the neighbor point in X-direction
C  JSX the distance in words to the neighbor point in Y-direction
C
      IF (NSTEP.EQ.1) THEN
      JSX=NLAT
      JSY=1
      NV=NLAT
      NF=NLON
      ELSEIF ( NSTEP.EQ.NLON) THEN
      JSX=1
      JSY=NLON
      NV=NLON
      NF=NLAT
      JSF=1
      JSV=NLAT
      ELSE
      JSX=1
      JSY=1
      NV=1
      NF=NLON
      ENDIF
C
C     DEFINE RHS
C
      CALL RHSINI(NV,NF,NSTEP,XY,RHS,EPS)
C
C      SOLVE FOR XDASH
C
      CALL INVLOW_V(NV,NF,RHS,XDASH,H,EPS,.FALSE.,
     + WRK(1,1),
     + WRK(1,2),
     + WRK(1,3),
     + WRK(1,4),
     + WRK(1,5),
     + WRK(1,6),
     + WRK(1,7),
     + WRK(1,8),
     + WRK(1,9),
     + WRK(1,10),
     + WRK(1,11),
     + WRK(1,12),
     + WRK(1,13)
     + )
C
C       ADD CORRECTION TO GET FILTERED VALUES.
C
      DO 30 J=2,NLAT-2
      DO 30 I=2,NLON-2
      K=1+(I-1)*JSX+(J-1)*JSY
      XY(I,J) = XY(I,J) + XDASH(K)
   30 CONTINUE
      RETURN
      END

      SUBROUTINE INVLOW_V(NP,N,BB,XANS,H,EP,INIT,
     + A,B,C,D,E,DELTA,BETA,W,GAM,AP,F,Z,PI
     + )
C
C       GAUSSIAN ELIMINATION FOR LOW-PASS FILTER.
C
C       SIXTH-ORDER LOW-PASS IMPLICIT TANGENT FILTER.
C       (REF: WILLIAM H RAYMOND, MWR, 116, 2132-2124)
C
      IMPLICIT NONE
      LOGICAL INIT
      REAL*8 H,PI,GAM,BB,BETA,DELTA,W,EP,XANS
      INTEGER I,JL,N,NP
      REAL*8 A,B,C,D,E,F,Z
      REAL*8 NSAVE,AP
C
C     PARAMETER(NNMAX=200)
C
      DIMENSION A(N),B(N),C(N),D(N),E(N),
     +            DELTA(N),BETA(N),W(N),GAM(N),
     +            AP(N),F(N),Z(N),PI(N)
      SAVE NSAVE
c     SAVE A,B,C,D,E,DELTA,BETA,W,GAM,PI,AP,F,Z,NSAVE
      DIMENSION  H(NP,N),XANS(NP,N),BB(NP,N)
      DATA NSAVE /0/
C
C---------------------------------------------------------------
C
C       SKIP INITIALIZATION OF MATRIX ON REPEAT CALLS.
      IF ( .NOT.INIT ) GO TO 100
      NSAVE = N
C
C       INITIALIZE THE MATRIX
C
      DO 10 I=4,N-3
      Z(I) = 1-EP
      A(I) = 6*(1+EP)
      B(I) = 15*(1-EP)
      C(I) = 20*(1+EP)
      D(I) = B(I)
      E(I) = A(I)
      F(I) = Z(I)
   10 CONTINUE
C
      Z(1) = 0
      Z(2) = 0
      Z(3) = 0
C
      A(1) = 0
      A(2) = 0
      A(3) = 1+EP
C
      B(1) = 0
      B(2) = 1-EP
      B(3) = 4*(1-EP)
C
      C(1) = 1
      C(2) = 2*(1+EP)
      C(3) = 6*(1+EP)
C
      D(1) = 0
      D(2) = 1-EP
      D(3) = 4*(1-EP)
C
      E(1) = 0
      E(2) = 0
      E(3) = 1+EP
C
      F(1) = 0
      F(2) = 0
      F(3) = 0
C
C
      Z(N-2) = 0
      Z(N-1) = 0
      Z(N) = 0
C
      A(N-2) = 1+EP
      A(N-1) = 0
      A(N) = 0
C
      B(N-2) = 4*(1-EP)
      B(N-1) = 1-EP
      B(N) = 0
C
      C(N-2) = 6*(1+EP)
      C(N-1) = 2*(1+EP)
      C(N) = 1
C
      D(N-2) = 4*(1-EP)
      D(N-1) = 1-EP
      D(N) = 0
C
      E(N-2) = 1+EP
      E(N-1) = 0
      E(N) = 0
C
      F(N-2) = 0
      F(N-1) = 0
      F(N) = 0
C
C       Step One.
C
      BETA(1) = D(1)/C(1)
      DELTA(2) = B(2)
      W(1) = 1./C(1)
      PI(1) = F(1)*W(1)
      AP(1) = 0
      AP(2) = 0
      AP(3) = A(3)
      W(2) = 1./(C(2)-DELTA(2)*BETA(1))
      GAM(1) = E(1)/C(1)
      BETA(2) = (D(2)-DELTA(2)*GAM(1))*W(2)
      GAM(2) = (E(2)-PI(1)*DELTA(2))*W(2)
      PI(2) = F(2)*W(2)
      DELTA(3) = (B(3)-AP(3)*BETA(1))
      W(3) =1./(C(3)-DELTA(3)*BETA(2)-AP(3)*GAM(1))
      BETA(3) = (D(3)-AP(3)*PI(1)-DELTA(3)*GAM(2))*W(3)
      GAM(3) = (E(3)-DELTA(3)*PI(2))*W(3)
      PI(3) = F(3)*W(3)
C
C       Step Two
      DO 20 I=4,N
      AP(I) = A(I)-Z(I)*BETA(I-3)
      DELTA(I) = B(I)-AP(I)*BETA(I-2)-Z(I)*GAM(I-3)
      W(I) = 1./(C(I)-AP(I)*GAM(I-2)-DELTA(I)*BETA(I-1)
     +           -Z(I)*PI(I-3))
      BETA(I) = (D(I)-AP(I)*PI(I-2)-DELTA(I)*GAM(I-1))*W(I)
      GAM(I) = (E(I)-DELTA(I)*PI(I-1))*W(I)
      PI(I) = F(I)*W(I)
   20 CONTINUE
C
      RETURN
C
  100 CONTINUE
C
C
C       Step Three
C
      DO 25 JL=1,NP
      H(JL,1) = BB(JL,1)*W(1)
      H(JL,2) = (BB(JL,2)-DELTA(2)*H(JL,1))*W(2)
      H(JL,3) = (BB(JL,3)-DELTA(3)*H(JL,2)-AP(3)*H(JL,1))*W(3)
   25 CONTINUE
      DO 30 I=4,N
      DO 30 JL=1,NP
      H(JL,I) = (BB(JL,I)-DELTA(I)*H(JL,I-1)-AP(I)*H(JL,I-2)
     +             -Z(I)*H(JL,I-3))*W(I)
   30 CONTINUE
C
C       Step Four
C
      DO 35 JL=1,NP
      XANS(JL,N) = H(JL,N)
      XANS(JL,N-1)=H(JL,N-1)-BETA(N-1)*XANS(JL,N)
      XANS(JL,N-2)=H(JL,N-2)-BETA(N-2)*XANS(JL,N-1)-GAM(N-2)*XANS(JL,N)
   35 CONTINUE
      DO 40 I=N-3,1,-1
      DO 40 JL=1,NP
      XANS(JL,I) = H(JL,I)-BETA(I)*XANS(JL,I+1)-GAM(I)*XANS(JL,I+2)
     +                                 -PI(I)*XANS(JL,I+3)
   40 CONTINUE
      RETURN
      END
      SUBROUTINE RHSINI(NHOR,N,ISTEP,XY,RHS,ZEPS)
      IMPLICIT NONE
      INTEGER ISTEP,NR1,NR2,NR3,NHOR,N,I,IL,JL,JS,N1,N2,N3,N4,N5
      INTEGER NM0,NM1,NM2,NM3,NM4
      REAL*8 RHS,XY,ZEPS
      DIMENSION RHS(NHOR,N),XY(*)
CF2PY INTENT(IN,OUT,COPY) RHS
C
C   RELATIVE POINTS
C
      NR1=1*ISTEP
      NR2=2*ISTEP
      NR3=3*ISTEP
C
C   4 LAST POINTS
C
      NM0 = (N-1)*ISTEP
      NM1 = (N-2)*ISTEP
      NM2 = (N-3)*ISTEP
      NM3 = (N-4)*ISTEP
      NM4 = (N-5)*ISTEP
C
C   5 FIRST POINTS
C
      N1 = 0*ISTEP
      N2 = 1*ISTEP
      N3 = 2*ISTEP
      N4 = 3*ISTEP
      N5 = 4*ISTEP
      if ( ISTEP.EQ.1) THEN
      JS=N
      else
      JS=1
      endif
C
      DO 10 JL=1,NHOR
      IL=(JL-1)*JS+1
      RHS(JL,1) = 0.
      RHS(JL,N) = 0.
      RHS(JL,2) = ZEPS*(XY(IL+  N1)-2*XY(IL+  N2)+XY(IL+  N3))
      RHS(JL,N-1) = ZEPS*(XY(IL+ NM2)-2*XY(IL+NM1)+XY(IL+  NM0))
      RHS(JL,  3) = ZEPS*(-1*(XY(IL+  N1)+XY(IL+  N5))
     +                   +4*(XY(IL+  N2)+XY(IL+  N4))
     +                   -6* XY(IL+  N3)         )
      RHS(JL,N-2) = ZEPS*(-1*(XY(IL+NM0)+XY(IL+NM4))
     +                    +4*(XY(IL+NM1)+XY(IL+NM3))
     +                    -6* XY(IL+NM2)         )
  10  CONTINUE
      DO 20 I=4,N-3
      DO 30 JL=1,NHOR
      IL=(I-1)*ISTEP+(JL-1)*JS+1
      RHS(JL,I) = ZEPS*( (XY(IL-NR3)+XY(IL+NR3))
     +                 - 6*(XY(IL-NR2)+XY(IL+NR2))
     +                 +15*(XY(IL-NR1)+XY(IL+NR1))
     +                 -20* XY(IL)       )
  30  CONTINUE
  20  CONTINUE
      END

