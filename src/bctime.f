C
C***********************************************************************
C                                                                      *
C            NACHOS II - A Finite Element Computer Program             *
C                        for Incompressible Flow Problems              *
C                                                                      *
C     Copyright (c) 1986,2019   National Technology & Engineering      *
C                               Solutions of Sandia, LLC (NTESS)       *
C                                                                      *
C                            All rights reserved.                      *
C                                                                      *
C     This software is distributed under the BSD 3-Clause License.     *
C                                                                      *
C***********************************************************************
C
      SUBROUTINE BCTIME (KNUM,IVAR,ICURVE,ISIDE,ITIME,TIME,KSTEP,KIND,  
     1X,Y,U,V,P,T,V1,V2,F)
C
C     ******************************************************************
C
C     SUBROUTINE TO EVALUATE RHS VECTOR AS A FUNCTION OF TIME
C
C     ******************************************************************
C
      COMMON /ELMDAT/ NNELM(6),NNCOR(6),NNMID(6),NNCTR(6),NSIDET(5),    
     1                NSIDEQ(7),NNSIDE(6,4,3)
C
      DIMENSION F(*), X(*), Y(*), U(*), V(*), P(*), T(*), V1(*), V2(*)
      DIMENSION XSURF(3), YSURF(3), USURF(3), VSURF(3), PSURF(3)
      DIMENSION TSURF(3), V1SURF(3), V2SURF(3), VALUE(3)
C
C     ******************************************************************
C
C     COLLECT DATA FOR APPROPRIATE ELEMENT SIDE
C
      IF (ICURVE.GT.6) CALL ERROR ('BCTIME','TIME CURVE NUMBER TOO LARGE
     1','CURVE NUMBER',ICURVE,'MAX CURVE NUMBER ALLOWED',6,' ',' ',1)
      DO 10 J=1,3
      KK=NNSIDE(KIND,ISIDE,J)
      XSURF(J)=X(KK)
      YSURF(J)=Y(KK)
      USURF(J)=U(KK)
      VSURF(J)=V(KK)
      TSURF(J)=T(KK)
      V1SURF(J)=V1(KK)
      V2SURF(J)=V2(KK)
      IF (J.EQ.2) GO TO 10
      PSURF(J)=P(KK)
   10 CONTINUE
      PSURF(2)=0.5*(PSURF(1)+PSURF(3))
C
C     FIND VALUE FOR TIME DEPENDENT VARAIABLE
C
      GO TO (20, 30, 40, 50, 60, 70), ICURVE
   20 CONTINUE
      CALL USRBC1 (VALUE,USURF,VSURF,PSURF,TSURF,V1SURF,V2SURF,XSURF,   
     1YSURF,3,KNUM,TIME,KSTEP)
      GO TO 80
   30 CONTINUE
      CALL USRBC2 (VALUE,USURF,VSURF,PSURF,TSURF,V1SURF,V2SURF,XSURF,   
     1YSURF,3,KNUM,TIME,KSTEP)
      GO TO 80
   40 CONTINUE
      CALL USRBC3 (VALUE,USURF,VSURF,PSURF,TSURF,V1SURF,V2SURF,XSURF,   
     1YSURF,3,KNUM,TIME,KSTEP)
      GO TO 80
   50 CONTINUE
      CALL USRBC4 (VALUE,USURF,VSURF,PSURF,TSURF,V1SURF,V2SURF,XSURF,   
     1YSURF,3,KNUM,TIME,KSTEP)
      GO TO 80
   60 CONTINUE
      CALL USRBC5 (VALUE,USURF,VSURF,PSURF,TSURF,V1SURF,V2SURF,XSURF,   
     1YSURF,3,KNUM,TIME,KSTEP)
      GO TO 80
   70 CONTINUE
      CALL USRBC6 (VALUE,USURF,VSURF,PSURF,TSURF,V1SURF,V2SURF,XSURF,   
     1YSURF,3,KNUM,TIME,KSTEP)
   80 CONTINUE
C
C     UPDATE RHS VECTOR
C
      IF (ITIME.GT.1) GO TO 100
C
C     TIME DEPENDENT VARIABLE
C
      IADD=0
      IF (IVAR.EQ.2) IADD=NNELM(KIND)
      DO 90 J=1,3
      KK=NNSIDE(KIND,ISIDE,J)+IADD
      F(KK)=F(KK)*VALUE(J)
   90 CONTINUE
      RETURN
C
C     TIME DEPENDENT FLUX
C
  100 CONTINUE
      DO 120 J=1,3
      KK=NNSIDE(KIND,ISIDE,J)
      IF (ABS(F(KK)).GT.1.0E10) GO TO 110
      F(KK)=F(KK)*VALUE(J)
  110 CONTINUE
      IF (ITIME.EQ.2) GO TO 120
      KK=KK+NNELM(KIND)
      IF (ABS(F(KK)).GT.1.0E10) GO TO 120
      F(KK)=F(KK)*VALUE(J)
  120 CONTINUE
C
      RETURN
      END
