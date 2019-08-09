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
      SUBROUTINE PREDCT (ITMDEP,KSTEP,IALGOR,DELTN,DELTNM,UN,UPNP1,ACCN,
     1SCRTCH)
C
C     ******************************************************************
C
C     SUBROUTINE TO PREDICT SOLUTION VECTORS USING AN EXPLICIT
C     ADAMS-BASHFORTH METHOD OR AN EXPLICIT FORWARD EULER METHOD
C
C     ******************************************************************
C
      COMMON /SZDAT/ NUMEL,NUMNOD,NUMVAR,NUMDOF,NODSOL
C
      DIMENSION UN(NUMNOD,*), UPNP1(NUMNOD,*), ACCN(NUMNOD,*)
      DIMENSION SCRTCH(NUMNOD,*)
C
C     ******************************************************************
C
      IF (ITMDEP.EQ.0) GO TO 70
      IF (KSTEP.EQ.1) GO TO 10
      IF (IALGOR.EQ.1) GO TO 30
      IF (KSTEP.EQ.2.AND.IALGOR.EQ.2) GO TO 30
      GO TO 50
C
C     FIRST STEP FORMULA (NO PREDICTOR)
C
   10 CONTINUE
      DO 20 I=1,NUMNOD
      DO 20 J=1,NUMVAR
      UPNP1(I,J)=UN(I,J)
   20 CONTINUE
      RETURN
C
C     STANDARD PREDICTION FORMULA - FORWARD EULER METHOD
C
   30 CONTINUE
      DO 40 I=1,NUMNOD
      DO 40 J=1,NUMVAR
      UPNP1(I,J)=UN(I,J)+DELTN*ACCN(I,J)
   40 CONTINUE
      RETURN
C
C     STANDARD PREDICTION FORMULA - ADAMS-BASHFORTH METHOD
C
   50 CONTINUE
      FACT1=0.5*DELTN*(2.+DELTN/DELTNM)
      FACT2=0.5*DELTN*(DELTN/DELTNM)
      DO 60 I=1,NUMNOD
      DO 60 J=1,NUMVAR
      UPNP1(I,J)=UN(I,J)+FACT1*ACCN(I,J)-FACT2*SCRTCH(I,J)
   60 CONTINUE
      RETURN
C
C     STEADY-STATE PREDICTOR
C
   70 CONTINUE
      DO 80 I=1,NUMNOD
      DO 80 J=1,NUMVAR
      UPNP1(I,J)=UN(I,J)
   80 CONTINUE
      RETURN
C
      END
