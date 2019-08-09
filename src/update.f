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
      SUBROUTINE UPDATE (KSTEP,IALGOR,DELTN,UN,ACCN,SCRTCH)
C
C     ******************************************************************
C
C     SUBROUTINE TO COMPUTE ACCELERATION VECTORS AND SHIFT SOLUTION
C     VECTORS FOR NEXT INTEGRATION/ITERATION STEP
C
C     ******************************************************************
C
      COMMON /SZDAT/  NUMEL,NUMNOD,NUMVAR,NUMDOF,NODSOL
      COMMON /PRBDAT/ IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,    
     1                IPNLTY,PNLTY
C
      DIMENSION UN(NUMNOD,*), ACCN(NUMNOD,*), SCRTCH(NUMNOD,*)
C
C     ******************************************************************
C
      IF (ITMDEP.EQ.1) GO TO 20
C
C     STEADY-STATE SOLUTIONS
C
      DO 10 I=1,NUMNOD
      DO 10 J=1,NUMVAR
      UN(I,J)=SCRTCH(I,J)
   10 CONTINUE
      RETURN
C
C     TRANSIENT SOLUTIONS
C
   20 CONTINUE
      IF (IALGOR.EQ.1) GO TO 30
      IF (KSTEP.LE.2.AND.IALGOR.EQ.2) GO TO 30
      GO TO 50
C
C     UPDATE FOR EULER FORMULA
C
   30 CONTINUE
      FACT=1./DELTN
      DO 40 I=1,NUMNOD
      DO 40 J=1,NUMVAR
      ACCEL=FACT*(SCRTCH(I,J)-UN(I,J))
      UN(I,J)=SCRTCH(I,J)
      SCRTCH(I,J)=ACCN(I,J)
      ACCN(I,J)=ACCEL
   40 CONTINUE
      RETURN
C
C     UPDATE TRAPEZOID RULE
C
   50 CONTINUE
      FACT=2./DELTN
      DO 60 I=1,NUMNOD
      DO 60 J=1,NUMVAR
      ACCEL=FACT*(SCRTCH(I,J)-UN(I,J))-ACCN(I,J)
      UN(I,J)=SCRTCH(I,J)
      SCRTCH(I,J)=ACCN(I,J)
      ACCN(I,J)=ACCEL
   60 CONTINUE
C
      RETURN
      END
