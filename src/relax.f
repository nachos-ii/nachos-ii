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
      SUBROUTINE RELAX (KSTEP,ACCF,UN,SCRTCH)
C
C     ******************************************************************
C
C     SUBROUTINE TO APPLY AN UNDER-RELAXATION FACTOR TO A
C     SOLUTION VECTOR
C
C     ******************************************************************
C
      COMMON /SZDAT/ NUMEL,NUMNOD,NUMVAR,NUMDOF,NODSOL
C
      DIMENSION UN(NUMNOD,*), SCRTCH(NUMNOD,*)
C
C     ******************************************************************
C
      IF (ACCF.EQ.0.) RETURN
      IF (KSTEP.EQ.1) THEN
      FACT1=0.0
      FACT2=1.0
       ELSE
      FACT1=ACCF
      FACT2=1.0-ACCF
      END IF
C
      DO 10 I=1,NUMNOD
      DO 10 J=1,NUMVAR
      SCRTCH(I,J)=FACT1*UN(I,J)+FACT2*SCRTCH(I,J)
   10 CONTINUE
C
      RETURN
      END
