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
      SUBROUTINE MAXMIN (LISTEL,ICON,IVAR,VAR,IWHERE,VMAX,VMIN,IELMX,   
     1INODMX,IELMN,INODMN)
C
C     ***************************************************************
C
C     SUBROUTINE TO FIND THE MAX AND MIN VALUES IN A SOLUTION VECTOR
C
C     ***************************************************************
C
      COMMON /SZDAT/  NUMEL,NUMNOD,NUMVAR,NUMDOF,NODSOL
      COMMON /ELMDAT/ NNELM(6),NNCOR(6),NNMID(6),NNCTR(6),NSIDET(5),    
     1                NSIDEQ(7),NNSIDE(6,4,3)
C
      DIMENSION LISTEL(*), ICON(NUMEL,*), VAR(*)
C
C     ***************************************************************
C
      VMAX=0.
      VMIN=1.0E20
      DO 20 N=1,NUMEL
      NN=ABS(LISTEL(N))
      KIND=MOD(NN,100)
      IF (IWHERE.EQ.1) THEN
      NPTS=NNCOR(KIND)
       ELSE
      NPTS=NNELM(KIND)
      END IF
      DO 10 I=1,NPTS
      INODE=ABS(ICON(N,I))
      IF (VAR(INODE).GT.VMAX) THEN
      VMAX=VAR(INODE)
      IELMX=N
      INODMX=I
      GO TO 10
      END IF
      IF (VAR(INODE).LT.VMIN) THEN
      VMIN=VAR(INODE)
      IELMN=N
      INODMN=I
      END IF
   10 CONTINUE
   20 CONTINUE
C
      RETURN
      END
