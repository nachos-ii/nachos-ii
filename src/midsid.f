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
      SUBROUTINE MIDSID (VAR,ICON,LISTEL)
C
C     ******************************************************************
C
C     SUBROUTINE TO INTERPOLATE LINEAR ELEMENT DATA TO THE MIDSIDE NODES
C     OF QUADRATIC ELEMENTS
C
C     ******************************************************************
C
      COMMON /SZDAT/  NUMEL,NUMNOD,NUMVAR,NUMDOF,NODSOL
      COMMON /ELMDAT/ NNELM(6),NNCOR(6),NNMID(6),NNCTR(6),NSIDET(5),    
     1                NSIDEQ(7),NNSIDE(6,4,3)
C
      DIMENSION VAR(*), ICON(NUMEL,*), LISTEL(*)
      DIMENSION VARLOC(9)
C
C     ******************************************************************
C
C     LOOP ON ELEMENTS
C
      DO 70 N=1,NUMEL
      NN=ABS(LISTEL(N))
      KIND=MOD(NN,100)
      NCOR=NNCOR(KIND)
      NMID=NNMID(KIND)
      NCNTR=NNCTR(KIND)
      DO 10 I=1,NCOR
      INODE=ABS(ICON(N,I))
      VARLOC(I)=VAR(INODE)
   10 CONTINUE
C
C     INTERPOLATE VARIABLE TO MIDSIDE NODES
C
      GO TO (20, 20, 30, 30, 40, 40), KIND
C
C     TRIANGULAR ELEMENTS (6 NODE)
C
   20 CONTINUE
      VARLOC(4)=0.5*(VARLOC(1)+VARLOC(2))
      VARLOC(5)=0.5*(VARLOC(2)+VARLOC(3))
      VARLOC(6)=0.5*(VARLOC(3)+VARLOC(1))
      GO TO 50
C
C     QUADRILATERAL ELEMENTS (8 NODE)
C
   30 CONTINUE
      VARLOC(5)=0.5*(VARLOC(1)+VARLOC(2))
      VARLOC(6)=0.5*(VARLOC(2)+VARLOC(3))
      VARLOC(7)=0.5*(VARLOC(3)+VARLOC(4))
      VARLOC(8)=0.5*(VARLOC(4)+VARLOC(1))
      GO TO 50
C
C     QUADRILATERAL ELEMENTS (9 NODE)
C
   40 CONTINUE
      VARLOC(5)=0.5*(VARLOC(1)+VARLOC(2))
      VARLOC(6)=0.5*(VARLOC(2)+VARLOC(3))
      VARLOC(7)=0.5*(VARLOC(3)+VARLOC(4))
      VARLOC(8)=0.5*(VARLOC(4)+VARLOC(1))
      VARLOC(9)=0.25*(VARLOC(1)+VARLOC(2)+VARLOC(3)+VARLOC(4))
C
   50 CONTINUE
      I1=NCOR+1
      I2=NCOR+NMID+NCNTR
      DO 60 I=I1,I2
      INODE=ABS(ICON(N,I))
      VAR(INODE)=VARLOC(I)
   60 CONTINUE
C
   70 CONTINUE
      RETURN
      END
