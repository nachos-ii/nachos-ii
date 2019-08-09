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
      SUBROUTINE SPOINT (ITYPE,LISTEL,ICON,UN)
C
C     ******************************************************************
C
C     SUBROUTINE TO EVALUATE THE SOLUTION AT THE SPECIAL OUTPUT POINTS
C
C     ******************************************************************
C
      COMMON /SZDAT/  NUMEL,NUMNOD,NUMVAR,NUMDOF,NODSOL
      COMMON /PRBDAT/ IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,    
     1                IPNLTY,PNLTY
      COMMON /ELMDAT/ NNELM(6),NNCOR(6),NNMID(6),NNCTR(6),NSIDET(5),    
     1                NSIDEQ(7),NNSIDE(6,4,3)
      COMMON /SPTDAT/ ISPT,NSPT,NELSPT(50),PTS(50,2),SPTVAL(50,6)
C
      DIMENSION ICON(NUMEL,*), LISTEL(*)
      DIMENSION UN(NUMNOD,*)
      DIMENSION UNELM(9,6)
C
C     ******************************************************************
C
C     SET UP POINTERS FOR SPECIAL POINTS AND ELEMENTS
C
      DO 30 I=1,NSPT
      NEL=NELSPT(I)
      SPT=PTS(I,1)
      TPT=PTS(I,2)
      NN=ABS(LISTEL(NEL))
      KIND=MOD(NN,100)
      NN=NNELM(KIND)
      NC=NNCOR(KIND)
      DO 10 J=1,NN
      INODE=ABS(ICON(NEL,J))
      DO 10 K=1,NUMVAR
      UNELM(J,K)=UN(INODE,K)
   10 CONTINUE
C
C     EVALUATE SOLUTION AT SPECIAL POINTS
C
      DO 20 K=1,NUMVAR
      CALL PTFUNC (KIND,K,SPT,TPT,UNELM(1,K),VALUE)
      SPTVAL(I,K)=VALUE
   20 CONTINUE
   30 CONTINUE
C
      RETURN
      END
