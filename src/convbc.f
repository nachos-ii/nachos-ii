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
      SUBROUTINE CONVBC (KIND,IN,VALUE,NGAUSS,GS,GSWT,X,Y,KTQ,FTQ)
C
C     ******************************************************************
C
C     SUBROUTINE TO CALCULATE FLUX TERMS FOR A CONVECTIVE BOUNDARY
C     CONDITION (CONSTANT H AND TREF)
C
C     ******************************************************************
C
      REAL KTQ
C
      COMMON /PRBDAT/ IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,    
     1                IPNLTY,PNLTY
      COMMON /ELMDAT/ NNELM(6),NNCOR(6),NNMID(6),NNCTR(6),NSIDET(5),    
     1                NSIDEQ(7),NNSIDE(6,4,3)
      COMMON /BCDAT/  BCSET(2,20,2),IPELM,IPNODE,PBC
C
      DIMENSION X(*), Y(*)
      DIMENSION GS(*), GSWT(*)
      DIMENSION XX(3), YY(3), INODE(3)
      DIMENSION KTQ(9,9), FTQ(9)
C
C     STATEMENT FUNCTION - ONE-DIMENSIONAL SHAPE FUNCTION
C
      FS(I,S)=(S**2-S)*.5*(I-2)*(I-3)*.5-(1.-S**2)*(I-1)*(I-3)+(S**2+S)*
     1.5*(I-1)*(I-2)*.5
C
C     ******************************************************************
C
C     LOCATE DATA FOR ELEMENT SIDE
C
      NN=NNELM(KIND)
      DO 10 J=1,3
      INODE(J)=NNSIDE(KIND,IN,J)
   10 CONTINUE
      DO 20 J=1,3
      I=INODE(J)
      XX(J)=X(I)
      YY(J)=Y(I)
   20 CONTINUE
C
C     COMPUTATION OF FLUX VECTOR
C
      IVALUE=INT(VALUE)
      VALUE=BCSET(1,IVALUE,1)*BCSET(1,IVALUE,2)
      DO 70 I=1,3
      FD=0.
      DO 60 K=1,NGAUSS
      S=GS(K)
      WW=GSWT(K)
      A=FS(I,S)
      GO TO (30, 40, 30, 40, 30, 40), KIND
   30 CONTINUE
      RAD=(1.-S)*XX(1)*.5+(1.+S)*XX(3)*.5
      C=(XX(3)-XX(1))*.5
      D=(YY(3)-YY(1))*.5
      GO TO 50
   40 CONTINUE
      RAD=(S**2-S)*XX(1)*.5+(1.-S**2)*XX(2)+(S**2+S)*XX(3)*.5
      C=(2.*S-1.)*XX(1)*.5+(-2.*S)*XX(2)+(2.*S+1.)*XX(3)*.5
      D=(2.*S-1.)*YY(1)*.5+(-2.*S)*YY(2)+(2.*S+1.)*YY(3)*.5
   50 CONTINUE
      IF (IAXSYM.EQ.0) RAD=1.0
      DELTA=SQRT(C**2+D**2)
      FD=FD+WW*A*RAD*DELTA
   60 CONTINUE
      II=INODE(I)
      FTQ(II)=FTQ(II)+FD*VALUE
   70 CONTINUE
C
C     COMPUTATION OF DIFFUSION CORRECTION MATRIX
C
      DO 120 J=1,3
      DO 120 I=1,3
      FD=0.
      DO 110 K=1,NGAUSS
      S=GS(K)
      WW=GSWT(K)
      A=FS(I,S)
      E=FS(J,S)
      GO TO (80, 90, 80, 90, 80, 90), KIND
   80 CONTINUE
      RAD=(1.-S)*XX(1)*.5+(1.+S)*XX(3)*.5
      C=(XX(3)-XX(1))*.5
      D=(YY(3)-YY(1))*.5
      GO TO 100
   90 CONTINUE
      RAD=(S**2-S)*XX(1)*.5+(1.-S**2)*XX(2)+(S**2+S)*XX(3)*.5
      C=(2.*S-1.)*XX(1)*.5+(-2.*S)*XX(2)+(2.*S+1.)*XX(3)*.5
      D=(2.*S-1.)*YY(1)*.5+(-2.*S)*YY(2)+(2.*S+1.)*YY(3)*.5
  100 CONTINUE
      IF (IAXSYM.EQ.0) RAD=1.0
      DELTA=SQRT(C**2+D**2)
      FD=FD+WW*A*RAD*DELTA*E
  110 CONTINUE
      II=INODE(I)
      JJ=INODE(J)
      KTQ(II,JJ)=KTQ(II,JJ)+FD*BCSET(1,IVALUE,1)
  120 CONTINUE
C
      RETURN
      END
