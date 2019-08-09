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
      SUBROUTINE FLUXBC (KIND,IN,VALUE,NGAUSS,GS,GSWT,X,Y,F,KINDBC)
C
C     ******************************************************************
C
C     SUBROUTINE TO CALCULATE CONSISTENT NODAL POINT FORCES AND FLUXES
C
C     ******************************************************************
C
      COMMON /PRBDAT/ IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,    
     1                IPNLTY,PNLTY
      COMMON /ELMDAT/ NNELM(6),NNCOR(6),NNMID(6),NNCTR(6),NSIDET(5),    
     1                NSIDEQ(7),NNSIDE(6,4,3)
C
      DIMENSION GS(*), GSWT(*)
      DIMENSION X(*), Y(*), F(*)
      DIMENSION FF(3), FX(3), FY(3)
      DIMENSION XX(3), YY(3), INODE(3)
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
      IF (KINDBC.GT.0) GO TO 110
      GO TO (30, 70, 30, 70, 30, 70), KIND
C
C     NODAL POINT FORCES - SUBPARAMETRIC ELEMENT
C
   30 CONTINUE
      XDIFF=(XX(3)-XX(1))*.5
      YDIFF=(YY(3)-YY(1))*.5
      DO 40 I=1,3
      FF(I)=0.
      DO 40 J=1,NGAUSS
      S=GS(J)
      WW=GSWT(J)
      RAD=(1.-S)*.5*XX(1)+(1.+S)*.5*XX(3)
      IF (IAXSYM.EQ.0) RAD=1.0
      B=FS(I,S)
      FF(I)=FF(I)+WW*RAD*B
   40 CONTINUE
      IF (KINDBC.EQ.-1) THEN
C
C     NORMAL STRESSES
C
      DO 50 J=1,3
      I=INODE(J)
      F(I)=F(I)-FF(J)*YDIFF*VALUE
      F(I+NN)=F(I+NN)+FF(J)*XDIFF*VALUE
   50 CONTINUE
      RETURN
       ELSE
C
C     TANGENTIAL STRESSES
C
      DO 60 J=1,3
      I=INODE(J)
      F(I)=F(I)+FF(J)*XDIFF*VALUE
      F(I+NN)=F(I+NN)-FF(J)*YDIFF*VALUE
   60 CONTINUE
      RETURN
      END IF
C
C     NODAL POINT FORCES - ISOPARAMETRIC ELEMENT
C
   70 CONTINUE
      DO 80 I=1,3
      FX(I)=0.
      FY(I)=0.
      DO 80 J=1,NGAUSS
      S=GS(J)
      WW=GSWT(J)
      RAD=(S**2-S)*.5*XX(1)+(1.-S**2)*XX(2)+(S**2+S)*.5*XX(3)
      IF (IAXSYM.EQ.0) RAD=1.0
      B=FS(I,S)
      C=(2.*S-1.)*.5*XX(1)+(-2.*S)*XX(2)+(2.*S+1.)*.5*XX(3)
      D=(2.*S-1.)*.5*YY(1)+(-2.*S)*YY(2)+(2.*S+1.)*.5*YY(3)
      FX(I)=FX(I)-WW*RAD*B*D
      FY(I)=FY(I)+WW*RAD*B*C
   80 CONTINUE
      IF (KINDBC.EQ.-1) THEN
C
C     NORMAL STRESSES
C
      DO 90 J=1,3
      I=INODE(J)
      F(I)=F(I)+FX(J)*VALUE
      F(I+NN)=F(I+NN)+FY(J)*VALUE
   90 CONTINUE
      RETURN
       ELSE
C
C     TANGENTIAL STRESSES
C
      DO 100 J=1,3
      I=INODE(J)
      F(I)=F(I)+FY(J)*VALUE
      F(I+NN)=F(I+NN)+FX(J)*VALUE
  100 CONTINUE
      RETURN
      END IF
C
  110 CONTINUE
      GO TO (120, 140, 120, 140, 120, 140), KIND
C
C     NODAL POINT FLUXES - SUBPARAMETRIC ELEMENT
C
  120 CONTINUE
      XDIFF=(XX(3)-XX(1))*.5
      YDIFF=(YY(3)-YY(1))*.5
      DELTA=SQRT(XDIFF**2+YDIFF**2)
      IF (IAXSYM.EQ.0) THEN
      FX(1)=DELTA/3.
      FX(2)=DELTA*4./3.
      FX(3)=DELTA/3.
       ELSE
      FX(1)=DELTA*XX(1)/3.
      FX(2)=DELTA*(XX(1)+XX(3))*2./3.
      FX(3)=DELTA*XX(3)/3.
      END IF
      DO 130 J=1,3
      I=INODE(J)
      F(I)=F(I)+VALUE*FX(J)
  130 CONTINUE
      RETURN
C
C     NODAL POINT FLUXES - ISOPARAMETRIC ELEMENT
C
  140 CONTINUE
      DO 150 I=1,3
      FF(I)=0.
      DO 150 J=1,NGAUSS
      S=GS(J)
      WW=GSWT(J)
      RAD=(S**2-S)*.5*XX(1)+(1.-S**2)*XX(2)+(S**2+S)*.5*XX(3)
      IF (IAXSYM.EQ.0) RAD=1.0
      B=FS(I,S)
      C=(2.*S-1.)*.5*XX(1)+(-2.*S)*XX(2)+(2.*S+1.)*.5*XX(3)
      D=(2.*S-1.)*.5*YY(1)+(-2.*S)*YY(2)+(2.*S+1.)*.5*YY(3)
      DELTA=SQRT(C**2+D**2)
      FF(I)=FF(I)+WW*RAD*B*DELTA
  150 CONTINUE
      DO 160 J=1,3
      I=INODE(J)
      F(I)=F(I)+VALUE*FF(J)
  160 CONTINUE
C
      RETURN
      END
