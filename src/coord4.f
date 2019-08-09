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
      SUBROUTINE COORD4 (XB,YB,XX,YY,S,T,INODE,ICONV)
C
C     ******************************************************************
C
C     SUBROUTINE TO LOCATE NORMALIZED COORDINATES S,T IN PARENT ELEMENT
C     SUBPARAMETRIC EIGHT- AND NINE-NODE QUADRILATERALS
C
C     ******************************************************************
C
      DIMENSION XX(*), YY(*)
      DIMENSION SS(4), TT(4)
C
      DATA SS/-1.,1.,1.,-1./
      DATA TT/-1.,-1.,1.,1./
      DATA TOL/.005/,IMAX/20/
C
C     STATEMENT FUNCTIONS FOR SHAPE FUNCTIONS AND DERIVATIVES
C
      FX(SP,TP)=.25*XX(1)*(1.-SP)*(1.-TP)+.25*XX(2)*(1.+SP)*(1.-TP)+.25*
     1XX(3)*(1.+SP)*(1.+TP)+.25*XX(4)*(1.-SP)*(1.+TP)
C
      FY(SP,TP)=.25*YY(1)*(1.-SP)*(1.-TP)+.25*YY(2)*(1.+SP)*(1.-TP)+.25*
     1YY(3)*(1.+SP)*(1.+TP)+.25*YY(4)*(1.-SP)*(1.+TP)
C
      DFXDS(SP,TP)=.25*XX(1)*(-1.)*(1.-TP)+.25*XX(2)*(1.-TP)+.25*XX(3)*(
     11.+TP)+.25*XX(4)*(-1.)*(1.+TP)
C
      DFYDS(SP,TP)=.25*YY(1)*(-1.)*(1.-TP)+.25*YY(2)*(1.-TP)+.25*YY(3)*(
     11.+TP)+.25*YY(4)*(-1.)*(1.+TP)
C
      DFXDT(SP,TP)=.25*XX(1)*(-1.)*(1.-SP)+.25*XX(2)*(-1.)*(1.+SP)+.25*X
     1X(3)*(1.+SP)+.25*XX(4)*(1.-SP)
C
      DFYDT(SP,TP)=.25*YY(1)*(-1.)*(1.-SP)+.25*YY(2)*(-1.)*(1.+SP)+.25*Y
     1Y(3)*(1.+SP)+.25*YY(4)*(1.-SP)
C
C     ******************************************************************
C
C     CHECK FOR COINCIDENT NODE
C
      IF (INODE.EQ.0) GO TO 10
C
      S=SS(INODE)
      T=TT(INODE)
      RETURN
C
C     USE NEWTON RAPHSON PROCEDURE TO COMPUTE S,T COORDINATES
C
   10 CONTINUE
      SG=0.
      TG=0.
      ICOUNT=0
   20 CONTINUE
C
C     BUILD JACOBIAN AND INVERT
C
      A11=DFXDS(SG,TG)
      A12=DFYDS(SG,TG)
      A21=DFXDT(SG,TG)
      A22=DFYDT(SG,TG)
      F1=FX(SG,TG)-XB
      F2=FY(SG,TG)-YB
      DETA=A11*A22-A12*A21
      SNEW=SG-(A22*F1-A21*F2)/DETA
      TNEW=TG-(-A12*F1+A11*F2)/DETA
C
C     CHECK CONVERGENCE
C
      ICOUNT=ICOUNT+1
      DS=ABS(SNEW-SG)
      DT=ABS(TNEW-TG)
      IF (DS.LT.TOL.AND.DT.LT.TOL) GO TO 30
      SG=SNEW
      TG=TNEW
      IF (ICOUNT.EQ.IMAX) GO TO 40
      GO TO 20
C
   30 CONTINUE
      S=SNEW
      T=TNEW
      RETURN
C
   40 CONTINUE
      ICONV=1
C
      RETURN
      END
