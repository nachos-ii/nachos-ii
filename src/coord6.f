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
      SUBROUTINE COORD6 (XB,YB,XX,YY,S,T,INODE,ICONV)
C
C     ******************************************************************
C
C     SUBROUTINE TO LOCATE NORMALIZED COORDINATES S,T IN PARENT ELEMENT
C     ISOPARAMETRIC SIX-NODE TRIANGLE
C
C     ******************************************************************
C
      DIMENSION XX(*), YY(*)
      DIMENSION SS(6), TT(6)
C
      DATA SS/1.,0.,0.,.5,0.,.5/
      DATA TT/0.,1.,0.,.5,.5,0./
      DATA TOL/.005/,IMAX/20/
C
C     STATEMENT FUNCTIONS FOR SHAPE FUNCTIONS AND DERIVATIVES
C
      FX(SP,TP)=XX(1)*SP*(2.*SP-1.)+XX(2)*TP*(2.*TP-1.)+XX(3)*(2.*(SP+TP
     1)**2-3.*(SP+TP)+1.)+XX(4)*4.*SP*TP+XX(5)*4.*(-TP*TP-SP*TP+TP)+XX(6
     2)*4.*(-SP*SP-SP*TP+SP)
C
      FY(SP,TP)=YY(1)*SP*(2.*SP-1.)+YY(2)*TP*(2.*TP-1.)+YY(3)*(2.*(SP+TP
     1)**2-3.*(SP+TP)+1.)+YY(4)*4.*SP*TP+YY(5)*4.*(-TP*TP-SP*TP+TP)+YY(6
     2)*4.*(-SP*SP-SP*TP+SP)
C
      DFXDS(SP,TP)=XX(1)*(4.*SP-1.)+XX(3)*(4.*SP+4.*TP-3.)+XX(4)*(4.*TP)
     1+XX(5)*(-4.*TP)+XX(6)*(-8.*SP-4.*TP+4.)
C
      DFYDS(SP,TP)=YY(1)*(4.*SP-1.)+YY(3)*(4.*SP+4.*TP-3.)+YY(4)*(4.*TP)
     1+YY(5)*(-4.*TP)+YY(6)*(-8.*SP-4.*TP+4.)
C
      DFXDT(SP,TP)=XX(2)*(4.*TP-1.)+XX(3)*(4.*TP+4.*SP-3.)+XX(4)*(4.*SP)
     1+XX(5)*(-8.*TP-4.*SP+4.)+XX(6)*(-4.*SP)
C
      DFYDT(SP,TP)=YY(2)*(4.*TP-1.)+YY(3)*(4.*TP+4.*SP-3.)+YY(4)*(4.*SP)
     1+YY(5)*(-8.*TP-4.*SP+4.)+YY(6)*(-4.*SP)
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
