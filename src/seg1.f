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
      SUBROUTINE SEG1 (X1,Y1,X2,Y2,GRAD,XS,YS,N)
C
C     ******************************************************************
C
C     SUBROUTINE TO DIVIDE A LINE WITH END POINTS X1,Y1 AND X2,Y2
C     INTO N-1 SEGMENTS WITH A GRADIENT=GRAD
C
C     ******************************************************************
C
      REAL LLS,LSS,L1,LF,NSEG
C
      DIMENSION XS(*), YS(*)
C
C     ******************************************************************
C
      IF (N.EQ.1) GO TO 30
      DEL=0.0
      NSEG=N-1
      NS=NSEG
      LF=SQRT((X2-X1)**2+(Y2-Y1)**2)
      LLS=(2.*LF/NSEG)*GRAD/(GRAD+1)
      LSS=LLS/GRAD
      IF (NSEG.EQ.1.) GO TO 10
      DEL=(LLS-LSS)/(NSEG-1.)
   10 SUM=0.0
      DO 20 I=1,NS
      RJ=I-1
      L1=LLS-RJ*DEL
      XS(I)=X1+(X2-X1)*SUM/LF
      YS(I)=Y1+(Y2-Y1)*SUM/LF
   20 SUM=SUM+L1
      XS(N)=X2
      YS(N)=Y2
      RETURN
   30 XS(1)=X1
      YS(1)=Y1
C
      RETURN
      END
