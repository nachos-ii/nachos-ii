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
      SUBROUTINE SRCHQ (NN,XX,YY,XB,YB,S,T,IFIND)
C
C     ******************************************************************
C
C     SUBROUTINE TO FIND THE LOCATION OF A POINT IN A MESH
C     ISOPARAMETRIC QUADRILATERALS AND TRIANGLES
C
C     ******************************************************************
C
      DIMENSION XX(*), YY(*), AL(9)
C
      DATA EPS/.01/,ETA/.20/,STEST/1.001/
C
C     ******************************************************************
C
      IFIND=0
      AL(7)=0.
      AL(8)=0.
      AL(9)=0.
C
C     CHECK FOR ELEMENT IN VICINITY OF NODE
C
      IF (NN.EQ.6) GO TO 10
      XC=-(XX(1)+XX(2)+XX(3)+XX(4))/4.+(XX(5)+XX(6)+XX(7)+XX(8))/2.
      YC=-(YY(1)+YY(2)+YY(3)+YY(4))/4.+(YY(5)+YY(6)+YY(7)+YY(8))/2.
      GO TO 20
   10 CONTINUE
      XC=-(XX(1)+XX(2)+XX(3))/9.+4.*(XX(4)+XX(5)+XX(6))/9.
      YC=-(YY(1)+YY(2)+YY(3))/9.+4.*(YY(4)+YY(5)+YY(6))/9.
   20 CONTINUE
      DO 30 I=1,NN
      A=XX(I)-XC
      B=YY(I)-YC
      AL(I)=SQRT(A**2+B**2)
   30 CONTINUE
      ALM=MAX(AL(1),AL(2),AL(3),AL(4),AL(5),AL(6),AL(7),AL(8))
      ALTEST=ALM*(1.+ETA)
      A=XB-XC
      B=YB-YC
      ALB=SQRT(A**2+B**2)
      IF (ALB.GT.ALTEST) RETURN
C
C     CHECK FOR COINCIDENT NODES
C
      INODE=0
      ALTEST=ALM*EPS
      DO 40 I=1,NN
      A=XX(I)-XB
      B=YY(I)-YB
      ALB=SQRT(A**2+B**2)
      INODE=I
      IF (ALB.LT.ALTEST) GO TO 50
   40 CONTINUE
      INODE=0
   50 CONTINUE
C
C     FIND NORMALIZED COORDINATES IN PARENT ELEMENT
C
      ICONV=0
      IF (NN.EQ.6) CALL COORD6 (XB,YB,XX,YY,S,T,INODE,ICONV)
      IF (NN.EQ.8) CALL COORD8 (XB,YB,XX,YY,S,T,INODE,ICONV)
      IF (NN.EQ.9) CALL COORD9 (XB,YB,XX,YY,S,T,INODE,ICONV)
C
C     FINAL CHECK ON LOCATION OF POINT
C
      IF (ICONV.EQ.1) RETURN
      A=S
      B=T
      IF (NN.GE.8) A=ABS(S)
      IF (NN.GE.8) B=ABS(T)
      IF (A.LE.STEST.AND.B.LE.STEST) GO TO 60
      RETURN
C
   60 CONTINUE
      IFIND=1
      RETURN
      END
