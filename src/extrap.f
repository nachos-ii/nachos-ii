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
      SUBROUTINE EXTRAP (KIND,GPFLUX)
C
C     ******************************************************************
C
C     SUBROUTINE TO EXTRAPOLATE INTEGRATION POINT VALUES TO THE NODES
C
C     ******************************************************************
C
      DIMENSION GPFLUX(*), DUM(4)
      DIMENSION TRANSQ(4,4), TRANST(3,3)
C
      DATA ((TRANSQ(I,J),I=1,4),J=1,4)/1.866025404,-.5,.133974596,-.5,  
     1     -.5,1.866025404,-.5,.133974596,.133974596,-.5,1.866025404,   
     2     -.5,-.5,.133974596,-.5,1.866025404/
      DATA ((TRANST(I,J),I=1,3),J=1,3)/1.444444444,-.222222222,         
     1     -.222222222,-.222222222,1.444444444,-.222222222,-.222222222, 
     2     -.222222222,1.444444444/
C
C     ******************************************************************
C
      GO TO (10, 10, 50, 50, 50, 50), KIND
C
C     TRIANGLE TRANSFORMATION
C
   10 CONTINUE
      DO 20 I=1,3
      DUM(I)=GPFLUX(I)
   20 CONTINUE
      DO 40 I=1,3
      SUM=0.
      DO 30 J=1,3
      SUM=SUM+TRANST(I,J)*DUM(J)
   30 CONTINUE
      GPFLUX(I)=SUM
   40 CONTINUE
      RETURN
C
C     QUADRILATERAL TRANSFORMATION
C
   50 CONTINUE
      DO 60 I=1,4
      DUM(I)=GPFLUX(I)
   60 CONTINUE
      DO 80 I=1,4
      SUM=0.
      DO 70 J=1,4
      SUM=SUM+TRANSQ(I,J)*DUM(J)
   70 CONTINUE
      GPFLUX(I)=SUM
   80 CONTINUE
C
      RETURN
      END
