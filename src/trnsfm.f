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
      SUBROUTINE TRNSFM (X,Y,XO,YO,NC)
C
C     ******************************************************************
C
C     SUBROUTINE TO TRANSFORM POLAR CCORDINATES TO CARTESIAN COORDINATES
C
C     ******************************************************************
C
      DIMENSION X(*), Y(*)
C
C     ******************************************************************
C
      PI=ATAN(1.0)/45.
C
      DO 10 I=1,NC
      R=X(I)
      ANG=Y(I)*PI
      IF (R.EQ.0.) GO TO 10
      X(I)=XO+R*COS(ANG)
      Y(I)=YO+R*SIN(ANG)
   10 CONTINUE
C
      RETURN
      END
