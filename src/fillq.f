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
      SUBROUTINE FILLQ (X,Y,IFLAG)
C
C     ******************************************************************
C
C     SUBROUTINE TO FILL IN ANY OMITTED BOUNDARY POINTS BY AVERAGING
C     FOR A QUADRILATERAL REGION
C
C     ******************************************************************
C
      DIMENSION X(*), Y(*)
C
C     ******************************************************************
C
      IF (X(9).NE.0.0) GO TO 20
      IF (X(10).NE.0.0) GO TO 20
      IF (X(11).NE.0.0) GO TO 20
      IF (X(12).NE.0.0) GO TO 20
      IF (X(5).NE.0.0) GO TO 10
      IF (X(6).NE.0.0) GO TO 10
      IF (X(7).NE.0.0) GO TO 10
      IF (X(8).NE.0.0) GO TO 10
      RETURN
C
   10 CONTINUE
      IF (X(8).EQ.0.0) X(8)=(X(1)+X(4))/2.
      IF (Y(8).EQ.0.0) Y(8)=(Y(1)+Y(4))/2.
      IF (X(7).EQ.0.0) X(7)=(X(4)+X(3))/2.
      IF (Y(7).EQ.0.0) Y(7)=(Y(4)+Y(3))/2.
      IF (X(6).EQ.0.0) X(6)=(X(3)+X(2))/2.
      IF (Y(6).EQ.0.0) Y(6)=(Y(3)+Y(2))/2.
      IF (X(5).EQ.0.0) X(5)=(X(2)+X(1))/2.
      IF (Y(5).EQ.0.0) Y(5)=(Y(2)+Y(1))/2.
      IFLAG=2
      RETURN
C
   20 CONTINUE
      IF (X(5).EQ.0.0) X(5)=X(1)+(X(2)-X(1))/3.
      IF (Y(5).EQ.0.0) Y(5)=Y(1)+(Y(2)-Y(1))/3.
      IF (X(6).EQ.0.0) X(6)=X(2)+(X(3)-X(2))/3.
      IF (Y(6).EQ.0.0) Y(6)=Y(2)+(Y(3)-Y(2))/3.
      IF (X(7).EQ.0.0) X(7)=X(3)+(X(4)-X(3))/3.
      IF (Y(7).EQ.0.0) Y(7)=Y(3)+(Y(4)-Y(3))/3.
      IF (X(8).EQ.0.0) X(8)=X(4)+(X(1)-X(4))/3.
      IF (Y(8).EQ.0.0) Y(8)=Y(4)+(Y(1)-Y(4))/3.
      IF (X(9).EQ.0.0) X(9)=X(1)+2.*(X(2)-X(1))/3.
      IF (Y(9).EQ.0.0) Y(9)=Y(1)+2.*(Y(2)-Y(1))/3.
      IF (X(10).EQ.0.0) X(10)=X(2)+2.*(X(3)-X(2))/3.
      IF (Y(10).EQ.0.0) Y(10)=Y(2)+2.*(Y(3)-Y(2))/3.
      IF (X(11).EQ.0.0) X(11)=X(3)+2.*(X(4)-X(3))/3.
      IF (Y(11).EQ.0.0) Y(11)=Y(3)+2.*(Y(4)-Y(3))/3.
      IF (X(12).EQ.0.0) X(12)=X(4)+2.*(X(1)-X(4))/3.
      IF (Y(12).EQ.0.0) Y(12)=Y(4)+2.*(Y(1)-Y(4))/3.
      IFLAG=3
C
      RETURN
      END
