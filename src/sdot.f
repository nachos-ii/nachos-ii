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
      FUNCTION SDOT (N,VX,INCRX,VY,INCRY)
C     ******************************************************************
C
C     FUNCTION TO COMPUTE THE INNER PRODUCT OF TWO VECTORS
C
C     THIS ROUTINE IS A SIMPLIFIED VERSION OF THE CRAY LIBRARY 
C     ROUTINE OF THE SAME NAME (ASSUMES INCRX=INCRY) (CFT 1.11)
C     
C     ******************************************************************
C
      DIMENSION VX(N), VY(N)
C
      SDOT=0.
      DO 10 I=1,N,INCRX
      SDOT=SDOT+VX(I)*VY(I)
   10 CONTINUE
C
      RETURN
      END
