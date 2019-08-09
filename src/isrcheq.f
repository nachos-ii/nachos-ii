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
      FUNCTION ISRCHEQ (N,IARRAY,INCR,ITARGET)
C
C     ******************************************************************
C
C     FUNCTION TO RETURN THE FIRST LOCATION IN AN INTEGER ARRAY 
C     THAT IS EQUAL TO THE INTEGER TARGET VALUE     
C
C     THIS ROUTINE IS EQUIVALENT TO THE CRAY LIBRARY ROUTINE
C     OF THE SAME NAME (CFT 1.14)
C
C     ******************************************************************
C
      DIMENSION IARRAY (N)
C
      DO 10 I=1,N,INCR
      IF (IARRAY(I).EQ.ITARGET) THEN
      ISRCHEQ=I
      RETURN
      END IF
   10 CONTINUE
      ISRCHEQ=N+1
C
      RETURN
      END
