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
      SUBROUTINE GATHER (N,A,B,INDEX)
C
C     ******************************************************************
C
C     SUBROUTINE TO GATHER A VECTOR FROM A SOURCE VECTOR
C
C     THIS ROUTINE IS EQUIVALENT TO THE CRAY LIBRARY ROUTINE
C     OF THE SAME NAME (CFT 1.11)
C     
C     ******************************************************************
C
      DIMENSION A(N), B(N), INDEX(N)
C
      DO 10 I=1,N
      A(I)=B(INDEX(I))
   10 CONTINUE
C
      RETURN 
      END
