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
      SUBROUTINE INVRT3 (A,DET)
C
C     ******************************************************************
C
C     SUBROUTINE TO INVERT A 3X3 MATRIX
C     INVERSE IS RETURNED IN THE ORIGINAL MATRIX
C
C     ******************************************************************
C
      DIMENSION A(3,3)
C
C     ******************************************************************
C
      F11=A(2,2)*A(3,3)-A(2,3)*A(3,2)
      F21=A(1,2)*A(3,3)-A(3,2)*A(1,3)
      F31=A(2,3)*A(1,2)-A(1,3)*A(2,2)
C
      F12=A(2,1)*A(3,3)-A(3,1)*A(2,3)
      F22=A(1,1)*A(3,3)-A(1,3)*A(3,1)
      F32=A(1,1)*A(2,3)-A(2,1)*A(1,3)
C
      F13=A(2,1)*A(3,2)-A(3,1)*A(2,2)
      F23=A(1,1)*A(3,2)-A(3,1)*A(1,2)
      F33=A(1,1)*A(2,2)-A(1,2)*A(2,1)
C
      DET=A(1,1)*F11-A(1,2)*F12+A(1,3)*F13
      IF (DET.EQ.0.0) RETURN
C
      A(1,1)=F11/DET
      A(2,1)=-F12/DET
      A(3,1)=F13/DET
      A(1,2)=-F21/DET
      A(2,2)=F22/DET
      A(3,2)=-F23/DET
      A(1,3)=F31/DET
      A(2,3)=-F32/DET
      A(3,3)=F33/DET
C
      RETURN
      END
