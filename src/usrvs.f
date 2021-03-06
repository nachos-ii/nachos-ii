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
      SUBROUTINE USRVS ( QVALUE, T, VAR1, VAR2, X, Y, NNODES, MAT,  
     *                   NELEM, TIME, ICYCLE, IEQN)
C
C     ******************************************************************
C
C     USER SUBROUTINE TO EVALUATE VOLUMETRIC SOURCE (AUXILIARY EQN)
C
C     ******************************************************************
C
      DIMENSION QVALUE(*), T(*), X(*), Y(*), VAR1(*), VAR2(*)
C
C     ******************************************************************
C
      RETURN
      END
