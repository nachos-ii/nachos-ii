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
      SUBROUTINE USRBC4 (VALUE, UBD, VBD, PBD, TBD, V1BD, V2BD, 
     *                   XBD, YBD, NNODES, NELEM, TIME, KCYCLE)
C
C     ******************************************************************
C
C     USER SUBROUTINE TO EVALUATE A TIME-DEPENDENT BOUNDARY CONDITION
C
C     ******************************************************************
C
      DIMENSION VALUE(*), UBD(*), VBD(*), PBD(*), TBD(*), V1BD(*),
     *          V2BD(*), XBD(*), YBD(*)
C
C     ******************************************************************
C
      RETURN
      END
