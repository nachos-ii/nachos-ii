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
      SUBROUTINE USRBC1 ( VALUE, UBD, VBD, PBD, TBD, V1BD, V2BD, 
     *                    XBD, YBD, NNODES, NELEM, TIME, KCYCLE)
C
C     ******************************************************************
C
C     USER SUBROUTINE TO EVALUATE A TIME-DEPENDENT BOUNDARY CONDITION
C
C     ******************************************************************
C
      DIMENSION VALUE(*), UBD(*), VBD(*), PBD(*), TBD(*), V1BD(*), V2BD(
     1*), XBD(*), YBD(*)
C
C     ******************************************************************
C
      RETURN
      END
      SUBROUTINE USRBC2 ( VALUE, UBD, VBD, PBD, TBD, V1BD, V2BD,
     *                    XBD, YBD, NNODES, NELEM, TIME, KCYCLE)
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
      SUBROUTINE USRBC3 ( VALUE, UBD, VBD, PBD, TBD, V1BD, V2BD,
     *                    XBD, YBD, NNODES, NELEM, TIME, KCYCLE)
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
      SUBROUTINE USRBC5 ( VALUE, UBD, VBD, PBD, TBD, V1BD, V2BD, 
     *                    XBD, YBD, NNODES, 1NELEM, TIME, KCYCLE)
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
      SUBROUTINE USRBC6 ( VALUE, UBD, VBD, PBD, TBD, V1BD, V2BD,
     *                    XBD, YBD, NNODES, NELEM, TIME, KCYCLE)
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
      SUBROUTINE USRCON ( AK, T, S, VAR1, VAR2, NNODES, MAT)
C
C     ******************************************************************
C
C     USER SUBROUTINE TO EVALUATE THERMAL CONDUCTIVITY
C
C     ******************************************************************
C
      DIMENSION AK(*), T(*), S(*), VAR1(*), VAR2(*)
C
C     ******************************************************************
C
      RETURN
      END
      SUBROUTINE  USRCP (RHOCP, T, VAR1, VAR2, NNODES, MAT, RHO)
C
C     ******************************************************************
C
C     USER SUBROUTINE TO EVALUATE SPECIFIC HEAT
C
C     ******************************************************************
C
      DIMENSION RHOCP(*), T(*), VAR1(*), VAR2(*)
C
C     ******************************************************************
C
      RETURN
      END
      SUBROUTINE USRDIF ( AD, T, VAR1, VAR2, NNODES, MAT, IEQN)
C
C     ******************************************************************
C
C     USER SUBROUTINE TO EVALUATE DIFFUSION COEFFICIENT (AUXILIARY EQN)
C
C     ******************************************************************
C
      DIMENSION AD(*), T(*), VAR1(*), VAR2(*)
C
C     ******************************************************************
C
      RETURN
      END
      SUBROUTINE USREX ( ZETA, T, VAR1, VAR2, NNODES, MAT, IEQN)
C
C     ******************************************************************
C
C     USER SUBROUTINE TO EVALUATE EXPANSION COEFFICIENT (AUXILIARY EQN)
C
C     ******************************************************************
C
      DIMENSION ZETA(*), T(*), VAR1(*), VAR2(*)
C
C     ******************************************************************
C
      RETURN
      END
      SUBROUTINE USREXT ( BETA, T, VAR1, VAR2, NNODES, MAT)
C
C     ******************************************************************
C
C     USER SUBROUTINE TO EVALUATE THERMAL EXPANSION COEFFICIENT
C
C     ******************************************************************
C
      DIMENSION BETA(*), T(*), VAR1(*), VAR2(*)
C
C     ******************************************************************
C
      RETURN
      END
      SUBROUTINE USRGRV ( GX, GY, X, Y, NNODES, NELEM)
C
C     ******************************************************************
C
C     USER SUBROUTINE TO EVALUATE VARIABLY ORIENTED GRAVITY VECTOR
C
C     ******************************************************************
C
      DIMENSION X(*), Y(*)
C
C     ******************************************************************
C
      RETURN
      END
      SUBROUTINE USRHTC ( HT, TSURF, TREF, XSURF, YSURF, TIME, NELEM,
     *                    IVALUE)
C
C     ******************************************************************
C
C     USER SUBROUTINE TO EVALUATE VARIABLE HEAT TRANSFER COEFFICIENT
C
C     ******************************************************************
C
      RETURN
      END
      SUBROUTINE USRMSH ( X, Y, MAXI, MAXJ)
C
C     ******************************************************************
C
C     USER SUBROUTINE TO GENERATE MESH POINT LOCATIONS
C
C     ******************************************************************
C
      DIMENSION X(MAXI,*), Y(MAXI,*)
C
C     ******************************************************************
C
      RETURN
      END
      SUBROUTINE USRVHS ( QVALUE, T, VAR1, VAR2, X, Y, NNODES, MAT, 
     *                    NELEM, TIME, ICYCLE)
C
C     ******************************************************************
C
C     USER SUBROUTINE TO EVALUATE VOLUMETRIC HEAT SOURCE
C
C     ******************************************************************
C
      DIMENSION QVALUE(*), T(*), X(*), Y(*), VAR1(*), VAR2(*)
C
C     ******************************************************************
C
      RETURN
      END
      SUBROUTINE USRVIS ( AMU, T, S, VAR1, VAR2, NNODES, MAT, ICYCLE)
      DIMENSION AMU(*), T(*), S(*), VAR1(*), VAR2(*)
C
C     ***************************************************************
C
C     USER SUBROUTINE TO EVALUATE DYNAMIC VISCOSITY
C
C     ***************************************************************
C
      RETURN
      END
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
