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
      SUBROUTINE FLUX6 (NSTRST,ST,TT)
C
C     ******************************************************************
C
C     SUBROUTINE TO CALCULATE THE VALUES OF THE SIX NODE SHAPE FUNCTIONS
C     AND THEIR DERIVATIVES AT THE FLUX EVALUATION POINTS.
C     THE POINTS MAY BE ON THE ELEMENT BOUNDARY OR AT INTERIOR POINTS
C
C     ******************************************************************
C
      COMMON /TRI6/ F6Q(6,7),F6L(3,7),DF6QDS(6,7),DF6QDT(6,7)
C
      DIMENSION ST(*), TT(*)
C
C     ******************************************************************
C
C     EVALUATE  FUNCTIONS
C
      DO 10 I=1,NSTRST
      A=ST(I)
      B=TT(I)
      F6Q(1,I)=2.*A*A-A
      F6Q(2,I)=2.*B*B-B
      F6Q(3,I)=2.*(A+B)**2-3.*(A+B)+1.
      F6Q(4,I)=4.*A*B
      F6Q(5,I)=-4.*(B**2+A*B-B)
      F6Q(6,I)=-4.*(A**2+A*B-A)
C
      F6L(1,I)=A
      F6L(2,I)=B
      F6L(3,I)=1.-A-B
C
      DF6QDS(1,I)=4.*A-1.
      DF6QDS(2,I)=0.
      DF6QDS(3,I)=4.*(A+B)-3.
      DF6QDS(4,I)=4.*B
      DF6QDS(5,I)=-4.*B
      DF6QDS(6,I)=4.*(-2.*A-B+1.)
C
      DF6QDT(1,I)=0.
      DF6QDT(2,I)=4.*B-1.
      DF6QDT(3,I)=4.*(A+B)-3.
      DF6QDT(4,I)=4.*A
      DF6QDT(5,I)=4.*(-2.*B-A+1.)
      DF6QDT(6,I)=-4.*A
   10 CONTINUE
C
      RETURN
      END
