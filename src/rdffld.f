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
      SUBROUTINE RDFFLD (NTAPE,X,IY,Z)
C
C     ******************************************************************
C
C     SUBROUTINE TO READ A LINE OF INPUT FROM THE NACHOS INPUT FILE
C
C     ******************************************************************
C
      CHARACTER*10 Z
C
      DIMENSION X(150), IY(150), Z(150)
C
C     ******************************************************************
C
C     READ DATA
C
      READ (NTAPE, END=10) X,IY,Z
      RETURN
   10 CONTINUE
      CALL ERROR ('RDFFLD','END OF DATA ENCOUNTERED',' ',0,' ',0,' ',   
     1' ',1)
      END
