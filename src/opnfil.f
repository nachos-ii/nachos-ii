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
      SUBROUTINE OPNFIL (IOPEN)
C
C     ******************************************************************
C
C     SUBROUTINE TO OPEN REQUIRED FILES
C
C     ******************************************************************
C
      CHARACTER*132 FILNAM
C
      COMMON /TAPES/  NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7, 
     1                NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
      COMMON /EXOID/  NTPX0,NTPX12
      COMMON /NTPDAT/ IFILES(16)
C
C     ******************************************************************
C
C     BRANCH TO NEEDED FILES
C
      GO TO (10, 20), IOPEN
C
C     OPEN FILES REQUIRED FOR STANDARD EXECUTION
C
   10 CONTINUE
C
C     OUTPUT FILE
C
      IFILES(2)=1
      IUNIT=NOUT
      CALL EXNAME (NOUT,FILNAM,LN)
      OPEN (UNIT=NOUT, FILE=FILNAM(1:LN), STATUS='unknown',FORM='formatt
     1ed', ERR=30)
C
C     IPUT FILE
C
      IFILES(1)=1
      IUNIT=NIN
      CALL EXNAME (NIN,FILNAM,LN)
      OPEN (UNIT=NIN, FILE=FILNAM(1:LN), STATUS='old', FORM='formatted',
     1 ERR=30)
C
C     EXTERNAL MESH DATA FILE ('GENESIS' FORMAT)
C
      IFILES(3)=1
      IUNIT=NTP0
      CALL EXNAME (NTP0,FILNAM,LN)
      OPEN (UNIT=NTP0, FILE=FILNAM(1:LN), STATUS='unknown', FORM='unform
     1atted', ERR=30)
C
C     FREE FIELD INPUT FILE
C
      IFILES(4)=1
      IUNIT=NTP1
      OPEN (UNIT=NTP1, STATUS='unknown', FORM='unformatted', ERR=30)
C
C     ELEMENT DATA FILE
C
      IFILES(5)=1
      IUNIT=NTP2
      OPEN (UNIT=NTP2, STATUS='unknown', FORM='unformatted', ERR=30)
C
C     ELEMENT MATRIX FILE (PENALTY FORMULATION)
C
      IFILES(6)=1
      IUNIT=NTP3
      OPEN (UNIT=NTP3, STATUS='unknown', FORM='unformatted', ERR=30)
C
C     ELEMENT MATRIX FILE
C
      IFILES(7)=1
      IUNIT=NTP4
      OPEN (UNIT=NTP4, STATUS='unknown', FORM='unformatted', ERR=30)
C
C     STREAM FUNCTION DATA FILE
C
      IFILES(8)=1
      IUNIT=NTP5
      OPEN (UNIT=NTP5, STATUS='unknown', FORM='unformatted', ERR=30)
C
C     SPECIAL POINTS DATA FILE
C
      IFILES(9)=1
      IUNIT=NTP6
      OPEN (UNIT=NTP6, STATUS='unknown', FORM='unformatted', ERR=30)
C
C     FLUX DATA FILE
C
      IFILES(10)=1
      IUNIT=NTP7
      OPEN (UNIT=NTP7, STATUS='unknown', FORM='unformatted', ERR=30)
C
C     STRESS DATA FILE
C
      IFILES(11)=1
      IUNIT=NTP8
      OPEN (UNIT=NTP8, STATUS='unknown', FORM='unformatted', ERR=30)
C
C     SOLUTION DATA FILE
C
      IFILES(12)=1
      IUNIT=NTP9
      CALL EXNAME (NTP9,FILNAM,LN)
      OPEN (UNIT=NTP9, FILE=FILNAM(1:LN), STATUS='new', FORM='unformatte
     1d', ERR=30)
C
C     MATRIX SOLVER SCRATCH FILE
C
      IFILES(13)=1
      IUNIT=NTP10
      OPEN (UNIT=NTP10, STATUS='unknown', FORM='unformatted', ERR=30)
C
C     SCRATCH FILE
C
      IFILES(14)=1
      IUNIT=NTP11
      OPEN (UNIT=NTP11, STATUS='unknown', FORM='unformatted', ERR=30)
C
C     POST-PROCESSING DATA FILE ('EXODUS' FORMAT)
C
      IFILES(15)=1
      IUNIT=NTP12
      CALL EXNAME (NTP12,FILNAM,LN)
      OPEN (UNIT=NTP12, FILE=FILNAM(1:LN), STATUS='new', FORM='unformatt
     1ed', ERR=30)
      RETURN
C
C     OPEN FILE WITH PREVIOUS SOLUTION DATA (RESTART)
C
   20 CONTINUE
      IFILES(16)=1
      IUNIT=NTP13
      CALL EXNAME (NTP13,FILNAM,LN)
      OPEN (UNIT=NTP13, FILE=FILNAM(1:LN), STATUS='old', FORM='unformatt
     1ed', ERR=30)
      RETURN
C
   30 CONTINUE
      CALL ERROR ('OPNFIL','ERROR OPENING FILE','UNIT NUMBER',IUNIT,    
     1' ',0,' ',' ',1)
      END
