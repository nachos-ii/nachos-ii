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
      SUBROUTINE PREVUE
C
C     ******************************************************************
C
C     SUBROUTINE TO READ INPUT DATA FILE, VERIFY COMMAND SEQUENCE AND
C     PRESET CRITICAL PARAMETERS FOR BLANK COMMON ALLOCATION
C
C     ******************************************************************
C
      CHARACTER*10 CWORD,CDATA
C
      COMMON /INDATR/ RDATA(150)
      COMMON /INDATI/ IDATA(150)
      COMMON /INDATC/ CDATA(150)
      COMMON /TAPES/  NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7, 
     1                NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
      COMMON /LMTDAT/ MAXELM,MAXNOD,MAXBC,MAXEDG,MAXMAT,MAXNPT,MAXPLT,  
     1                MAXBRY,MAXBLK
      COMMON /PRBDAT/ IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,    
     1                IPNLTY,PNLTY
      COMMON /CMMNDS/ CWORD(10)
C
      DIMENSION ICMMND(7)
C
      DATA (ICMMND(I),I=1,7)/7*0/
      DATA ISTOP/0/
C
C     ******************************************************************
C
C     READ INPUT FILE
C
   10 CONTINUE
      CALL RDFFLD (NIN,RDATA,IDATA,CDATA)
C
C     BRANCH TO SPECIFIED COMMAND
C
      DO 20 I=1,10
      IF (CDATA(1).EQ.CWORD(I)) GO TO (30, 40, 50, 60, 70, 80, 90,      
     1100, 180, 190), I
   20 CONTINUE
      GO TO 10
C
C     MATERIALS COMMAND
C
   30 CONTINUE
      ICMMND(1)=1
      GO TO 10
C
C     MESH COMMAND
C
   40 CONTINUE
      ICMMND(2)=1
      GO TO 10
C
C     ELEMENTS COMMAND
C
   50 CONTINUE
      IF (ICMMND(2).EQ.0) THEN
      CALL ERROR ('PREVUE','MESH COMMAND MUST PRECEDE ELEMENTS COMMAND',
     1' ',0,' ',0,' ',' ',0)
      ISTOP=1
      END IF
      ICMMND(3)=1
      GO TO 10
C
C     FORMKF COMMAND
C
   60 CONTINUE
      IF (ICMMND(1).EQ.0) THEN
      CALL ERROR ('PREVUE','MATERIALS COMMAND MUST PRECEDE FORMKF COMMAN
     1D',' ',0,' ',0,' ',' ',0)
      ISTOP=1
      END IF
      IF (ICMMND(2).EQ.0) THEN
      CALL ERROR ('PREVUE','MESH COMMAND MUST PRECEDE FORMKF COMMAND',  
     1' ',0,' ',0,' ',' ',0)
      ISTOP=1
      END IF
      IF (ICMMND(3).EQ.0) THEN
      CALL ERROR ('PREVUE','ELEMENTS COMMAND MUST PRECEDE FORMKF COMMAND
     1',' ',0,' ',0,' ',' ',0)
      ISTOP=1
      END IF
      ICMMND(4)=1
      IF (CDATA(6).EQ.'LINEAR    ') IPFUNC=1
      IF (CDATA(4).EQ.'PENALTY   ') IPFUNC=0
      GO TO 10
C
C     OUTPUT COMMAND
C
   70 CONTINUE
      IF (ICMMND(3).EQ.0) THEN
      CALL ERROR ('PREVUE','MESH AND ELEMENTS COMMANDS MUST PRECEDE OUTP
     1UT COMMAND',' ',0,' ',0,' ',' ',0)
      ISTOP=1
      END IF
      GO TO 10
C
C     SOLVE COMMAND
C
   80 CONTINUE
      IF (ICMMND(1).EQ.0) THEN
      CALL ERROR ('PREVUE','MATERIALS COMMAND MUST PRECEDE SOLVE COMMAND
     1',' ',0,' ',0,' ',' ',0)
      ISTOP=1
      END IF
      IF (ICMMND(2).EQ.0) THEN
      CALL ERROR ('PREVUE','MESH COMMAND MUST PRECEDE SOLVE COMMAND',   
     1' ',0,' ',0,' ',' ',0)
      ISTOP=1
      END IF
      IF (ICMMND(3).EQ.0) THEN
      CALL ERROR ('PREVUE','ELEMENTS COMMAND MUST PRECEDE SOLVE COMMAND'
     1,' ',0,' ',0,' ',' ',0)
      ISTOP=1
      END IF
      IF (ICMMND(4).EQ.0) THEN
      CALL ERROR ('PREVUE','FORMKF COMMAND MUST PRECEDE SOLVE COMMAND', 
     1' ',0,' ',0,' ',' ',0)
      ISTOP=1
      END IF
      ICMMND(5)=1
      GO TO 10
C
C     QUASI-NEWTON CHECK (NOT OPERATIONAL IN THIS VERSION)
C
C  90 CONTINUE
C
C     STREAM COMMAND
C
   90 CONTINUE
      IF (ICMMND(5).EQ.0) THEN
      CALL ERROR ('PREVUE','SOLVE COMMAND MUST PRECEDE STREAM COMMAND', 
     1' ',0,' ',0,' ',' ',0)
      ISTOP=1
      END IF
      ICMMND(6)=1
      GO TO 10
C
C     FLUX COMMAND
C
  100 CONTINUE
      IF (ICMMND(5).EQ.0) THEN
      CALL ERROR ('PREVUE','SOLVE COMMAND MUST PRECEDE FLUX COMMAND',   
     1' ',0,' ',0,' ',' ',0)
      ISTOP=1
      END IF
      IF (CDATA(4).EQ.'SAVE      ') ICMMND(7)=1
      GO TO 10
C
C     POST COMMAND
C
  180 CONTINUE
      IF (ICMMND(5).EQ.0) THEN
      CALL ERROR ('PREVUE','SOLVE COMMAND MUST PRECEDE POST COMMAND',   
     1' ',0,' ',0,' ',' ',0)
      ISTOP=1
      END IF
      GO TO 10
C
C     STOP COMMAND
C
  190 CONTINUE
      IF (ISTOP.EQ.1) CALL ERROR ('PREVUE','ERROR IN SEQUENCE OF INPUT C
     1OMMANDS',' ',0,' ',0,' ',' ',1)
C
      RETURN
      END
