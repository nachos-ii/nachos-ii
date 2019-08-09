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
      SUBROUTINE RDINPT
C
C     ******************************************************************
C
C     SUBROUTINE TO READ,CHECK AND PRINT INPUT DATA
C     INPUT IS READ UNDER A FREE FIELD FORMAT IN SUBROUTINES FREFLD AND
C     GETINP, DECODED AND PLACED ON AN INTERNAL WORKING FILE
C     SUBROUTINES FREFLD AND GETINP ARE PART OF THE EXTERNAL "SUPES"
C     LIBRARY (SAND86-0911)
C
C     ******************************************************************
C
      CHARACTER*8  CODNAM,VERSN,RDATE,RTIME,HRDWRE,SFTWRE
      CHARACTER*10 CVALUE,NAME
      CHARACTER*80 HED,CMMNT,LABEL
C
      COMMON /HEADER/ HED,CMMNT(10)
      COMMON /RUNDAT/ CODNAM,VERSN,RDATE,RTIME,HRDWRE,SFTWRE
      COMMON /INDATR/ RVALUE(150)
      COMMON /INDATI/ IVALUE(150)
      COMMON /INDATC/ CVALUE(150)
      COMMON /TAPES/  NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7, 
     1                NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
C
      DIMENSION KVALUE(150)
C
C     ******************************************************************
C
C     PRINT RUN-TIME DATA
C
      WRITE (NOUT, 50)
      CALL BANNER (84,'NACHOS II',NOUT)
      WRITE (NOUT, 60)
      WRITE (NOUT, 70) VERSN
      WRITE (NOUT, 80) RDATE
      WRITE (NOUT, 90) RTIME
      WRITE (NOUT, 100) HRDWRE,SFTWRE
C
C     PUT LARGE FILES ON SOLID STATE DISK, IF POSSIBLE
C     NOT USED ON CRAY UNICOS
C
C      CALL OPNSSD
C
      WRITE (NOUT, 60)
      WRITE (NOUT, 110)
C
C     EXTRACT TITLE AND COMMENTS FROM INPUT FILE
C
      NLINE=0
      NCMMNT=0
   10 CONTINUE
      CALL GETINP (NIN,0,' ',LABEL,IOSTAT)
      IF (IOSTAT.LT.0) CALL ERROR ('RDINPT','EOF ENCOUNTERED','LINE NUMB
     1ER',NLINE,' ',0,'AT SUBROUTINE','GETINP',1)
      IF (IOSTAT.GT.0) CALL ERROR ('RDINPT','READ ERROR','LINE NUMBER', 
     1NLINE,' ',0,'AT SUBROUTINE','GETINP',1)
      IF (LABEL(1:1).NE.'$') GO TO 20
      NLINE=NLINE+1
      IF (NLINE.EQ.1) THEN
      HED=LABEL(2:80)
       ELSE
      NCMMNT=NCMMNT+1
      CMMNT(NCMMNT)=LABEL(2:80)
      END IF
      GO TO 10
C
C     LOOP OVER EACH LINE OF INPUT
C
   20 CONTINUE
      BACKSPACE (NIN)
      NLINE=0
   30 CONTINUE
      NLINE=NLINE+1
      CALL FREFLD (NIN,NOUT,'AUTO',150,IOSTAT,NFIELD,KVALUE,CVALUE,     
     1IVALUE,RVALUE)
      IF (IOSTAT.GT.0) CALL ERROR ('RDINPT','READ ERROR','LINE NUMBER', 
     1NLINE,' ',0,'AT SUBROUTINE','FREFLD',1)
      IF (IOSTAT.LT.0) GO TO 40
C
C     WRITE DECODED DATA TO INTERNAL FILE
C
      IF (KVALUE(1).EQ.-1) GO TO 30
      WRITE (NTP1) RVALUE, IVALUE,CVALUE
      GO TO 30
C
   40 CONTINUE
      RETURN
C
   50 FORMAT ('1',/////)
   60 FORMAT (//,5X,'***************************************************
     1********',//)
   70 FORMAT (5X,'VERSION -- ',A8,//)
   80 FORMAT (5X,'DATE OF EXECUTION -- ',A8,//)
   90 FORMAT (5X,'TIME OF EXECUTION -- ',A8,//)
  100 FORMAT (5X,'SYSTEM HARDWARE -- ',A8,//,5X,'SYSTEM SOFTWARE -- ',  
     1A8)
  110 FORMAT (2X,'LINE:  DIRECT LIST OF INPUT DATA' ,/)
      END
