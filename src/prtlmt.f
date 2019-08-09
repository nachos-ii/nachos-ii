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
      SUBROUTINE PRTLMT (LISTEL)
C
C     ******************************************************************
C
C     SUBROUTINE TO SET LIMITS FOR SELECTIVE PRINTING OF OUTPUT
C
C     ******************************************************************
C
      CHARACTER*10 CDATA
C
      COMMON /INDATR/ RDATA(150)
      COMMON /INDATI/ IDATA(150)
      COMMON /INDATC/ CDATA(150)
      COMMON /TAPES/  NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7, 
     1                NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
      COMMON /CONTRL/ IEDIT,IPRINT
C
      DIMENSION LISTEL(*)
C
C     ******************************************************************
C
C     SET OUTPUT LIMITS
C
   10 CONTINUE
      CALL RDFFLD (NIN,RDATA,IDATA,CDATA)
      IF (CDATA(1).EQ.'END       ') GO TO 70
      IF (CDATA(1).EQ.'SINGLE    ') GO TO 20
      IF (CDATA(1).EQ.'STRING    ') GO TO 40
      CALL ERROR ('PRTLMT','UNRECOGNIZED OUTPUT COMMAND',' ',0,' ',0,   
     1'WORD',CDATA(3),1)
C
C     SINGLE ELEMENT SPECIFICATION
C
   20 CONTINUE
      DO 30 I=1,50
      II=IDATA(I+1)
      IF (II.EQ.0) GO TO 10
      ITEMP=ABS(LISTEL(II))
      LISTEL(II)=-ITEMP
   30 CONTINUE
      GO TO 10
C
C     ELEMENT STRING SPECIFICATION
C
   40 CONTINUE
      DO 60 I=1,25
      II=IDATA(2*I)
      JJ=IDATA(2*I+1)
      IF (II.EQ.0) GO TO 10
      DO 50 J=II,JJ
      ITEMP=ABS(LISTEL(J))
      LISTEL(J)=-ITEMP
   50 CONTINUE
   60 CONTINUE
      GO TO 10
C
   70 CONTINUE
      IEDIT=1
C
      RETURN
      END
