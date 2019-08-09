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
      SUBROUTINE BCKSUB (IRESOL,NBLKS,IPTECV,KROWV,ILHED,IDISK,QQ,TEMP, 
     1SK,R1,ECV)
C
C     ******************************************************************
C
C     SUBROUTINE TO PERFORM BACKSUBTITUTION FOR THE EQUATIONS REDUCED
C     BY SUBROUTINE FRONT
C
C     ******************************************************************
C
      COMMON /TAPES/ NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7,  
     1               NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
      COMMON /SZDAT/ NUMEL,NUMNOD,NUMVAR,NUMDOF,NODSOL
C
      DIMENSION ILHED(*), IDISK(*), KROWV(*)
      DIMENSION QQ(*), ECV(*), SK(*), TEMP(*), R1(*)
C
C     ******************************************************************
C
C     INITIALIZE SOLUTION VECTOR
C
      DO 10 I=1,NUMDOF
      SK(I)=0.
   10 CONTINUE
C
C     BACKSUBSTITUTION
C
      DO 30 I=1,NUMDOF
C
C     RETRIEVE A PORTION OF "ECV" FROM DISK IF NECESSARY
C
      IF (IPTECV.EQ.1) THEN
      IPTECV=IDISK(NBLKS)+1
C
C     THE FOLLOWING BUFFERED INPUT STATEMENT MAY BE REPLACED WITH A
C     READ STATEMENT OF THE FORM
C     READ (NTP10) (ECV(II),II=1,IPTECV-1)
C     ALSO, THE FOLLOWING UNIT STATEMENT SHOULD BE REMOVED IF BUFFERED
C     INPUT IS NOT USED
C
      BACKSPACE (NTP10)
c      BUFFER IN (NTP10,0) (ECV(1),ECV(IPTECV-1))
      READ (NTP10) (ECV(II),II=1,IPTECV-1)
      BACKSPACE (NTP10)
      NBLKS=NBLKS-1
c      IF (UNIT(NTP10).GT.0.) CALL ERROR ('BCKSUB','ERROR IN BUFFER IN OP
c     1ERATION',' ',0,' ',0,' ',' ',1)
      END IF
C
C     RETRIEVE ELIMINATED EQUATIONS FROM "ECV" ONE-BY-ONE
C
      IF (IRESOL.EQ.0) THEN
      IPTECV=IPTECV-2
       ELSE
      IPTECV=IPTECV-4
      END IF
      LCOL=KROWV(NUMDOF+1-I)
      KRO1=ECV(IPTECV)
      LPVCO1=ECV(IPTECV+1)
      IF (IRESOL.EQ.0) THEN
      IPTECV=IPTECV-2*LCOL
       ELSE
      IPTECV=IPTECV-4*LCOL
      END IF
      LCO=ECV(IPTECV+LCOL-1+LPVCO1)
C
C     SOLVE CURRENT EQUATION FOR THE "LCO"-TH D.O.F.
C
      DO 20 J=1,LCOL
      QQ(J)=ECV(IPTECV-1+J)
   20 ILHED(J)=ECV(IPTECV+LCOL-1+J)
      QQ(LPVCO1)=0.0
      CALL GATHER (LCOL,TEMP,SK,ILHED)
      SK(LCO)=R1(KRO1)-SDOT(LCOL,QQ,1,TEMP,1)
   30 CONTINUE
      RETURN
C
      END
