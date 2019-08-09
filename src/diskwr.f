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
      SUBROUTINE DISKWR (IRESOL,MXBLKS,MXFRNT,                          
     1 MXBUFF,NBLKS,KROW,LCOL,KRO,KPVRO,LPVCO,IPTECV,PIVOT,ILHED,IKHED, 
     2 LHED,KHED,IDISK,QQ,PVKOL,EQ,ECV)
C
C     ******************************************************************
C
C     SUBROUTINE TO BUFFER THE REDUCED EQUATIONS TO AND FROM A DISK
C     FILE. STORAGE OF THE EQUATIONS IS INITIALLY IN VECTOR ECV.
C
C     ******************************************************************
C
      COMMON /TAPES/ NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7,  
     1               NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
C
      DIMENSION ILHED(*), IKHED(*), LHED(*), KHED(*), IDISK(*)
      DIMENSION QQ(*), PVKOL(*), ECV(*)
      DIMENSION EQ(MXFRNT,*)
C
C     ******************************************************************
C
      ITEST=IPTECV+2*KROW+3
      IF (IRESOL.GE.1) ITEST=ITEST+2*KROW+2
      IF (ITEST.GT.MXBUFF) THEN
      NBLKS=NBLKS+1
      IF (NBLKS.GT.MXBLKS) THEN
      WRITE (NOUT, 40) MXBLKS
      CALL ERROR ('DISKWR','SEE PREVIOUS MESSAGE',' ',0,' ',0,' ',' ',1)
      END IF
C
C     THE FOLLOWING BUFFERED OUTPUT STATEMENT MAY BE REPLACED WITH A
C     WRITE STATEMENT OF THE FORM
C     WRITE (NTP10) (ECV(II),II=1,IPTECV-1)
C     ALSO, THE FOLLOWING UNIT STATEMENT SHOULD BE REMOVED IF BUFFERED
C     OUTPUT IS NOT USED
C
c      BUFFER OUT (NTP10,0) (ECV(1),ECV(IPTECV-1))
      WRITE (NTP10) (ECV(II),II=1,IPTECV-1)
      IDISK(NBLKS)=IPTECV-1
      IPTECV=1
c      IF (UNIT(NTP10).GT.0.) CALL ERROR ('DISKWR','ERROR IN BUFFER OUT O
c     1PERATION',' ',0,' ',0,' ',' ',1)
      END IF
C
C     SAVE THE ELIMINATED EQUATION AND HEADING VECTORS IN "ECV"
C
      IF (IRESOL.EQ.0) THEN
CDIR$ IVDEP
      DO 10 I=1,KROW
      ECV(IPTECV-1+I)=QQ(I)
      ECV(IPTECV+KROW-1+I)=ILHED(I)
   10 CONTINUE
      IPTECV=IPTECV+2*KROW+2
      ECV(IPTECV-2)=KRO
      ECV(IPTECV-1)=LPVCO
       ELSE
CDIR$ IVDEP
      DO 20 I=1,KROW
      ECV(IPTECV-1+I)=QQ(I)
      ECV(IPTECV+KROW-1+I)=ILHED(I)
      ECV(IPTECV+2*KROW-1+I)=PVKOL(I)
      ECV(IPTECV+3*KROW-1+I)=KHED(I)
   20 CONTINUE
      IPTECV=IPTECV+4*KROW+4
      ECV(IPTECV-4)=KRO
      ECV(IPTECV-3)=LPVCO
      ECV(IPTECV-2)=PIVOT
      ECV(IPTECV-1)=KPVRO
      END IF
C
C     TEST FOR OUT-OF-BOUNDS SUBSCRIPTS IN "ECV" VECTOR
C
      IF (IPTECV.GT.MXBUFF) THEN
      WRITE (NOUT, 50) MXBUFF,IPTECV
      CALL ERROR ('DISKWR','SEE PREVIOUS MESSAGE',' ',0,' ',0,' ',' ',1)
      END IF
C
C     DELETE THE FIRST PIVOTAL ROW AND COLUMN FROM THE FRONTAL MATRIX
C
      TEMPS=EQ(KROW,KROW)
CDIR$ IVDEP
      DO 30 K=1,KROW
      EQ(K,LPVCO)=EQ(K,KROW)
      EQ(KPVRO,K)=EQ(KROW,K)
      EQ(K,KROW)=0.
      EQ(KROW,K)=0.
   30 CONTINUE
      EQ(KPVRO,LPVCO)=TEMPS
      EQ(KPVRO,KROW)=0.
      EQ(KROW,LPVCO)=0.
C
C     REARRANGE THE HEADING VECTORS
C
      KHED(KPVRO)=KHED(KROW)
      LHED(LPVCO)=LHED(KROW)
      IKHED(KPVRO)=IKHED(KROW)
      ILHED(LPVCO)=ILHED(KROW)
      RETURN
C
   40 FORMAT (/,10X,'MAXIMUM NUMBER OF DISK BLOCK TRANSFERS EXCEEDED',  
     1/,10X,'MXBLKS = ',I5)
   50 FORMAT (/,10X,'MAXIMUM ALLOWED BUFFER SIZE EXCEEDED',/,10X,       
     1'MXBUFF = ',I10,5X,'IPTECV AFTER ASSEMBLY =',I10)
      END
