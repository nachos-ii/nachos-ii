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
      SUBROUTINE TSTEP (ITYPE,KSTEP,IALGOR,TOLT,DELTN,DELTNM,DELTMN,    
     1DELTMX)
C
C     ******************************************************************
C
C     SUBROUTINE TO COMPUTE A NEW TIME STEP
C
C     ******************************************************************
C
      CHARACTER*15 ETYPE
C
      COMMON /TAPES/  NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7, 
     1                NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
      COMMON /PRBDAT/ IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,    
     1                IPNLTY,PNLTY
      COMMON /BNORM/  BNRMU,BNRMV,BNRMUV,BNRMT,BNRMX1,BNRMX2
C
      DIMENSION ETYPE(4)
C
      DATA (ETYPE(I),I=1,4)/'MOMENTUM','ENERGY','1ST AUXILIARY',        
     1     '2ND AUXILIARY'/
C
C     ******************************************************************
C
      IF (KSTEP.LE.2) RETURN
      IF (IALGOR.EQ.1) THEN
      POWER=1./2.
      FACTOR=2.0
       ELSE
      POWER=1./3.
      FACTOR=3.0*(1.0+DELTNM/DELTN)
      END IF
      IF (ITYPE.GE.2) GO TO 10
C
C     ISOTHERMAL CASE
C
      DELT=DELTN*(TOLT*FACTOR/BNRMUV)**POWER
      IF (DELT.LT.DELTMN) THEN
      DELT=DELTMN
      WRITE (NOUT, 60)
      END IF
      IF (DELT.GT.DELTMX) THEN
      DELT=DELTMX
      WRITE (NOUT, 70)
      END IF
      IF (DELT.LE.DELTN/2.) WRITE (NOUT, 40)
      DELTNM=DELTN
      DELTN=DELT
      RETURN
C
C     NON-ISOTHERMAL, MULTIPLE EQUATION CASE
C
   10 CONTINUE
      IF (BNRMUV.LT.1.0E-10) THEN
      UDELT=0.0
       ELSE
      UDELT=DELTN*(TOLT*FACTOR/BNRMUV)**POWER
      END IF
      IF (BNRMT.LT.1.0E-10) THEN
      TDELT=0.0
       ELSE
      TDELT=DELTN*(TOLT*FACTOR/BNRMT)**POWER
      END IF
      IF (UDELT.LE.TDELT) THEN
      IEQN=1
      DELT=UDELT
       ELSE
      IEQN=2
      DELT=TDELT
      END IF
C
      IF (IVAR1.EQ.0) GO TO 30
      IF (BNRMX1.LT.1.0E-10) THEN
      GO TO 20
       ELSE
      X1DELT=DELTN*(TOLT*FACTOR/BNRMX1)**POWER
      END IF
      IF (X1DELT.LT.DELT) THEN
      IEQN=3
      DELT=X1DELT
      END IF
C
   20 CONTINUE
      IF (IVAR2.EQ.0) GO TO 30
      IF (BNRMX2.LT.1.0E-10) THEN
      GO TO 30
       ELSE
      X2DELT=DELTN*(TOLT*FACTOR/BNRMX2)**POWER
      END IF
      IF (X2DELT.LT.DELT) THEN
      IEQN=4
      DELT=X2DELT
      END IF
C
   30 CONTINUE
      WRITE (NOUT, 50) ETYPE(IEQN)
      IF (DELT.LT.DELTMN) THEN
      DELT=DELTMN
      WRITE (NOUT, 60)
      END IF
      IF (DELT.GT.DELTMX) THEN
      DELT=DELTMX
      WRITE (NOUT, 70)
      END IF
      IF (DELT.LE.DELTN/2.) WRITE (NOUT, 40)
      DELTNM=DELTN
      DELTN=DELT
      RETURN
C
   40 FORMAT (/,10X,'*** WARNING - TIME STEP REDUCED BY MORE THAN 1/2 **
     1*')
   50 FORMAT (/,10X,A,' EQUATION CONTROLS THE TIME STEP')
   60 FORMAT (/,10X,'*** WARNING - TIME STEP RESET TO DELT MIN ***')
   70 FORMAT (/,10X,'*** WARNING - TIME STEP RESET TO DELT MAX ***')
      END
