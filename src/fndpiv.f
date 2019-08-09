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
      SUBROUTINE FNDPIV (IPIVOT,IFULPV,NELL,IPIV,NPIV,LC,KR,MXFRNT,     
     1 NUMDOF,LPIV,KPIV,ILHED,IKHED,NRP,NCP,RPIVOT,SPIVOT,LPVCO,KPVRO,  
     2 KRO,PIVOT,DET,PVKOL1,QQ1,PVKOL2,QQ2,PVKOL3,QQ3,EQ)
C
C     ******************************************************************
C
C     SUBROUTINE TO COMPUTE AND TEST PIVOT VALUES USED IN SUBROUTINE
C     FRONT. PIVOT SELECTION IS MADE BY FULL PIVOTING STRATEGY AMOUNG
C     THE FULLY SUMMED EQUATIONS IN THE FRONTAL MATRIX.
C
C     ******************************************************************
C
      COMMON /TAPES/ NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7,  
     1               NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
C
      DIMENSION LPIV(*), KPIV(*), ILHED(*), IKHED(*), NRP(*), NCP(*)
      DIMENSION PVKOL1(*), QQ1(*), PVKOL2(*), QQ2(*), PVKOL3(*), QQ3(*)
      DIMENSION EQ(MXFRNT,*)
C
C     ******************************************************************
C
C     IF DESIRED, CHECK THE FIRST PIVOT
C
      IF (IPIVOT.EQ.0) THEN
      LPVCO=LPIV(1)
      KPVRO=KPIV(1)
      L0=1
      K0=1
      PIVOT=EQ(KPVRO,LPVCO)
      IF (IPIV.GE.1) THEN
      PIVOT=PIVOT-PVKOL1(KPVRO)*QQ1(LPVCO)
      IF (IPIV.GE.2) THEN
      PIVOT=PIVOT-PVKOL2(KPVRO)*QQ2(LPVCO)
      IF (IPIV.GE.3) THEN
      PIVOT=PIVOT-PVKOL3(KPVRO)*QQ3(LPVCO)
      END IF
      END IF
      END IF
C
C     PIVOTAL CHOICE OK
C
      IF (ABS(PIVOT).GT.RPIVOT) GO TO 40
C
C     CHECK NEXT DIAGONAL PIVOT,
C     IF NO OFF-DIAGONAL PIVOTS HAVE BEEN USED
C
      IF (IFULPV.EQ.1) GO TO 20
      DO 10 L=2,LC
      LPVCO=LPIV(L)
      KPVRO=KPIV(L)
      K0=L
      L0=L
      PIVOT=EQ(KPVRO,LPVCO)
      IF (ABS(PIVOT).LT.RPIVOT) GO TO 10
      IF (IPIV.GE.1) THEN
      PIVOT=PIVOT-PVKOL1(KPVRO)*QQ1(LPVCO)
      IF (IPIV.GE.2) THEN
      PIVOT=PIVOT-PVKOL2(KPVRO)*QQ2(LPVCO)
      IF (IPIV.GE.3) THEN
      PIVOT=PIVOT-PVKOL3(KPVRO)*QQ3(LPVCO)
      END IF
      END IF
      END IF
C
C     NEW PIVOT OK?
C
      IF (ABS(PIVOT).GT.RPIVOT) THEN
      NPIV=NPIV+1
      GO TO 40
      END IF
   10 CONTINUE
      END IF
C
C     OTHERWISE FIND A NEW PIVOT
C
   20 CONTINUE
      PIVOT=0.
      DO 30 L=1,LC
      DO 30 K=1,KR
      PIVA=EQ(KPIV(K),LPIV(L))
      IF (IPIV.GE.1) THEN
      PIVA=PIVA-PVKOL1(KPIV(K))*QQ1(LPIV(L))
      IF (IPIV.GE.2) THEN
      PIVA=PIVA-PVKOL2(KPIV(K))*QQ2(LPIV(L))
      IF (IPIV.GE.3) THEN
      PIVA=PIVA-PVKOL3(KPIV(K))*QQ3(LPIV(L))
      END IF
      END IF
      END IF
      IF (ABS(PIVA).GE.ABS(PIVOT)) THEN
      PIVOT=PIVA
      L0=L
      K0=K
      END IF
   30 CONTINUE
      LPVCO=LPIV(L0)
      KPVRO=KPIV(K0)
C
C     TEST THE NEW PIVOT
C
      IF (ABS(PIVOT).LT.SPIVOT) THEN
      WRITE (NOUT, 90) NELL,IKHED(KPVRO),ILHED(LPVCO),PIVOT
      CALL ERROR ('FNDPIV','SEE PREVIOUS MESSAGE',' ',0,' ',0,' ',' ',1)
      END IF
C
C     TAKE NOTE OF THE PIVOTING OPERATION
C
      IF (L0.NE.1 .OR. K0.NE.1) THEN
      NPIV=NPIV+1
      IFULPV=1
      END IF
C
C     REARRANGE THE POINTERS
C
   40 KRO=IKHED(KPVRO)
      LCO=ILHED(LPVCO)
      DO 50 L=L0+1,LC
   50 LPIV(L-1)=LPIV(L)
      DO 60 K=K0+1,KR
   60 KPIV(K-1)=KPIV(K)
      LC=LC-1
      KR=KR-1
C
C     UPDATE THE NUMBER OF EQUATIONS TO BE ELIMINATED
C
      IPIV=IPIV+1
C
C     UPDATE THE DETERMINANT
C
C      PIVTEM=PIVOT
C      ABSPIV=ABS(PIVOT)
C      IF (ABSPIV.GT.1.0E5) PIVTEM=PIVOT/ABSPIV
C      DET=DET*PIVTEM*(-1)**(KRO+LCO+NRP(KRO)+NCP(LCO))
C
      DO 70 IPERM=KRO+1,NUMDOF
   70 NRP(IPERM)=NRP(IPERM)-1
      DO 80 IPERM=LCO+1,NUMDOF
   80 NCP(IPERM)=NCP(IPERM)-1
      RETURN
C
   90 FORMAT (/,10X,'MATRIX IS SINGULAR OR ILL CONDITIONED ',/,10X,     
     1'ELEMENT = ',I5,5X,'KRO = ',I5,5X,'LCO = ',I5,5X,'PIVOT = ',E15.7)
      END
