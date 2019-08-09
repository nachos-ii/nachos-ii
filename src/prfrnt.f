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
      SUBROUTINE PRFRNT (NOP,LISTEL,NCN,MDF,NOPP,NCHECK)
C
C     ******************************************************************
C
C     SUBROUTINE TO SET UP PARAMETERS AND POINTER VECTORS FOR FRONTAL
C     ELIMINATION PROCEDURE (PRE-FRONT PROCEDURE)
C
C     ORIGINAL CODE BY P.HOOD, UNIVERSITY OF WALES, SWANSEA,WALES, 1976.
C     PRESENT VERSION BY R.E.BENNER, SANDIA NATIONAL LABS, 1985.
C
C     ARRAY USAGE - NCN(NUMEL):        NUMBER OF DOF IN EACH ELEMENT
C                   NOP(NUMEL,MAXNPT): CONNECTIVITY, SIGN DETERMINES
C                                      LAST APPEARANCE OF NODE
C                   MDF(NUMNOD):       NUMBER OF DOF AT EACH NODE
C                   NOPP (NUMNOD):     NUMBER OF FIRST DOF AT EACH NODE
C                   NCHECK(NUMDOF):    SCRATCH ARRAY TO CHECK STATUS OF
C                                      EACH DOF
C     ******************************************************************
C
      COMMON /TAPES/  NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7, 
     1                NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
      COMMON /SZDAT/  NUMEL,NUMNOD,NUMVAR,NUMDOF,NODSOL
      COMMON /PRBDAT/ IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,    
     1                IPNLTY,PNLTY
      COMMON /ELMDAT/ NNELM(6),NNCOR(6),NNMID(6),NNCTR(6),NSIDET(5),    
     1                NSIDEQ(7),NNSIDE(6,4,3)
      COMMON /FRNT/   RPIVOT,SPIVOT,IRESOL,ISOLVE,NSUM,IPIVOT,MWGA
      COMMON /FRNTSZ/ MXFRNT,MXDOFE,MXBLKS,MXBUFF,LIWK,LRWK
C
      DIMENSION NOP(NUMEL,*), NCN(*), MDF(*), NOPP(*), NCHECK(*)
      DIMENSION LISTEL(*)
C
C     ******************************************************************
C
C     SET UP DATA VECTORS (NOP,MDF,NOPP)
C
      MXDOFE=0
      DO 30 N=1,NUMEL
      NN=ABS(LISTEL(N))
      KIND=MOD(NN,100)
      NCOR=NNCOR(KIND)
      NMID=NNMID(KIND)
      NCNTR=NNCTR(KIND)
      NN=NNELM(KIND)
      NDOFEL=0
      NVEL=2
      NPRESS=1
      IF (IPNLTY.EQ.1.OR.IPFUNC.EQ.1) NPRESS=0
      NTEMP=IFORCE+IFREE
      NVAR=IVAR1+IVAR2
C
      KDOF=NVEL+NPRESS+NTEMP+NVAR
      DO 10 I=1,NCOR
      INODE=NOP(N,I)
      MDF(INODE)=KDOF
      NDOFEL=NDOFEL+KDOF
   10 CONTINUE
      KDOF=NVEL+NTEMP+NVAR
      DO 20 I=1,NMID
      INODE=NOP(N,I+NCOR)
      MDF(INODE)=KDOF
      NDOFEL=NDOFEL+KDOF
   20 CONTINUE
      IF (NCNTR.EQ.1) THEN
      KDOF=NVEL+NTEMP+NVAR
      INODE=NOP(N,NN)
      MDF(INODE)=KDOF
      NDOFEL=NDOFEL+KDOF
      END IF
      IF (IPFUNC.EQ.1) THEN
      KDOF=3
      INODE=NOP(N,NN+1)
      MDF(INODE)=KDOF
      NDOFEL=NDOFEL+KDOF
      END IF
      NCN(N)=NDOFEL
      IF (NDOFEL.GT.MXDOFE) MXDOFE=NDOFEL
   30 CONTINUE
C
      KDOF=1
      DO 40 I=1,NODSOL
      NOPP(I)=KDOF
      KDOF=KDOF+MDF(I)
   40 CONTINUE
      NUMDOF=KDOF-1
C
C     CHECK UNIQUENESS OF EACH DEGREE OF FREEDOM (D.O.F.) IN
C     EACH ELEMENT
C
      DO 50 I=1,NUMEL
      NN=ABS(LISTEL(I))
      KIND=MOD(NN,100)
      NBN=NNELM(KIND)+IPFUNC
      DO 50 J=1,NBN-1
      J0=NOP(I,J)
      DO 50 K=1,MDF(J0)
      DO 50 L=J+1,NBN
      L0=NOP(I,L)
      DO 50 M=1,MDF(L0)
      IF (NOPP(J0)+K .EQ. NOPP(L0)+M) THEN
      WRITE (NOUT, 110) I,J,L,NOP(I,J),NOP(I,L),NOPP(J0),NOPP(L0),K,M
      CALL ERROR ('PRFRNT','TWO DOF WITH SAME NODE NO.',' ',0,' ',0,    
     1' ',' ',1)
      END IF
   50 CONTINUE
C
C     FIND LAST APPEARANCE OF EACH NODE.  SIGNAL LAST APPEARANCE
C     OF A NODE BY CHANGING THE SIGN OF THE CORRESPONDING ENTRY
C     IN "NOP".  SEARCH BEGINS WITH THE LAST ELEMENT, AND THE FIRST
C     APPEARANCE OF EACH NODE IS SOUGHT.
C
      DO 60 I=1,NUMDOF
      NCHECK(I)=0
   60 CONTINUE
      DO 70 I=1,NUMEL
      NUMELP=NUMEL+1-I
      NN=ABS(LISTEL(NUMELP))
      KIND=MOD(NN,100)
      NBN=NNELM(KIND)+IPFUNC
      DO 70 J=1,NBN
      K=NOP(NUMELP,J)
      IF (NCHECK(K).EQ.0) THEN
      NCHECK(K)=1
      NOP(NUMELP,J)=-NOP(NUMELP,J)
      END IF
   70 CONTINUE
C
C     CALCULATE FRONTWIDTH "MXFRNT"
C
      DO 80 I=1,NUMDOF
      NCHECK(I)=0
   80 CONTINUE
      NELIM=0
      NAPP=0
      MXFRNT=0
      NMAXLC=0
      DO 100 I=1,NUMEL
      NNEW=0
      NNEWEL=0
      NN=ABS(LISTEL(I))
      KIND=MOD(NN,100)
      NBN=NNELM(KIND)+IPFUNC
      DO 90 J=1,NBN
      K=NOP(I,J)
      KABS=ABS(K)
      IF (NCHECK(KABS).EQ.1) THEN
      IF (K.LE.0) THEN
      NCHECK(KABS)=0
      NNEWEL=NNEWEL+MDF(KABS)
      END IF
       ELSE
      NNEW=NNEW+MDF(KABS)
      IF (K.LE.0) THEN
      NCHECK(KABS)=0
      NNEWEL=NNEWEL+MDF(KABS)
       ELSE
      NCHECK(KABS)=1
      END IF
      END IF
   90 CONTINUE
      NMAXLC=NAPP+NNEW-NELIM+NSUM
      NAPP=NAPP+NNEW
      NELIM=NELIM+NNEWEL
      IF (MXFRNT.LT.NMAXLC) THEN
      MXFRNT=NMAXLC
      NELM=I
      END IF
  100 CONTINUE
C
      RETURN
C
  110 FORMAT (/2X,'PREFRONT ERROR:  TWO D.O.F. WITH SAME NODE NO.',     
     1  2X,'I=',I4,2X,'J=',I4,2X,'L=',I4,/,50X,'NOP(I,J)=',I5,2X,       
     2  'NOP(I,L)=',I5,2X,'NOPP(J0)=',I5,2X,'NOPP(L0)=',I5,2X,'K=',I3,  
     3  2X,'M=',I3)
      END
