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
      SUBROUTINE CNNCT (ICON,IJK)
C
C     ******************************************************************
C
C     SUBROUTINE TO CONSTRUCT A CONNECTIVITY ARRAY AND ASSIGN GLOBAL
C     NODAL POINT NUMBERS
C
C     ******************************************************************
C
      COMMON /TAPES/  NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7, 
     1                NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
      COMMON /SZDAT/  NUMEL,NUMNOD,NUMVAR,NUMDOF,NODSOL
      COMMON /LMTDAT/ MAXELM,MAXNOD,MAXBC,MAXEDG,MAXMAT,MAXNPT,MAXPLT,  
     1                MAXBRY,MAXBLK
C
      DIMENSION ICON(NUMEL,*), IJK(*)
      DIMENSION X(9), Y(9), IDUM(9)
      DIMENSION NBC(30), IBC(30), BCVAL(30)
C
C     ******************************************************************
C
C     INITIALIZE ARRAYS
C
      REWIND (NTP2)
      REWIND (NTP11)
      DO 10 I=1,MAXNOD
      IJK(I)=0
   10 CONTINUE
      DO 20 I=1,NUMEL
      DO 20 J=1,MAXNPT
      ICON(I,J)=0
   20 CONTINUE
C
C     LOOP ON ELEMENTS, PACK NODE IDENTIFIERS (NICNAMES)
C
      DO 30 N=1,NUMEL
      READ (NTP2) KIND,MAT,NN,(X(I),I=1,NN),(Y(I),I=1,NN),              
     1(ICON(N,I),I=1,NN),NUMBC,(NBC(I),IBC(I),BCVAL(I),I=1,NUMBC)
      WRITE (NTP11) KIND,MAT,NN,(X(I),I=1,NN),(Y(I),I=1,NN),            
     1(ICON(N,I),I=1,NN),NUMBC,(NBC(I),IBC(I),BCVAL(I),I=1,NUMBC)
   30 CONTINUE
C
      NNODES=1
      DO 60 I=1,NUMEL
      DO 50 J=1,MAXNPT
      IJCHCK=ICON(I,J)
      IF (IJCHCK.EQ.0) GO TO 50
      DO 40 K=1,NNODES
      IJTEST=IJK(K)
      IF (IJCHCK.EQ.IJTEST) GO TO 50
   40 CONTINUE
      IF (NNODES.GT.MAXNOD) CALL ERROR ('CNNCT','MAXIMUM NUMBER OF NODES
     1 EXCEEDED','CURRENT ELEMENT NUMBER',I,'NODES ALLOWED',MAXNOD,' ', 
     2' ',1)
      IJK(NNODES)=IJCHCK
      NNODES=NNODES+1
   50 CONTINUE
   60 CONTINUE
      NNODES=NNODES-1
      NUMNOD=NNODES
C
C     REPLACE NODE IDENTIFIERS WITH GLOBAL NODE NUMBERS
C
      DO 100 I=1,NUMEL
      DO 90 J=1,MAXNPT
      IJCHCK=ICON(I,J)
      IF (IJCHCK.EQ.0) GO TO 90
      DO 70 K=1,NNODES
      IJTEST=IJK(K)
      IF (IJTEST.EQ.IJCHCK) GO TO 80
   70 CONTINUE
      CALL ERROR ('CNNCT','INVALID NODAL POINT NUMBER','CURRENT ELEMENT 
     1NUMBER',I,' ',0,' ',' ',1)
   80 CONTINUE
      ICON(I,J)=K
   90 CONTINUE
  100 CONTINUE
C
C     REWRITE NTP2 WITH NEW CONNECTIVITY
C
      REWIND (NTP2)
      REWIND (NTP11)
      DO 110 N=1,NUMEL
      READ (NTP11) KIND,MAT,NN,(X(I),I=1,NN),(Y(I),I=1,NN),             
     1(IDUM(I),I=1,NN),NUMBC,(NBC(I),IBC(I),BCVAL(I),I=1,NUMBC)
      WRITE (NTP2) KIND,MAT,NN,(X(I),I=1,NN),(Y(I),I=1,NN),             
     1(ICON(N,I),I=1,NN),NUMBC,(NBC(I),IBC(I),BCVAL(I),I=1,NUMBC)
  110 CONTINUE
      WRITE (NOUT, 120) NUMNOD
      RETURN
C
  120 FORMAT (//,3X,'CONNECTIVITY ARRAY CONSTRUCTED FOR ',I7,' GLOBAL NO
     1DES',//)
      END
