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
      SUBROUTINE GENWRT (X,Y,ICON,LISTEL,MAP,ICONB)
C
C     *****************************************************************
C
C     SUBROUTINE TO CREATE THE 'GENESIS' FILE FOR AN INTERNALLY
C     GENERATED MESH
C
C     *****************************************************************
C
      CHARACTER*80 HED,CMMNT
C
      COMMON /HEADER/ HED,CMMNT(10)
      COMMON /TAPES/  NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7, 
     1                NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
      COMMON /SZDAT/  NUMEL,NUMNOD,NUMVAR,NUMDOF
      COMMON /GENDAT/ NELBLK,NUMNPS,LNPSNL,NUMESS,LESSEL,LESSNL,MXEBLK, 
     1                IBLK(100)
      COMMON /ELMDAT/ NNELM(6),NNCOR(6),NNMID(6),NNCTR(6),NSIDET(5),    
     1                NSIDEQ(7),NNSIDE(6,4,3)
C
      DIMENSION X(*), Y(*), ICON(NUMEL,*), LISTEL(*)
      DIMENSION MAP(*), ICONB(NUMEL,*)
C
      DATA NDIM,NVERSN,NATRIB/2,1,0/
      DATA DUMWRT/0./
C
C     ****************************************************************
C
C     TITLE RECORD
C
      WRITE (NTP12) HED
C
C     SIZING PARAMETER RECORD
C
      WRITE (NTP12) NUMNOD,NDIM,NUMEL,NELBLK,NUMNPS,LNPSNL,NUMESS,      
     1LESSEL,LESSNL,NVERSN
C
C     NODAL POINT COORDINATES RECORD
C
      WRITE (NTP12) (X(I),I=1,NUMNOD),(Y(I),I=1,NUMNOD)
C
C     ELEMENT ORDERING RECORD
C
      DO 10 I=1,NUMEL
      MAP(I)=I
   10 CONTINUE
      WRITE (NTP12) (MAP(I),I=1,NUMEL)
C
C     ELEMENT BLOCKS RECORD
C
      DO 40 K=1,NELBLK
      IDEBLK=K
      ITEST=IBLK(K)
      KIND=MOD(ITEST,100)
      NELNOD=NNELM(KIND)
      NUM=0
      DO 30 N=1,NUMEL
      IF (ITEST.EQ.LISTEL(N)) THEN
      NUM=NUM+1
      DO 20 I=1,NELNOD
      ICONB(NUM,I)=ICON(N,I)
   20 CONTINUE
      END IF
   30 CONTINUE
      NUMELB=NUM
C
      WRITE (NTP12) IDEBLK,NUMELB,NELNOD,NATRIB
C
      WRITE (NTP12) ((ICONB(N,I),I=1,NELNOD),N=1,NUMELB)
C
      WRITE (NTP12) DUMWRT
   40 CONTINUE
C
C     NODAL POINT BOUNDARY CONDITION RECORDS
C
      WRITE (NTP12) DUMWRT
      WRITE (NTP12) DUMWRT
      WRITE (NTP12) DUMWRT
C
      WRITE (NTP12) DUMWRT
      WRITE (NTP12) DUMWRT
C
C     ELEMENT SIDE BOUNDARY CONDITION DATA RECORDS
C
      WRITE (NTP12) DUMWRT
      WRITE (NTP12) DUMWRT
      WRITE (NTP12) DUMWRT
      WRITE (NTP12) DUMWRT
      WRITE (NTP12) DUMWRT
C
      WRITE (NTP12) DUMWRT
      WRITE (NTP12) DUMWRT
      WRITE (NTP12) DUMWRT
C
      RETURN
      END
