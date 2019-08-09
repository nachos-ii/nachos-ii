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
      SUBROUTINE GENRD2 (X,Y,ICON,LISTEL,MAP,ICONB)
C
C     ******************************************************************
C
C     SUBROUTINE TO READ AND STORE MESH DATA FROM THE 'GENESIS' FILE
C     'GENESIS' DATA IS REWRITTEN TO THE 'EXODUS' FILE
C
C     ******************************************************************
C
      COMMON /TAPES/  NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7, 
     1                NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
      COMMON /SZDAT/  NUMEL,NUMNOD,NUMVAR,NUMDOF,NODSOL
      COMMON /LMTDAT/ MAXELM,MAXNOD,MAXBC,MAXEDG,MAXMAT,MAXNPT,MAXPLT,  
     1                MAXBRY,MAXBLK
      COMMON /GENDAT/ NELBLK,NUMNPS,LNPSNL,NUMESS,LESSEL,LESSNL,MXEBLK, 
     1                IBLK(100)
C
      DIMENSION X(*), Y(*), ICON(NUMEL,*)
      DIMENSION LISTEL(*), MAP(*), ICONB(MAXNPT,*)
      DIMENSION KNDLST(9)
C
      DATA (KNDLST(I),I=1,9)/0,0,0,0,0,1,0,3,5/
      DATA DUMWRT/0./
C
C     ******************************************************************
C
C     TITLE RECORD
C
      READ (NTP0)
C
C     SIZING PARAMETERS RECORD
C
      READ (NTP0)
C
C     NODAL POINT COORDINATES RECORD
C
      READ (NTP0) (X(I),I=1,NUMNOD),(Y(I),I=1,NUMNOD)
      WRITE (NTP12) (X(I),I=1,NUMNOD),(Y(I),I=1,NUMNOD)
C
C     ELEMENT ORDERING RECORD, INVERT MAPPING
C
      READ (NTP0) (LISTEL(I),I=1,NUMEL)
      WRITE (NTP12) (LISTEL(I),I=1,NUMEL)
C
      DO 10 I=1,NUMEL
      MAP(LISTEL(I))=I
   10 CONTINUE
C
C     ELEMENT BLOCKS RECORDS
C
      KUREL=1
      DO 30 K=1,NELBLK
      READ (NTP0) IDEBLK,NUMELB,NELNOD,NATRIB
      WRITE (NTP12) IDEBLK,NUMELB,NELNOD,NATRIB
C
      READ (NTP0) ((ICONB(J,I),J=1,NELNOD),I=1,NUMELB)
      WRITE (NTP12) ((ICONB(J,I),J=1,NELNOD),I=1,NUMELB)
C
      READ (NTP0)
      WRITE (NTP12) DUMWRT
C
C     CHECK ELEMENT TYPE AND STORE NODE ORDERING FOR CONNECTIVITY ARRAY
C
      KIND=KNDLST(NELNOD)
      IF (KIND.EQ.0) CALL ERROR ('GENRD2','SPECIFIED ELEMENT TYPE NOT S 
     1UPPORTED','NUMBER OF NODES',NELNOD,'ELEMENT BLOCK ID',IDEBLK,' ', 
     2' ',1)
C
      DO 30 I=1,NUMELB
      IEL=MAP(KUREL)
      LISTEL(IEL)=-IDEBLK
      DO 20 J=1,NELNOD
      ICON(IEL,J)=ICONB(J,I)
   20 CONTINUE
      KUREL=KUREL+1
   30 CONTINUE
C
      RETURN
      END
