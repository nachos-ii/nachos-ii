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
      SUBROUTINE GENRD1
C
C     ******************************************************************
C
C     SUBROUTINE TO READ THE CRITICAL INPUT AND SIZING PARAMETERS
C     FROM THE 'GENESIS' FILE
C     'GENESIS' DATA IS REWRITTEN TO THE 'EXODUS' FILE
C
C     ******************************************************************
C
      CHARACTER*80 TITLE
C
      COMMON /TAPES/  NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7, 
     1                NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
      COMMON /SZDAT/  NUMEL,NUMNOD,NUMVAR,NUMDOF,NODSOL
      COMMON /LMTDAT/ MAXELM,MAXNOD,MAXBC,MAXEDG,MAXMAT,MAXNPT,MAXPLT,  
     1                MAXBRY,MAXBLK
      COMMON /GENDAT/ NELBLK,NUMNPS,LNPSNL,NUMESS,LESSEL,LESSNL,MXEBLK, 
     1                IBLK(100)
C
C     ******************************************************************
C
C     TITLE RECORD
C
      READ (NTP0) TITLE
      WRITE (NOUT, 20) TITLE
      WRITE (NTP12) TITLE
C
C     SIZING PARAMETER RECORD
C
      READ (NTP0) NUMNOD,NDIM,NUMEL,NELBLK,NUMNPS,LNPSNL,NUMESS,LESSEL, 
     1LESSNL,NVERSN
      WRITE (NTP12) NUMNOD,NDIM,NUMEL,NELBLK,NUMNPS,LNPSNL,NUMESS,      
     1LESSEL,LESSNL,NVERSN
      IF (NDIM.NE.2) CALL ERROR ('GENRD1','MESH DATA NOT SET FOR 2-D GEO
     1METRY','NO. OF DIMENSIONS',NDIM,' ',0,' ',' ',1)
      WRITE (NOUT, 30) NVERSN,NUMEL,NUMNOD,NELBLK,NUMNPS,NUMESS
C
C     NODAL COORDINATES RECORD
C
      READ (NTP0)
C
C     ELEMENT ORDER RECORD
C
      READ (NTP0)
C
C     ELEMENT BLOCKS RECORDS
C
      DO 10 K=1,NELBLK
      READ (NTP0) IDEBLK,NUMELB,NELNOD,NATRIB
      IF (NELNOD.GT.MAXNPT) MAXNPT=NELNOD
      IF (NUMELB.GT.MXEBLK) MXEBLK=NUMELB
      READ (NTP0)
      READ (NTP0)
   10 CONTINUE
      WRITE (NOUT, 40) MAXNPT
C
      RETURN
C
   20 FORMAT (3X,'DATA FROM MESH FILE - ',//,10X,'HEADING - ',A,/)
   30 FORMAT (10X,'VERSION NUMBER FOR GENESIS FORMAT -     ',I7,/,10X,  
     1'NUMBER OF ELEMENTS IN MESH -            ',I7,/,10X,              
     2'NUMBER OF NODES IN MESH -               ',I7,/,10X,              
     3'NUMBER OF ELEMENT BLOCKS -              ',I7,/,10X,              
     4'NUMBER OF BOUNDARY CONDITIONS (NODES) - ',I7,/,10X,              
     5'NUMBER OF BOUNDARY CONDITIONS (SIDES) - ',I7,/)
   40 FORMAT (10X,'MAXIMUM NUMBER OF NODES/ELEMENT -       ',I7,/)
      END
