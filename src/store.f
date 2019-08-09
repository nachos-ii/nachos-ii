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
      SUBROUTINE STORE (X,Y,ICON,LISTEL,LISTND)
C
C     ******************************************************************
C
C     SUBROUTINE TO STORE BASIC ELEMENT AND NODE DATA
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
      COMMON /MATDAT/ NMAT,NPROP,PROP(15,10),NXMAT,NXPROP,XPROP(15,10)
      COMMON /PRBDAT/ IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,    
     1                IPNLTY,PNLTY
C
      DIMENSION X(*), Y(*), ICON(NUMEL,*), LISTEL(*), LISTND(*)
      DIMENSION XX(9), YY(9)
C
C     ******************************************************************
C
      REWIND (NTP2)
C
C     STORE COORDINATES AND ELEMENT TYPE
C
      DO 10 I=1,NUMNOD
      LISTND(I)=0
   10 CONTINUE
      DO 50 N=1,NUMEL
      READ (NTP2) KIND,MAT,NN,(XX(I),I=1,NN),(YY(I),I=1,NN),            
     1(ICON(N,I),I=1,NN)
      IF (IPFUNC.EQ.1) ICON(N,NN+1)=NUMNOD+N
      DO 20 I=1,NN
      INODE=ICON(N,I)
      X(INODE)=XX(I)
      Y(INODE)=YY(I)
   20 CONTINUE
      IF (PROP(15,MAT).EQ.3.0) GO TO 40
      DO 30 I=1,NN
      INODE=ICON(N,I)
      LISTND(INODE)=LISTND(INODE)+1
   30 CONTINUE
   40 CONTINUE
      LISTEL(N)=MAT*100+KIND
   50 CONTINUE
C
C     RECORD ELEMENT BLOCK DATA FOR USE IN 'GENESIS' FILE
C
      NBLK=1
      IBLK(1)=LISTEL(1)
      DO 80 N=1,NUMEL
      ITEST=LISTEL(N)
      NB=NBLK
      DO 60 I=1,NB
      IF (ITEST.EQ.IBLK(I)) GO TO 70
   60 CONTINUE
      NBLK=NBLK+1
      IF (NBLK.GT.MAXBLK) CALL ERROR ('STORE','MAXIMUM NUMBER OF ELEMENT
     1 BLOCKS EXCEEDED','MAX BLOCKS',MAXBLK,'NUMBER OF BLOCKS',NBLK,    
     2'INCREASE MAXBLK','AND IBLK(MAXBLK)',1)
      IBLK(NBLK)=ITEST
   70 CONTINUE
   80 CONTINUE
      NELBLK=NBLK
C
      RETURN
      END
