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
      SUBROUTINE FORMKF
C     ******************************************************************
C
C     SUBROUTINE TO DIRECT THE COMPUTATION AND STORAGE OF THE
C     ELEMENT LEVEL MATRICES
C
C     ******************************************************************
C
      REAL MASS,KUV,KP,KPU,KT,KTQ,KTR,KV1,KV2
C
      COMMON /TAPES/  NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7, 
     1                NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
      COMMON /SZDAT/  NUMEL,NUMNOD,NUMVAR,NUMDOF,NODSOL
      COMMON /PRBDAT/ IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,    
     1                IPNLTY,PNLTY
      COMMON /GAUSS2/ GSPT3(3),GSWT3(3)
      COMMON /HAMMR2/ HMPT1(7),HMPT2(7),HMWT(7)
      COMMON /TRI6/   F6Q(6,7),F6L(3,7),DF6QDS(6,7),DF6QDT(6,7)
      COMMON /QUAD8/  F8Q(8,9),F8L(4,9),DF8QDS(8,9),DF8QDT(8,9)
      COMMON /QUAD9/  F9Q(9,9),F9L(4,9),DF9QDS(9,9),DF9QDT(9,9)
C
      DIMENSION X(9), Y(9), ICONEL(9), NBC(30), IBC(30), BCVAL(30)
      DIMENSION KUV(22,22,4), CX(9,9,9), CY(9,9,9), C(9,9,9), MASS(9,9)
      DIMENSION FUV(22), FBDY(9,9,4), KP(18,18), KPU(3,18)
      DIMENSION KT(9,9,4), KTQ(9,9), KTR(9,3,9)
      DIMENSION FT(9), FTQ(9), FTR(9,3), FTS(9,4)
      DIMENSION KV1(9,9,4), KV2(9,9,4)
      DIMENSION FV1(9), FV2(9)
C
      DATA NGAUSS/3/,NHAM/7/
C
C     ******************************************************************
C
C     EVALUATE SHAPE FUNCTIONS AND THEIR DERIVATIVES AT THE INTEGRATION
C     POINTS
C
      CALL SHAPE6 (NHAM,HMPT1,HMPT2)
      CALL SHAPE8 (NGAUSS,GSPT3,GSPT3)
      CALL SHAPE9 (NGAUSS,GSPT3,GSPT3)
C
C     LOOP ON ELEMENTS
C
      ISTOP=0
      DO 50 LL=1,NUMEL
      READ (NTP2) KIND,MAT,NN,(X(I),I=1,NN),(Y(I),I=1,NN),(ICONEL(I),   
     1I=1,NN),NUMBC,(NBC(I),IBC(I),BCVAL(I),I=1,NUMBC)
C
C     CALCULATE ELEMENT DIFFUSION, CONVECTION AND MASS MATRICES
C
      IMAP=MOD(KIND,2)
      GO TO (10, 10, 20, 20, 30, 30), KIND
C
C     SUB AND ISOPARAMETRIC (SIX NODE) TRIANGLE (TRI6/3-TRI6/6)
C
   10 CONTINUE
      CALL TRIB (IMAP,MAT,NHAM,HMPT1,HMPT2,HMWT,X,Y,KUV,CX,CY,MASS,FBDY,
     1KT,FTS,KV1,KV2,KP,KPU,C)
      NSZ1=15
      NSZ2=6
      NSZ3=3
      NSZ4=3
      NSZ5=12
      GO TO 40
C
C     SUB AND ISOPARAMETRIC (EIGHT NODE) QUADRILATERAL (QUAD8/4-QUAD8/8)
C
   20 CONTINUE
      CALL QUADB (IMAP,MAT,NGAUSS,GSPT3,GSPT3,GSWT3,X,Y,KUV,CX,CY,MASS, 
     1FBDY,KT,FTS,KV1,KV2,KP,KPU,C)
      NSZ1=20
      NSZ2=8
      NSZ3=4
      NSZ4=3
      NSZ5=16
      GO TO 40
C
C     SUB AND ISOPARAMETRIC (NINE NODE) QUADRILATERAL (QUAD9/4-QUAD9/9)
C
   30 CONTINUE
      CALL QUADC (IMAP,MAT,NGAUSS,GSPT3,GSPT3,GSWT3,X,Y,KUV,CX,CY,MASS, 
     1FBDY,KT,FTS,KV1,KV2,KP,KPU,C)
      NSZ1=22
      NSZ2=9
      NSZ3=4
      NSZ4=3
      NSZ5=18
C
   40 CONTINUE
      IF (IMAP.LT.0) THEN
      WRITE (NOUT, 60) LL
      ISTOP=1
      END IF
C
C     APPLY BOUNDARY CONDITIONS
C
      CALL BNDRY (LL,KIND,NN,NUMBC,NBC,IBC,BCVAL,NGAUSS,GSPT3,GSWT3,    
     1X,Y,KUV,FUV,KT,FT,KTQ,FTQ,KTR,FTR,KV1,FV1,KV2,FV2,IBCPT1,IBCPT2,  
     2IBCPT3)
C
C     WRITE ELEMENT DATA TO TAPE
C
      WRITE (NTP4) NSZ1,NSZ2,NSZ3,NSZ4,NSZ5,IBCPT1,IBCPT2,IBCPT3,       
     1KUV,CX,CY,MASS,KP,C,FUV,FBDY,KT,KTQ,KTR,FT,FTQ,FTR,FTS,           
     2KV1,FV1,KV2,FV2
C
      IF (IPNLTY.EQ.1) THEN
      WRITE (NTP3) NSZ4,NSZ5,KPU
      END IF
   50 CONTINUE
      IF (ISTOP.GT.0) CALL ERROR ('FORMKF','BAD ELEMENT JACOBIAN',      
     1' ',0,' ',0,' ',' ',1)
C
   60 FORMAT (//,10X,'ELEMENT NO.' I5,' HAS A NEGATIVE JACOBIAN' )
      END
