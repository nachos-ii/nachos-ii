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
      SUBROUTINE MESH (XN,YN,IBND,CC,HH,CR,HR,CA,HA,NN,NM,IPRINT)
C
C     ******************************************************************
C
C     PROGRAM TO GENERATE AN ARRAY OF NODAL POINTS USING THE
C     ISOPARAMETRIC SHAPE FUNCTIONS
C
C     ******************************************************************
C
      CHARACTER*10 CDATA,WORD
C
      COMMON /INDATR/ RDATA(150)
      COMMON /INDATI/ IDATA(150)
      COMMON /INDATC/ CDATA(150)
      COMMON /TAPES/  NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7, 
     1                NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
      COMMON /LMTDAT/ MAXELM,MAXNOD,MAXBC,MAXEDG,MAXMAT,MAXNPT,MAXPLT,  
     1                MAXBRY,MAXBLK
C
      DIMENSION XN(NN,*), YN(NN,*), IBND(NN,*)
      DIMENSION CC(*), HH(*), CR(*), HR(*), CA(*), HA(*)
      DIMENSION AN(12), X(12), Y(12), XI(9), YI(9), A(5), DX(4), DY(4)
C
      DATA DUM/-1.23E-21/
C
C     ******************************************************************
C
      NP=0
      IEXT=0
      DO 10 I=1,NN
      DO 10 J=1,NM
      XN(I,J)=DUM
      YN(I,J)=DUM
   10 CONTINUE
   20 CONTINUE
      NP=NP+1
      IREF=0
      ILAP=0
      CALL RDFFLD (NIN,RDATA,IDATA,CDATA)
      IF (CDATA(1).EQ.'END       ') RETURN
      IF (CDATA(1).EQ.'TBLOCK    ') GO TO 220
      IF (CDATA(1).EQ.'QBLOCK    ') GO TO 120
      IF (CDATA(1).EQ.'POINT     ') GO TO 30
      IF (CDATA(1).EQ.'ARC       ') GO TO 40
      IF (CDATA(1).EQ.'COPY      ') GO TO 330
      IF (CDATA(1).EQ.'REFLECT   ') GO TO 330
      IF (CDATA(1).EQ.'FILLIN    ') GO TO 400
      IF (CDATA(1).EQ.'EXTDEF    ') GO TO 390
      CALL ERROR ('MESH','UNRECOGNIZED MESH COMMAND',' ',0,' ',0,'WORD',
     1CDATA(1),1)
C
C     READ DATA FOR A POINT
C
   30 CONTINUE
      NCOORD=1
      IMIN=IDATA(2)
      IMAX=IMIN
      JMIN=IDATA(3)
      JMAX=JMIN
      WORD=CDATA(6)
      XO=RDATA(7)
      YO=RDATA(8)
      X(1)=RDATA(4)
      Y(1)=RDATA(5)
      IF (IMIN.GT.NN) CALL ERROR ('MESH','INPUT I VALUE EXCEEDS SPECIFIE
     1D MESH IMAX','MESH PART',NP,'I VALUE',IMIN,'MESH COMMAND','POINT',
     21)
      IF (JMIN.GT.NM) CALL ERROR ('MESH','INPUT J VALUE EXCEEDS SPECIFIE
     1D MESH JMAX','MESH PART',NP,'J VALUE',JMIN,'MESH COMMAND','POINT',
     21)
      IF (WORD.EQ.'POLAR') CALL TRNSFM (X,Y,XO,YO,NCOORD)
      XN(IMIN,JMIN)=X(1)
      YN(IMIN,JMIN)=Y(1)
      GO TO 520
C
C     READ DATA FOR AN ARC
C
   40 CONTINUE
      NCOORD=4
      IMIN=IDATA(2)
      JMIN=IDATA(3)
      IMAX=IDATA(4)
      JMAX=IDATA(5)
      G1=RDATA(6)
      WORD=CDATA(7)
      XO=RDATA(8)
      YO=RDATA(9)
      NI=IMAX-IMIN
      NJ=JMAX-JMIN
      IF (ABS(NI-NJ).NE.ABS(NI+NJ)) CALL ERROR ('MESH','INCORRECT I,J SP
     1ECIFICATION FOR AN ARC','MESH PART',NP,' ',0,' ',' ',1)
      IF (IMIN.GT.NN.OR.IMAX.GT.NN) CALL ERROR ('MESH','INPUT I VALUE EX
     1CEEDS SPECIFIED MESH IMAX','MESH PART',NP,' ',0,'MESH COMMAND',   
     2'ARC',1)
      IF (JMIN.GT.NM.OR.JMAX.GT.NM) CALL ERROR ('MESH','INPUT J VALUE EX
     1CEEDS SPECIFIED MESH JMAX','MESH PART',NP,' ',0,'MESH COMMAND',   
     2'ARC',1)
      IF (IMIN.EQ.IMAX) N1=JMAX-JMIN+1
      IF (JMIN.EQ.JMAX) N1=IMAX-IMIN+1
      IF (G1.LE.0.) G1=1.0
      IF (N1.GE.MAXEDG) CALL ERROR ('MESH','TOO MANY NODES SPECIFIED ON 
     1ONE ARC','NODES SPECIFIED',N1,'NODES ALLOWED',MAXEDG,' ',' ',1)
      CALL SEG1 (-1.,1.,1.,1.,G1,CR,HR,N1)
      CALL RDFFLD (NIN,RDATA,IDATA,CDATA)
      DO 50 I=1,NCOORD
      X(I)=RDATA(I)
   50 CONTINUE
      CALL RDFFLD (NIN,RDATA,IDATA,CDATA)
      DO 60 I=1,NCOORD
      Y(I)=RDATA(I)
   60 CONTINUE
      IF (WORD.EQ.'POLAR') CALL TRNSFM (X,Y,XO,YO,NCOORD)
      N1=1
      DO 110 I=IMIN,IMAX
      DO 110 J=JMIN,JMAX
      XN(I,J)=0.
      YN(I,J)=0.
      S=CR(N1)
      IF (X(3).EQ.0.) GO TO 80
      IF (X(4).EQ.0.) GO TO 70
      AN(1)=1./16.*(1.-S)*(9.*S*S-1.)
      AN(2)=1./16.*(1.+S)*(9.*S*S-1.)
      AN(3)=9./16.*(1.-S*S)*(1.-3.*S)
      AN(4)=9./16.*(1.-S*S)*(1.+3.*S)
      IFLAG=4
      GO TO 90
   70 CONTINUE
      AN(1)=.5*S*(S-1.)
      AN(2)=.5*S*(S+1.)
      AN(3)=1.-S*S
      IFLAG=3
      GO TO 90
   80 CONTINUE
      AN(1)=.5*(1.-S)
      AN(2)=.5*(1.+S)
      IFLAG=2
   90 CONTINUE
      DO 100 K=1,IFLAG
      XN(I,J)=XN(I,J)+AN(K)*X(K)
      YN(I,J)=YN(I,J)+AN(K)*Y(K)
  100 CONTINUE
      N1=N1+1
  110 CONTINUE
      GO TO 520
C
C     READ DATA FOR A QUADRILATERAL BLOCK
C
  120 CONTINUE
      NCOORD=12
      IMIN=IDATA(2)
      JMIN=IDATA(3)
      IMAX=IDATA(4)
      JMAX=IDATA(5)
      G1=RDATA(6)
      G2=RDATA(7)
      G3=RDATA(8)
      G4=RDATA(9)
      WORD=CDATA(10)
      XO=RDATA(11)
      YO=RDATA(12)
      IF (G1.EQ.0.) G1=1.
      IF (G2.EQ.0.) G2=1.
      IF (G3.EQ.0.) G3=1.
      IF (G4.EQ.0.) G4=1.
      IF (IMIN.GT.NN.OR.IMAX.GT.NN) CALL ERROR ('MESH','INPUT I VALUE EX
     1CEEDS SPECIFIED MESH IMAX','MESH PART',NP,' ',0,'MESH COMMAND',   
     2'QBLOCK',1)
      IF (JMIN.GT.NM.OR.JMAX.GT.NM) CALL ERROR ('MESH','INPUT J VALUE EX
     1CEEDS SPECIFIED MESH JMAX','MESH PART',NP,' ',0,'MESH COMMAND',   
     2'QBLOCK',1)
      CALL RDFFLD (NIN,RDATA,IDATA,CDATA)
      DO 130 I=1,NCOORD
      X(I)=RDATA(I)
  130 CONTINUE
      CALL RDFFLD (NIN,RDATA,IDATA,CDATA)
      DO 140 I=1,NCOORD
      Y(I)=RDATA(I)
  140 CONTINUE
      IF (WORD.EQ.'POLAR') CALL TRNSFM (X,Y,XO,YO,NCOORD)
      XN(IMIN,JMIN)=X(1)
      YN(IMIN,JMIN)=Y(1)
      IF (IMAX.LE.0) GO TO 520
C
C     CALCULATE THE MESH POINT COORDINATES ALONG SIDES 1 AND 3
C     IN THE PSI-ETA PLANE, STORE IN C AND H VECTORS
C
      NI=IMAX-IMIN+1
      NJ=JMAX-JMIN+1
      IF (NI.GE.MAXEDG) CALL ERROR ('MESH','TOO MANY NODES ON EDGE OF O 
     1NE QBLOCK','NODES SPECIFIED',NI,'NODES ALLOWED',MAXEDG,' ',' ',1)
      IF (NJ.GE.MAXEDG) CALL ERROR ('MESH','TOO MANY NODES ON EDGE OF O 
     1NE QBLOCK','NODES SPECIFIED',NJ,'NODES ALLOWED',MAXEDG,' ',' ',1)
      CALL SEG1 (-1.,-1.,1.,-1.,G1,CC,HH,NI)
      CALL SEG1 (-1.,1.,1.,1.,G3,CR,HR,NI)
C
C     GIVEN C AND H,  GENERATE J=CONSTANT LINES FROM IMIN TO IMAX
C
      DO 210 I=IMIN,IMAX
      KI=I-IMIN+1
      G=G1
      IF (NI.GT.1) G=G4-(G4-G2)/(NI-1)*(KI-1)
      CALL SEG1 (CC(KI),HH(KI),CR(KI),HR(KI),G,CA,HA,NJ)
      DO 210 J=JMIN,JMAX
      KJ=J-JMIN+1
      PSI=CA(KJ)
      ETA=HA(KJ)
C
C     CALCULATE XN AND YN
C
      XN(I,J)=0.0
      YN(I,J)=0.0
      AN(1)=1./4.*(1.-PSI)*(1.-ETA)
      AN(2)=1./4.*(1.+PSI)*(1.-ETA)
      AN(3)=1./4.*(1.+PSI)*(1.+ETA)
      AN(4)=1./4.*(1.-PSI)*(1.+ETA)
      IFLAG=1
      DO 150 KK=5,12
      II=17-KK
      AN(II)=0.0
      IF (X(II).NE.0.0) GO TO 160
  150 CONTINUE
      GO TO 170
  160 CONTINUE
      CALL FILLQ (X,Y,IFLAG)
  170 CONTINUE
      IF (IFLAG-2) 200, 180, 190
  180 CONTINUE
      AN(1)=AN(1)*(-PSI-ETA-1.)
      AN(2)=AN(2)*(PSI-ETA-1.)
      AN(3)=AN(3)*(PSI+ETA-1.)
      AN(4)=AN(4)*(-PSI+ETA-1.)
      AN(5)=1./2.*(1.-PSI**2)*(1.-ETA)
      AN(6)=1./2.*(1.-ETA**2)*(1.+PSI)
      AN(7)=1./2.*(1.-PSI**2)*(1.+ETA)
      AN(8)=1./2.*(1.-ETA**2)*(1.-PSI)
      GO TO 200
  190 CONTINUE
      AN(1)=AN(1)*(1./8.*(-10.+9.*(PSI**2+ETA**2)))
      AN(2)=AN(2)*(1./8.*(-10.+9.*(PSI**2+ETA**2)))
      AN(3)=AN(3)*(1./8.*(-10.+9.*(PSI**2+ETA**2)))
      AN(4)=AN(4)*(1./8.*(-10.+9.*(PSI**2+ETA**2)))
      AN(5)=9./32.*(1.-PSI**2)*(1.-ETA)*(1.-3.*PSI)
      AN(6)=9./32.*(1.-ETA**2)*(1.+PSI)*(1.-3.*ETA)
      AN(7)=9./32.*(1.-PSI**2)*(1.+ETA)*(1.+3.*PSI)
      AN(8)=9./32.*(1.-ETA**2)*(1.-PSI)*(1.+3.*ETA)
      AN(9)=9./32.*(1.-PSI**2)*(1.-ETA)*(1.+3.*PSI)
      AN(10)=9./32.*(1.-ETA**2)*(1.+PSI)*(1.+3.*ETA)
      AN(11)=9./32.*(1.-PSI**2)*(1.+ETA)*(1.-3.*PSI)
      AN(12)=9./32.*(1.-ETA**2)*(1.-PSI)*(1.-3.*ETA)
  200 CONTINUE
      DO 210 II=1,12
      XN(I,J)=XN(I,J)+AN(II)*X(II)
      YN(I,J)=YN(I,J)+AN(II)*Y(II)
  210 CONTINUE
      GO TO 520
C
C     READ DATA FOR A TRIANGULAR BLOCK
C
  220 CONTINUE
      NCOORD=9
      I1=IDATA(2)
      J1=IDATA(3)
      I2=IDATA(4)
      J2=IDATA(5)
      I3=IDATA(6)
      J3=IDATA(7)
      G1=RDATA(8)
      G2=RDATA(9)
      WORD=CDATA(10)
      XO=RDATA(11)
      YO=RDATA(12)
      IF (G1.EQ.0.) G1=1.0
      IF (G2.EQ.0.) G2=1.0
      IDIFF=ABS(I1-I2)
      JDIFF=ABS(J1-J2)
      NI=MAX(IDIFF,JDIFF)+1
      IDIFF=ABS(I1-I3)
      JDIFF=ABS(J1-J3)
      NJ=MAX(IDIFF,JDIFF)+1
      IF (NI.NE.NJ) CALL ERROR ('MESH','INCORRECT I,J SPECIFICATION FOR 
     1TBLOCK','MESH PART',NP,' ',0,' ',' ',1)
      IMIN=MIN(I1,I2,I3)
      JMIN=MIN(J1,J2,I3)
      IMAX=MAX(I1,I2,I3)
      JMAX=MAX(J1,J2,J3)
      IF (I1.GT.NN.OR.I2.GT.NN.OR.I3.GT.NN) CALL ERROR ('MESH','INPUT I 
     1VALUE EXCEEDS SPECIFIED MESH IMAX','MESH PART',NP,' ',0,'MESH COMM
     2AND','TBLOCK',1)
      IF (J1.GT.NM.OR.J2.GT.NM.OR.J3.GT.NM) CALL ERROR ('MESH','INPUT J 
     1VALUE EXCEEDS SPECIFIED MESH JMAX','MESH PART',NP,' ',0,'MESH COMM
     2AND','TBLOCK',1)
      CALL RDFFLD (NIN,RDATA,IDATA,CDATA)
      DO 230 I=1,NCOORD
      X(I)=RDATA(I)
  230 CONTINUE
      CALL RDFFLD (NIN,RDATA,IDATA,CDATA)
      DO 240 I=1,NCOORD
      Y(I)=RDATA(I)
  240 CONTINUE
      IF (WORD.EQ.'POLAR') CALL TRNSFM (X,Y,XO,YO,NCOORD)
C
C     CALCULATE THE MESH POINT COORDINATES ALONG SIDES 1 AND 3 OF THE
C     TRIANGLE IN THE PSI-ETA PLANE, STORE IN THE C AND H VECTORS
C
      IF (NI.GE.MAXEDG) CALL ERROR ('MESH','TOO MANY NODES ON EDGE OF O 
     1NE TBLOCK','NODES SPECIFIED',NI,'NODES ALLOWED',MAXEDG,' ',' ',1)
      IF (NJ.GE.MAXEDG) CALL ERROR ('MESH','TOO MANY NODES ON EDGE OF O 
     1NE TBLOCK','NODES SPECIFIED',NJ,'NODES ALLOWED',MAXEDG,' ',' ',1)
      CALL SEG1 (0.,0.,1.,0.,G1,CC,HH,NI)
      CALL SEG1 (0.,0.,0.,1.,G2,CA,HA,NJ)
C
C     CALCULATE XN AND YN
C
      DO 320 I=1,NI
      PSI=CC(I)
      NJ=NI-(I-1)
      DO 320 J=1,NJ
      ETA=HA(J)
      BETA=1.0-PSI-ETA
      AN(1)=BETA
      AN(2)=PSI
      AN(3)=ETA
      IFLAG=1
      DO 250 KK=4,10
      KI=14-KK
      AN(KI)=0.
      IF (X(KI).NE.0.0) GO TO 260
  250 CONTINUE
      GO TO 270
  260 CONTINUE
      CALL FILLT (X,Y,IFLAG)
  270 CONTINUE
      IF (IFLAG-2) 300, 280, 290
  280 CONTINUE
      AN(1)=AN(1)*(2.*BETA-1.0)
      AN(2)=AN(2)*(2.*PSI-1.0)
      AN(3)=AN(3)*(2.*ETA-1.0)
      AN(4)=4.0*BETA*PSI
      AN(5)=4.0*PSI*ETA
      AN(6)=4.0*ETA*BETA
      GO TO 300
  290 CONTINUE
      AN(1)=AN(1)*.5*(3.0*BETA-1.0)*(3.0*BETA-2.0)
      AN(2)=AN(2)*.5*(3.0*PSI-1.0)*(3.0*PSI-2.0)
      AN(3)=AN(3)*.5*(3.0*ETA-1.0)*(3.0*ETA-2.0)
      AN(4)=9.0*BETA*PSI*(3.0*BETA-1.0)/2.0
      AN(5)=9.0*PSI*ETA*(3.0*PSI-1.0)/2.0
      AN(6)=9.0*ETA*BETA*(3.0*ETA-1.0)/2.0
      AN(7)=9.0*BETA*PSI*(3.0*PSI-1.0)/2.0
      AN(8)=9.0*PSI*ETA*(3.0*ETA-1.0)/2.0
      AN(9)=9.0*ETA*BETA*(3.0*BETA-1.0)/2.0
      AN(10)=27.0*BETA*PSI*ETA
  300 CONTINUE
      II=(I2*(I-1)+I3*(J-1)+I1*(NI-1)-I1*(I-1)-I1*(J-1))/(NI-1)
      JJ=(J2*(I-1)+J3*(J-1)+J1*(NI-1)-J1*(I-1)-J1*(J-1))/(NI-1)
      XN(II,JJ)=0.
      YN(II,JJ)=0.
      DO 310 KK=1,10
      XN(II,JJ)=XN(II,JJ)+AN(KK)*X(KK)
      YN(II,JJ)=YN(II,JJ)+AN(KK)*Y(KK)
  310 CONTINUE
  320 CONTINUE
      GO TO 520
C
C     COPY/REFLECT EXISTING DATA BLOCK
C
  330 CONTINUE
      NCOORD=2
      IREF=-1
      IF (CDATA(1).EQ.'REFLECT   ') IREF=1
      IMIN=IDATA(2)
      JMIN=IDATA(3)
      IMAX=IDATA(4)
      JMAX=IDATA(5)
      INEW1=IDATA(6)
      JNEW1=IDATA(7)
      INEW2=IDATA(8)
      JNEW2=IDATA(9)
      IF (INEW1.GT.NN.OR.INEW2.GT.NN) CALL ERROR ('MESH','INPUT I VALUE 
     1EXCEEDS SPECIFIED MESH IMAX','MESH PART',NP,' ',0,'MESH COMMAND', 
     2'COPY OR REFLECT',1)
      IF (JNEW1.GT.NM.OR.JNEW2.GT.NM) CALL ERROR ('MESH','INPUT J VALUE 
     1EXCEEDS SPECIFIED MESH JMAX','MESH PART',NP,' ',0,'MESH COMMAND', 
     2'COPY OR REFLECT',1)
      NI=IMAX-IMIN
      NINEW=INEW2-INEW1
      ISIGN=1
      IF (NINEW.LT.0) ISIGN=-1
      IF (ABS(NINEW).NE.NI) INEW2=INEW1+ISIGN*NI
      NI=INEW1-IMIN
      NJ=JMAX-JMIN
      NJNEW=JNEW2-JNEW1
      JSIGN=1
      IF (NJNEW.LT.0) JSIGN=-1
      IF (ABS(NJNEW).NE.NJ) JNEW2=JNEW1+JSIGN*NJ
      NJ=JNEW1-JMIN
      CALL RDFFLD (NIN,RDATA,IDATA,CDATA)
      X(1)=RDATA(1)
      Y(1)=RDATA(2)
      X(2)=RDATA(3)
      Y(2)=RDATA(4)
      DO 340 I=3,12
      X(I)=0.
      Y(I)=0.
  340 CONTINUE
      G1=1.0
      G2=1.0
      G3=1.0
      G4=1.0
      WORD=CDATA(5)
      XO=RDATA(6)
      YO=RDATA(7)
      IF (WORD.EQ.'POLAR') CALL TRNSFM (X,Y,XO,YO,NCOORD)
      IF (IREF.EQ.1) GO TO 370
C
C     TRANSLATE DATA BLOCK
C
      DELX=X(1)-XN(IMIN,JMIN)
      DELY=Y(1)-YN(IMIN,JMIN)
      DO 350 I=IMIN,IMAX
      DO 350 J=JMIN,JMAX
      XN(I+NI,J+NJ)=XN(I,J)+DELX
      YN(I+NI,J+NJ)=YN(I,J)+DELY
  350 CONTINUE
C
C     ROTATE DATA BLOCK
C
      DELX=XN(IMAX,JMAX)-XN(IMIN,JMIN)
      DELY=YN(IMAX,JMAX)-YN(IMIN,JMIN)
      IMIN=INEW1
      JMIN=JNEW1
      IMAX=INEW2
      JMAX=JNEW2
      V=SQRT(DELX**2+DELY**2)
      VX=DELX/V
      ANGO=ACOS(VX)
      DELX=X(2)-X(1)
      DELY=Y(2)-Y(1)
      V=SQRT(DELX**2+DELY**2)
      VX=DELX/V
      ANGN=ACOS(VX)
      ANG=ANGN-ANGO
      IF (ABS(ANG).LE.1.0E-2) GO TO 520
      XO=XN(INEW1,JNEW1)
      YO=YN(INEW1,JNEW1)
      DO 360 I=INEW1,INEW2
      DO 360 J=JNEW1,JNEW2
      DELX=XN(I,J)-XO
      DELY=YN(I,J)-YO
      XN(I,J)=XO+DELX*COS(ANG)-DELY*SIN(ANG)
      YN(I,J)=YO+DELX*SIN(ANG)+DELY*COS(ANG)
  360 CONTINUE
      GO TO 520
C
C     REFLECT DATA BLOCK
C
  370 CONTINUE
      DELX=X(2)-X(1)
      DELY=Y(2)-Y(1)
      IF (ABS(DELX).LT.1.0E-10) DELX=1.0E-10
      IF (ABS(DELY).LT.1.0E-10) DELY=1.0E-10
      SLOPE=DELY/DELX
      PSLOPE=-1.0/SLOPE
      V=PSLOPE-SLOPE
      DO 380 I=IMIN,IMAX
      DO 380 J=JMIN,JMAX
      XP=XN(I,J)
      YP=YN(I,J)
      XS=(Y(1)-YP)/V-(SLOPE*X(1)-PSLOPE*XP)/V
      YS=YP+PSLOPE*(XS-XP)
      DELX=XS-XP
      DELY=YS-YP
      IINC=ISIGN*(I-IMIN)+INEW1
      JINC=JSIGN*(J-JMIN)+JNEW1
      XN(IINC,JINC)=XS+DELX
      YN(IINC,JINC)=YS+DELY
  380 CONTINUE
      IMIN=MIN(INEW1,INEW2)
      JMIN=MIN(JNEW1,JNEW2)
      IMAX=MAX(INEW1,INEW2)
      JMAX=MAX(JNEW1,JNEW2)
      GO TO 520
C
C     USER SUBROUTINE TO GENERATE COORDINATE DATA
C
  390 CONTINUE
      IEXT=IEXT+1
      MAXI=NN
      MAXJ=NM
      CALL USRMSH (XN,YN,MAXI,MAXJ)
      GO TO 520
C
C     FILL IN COORDINATES OF UNDEFINED POINTS (VIA LAPLACE EQUATION)
C
  400 CONTINUE
      ILAP=1
      DO 410 I=1,12
      X(I)=0.
      Y(I)=0.
  410 CONTINUE
      G1=1.0
      G2=1.0
      G3=1.0
      G4=1.0
      DO 420 I=1,NN
      DO 420 J=1,NM
      IF (XN(I,J).NE.DUM) IBND(I,J)=1
  420 CONTINUE
      INC=IDATA(6)
      IF (INC.EQ.0) INC=1
      ALPHA=RDATA(7)
      IF (ALPHA.EQ.0.) ALPHA=.5
      OMEGA=RDATA(8)
      IF (OMEGA.EQ.0.) OMEGA=1.86
      ITER=IDATA(9)
      IF (ITER.EQ.0) ITER=100
      TOL=RDATA(10)
      IF (TOL.EQ.0.) TOL=.0001
      IMIN=IDATA(2)+INC
      JMIN=IDATA(3)+INC
      IMAX=IDATA(4)-INC
      JMAX=IDATA(5)-INC
      IF (IMAX.LE.0) IMAX=NN-INC
      IF (JMAX.LE.0) JMAX=NM-INC
      ALFA=0.
      ICONV=0
      DO 480 IT=1,ITER
      IF (IT.GT.30) ALFA=ALPHA
      RESID=0.
      DO 470 I=IMIN,IMAX,INC
      DO 470 J=JMIN,JMAX,INC
      IF (IBND(I,J).EQ.1) GO TO 470
      XI(1)=XN(I,J-INC)
      XI(2)=XN(I+INC,J-INC)
      XI(3)=XN(I+INC,J)
      XI(4)=XN(I+INC,J+INC)
      XI(5)=XN(I,J+INC)
      XI(6)=XN(I-INC,J+INC)
      XI(7)=XN(I-INC,J)
      XI(8)=XN(I-INC,J-INC)
      XI(9)=XI(1)
      YI(1)=YN(I,J-INC)
      YI(2)=YN(I+INC,J-INC)
      YI(3)=YN(I+INC,J)
      YI(4)=YN(I+INC,J+INC)
      YI(5)=YN(I,J+INC)
      YI(6)=YN(I-INC,J+INC)
      YI(7)=YN(I-INC,J)
      YI(8)=YN(I-INC,J-INC)
      YI(9)=YI(1)
      DO 430 L=1,4
      A(L)=0.5*((XI(2*L)-XN(I,J))*(YI(2*L+1)-YI(2*L-1))-                
     1(YI(2*L)-YN(I,J))*(XI(2*L+1)-XI(2*L-1)))
      DX(L)=XI(2*L+1)-XN(I,J)
      DY(L)=YI(2*L+1)-YN(I,J)
  430 CONTINUE
      A(5)=A(1)
      IF (IT.GT.1) GO TO 440
      A(1)=1.0
      A(2)=1.0
      A(3)=1.0
      A(4)=1.0
      A(5)=1.0
  440 CONTINUE
      DXA=0.
      DYA=0.
      DO 450 L=1,5
      A(L)=ABS(A(L))
  450 CONTINUE
      DO 460 L=1,4
      DXA=ALFA*(A(L+1)-A(L))/(A(L+1)+A(L))*(-DY(L))+                    
     1(1.0-ALFA)*DX(L)+DXA
      DYA=ALFA*(A(L+1)-A(L))/(A(L+1)+A(L))*DX(L)+                       
     1(1.0-ALFA)*DY(L)+DYA
  460 CONTINUE
      RESID=RESID+ABS(DXA)+ABS(DYA)
      XN(I,J)=XN(I,J)+OMEGA*DXA/4.
      YN(I,J)=YN(I,J)+OMEGA*DYA/4.
  470 CONTINUE
      IF (RESID.LT.TOL) THEN
      ICONV=1
      GO TO 490
      END IF
  480 CONTINUE
      WRITE (NOUT, 640) IT,RESID
  490 CONTINUE
      IF (INC.EQ.1) GO TO 520
      MINI=IMIN-INC
      MAXI=IMAX+INC
      MINJ=JMIN-INC
      MAXJ=JMAX+INC
      DO 510 I=MINI,MAXI,INC
      DO 510 J=MINJ,MAXJ,INC
      IINC=I+INC
      JINC=J+INC
      DELX=(XN(IINC,JINC)-XN(I,J))/INC
      DELY=(YN(IINC,JINC)-YN(I,J))/INC
      INCC=INC-1
      DO 500 IL=1,INCC
      IINC1=I+IL
      JINC1=J+IL
      XN(IINC1,JINC1)=XN(I,J)+DELX*IL
      YN(IINC1,JINC1)=YN(I,J)+DELY*IL
  500 CONTINUE
  510 CONTINUE
C
C     PRINT COORDINATE DATA
C
  520 CONTINUE
      IF (NP.EQ.1.AND.IPRINT.GT.2) THEN
      WRITE (NOUT, 580)
      WRITE (NOUT, 650)
      END IF
      IF (IPRINT.LT.3) GO TO 550
      IF (IEXT.GT.0) THEN
      WRITE (NOUT, 590) NP
      WRITE (NOUT, 650)
      GO TO 20
      END IF
      WRITE (NOUT, 600) NP,IMIN,IMAX,JMIN,JMAX,G1,G2,G3,G4
      WRITE (NOUT, 620)
      DO 540 I=1,12
      IF (I.LT.5) GO TO 530
      IF (X(I).EQ.0.0.AND.Y(I).EQ.0.0) GO TO 540
  530 CONTINUE
      WRITE (NOUT, 610) I,X(I),I,Y(I)
  540 CONTINUE
      IF (IREF.EQ.-1) WRITE (NOUT, 660)
      IF (IREF.EQ.1) WRITE (NOUT, 670)
      IF (ILAP.EQ.1) WRITE (NOUT, 680) INC,ALPHA,OMEGA,ITER,TOL
      IF (ILAP.EQ.1.AND.ICONV.EQ.1) WRITE (NOUT, 630) IT
  550 CONTINUE
      IF (IPRINT.LT.4) GO TO 20
      WRITE (NOUT, 560) NP
      WRITE (NOUT, 570) ((I,J,XN(I,J),YN(I,J),I=IMIN,IMAX),J=JMIN,JMAX)
      WRITE (NOUT, 650)
      GO TO 20
C
  560 FORMAT (////,10X,'MESH POINT COORDINATES FOR PART NO.', I2/16X,
     1'I',14X,'J',10X,'X',14X,'Y'//)
  570 FORMAT (12X,2I5,2E15.4)
  580 FORMAT (///,3X,'MESH POINT DATA FOR PROBLEM GEOMETRY' ,//)
  590 FORMAT (///,10X,'INPUT DATA FOR PART NO.' ,I3,//,10X,'DATA GENERAT
     1ED BY EXTERNAL SUBROUTINE' ,/)
  600 FORMAT (///10X,'INPUT DATA FOR PART NO.',I3,//16X,'IMIN=',I3,5X,
     1'IMAX=',I3,/16X,'JMIN=',I3,5X,'JMAX=',I3,//10X,'GRADIENTS FOR SIZE
     2VARIATION' ,//16X,'G1=',F6.2,5X,'G2=',F6.2,/16X,'G3=',F6.2,5X,
     3'G4=',F6.2/)
  610 FORMAT (16X,'X(',I2,')=',F10.4,5X,'Y(',I2,')=',F10.4)
  620 FORMAT (10X,'DATA POINTS' /)
  630 FORMAT (/,10X,I4,' ITERATIONS REQUIRED FOR CONVERGENCE OF' ,/,10X,
     1'LAPLACE MESH GENERATOR' )
  640 FORMAT (10X,'***WARNING-' ,/,10X,'NO CONVERGENCE OF LAPLACE MESH G
     1ENERATOR' ,/,10X,'IN ',I3,' ITERATIONS -- RESIDUAL= '  ,E15.7)
  650 FORMAT (///,2X,'**************************************************
     1******************************************************************
     2**************'   ,//)
  660 FORMAT (//,10X,'MESH POINTS IN THIS PART WERE' ,/,10X,'GENERATED V
     1IA A COPY COMMAND' )
  670 FORMAT (//,10X,'MESH POINTS IN THIS PART WERE' ,/,10X,'GENERATED V
     1IA A REFLECTION COMMAND' )
  680 FORMAT (//,10X,'MESH POINTS IN THIS PART WERE GENERATED' ,/,10X,'V
     1IA AN ITERATIVE SOLUTION OF LAPLACES EQUATION'   ,//,10X,'PARAMETE
     2RS-'          ,/,16X,'IJ INCREMENT - ' ,I5,/,16X,'ALPHA - ',E15.7,
     3/,16X,'OMEGA - ',E15.7,/,16X,'MAX ITERATION - ' ,I5,/,16X,'CONVERG
     4ENCE TOL - '         ,E15.7)
      END
