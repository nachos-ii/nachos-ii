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
      SUBROUTINE ASSMT (KNUM,KIND,MAT,NN,NCOR,IBCPT2,IBCPT3,CP,CX,CY,   
     1KT,KTQ,KTR,FT,FTQ,FTR,FTS,X,Y,UP,VP,PP,TP,V1P,V2P,CXYT,CXYTP)
C
C     ******************************************************************
C
C     SUBROUTINE TO ASSEMBLE ELEMENT MATRICES FOR ENERGY EQUATION
C
C     ******************************************************************
C
      REAL KT,KTQ,KTR
C
      COMMON /RSTART/ IRSTRT,NSTEPS
      COMMON /MATDAT/ NMAT,NPROP,PROP(15,10),NXMAT,NXPROP,XPROP(15,10)
      COMMON /PRBDAT/ IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,    
     1                IPNLTY,PNLTY
      COMMON /SLNDAT/ IALGOR,IBLDUV,IBLDT,IJACOB,IAUTO,NCYCLM,NPRNT,    
     1                KSTEP,KPRNT,TIME,DELTN,DELTNM
C
      DIMENSION KT(9,9,4), CX(9,9,9), CY(9,9,9), FT(9), FTS(9,4)
      DIMENSION CP(9,9)
      DIMENSION KTQ(9,9), FTQ(9), KTR(9,3,9), FTR(9,3)
      DIMENSION CXYT(9,9), CXYTP(9,18)
      DIMENSION UP(9), VP(9), PP(9), TP(9), V1P(9), V2P(9), X(9), Y(9)
      DIMENSION AC(4), S(4)
C
      PARAMETER (BIG=1.0E40)
C
C     ******************************************************************
C
C     EVALUATE THERMAL CONDUCTIVITY
C
      COND=PROP(4,MAT)
      IF (PROP(15,MAT).EQ.4.0) COND=PROP(7,MAT)
      IF (PROP(8,MAT).EQ.1.0) GO TO 30
      IF (KSTEP.EQ.1.AND.ITMDEP.EQ.0.AND.IRSTRT.EQ.0) GO TO 10
      CALL DEFORM (KIND,X,Y,UP,VP,S)
      CALL USRCON (AC,TP,S,V1P,V2P,NCOR,MAT)
      KLEND=NCOR
      GO TO 50
   10 CONTINUE
      DO 20 I=1,NCOR
      AC(I)=COND
   20 CONTINUE
      KLEND=NCOR
      GO TO 50
   30 CONTINUE
      DO 40 I=1,NCOR
      AC(I)=0.
   40 CONTINUE
      AC(1)=COND
      KLEND=1
   50 CONTINUE
C
C     CONSTRUCT DIFFUSION MATRIX
C
      DO 70 I=1,NN
      DO 70 J=1,NN
      FD=0.
      DO 60 KL=1,KLEND
      FD=FD+AC(KL)*KT(I,J,KL)
   60 CONTINUE
      KT(I,J,1)=FD+KTQ(I,J)
   70 CONTINUE
C
C     MODIFY THE DIFFUSION MATRIX FOR RADIATION, IF REQUIRED
C
      MM=IBCPT3
      IQRAD=MM/10
      IF (IQRAD.EQ.0) GO TO 100
      ISIDE=MOD(MM,10)/1
      CALL RADMOD (KIND,IQRAD,ISIDE,KTR,FTR,TP,NN,TIME,X,Y,KNUM)
      DO 80 I=1,NN
      DO 80 J=1,NN
      KT(I,J,1)=KT(I,J,1)+KTR(I,1,J)
   80 CONTINUE
      DO 90 J=1,NN
      FT(J)=FT(J)+FTR(J,1)
   90 CONTINUE
  100 CONTINUE
C
C     CHECK BOUNDARY CONDITIONS
C
      DO 110 I=1,NN
      FD=ABS(KT(I,I,1))
      IF (FD.GT.1.0E15) KT(I,I,1)=BIG
  110 CONTINUE
C
C     CONSTRUCT ADVECTION MATRIX
C
      RHOCP=PROP(1,MAT)*PROP(3,MAT)
      IF (PROP(15,MAT).EQ.4.0) THEN
      MATF=INT(PROP(1,MAT))
      RHOCP=PROP(1,MATF)*PROP(3,MATF)
      END IF
      DO 130 I=1,NN
      DO 130 K=1,NN
      FDX=0.
      FDY=0.
      DO 120 J=1,NN
      FDX=FDX+CX(I,J,K)*UP(J)
      FDY=FDY+CY(I,J,K)*VP(J)
  120 CONTINUE
      CXYT(I,K)=RHOCP*(FDX+FDY)
  130 CONTINUE
C
C     CONSTRUCT JACOBIAN TERMS, IF REQUIRED
C
      DO 140 I=1,9
      DO 140 J=1,18
      CXYTP(I,J)=0.
  140 CONTINUE
      IF (IJACOB.EQ.0) GO TO 170
      DO 160 I=1,NN
      DO 160 K=1,NN
      FDX=0.
      FDY=0.
      DO 150 J=1,NN
      FDX=FDX+CX(I,K,J)*TP(J)
      FDY=FDY+CY(I,K,J)*TP(J)
  150 CONTINUE
      CXYTP(I,K)=RHOCP*FDX
      CXYTP(I,K+NN)=RHOCP*FDY
  160 CONTINUE
  170 CONTINUE
C
C     CONSTRUCT CAPACITANCE MATRIX, IF REQUIRED
C
      IF (ITMDEP.EQ.0) GO TO 190
      IF (PROP(15,MAT).EQ.4.0) RHOCP=PROP(6,MAT)
      DO 180 I=1,NN
      DO 180 J=1,NN
      CP(I,J)=RHOCP*CP(I,J)
  180 CONTINUE
  190 CONTINUE
C
C     UPDATE VARIABLE RHS, IF REQUIRED
C
      MM=IBCPT2
      ITIME=MM/100000000
      IF (ITIME.EQ.0) GO TO 200
      ICURVE=MOD(MM,100000000)/10000000
      ISIDE=MOD(MM,10000000)/1000000
      CALL BCTIME (KNUM,4,ICURVE,ISIDE,ITIME,TIME,KSTEP,KIND,X,Y,       
     1UP,VP,PP,TP,V1P,V2P,FT)
  200 CONTINUE
C
C     CONSTRUCT SOURCE TERMS
C
      DO 210 I=1,NCOR
      AC(I)=PROP(9,MAT)
  210 CONTINUE
      IF (PROP(9,MAT).EQ.-1.0) CALL USRVHS (AC,TP,V1P,V2P,X,Y,NCOR,     
     1MAT,KNUM,TIME,KSTEP)
      DO 230 I=1,NN
      FD=0.
      DO 220 J=1,NCOR
      FD=FD+FTS(I,J)*AC(J)
  220 CONTINUE
      FT(I)=FT(I)+FD+FTQ(I)
  230 CONTINUE
      IF (PROP(10,MAT).EQ.-1.0) THEN
      CALL DEFORM (KIND,X,Y,UP,VP,S)
      DO 240 I=1,NCOR
      AC(I)=PROP(2,MAT)
  240 CONTINUE
      IF (PROP(8,MAT).EQ.-1.0) CALL USRVIS (AC,TP,S,V1P,V2P,NCOR,MAT,   
     1KSTEP)
      AVG=0.
      DO 250 I=1,NCOR
      AVG=AVG+AC(I)*S(I)
  250 CONTINUE
      AVG=AVG/NCOR
      DO 270 I=1,NN
      FD=0.
      DO 260 J=1,NCOR
      FD=FD+FTS(I,J)*AVG
  260 CONTINUE
      FT(I)=FT(I)+FD
  270 CONTINUE
      END IF
C
      RETURN
      END
