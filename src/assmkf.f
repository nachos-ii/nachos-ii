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
      SUBROUTINE ASSMKF (KNUM,KIND,NN,NCOR,MASS,CXY,CXYP,KUV,FUV,FBDY,  
     1CP,CXYT,CXYTP,KT,FT,CPV1,CPV2,CXYV1,CXYV2,CXYV1P,CXYV2P,KV1,KV2,  
     2FV1,FV2,U,V,P,T,V1,V2,UP,VP,PP,TP,V1P,V2P,UDOT,VDOT,PDOT,TDOT,    
     3V1DOT,V2DOT,GX,GY,ELSTIF,ELFV)
C
C     ******************************************************************
C
C     SUBROUTINE TO ASSEMBLE COMPONENT ELEMENT MATRICES FOR VARIOUS
C     STEADY STATE AND TRANSIENT SOLUTION ALGORITHMS
C
C     ******************************************************************
C
      REAL MASS,KUV,KT,KV1,KV2
C
      COMMON /PRBDAT/ IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,    
     1                IPNLTY,PNLTY
      COMMON /SLNDAT/ IALGOR,IBLDUV,IBLDT,IJACOB,IAUTO,NCYCLM,          
     1                NPRNT,KSTEP,KPRNT,TIME,DELTN,DELTNM
C
      DIMENSION KUV(22,22,4), CXY(9,9), CXYP(18,18), MASS(9,9), FUV(22)
      DIMENSION FBDY(9,9,4)
      DIMENSION KT(9,9,4), CXYT(9,9), CXYTP(9,18), CP(9,9), FT(9)
      DIMENSION KV1(9,9,4), KV2(9,9,4), CXYV1(9,9), CXYV2(9,9)
      DIMENSION CXYV1P(9,18), CXYV2P(9,18), CPV1(9,9), CPV2(9,9)
      DIMENSION FV1(9), FV2(9)
      DIMENSION ELSTIF(50,50), ELFV(50)
      DIMENSION U(9), V(9), P(9), T(9), V1(9), V2(9)
      DIMENSION UP(9), VP(9), PP(9), TP(9), V1P(9), V2P(9)
      DIMENSION UDOT(9), VDOT(9), PDOT(9), TDOT(9), V1DOT(9), V2DOT(9)
C
C     ******************************************************************
C
      DO 10 I=1,50
      ELFV(I)=0.
      DO 10 J=1,50
      ELSTIF(I,J)=0.
   10 CONTINUE
      IF (ITMDEP.EQ.1) GO TO 180
C
C     **********   STEADY-STATE SOLUTION PROCEDURES   **********
C
C     MOMENTUM AND CONTINUITY EQUATIONS
C
      LOOP=2*NN+NCOR
      DO 20 I=1,LOOP
      ELFV(I)=FUV(I)
      DO 20 J=1,LOOP
      ELSTIF(I,J)=KUV(I,J,1)
   20 CONTINUE
      IF (IALGOR.EQ.2.OR.IJACOB.EQ.1) THEN
      DO 40 I=1,NN
      FDX=0.
      FDY=0.
      FDXY=0.
      FDYX=0.
      DO 30 J=1,NN
      FDX=FDX+CXYP(I,J)*UP(J)
      FDXY=FDXY+CXYP(I,J+NN)*VP(J)
      FDYX=FDYX+CXYP(I+NN,J)*UP(J)
      FDY=FDY+CXYP(I+NN,J+NN)*VP(J)
   30 CONTINUE
      ELFV(I)=ELFV(I)+FDX+FDXY
      ELFV(I+NN)=ELFV(I+NN)+FDYX+FDY
   40 CONTINUE
      END IF
      DO 50 I=1,NN
      DO 50 J=1,NN
      ELSTIF(I,J)=ELSTIF(I,J)+CXY(I,J)+CXYP(I,J)
      ELSTIF(I,J+NN)=ELSTIF(I,J+NN)+CXYP(I,J+NN)
      ELSTIF(I+NN,J)=ELSTIF(I+NN,J)+CXYP(I+NN,J)
      ELSTIF(I+NN,J+NN)=ELSTIF(I+NN,J+NN)+CXY(I,J)+CXYP(I+NN,J+NN)
   50 CONTINUE
C
C     COUPLED TRANSPORT EQUATIONS
C
      IF (IBLDT.EQ.0) RETURN
      IADD=2*NN+NCOR
      DO 60 I=1,NN
      ELFV(I+IADD)=FT(I)
      DO 60 J=1,NN
      ELSTIF(I,J+IADD)=-FBDY(I,J,1)*GX
      ELSTIF(I+NN,J+IADD)=-FBDY(I,J,1)*GY
      ELSTIF(I+IADD,J+IADD)=CXYT(I,J)+KT(I,J,1)
   60 CONTINUE
      IF (IALGOR.EQ.2.OR.IJACOB.EQ.1) THEN
      DO 80 I=1,NN
      FDX=0.
      FDY=0.
      DO 70 J=1,NN
      FDX=FDX+CXYTP(I,J)*UP(J)
      FDY=FDY+CXYTP(I,J+NN)*VP(J)
   70 CONTINUE
      ELFV(I+IADD)=ELFV(I+IADD)+FDX+FDY
   80 CONTINUE
      DO 90 I=1,NN
      DO 90 J=1,NN
      ELSTIF(I+IADD,J)=ELSTIF(I+IADD,J)+CXYTP(I,J)
      ELSTIF(I+IADD,J+NN)=ELSTIF(I+IADD,J+NN)+CXYTP(I,J+NN)
   90 CONTINUE
      END IF
C
C     COUPLED AUXILIARY TRANSPORT EQUATIONS
C
      IF (IVAR1.EQ.0) RETURN
      IADD=3*NN+NCOR
      DO 100 I=1,NN
      ELFV(I+IADD)=FV1(I)
      DO 100 J=1,NN
      ELSTIF(I,J+IADD)=-FBDY(I,J,2)*GX
      ELSTIF(I+NN,J+IADD)=-FBDY(I,J,2)*GY
      ELSTIF(I+IADD,J+IADD)=CXYV1(I,J)+KV1(I,J,1)
  100 CONTINUE
      IF (IALGOR.EQ.2.OR.IJACOB.EQ.1) THEN
      DO 120 I=1,NN
      FDX=0.
      FDY=0.
      DO 110 J=1,NN
      FDX=FDX+CXYV1P(I,J)*UP(J)
      FDY=FDY+CXYV1P(I,J+NN)*VP(J)
  110 CONTINUE
      ELFV(I+IADD)=ELFV(I+IADD)+FDX+FDY
  120 CONTINUE
      DO 130 I=1,NN
      DO 130 J=1,NN
      ELSTIF(I+IADD,J)=ELSTIF(I+IADD,J)+CXYV1P(I,J)
      ELSTIF(I+IADD,J+NN)=ELSTIF(I+IADD,J+NN)+CXYV1P(I,J+NN)
  130 CONTINUE
      END IF
C
      IF (IVAR2.EQ.0) RETURN
      IADD=4*NN+NCOR
      DO 140 I=1,NN
      ELFV(I+IADD)=FV2(I)
      DO 140 J=1,NN
      ELSTIF(I,J+IADD)=-FBDY(I,J,3)*GX
      ELSTIF(I+NN,J+IADD)=-FBDY(I,J,3)*GY
      ELSTIF(I+IADD,J+IADD)=CXYV2(I,J)+KV2(I,J,1)
  140 CONTINUE
      IF (IALGOR.EQ.2.OR.IJACOB.EQ.1) THEN
      DO 160 I=1,NN
      FDX=0.
      FDY=0.
      DO 150 J=1,NN
      FDX=FDX+CXYV2P(I,J)*UP(J)
      FDY=FDY+CXYV2P(I,J+NN)*VP(J)
  150 CONTINUE
      ELFV(I+IADD)=ELFV(I+IADD)+FDX+FDY
  160 CONTINUE
      DO 170 I=1,NN
      DO 170 J=1,NN
      ELSTIF(I+IADD,J)=ELSTIF(I+IADD,J)+CXYV2P(I,J)
      ELSTIF(I+IADD,J+NN)=ELSTIF(I+IADD,J+NN)+CXYV2P(I,J+NN)
  170 CONTINUE
      END IF
      RETURN
C
C     **********   TRANSIENT SOLUTION PROCEDURES   **********
C
C     MOMENTUM AND CONTINUITY EQUATIONS
C
  180 CONTINUE
      FACT=1./DELTN
      IF (IALGOR.EQ.2.AND.KSTEP.GE.2) FACT=2./DELTN
C
      LOOP=2*NN+NCOR
      DO 190 I=1,LOOP
      ELFV(I)=FUV(I)
      DO 190 J=1,LOOP
      ELSTIF(I,J)=KUV(I,J,1)
  190 CONTINUE
      IF (IJACOB.EQ.1) THEN
      DO 210 I=1,NN
      FDX=0.
      FDY=0.
      FDXY=0.
      FDYX=0.
      DO 200 J=1,NN
      FDX=FDX+CXYP(I,J)*UP(J)
      FDXY=FDXY+CXYP(I,J+NN)*VP(J)
      FDYX=FDYX+CXYP(I+NN,J)*UP(J)
      FDY=FDY+CXYP(I+NN,J+NN)*VP(J)
  200 CONTINUE
      ELFV(I)=ELFV(I)+FDX+FDXY
      ELFV(I+NN)=ELFV(I+NN)+FDYX+FDY
  210 CONTINUE
      END IF
      DO 230 I=1,NN
      FDX=0.
      FDY=0.
      FDXY=0.
      FDYX=0.
      DO 220 J=1,NN
      FDX=FDX+MASS(I,J)*U(J)
      FDY=FDY+MASS(I,J)*V(J)
      FDXY=FDXY+MASS(I,J)*UDOT(J)
      FDYX=FDYX+MASS(I,J)*VDOT(J)
  220 CONTINUE
      ELFV(I)=ELFV(I)+FDX*FACT+FDXY*(IALGOR-1)
      ELFV(I+NN)=ELFV(I+NN)+FDY*FACT+FDYX*(IALGOR-1)
  230 CONTINUE
      DO 240 I=1,NN
      DO 240 J=1,NN
      ELSTIF(I,J)=ELSTIF(I,J)+CXY(I,J)+CXYP(I,J)+FACT*MASS(I,J)
      ELSTIF(I,J+NN)=ELSTIF(I,J+NN)+CXYP(I,J+NN)
      ELSTIF(I+NN,J)=ELSTIF(I+NN,J)+CXYP(I+NN,J)
      ELSTIF(I+NN,J+NN)=ELSTIF(I+NN,J+NN)+CXY(I,J)+CXYP(I+NN,J+NN)+     
     1FACT*MASS(I,J)
  240 CONTINUE
C
C     COUPLED TRANSPORT EQUATIONS
C
      IF (IBLDT.EQ.0) RETURN
      IADD=2*NN+NCOR
      DO 260 I=1,NN
      FDX=0.
      FDY=0.
      DO 250 J=1,NN
      FDX=FDX+CP(I,J)*T(J)
      FDY=FDY+CP(I,J)*TDOT(J)
  250 CONTINUE
      ELFV(I+IADD)=FT(I)+FDX*FACT+FDY*(IALGOR-1)
  260 CONTINUE
      DO 270 I=1,NN
      DO 270 J=1,NN
      ELSTIF(I,J+IADD)=-FBDY(I,J,1)*GX
      ELSTIF(I+NN,J+IADD)=-FBDY(I,J,1)*GY
      ELSTIF(I+IADD,J+IADD)=CXYT(I,J)+KT(I,J,1)+CP(I,J)*FACT
  270 CONTINUE
      IF (IJACOB.EQ.1) THEN
      DO 290 I=1,NN
      FDX=0.
      FDY=0.
      DO 280 J=1,NN
      FDX=FDX+CXYTP(I,J)*UP(J)
      FDY=FDY+CXYTP(I,J+NN)*VP(J)
  280 CONTINUE
      ELFV(I+IADD)=ELFV(I+IADD)+FDX+FDY
  290 CONTINUE
      DO 300 I=1,NN
      DO 300 J=1,NN
      ELSTIF(I+IADD,J)=ELSTIF(I+IADD,J)+CXYTP(I,J)
      ELSTIF(I+IADD,J+NN)=ELSTIF(I+IADD,J+NN)+CXYTP(I,J+NN)
  300 CONTINUE
      END IF
C
C     COUPLED AUXILIARY TRANSPORT EQUATIONS
C
      IF (IVAR1.EQ.0) RETURN
      IADD=3*NN+NCOR
      DO 320 I=1,NN
      FDX=0.
      FDY=0.
      DO 310 J=1,NN
      FDX=FDX+CPV1(I,J)*V1(J)
      FDY=FDY+CPV1(I,J)*V1DOT(J)
  310 CONTINUE
      ELFV(I+IADD)=FV1(I)+FDX*FACT+FDY*(IALGOR-1)
  320 CONTINUE
      DO 330 I=1,NN
      DO 330 J=1,NN
      ELSTIF(I,J+IADD)=-FBDY(I,J,2)*GX
      ELSTIF(I+NN,J+IADD)=-FBDY(I,J,2)*GY
      ELSTIF(I+IADD,J+IADD)=CXYV1(I,J)+KV1(I,J,1)+CPV1(I,J)*FACT
  330 CONTINUE
      IF (IJACOB.EQ.1) THEN
      DO 350 I=1,NN
      FDX=0.
      FDY=0.
      DO 340 J=1,NN
      FDX=FDX+CXYV1P(I,J)*UP(J)
      FDY=FDY+CXYV1P(I,J+NN)*VP(J)
  340 CONTINUE
      ELFV(I+IADD)=ELFV(I+IADD)+FDX+FDY
  350 CONTINUE
      DO 360 I=1,NN
      DO 360 J=1,NN
      ELSTIF(I+IADD,J)=ELSTIF(I+IADD,J)+CXYV1P(I,J)
      ELSTIF(I+IADD,J+NN)=ELSTIF(I+IADD,J+NN)+CXYV1P(I,J+NN)
  360 CONTINUE
      END IF
C
      IF (IVAR2.EQ.0) RETURN
      IADD=4*NN+NCOR
      DO 380 I=1,NN
      FDX=0.
      FDY=0.
      DO 370 J=1,NN
      FDX=FDX+CPV2(I,J)*V2(J)
      FDY=FDY+CPV2(I,J)*V2DOT(J)
  370 CONTINUE
      ELFV(I+IADD)=FV2(I)+FDX*FACT+FDY*(IALGOR-1)
  380 CONTINUE
      DO 390 I=1,NN
      DO 390 J=1,NN
      ELSTIF(I,J+IADD)=-FBDY(I,J,3)*GX
      ELSTIF(I+NN,J+IADD)=-FBDY(I,J,3)*GY
      ELSTIF(I+IADD,J+IADD)=CXYV2(I,J)+KV2(I,J,1)+CPV2(I,J)*FACT
  390 CONTINUE
      IF (IJACOB.EQ.1) THEN
      DO 410 I=1,NN
      FDX=0.
      FDY=0.
      DO 400 J=1,NN
      FDX=FDX+CXYV2P(I,J)*UP(J)
      FDY=FDY+CXYV2P(I,J+NN)*VP(J)
  400 CONTINUE
      ELFV(I+IADD)=ELFV(I+IADD)+FDX+FDY
  410 CONTINUE
      DO 420 I=1,NN
      DO 420 J=1,NN
      ELSTIF(I+IADD,J)=ELSTIF(I+IADD,J)+CXYV2P(I,J)
      ELSTIF(I+IADD,J+NN)=ELSTIF(I+IADD,J+NN)+CXYV2P(I,J+NN)
  420 CONTINUE
      END IF
C
      RETURN
      END
