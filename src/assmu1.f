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
      SUBROUTINE ASSMU1 (KNUM,KIND,MAT,NN,NCOR,IBCPT1,MASS,CX,CY,KUV,   
     1FUV,KP,FBDY,X,Y,UP,VP,PP,TP,V1P,V2P,CXY,CXYP,GX,GY)
C
C     ******************************************************************
C
C     SUBROUTINE TO ASSEMBLE ELEMENT MATRICES FOR MOMENTUM AND
C     CONTINUITY EQUATIONS (NAVIER-STOKES EQUATIONS)
C
C     ******************************************************************
C
      REAL MASS,KUV,KP
C
      COMMON /RSTART/ IRSTRT,NSTEPS
      COMMON /MATDAT/ NMAT,NPROP,PROP(15,10),NXMAT,NXPROP,XPROP(15,10)
      COMMON /PRBDAT/ IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,    
     1                IPNLTY,PNLTY
      COMMON /SLNDAT/ IALGOR,IBLDUV,IBLDT,IJACOB,IAUTO,NCYCLM,NPRNT,    
     1                KSTEP,KPRNT,TIME,DELTN,DELTNM
C
      DIMENSION KUV(22,22,4), CX(9,9,9), CY(9,9,9), FUV(22)
      DIMENSION FBDY(9,9,4), MASS(9,9), KP(18,18)
      DIMENSION CXY(9,9), CXYP(18,18)
      DIMENSION UP(9), VP(9), PP(9), TP(9), V1P(9), V2P(9), X(9), Y(9)
      DIMENSION AC(4), BC(4), CC(4), S(4)
C
      PARAMETER (BIG=1.0E40)
C
C     ******************************************************************
C
C     EVALUATE FLUID VISCOSITY
C
      IF (PROP(8,MAT).EQ.1.0) GO TO 30
      IF (KSTEP.EQ.1.AND.ITMDEP.EQ.0.AND.IRSTRT.EQ.0) GO TO 10
      IF (PROP(15,MAT).EQ.3.0) GO TO 10
      IF (PROP(15,MAT).EQ.2.0) CALL DEFORM (KIND,X,Y,UP,VP,S)
      CALL USRVIS (AC,TP,S,V1P,V2P,NCOR,MAT,KSTEP)
      KLEND=NCOR
      GO TO 50
   10 CONTINUE
      DO 20 I=1,NCOR
      AC(I)=PROP(2,MAT)
   20 CONTINUE
      KLEND=NCOR
      GO TO 50
   30 CONTINUE
      DO 40 I=1,NCOR
      AC(I)=0.
   40 CONTINUE
      AC(1)=PROP(2,MAT)
      KLEND=1
   50 CONTINUE
C
C     CONSTRUCT DIFFUSION MATRIX
C
      NROW=2*NN
      DO 70 I=1,NROW
      DO 70 J=1,NROW
      FD=0.
      DO 60 KL=1,KLEND
      FD=FD+AC(KL)*KUV(I,J,KL)
   60 CONTINUE
      KUV(I,J,1)=FD+(1.0/PNLTY)*KP(I,J)*IPNLTY
   70 CONTINUE
C
C     CHECK BOUNDARY CONDITIONS
C
      NROW=2*NN+NCOR
      DO 80 I=1,NROW
      FD=ABS(KUV(I,I,1))
      IF (FD.GT.1.0E15) KUV(I,I,1)=BIG
   80 CONTINUE
C
C     CONSTRUCT ADVECTION MATRIX
C
      RHO=PROP(1,MAT)
      DO 100 I=1,NN
      DO 100 K=1,NN
      FDX=0.
      FDY=0.
      DO 90 J=1,NN
      FDX=FDX+CX(I,J,K)*UP(J)
      FDY=FDY+CY(I,J,K)*VP(J)
   90 CONTINUE
      CXY(I,K)=RHO*(FDX+FDY)
  100 CONTINUE
C
C     CONSTRUCT JACOBIAN TERMS, IF REQUIRED
C
      DO 110 I=1,18
      DO 110 K=1,18
      CXYP(I,K)=0.
  110 CONTINUE
      IF (IJACOB.EQ.0) GO TO 140
      DO 130 I=1,NN
      DO 130 K=1,NN
      FDX=0.
      FDY=0.
      FDXY=0.
      FDYX=0.
      DO 120 J=1,NN
      FDX=FDX+CX(I,K,J)*UP(J)
      FDY=FDY+CY(I,K,J)*VP(J)
      FDXY=FDXY+CX(I,K,J)*VP(J)
      FDYX=FDYX+CY(I,K,J)*UP(J)
  120 CONTINUE
      CXYP(I,K)=RHO*FDX
      CXYP(I+NN,K+NN)=RHO*FDY
      CXYP(I,K+NN)=RHO*FDYX
      CXYP(I+NN,K)=RHO*FDXY
  130 CONTINUE
  140 CONTINUE
C
C     CONSTRUCT MASS MATRIX, IF REQUIRED
C
      IF (ITMDEP.EQ.0) GO TO 160
      DO 150 I=1,NN
      DO 150 J=1,NN
      MASS(I,J)=RHO*MASS(I,J)
  150 CONTINUE
  160 CONTINUE
C
C     UPDATE VARIABLE RHS, IF REQUIRED
C
      MM=IBCPT1
      ITIME=MM/100000000
      IF (ITIME.EQ.0) GO TO 170
      ICURVE=MOD(MM,100000000)/10000000
      ISIDE=MOD(MM,10000000)/1000000
      CALL BCTIME (KNUM,1,ICURVE,ISIDE,ITIME,TIME,KSTEP,KIND,X,Y,       
     1UP,VP,PP,TP,V1P,V2P,FUV)
  170 CONTINUE
      ITIME=MOD(MM,1000000)/100000
      IF (ITIME.EQ.0) GO TO 180
      ICURVE=MOD(MM,100000)/10000
      ISIDE=MOD(MM,10000)/1000
      CALL BCTIME (KNUM,2,ICURVE,ISIDE,ITIME,TIME,KSTEP,KIND,X,Y,       
     1UP,VP,PP,TP,V1P,V2P,FUV)
  180 CONTINUE
      ITIME=MOD(MM,1000)/100
      IF (ITIME.EQ.0) GO TO 190
      ICURVE=MOD(MM,100)/10
      ISIDE=MOD(MM,10)/1
      CALL BCTIME (KNUM,3,ICURVE,ISIDE,ITIME,TIME,KSTEP,KIND,X,Y,       
     1UP,VP,PP,TP,V1P,V2P,FUV)
  190 CONTINUE
C
C     CONSTRUCT BUOYANCY TERMS, IF REQUIRED
C
      IF (IFREE.EQ.0) GO TO 280
      TREF=PROP(11,MAT)
      V1REF=XPROP(5,MAT)
      V2REF=XPROP(11,MAT)
      GX=PROP(6,MAT)
      GY=PROP(7,MAT)
      IF (GX.EQ.-1.0) CALL USRGRV (GX,GY,X,Y,NCOR,KNUM)
      AC(1)=PROP(5,MAT)
      BC(1)=XPROP(3,MAT)
      CC(1)=XPROP(9,MAT)
      IF (PROP(8,MAT).EQ.-1.0) CALL USREXT (AC,TP,V1P,V2P,NCOR,MAT)
      IF (PROP(8,MAT).EQ.-1.0.AND.IVAR1.EQ.1) CALL USREX (BC,TP,V1P,    
     1V2P,NCOR,MAT,1)
      IF (PROP(8,MAT).EQ.-1.0.AND.IVAR2.EQ.1) CALL USREX (CC,TP,V1P,    
     1V2P,NCOR,MAT,2)
      DO 210 I=1,NN
      DO 210 K=1,NN
      FD=0.
      FDX=0.
      FDY=0.
      DO 200 KL=1,KLEND
      FD=FD+AC(KL)*FBDY(I,K,KL)
      FDX=FDX+BC(KL)*FBDY(I,K,KL)
      FDY=FDY+CC(KL)*FBDY(I,K,KL)
  200 CONTINUE
      FBDY(I,K,1)=FD*RHO
      FBDY(I,K,2)=FDX*RHO
      FBDY(I,K,3)=FDY*RHO
  210 CONTINUE
      DO 230 I=1,NN
      FDX=0.
      FDY=0.
      DO 220 J=1,NN
      FDX=FDX-FBDY(I,J,1)*GX*TREF
      FDY=FDY-FBDY(I,J,1)*GY*TREF
  220 CONTINUE
      FUV(I)=FUV(I)+FDX
      FUV(I+NN)=FUV(I+NN)+FDY
  230 CONTINUE
      IF (IVAR1.EQ.0) GO TO 280
      DO 250 I=1,NN
      FDX=0.
      FDY=0.
      DO 240 J=1,NN
      FDX=FDX-FBDY(I,J,2)*GX*V1REF
      FDY=FDY-FBDY(I,J,2)*GY*V1REF
  240 CONTINUE
      FUV(I)=FUV(I)+FDX
      FUV(I+NN)=FUV(I+NN)+FDY
  250 CONTINUE
      IF (IVAR2.EQ.0) GO TO 280
      DO 270 I=1,NN
      FDX=0.
      FDY=0.
      DO 260 J=1,NN
      FDX=FDX-FBDY(I,J,3)*GX*V2REF
      FDY=FDY-FBDY(I,J,3)*GY*V2REF
  260 CONTINUE
      FUV(I)=FUV(I)+FDX
      FUV(I+NN)=FUV(I+NN)+FDY
  270 CONTINUE
  280 CONTINUE
C
      RETURN
      END
