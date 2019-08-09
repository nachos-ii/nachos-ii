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
      SUBROUTINE ASSMU2 (KNUM,KIND,MAT,NN,NCOR,IBCPT1,MASS,C,KUV,       
     1FUV,KP,FBDY,X,Y,UP,VP,PP,TP,V1P,V2P,CXY,CXYP,GX,GY)
C
C     ******************************************************************
C
C     SUBROUTINE TO ASSEMBLE ELEMENT MATRICES FOR MOMENTUM AND
C     CONTINUITY EQUATIONS (DARCY-BRINKMAN EQUATIONS)
C
C     ******************************************************************
C
      REAL MASS,KUV,KP,KDIAG
C
      COMMON /RSTART/ IRSTRT,NSTEPS
      COMMON /MATDAT/ NMAT,NPROP,PROP(15,10),NXMAT,NXPROP,XPROP(15,10)
      COMMON /PRBDAT/ IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,    
     1                IPNLTY,PNLTY
      COMMON /SLNDAT/ IALGOR,IBLDUV,IBLDT,IJACOB,IAUTO,NCYCLM,NPRNT,    
     1                KSTEP,KPRNT,TIME,DELTN,DELTNM
C
      DIMENSION KUV(22,22,4), C(9,9,9), FUV(22)
      DIMENSION FBDY(9,9,4), MASS(9,9), KP(18,18)
      DIMENSION CXY(9,9), CXYP(18,18)
      DIMENSION KDIAG(22)
      DIMENSION UP(9), VP(9), PP(9), TP(9), V1P(9), V2P(9), X(9), Y(9)
      DIMENSION AC(4), BC(4), CC(4), S(4)
C
      PARAMETER (BIG=1.0E40)
      PARAMETER (EPS=1.0E-25)
C
C     ******************************************************************
C
C     EVALUATE EFFECTIVE FLUID VISCOSITY
C
      IF (PROP(8,MAT).EQ.1.0) GO TO 30
      IF (KSTEP.EQ.1.AND.ITMDEP.EQ.0.AND.IRSTRT.EQ.0) GO TO 10
      IF (PROP(15,MAT).EQ.3.0) GO TO 10
      CALL USRVIS (AC,TP,S,V1P,V2P,NCOR,MAT,KSTEP)
      KLEND=NCOR
      GO TO 50
   10 CONTINUE
      DO 20 I=1,NCOR
      AC(I)=PROP(3,MAT)
   20 CONTINUE
      KLEND=NCOR
      GO TO 50
   30 CONTINUE
      DO 40 I=1,NCOR
      AC(I)=0.
   40 CONTINUE
      AC(1)=PROP(3,MAT)
      KLEND=1
   50 CONTINUE
C
C     CONSTRUCT DIFFUSION MATRIX
C
      NROW=2*NN
      DO 70 I=1,NROW
      KDIAG(I)=KUV(I,I,1)
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
      FD=ABS(KDIAG(I))
      IF (FD.GT.1.0E15) KUV(I,I,1)=BIG
   80 CONTINUE
C
C     CONSTRUCT ADVECTION MATRICES
C
      MATF=INT(PROP(1,MAT))
      RHOC=PROP(1,MATF)*PROP(5,MAT)/SQRT(PROP(2,MAT))
      IF (PROP(8,MATF).EQ.-1.0) THEN
      IF (PROP(15,MATF).EQ.2.0) CALL DEFORM (KIND,X,Y,UP,VP,S)
      CALL USRVIS (AC,TP,S,V1P,V2P,NCOR,MATF,KSTEP)
      SUM=0.
      DO 90 I=1,NCOR
      SUM=SUM+AC(I)
   90 CONTINUE
      AMU=(SUM/NCOR)/PROP(2,MAT)
       ELSE
      AMU=PROP(2,MATF)/PROP(2,MAT)
      END IF
      DO 100 I=1,NN
      DO 100 J=1,NN
      KUV(I,J,1)=KUV(I,J,1)+AMU*MASS(I,J)
      KUV(I+NN,J+NN,1)=KUV(I+NN,J+NN,1)+AMU*MASS(I,J)
  100 CONTINUE
      DO 120 I=1,NN
      DO 120 K=1,NN
      FDX=0.
      DO 110 J=1,NN
      FD=SQRT(UP(J)**2+VP(J)**2)
      FDX=FDX+C(I,J,K)*FD
  110 CONTINUE
      CXY(I,K)=RHOC*FDX
  120 CONTINUE
C
C     CONSTRUCT JACOBIAN TERMS, IF REQUIRED
C
      DO 130 I=1,18
      DO 130 K=1,18
      CXYP(I,K)=0.
  130 CONTINUE
      IF (IJACOB.EQ.0) GO TO 160
      IF (KSTEP.EQ.1) GO TO 160
      DO 150 I=1,NN
      DO 150 K=1,NN
      FDX=0.
      FDY=0.
      FDXY=0.
      DO 140 J=1,NN
      FD=SQRT(UP(J)**2+VP(J)**2+EPS**2)
      FDX=FDX+C(I,K,J)*UP(J)*UP(J)/FD
      FDY=FDY+C(I,K,J)*VP(J)*VP(J)/FD
      FDXY=FDXY+C(I,K,J)*UP(J)*VP(J)/FD
  140 CONTINUE
      CXYP(I,K)=RHOC*FDX
      CXYP(I+NN,K+NN)=RHOC*FDY
      CXYP(I,K+NN)=RHOC*FDXY
      CXYP(I+NN,K)=RHOC*FDXY
  150 CONTINUE
  160 CONTINUE
C
C     CONSTRUCT MASS MATRIX, IF REQUIRED
C
      IF (ITMDEP.EQ.0) GO TO 180
      PHI=1.0/PROP(4,MAT)
      DO 170 I=1,NN
      DO 170 J=1,NN
      MASS(I,J)=PHI*MASS(I,J)
  170 CONTINUE
  180 CONTINUE
C
C     UPDATE VARIABLE RHS, IF REQUIRED
C
      MM=IBCPT1
      ITIME=MM/100000000
      IF (ITIME.EQ.0) GO TO 190
      ICURVE=MOD(MM,100000000)/10000000
      ISIDE=MOD(MM,10000000)/1000000
      CALL BCTIME (KNUM,1,ICURVE,ISIDE,ITIME,TIME,KSTEP,KIND,X,Y,       
     1UP,VP,PP,TP,V1P,V2P,FUV)
  190 CONTINUE
      ITIME=MOD(MM,1000000)/100000
      IF (ITIME.EQ.0) GO TO 200
      ICURVE=MOD(MM,100000)/10000
      ISIDE=MOD(MM,10000)/1000
      CALL BCTIME (KNUM,2,ICURVE,ISIDE,ITIME,TIME,KSTEP,KIND,X,Y,       
     1UP,VP,PP,TP,V1P,V2P,FUV)
  200 CONTINUE
      ITIME=MOD(MM,1000)/100
      IF (ITIME.EQ.0) GO TO 210
      ICURVE=MOD(MM,100)/10
      ISIDE=MOD(MM,10)/1
      CALL BCTIME (KNUM,3,ICURVE,ISIDE,ITIME,TIME,KSTEP,KIND,X,Y,       
     1UP,VP,PP,TP,V1P,V2P,FUV)
  210 CONTINUE
C
C     CONSTRUCT BUOYANCY TERMS, IF REQUIRED
C
      RHO=PROP(1,MATF)
      IF (IFREE.EQ.0) GO TO 300
      TREF=PROP(11,MATF)
      V1REF=XPROP(5,MATF)
      V2REF=XPROP(11,MATF)
      GX=PROP(6,MATF)
      GY=PROP(7,MATF)
      IF (GX.EQ.-1.0) CALL USRGRV (GX,GY,X,Y,NCOR,KNUM)
      AC(1)=PROP(5,MATF)
      BC(1)=XPROP(3,MATF)
      CC(1)=XPROP(9,MATF)
      IF (PROP(8,MATF).EQ.-1.0) CALL USREXT (AC,TP,V1P,V2P,NCOR,MATF)
      IF (PROP(8,MATF).EQ.-1.0.AND.IVAR1.EQ.1) CALL USREX (BC,TP,V1P,   
     1V2P,NCOR,MATF,1)
      IF (PROP(8,MATF).EQ.-1.0.AND.IVAR2.EQ.1) CALL USREX (CC,TP,V1P,   
     1V2P,NCOR,MATF,2)
      DO 230 I=1,NN
      DO 230 K=1,NN
      FD=0.
      FDX=0.
      FDY=0.
      DO 220 KL=1,KLEND
      FD=FD+AC(KL)*FBDY(I,K,KL)
      FDX=FDX+BC(KL)*FBDY(I,K,KL)
      FDY=FDY+CC(KL)*FBDY(I,K,KL)
  220 CONTINUE
      FBDY(I,K,1)=FD*RHO
      FBDY(I,K,2)=FDX*RHO
      FBDY(I,K,3)=FDY*RHO
  230 CONTINUE
      DO 250 I=1,NN
      FDX=0.
      FDY=0.
      DO 240 J=1,NN
      FDX=FDX-FBDY(I,J,1)*GX*TREF
      FDY=FDY-FBDY(I,J,1)*GY*TREF
  240 CONTINUE
      FUV(I)=FUV(I)+FDX
      FUV(I+NN)=FUV(I+NN)+FDY
  250 CONTINUE
      IF (IVAR1.EQ.0) GO TO 300
      DO 270 I=1,NN
      FDX=0.
      FDY=0.
      DO 260 J=1,NN
      FDX=FDX-FBDY(I,J,2)*GX*V1REF
      FDY=FDY-FBDY(I,J,2)*GY*V1REF
  260 CONTINUE
      FUV(I)=FUV(I)+FDX
      FUV(I+NN)=FUV(I+NN)+FDY
  270 CONTINUE
      IF (IVAR2.EQ.0) GO TO 300
      DO 290 I=1,NN
      FDX=0.
      FDY=0.
      DO 280 J=1,NN
      FDX=FDX-FBDY(I,J,3)*GX*V2REF
      FDY=FDY-FBDY(I,J,3)*GY*V2REF
  280 CONTINUE
      FUV(I)=FUV(I)+FDX
      FUV(I+NN)=FUV(I+NN)+FDY
  290 CONTINUE
  300 CONTINUE
C
      RETURN
      END
