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
      SUBROUTINE ASSMVI (KNUM,IEQN,KIND,MAT,NN,NCOR,IBCPT2,CV,CX,CY,    
     1KV,FV,FVS,X,Y,UP,VP,PP,TP,V1P,V2P,CXYV,CXYVP)
C
C     ******************************************************************
C
C     SUBROUTINE TO ASSEMBLE ELEMENT MATRICES FOR AUXILIARY
C     TRANSPORT EQUATIONS
C
C     ******************************************************************
C
      REAL KV
C
      COMMON /RSTART/ IRSTRT,NSTEPS
      COMMON /MATDAT/ NMAT,NPROP,PROP(15,10),NXMAT,NXPROP,XPROP(15,10)
      COMMON /PRBDAT/ IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,    
     1                IPNLTY,PNLTY
      COMMON /SLNDAT/ IALGOR,IBLDUV,IBLDT,IJACOB,IAUTO,NCYCLM,NPRNT,    
     1                KSTEP,KPRNT,TIME,DELTN,DELTNM
C
      DIMENSION KV(9,9,4), CX(9,9,9), CY(9,9,9), FV(9), FVS(9,4)
      DIMENSION CV(9,9)
      DIMENSION CXYV(9,9), CXYVP(9,18)
      DIMENSION UP(9), VP(9), PP(9), TP(9), V1P(9), V2P(9), X(9), Y(9)
      DIMENSION VTEMP(9)
      DIMENSION AC(4), S(4)
C
      PARAMETER (BIG=1.0E40)
C
C     ******************************************************************
C
C     SELECT AUXILIARY VARIABLE
C
      DO 10 I=1,NN
      IF (IEQN.EQ.1) THEN
      VTEMP(I)=V1P(I)
       ELSE
      VTEMP(I)=V2P(I)
      END IF
   10 CONTINUE
C
C     EVALUATE DIFFUSIVITY
C
      IPT=2+(IEQN-1)*6
      IF (PROP(8,MAT).EQ.1.0) GO TO 40
      IF (KSTEP.EQ.1.AND.ITMDEP.EQ.0) GO TO 20
      CALL USRDIF (AC,TP,V1P,V2P,NCOR,MAT,IEQN)
      KLEND=NCOR
      GO TO 60
   20 CONTINUE
      DO 30 I=1,NCOR
      AC(I)=XPROP(IPT,MAT)
   30 CONTINUE
      KLEND=NCOR
      GO TO 60
   40 CONTINUE
      DO 50 I=1,NCOR
      AC(I)=0.
   50 CONTINUE
      AC(1)=XPROP(IPT,MAT)
      KLEND=1
   60 CONTINUE
C
C     CONSTRUCT DIFFUSION MATRIX
C
      DO 80 I=1,NN
      DO 80 J=1,NN
      FD=0.
      DO 70 KL=1,KLEND
      FD=FD+AC(KL)*KV(I,J,KL)
   70 CONTINUE
      KV(I,J,1)=FD
   80 CONTINUE
C
C     CHECK BOUNDARY CONDITIONS
C
      DO 90 I=1,NN
      FD=ABS(KV(I,I,1))
      IF (FD.GT.1.0E15) KV(I,I,1)=BIG
   90 CONTINUE
C
C     CONSTRUCT ADVECTION MATRIX
C
      IPT=1+(IEQN-1)*6
      RHOCP=XPROP(IPT,MAT)
      DO 110 I=1,NN
      DO 110 K=1,NN
      FDX=0.
      FDY=0.
      DO 100 J=1,NN
      FDX=FDX+CX(I,J,K)*UP(J)
      FDY=FDY+CY(I,J,K)*VP(J)
  100 CONTINUE
      CXYV(I,K)=RHOCP*(FDX+FDY)
  110 CONTINUE
C
C     CONSTRUCT JACOBIAN TERMS, IF REQUIRED
C
      DO 120 I=1,9
      DO 120 J=1,18
      CXYVP(I,J)=0.
  120 CONTINUE
      IF (IJACOB.EQ.0) GO TO 150
      DO 140 I=1,NN
      DO 140 K=1,NN
      FDX=0.
      FDY=0.
      DO 130 J=1,NN
      FDX=FDX+CX(I,K,J)*VTEMP(J)
      FDY=FDY+CY(I,K,J)*VTEMP(J)
  130 CONTINUE
      CXYVP(I,K)=RHOCP*FDX
      CXYVP(I,K+NN)=RHOCP*FDY
  140 CONTINUE
  150 CONTINUE
C
C     CONSTRUCT CAPACITANCE MATRIX, IF REQUIRED
C
      IF (ITMDEP.EQ.0) GO TO 170
      DO 160 I=1,NN
      DO 160 J=1,NN
      CV(I,J)=RHOCP*CV(I,J)
  160 CONTINUE
  170 CONTINUE
C
C     UPDATE VARIABLE RHS, IF REQUIRED
C
      MM=IBCPT2
      IF (IEQN.EQ.1) THEN
      ITIME=MOD(MM,1000000)/100000
      IF (ITIME.EQ.0) GO TO 180
      ICURVE=MOD(MM,100000)/10000
      ISIDE=MOD(MM,10000)/1000
      CALL BCTIME (KNUM,5,ICURVE,ISIDE,ITIME,TIME,KSTEP,KIND,X,Y,       
     1UP,VP,PP,TP,V1P,V2P,FV)
  180 CONTINUE
       ELSE
      ITIME=MOD(MM,1000)/100
      IF (ITIME.EQ.0) GO TO 190
      ICURVE=MOD(MM,100)/10
      ISIDE=MOD(MM,10)/1
      CALL BCTIME (KNUM,6,ICURVE,ISIDE,ITIME,TIME,KSTEP,KIND,X,Y,       
     1UP,VP,PP,TP,V1P,V2P,FV)
  190 CONTINUE
      END IF
C
C     CONSTRUCT SOURCE TERMS
C
      IPT=4+(IEQN-1)*6
      DO 200 I=1,NCOR
      AC(I)=XPROP(IPT,MAT)
  200 CONTINUE
      IF (XPROP(IPT,MAT).EQ.-1.0) CALL USRVS (AC,TP,V1P,V2P,X,Y,        
     1NCOR,MAT,KNUM,TIME,KSTEP,IEQN)
      DO 220 I=1,NN
      FD=0.
      DO 210 J=1,NCOR
      FD=FD+FVS(I,J)*AC(J)
  210 CONTINUE
      FV(I)=FV(I)+FD
  220 CONTINUE
C
      RETURN
      END
