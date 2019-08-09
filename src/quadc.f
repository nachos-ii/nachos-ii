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
      SUBROUTINE QUADC (IMAP,MAT,NGAUSS,GS1,GS2,GSWT,X,Y,KUV,CX,CY,MASS,
     1FBDY,KT,FTS,KV1,KV2,KP,KPU,CPOR)
C
C     ******************************************************************
C
C     SUBROUTINE FOR THE ASSEMBLY OF THE DIFFUSION,CONVECTION,MASS,
C     PRESSURE AND SOURCE MATRICES FOR SUB AND ISOPARAMETRIC NINE-NODE
C     (LAGRANGE) QUADS
C
C     ******************************************************************
C
      REAL MASS,KUV,KP,KPU,KT,KV1,KV2,MASSP
C
      COMMON /MATDAT/ NMAT,NPROP,PROP(15,10),NXMAT,NXPROP,XPROP(15,10)
      COMMON /PRBDAT/ IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,    
     1                IPNLTY,PNLTY
      COMMON /QUAD9/  F9Q(9,9),F9L(4,9),DF9QDS(9,9),DF9QDT(9,9)
C
      DIMENSION GS1(*), GS2(*), GSWT(*)
      DIMENSION X(*), Y(*)
      DIMENSION KUV(22,22,4), CX(9,9,9), CY(9,9,9), CPOR(9,9,9)
      DIMENSION KP(18,18), KPU(3,18), MASS(9,9)
      DIMENSION KT(9,9,4), FTS(9,4)
      DIMENSION KV1(9,9,4), KV2(9,9,4)
      DIMENSION FBDY(9,9,4), FF(9,9,4)
      DIMENSION FL(4,9), MASSP(3,3), QT(3,18)
C
C     STATEMENT FUNCTIONS FOR THE JACOBIAN OF THE PARAMETRIC MAPPING
C
      F11Q(S,T)=.25*(D1+S*D2+T*D3+S*T*D4+S**2*D5+(S**2)*T*D6)
      F12Q(S,T)=-.25*(B1+T*B2+S*B3+S*T*B4+T**2*B5+S*(T**2)*B6)
      F21Q(S,T)=-.25*(C1+S*C2+T*C3+S*T*C4+S**2*C5+(S**2)*T*C6)
      F22Q(S,T)=.25*(A1+T*A2+S*A3+S*T*A4+T**2*A5+S*(T**2)*A6)
      DETJQ(S,T)=F11Q(S,T)*F22Q(S,T)-F12Q(S,T)*F21Q(S,T)
C
      RAD9L(L)=F9L(1,L)*X(1)+F9L(2,L)*X(2)+F9L(3,L)*X(3)+F9L(4,L)*X(4)
      RAD9Q(L)=F9Q(1,L)*X(1)+F9Q(2,L)*X(2)+F9Q(3,L)*X(3)+F9Q(4,L)*X(4)+ 
     1F9Q(5,L)*X(5)+F9Q(6,L)*X(6)+F9Q(7,L)*X(7)+F9Q(8,L)*X(8)+          
     2F9Q(9,L)*X(9)
C
      X9L(L)=F9L(1,L)*X(1)+F9L(2,L)*X(2)+F9L(3,L)*X(3)+F9L(4,L)*X(4)
      Y9L(L)=F9L(1,L)*Y(1)+F9L(2,L)*Y(2)+F9L(3,L)*Y(3)+F9L(4,L)*Y(4)
      X9Q(L)=F9Q(1,L)*X(1)+F9Q(2,L)*X(2)+F9Q(3,L)*X(3)+F9Q(4,L)*X(4)+   
     1F9Q(5,L)*X(5)+F9Q(6,L)*X(6)+F9Q(7,L)*X(7)+F9Q(8,L)*X(8)+          
     2F9Q(9,L)*X(9)
      Y9Q(L)=F9Q(1,L)*Y(1)+F9Q(2,L)*Y(2)+F9Q(3,L)*Y(3)+F9Q(4,L)*Y(4)+   
     1F9Q(5,L)*Y(5)+F9Q(6,L)*Y(6)+F9Q(7,L)*Y(7)+F9Q(8,L)*Y(8)+          
     2F9Q(9,L)*Y(9)
C
C     ******************************************************************
C
C     INITIALIZE MATRICES
C
      DO 10 I=1,3
      DO 10 J=1,18
      KPU(I,J)=0.
      QT(I,J)=0.
   10 CONTINUE
      DO 20 I=1,9
      DO 20 J=1,4
      FTS(I,J)=0.
   20 CONTINUE
      DO 40 I=1,9
      DO 40 J=1,9
      MASS(I,J)=0.
      DO 30 K=1,9
      CX(I,J,K)=0.
      CY(I,J,K)=0.
      CPOR(I,J,K)=0.
   30 CONTINUE
      DO 40 K=1,4
      KT(I,J,K)=0.
      FBDY(I,J,K)=0.
      FF(I,J,K)=0.
      KV1(I,J,K)=0.
      KV2(I,J,K)=0.
   40 CONTINUE
      DO 50 I=1,18
      DO 50 J=1,18
      KP(I,J)=0.
   50 CONTINUE
      DO 60 I=1,22
      DO 60 J=1,22
      DO 60 K=1,4
      KUV(I,J,K)=0.
   60 CONTINUE
C
C     SET GEOMETRIC PARAMETERS FOR PARAMETRIC MAPPING
C     EVALUATE DISCONTINUOUS PRESSURE FUNCTIONS
C
      IM=IMAP
      IF (IMAP.EQ.1) THEN
      A1=-X(1)+X(2)+X(3)-X(4)
      A2=X(1)-X(2)+X(3)-X(4)
      A3=0.
      A4=0.
      A5=0.
      A6=0.
      B1=-Y(1)+Y(2)+Y(3)-Y(4)
      B2=Y(1)-Y(2)+Y(3)-Y(4)
      B3=0.
      B4=0.
      B5=0.
      B6=0.
      C1=-X(1)-X(2)+X(3)+X(4)
      C2=X(1)-X(2)+X(3)-X(4)
      C3=0.
      C4=0.
      C5=0.
      C6=0.
      D1=-Y(1)-Y(2)+Y(3)+Y(4)
      D2=Y(1)-Y(2)+Y(3)-Y(4)
      D3=0.
      D4=0.
      D5=0.
      D6=0.
      DO 70 I=1,9
      FL(1,I)=1.0
      FL(2,I)=X9L(I)
      FL(3,I)=Y9L(I)
      FL(4,I)=0.
   70 CONTINUE
       ELSE
      A1=2.*(X(6)-X(8))
      A2=X(1)-X(2)+X(3)-X(4)
      A3=4.*(X(6)+X(8)-2.*X(9))
      A4=2.*(-X(1)-X(2)+X(3)+X(4)+2.*X(5)-2.*X(7))
      A5=-X(1)+X(2)+X(3)-X(4)-2.*X(6)+2.*X(8)
      A6=2.*(X(1)+X(2)+X(3)+X(4)-2.*(X(5)+X(6)+X(7)+X(8))+4.*X(9))
      B1=2.*(Y(6)-Y(8))
      B2=Y(1)-Y(2)+Y(3)-Y(4)
      B3=4.*(Y(6)+Y(8)-2.*Y(9))
      B4=2.*(-Y(1)-Y(2)+Y(3)+Y(4)+2.*Y(5)-2.*Y(7))
      B5=-Y(1)+Y(2)+Y(3)-Y(4)-2.*Y(6)+2.*Y(8)
      B6=2.*(Y(1)+Y(2)+Y(3)+Y(4)-2.*(Y(5)+Y(6)+Y(7)+Y(8))+4.*Y(9))
      C1=2.*(X(7)-X(5))
      C2=X(1)-X(2)+X(3)-X(4)
      C3=4.*(X(5)+X(7)-2.*X(9))
      C4=2.*(-X(1)+X(2)+X(3)-X(4)-2.*X(6)+2.*X(8))
      C5=-X(1)-X(2)+X(3)+X(4)+2.*X(5)-2.*X(7)
      C6=2.*(X(1)+X(2)+X(3)+X(4)-2.*(X(5)+X(6)+X(7)+X(8))+4.*X(9))
      D1=2.*(Y(7)-Y(5))
      D2=Y(1)-Y(2)+Y(3)-Y(4)
      D3=4.*(Y(5)+Y(7)-2.*Y(9))
      D4=2.*(-Y(1)+Y(2)+Y(3)-Y(4)-2.*Y(6)+2.*Y(8))
      D5=-Y(1)-Y(2)+Y(3)+Y(4)+2.*Y(5)-2.*Y(7)
      D6=2.*(Y(1)+Y(2)+Y(3)+Y(4)-2.*(Y(5)+Y(6)+Y(7)+Y(8))+4.*Y(9))
      DO 80 I=1,9
      FL(1,I)=1.0
      FL(2,I)=X9Q(I)
      FL(3,I)=Y9Q(I)
      FL(4,I)=0.
   80 CONTINUE
      END IF
C
C     LOOP ON INTEGRATION POINTS
C
      L=1
      KLEND=1
      IF (PROP(8,MAT).EQ.-1.0) KLEND=4
      DO 230 K=1,NGAUSS
      DO 230 J=1,NGAUSS
      A=F11Q(GS1(J),GS2(K))
      B=F21Q(GS1(J),GS2(K))
      C=F12Q(GS1(J),GS2(K))
      D=F22Q(GS1(J),GS2(K))
      E=DETJQ(GS1(J),GS2(K))
      IF (ABS(E).LE.1.0E-15) CALL ERROR ('QUADC','ZERO JACOBIAN' ,' ',0,
     1' ',0,' ',' ',1)
      EI=1./E
      IF (IMAP.EQ.1) THEN
      R=RAD9L(L)
       ELSE
      R=RAD9Q(L)
      END IF
      IF (EI.LT.0) IM=-1
      IF (IAXSYM.EQ.0) R=1.0
      WW=GSWT(J)*GSWT(K)
C
C     SOURCE TERM
C
      DO 90 I=1,9
      DO 90 M=1,4
      FTS(I,M)=FTS(I,M)+WW*F9Q(I,L)*F9L(M,L)*R*E
   90 CONTINUE
C
C     VISCOUS AND BODY FORCE TERMS
C
      DO 100 KL=1,KLEND
      CC=1.
      IF (PROP(8,MAT).EQ.-1.0) CC=F9L(KL,L)
      DO 100 I=1,9
      DO 100 M=I,9
      AA=A*DF9QDS(I,L)+C*DF9QDT(I,L)
      BB=B*DF9QDS(I,L)+D*DF9QDT(I,L)
      KUV(I,M,KL)=KUV(I,M,KL)+WW*CC*((AA*(A*DF9QDS(M,L)+                
     1C*DF9QDT(M,L)))*EI)*R
      KUV(I+9,M,KL)=KUV(I+9,M,KL)+WW*CC*((AA*(B*DF9QDS(M,L)+            
     1D*DF9QDT(M,L)))*EI)*R
      KUV(I,M+9,KL)=KUV(I,M+9,KL)+WW*CC*((BB*(A*DF9QDS(M,L)+            
     1C*DF9QDT(M,L)))*EI)*R
      KUV(I+9,M+9,KL)=KUV(I+9,M+9,KL)+WW*CC*((BB*(B*DF9QDS(M,L)+        
     1D*DF9QDT(M,L)))*EI)*R
      FF(I,M,KL)=FF(I,M,KL)+WW*CC*(F9Q(I,L)*F9Q(M,L))*E/R
      FBDY(I,M,KL)=FBDY(I,M,KL)+WW*CC*(F9Q(I,L)*F9Q(M,L))*E*R
  100 CONTINUE
C
C     PRESSURE TERMS
C
      IF (IPNLTY.EQ.1) GO TO 130
      IF (IPFUNC.EQ.0) THEN
      DO 110 I=1,4
      DO 110 M=1,9
      FL(I,M)=F9L(I,M)
  110 CONTINUE
      END IF
      DO 120 I=1,9
      DO 120 M=1,4
      AA=FL(M,L)
      KUV(I,M+18,1)=KUV(I,M+18,1)-WW*((A*DF9QDS(I,L)+                   
     1C*DF9QDT(I,L))*AA)*R-WW*(F9Q(I,L)*FL(M,L))*E*IAXSYM
      KUV(I+9,M+18,1)=KUV(I+9,M+18,1)-WW*((B*DF9QDS(I,L)+               
     1D*DF9QDT(I,L))*AA)*R
      KUV(M+18,I,1)=KUV(I,M+18,1)
      KUV(M+18,I+9,1)=KUV(I+9,M+18,1)
  120 CONTINUE
  130 CONTINUE
C
C     CONVECTIVE TERMS
C
      DO 140 I=1,9
      DO 140 M=1,9
      DO 140 N=1,9
      AA=F9Q(I,L)*F9Q(M,L)
      CX(I,M,N)=CX(I,M,N)+WW*AA*(A*DF9QDS(N,L)+C*DF9QDT(N,L))*R
      CY(I,M,N)=CY(I,M,N)+WW*AA*(B*DF9QDS(N,L)+D*DF9QDT(N,L))*R
      CPOR(I,M,N)=CPOR(I,M,N)+WW*AA*F9Q(N,L)*E*R
  140 CONTINUE
C
C     MASS (CAPACITANCE) TERMS
C
      DO 150 I=1,9
      DO 150 M=I,9
      MASS(I,M)=MASS(I,M)+WW*(F9Q(I,L)*F9Q(M,L))*E*R
  150 CONTINUE
C
C     PENALTY TERMS
C
      IF (IPNLTY.EQ.0) GO TO 220
      DO 160 I=1,3
      DO 160 M=1,3
      MASSP(I,M)=MASSP(I,M)+WW*(FL(I,L)*FL(M,L))*E*R
  160 CONTINUE
      CALL INVRT3 (MASSP,DET)
      DO 170 M=1,3
      DO 170 I=1,9
      AA=FL(M,L)
      QT(M,I)=QT(M,I)+WW*AA*(A*DF9QDS(I,L)+C*DF9QDT(I,L))*R+            
     1WW*(F9Q(I,L)*FL(M,L))*E*IAXSYM
      QT(M,I+9)=QT(M,I+9)+WW*AA*(B*DF9QDS(I,L)+D*DF9QDT(I,L))*R
  170 CONTINUE
      DO 190 I=1,3
      DO 190 N=1,18
      AA=0.
      DO 180 M=1,3
      AA=AA+MASSP(I,M)*QT(M,N)
  180 CONTINUE
      KPU(I,N)=AA
  190 CONTINUE
      DO 210 I=1,18
      DO 210 N=1,18
      AA=0.
      DO 200 M=1,3
      AA=AA+QT(M,I)*KPU(M,N)
  200 CONTINUE
      KP(I,N)=AA
  210 CONTINUE
C
  220 CONTINUE
      L=L+1
  230 CONTINUE
C
C     CONSTRUCT REMAINDER OF SYMMETRIC MATRICES
C
      DO 240 I=1,9
      DO 240 M=I,9
      IF (I.EQ.M) GO TO 240
      MASS(M,I)=MASS(I,M)
  240 CONTINUE
      DO 250 KL=1,KLEND
      DO 250 I=1,9
      DO 250 M=I,9
      IF (M.EQ.I) GO TO 250
      KUV(M,I,KL)=KUV(I,M,KL)
      KUV(M+9,I,KL)=KUV(I,M+9,KL)
      KUV(M,I+9,KL)=KUV(I+9,M,KL)
      KUV(M+9,I+9,KL)=KUV(I+9,M+9,KL)
      FF(M,I,KL)=FF(I,M,KL)
      FBDY(M,I,KL)=FBDY(I,M,KL)
  250 CONTINUE
C
C     RE-ARRANGE VISCOUS TERMS AND CONSTRUCT REMAINING DIFFUSION TERMS
C
      DO 260 KL=1,KLEND
      DO 260 I=1,9
      DO 260 M=1,9
      AA=KUV(I,M,KL)
      BB=KUV(I+9,M+9,KL)
      KUV(I,M,KL)=2.*AA+BB+2.*FF(I,M,KL)*IAXSYM
      KUV(I+9,M+9,KL)=2.*BB+AA
C
      KT(I,M,KL)=AA+BB
C
      KV1(I,M,KL)=AA+BB
C
      KV2(I,M,KL)=AA+BB
  260 CONTINUE
C
      IMAP=IM
      RETURN
      END
