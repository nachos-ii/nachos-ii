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
      SUBROUTINE QUADB (IMAP,MAT,NGAUSS,GS1,GS2,GSWT,X,Y,KUV,CX,CY,MASS,
     1FBDY,KT,FTS,KV1,KV2,KP,KPU,CPOR)
C
C     ******************************************************************
C
C     SUBROUTINE FOR THE ASSEMBLY OF THE DIFFUSION,CONVECTION,MASS,
C     PRESSURE AND SOURCE MATRICES FOR SUB AND ISOPARAMETRIC EIGHT-NODE
C     (SERENDIPITY) QUADS
C
C     ******************************************************************
C
      REAL MASS,KUV,KP,KPU,KT,KV1,KV2,MASSP
C
      COMMON /MATDAT/ NMAT,NPROP,PROP(15,10),NXMAT,NXPROP,XPROP(15,10)
      COMMON /PRBDAT/ IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,    
     1                IPNLTY,PNLTY
      COMMON /QUAD8/  F8Q(8,9),F8L(4,9),DF8QDS(8,9),DF8QDT(8,9)
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
      F11Q(S,T)=.25*(D1+S*D2+T*D3+S*T*D4+S**2*D5)
      F12Q(S,T)=-.25*(B1+T*B2+S*B3+S*T*B4+T**2*B5)
      F21Q(S,T)=-.25*(C1+S*C2+T*C3+S*T*C4+S**2*C5)
      F22Q(S,T)=.25*(A1+T*A2+S*A3+S*T*A4+T**2*A5)
      DETJQ(S,T)=F11Q(S,T)*F22Q(S,T)-F12Q(S,T)*F21Q(S,T)
C
      RAD8L(L)=F8L(1,L)*X(1)+F8L(2,L)*X(2)+F8L(3,L)*X(3)+F8L(4,L)*X(4)
      RAD8Q(L)=F8Q(1,L)*X(1)+F8Q(2,L)*X(2)+F8Q(3,L)*X(3)+F8Q(4,L)*X(4)+ 
     1F8Q(5,L)*X(5)+F8Q(6,L)*X(6)+F8Q(7,L)*X(7)+F8Q(8,L)*X(8)
C
      X8L(L)=F8L(1,L)*X(1)+F8L(2,L)*X(2)+F8L(3,L)*X(3)+F8L(4,L)*X(4)
      Y8L(L)=F8L(1,L)*Y(1)+F8L(2,L)*Y(2)+F8L(3,L)*Y(3)+F8L(4,L)*Y(4)
      X8Q(L)=F8Q(1,L)*X(1)+F8Q(2,L)*X(2)+F8Q(3,L)*X(3)+F8Q(4,L)*X(4)+   
     1F8Q(5,L)*X(5)+F8Q(6,L)*X(6)+F8Q(7,L)*X(7)+F8Q(8,L)*X(8)
      Y8Q(L)=F8Q(1,L)*Y(1)+F8Q(2,L)*Y(2)+F8Q(3,L)*Y(3)+F8Q(4,L)*Y(4)+   
     1F8Q(5,L)*Y(5)+F8Q(6,L)*Y(6)+F8Q(7,L)*Y(7)+F8Q(8,L)*Y(8)
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
C
      IM=IMAP
      IF (IMAP.EQ.1) THEN
      A1=-X(1)+X(2)+X(3)-X(4)
      A2=X(1)-X(2)+X(3)-X(4)
      A3=0.
      A4=0.
      A5=0.
      B1=-Y(1)+Y(2)+Y(3)-Y(4)
      B2=Y(1)-Y(2)+Y(3)-Y(4)
      B3=0.
      B4=0.
      B5=0.
      C1=-X(1)-X(2)+X(3)+X(4)
      C2=X(1)-X(2)+X(3)-X(4)
      C3=0.
      C4=0.
      C5=0.
      D1=-Y(1)-Y(2)+Y(3)+Y(4)
      D2=Y(1)-Y(2)+Y(3)-Y(4)
      D3=0.
      D4=0.
      D5=0.
      DO 70 I=1,9
      FL(1,I)=1.0
      FL(2,I)=X8L(I)
      FL(3,I)=Y8L(I)
      FL(4,I)=0.
   70 CONTINUE
       ELSE
      A1=2.*(X(6)-X(8))
      A2=X(1)-X(2)+X(3)-X(4)
      A3=2.*(X(1)+X(2)+X(3)+X(4)-2.*X(5)-2.*X(7))
      A4=2.*(-X(1)-X(2)+X(3)+X(4)+2.*X(5)-2.*X(7))
      A5=-X(1)+X(2)+X(3)-X(4)-2.*X(6)+2.*X(8)
      B1=2.*(Y(6)-Y(8))
      B2=Y(1)-Y(2)+Y(3)-Y(4)
      B3=2.*(Y(1)+Y(2)+Y(3)+Y(4)-2.*Y(5)-2.*Y(7))
      B4=2.*(-Y(1)-Y(2)+Y(3)+Y(4)+2.*Y(5)-2.*Y(7))
      B5=-Y(1)+Y(2)+Y(3)-Y(4)-2.*Y(6)+2.*Y(8)
      C1=2.*(X(7)-X(5))
      C2=X(1)-X(2)+X(3)-X(4)
      C3=2.*(X(1)+X(2)+X(3)+X(4)-2.*X(6)-2.*X(8))
      C4=2.*(-X(1)+X(2)+X(3)-X(4)-2.*X(6)+2.*X(8))
      C5=-X(1)-X(2)+X(3)+X(4)+2.*X(5)-2.*X(7)
      D1=2.*(Y(7)-Y(5))
      D2=Y(1)-Y(2)+Y(3)-Y(4)
      D3=2.*(Y(1)+Y(2)+Y(3)+Y(4)-2.*Y(6)-2.*Y(8))
      D4=2.*(-Y(1)+Y(2)+Y(3)-Y(4)-2.*Y(6)+2.*Y(8))
      D5=-Y(1)-Y(2)+Y(3)+Y(4)+2.*Y(5)-2.*Y(7)
      DO 80 I=1,9
      FL(1,I)=1.0
      FL(2,I)=X8Q(I)
      FL(3,I)=Y8Q(I)
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
      IF (ABS(E).LE.1.0E-15) CALL ERROR ('QUADB','ZERO JACOBIAN' ,' ',0,
     1' ',0,' ',' ',1)
      EI=1./E
      IF (IMAP.EQ.1) THEN
      R=RAD8L(L)
       ELSE
      R=RAD8Q(L)
      END IF
      IF (EI.LT.0) IM=-1
      IF (IAXSYM.EQ.0) R=1.0
      WW=GSWT(J)*GSWT(K)
C
C     SOURCE TERM
C
      DO 90 I=1,8
      DO 90 M=1,4
      FTS(I,M)=FTS(I,M)+WW*F8Q(I,L)*F8L(M,L)*R*E
   90 CONTINUE
C
C     VISCOUS AND BODY FORCE TERMS
C
      DO 100 KL=1,KLEND
      CC=1.
      IF (PROP(8,MAT).EQ.-1.0) CC=F8L(KL,L)
      DO 100 I=1,8
      DO 100 M=I,8
      AA=A*DF8QDS(I,L)+C*DF8QDT(I,L)
      BB=B*DF8QDS(I,L)+D*DF8QDT(I,L)
      KUV(I,M,KL)=KUV(I,M,KL)+WW*CC*((AA*(A*DF8QDS(M,L)+                
     1C*DF8QDT(M,L)))*EI)*R
      KUV(I+8,M,KL)=KUV(I+8,M,KL)+WW*CC*((AA*(B*DF8QDS(M,L)+            
     1D*DF8QDT(M,L)))*EI)*R
      KUV(I,M+8,KL)=KUV(I,M+8,KL)+WW*CC*((BB*(A*DF8QDS(M,L)+            
     1C*DF8QDT(M,L)))*EI)*R
      KUV(I+8,M+8,KL)=KUV(I+8,M+8,KL)+WW*CC*((BB*(B*DF8QDS(M,L)+        
     1D*DF8QDT(M,L)))*EI)*R
      FF(I,M,KL)=FF(I,M,KL)+WW*CC*(F8Q(I,L)*F8Q(M,L))*E/R
      FBDY(I,M,KL)=FBDY(I,M,KL)+WW*CC*(F8Q(I,L)*F8Q(M,L))*E*R
  100 CONTINUE
C
C     PRESSURE TERMS
C
      IF (IPNLTY.EQ.1) GO TO 130
      IF (IPFUNC.EQ.0) THEN
      DO 110 I=1,4
      DO 110 M=1,9
      FL(I,M)=F8L(I,M)
  110 CONTINUE
      END IF
      DO 120 I=1,8
      DO 120 M=1,4
      AA=FL(M,L)
      KUV(I,M+16,1)=KUV(I,M+16,1)-WW*((A*DF8QDS(I,L)+C*DF8QDT(I,L))     
     1*AA)*R-WW*(F8Q(I,L)*FL(M,L))*E*IAXSYM
      KUV(I+8,M+16,1)=KUV(I+8,M+16,1)-WW*((B*DF8QDS(I,L)+D*DF8QDT(I,L)) 
     1*AA)*R
      KUV(M+16,I,1)=KUV(I,M+16,1)
      KUV(M+16,I+8,1)=KUV(I+8,M+16,1)
  120 CONTINUE
  130 CONTINUE
C
C     CONVECTIVE TERMS
C
      DO 140 I=1,8
      DO 140 M=1,8
      DO 140 N=1,8
      AA=F8Q(I,L)*F8Q(M,L)
      CX(I,M,N)=CX(I,M,N)+WW*AA*(A*DF8QDS(N,L)+C*DF8QDT(N,L))*R
      CY(I,M,N)=CY(I,M,N)+WW*AA*(B*DF8QDS(N,L)+D*DF8QDT(N,L))*R
      CPOR(I,M,N)=CPOR(I,M,N)+WW*AA*F8Q(N,L)*E*R
  140 CONTINUE
C
C     MASS (CAPACITANCE) TERMS
C
      DO 150 I=1,8
      DO 150 M=I,8
      MASS(I,M)=MASS(I,M)+WW*(F8Q(I,L)*F8Q(M,L))*E*R
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
      DO 170 I=1,8
      AA=FL(M,L)
      QT(M,I)=QT(M,I)+WW*AA*(A*DF8QDS(I,L)+C*DF8QDT(I,L))*R+            
     1WW*(F8Q(I,L)*FL(M,L))*E*IAXSYM
      QT(M,I+8)=QT(M,I+8)+WW*AA*(B*DF8QDS(I,L)+D*DF8QDT(I,L))*R
  170 CONTINUE
      DO 190 I=1,3
      DO 190 N=1,16
      AA=0.
      DO 180 M=1,3
      AA=AA+MASSP(I,M)*QT(M,N)
  180 CONTINUE
      KPU(I,N)=AA
  190 CONTINUE
      DO 210 I=1,16
      DO 210 N=1,16
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
      DO 240 I=1,8
      DO 240 M=I,8
      IF (I.EQ.M) GO TO 240
      MASS(M,I)=MASS(I,M)
  240 CONTINUE
      DO 250 KL=1,KLEND
      DO 250 I=1,8
      DO 250 M=I,8
      IF (M.EQ.I) GO TO 250
      KUV(M,I,KL)=KUV(I,M,KL)
      KUV(M+8,I,KL)=KUV(I,M+8,KL)
      KUV(M,I+8,KL)=KUV(I+8,M,KL)
      KUV(M+8,I+8,KL)=KUV(I+8,M+8,KL)
      FF(M,I,KL)=FF(I,M,KL)
      FBDY(M,I,KL)=FBDY(I,M,KL)
  250 CONTINUE
C
C     RE-ARRANGE VISCOUS TERMS AND CONSTRUCT REMAINING DIFFUSION TERMS
C
      DO 260 KL=1,KLEND
      DO 260 I=1,8
      DO 260 M=1,8
      AA=KUV(I,M,KL)
      BB=KUV(I+8,M+8,KL)
      KUV(I,M,KL)=2.*AA+BB+2.*FF(I,M,KL)*IAXSYM
      KUV(I+8,M+8,KL)=2.*BB+AA
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
