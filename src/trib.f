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
      SUBROUTINE TRIB (IMAP,MAT,NHAM,HM1,HM2,HMWT,X,Y,KUV,CX,CY,MASS,   
     1FBDY,KT,FTS,KV1,KV2,KP,KPU,CPOR)
C
C     ******************************************************************
C
C     SUBROUTINE FOR THE ASSEMBLY OF THE DIFFUSION,CONVECTION,MASS,
C     PRESSURE AND SOURCE MATRICES FOR SUB AND ISOPARAMETRIC TRIANGLES
C
C     ******************************************************************
C
      REAL MASS,KUV,KP,KPU,KT,KV1,KV2,MASSP
C
      COMMON /MATDAT/ NMAT,NPROP,PROP(15,10),NXMAT,NXPROP,XPROP(15,10)
      COMMON /PRBDAT/ IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,    
     1                IPNLTY,PNLTY
      COMMON /TRI6/   F6Q(6,7),F6L(3,7),DF6QDS(6,7),DF6QDT(6,7)
C
      DIMENSION HM1(*), HM2(*), HMWT(*)
      DIMENSION X(*), Y(*)
      DIMENSION KUV(22,22,4), CX(9,9,9), CY(9,9,9), CPOR(9,9,9)
      DIMENSION KP(18,18), KPU(3,18), MASS(9,9)
      DIMENSION KT(9,9,4), FTS(9,4)
      DIMENSION KV1(9,9,4), KV2(9,9,4)
      DIMENSION FBDY(9,9,4), FF(9,9,4)
      DIMENSION FL(3,7), MASSP(3,3), QT(3,18)
C
C     STATEMENT FUNCTIONS FOR THE JACOBIAN OF THE PARAMETRIC MAPPING
C
      F11T(S,T)=S*A1+T*A2+A3
      F12T(S,T)=-(S*C1+T*C2+C3)
      F21T(S,T)=-(S*B1+T*B2+B3)
      F22T(S,T)=S*D1+T*D2+D3
      DETJT(S,T)=F11T(S,T)*F22T(S,T)-F12T(S,T)*F21T(S,T)
C
      RAD6L(L)=F6L(1,L)*X(1)+F6L(2,L)*X(2)+F6L(3,L)*X(3)
      RAD6Q(L)=F6Q(1,L)*X(1)+F6Q(2,L)*X(2)+F6Q(3,L)*X(3)+F6Q(4,L)*X(4)+ 
     1F6Q(5,L)*X(5)+F6Q(6,L)*X(6)
C
      X6L(L)=F6L(1,L)*X(1)+F6L(2,L)*X(2)+F6L(3,L)*X(3)
      Y6L(L)=F6L(1,L)*Y(1)+F6L(2,L)*Y(2)+F6L(3,L)*Y(3)
      X6Q(L)=F6Q(1,L)*X(1)+F6Q(2,L)*X(2)+F6Q(3,L)*X(3)+F6Q(4,L)*X(4)+   
     1F6Q(5,L)*X(5)+F6Q(6,L)*X(6)
      Y6Q(L)=F6Q(1,L)*Y(1)+F6Q(2,L)*Y(2)+F6Q(3,L)*Y(3)+F6Q(4,L)*Y(4)+   
     1F6Q(5,L)*Y(5)+F6Q(6,L)*Y(6)
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
      A1=0.
      A2=0.
      A3=Y(2)-Y(3)
      B1=0.
      B2=0.
      B3=X(2)-X(3)
      C1=0.
      C2=0.
      C3=Y(1)-Y(3)
      D1=0.
      D2=0.
      D3=X(1)-X(3)
      DO 70 I=1,7
      FL(1,I)=1.0
      FL(2,I)=X6L(I)
      FL(3,I)=Y6L(I)
   70 CONTINUE
       ELSE
      A1=4.*(Y(3)+Y(4)-Y(5)-Y(6))
      A2=4.*(Y(2)+Y(3)-2.*Y(5))
      A3=4.*Y(5)-3.*Y(3)-Y(2)
      B1=4.*(X(3)+X(4)-X(5)-X(6))
      B2=4.*(X(2)+X(3)-2.*X(5))
      B3=4.*X(5)-3.*X(3)-X(2)
      C1=4.*(Y(1)+Y(3)-2.*Y(6))
      C2=4.*(Y(3)+Y(4)-Y(5)-Y(6))
      C3=4.*Y(6)-3.*Y(3)-Y(1)
      D1=4.*(X(1)+X(3)-2.*X(6))
      D2=4.*(X(3)+X(4)-X(5)-X(6))
      D3=4.*X(6)-3.*X(3)-X(1)
      DO 80 I=1,7
      FL(1,I)=1.0
      FL(2,I)=X6Q(I)
      FL(3,I)=Y6Q(I)
   80 CONTINUE
      END IF
C
C     LOOP ON INTEGRATION POINTS
C
      L=1
      KLEND=1
      IF (PROP(8,MAT).EQ.-1.0) KLEND=3
      DO 230 J=1,NHAM
      A=F11T(HM1(J),HM2(J))
      B=F21T(HM1(J),HM2(J))
      C=F12T(HM1(J),HM2(J))
      D=F22T(HM1(J),HM2(J))
      E=DETJT(HM1(J),HM2(J))
      IF (ABS(E).LE.1.0E-15) CALL ERROR ('TRIB','ZERO JACOBIAN' ,' ',0, 
     1' ',0,' ',' ',1)
      EI=1./E
      IF (IMAP.EQ.1) THEN
      R=RAD6L(L)
       ELSE
      R=RAD6Q(L)
      END IF
      IF (EI.LT.0.) IM=-1
      IF (IAXSYM.EQ.0) R=1.0
      WW=HMWT(J)
C
C     SOURCE TERM
C
      DO 90 I=1,6
      DO 90 M=1,3
      FTS(I,M)=FTS(I,M)+WW*F6Q(I,L)*F6L(M,L)*R*E
   90 CONTINUE
C
C     VISCOUS AND BODY FORCE TERM
C
      DO 100 KL=1,KLEND
      CC=1.
      IF (PROP(8,MAT).EQ.-1.0) CC=F6L(KL,L)
      DO 100 I=1,6
      DO 100 M=I,6
      AA=A*DF6QDS(I,L)+C*DF6QDT(I,L)
      BB=B*DF6QDS(I,L)+D*DF6QDT(I,L)
      KUV(I,M,KL)=KUV(I,M,KL)+WW*CC*((AA*(A*DF6QDS(M,L)+C*DF6QDT(M,L))) 
     1*EI)*R
      KUV(I+6,M,KL)=KUV(I+6,M,KL)+WW*CC*((AA*(B*DF6QDS(M,L)+            
     1D*DF6QDT(M,L)))*EI)*R
      KUV(I,M+6,KL)=KUV(I,M+6,KL)+WW*CC*((BB*(A*DF6QDS(M,L)+            
     1C*DF6QDT(M,L)))*EI)*R
      KUV(I+6,M+6,KL)=KUV(I+6,M+6,KL)+WW*CC*((BB*(B*DF6QDS(M,L)+        
     1D*DF6QDT(M,L)))*EI)*R
      FF(I,M,KL)=FF(I,M,KL)+WW*CC*(F6Q(I,L)*F6Q(M,L))*E/R
      FBDY(I,M,KL)=FBDY(I,M,KL)+WW*CC*(F6Q(I,L)*F6Q(M,L))*E*R
  100 CONTINUE
C
C     PRESSURE TERMS
C
      IF (IPNLTY.EQ.1) GO TO 130
      IF (IPFUNC.EQ.0) THEN
      DO 110 I=1,3
      DO 110 M=1,7
      FL(I,M)=F6L(I,M)
  110 CONTINUE
      END IF
      DO 120 I=1,6
      DO 120 M=1,3
      AA=FL(M,L)
      KUV(I,M+12,1)=KUV(I,M+12,1)-WW*((A*DF6QDS(I,L)+C*DF6QDT(I,L))     
     1*AA)*R-WW*(F6Q(I,L)*FL(M,L))*E*IAXSYM
      KUV(I+6,M+12,1)=KUV(I+6,M+12,1)-WW*((B*DF6QDS(I,L)+D*DF6QDT(I,L)) 
     1*AA)*R
      KUV(M+12,I,1)=KUV(I,M+12,1)
      KUV(M+12,I+6,1)=KUV(I+6,M+12,1)
  120 CONTINUE
  130 CONTINUE
C
C     CONVECTIVE TERMS
C
      DO 140 I=1,6
      DO 140 M=1,6
      DO 140 N=1,6
      AA=F6Q(I,L)*F6Q(M,L)
      CX(I,M,N)=CX(I,M,N)+WW*AA*(A*DF6QDS(N,L)+C*DF6QDT(N,L))*R
      CY(I,M,N)=CY(I,M,N)+WW*AA*(B*DF6QDS(N,L)+D*DF6QDT(N,L))*R
      CPOR(I,M,N)=CPOR(I,M,N)+WW*AA*F6Q(N,L)*E*R
  140 CONTINUE
C
C     MASS (CAPACITANCE) TERMS
C
      DO 150 I=1,6
      DO 150 M=I,6
      MASS(I,M)=MASS(I,M)+WW*(F6Q(I,L)*F6Q(M,L))*E*R
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
      DO 170 I=1,6
      AA=FL(M,L)
      QT(M,I)=QT(M,I)+WW*AA*(A*DF6QDS(I,L)+C*DF6QDT(I,L))*R+            
     1WW*(F6Q(I,L)*FL(M,L))*E*IAXSYM
      QT(M,I+6)=QT(M,I+6)+WW*AA*(B*DF6QDS(I,L)+D*DF6QDT(I,L))*R
  170 CONTINUE
      DO 190 I=1,3
      DO 190 N=1,12
      AA=0.
      DO 180 M=1,3
      AA=AA+MASSP(I,M)*QT(M,N)
  180 CONTINUE
      KPU(I,N)=AA
  190 CONTINUE
      DO 210 I=1,12
      DO 210 N=1,12
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
      DO 240 I=1,6
      DO 240 M=I,6
      IF (I.EQ.M) GO TO 240
      MASS(M,I)=MASS(I,M)
  240 CONTINUE
      DO 250 KL=1,KLEND
      DO 250 I=1,6
      DO 250 M=I,6
      IF (I.EQ.M) GO TO 250
      KUV(M,I,KL)=KUV(I,M,KL)
      KUV(M+6,I,KL)=KUV(I,M+6,KL)
      KUV(M,I+6,KL)=KUV(I+6,M,KL)
      KUV(M+6,I+6,KL)=KUV(I+6,M+6,KL)
      FF(M,I,KL)=FF(I,M,KL)
      FBDY(M,I,KL)=FBDY(I,M,KL)
  250 CONTINUE
C
C     RE-ARRANGE VISCOUS TERMS AND CONSTRUCT REMAINING DIFFUSION TERMS
C
      DO 260 KL=1,KLEND
      DO 260 I=1,6
      DO 260 M=1,6
      AA=KUV(I,M,KL)
      BB=KUV(I+6,M+6,KL)
      KUV(I,M,KL)=2.*AA+BB+FF(I,M,KL)*IAXSYM
      KUV(I+6,M+6,KL)=2.*BB+AA
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
