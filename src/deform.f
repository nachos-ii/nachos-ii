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
      SUBROUTINE DEFORM (KIND,X,Y,U,V,SCNDII)
C
C     ******************************************************************
C
C     SUBROUTINE TO EVALUATE THE INVARIANTS OF THE RATE OF DEFORMATION
C     TENSOR AT THE ELEMENT INTEGRATION POINTS
C
C     ******************************************************************
C
      COMMON /PRBDAT/ IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,    
     1                IPNLTY ,PNLTY
      COMMON /GAUSS1/ GSPT1(4),GSPT2(4)
      COMMON /HAMMR1/ HAMPT1(3),HAMPT2(3)
      COMMON /TRI6/   F6Q(6,7),F6L(3,7),DF6QDS(6,7),DF6QDT(6,7)
      COMMON /QUAD8/  F8Q(8,9),F8L(4,9),DF8QDS(8,9),DF8QDT(8,9)
      COMMON /QUAD9/  F9Q(9,9),F9L(4,9),DF9QDS(9,9),DF9QDT(9,9)
C
      DIMENSION X(9), Y(9), U(9), V(9), SCNDII(4)
C
      DATA NGAUSS/4/,NHAM/3/
C
C     STATEMENT FUNCTIONS FOR THE JACOBIAN OF THE PARAMETRIC MAPPING
C     AND OTHER ELEMENT FUNCTIONS
C
      F11Q(S,T)=.25*(D1+S*D2+T*D3+S*T*D4+S**2*D5+S**2*T*D6)
      F12Q(S,T)=-.25*(B1+T*B2+S*B3+S*T*B4+T**2*B5+S*T**2*C6)
      F21Q(S,T)=-.25*(C1+S*C2+T*C3+S*T*C4+S**2*C5+S**2*T*B6)
      F22Q(S,T)=.25*(A1+T*A2+S*A3+S*T*A4+T**2*A5+S*T**2*A6)
      DETJQ(S,T)=F11Q(S,T)*F22Q(S,T)-F12Q(S,T)*F21Q(S,T)
C
      F11T(S,T)=S*A1+T*A2+A3
      F12T(S,T)=-(S*C1+T*C2+C3)
      F21T(S,T)=-(S*B1+T*B2+B3)
      F22T(S,T)=S*D1+T*D2+D3
      DETJT(S,T)=F11T(S,T)*F22T(S,T)-F12T(S,T)*F21T(S,T)
C
      U6Q(L)=F6Q(1,L)*U(1)+F6Q(2,L)*U(2)+F6Q(3,L)*U(3)+F6Q(4,L)*U(4)+   
     1F6Q(5,L)*U(5)+F6Q(6,L)*U(6)
      U8Q(L)=F8Q(1,L)*U(1)+F8Q(2,L)*U(2)+F8Q(3,L)*U(3)+F8Q(4,L)*U(4)+   
     1F8Q(5,L)*U(5)+F8Q(6,L)*U(6)+F8Q(7,L)*U(7)+F8Q(8,L)*U(8)
      U9Q(L)=F9Q(1,L)*U(1)+F9Q(2,L)*U(2)+F9Q(3,L)*U(3)+F9Q(4,L)*U(4)+   
     1F9Q(5,L)*U(5)+F9Q(6,L)*U(6)+F9Q(7,L)*U(7)+F9Q(8,L)*U(8)+          
     2F9Q(9,L)*U(9)
C
      RAD6L(L)=F6L(1,L)*X(1)+F6L(2,L)*X(2)+F6L(3,L)*X(3)
      RAD8L(L)=F8L(1,L)*X(1)+F8L(2,L)*X(2)+F8L(3,L)*X(3)+F8L(4,L)*X(4)
      RAD9L(L)=F9L(1,L)*X(1)+F9L(2,L)*X(2)+F9L(3,L)*X(3)+F9L(4,L)*X(4)
C
      RAD6Q(L)=F6Q(1,L)*X(1)+F6Q(2,L)*X(2)+F6Q(3,L)*X(3)+F6Q(4,L)*X(4)+ 
     1F6Q(5,L)*X(5)+F6Q(6,L)*X(6)
      RAD8Q(L)=F8Q(1,L)*X(1)+F8Q(2,L)*X(2)+F8Q(3,L)*X(3)+F8Q(4,L)*X(4)+ 
     1F8Q(5,L)*X(5)+F8Q(6,L)*X(6)+F8Q(7,L)*X(7)+F8Q(8,L)*X(8)
      RAD9Q(L)=F9Q(1,L)*X(1)+F9Q(2,L)*X(2)+F9Q(3,L)*X(3)+F9Q(4,L)*X(4)+ 
     1F9Q(5,L)*X(5)+F9Q(6,L)*X(6)+F9Q(7,L)*X(7)+F9Q(8,L)*X(8)+          
     2F9Q(9,L)*X(9)
C
C     ******************************************************************
C
C     SELECT ELEMENT TYPE
C
      GO TO (110, 120, 10, 20, 10, 30), KIND
C
C     QUADRILATERAL ELEMENTS
C
   10 CONTINUE
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
      GO TO 40
   20 CONTINUE
      A1=2.*(X(6)-X(8))
      A2=X(1)-X(2)+X(3)-X(4)
      A3=2.*(X(1)+X(2)+X(3)+X(4)-2.*X(5)-2.*X(7))
      A4=2.*(-X(1)-X(2)+X(3)+X(4)+2.*X(5)-2.*X(7))
      A5=-X(1)+X(2)+X(3)-X(4)-2.*X(6)+2.*X(8)
      A6=0.
      B1=2.*(Y(6)-Y(8))
      B2=Y(1)-Y(2)+Y(3)-Y(4)
      B3=2.*(Y(1)+Y(2)+Y(3)+Y(4)-2.*Y(5)-2.*Y(7))
      B4=2.*(-Y(1)-Y(2)+Y(3)+Y(4)+2.*Y(5)-2.*Y(7))
      B5=-Y(1)+Y(2)+Y(3)-Y(4)-2.*Y(6)+2.*Y(8)
      B6=0.
      C1=2.*(X(7)-X(5))
      C2=X(1)-X(2)+X(3)-X(4)
      C3=2.*(X(1)+X(2)+X(3)+X(4)-2.*X(6)-2.*X(8))
      C4=2.*(-X(1)+X(2)+X(3)-X(4)-2.*X(6)+2.*X(8))
      C5=-X(1)-X(2)+X(3)+X(4)+2.*X(5)-2.*X(7)
      C6=0.
      D1=2.*(Y(7)-Y(5))
      D2=Y(1)-Y(2)+Y(3)-Y(4)
      D3=2.*(Y(1)+Y(2)+Y(3)+Y(4)-2.*Y(6)-2.*Y(8))
      D4=2.*(-Y(1)+Y(2)+Y(3)-Y(4)-2.*Y(6)+2.*Y(8))
      D5=-Y(1)-Y(2)+Y(3)+Y(4)+2.*Y(5)-2.*Y(7)
      D6=0.
      GO TO 40
   30 CONTINUE
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
   40 CONTINUE
C
C     COMPUTE STRAIN RATES AT GAUSS POINTS
C
      DO 100 I=1,NGAUSS
      EXX=0.
      EYY=0.
      EXY=0.
      A=F11Q(GSPT1(I),GSPT2(I))
      B=F21Q(GSPT1(I),GSPT2(I))
      C=F12Q(GSPT1(I),GSPT2(I))
      D=F22Q(GSPT1(I),GSPT2(I))
      EI=1./DETJQ(GSPT1(I),GSPT2(I))
      GO TO (50, 50, 50, 50, 70, 70), KIND
   50 CONTINUE
      RG=RAD8L(I)
      IF (KIND.EQ.4) RG=RAD8Q(I)
      UG=U8Q(I)
      DO 60 J=1,8
      AA=DF8QDS(J,I)*A+DF8QDT(J,I)*C
      BB=DF8QDS(J,I)*B+DF8QDT(J,I)*D
      EXX=EXX+AA*U(J)*EI
      EYY=EYY+BB*V(J)*EI
      EXY=EXY+AA*V(J)*EI+BB*U(J)*EI
   60 CONTINUE
      GO TO 90
   70 CONTINUE
      RG=RAD9L(I)
      IF (KIND.EQ.6) RG=RAD9Q(I)
      UG=U9Q(I)
      DO 80 J=1,9
      AA=DF9QDS(J,I)*A+DF9QDT(J,I)*C
      BB=DF9QDS(J,I)*B+DF9QDT(J,I)*D
      EXX=EXX+AA*U(J)*EI
      EYY=EYY+BB*V(J)*EI
      EXY=EXY+AA*V(J)*EI+BB*U(J)*EI
   80 CONTINUE
   90 CONTINUE
      SCNDII(I)=2.*(EXX**2+EYY**2)+EXY**2+IAXSYM*2.*(UG/RG)**2
  100 CONTINUE
C
C     EXTRAPOLATE GAUSS POINT VALUES TO NODAL POINTS
C
      CALL EXTRAP (KIND,SCNDII)
      RETURN
C
C     TRIANGULAR ELEMENTS
C
  110 CONTINUE
      A1=0.
      A2=0.
      A3=Y(2)-Y(3)
      B1=0.
      B2=0.
      B3=-X(2)+X(3)
      C1=0.
      C2=0.
      C3=-Y(1)+Y(3)
      D1=0.
      D2=0.
      D3=X(1)-X(3)
      GO TO 130
  120 CONTINUE
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
  130 CONTINUE
C
C     COMPUTE STRAIN RATES AT INTEGRATION POINTS
C
      DO 150 I=1,NHAM
      EXX=0.
      EYY=0.
      EXY=0.
      A=F11T(HAMPT1(I),HAMPT2(I))
      B=F12T(HAMPT1(I),HAMPT2(I))
      C=F21T(HAMPT1(I),HAMPT2(I))
      D=F22T(HAMPT1(I),HAMPT2(I))
      EI=1./DETJT(HAMPT1(I),HAMPT2(I))
      RG=RAD6L(I)
      IF (KIND.EQ.2) RG=RAD6Q(I)
      UG=U6Q(I)
      DO 140 J=1,6
      AA=DF6QDS(J,I)*A+DF6QDT(J,I)*B
      BB=DF6QDS(J,I)*C+DF6QDT(J,I)*D
      EXX=EXX+AA*U(J)*EI
      EYY=EYY+BB*V(J)*EI
      EXY=EXY+AA*V(J)*EI+BB*U(J)*EI
  140 CONTINUE
      SCNDII(I)=2.*(EXX**2+EYY**2)+EXY**2+IAXSYM*2.*(UG/RG)**2
  150 CONTINUE
C
C     EXTRAPOLATE INTEGRATION POINT VALUES TO NODAL POINTS
C
      CALL EXTRAP (KIND,SCNDII)
C
      RETURN
      END
