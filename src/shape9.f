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
      SUBROUTINE SHAPE9 (NGAUSS,GS1,GS2)
C
C     ******************************************************************
C
C     SUBROUTINE TO CALCULATE THE VALUES OF THE NINE NODE 'LAGRANGE'
C     SHAPE FUNCTIONS AND THEIR DERIVATIVES AT THE INTEGRATION POINTS
C     (N X N GAUSS INTEGRATION RULE)
C
C     ******************************************************************
C
      COMMON /QUAD9/ F9Q(9,9),F9L(4,9),DF9QDS(9,9),DF9QDT(9,9)
C
      DIMENSION GS1(*), GS2(*)
      DIMENSION SI(8), TI(8)
C
      DATA (SI(I),I=1,8)/-1.,1.,1.,-1.,0.,1.,0.,-1./
      DATA (TI(I),I=1,8)/-1.,-1.,1.,1.,-1.,0.,1.,0./
C
C     STATEMENT FUNCTIONS FOR THE SHAPE FUNCTIONS AND THEIR DERIVATIVES
C
      FC(I,S,T)=.25*(1.+S*SI(I))*(1.+T*TI(I))*(S*SI(I)+T*TI(I)-1.)
      FM(I,S,T)=.50*((1.-S**2)*(1.+T*TI(I))*TI(I)**2+(1.-T**2)*(1.+S*SI 
     1(I))*SI(I)**2)
      BUBBLE(S,T)=(1.-S**2)*(1.-T**2)
      G(I,S,T)=.25*(1.+S*SI(I))*(1.+T*TI(I))
      DFCDS(I,S,T)=.25*(1.+T*TI(I))*(2.*S*SI(I)+T*TI(I))*SI(I)
      DFMDS(I,S,T)=.50*((1.+T*TI(I))*(-2.*S)*TI(I)**2+(1.-T**2)*SI(I))
      DFCDT(I,S,T)=.25*(1.+S*SI(I))*(2.*T*TI(I)+S*SI(I))*TI(I)
      DFMDT(I,S,T)=.50*((1.+S*SI(I))*(-2.*T)*SI(I)**2+(1.-S**2)*TI(I))
      DBUBDS(S,T)=-2.*S*(1.-T**2)
      DBUBDT(S,T)=-2.*T*(1.-S**2)
C
C     ******************************************************************
C
C     EVALUATE FUNCTIONS
C
      DO 10 I=1,4
      L=1
      DO 10 K=1,NGAUSS
      DO 10 J=1,NGAUSS
      DF9QDS(I,L)=DFCDS(I,GS1(J),GS2(K))+.25*DBUBDS(GS1(J),GS2(K))
      DF9QDT(I,L)=DFCDT(I,GS1(J),GS2(K))+.25*DBUBDT(GS1(J),GS2(K))
      F9Q(I,L)=FC(I,GS1(J),GS2(K))+.25*BUBBLE(GS1(J),GS2(K))
      F9L(I,L)=G(I,GS1(J),GS2(K))
      L=L+1
   10 CONTINUE
      DO 20 I=5,8
      L=1
      DO 20 K=1,NGAUSS
      DO 20 J=1,NGAUSS
      DF9QDS(I,L)=DFMDS(I,GS1(J),GS2(K))-.5*DBUBDS(GS1(J),GS2(K))
      DF9QDT(I,L)=DFMDT(I,GS1(J),GS2(K))-.5*DBUBDT(GS1(J),GS2(K))
      F9Q(I,L)=FM(I,GS1(J),GS2(K))-.5*BUBBLE(GS1(J),GS2(K))
      L=L+1
   20 CONTINUE
      L=1
      DO 30 K=1,NGAUSS
      DO 30 J=1,NGAUSS
      DF9QDS(9,L)=DBUBDS(GS1(J),GS2(K))
      DF9QDT(9,L)=DBUBDT(GS1(J),GS2(K))
      F9Q(9,L)=BUBBLE(GS1(J),GS2(K))
      L=L+1
   30 CONTINUE
C
      RETURN
      END
