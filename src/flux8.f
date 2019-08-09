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
      SUBROUTINE FLUX8 (NSTRSQ,SQ,TQ)
C
C     ******************************************************************
C
C     SUBROUTINE TO CALCULATE THE VALUES OF THE EIGHT NODE 'SERENDIPITY'
C     SHAPE FUNCTIONS AND THEIR DERIVATIVES AT THE FLUX EVALUATION
C     POINTS.
C     THE POINTS MAY BE ON THE ELEMENT BOUNDARY OR AT INTERIOR POINTS
C
C     ******************************************************************
C
      COMMON /QUAD8/ F8Q(8,9),F8L(4,9),DF8QDS(8,9),DF8QDT(8,9)
C
      DIMENSION SQ(*), TQ(*)
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
      G(I,S,T)=.25*(1.+S*SI(I))*(1.+T*TI(I))
      DFCDS(I,S,T)=.25*(1.+T*TI(I))*(2.*S*SI(I)+T*TI(I))*SI(I)
      DFMDS(I,S,T)=.50*((1.+T*TI(I))*(-2.*S)*TI(I)**2+(1.-T**2)*SI(I))
      DFCDT(I,S,T)=.25*(1.+S*SI(I))*(2.*T*TI(I)+S*SI(I))*TI(I)
      DFMDT(I,S,T)=.50*((1.+S*SI(I))*(-2.*T)*SI(I)**2+(1.-S**2)*TI(I))
C
C     ******************************************************************
C
C     EVALUATE CORNER FUNCTIONS
C
      DO 10 I=1,4
C
      DO 10 J=1,NSTRSQ
      DF8QDS(I,J)=DFCDS(I,SQ(J),TQ(J))
      DF8QDT(I,J)=DFCDT(I,SQ(J),TQ(J))
      F8Q(I,J)=FC(I,SQ(J),TQ(J))
      F8L(I,J)=G(I,SQ(J),TQ(J))
   10 CONTINUE
C
C     EVALUATE MID-SIDE FUNCTIONS
C
      DO 20 I=5,8
C
      DO 20 J=1,NSTRSQ
      DF8QDS(I,J)=DFMDS(I,SQ(J),TQ(J))
      DF8QDT(I,J)=DFMDT(I,SQ(J),TQ(J))
      F8Q(I,J)=FM(I,SQ(J),TQ(J))
   20 CONTINUE
C
      RETURN
      END
