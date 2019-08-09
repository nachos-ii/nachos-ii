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
      SUBROUTINE PTFUNC (KIND,IVAR,S,T,VAR,V)
C
C     ******************************************************************
C
C     SUBROUTINE TO EVALUATE A DEPENDENT VARIABLE AT A GIVEN POINT (S,T)
C     WITHIN AN ELEMENT
C
C     ******************************************************************
C
      DIMENSION VAR(9)
C
C     ******************************************************************
C
      GO TO (10, 10, 20, 20, 30, 30), KIND
C
C     TRIANGULAR ELEMENTS (TRI6/3,TRI6/6)
C
   10 CONTINUE
      IF (IVAR.EQ.3) THEN
      V=S*VAR(1)+T*VAR(2)+(1.-S-T)*VAR(3)
       ELSE
      V=(2.*S*S-S)*VAR(1)+(2.*T*T-T)*VAR(2)+(2.*(S+T)**2-3.*(S+T)+1.)*  
     1VAR(3)+(4.*S*T)*VAR(4)+(-4.*T*T-4.*S*T+4.*T)*VAR(5)+(-4.*S*S-4.*S*
     2T+4.*S)*VAR(6)
      END IF
      RETURN
C
C     QUADRILATERAL ELEMENTS (QUAD8/4,QUAD8/8)
C
   20 CONTINUE
      IF (IVAR.EQ.3) THEN
      V=.25*(1.-S)*(1.-T)*VAR(1)+.25*(1.+S)*(1.-T)*VAR(2)+.25*(1.+S)*   
     1(1.+T)*VAR(3)+.25*(1.-S)*(1.+T)*VAR(4)
       ELSE
      V=.25*(1.-S)*(1.-T)*(-S-T-1.)*VAR(1)+.25*(1.+S)*(1.-T)*(S-T-1.)*  
     1VAR(2)+.25*(1.+S)*(1.+T)*(S+T-1.)*VAR(3)+.25*(1.-S)*(1.+T)*       
     2(-S+T-1.)*VAR(4)+.50*(1.-S*S)*(1.-T)*VAR(5)+.50*(1.+S)*(1.-T*T)*  
     3VAR(6)+.50*(1.-S*S)*(1.+T)*VAR(7)+.50*(1.-S)*(1.-T*T)*VAR(8)
      END IF
      RETURN
C
C     QUADRILATERAL ELEMENTS (QUAD9/4,QUAD9/9)
C
   30 CONTINUE
      IF (IVAR.EQ.3) THEN
      V=.25*(1.-S)*(1.-T)*VAR(1)+.25*(1.+S)*(1.-T)*VAR(2)+.25*(1.+S)*   
     1(1.+T)*VAR(3)+.25*(1.-S)*(1.+T)*VAR(4)
       ELSE
      V=.25*(1.-S)*(1.-T)*(S*T)*VAR(1)+.25*(1.+S)*(1.-T)*(-S*T)*VAR(2)+ 
     1.25*(1.+S)*(1.+T)*(S*T)*VAR(3)+.25*(1.-S)*(1.+T)*(-S*T)*VAR(4)+   
     2.50*(1.-S*S)*(1.-T)*(-T)*VAR(5)+.50*(1.+S)*(1.-T*T)*(S)*VAR(6)+   
     3.50*(1.-S*S)*(1.+T)*(T)*VAR(7)+.50*(1.-S)*(1.-T*T)*(-S)*VAR(8)+   
     4(1.-S*S)*(1.-T*T)*VAR(9)
      END IF
C
      RETURN
      END
