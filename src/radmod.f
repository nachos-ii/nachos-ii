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
      SUBROUTINE RADMOD (KIND,IVALUE,ISIDE,KTR,FR,T,NN,TIME,X,Y,KNUM)
C
C     ******************************************************************
C
C     SUBROUTINE TO CONSTRUCT THE MATRICES FOR A RADIATION
C     BOUNDARY CONDITION
C
C     ******************************************************************
C
      REAL KTR
C
      COMMON /ELMDAT/ NNELM(6),NNCOR(6),NNMID(6),NNCTR(6),NSIDET(5),    
     1                NSIDEQ(7),NNSIDE(6,4,3)
      COMMON /BCDAT/  BCSET(2,20,2),IPELM,IPNODE,PBC
C
      DIMENSION KTR(9,3,9), FR(9,3), H(3), II(3)
      DIMENSION T(9), X(9), Y(9)
C
C     ******************************************************************
C
C     COLLECT DATA FOR APPROPRIATE SIDE OF ELEMENT
C
      DO 10 I=1,3
      II(I)=NNSIDE(KIND,ISIDE,I)
   10 CONTINUE
C
C     COMPUTE EFFECTIVE HEAT TRANSFER COEFFICIENT
C
      DO 20 I=1,3
      IPT=II(I)
      IF (BCSET(2,IVALUE,1).LT.1.0E20) THEN
      H(I)=BCSET(2,IVALUE,1)*(BCSET(2,IVALUE,2)+T(IPT))*                
     1(BCSET(2,IVALUE,2)**2+T(IPT)**2)
       ELSE
      TSURF=T(IPT)
      XSURF=X(IPT)
      YSURF=Y(IPT)
      TREF=BCSET(2,IVALUE,2)
      CALL USRHTC (HT,TSURF,TREF,XSURF,YSURF,TIME,KNUM,IVALUE)
      BCSET(2,IVALUE,2)=TREF
      H(I)=HT
      END IF
   20 CONTINUE
C
C     CONSTRUCT RADIATION MATRIX AND VECTOR
C
      DO 40 I=1,NN
      DO 40 J=1,NN
      FD=0.
      DO 30 K=1,3
   30 FD=FD+KTR(I,K,J)*H(K)
      KTR(I,1,J)=FD
   40 CONTINUE
      DO 60 I=1,NN
      FD=0.
      DO 50 K=1,3
   50 FD=FD+FR(I,K)*H(K)
      FR(I,1)=FD*BCSET(2,IVALUE,2)
   60 CONTINUE
C
      RETURN
      END
