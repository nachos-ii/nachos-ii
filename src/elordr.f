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
      SUBROUTINE ELORDR (KNUM,ELSTIF,ELFV,LISTEL)
C
C     ******************************************************************
C
C     SUBROUTINE TO  RE-ORDER ELEMENT EQUATIONS TO AGREE WITH
C     ORDERING IN FRONTAL PROGRAM
C
C     ******************************************************************
C
      COMMON /PRBDAT/ IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,    
     1                IPNLTY,PNLTY
      COMMON /ELMDAT/ NNELM(6),NNCOR(6),NNMID(6),NNCTR(6),NSIDET(5),    
     1                NSIDEQ(7),NNSIDE(6,4,3)
C
      DIMENSION LISTEL(*)
      DIMENSION ELSTIF(50,50), ELFV(50), STIF(50,50), F(50)
      DIMENSION LIST(50)
C
C     ******************************************************************
C
C     CONSTRUCT POINTER VECTOR FOR LOCATION OF DIFFERENT VARIABLES
C
      NN=ABS(LISTEL(KNUM))
      KIND=MOD(NN,100)
      NCOR=NNCOR(KIND)
      NMID=NNMID(KIND)
      NCNTR=NNCTR(KIND)
      NN=NNELM(KIND)
      NVEL=2
      NP=1
      NTEMP=IFORCE+IFREE
      NVAR=IVAR1+IVAR2
      NPRESS=1
      IF (IPFUNC.EQ.1.OR.IPNLTY.EQ.1) NPRESS=0
      DO 10 I=1,50
      LIST(I)=0
   10 CONTINUE
C
C     CORNER NODES
C
      NSKIP=NVEL+NPRESS+NTEMP+NVAR
      DO 20 I=1,NCOR
      II=1+(I-1)*NSKIP
      LIST(II)=I
      II=II+1
      LIST(II)=I+NN
   20 CONTINUE
      IF (NPRESS.EQ.1) THEN
      DO 30 I=1,NCOR
      II=NVEL+1+(I-1)*NSKIP
      LIST(II)=I+NVEL*NN
   30 CONTINUE
      END IF
      IF (NTEMP.EQ.1) THEN
      DO 40 I=1,NCOR
      II=NVEL+NPRESS+1+(I-1)*NSKIP
      LIST(II)=I+NVEL*NN+NP*NCOR
   40 CONTINUE
      END IF
      IF (IVAR1.EQ.1) THEN
      DO 50 I=1,NCOR
      II=NVEL+NPRESS+NTEMP+1+(I-1)*NSKIP
      LIST(II)=I+NVEL*NN+NP*NCOR+NTEMP*NN
   50 CONTINUE
      END IF
      IF (IVAR2.EQ.1) THEN
      DO 60 I=1,NCOR
      II=NVEL+NPRESS+NTEMP+IVAR1+1+(I-1)*NSKIP
      LIST(II)=I+NVEL*NN+NP*NCOR+NTEMP*NN+IVAR1*NN
   60 CONTINUE
      END IF
C
C     MIDSIDE NODES
C
      III=NSKIP*NCOR
      NSKIP=NVEL+NTEMP+NVAR
      DO 70 I=1,NMID
      II=1+III+(I-1)*NSKIP
      LIST(II)=I+NCOR
      II=II+1
      LIST(II)=I+NCOR+NN
   70 CONTINUE
      IF (NTEMP.EQ.1) THEN
      DO 80 I=1,NMID
      II=NVEL+1+III+(I-1)*NSKIP
      LIST(II)=I+NVEL*NN+NP*NCOR+NCOR
   80 CONTINUE
      END IF
      IF (IVAR1.EQ.1) THEN
      DO 90 I=1,NMID
      II=NVEL+NTEMP+1+III+(I-1)*NSKIP
      LIST(II)=I+NVEL*NN+NP*NCOR+NTEMP*NN+NCOR
   90 CONTINUE
      END IF
      IF (IVAR2.EQ.1) THEN
      DO 100 I=1,NMID
      II=NVEL+NTEMP+IVAR1+1+III+(I-1)*NSKIP
      LIST(II)=I+NVEL*NN+NP*NCOR+NTEMP*NN+IVAR1*NN+NCOR
  100 CONTINUE
      END IF
C
C     CENTER NODE
C
      III=III+NSKIP*NMID
      IF (NCNTR.EQ.0) GO TO 110
      II=III+1
      LIST(II)=NN
      II=II+1
      LIST(II)=NVEL*NN
      IF (NTEMP.EQ.1) THEN
      II=II+1
      LIST(II)=NVEL*NN+NP*NCOR+NN
      END IF
      IF (IVAR1.EQ.1) THEN
      II=II+1
      LIST(II)=NVEL*NN+NP*NCOR+NTEMP*NN+NN
      END IF
      IF (IVAR2.EQ.1) THEN
      II=II+1
      LIST(II)=NVEL*NN+NP*NCOR+NTEMP*NN+IVAR1*NN+NN
      END IF
      III=III+NVEL+NTEMP+IVAR1+IVAR2
C
C     DISCONTINUOUS PRESSURE
C
  110 CONTINUE
      IF (IPFUNC.EQ.1) THEN
      II=III+1
      LIST(II)=NVEL*NN+1
      II=II+1
      LIST(II)=NVEL*NN+2
      II=II+1
      LIST(II)=NVEL*NN+3
      END IF
C
      KDOF=NVEL*NN+NTEMP*NN+NVAR*NN
      IF (IPFUNC.EQ.0.AND.IPNLTY.EQ.0) KDOF=KDOF+NP*NCOR
      IF (IPFUNC.EQ.1) KDOF=KDOF+3
C
C     COPY MATRIX AND VECTOR TO WORK SPACE
C
      DO 120 I=1,50
      F(I)=ELFV(I)
      DO 120 J=1,50
      STIF(I,J)=ELSTIF(I,J)
  120 CONTINUE
C
C     RE-ORDER MATRIX AND RHS VECTOR
C
      DO 140 I=1,KDOF
      II=LIST(I)
      DO 130 J=1,KDOF
      JJ=LIST(J)
      ELSTIF(I,J)=STIF(II,JJ)
  130 CONTINUE
      ELFV(I)=F(II)
  140 CONTINUE
C
      RETURN
      END
