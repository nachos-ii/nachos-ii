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
      SUBROUTINE EXTRCT (X,Y,MDF,NOPP,LISTEL,LISTND,ICON,SK,SCRTCH)
C
C     ******************************************************************
C
C     SUBROUTINE TO EXTRACT THE SOLUTION FROM THE GLOBAL DOF VECTOR
C     USED IN THE FRONTAL PROCEDURE
C     ALSO GENERATES PRESSURE DATA AT ELEMENT MIDSIDE NODES
C
C     ******************************************************************
C
      REAL KPU
C
      COMMON /TAPES/  NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7, 
     1                NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
      COMMON /SZDAT/  NUMEL,NUMNOD,NUMVAR,NUMDOF,NODSOL
      COMMON /MATDAT/ NMAT,NPROP,PROP(15,10),NXMAT,NXPROP,XPROP(15,10)
      COMMON /PRBDAT/ IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,    
     1                IPNLTY,PNLTY
      COMMON /ELMDAT/ NNELM(6),NNCOR(6),NNMID(6),NNCTR(6),NSIDET(5),    
     1                NSIDEQ(7),NNSIDE(6,4,3)
      COMMON /BCDAT/  BCSET(2,20,2),IPELM,IPNODE,PBC
C
      DIMENSION X(*), Y(*)
      DIMENSION MDF(*), NOPP(*)
      DIMENSION ICON(NUMEL,*), LISTEL(*), LISTND(*)
      DIMENSION SCRTCH(NUMNOD,*), SK(*)
      DIMENSION KPU(3,18)
      DIMENSION U(9), V(9), P(3), IPLACE(5)
C
      DATA (IPLACE(I),I=1,5)/1,2,4,5,6/
C
C     ******************************************************************
C
      DO 10 I=1,NUMNOD
      SCRTCH(I,3)=0.
   10 CONTINUE
C
C     LOOP ON ELEMENTS
C
      DO 130 N=1,NUMEL
      NN=ABS(LISTEL(N))
      MAT=NN/100
      ISOLID=0
      IF (PROP(15,MAT).EQ.3.0) ISOLID=1
      KIND=MOD(NN,100)
      NCOR=NNCOR(KIND)
      NMID=NNMID(KIND)
      NCNTR=NNCTR(KIND)
      NN=NNELM(KIND)
C
C     LOOP ON CORNER NODES
C
      DO 20 I=1,NCOR
      INODE=ABS(ICON(N,I))
      KDOF=NOPP(INODE)
      NLDOF=MDF(INODE)
      DO 20 J=1,NLDOF
      JJ=J
      IF (IPFUNC.EQ.1.OR.IPNLTY.EQ.1) JJ=IPLACE(J)
      SCRTCH(INODE,JJ)=SK(KDOF+J-1)
   20 CONTINUE
C
C     LOOP ON MIDSIDE NODES
C
      DO 30 I=1,NMID
      INODE=ABS(ICON(N,I+NCOR))
      KDOF=NOPP(INODE)
      NLDOF=MDF(INODE)
      DO 30 J=1,NLDOF
      JJ=IPLACE(J)
      SCRTCH(INODE,JJ)=SK(KDOF+J-1)
   30 CONTINUE
C
C     CENTER NODE
C
      IF (NCNTR.EQ.1) THEN
      INODE=ABS(ICON(N,NN))
      KDOF=NOPP(INODE)
      NLDOF=MDF(INODE)
      DO 40 J=1,NLDOF
      JJ=IPLACE(J)
      SCRTCH(INODE,JJ)=SK(KDOF+J-1)
   40 CONTINUE
      END IF
C
C     DISCONTINUOUS PRESSURE FIELD
C
      IF (IPFUNC.EQ.1) THEN
      INODE=ABS(ICON(N,NN+1))
      KDOF=NOPP(INODE)
      NLDOF=MDF(INODE)
      DO 50 J=1,NLDOF
      P(J)=SK(KDOF+J-1)
   50 CONTINUE
      IF (ISOLID.EQ.1) GO TO 70
      DO 60 I=1,NCOR
      INODE=ABS(ICON(N,I))
      PNODE=P(1)+P(2)*X(INODE)+P(3)*Y(INODE)
      SCRTCH(INODE,3)=SCRTCH(INODE,3)+PNODE/LISTND(INODE)
   60 CONTINUE
   70 CONTINUE
      END IF
C
C     PENALTY PRESSURE FIELD
C
      IF (IPNLTY.EQ.1) THEN
      READ (NTP3) NSZ1,NSZ2,KPU
      DO 80 I=1,NN
      INODE=ABS(ICON(N,I))
      U(I)=SCRTCH(INODE,1)
      V(I)=SCRTCH(INODE,2)
   80 CONTINUE
      DO 100 I=1,3
      SUMX=0.
      SUMY=0.
      DO 90 J=1,NN
      SUMX=SUMX+KPU(I,J)*U(J)
      SUMY=SUMY+KPU(I,J+NN)*V(J)
   90 CONTINUE
      P(I)=-(SUMX+SUMY)*(1.0/PNLTY)
  100 CONTINUE
      IF (ISOLID.EQ.1) GO TO 120
      DO 110 I=1,NCOR
      INODE=ABS(ICON(N,I))
      PNODE=P(1)+P(2)*X(INODE)+P(3)*Y(INODE)
      SCRTCH(INODE,3)=SCRTCH(INODE,3)+PNODE/LISTND(INODE)
  110 CONTINUE
  120 CONTINUE
      END IF
  130 CONTINUE
C
C     GENERATE PRESSURE DATA AT MIDSIDE NODES
C
      CALL MIDSID (SCRTCH(1,3),ICON,LISTEL)
C
      IF (IPFUNC.EQ.1.OR.IPNLTY.EQ.1) THEN
      IF (IPELM.EQ.0) GO TO 150
      INODE=ABS(ICON(IPELM,IPNODE))
      PNODE=SCRTCH(INODE,3)
      PDIFF=PNODE-PBC
      DO 140 I=1,NUMNOD
      SCRTCH(I,3)=SCRTCH(I,3)-PDIFF
  140 CONTINUE
      END IF
  150 CONTINUE
C
      RETURN
      END
