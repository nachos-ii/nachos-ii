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
      SUBROUTINE POINT (X,Y,ICON,LISTEL)
C
C     ******************************************************************
C
C     SUBROUTINE TO LOCATE SPECIAL OUTPUT POINTS IN THE MESH
C
C     ******************************************************************
C
      CHARACTER*10 CDATA
C
      COMMON /INDATR/ RDATA(150)
      COMMON /INDATI/ IDATA(150)
      COMMON /INDATC/ CDATA(150)
      COMMON /TAPES/  NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7, 
     1                NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
      COMMON /SZDAT/  NUMEL,NUMNOD,NUMVAR,NUMDOF,NODSOL
      COMMON /ELMDAT/ NNELM(6),NNCOR(6),NNMID(6),NNCTR(6),NSIDET(5),    
     1                NSIDEQ(7),NNSIDE(6,4,3)
      COMMON /SPTDAT/ ISPT,NSPT,NELSPT(50),PTS(50,2),SPTVAL(50,6)
C
      DIMENSION X(*), Y(*), ICON(NUMEL,*), LISTEL(*)
      DIMENSION XE(9), YE(9), XX(50), YY(50)
C
C     ******************************************************************
C
C     READ AND STORE INPUT DATA
C
   10 CONTINUE
      CALL RDFFLD (NIN,RDATA,IDATA,CDATA)
      IF (CDATA(1).EQ.'END       ') GO TO 120
C
      J=1
      DO 20 I=1,50
      XX(I)=RDATA(J)
      YY(I)=RDATA(J+1)
      J=J+2
      NPT=I
      A=ABS(XX(I)+YY(I))
      IF (A.LE.1.0E-5) GO TO 30
   20 CONTINUE
      NPT=50
      GO TO 40
   30 CONTINUE
      NPT=NPT-1
   40 CONTINUE
      WRITE (NOUT, 130) NPT,(I,XX(I),YY(I),I=1,NPT)
C
C     LOOP ON SPECIAL POINTS
C
      DO 110 I=1,NPT
      II=NSPT+I
      XB=XX(I)
      YB=YY(I)
C
C     LOOP ON ELEMENTS
C
      DO 90 N=1,NUMEL
      NN=ABS(LISTEL(N))
      KIND=MOD(NN,100)
      NN=NNELM(KIND)
      DO 50 J=1,NN
      INODE=ICON(N,J)
      XE(J)=X(INODE)
      YE(J)=Y(INODE)
   50 CONTINUE
      GO TO (60, 70, 60, 70, 60, 70), KIND
C
C     LOCATE ELEMENT FOR EACH SPECIAL POINT
C
   60 CONTINUE
      NN=NNCOR(KIND)
      CALL SRCHL (NN,XE,YE,XB,YB,S,T,IFIND)
      GO TO 80
   70 CONTINUE
      NN=NNELM(KIND)
      CALL SRCHQ (NN,XE,YE,XB,YB,S,T,IFIND)
   80 CONTINUE
      IF (IFIND.EQ.0) GO TO 90
      NELSPT(II)=N
      PTS(II,1)=S
      PTS(II,2)=T
      GO TO 100
   90 CONTINUE
      NELSPT(II)=0
      WRITE (NOUT, 160) I,XB,YB
  100 CONTINUE
  110 CONTINUE
      I1=NSPT+1
      I2=NSPT+NPT
      WRITE (NOUT, 140)
      WRITE (NOUT, 150) (I,NELSPT(I),PTS(I,1),PTS(I,2),I=I1,I2)
      WRITE (NOUT, 170)
      NSPT=NSPT+NPT
      GO TO 10
C
  120 CONTINUE
      RETURN
C
  130 FORMAT (///,3X,'DEPENDENT VARIABLE DATA (I.E. VELOCITY,PRESSURE,TE
     1MPERATURE,ETC)' ,/,3X,'FOR THE FOLLOWING' ,I3,' POINTS HAVE BEEN R
     2EQUESTED--' ,//,(10X,'PT',I5,2X,'X= 'E15.7,2X,'Y= 'E15.7,:,5X,'PT'
     3    ,I5,2X,'X= 'E15.7,2X,'Y= 'E15.7))
  140 FORMAT (//,3X,'SPECIAL OUTPUT POINTS WERE LOCATED IN MESH AS FOLLO
     1WS--' ,/)
  150 FORMAT (10X,'PT'I5,' FOUND IN ELEMENT' I5,'  S='E15.7,2X,'  T='E15
     1.7)
  160 FORMAT (//,10X,'--SPECIAL OUTPUT POINT' I5,' XB= 'E15.7,' YB= 'E15
     1.7,' WAS NOT LOCATED IN MESH' )
  170 FORMAT (///,2X,'**************************************************
     1******************************************************************
     2**************'   ,//)
      END
