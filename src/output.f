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
      SUBROUTINE OUTPUT (ITYPE,LISTEL,ICON,UN)
C
C     ******************************************************************
C
C     SUBROUTINE TO OUTPUT SOLUTION FIELDS
C
C     ******************************************************************
C
      COMMON /TAPES/  NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7, 
     1                NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
      COMMON /CONTRL/ IEDIT,IPRINT
      COMMON /SZDAT/  NUMEL,NUMNOD,NUMVAR,NUMDOF,NODSOL
      COMMON /MATDAT/ NMAT,NPROP,PROP(15,10),NXMAT,NXPROP,XPROP(15,10)
      COMMON /PRBDAT/ IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,    
     1                IPNLTY,PNLTY
      COMMON /ELMDAT/ NNELM(6),NNCOR(6),NNMID(6),NNCTR(6),NSIDET(5),    
     1                NSIDEQ(7),NNSIDE(6,4,3)
      COMMON /SPTDAT/ ISPT,NSPT,NELSPT(50),PTS(50,2),SPTVAL(50,6)
C
      DIMENSION LISTEL(*), ICON(NUMEL,*), UN(NUMNOD,*)
      DIMENSION UNELM(9,6)
C
C     ******************************************************************
C
C     PRINT SOLUTION FIELDS
C
      WRITE (NOUT, 50)
      WRITE (NOUT, 60)
      IF (ITYPE.GE.2) WRITE (NOUT, 70)
      IF (ITYPE.GE.2.AND.(IVAR1+IVAR2).GT.0) WRITE (NOUT, 80)
      WRITE (NOUT, 90)
      DO 40 N=1,NUMEL
      IF (IEDIT.EQ.0) GO TO 10
      IF (LISTEL(N).GT.0) GO TO 40
   10 CONTINUE
      NN=ABS(LISTEL(N))
      KIND=MOD(NN,100)
      MAT=NN/100
      NN=NNELM(KIND)
      NCOR=NNCOR(KIND)
      DO 20 I=1,NN
      INODE=ABS(ICON(N,I))
      DO 20 J=1,NUMVAR
      UNELM(I,J)=UN(INODE,J)
   20 CONTINUE
      IF (PROP(15,MAT).EQ.3.0) GO TO 30
      WRITE (NOUT, 100) N,(UNELM(I,1),I=1,NN)
      WRITE (NOUT, 110) (UNELM(I,2),I=1,NN)
      WRITE (NOUT, 120) (UNELM(I,3),I=1,NCOR)
   30 CONTINUE
      IF (ITYPE.GE.2) WRITE (NOUT, 130) N,(UNELM(I,4),I=1,NN)
      IF (ITYPE.GE.2.AND.IVAR1.EQ.1) WRITE (NOUT, 160) (UNELM(I,5),I=1, 
     1NN)
      IF (ITYPE.GE.2.AND.IVAR2.EQ.1) WRITE (NOUT, 170) (UNELM(I,6),I=1, 
     1NN)
   40 CONTINUE
C
C     PRINT SOLUTION FOR SPECIAL POINTS
C
      IF (ISPT.EQ.0) RETURN
      WRITE (NOUT, 140)
      WRITE (NOUT, 150) (I,(SPTVAL(I,J),J=1,6),I=1,NSPT)
      RETURN
C
   50 FORMAT (////,2X,'*************************************************
     1******************************************************************
     2***************',//,3X,'CURRENT VALUES FOR :',/)
   60 FORMAT (10X,'VELOCITY AND PRESSURE FIELDS' )
   70 FORMAT (10X,'TEMPERATURE FIELD' )
   80 FORMAT (10X,'AUXILIARY VARIABLE FIELDS')
   90 FORMAT (////,67X,'--NODES--',//,2X,'ELEMENTS',7X,'-1-',12X,'-2-',1
     12X,'-3-',12X,'-4-',12X,'-5-',12X,'-6-',12X,'-7-',12X,'-8-',/)
  100 FORMAT (1X,I4,1X,'U',2X,8E15.7,:,/,9X,E15.7)
  110 FORMAT (6X,'V',2X,8E15.7,:,/,9X,E15.7)
  120 FORMAT (6X,'P',2X,4E15.7)
  130 FORMAT (1X,I4,1X,'T',2X,8E15.7,:,/,9X,E15.7)
  140 FORMAT (//,47X,'--SOLUTION AT REQUESTED OUTPUT POINTS--' ,/)
  150 FORMAT (10X,I5,2X,' U = ',E15.7,2X,' V = ',E15.7,2X,' P = ',E15.7,
     1/,17X,' T = ',E15.7,2X,' V1 = ',E15.7,2X,' V2 = ',E15.7)
  160 FORMAT (5X,'V1',2X,8E15.7,:,/,9X,E15.7)
  170 FORMAT (5X,'V2',2X,8E15.7,:,/,9X,E15.7)
      END
