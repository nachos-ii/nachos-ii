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
      SUBROUTINE ELEMNT (X,Y,IJK,IBC,BCVAL,NI,NJ,ORDER,IPRINT)
C
C     ******************************************************************
C
C     SUBROUTINE TO READ, STORE AND ORGANIZE ELEMENT AND BOUNDARY
C     CONDITION DATA
C
C     ******************************************************************
C
      CHARACTER*10 CDATA,CODES,BCCODE,TYPE1,TYPE2,TYPE,ORDER
C
      COMMON /INDATR/ RDATA(150)
      COMMON /INDATI/ IDATA(150)
      COMMON /INDATC/ CDATA(150)
      COMMON /TAPES/  NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7, 
     1                NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
      COMMON /SZDAT/  NUMEL,NUMNOD,NUMVAR,NUMDOF,NODSOL
      COMMON /LMTDAT/ MAXELM,MAXNOD,MAXBC,MAXEDG,MAXMAT,MAXNPT,MAXPLT,  
     1                MAXBRY,MAXBLK
      COMMON /MATDAT/ NMAT,NPROP,PROP(15,10),NXMAT,NXPROP,XPROP(15,10)
      COMMON /ELMDAT/ NNELM(6),NNCOR(6),NNMID(6),NNCTR(6),NSIDET(5),    
     1                NSIDEQ(7),NNSIDE(6,4,3)
      COMMON /BCDAT/  BCSET(2,20,2),IPELM,IPNODE,PBC
C
      DIMENSION X(NI,*), Y(NI,*), IBC(*), BCVAL(*), IJK(10,*)
      DIMENSION CODES(15), BCCODE(30), TYPE(30)
      DIMENSION NBC(30), NNBC(30), VAL(30)
      DIMENSION XX(9), YY(9), II(9), JJ(9)
C
      DATA (CODES(I),I=1,15)/'TRI6/3','TRI6/6','QUAD8/4','QUAD8/8',     
     1     'QUAD9/4','QUAD9/9',' ',' ','BC','ILOOP','JLOOP','IEND',     
     2     'JEND','SET','END'/
      DATA (BCCODE(I),I=1,30)/'U','V','P','T','V1','V2','STICK','USIDE',
     1     'VSIDE','TSIDE','V1SIDE','V2SIDE','UVARY','VVARY','TVARY',   
     2     'V1VARY','V2VARY','TNRMLSIDE','QSIDE','FLUX1SIDE',           
     3     'FLUX2SIDE','TNRMLVARY','QVARY','FLUX1VARY','FLUX2VARY',     
     4     'QCONV','QRAD','TSHRSIDE','TSHRVARY','SLIP'/
C
C     ******************************************************************
C
      LL=0
      LC=0
      LINE=0
      MAXLN1=45
      MAXLN2=55
      IBCSET=0
      KURBCN=0
C
C     READ ELEMENT AND BC DATA,  SET UP LOOPING DATA
C
      LC1=-1
      LC2=-1
      IINC=0
      JINC=0
      TYPE1='BLANK'
      TYPE2='BLANK'
   10 CONTINUE
      CALL RDFFLD (NIN,RDATA,IDATA,CDATA)
      IF (LC2.GE.0) LC2=LC2+1
      IF (LC2.LT.0.AND.LC1.GE.0) LC1=LC1+1
      INST=0
   20 CONTINUE
      INST=INST+1
      IF (CDATA(1).NE.CODES(INST).AND.INST.LE.14) GO TO 20
      IF (INST.EQ.15.AND.CDATA(1).NE.CODES(15)) CALL ERROR ('ELEMNT',   
     1'UNRECOGNIZED ELEMENTS COMMAND',' ',0,' ',0,'WORD',CDATA(1),1)
      GO TO (120, 110, 130, 130, 140, 140, 140, 140, 200, 30, 30, 50, 50
     1, 190, 250), INST
   30 CONTINUE
      IF (LC2.GE.0) CALL ERROR ('ELEMNT','I/J LOOP PREVIOUSLY DEFINED', 
     1' ',0,' ',0,'LOOP COMMAND',CDATA(1),1)
      IF (LC1.LT.0) GO TO 40
      LC2=0
      NL2=IDATA(2)
      INC2=IDATA(3)
      TYPE2=CDATA(1)
      GO TO 10
   40 CONTINUE
      LC1=0
      NL1=IDATA(2)
      INC1=IDATA(3)
      TYPE1=CDATA(1)
      GO TO 10
C
   50 CONTINUE
      IF (LC2.GE.0) GO TO 80
      NL1=NL1-1
      IF (NL1.LE.0) GO TO 70
      DO 60 N=1,LC1
      BACKSPACE (NIN)
   60 CONTINUE
      IF (TYPE1.EQ.CODES(10)) IINC=IINC+INC1
      IF (TYPE1.EQ.CODES(11)) JINC=JINC+INC1
      LC1=0
      GO TO 10
   70 CONTINUE
      IINC=0
      JINC=0
      LC1=-1
      LC2=-1
      TYPE1='BLANK'
      TYPE2='BLANK'
      GO TO 10
C
   80 CONTINUE
      NL2=NL2-1
      IF (NL2.LE.0) GO TO 100
      DO 90 N=1,LC2
      BACKSPACE (NIN)
   90 CONTINUE
      IF (TYPE2.EQ.CODES(10)) IINC=IINC+INC2
      IF (TYPE2.EQ.CODES(11)) JINC=JINC+INC2
      LC2=0
      GO TO 10
  100 CONTINUE
      LC1=LC1+LC2
      LC2=-1
      IF (TYPE2.EQ.CODES(10)) IINC=0
      IF (TYPE2.EQ.CODES(11)) JINC=0
      TYPE2='BLANK'
      GO TO 10
C
C     STORE PRELIMINARY ELEMENT DATA IN IJK ARRAY
C
C     TRIANGULAR ELEMENTS
C
  110 CONTINUE
      NN=6
      IJINC=0
      GO TO 160
  120 CONTINUE
      NN=6
      IJINC=0
      IF (IDATA(9).NE.0) GO TO 160
      IDATA (9)=(IDATA(5)-IDATA(3))/2+IDATA(3)
      IDATA (10)=(IDATA(6)-IDATA(4))/2+IDATA(4)
      IDATA (11)=(IDATA(7)-IDATA(5))/2+IDATA(5)
      IDATA (12)=(IDATA(8)-IDATA(6))/2+IDATA(6)
      IDATA (13)=(IDATA(3)-IDATA(7))/2+IDATA(7)
      IDATA (14)=(IDATA(4)-IDATA(8))/2+IDATA(8)
      GO TO 160
C
C     QUADRILATERAL ELEMENTS
C
  130 CONTINUE
      NN=8
      IJINC=2
      GO TO 150
  140 CONTINUE
      NN=9
      IJINC=2
  150 CONTINUE
      IF (IDATA(5).EQ.0.OR.IDATA(11).NE.0) GO TO 160
      IDATA (11)=(IDATA(5)-IDATA(3))/2+IDATA(3)
      IDATA (12)=(IDATA(6)-IDATA(4))/2+IDATA(4)
      IDATA (13)=(IDATA(7)-IDATA(5))/2+IDATA(5)
      IDATA (14)=(IDATA(8)-IDATA(6))/2+IDATA(6)
      IDATA (15)=(IDATA(7)-IDATA(9))/2+IDATA(9)
      IDATA (16)=(IDATA(8)-IDATA(10))/2+IDATA(10)
      IDATA (17)=(IDATA(9)-IDATA(3))/2+IDATA(3)
      IDATA (18)=(IDATA(10)-IDATA(4))/2+IDATA(4)
      IDATA (19)=(IDATA(13)-IDATA(17))/2+IDATA(17)
      IDATA (20)=(IDATA(14)-IDATA(18))/2+IDATA(18)
  160 CONTINUE
      NUMEL=NUMEL+1
      I=IDATA(3)
      J=IDATA(4)
      IF (NUMEL.GT.MAXELM) CALL ERROR ('ELEMNT','MAXIMUM NUMBER OF ELEME
     1NTS EXCEEDED','ELEMENTS SPECIFIED',NUMEL,'ELEMENTS ALLOWED',      
     2MAXELM,' ',' ',1)
      MAT=IDATA(2)
      IJK(1,NUMEL)=100*NN+10000*MAT+1000000*INST
      IJK(2,NUMEL)=10000*(J+JINC)+I+IINC
      IF (IDATA(5).EQ.0) GO TO 180
      DO 170 KK=2,NN
      KK2=2*KK+2
      IJK(1+KK,NUMEL)=10000*(IDATA(KK2)+JINC)+IDATA(KK2-1)+IINC
  170 CONTINUE
      GO TO 10
  180 CONTINUE
      IJK(3,NUMEL)=10000*(J+JINC)+I+IJINC+IINC
      IJK(4,NUMEL)=10000*(J+JINC+IJINC)+I+IJINC+IINC
      IJK(5,NUMEL)=10000*(J+JINC+IJINC)+I+IINC
      IF (NN.LE.4) GO TO 10
      IJK(6,NUMEL)=10000*(J+JINC)+I+IINC+1
      IJK(7,NUMEL)=10000*(J+JINC+1)+I+IJINC+IINC
      IJK(8,NUMEL)=10000*(J+JINC+IJINC)+I+IJINC-1+IINC
      IJK(9,NUMEL)=10000*(J+JINC+IJINC-1)+I+IINC
      IF (NN.LE.8) GO TO 10
      IJK(10,NUMEL)=10000*(J+JINC+1)+I+IINC+1
      GO TO 10
C
C     STORE BOUNDARY CONDITION DATA
C
  190 CONTINUE
      LL=0
      IF (CDATA(2).EQ.BCCODE(26)) LL=1
      IF (CDATA(2).EQ.BCCODE(27)) LL=2
      IF (LL.EQ.0) CALL ERROR ('ELEMNT','UNRECOGNIZED DATA ON SET CARD',
     1' ',0,' ',0,'WORD',CDATA(2),1)
      N=IDATA(3)
      IF (CDATA(4).EQ.'VARIABLE  ') RDATA(4)=1.0E25
      BCSET(LL,N,1)=RDATA(4)
      BCSET(LL,N,2)=RDATA(5)
      IBCSET=IBCSET+1
      GO TO 10
  200 CONTINUE
      I1=IDATA(3)+IINC
      J1=IDATA(4)+JINC
      NCHECK=10000*J1+I1
      DO 210 N=1,NUMEL
      IF (NCHECK.EQ.IJK(2,N)) GO TO 220
  210 CONTINUE
      CALL ERROR ('ELEMNT','BC APPLIED TO AN UNDEFINED ELEMENT','I',I1, 
     1'J',J1,'BC TYPE',CDATA(2),1)
  220 CONTINUE
      IF (KURBCN.EQ.MAXBC) CALL ERROR ('ELEMNT','MAXIMUM NUMBER OF BOUND
     1ARY CONDITIONS EXCEEDED','BC NUMBER',KURBCN,'BC ALLOWED',MAXBC,   
     2' ',' ',1)
      IJK(1,N)=IJK(1,N)+1
      DO 230 L=1,30
      IF (BCCODE(L).EQ.CDATA(2)) GO TO 240
  230 CONTINUE
      CALL ERROR ('ELEMNT','UNRECOGNIZED BOUNDARY CONDITION','I',I1,'J',
     1J1,'BC WORD',CDATA(2),1)
  240 CONTINUE
      IF (L.GE.7.AND.IDATA(5).GT.4) CALL ERROR ('ELEMNT','BC ON IMPROPER
     1 SIDE OF ELEMENT','I',I1,'J',J1,'BC TYPE',CDATA(2),1)
      KURBCN=KURBCN+1
      IBC(KURBCN)=10000000*IDATA(5)+100000*L+N
      BCVAL(KURBCN)=RDATA(6)
      NUMBC=KURBCN
      GO TO 10
C
C     DECODE AND PRINT ELEMENT DATA
C
  250 CONTINUE
      WRITE (NOUT, 450) NUMEL
      IF (IPRINT.GE.2) WRITE (NOUT, 440)
      IF (NUMEL.EQ.0) RETURN
C
C     DETERMINE NODAL COORDINATES
C
      NL1=1
      LL=0
  260 CONTINUE
      IF (ORDER.EQ.'PRESCRIBED') GO TO 280
      LL=0
      NEW=20000000
      DO 270 N=1,NUMEL
      IJ=IJK(2,N)
      IF (IJK(2,N).GT.20000000) GO TO 270
      IF (NEW.LE.IJ) GO TO 270
      LL=N
      NEW=IJ
  270 CONTINUE
      IF (LL.EQ.0) GO TO 380
      GO TO 290
  280 CONTINUE
      LL=LL+1
      IF (LL.GT.NUMEL) GO TO 380
  290 CONTINUE
      KIND=IJK(1,LL)/1000000
      NN=MOD(IJK(1,LL),10000)/100
      MAT=MOD(IJK(1,LL),100000)/10000
      DO 300 N=1,NN
      II(N)=0
      JJ(N)=0
      XX(N)=0.
      YY(N)=0.
  300 CONTINUE
      DO 310 N=1,NN
      I=MOD(IJK(1+N,LL),10000)
      II(N)=I
      J=IJK(1+N,LL)/10000
      JJ(N)=J
      XX(N)=X(I,J)
      YY(N)=Y(I,J)
  310 CONTINUE
      GO TO (320, 340, 330, 340, 330, 340, 330, 340), KIND
  320 CONTINUE
      XX(4)=(XX(1)+XX(2))*.5
      XX(5)=(XX(2)+XX(3))*.5
      XX(6)=(XX(3)+XX(1))*.5
      YY(4)=(YY(1)+YY(2))*.5
      YY(5)=(YY(2)+YY(3))*.5
      YY(6)=(YY(3)+YY(1))*.5
      GO TO 340
  330 CONTINUE
      XX(5)=(XX(1)+XX(2))*.5
      XX(6)=(XX(2)+XX(3))*.5
      XX(7)=(XX(3)+XX(4))*.5
      XX(8)=(XX(4)+XX(1))*.5
      YY(5)=(YY(1)+YY(2))*.5
      YY(6)=(YY(2)+YY(3))*.5
      YY(7)=(YY(3)+YY(4))*.5
      YY(8)=(YY(4)+YY(1))*.5
      XX(9)=-(XX(1)+XX(2)+XX(3)+XX(4))*.25+(XX(5)+XX(6)+XX(7)+XX(8))*.5
      YY(9)=-(YY(1)+YY(2)+YY(3)+YY(4))*.25+(YY(5)+YY(6)+YY(7)+YY(8))*.5
  340 CONTINUE
C
C     COLLECT BOUNDARY CONDITIONS FOR ELEMENT
C
      NL=NN+1
      NUMBCL=MOD(IJK(1,LL),100)
      IF (NUMBCL.EQ.0) GO TO 360
      KBC=1
      DO 350 I=1,NUMBC
      IEL=MOD(IBC(I),100000)
      IF (IEL.NE.LL) GO TO 350
      NBC(KBC)=MOD(IBC(I),10000000)/100000
      NNBC(KBC)=IBC(I)/10000000
      VAL(KBC)=BCVAL(I)
      KBC=KBC+1
  350 CONTINUE
  360 CONTINUE
C
C     WRITE ELEMENT DATA ON NTAPE2
C
      WRITE (NTP2) KIND,MAT,NN,(XX(I),I=1,NN),(YY(I),I=1,NN),(IJK(I,LL),
     1I=2,NL),NUMBCL,(NNBC(I),NBC(I),VAL(I),I=1,NUMBCL)
      NNEL=NNELM(KIND)
      IF (NNEL.GT.MAXNPT) MAXNPT=NNEL
C
C     PRINT ELEMENT TYPE,MATERIAL AND NODAL POINT NUMBERS
C
      IF (IPRINT.LE.2) GO TO 370
      LINE=LINE+4
      IF (LINE.GE.MAXLN1) LINE=4
      IF (NL1.EQ.1) THEN
      LINE=LINE+4
      WRITE (NOUT, 460)
      END IF
      IF (LINE.EQ.4.AND.IPRINT.GE.2) WRITE (NOUT, 470)
      NL=12-NN
      NNL=6
      IF (NL.GE.6) NNL=NN
      WRITE (NOUT, 480) CODES(KIND),MAT,(II(I),JJ(I),I=1,NNL)
      WRITE (NOUT, 490) NL1,(XX(I),YY(I),I=1,NNL)
      IF (NL.GE.6) GO TO 370
      WRITE (NOUT, 500) (II(I),JJ(I),I=7,NN)
      WRITE (NOUT, 510) (XX(I),YY(I),I=7,NN)
  370 CONTINUE
      NL1=NL1+1
      IJK(2,LL)=IJK(2,LL)+20000000
      GO TO 260
C
C     PRINT BOUNDARY CONDITION DATA
C
  380 CONTINUE
      IF (IPRINT.GE.2) WRITE (NOUT, 440)
      IF (IPRINT.EQ.1) GO TO 430
      LINE=6
      DO 400 N=1,NUMEL
      I=MOD(IJK(2,N),10000)
      J=IJK(2,N)/10000
      J=J-2000
      NN=0
      NUMBCL=MOD(IJK(1,N),100)
      IF (NUMBCL.EQ.0) GO TO 400
      KBC=1
      DO 390 K=1,NUMBC
      IEL=MOD(IBC(K),100000)
      IF (IEL.NE.N) GO TO 390
      L=MOD(IBC(K),10000000)/100000
      TYPE(KBC)=BCCODE(L)
      NNBC(KBC)=IBC(K)/10000000
      VAL(KBC)=BCVAL(K)
      KBC=KBC+1
  390 CONTINUE
      NN=NUMBCL
      LINE=LINE+NN+1
      IF (LINE.GE.MAXLN2) LINE=NN+7
      IF (LINE.LE.NN+7) WRITE (NOUT, 520)
      WRITE (NOUT, 530) N,I,J,(TYPE(K),NNBC(K),VAL(K),K=1,NUMBCL)
  400 CONTINUE
      IF (IBCSET.EQ.0) GO TO 430
      WRITE (NOUT, 540)
      DO 420 KK=1,2
      DO 410 K=1,20
      IF (BCSET(KK,K,1).LT.0.) GO TO 420
      IF (BCSET(KK,K,1).GT.1.0E24) THEN
      WRITE (NOUT, 560) KK,K,'VARIABLE',BCSET(KK,K,2)
       ELSE
      WRITE (NOUT, 550) KK,K,BCSET(KK,K,1),BCSET(KK,K,2)
      END IF
  410 CONTINUE
  420 CONTINUE
  430 CONTINUE
      RETURN
C
  440 FORMAT (/////,2X,'************************************************
     1******************************************************************
     2****************'   ,//)
  450 FORMAT (///,3X,'DATA FOR 'I6,' ELEMENTS HAVE BEEN READ' )
  460 FORMAT ('1',/,53X,'--ELEMENT COORDINATE DATA--' ,///,' ELEMENT   M
     1ATR'              ,8X,'NODE 1',14X,'NODE 2',14X,'NODE 3',         
     214X,'NODE 4',14X,'NODE 5',14X,'NODE 6',/,' TYPE      TYPE' ,5X,'I/
     3X'   ,6X,'J/Y',8X,'I/X',6X,'J/Y',8X,'I/X',6X,'J/Y',8X,'I/X',6X,'J/
     4Y'   ,8X,'I/X',6X,'J/Y',8X,'I/X',6X,'J/Y',/,'   NUM   ',14X,'NODE 
     57'      ,14X,'NODE 8',14X,'NODE 9',14X,'NODE 10',13X,'NODE 11',13X
     6,'NODE 12',/,20X,'I/X',6X,'J/Y',8X,'I/X',6X,'J/Y',8X,'I/X',6X,'J/Y
     7'    ,8X,'I/X',6X,'J/Y',8X,'I/X',6X,'J/Y',8X,'I/X',6X,'J/Y'/)
  470 FORMAT ('1ELEMENT   MATR' ,8X,'NODE 1',14X,'NODE 2',14X,'NODE 3', 
     114X,'NODE 4',14X,'NODE 5',14X,'NODE 6',/,' TYPE      TYPE' ,5X,'I/
     2X'   ,6X,'J/Y',8X,'I/X',6X,'J/Y',8X,'I/X',6X,'J/Y',8X,'I/X',6X,'J/
     3Y'   ,8X,'I/X',6X,'J/Y',8X,'I/X',6X,'J/Y',/,'   NUM   ',14X,'NODE 
     47'      ,14X,'NODE 8',14X,'NODE 9',14X,'NODE 10',13X,'NODE 11',13X
     5,'NODE 12',/,20X,'I/X',6X,'J/Y',8X,'I/X',6X,'J/Y',8X,'I/X',6X,'J/Y
     6'    ,8X,'I/X',6X,'J/Y',8X,'I/X',6X,'J/Y',8X,'I/X',6X,'J/Y'/)
  480 FORMAT (/,1X,A8,1X,I4,I9,I9,5(I11,I9))
  490 FORMAT (1X,I5,9X,F9.3,F9.3,5(F11.3,F9.3))
  500 FORMAT (14X,I9,I9,5(I11,I9))
  510 FORMAT (15X,F9.3,F9.3,5(F11.3,F9.3))
  520 FORMAT ('1',///,3X,'BOUNDARY  CONDITIONS' ///,10X,'ELEMENT',4X,   
     1'ELEMENT',4X,'TYPE',5X,'NODE/',5X,'VALUE/',/,11X,'NUMBER',5X,     
     2'I   J',14X,'SIDE'6X,'SET NO.'//)
  530 FORMAT (9X,I5,4X,I5,I4,3X,A10,I4,4X,E12.4/14(30X,A10,I4,4X,E12.4  
     1/))
  540 FORMAT (///,10X,'QCONV=1,QRAD=2' ,7X,'SET NO.',8X,'--H--',8X,'--TB
     1ULK--')
  550 FORMAT (//,13X,I5,11X,I5,6X,2E15.7)
  560 FORMAT (//,13X,I5,11X,I5,13X,A,E15.7)
      END
