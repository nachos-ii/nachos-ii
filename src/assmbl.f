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
      SUBROUTINE ASSMBL (KNUM,ELSTIF,ELFV,ISOLVE,ICON,LISTEL,UN,UPNP1,  
     1ACCN,SCRTCH,XG,YG)
C
C     ******************************************************************
C
C     SUBROUTINE TO DIRECT THE ASSEMBLY OF THE ELEMENT MATRICES AND
C     R.H.S. VECTORS
C
C     ******************************************************************
C
      REAL MASS,KUV,KP,KT,KV1,KV2,KTQ,KTR
C
      COMMON /TAPES/  NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7, 
     1                NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
      COMMON /SZDAT/  NUMEL,NUMNOD,NUMVAR,NUMDOF,NODSOL
      COMMON /MATDAT/ NMAT,NPROP,PROP(15,10),NXMAT,NXPROP,XPROP(15,10)
      COMMON /PRBDAT/ IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,    
     1                IPNLTY,PNLTY
      COMMON /ELMDAT/ NNELM(6),NNCOR(6),NNMID(6),NNCTR(6),NSIDET(5),    
     1                NSIDEQ(7),NNSIDE(6,4,3)
      COMMON /SLNDAT/ IALGOR,IBLDUV,IBLDT,IJACOB,IAUTO,NCYCLM,NPRNT,    
     1                KSTEP,KPRNT,TIME,DELTN,DELTNM
C
      DIMENSION ICON(NUMEL,*), LISTEL(*), XG(*), YG(*)
      DIMENSION UN(NUMNOD,*), UPNP1(NUMNOD,*), ACCN(NUMNOD,*)
      DIMENSION SCRTCH(NUMNOD,*)
      DIMENSION ELSTIF(50,50), ELFV(50)
      DIMENSION KUV(22,22,4), CX(9,9,9), CY(9,9,9), C(9,9,9), MASS(9,9)
      DIMENSION FUV(22), FBDY(9,9,4), KP(18,18)
      DIMENSION KT(9,9,4), KTQ(9,9), KTR(9,3,9)
      DIMENSION FT(9), FTQ(9), FTR(9,3), FTS(9,4)
      DIMENSION KV1(9,9,4), KV2(9,9,4), FV1(9), FV2(9)
      DIMENSION CP(9,9), CPV1(9,9), CPV2(9,9)
      DIMENSION CXY(9,9), CXYP(18,18), CXYT(9,9), CXYTP(9,18)
      DIMENSION CXYV1(9,9), CXYV2(9,9), CXYV1P(9,18), CXYV2P(9,18)
      DIMENSION X(9), Y(9)
      DIMENSION UNELM(9,6), UPELM(9,6), UDOT(9,6)
C
C     ******************************************************************
C
C     READ COMPONENT MATRICES FROM TAPE
C
      READ (NTP4) NSZ1,NSZ2,NSZ3,NSZ4,NSZ5,IBCPT1,IBCPT2,IBCPT3,        
     1KUV,CX,CY,MASS,KP,C,FUV,FBDY,KT,KTQ,KTR,FT,FTQ,FTR,FTS,           
     2KV1,FV1,KV2,FV2
C
      DO 10 I=1,NSZ2
      DO 10 J=1,NSZ2
      CP(I,J)=MASS(I,J)
      CPV1(I,J)=MASS(I,J)
      CPV2(I,J)=MASS(I,J)
   10 CONTINUE
C
C     RETRIEVE CURRENT VALUES FOR ELEMENT UNKNOWNS
C
      NN=ABS(LISTEL(KNUM))
      MAT=NN/100
      KIND=MOD(NN,100)
      NN=NNELM(KIND)
      NCOR=NNCOR(KIND)
      DO 20 I=1,NN
      INODE=ABS(ICON(KNUM,I))
      X(I)=XG(INODE)
      Y(I)=YG(INODE)
      DO 20 J=1,NUMVAR
      UNELM(I,J)=UN(INODE,J)
      UPELM(I,J)=UPNP1(INODE,J)
      UDOT(I,J)=ACCN(INODE,J)
   20 CONTINUE
      IDARCY=0
      IF (PROP(15,MAT).EQ.4.0) IDARCY=1
C
C     CALL ASSEMBLY ROUTINES FOR THE APPROPRIATE EQUATIONS
C
      IF (IBLDUV.EQ.1.AND.IDARCY.EQ.0) THEN
      CALL ASSMU1 (KNUM,KIND,MAT,NN,NCOR,IBCPT1,MASS,CX,CY,KUV,FUV,KP,  
     1FBDY,X,Y,UPELM(1,1),UPELM(1,2),UPELM(1,3),UPELM(1,4),UPELM(1,5),  
     2UPELM(1,6),CXY,CXYP,GX,GY)
      END IF
C
      IF (IBLDUV.EQ.1.AND.IDARCY.EQ.1) THEN
      CALL ASSMU2 (KNUM,KIND,MAT,NN,NCOR,IBCPT1,MASS,C,KUV,FUV,KP,      
     1FBDY,X,Y,UPELM(1,1),UPELM(1,2),UPELM(1,3),UPELM(1,4),UPELM(1,5),  
     2UPELM(1,6),CXY,CXYP,GX,GY)
      END IF
C
      IF (IBLDT.EQ.1) THEN
      CALL ASSMT (KNUM,KIND,MAT,NN,NCOR,IBCPT2,IBCPT3,CP,CX,CY,KT,KTQ,  
     1KTR,FT,FTQ,FTR,FTS,X,Y,UPELM(1,1),UPELM(1,2),UPELM(1,3),          
     2UPELM(1,4),UPELM(1,5),UPELM(1,6),CXYT,CXYTP)
      END IF
C
      IF (IBLDT.EQ.1.AND.IVAR1.EQ.1) THEN
      CALL ASSMVI (KNUM,1,KIND,MAT,NN,NCOR,IBCPT2,CPV1,CX,CY,KV1,FV1,   
     1FTS,X,Y,UPELM(1,1),UPELM(1,2),UPELM(1,3),UPELM(1,4),UPELM(1,5),   
     2UPELM(1,6),CXYV1,CXYV1P)
      END IF
C
      IF (IBLDT.EQ.1.AND.IVAR2.EQ.1) THEN
      CALL ASSMVI (KNUM,2,KIND,MAT,NN,NCOR,IBCPT2,CPV2,CX,CY,KV2,FV2,   
     1FTS,X,Y,UPELM(1,1),UPELM(1,2),UPELM(1,3),UPELM(1,4),UPELM(1,5),   
     2UPELM(1,6),CXYV2,CXYV2P)
      END IF
C
C     ASSEMBLE EQUATIONS FOR APPROPRIATE SOLUTION OPTION
C
      CALL ASSMKF (KNUM,KIND,NN,NCOR,MASS,CXY,CXYP,KUV,FUV,FBDY,CP,     
     1CXYT,CXYTP,KT,FT,CPV1,CPV2,CXYV1,CXYV2,CXYV1P,CXYV2P,KV1,KV2,     
     2FV1,FV2,UNELM(1,1),UNELM(1,2),UNELM(1,3),UNELM(1,4),UNELM(1,5),   
     3UNELM(1,6),UPELM(1,1),UPELM(1,2),UPELM(1,3),UPELM(1,4),UPELM(1,5),
     4UPELM(1,6),UDOT(1,1),UDOT(1,2),UDOT(1,3),UDOT(1,4),UDOT(1,5),     
     5UDOT(1,6),GX,GY,ELSTIF,ELFV)
C
C     RE-ORDER ELEMENT MATRIX FOR USE IN FRONTAL SOLVER
C
      CALL ELORDR (KNUM,ELSTIF,ELFV,LISTEL)
C
      RETURN
      END
