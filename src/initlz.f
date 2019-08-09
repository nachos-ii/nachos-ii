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
      BLOCK DATA INITLZ
C
C     ******************************************************************
C
C     BLOCK DATA SUBROUTINE TO INITIALIZE VARIABLES STORED IN
C     NAMED COMMON BLOCKS
C
C     ******************************************************************
C
      CHARACTER*8  CODNAM,VERSN,RDATE,RTIME,HRDWRE,SFTWRE
      CHARACTER*10 CWORD,CDATA
      CHARACTER*80 HED,CMMNT
C
      COMMON /HEADER/ HED,CMMNT(10)
      COMMON /RUNDAT/ CODNAM,VERSN,RDATE,RTIME,HRDWRE,SFTWRE
      COMMON /INDATR/ RDATA(150)
      COMMON /INDATI/ IDATA(150)
      COMMON /INDATC/ CDATA(150)
      COMMON /TAPES/  NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7, 
     1                NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
      COMMON /NTPDAT/ IFILES(16)
      COMMON /CONTRL/ IEDIT,IPRINT
      COMMON /RSTART/ IRSTRT,NSTEPS
      COMMON /MEMDAT/ MEMTOT,MAXMEM
      COMMON /SZDAT/  NUMEL,NUMNOD,NUMVAR,NUMDOF,NODSOL
      COMMON /LMTDAT/ MAXELM,MAXNOD,MAXBC,MAXEDG,MAXMAT,MAXNPT,MAXPLT,  
     1                MAXBRY,MAXBLK
      COMMON /GENDAT/ NELBLK,NUMNPS,LNPSNL,NUMESS,LESSEL,LESSNL,MXEBLK, 
     1                IBLK(100)
      COMMON /MATDAT/ NMAT,NPROP,PROP(15,10),NXMAT,NXPROP,XPROP(15,10)
      COMMON /PRBDAT/ IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,    
     1                IPNLTY,PNLTY
      COMMON /ELMDAT/ NNELM(6),NNCOR(6),NNMID(6),NNCTR(6),NSIDET(5),    
     1                NSIDEQ(7),NNSIDE(6,4,3)
      COMMON /SPTDAT/ ISPT,NSPT,NELSPT(50),PTS(50,2),SPTVAL(50,6)
      COMMON /BCDAT/  BCSET(2,20,2),IPELM,IPNODE,PBC
      COMMON /GAUSS1/ GSPT1(4),GSPT2(4)
      COMMON /GAUSS2/ GSPT3(3),GSWT3(3)
      COMMON /HAMMR1/ HAMPT1(3),HAMPT2(3)
      COMMON /HAMMR2/ HMPT1(7),HMPT2(7),HMWT(7)
      COMMON /SLNDAT/ IALGOR,IBLDUV,IBLDT,IJACOB,IAUTO,NCYCLM,NPRNT,    
     1                KSTEP,KPRNT,TIME,DELTN,DELTNM
      COMMON /FRNT/   RPIVOT,SPIVOT,IRESOL,ISOLVE,NSUM,IPIVOT,MWGA
      COMMON /FRNTSZ/ MXFRNT,MXDOFE,MXBLKS,MXBUFF,LIWK,LRWK
      COMMON /NORM/   UMAX,VMAX,TMAX,V1MAX,V2MAX
      COMMON /ANORM/  ANRMU,ANRMV,ANRMUV,ANRMT,ANRMX1,ANRMX2
      COMMON /BNORM/  BNRMU,BNRMV,BNRMUV,BNRMT,BNRMX1,BNRMX2
      COMMON /FLXDAT/ NFLUX
      COMMON /CMMNDS/ CWORD(10)
C
      DATA HED/' '/
      DATA (CMMNT(I),I=1,10)/10*'END'/
      DATA CODNAM,VERSN,RDATE,RTIME,HRDWRE,SFTWRE/'NACHOS','V 03.00C',  
     1     ' ',' ',' ',' '/
      DATA (RDATA(I),I=1,150)/150*0.0/
      DATA (IDATA(I),I=1,150)/150*0/
      DATA (CDATA(I),I=1,150)/150*' '/
      DATA NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7,NTP8,NTP9,  
     1     NTP10,NTP11,NTP12,NTP13/5,6,10,11,12,13,14,15,16,17,18,19,   
     2     20,21,22,23/
      DATA (IFILES(I),I=1,16)/16*0/
      DATA IEDIT,IPRINT/0,0/
      DATA IRSTRT,NSTEPS/0,0/
      DATA MEMTOT,MAXMEM/0,120000/
      DATA NUMEL,NUMNOD,NUMVAR,NUMDOF,NODSOL/0,0,3,0,0/
      DATA MAXELM,MAXNOD,MAXBC,MAXEDG,MAXMAT,MAXNPT,MAXPLT,MAXBRY,      
     1     MAXBLK/0,0,3000,1000,10,0,1000,0,100/
      DATA NELBLK,NUMNPS,LNPSNL,NUMESS,LESSEL,LESSNL,MXEBLK/7*0/
      DATA (IBLK(I),I=1,100)/100*0/
      DATA NMAT,NPROP,NXMAT,NXPROP/0,0,0,0/
      DATA ((PROP(I,J),I=1,15),J=1,10)/150*0./
      DATA ((XPROP(I,J),I=1,15),J=1,10)/150*0./
      DATA IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,IPNLTY/        
     1     0,0,0,0,0,0,0,0/
      DATA PNLTY/0.0/
      DATA (NNELM(I),I=1,6)/6,6,8,8,9,9/
      DATA (NNCOR(I),I=1,6)/3,3,4,4,4,4/
      DATA (NNMID(I),I=1,6)/3,3,4,4,4,4/
      DATA (NNCTR(I),I=1,6)/0,0,0,0,1,1/
      DATA (NSIDET(I),I=1,5)/1,2,3,1,2/
      DATA (NSIDEQ(I),I=1,7)/1,2,3,4,1,2,3/
      DATA (((NNSIDE(I,J,K),I=1,6),J=1,4),K=1,3)/1,1,1,1,1,1,           
     1     2,2,2,2,2,2, 3,3,3,3,3,3, 0,0,4,4,4,4, 4,4,5,5,5,5,          
     2     5,5,6,6,6,6, 6,6,7,7,7,7, 0,0,8,8,8,8, 2,2,2,2,2,2,          
     3     3,3,3,3,3,3, 1,1,4,4,4,4, 0,0,1,1,1,1/
      DATA ISPT,NSPT,(NELSPT(I),I=1,50),((PTS(I,J),I=1,50),J=1,2),      
     1     ((SPTVAL(I,J),I=1,50),J=1,6)/52*0,400*0./
      DATA (((BCSET(I,J,K),I=1,2),J=1,20),K=1,2)/80*-1./
      DATA IPELM,IPNODE,PBC/0,0,0.0/
      DATA (GSPT1(I),I=1,4)/-.5773502691,.5773502691,.5773502691,       
     1     -.5773502691/
      DATA (GSPT2(I),I=1,4)/-.5773502691,-.5773502691,.5773502691,      
     1     .5773502691/
      DATA (GSPT3(I),I=1,3)/-.7745966692,0.,.7745966692/
      DATA (GSWT3(I),I=1,3)/.5555555556,.8888888889,.5555555556/
      DATA (HAMPT1(I),I=1,3)/.7333333333,.1333333333,.1333333333/
      DATA (HAMPT2(I),I=1,3)/.1333333333,.7333333333,.1333333333/
      DATA (HMPT1(I),I=1,7)/.3333333333,.0597158717,.4701420641,        
     1     .4701420641,.7974269853,.1012865073,.1012865073/
      DATA (HMPT2(I),I=1,7)/.3333333333,.4701420641,.0597158717,        
     1     .4701420641,.1012865073,.7974269853,.1012865073/
      DATA (HMWT(I),I=1,7)/.11250000,3*.06619707,3*.06296959/
      DATA IALGOR,IBLDUV,IBLDT,IJACOB,IAUTO,NCYCLM,NPRNT,KSTEP,KPRNT,   
     1     TIME,DELTN,DELTNM/9*0,3*0./
      DATA RPIVOT,SPIVOT,IRESOL,ISOLVE,NSUM,IPIVOT,MWGA/1.0E-4,1.0E-25, 
     1     0,0,20,0,0/
      DATA MXFRNT,MXDOFE,MXBLKS,MXBUFF,LIWK,LRWK/0,0,1000,0,0,0/
      DATA UMAX,VMAX,TMAX,V1MAX,V2MAX/5*0./
      DATA ANRMU,ANRMV,ANRMUV,ANRMT,ANRMX1,ANRMX2/6*0./
      DATA BNRMU,BNRMV,BNRMUV,BNRMT,BNRMX1,BNRMX2/6*0./
      DATA NFLUX/0/
      DATA (CWORD(I),I=1,10)/'MATERIALS','MESH','ELEMENTS','FORMKF',    
     1     'OUTPUT','SOLVE','STREAM','FLUX','POST','STOP'/
C
      END
