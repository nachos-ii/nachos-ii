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
      SUBROUTINE TNORM (ITYPE,UPNP1,SCRTCH)
C
C     *****************************************************************
C
C     SUBROUTINE TO COMPUTE NORMS FOR COMPUTATION OF NEW TIME STEP
C
C     *****************************************************************
C
      COMMON /TAPES/  NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7, 
     1                NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
      COMMON /SZDAT/  NUMEL,NUMNOD,NUMVAR,NUMDOF,NODSOL
      COMMON /PRBDAT/ IAXSYM,ITMDEP,IFORCE,IFREE,IVAR1,IVAR2,IPFUNC,    
     1                IPNLTY,PNLTY
      COMMON /NORM/   UMAX,VMAX,TMAX,V1MAX,V2MAX
      COMMON /BNORM/  BNRMU,BNRMV,BNRMUV,BNRMT,BNRMX1,BNRMX2
C
      DIMENSION UPNP1(NUMNOD,*), SCRTCH(NUMNOD,*)
C
C     *****************************************************************
C
C     COMPUTE VELOCITY NORMS
C
      USUM=0.
      VSUM=0.
      DO 10 I=1,NUMNOD
      A=(SCRTCH(I,1)-UPNP1(I,1))**2
      B=(SCRTCH(I,2)-UPNP1(I,2))**2
      USUM=USUM+A
      VSUM=VSUM+B
   10 CONTINUE
      UVMAX=SQRT(UMAX**2+VMAX**2)
      IF (UVMAX.LT.1.0E-10) UVMAX=1.0
      BNRMU=SQRT(USUM)/UVMAX
      BNRMV=SQRT(VSUM)/UVMAX
      BNRMUV=SQRT(BNRMU**2+BNRMV**2)
      WRITE (NOUT, 50) BNRMU,BNRMV,BNRMUV
C
C     COMPUTE TEMPERATURE NORM
C
      IF (ITYPE.EQ.1) RETURN
      TSUM=0.
      DO 20 I=1,NUMNOD
      A=(SCRTCH(I,4)-UPNP1(I,4))**2
      TSUM=TSUM+A
   20 CONTINUE
      BNRMT=SQRT(TSUM)/ABS(TMAX)
      WRITE (NOUT, 60) BNRMT
C
C     COMPUTE AUXILIARY VARIABLE NORMS
C
      BNRMX1=0.
      BNRMX2=0.
      IF (IVAR1.EQ.0) RETURN
      XSUM1=0.
      DO 30 I=1,NUMNOD
      A=(SCRTCH(I,5)-UPNP1(I,5))**2
      XSUM1=XSUM1+A
   30 CONTINUE
      BNRMX1=SQRT(XSUM1)/ABS(V1MAX)
      WRITE (NOUT, 70) BNRMX1
C
      IF (IVAR2.EQ.0) RETURN
      XSUM2=0.
      DO 40 I=1,NUMNOD
      A=(SCRTCH(I,6)-UPNP1(I,6))**2
      XSUM2=XSUM2+A
   40 CONTINUE
      BNRMX2=SQRT(XSUM2)/ABS(V2MAX)
      WRITE (NOUT, 80) BNRMX2
C
      RETURN
C
   50 FORMAT (/,10X,'NORMS ON INTEGRATION ERROR -',2X,' U NORM='        
     1E15.7,2X,'V NORM=',E15.7,3X,'UV NORM=',E15.7)
   60 FORMAT (/,10X,'NORM ON INTEGRATION ERROR -',2X,' T NORM=',        
     1E15.7)
   70 FORMAT (/,10X,'NORM ON INTEGRATION ERROR -',2X,' X1 NORM='        
     1E15.7)
   80 FORMAT (/,10X,'NORM ON INTEGRATION ERROR -',2X,' X2 NORM='        
     1E15.7)
      END
