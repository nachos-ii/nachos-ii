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
      SUBROUTINE MATRL
C
C     ******************************************************************
C
C     SUBROUTINE TO READ AND STORE MATERIAL PROPERTY DATA
C
C     ******************************************************************
C
      CHARACTER*10 CDATA,XNAME
      CHARACTER*20 A1,A2,A3,A4
C
      COMMON /INDATR/ RDATA(150)
      COMMON /INDATI/ IDATA(150)
      COMMON /INDATC/ CDATA(150)
      COMMON /TAPES/  NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7, 
     1                NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
      COMMON /LMTDAT/ MAXELM,MAXNOD,MAXBC,MAXEDG,MAXMAT,MAXNPT,MAXPLT,  
     1                MAXBRY,MAXBLK
      COMMON /MATDAT/ NMAT,NPROP,PROP(15,10),NXMAT,NXPROP,XPROP(15,10)
C
      DIMENSION XNAME(10)
C
C     ******************************************************************
C
      N=0
      M=0
   10 CONTINUE
      CALL RDFFLD (NIN,RDATA,IDATA,CDATA)
      IF (CDATA(1).EQ.'END       ') GO TO 50
      IF (CDATA(2).EQ.'AUXDATA   ') GO TO 30
C
C     READ FLUID, SOLID AND POROUS MATERIAL DATA
C
      N=N+1
      IF (N.GT.MAXMAT) CALL ERROR ('MATRL','TOO MANY MATERIALS SPECIFIED
     1','MAT NUMBER',N,'LIMIT',10,' ',' ',1)
      NCOUNT=N
      MAT=IDATA(3)
      XNAME(MAT)=CDATA(1)
      DO 20 I=1,12
      PROP(I,MAT)=RDATA(I+3)
   20 CONTINUE
      PROP(15,MAT)=1.0
      IF (CDATA(2).EQ.'NNEWTONIAN') PROP(15,MAT)=2.0
      IF (CDATA(2).EQ.'SOLID     ') PROP(15,MAT)=3.0
      IF (CDATA(2).EQ.'POROUS    ') PROP(15,MAT)=4.0
      IF (CDATA(9).EQ.'VARIABLE  ') PROP(6,MAT)=-1.0
      PROP(8,MAT)=1.0
      IF (CDATA(11).EQ.'VARIABLE  ') PROP(8,MAT)=-1.0
      IF (CDATA(2).EQ.'NNEWTONIAN') PROP(8,MAT)=-1.0
      IF (CDATA(12).EQ.'VARIABLE  ') PROP(9,MAT)=-1.0
      IF (PROP(15,MAT).EQ.4.0) GO TO 10
      PROP(10,MAT)=1.0
      IF (CDATA(13).EQ.'DISSIPATE ') PROP(10,MAT)=-1.0
      GO TO 10
C
C     READ AUXILIARY EQUATION DATA
C
   30 CONTINUE
      M=M+1
      MAT=IDATA(3)
      PROP(14,MAT)=1.0
      DO 40 I=1,12
      XPROP(I,MAT)=RDATA(I+3)
   40 CONTINUE
      IF (CDATA(7).EQ.'VARIABLE  ') XPROP(4,MAT)=-1.0
      IF (CDATA(13).EQ.'VARIABLE  ') XPROP(10,MAT)=-1.0
      GO TO 10
   50 CONTINUE
      NMAT=N
      NPROP=15
      NXMAT=M
      NXPROP=15
C
C     PRINT MATERIAL DATA
C
      WRITE (NOUT, 150) NMAT
      DO 140 I=1,N
      MAT=I
      PR=0.
      IF (PROP(4,MAT).NE.0) PR=PROP(2,MAT)*PROP(3,MAT)/PROP(4,MAT)
      IF (PROP(15,MAT).EQ.3.0) PR=0.
      A1='CONSTANT'
      IF (PROP(8,MAT).EQ.-1.0) A1='VARIABLE'
      A2='NEWTONIAN FLUID'
      IF (PROP(15,MAT).EQ.2.0) A2='NON-NEWTONIAN FLUID'
      IF (PROP(15,MAT).EQ.3.0) A2='SOLID'
      IF (PROP(15,MAT).EQ.4.0) GO TO 60
      A3='OMITTED'
      IF (PROP(10,MAT).LT.0.) A3='INCLUDED'
      IF (PROP(9,MAT).LT.0.) THEN
      A4='VARIABLE'
      WRITE (NOUT, 170) I,XNAME(I),A2,(PROP(J,I),J=1,5),A1,A4,A3
       ELSE
      WRITE (NOUT, 160) I,XNAME(I),A2,(PROP(J,I),J=1,5),A1,PROP(9,I),A3
      END IF
      GO TO 70
   60 CONTINUE
      A2='POROUS MEDIUM'
      PROP(12,MAT)=PROP(10,MAT)
      IF (PROP(9,MAT).LT.0.) THEN
      A4='VARIABLE'
      WRITE (NOUT, 270) I,XNAME(I),A2,INT(PROP(1,I)),(PROP(J,I),J=2,7), 
     1 A1,A4
       ELSE
      WRITE (NOUT, 280) I,XNAME(I),A2,INT(PROP(1,I)),(PROP(J,I),J=2,7), 
     1 A1,PROP(9,MAT)
      END IF
   70 CONTINUE
      IF (PROP(14,MAT).GT.0.) GO TO 80
      GO TO 90
   80 CONTINUE
      M=1
      IF (XPROP(7,MAT)+XPROP(8,MAT).NE.0.) M=2
      WRITE (NOUT, 200) M
      IF (XPROP(4,MAT).GE.0.) THEN
      WRITE (NOUT, 210) (XPROP(J,I),J=1,4)
       ELSE
      A1='VARIABLE'
      WRITE (NOUT, 220) (XPROP(J,I),J=1,3),A1
      END IF
      IF (M.EQ.1) GO TO 90
      IF (XPROP(10,MAT).GE.0.) THEN
      WRITE (NOUT, 230) (XPROP(J,I),J=7,10)
       ELSE
      A1='VARIABLE'
      WRITE (NOUT, 240) (XPROP(J,I),J=7,9),A1
      END IF
   90 CONTINUE
      IF (PROP(15,I).EQ.4.0) GO TO 130
      IF (PROP(6,I).EQ.-1.0) GO TO 100
      WRITE (NOUT, 180) PROP(6,I),PROP(7,I),PROP(11,I),PR
      GO TO 110
  100 CONTINUE
      WRITE (NOUT, 190) PROP(11,I),PR
  110 CONTINUE
      IF (PROP(14,MAT).GT.0.0) GO TO 120
      GO TO 130
  120 CONTINUE
      IF (M.EQ.1) THEN
      WRITE (NOUT, 250) XPROP(5,MAT)
       ELSE
      WRITE (NOUT, 250) XPROP(5,MAT)
      WRITE (NOUT, 260) XPROP(11,MAT)
      END IF
  130 CONTINUE
      NCOUNT=NCOUNT-1
      IF (PROP(15,MAT).EQ.3.0) PROP(2,MAT)=1.0E10
  140 CONTINUE
      RETURN
C
  150 FORMAT (///,5X,'MATERIAL PROPERTY DATA FOR ',I2,
     1' MATERIAL(S) -',//)
  160 FORMAT (12X,'MATERIAL NUMBER ',I2,//,15X,'MATERIAL NAME ..........
     1..................... ',A10,/,15X,'MATERIAL TYPE .................
     2.............. ',A20,/,15X,'DENSITY ..............................
     3....... ',E15.7,/,15X,'VISCOSITY .................................
     4.. ',E15.7,/,15X,'SPECIFIC HEAT ............................... ',
     5E15.7,/,15X,'CONDUCTIVITY ................................ ',     
     6E15.7,/,15X,'VOLUMETRIC EXPANSION COEFFICIENT ............ ',E15.7
     7,/,15X,'PROPERTY VARIATION .......................... ',A20,/,15X,
     8'VOLUMETRIC HEAT SOURCE ...................... ',E15.7,/,15X,'VISC
     9OUS DISSIPATION ......................... ',A20)
  170 FORMAT (12X,'MATERIAL NUMBER ',I2,//,15X,'MATERIAL NAME ..........
     1..................... ',A10,/,15X,'MATERIAL TYPE .................
     2.............. ',A20,/,15X,'DENSITY ..............................
     3....... ',E15.7,/,15X,'VISCOSITY .................................
     4.. ',E15.7,/,15X,'SPECIFIC HEAT ............................... ',
     5E15.7,/,15X,'CONDUCTIVITY ................................ ',     
     6E15.7,/,15X,'VOLUMETRIC EXPANSION COEFFICIENT ............ ',E15.7
     7,/,15X,'PROPERTY VARIATION .......................... ',A20,/,15X,
     8'VOLUMETRIC HEAT SOURCE ...................... ',A20,/,15X,       
     9'VISCOUS DISSIPATION ......................... ',A20)
  180 FORMAT (//,15X,'REFERENCE QUANTITIES -',//,15X,'GRAVITATIONAL VECT
     1OR :',/,18X,'X COMPONENT .............................. ',E15.7,/,
     218X,'Y COMPONENT .............................. ',E15.7,/,15X,    
     3'REFERENCE TEMPERATURE ....................... ',E15.7,/,15X,'PRAN
     4DTL NUMBER .............................. ',E15.7,//)
  190 FORMAT (//,15X,'REFERENCE QUANTITIES -',//,15X,'GRAVITATIONAL VECT
     1OR :',/,18X,'MAGNITUDE ................................... VARIABL
     2E',/,18X,'ORIENTATION ................................. VARIABLE',
     3/,15X,'REFERENCE TEMPERATURE ....................... ',E15.7,/,   
     415X,'PRANDTL NUMBER .............................. ',E15.7,//)
  200 FORMAT (//,10X,'PROPERTY DATA FOR ',I1,' AUXILIARY TRANSPORT EQUAT
     1ION(S)',//)
  210 FORMAT (10X,'VARIABLE 1 DATA -',/,15X,'CAPACITANCE ...............
     1.................. ',E15.7,/,15X,'DIFFUSION COEFFICIENT ..........
     2............. ',E15.7,/,15X,'VOLUMETRIC EXPANSION COEFFICIENT ....
     3........ ',E15.7,/,15X,'VOLUMETRIC SOURCE ........................
     4... ',E15.7)
  220 FORMAT (10X,'VARIABLE 1 DATA -',/,15X,'CAPACITANCE ...............
     1.................. ',E15.7,/,15X,'DIFFUSION COEFFICIENT ..........
     2............. ',E15.7,/,15X,'VOLUMETRIC EXPANSION COEFFICIENT ....
     3........ ',E15.7,/,15X,'VOLUMETRIC SOURCE ........................
     4... ',A20)
  230 FORMAT (10X,'VARIABLE 2 DATA -',/,15X,'CAPACITANCE ...............
     1.................. ',E15.7,/,15X,'DIFFUSION COEFFICIENT ..........
     2............. ',E15.7,/,15X,'VOLUMETRIC EXPANSION COEFFICIENT ....
     3........ ',E15.7,/,15X,'VOLUMETRIC SOURCE ........................
     4... ',E15.7)
  240 FORMAT (10X,'VARIABLE 2 DATA -',/,15X,'CAPACITANCE ...............
     1.................. ',E15.7,/,15X,'DIFFUSION COEFFICIENT ..........
     2............. ',E15.7,/,15X,'VOLUMETRIC EXPANSION COEFFICIENT ....
     3........ ',E15.7,/,15X,'VOLUMETRIC SOURCE ........................
     4... ',A20)
  250 FORMAT (15X,'REFERENCE VARIABLE 1 ........................ ',     
     1E15.7)
  260 FORMAT (15X,'REFERENCE VARIABLE 2 ........................ ',     
     1E15.7)
  270 FORMAT (12X,'MATERIAL NUMBER ',I2,//,15X,'MATERIAL NAME ..........
     1..................... ',A10,/,15X,'MATERIAL TYPE .................
     2.............. ',A20,/,15X,'SATURATING FLUID (NUMBER) ............
     3.......  ',I2,/,15X,'PERMEABLILITY ...............................
     4 ',E15.7,/,15X,'EFFECTIVE (BRINKMAN) VISCOSITY .............. ',  
     5E15.7,/,15X,'POROSITY .................................... ',E15.7
     6,/,15X,'INERTIA COEFFICIENT ......................... ',E15.7,/,  
     715X,'EFFECTIVE CAPACITANCE ....................... ',E15.7,/,15X, 
     8'EFFECTIVE CONDUCTIVITY ...................... ',E15.7,/,15X,     
     9'PROPERTY VARIATION ........................... ',A20,/,15X,'VOLUM
     $ETRIC HEAT SOURCE ...................... ',A20,//,15X,'REFERENCE Q
     $UANTITIES -',//,15X,'SEE PARAMETERS FOR SATURATING FLUID',//)
  280 FORMAT (12X,'MATERIAL NUMBER ',I2,//,15X,'MATERIAL NAME ..........
     1..................... ',A10,/,15X,'MATERIAL TYPE .................
     2.............. ',A20,/,15X,'SATURATING FLUID (NUMBER) ............
     3.......  ',I2,/,15X,'PERMEABLILITY ...............................
     4 ',E15.7,/,15X,'EFFECTIVE (BRINKMAN) VISCOSITY .............. ',  
     5E15.7,/,15X,'POROSITY .................................... ',E15.7
     6,/,15X,'INERTIA COEFFICIENT ......................... ',E15.7,/,  
     715X,'EFFECTIVE CAPACITANCE ....................... ',E15.7,/,15X, 
     8'EFFECTIVE CONDUCTIVITY ...................... ',E15.7,/,15X,     
     9'PROPERTY VARIATION ........................... ',A20,/,15X,'VOLUM
     $ETRIC HEAT SOURCE ...................... ',E15.7,//,15X,'REFERENCE
     $ QUANTITIES -',//,15X,'SEE PARAMETERS FOR SATURATING FLUID',//)
      END
