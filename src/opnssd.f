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
      SUBROUTINE OPNSSD
C
C     ******************************************************************
C
C     SUBROUTINE TO PUT REQUIRED FILES ON A SOLID STATE DISK
C     THIS ROUTINE USES CFTLIB ROUTINES IN THE CTSS OPERATING SYSTEM
C
C     ******************************************************************
C
      COMMON /TAPES/  NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7, 
     1                NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
C
      DATA DUMWRD/1.23E4/
C
C     ******************************************************************
C
C     PLACE ELEMENT MATRIX FILE ON SSD
C
C      CALL SETSSD (NTP4)
C
C     PLACE SOLUTION FILE ON SSD
C
C      CALL SETSSD (NTP9)
C
C     PLACE SCRATCH FILE FOR FRONTAL SOLVER ON SSD
C
C      CALL SETSSD (NTP10)
C
C     CHECK FOR FILE STATUS
C     WRITE TO UNIT TO FORCE FILE TO BE ASSIGNED
C
C      REWIND (NTP4)
C      WRITE (NTP4) DUMWRD
C      REWIND (NTP4)
C      CALL GETDU (NTP4,IUNIT1)
C      IF (IUNIT1.EQ.37B) THEN
C      WRITE (NOUT, 10)
C       ELSE
C      WRITE (NOUT, 20)
C      END IF
C
C      REWIND (NTP9)
C      WRITE (NTP9) DUMWRD
C      REWIND (NTP9)
C      CALL GETDU (NTP9,IUNIT2)
C      IF (IUNIT2.EQ.37B) THEN
C      WRITE (NOUT, 30)
C       ELSE
C      WRITE (NOUT, 40)
C      END IF
C
C      REWIND (NTP10)
C      WRITE (NTP10) DUMWRD
C      REWIND (NTP10)
C      CALL GETDU (NTP10,IUNIT3)
C      IF (IUNIT3.EQ.37B) THEN
C      WRITE (NOUT, 50)
C       ELSE
C      WRITE (NOUT, 60)
C      END IF
C
C   10 FORMAT (/,10X,'THE ELEMENT MATRIX FILE IS ON THE SSD')
C   20 FORMAT (/,10X,'THE ELEMENT MATRIX FILE IS NOT ON THE SSD')
C   30 FORMAT (/,10X,'THE SOLUTION FILE IS ON THE SSD')
C   40 FORMAT (/,10X,'THE SOLUTION FILE IS NOT ON THE SSD')
C   50 FORMAT (/,10X,'THE SOLUTION SCRATCH FILE IS ON THE SSD')
C   60 FORMAT (/,10X,'THE SOLUTION SCRATCH FILE IS NOT ON THE SSD')
C
      RETURN
      END
