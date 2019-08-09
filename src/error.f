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
      SUBROUTINE ERROR (SUBNAM,MESSAG,LABEL1,I,LABEL2,J,LABEL3,WORD,    
     1ISTOP)
C
C     ******************************************************************
C
C     SUBROUTINE TO PRINT ERROR MESSAGE AND TERMINATE EXECUTION
C
C     ******************************************************************
C
      CHARACTER*(*) SUBNAM,MESSAG,LABEL1,LABEL2,LABEL3,WORD
C
      COMMON /TAPES/ NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7,  
     1               NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
C
C     ******************************************************************
C
      WRITE (NOUT, 60)
      WRITE (NOUT, 10) SUBNAM
      WRITE (NOUT, 20) MESSAG
      WRITE (NOUT, 30)
      IF (LABEL1.NE.' ') WRITE (NOUT, 40) LABEL1,I
      IF (LABEL2.NE.' ') WRITE (NOUT, 40) LABEL2,J
      IF (LABEL3.NE.' ') WRITE (NOUT, 50) LABEL3,WORD
      WRITE (NOUT, 60)
C
      IF (ISTOP.EQ.0) RETURN
C
      CALL CLSFIL
C
      STOP 'ERROR'
C
   10 FORMAT (/,10X,' ERROR FOUND IN - ' ,A)
   20 FORMAT (/,10X,' DESCRIPTION - ' ,A)
   30 FORMAT (/,10X,' RELEVANT PARAMETERS - ')
   40 FORMAT (/,15X, A, ' = ' ,I10)
   50 FORMAT (/,15X, A, ' = ', A)
   60 FORMAT (/,10X,'* * * * * * * * * * * * * * * * * * * * * * * * * *
     1 * * * * * * * * * * * * * * * * * * * *  ',/,10X,'* * * * * * * *
     2 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
     3 * * * * *  ',/)
      END
