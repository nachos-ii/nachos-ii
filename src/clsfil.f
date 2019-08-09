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
      SUBROUTINE CLSFIL
C
C     ******************************************************************
C
C     SUBROUTINE TO CLOSE PREVIOUSLY OPENED DISK FILES
C
C     ******************************************************************
C
      COMMON /TAPES/  NIN,NOUT,NTP0,NTP1,NTP2,NTP3,NTP4,NTP5,NTP6,NTP7, 
     1                NTP8,NTP9,NTP10,NTP11,NTP12,NTP13
      COMMON /NTPDAT/ IFILES(16)
C
C     ******************************************************************
C
      IF (IFILES(1).EQ.1) CLOSE (UNIT=5, STATUS='keep')
      IF (IFILES(2).EQ.1) CLOSE (UNIT=6, STATUS='keep')
      IF (IFILES(3).EQ.1) CLOSE (UNIT=NTP0, STATUS='keep')
      IF (IFILES(4).EQ.1) CLOSE (UNIT=NTP1, STATUS='delete')
      IF (IFILES(5).EQ.1) CLOSE (UNIT=NTP2, STATUS='delete')
      IF (IFILES(6).EQ.1) CLOSE (UNIT=NTP3, STATUS='delete')
      IF (IFILES(7).EQ.1) CLOSE (UNIT=NTP4, STATUS='delete')
      IF (IFILES(8).EQ.1) CLOSE (UNIT=NTP5, STATUS='delete')
      IF (IFILES(9).EQ.1) CLOSE (UNIT=NTP6, STATUS='delete')
      IF (IFILES(10).EQ.1) CLOSE (UNIT=NTP7, STATUS='delete')
      IF (IFILES(11).EQ.1) CLOSE (UNIT=NTP8, STATUS='delete')
      IF (IFILES(12).EQ.1) CLOSE (UNIT=NTP9, STATUS='keep')
      IF (IFILES(13).EQ.1) CLOSE (UNIT=NTP10, STATUS='delete')
      IF (IFILES(14).EQ.1) CLOSE (UNIT=NTP11, STATUS='delete')
      IF (IFILES(15).EQ.1) CLOSE (UNIT=NTP12, STATUS='keep')
      IF (IFILES(16).EQ.1) CLOSE (UNIT=NTP13, STATUS='keep')
C
      RETURN
      END
