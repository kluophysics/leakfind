C*==outputs.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE OUTPUTS(OUTFILE,IBLOCH,NT,NP,NEW,IROTMAX,NPOL,MCD,
     &                   EPHOTO,PHOTOC,XESINT)
C     /****************************************************************/
C     # purpose       :  manage output of data into files
C     # notes         :
C     !!!! not yet completely finished !!!!
C
C         this is compressed output, detailed spectra are in files
C         fort.## if ip>=0, fort.4# contains 1st kind and fort.5#
C         the 2nd kind of photon polarization output (for mcd additional
C         files are 6#, or 7# for opposite magnetisation.
C         with open statements in the subroutines these files are
C         named xxx##.dat.
C         output should be improved to avoid overwriting of files
C         in succesive runs.
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:MEW,MPW
      USE MOD_SPEC_OUTPUT,ONLY:NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IBLOCH,IROTMAX,MCD,NEW,NP,NPOL,NT
      CHARACTER*20 OUTFILE
      REAL*8 EPHOTO(MEW),PHOTOC(MEW,MPW,MPW,2,2),XESINT(MEW,2)
C
C Local variables
C
      INTEGER JE,JP,JT
      REAL*8 PHOTO1,PHOTO2,PHOTOA,SPINA
C
C*** End of declarations rewritten by SPAG
C
      SELECT CASE (IBLOCH)
      CASE (0,1,2,4,7,8,9)
C
         SELECT CASE (NPOL)
         CASE (0)
C             unpolarized photons and p-s-dichroism
C             or magnetic dichroism
C           WRITE (NOUT1,99013)
            DO JT = 1,NT
               DO JP = 1,NP
                  DO JE = 1,NEW
                     PHOTO1 = PHOTOC(JE,JT,JP,1,1)
     &                        + PHOTOC(JE,JT,JP,2,1)
                     PHOTO2 = PHOTOC(JE,JT,JP,1,2)
     &                        + PHOTOC(JE,JT,JP,2,2)
                     PHOTOA = (PHOTO1+PHOTO2)/DBLE(2*IROTMAX)
                     SPINA = ((PHOTOC(JE,JT,JP,1,1)-PHOTOC(JE,JT,JP,2,1)
     &                       )
     &                       +(PHOTOC(JE,JT,JP,1,2)-PHOTOC(JE,JT,JP,2,2)
     &                       ))/DBLE(IROTMAX)
C
C                    WRITE (NOUT1,99010) EBIND(JE),JT,JP,PHOTOA,SPINA,
C    &                                PHOTOD
                  END DO
               END DO
            END DO
         CASE (3)
C             other orthogonal photon polarisation and dichroism
C             or magnetic dichroism
C           WRITE (NOUT1,99012)
            DO JT = 1,NT
               DO JP = 1,NP
                  DO JE = 1,NEW
                     PHOTO1 = (PHOTOC(JE,JT,JP,1,1)+PHOTOC(JE,JT,JP,2,1)
     &                        )/DBLE(2*IROTMAX)
C
                     PHOTO2 = (PHOTOC(JE,JT,JP,1,2)+PHOTOC(JE,JT,JP,2,2)
     &                        )/DBLE(2*IROTMAX)
C
                     PHOTOA = PHOTO1 + PHOTO2
                     SPINA = ((PHOTOC(JE,JT,JP,1,1)-PHOTOC(JE,JT,JP,2,1)
     &                       )
     &                       +(PHOTOC(JE,JT,JP,1,2)-PHOTOC(JE,JT,JP,2,2)
     &                       ))/DBLE(IROTMAX)
C
C                    WRITE (NOUT1,99009) EBIND(JE),JT,JP,PHOTO1,SPIN1,
C    &                                PHOTO2,SPIN2,PHOTOA,PHOTOD,SPINA
                  END DO
               END DO
            END DO
C
         CASE DEFAULT
C         single photon poalrization
            IF ( MCD.EQ.0 ) THEN
C             single photon polarisation, or spleed
C              WRITE (NOUT1,99011)
               DO JT = 1,NT
                  DO JP = 1,NP
                     DO JE = 1,NEW
                        PHOTOA = (PHOTOC(JE,JT,JP,1,1)
     &                           +PHOTOC(JE,JT,JP,2,1))/DBLE(2*IROTMAX)
                        SPINA = (PHOTOC(JE,JT,JP,1,1)
     &                          -PHOTOC(JE,JT,JP,2,1))/DBLE(IROTMAX)
C
C                       WRITE (NOUT1,99007) EBIND(JE),JT,JP,PHOTOA,SPINA
                     END DO
                  END DO
               END DO
            ELSE
C             mcd with single photon polarisation!
C              WRITE (NOUT1,99012)
               DO JT = 1,NT
                  DO JP = 1,NP
                     DO JE = 1,NEW
                        PHOTO1 = (PHOTOC(JE,JT,JP,1,1)
     &                           +PHOTOC(JE,JT,JP,2,1))/DBLE(2*IROTMAX)
C
                        PHOTO2 = (PHOTOC(JE,JT,JP,1,2)
     &                           +PHOTOC(JE,JT,JP,2,2))/DBLE(2*IROTMAX)
C
                        PHOTOA = PHOTO1 + PHOTO2
                        SPINA = ((PHOTOC(JE,JT,JP,1,1)-PHOTOC(JE,JT,JP,2
     &                          ,1))
     &                          +(PHOTOC(JE,JT,JP,1,2)-PHOTOC(JE,JT,JP,
     &                          2,2)))/DBLE(IROTMAX)
C
C                       WRITE (NOUT1,99009) EBIND(JE),JT,JP,PHOTO1,
C    &                         SPIN1,PHOTO2,SPIN2,PHOTOA,PHOTOD,SPINA
                     END DO
                  END DO
               END DO
            END IF
C
         END SELECT
C
      CASE (3)
         WRITE (*,99001) IBLOCH,OUTFILE
         WRITE (NOUT1,99001) IBLOCH,OUTFILE
C
      CASE (5)
         WRITE (NOUT1,99004)
         DO JE = 1,NEW
            PHOTOA = (XESINT(JE,1)+XESINT(JE,2))/DBLE(2*IROTMAX)
            SPINA = (XESINT(JE,1)-XESINT(JE,2))/DBLE(IROTMAX)
            WRITE (38,99003) EPHOTO(JE),XESINT(JE,1),XESINT(JE,2),
     &                       PHOTOA,SPINA
         END DO
         WRITE (*,99002) IBLOCH,OUTFILE
         WRITE (NOUT1,99002) IBLOCH,OUTFILE
C
      CASE (6)
         WRITE (NOUT1,99005)
         DO JE = 1,NEW
            PHOTOA = (XESINT(JE,1)+XESINT(JE,2))/DBLE(2*IROTMAX)
            SPINA = (XESINT(JE,1)-XESINT(JE,2))/DBLE(IROTMAX)
            WRITE (38,99003) EPHOTO(JE),XESINT(JE,1),XESINT(JE,2),
     &                       PHOTOA,SPINA
         END DO
         WRITE (*,99002) IBLOCH,OUTFILE
         WRITE (NOUT1,99002) IBLOCH,OUTFILE
C
      END SELECT
C
      RETURN
C
99001 FORMAT (1x,'ibloch =',i2,2x,'aes output written to :',a20)
99002 FORMAT (1x,'ibloch =',i2,2x,'xes, xas output written to :',a20)
C
99003 FORMAT (1x,e14.7,4(2x,e14.7))
C
99004 FORMAT (3x,'eph',4x,'ircp',3x,'ilcp',3x,'i0',5x,'polari')
99005 FORMAT (3x,'eph',4x,'ircp',3x,'ilcp',3x,'i0',5x,'dichro')
C
      END
C*==gcoutputs.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE GCOUTPUTS(OUTFILE,IBLOCH,IGCONV,NT,NP,NEW,IROTMAX,NPOL,
     &                     MCD,EBIND,EPHOTO,PHOTOC,XESINT,FTEMP,FWHM)
C     /****************************************************************/
C     # purpose       :  convolution by gaussian and fermi distribution
C                        manage output of data into files
C     # notes         :
C                         see also outputs
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:MEW,MPW
      USE MOD_SPEC_OUTPUT,ONLY:NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 FTEMP,FWHM
      INTEGER IBLOCH,IGCONV,IROTMAX,MCD,NEW,NP,NPOL,NT
      CHARACTER*20 OUTFILE
      REAL*8 EBIND(MEW),EPHOTO(MEW),PHOTOC(MEW,MPW,MPW,2,2),
     &       XESINT(MEW,2)
C
C Local variables
C
      INTEGER IS,JE,JP,JS,JT
      REAL*8 PHOTO1,PHOTO2,PHOTOA,SPINA,XI(:),YC(:),YI(:)
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE YC,XI,YI
      ALLOCATE (YC(MEW),XI(MEW),YI(MEW))
C
C     if there is no energy dependence, we have nothing to do
      IF ( NEW.LE.1 ) RETURN
C
      SELECT CASE (IBLOCH)
      CASE (0,2,4,7,8,9)
C        OPEN (36,FILE=OUTFILE,STATUS='unknown')
C
         SELECT CASE (NPOL)
         CASE (0)
C             unpolarized photons and p-s-dichroism
C             or magnetic dichroism
C           WRITE (36,99011)
            DO JT = 1,NT
               DO JP = 1,NP
C
                  DO IS = 1,2
                     DO JS = 1,2
                        DO JE = 1,NEW
                           XI(JE) = EBIND(JE)
                           YI(JE) = PHOTOC(JE,JT,JP,IS,JS)
                        END DO
                        CALL MAKECONV(XI,YI,NEW,IGCONV,FTEMP,FWHM,YC)
                        DO JE = 1,NEW
                           PHOTOC(JE,JT,JP,IS,JS) = YC(JE)
                        END DO
                     END DO
                  END DO
C
                  DO JE = 1,NEW
                     PHOTO1 = PHOTOC(JE,JT,JP,1,1)
     &                        + PHOTOC(JE,JT,JP,2,1)
                     PHOTO2 = PHOTOC(JE,JT,JP,1,2)
     &                        + PHOTOC(JE,JT,JP,2,2)
                     PHOTOA = (PHOTO1+PHOTO2)/DBLE(2*IROTMAX)
                     SPINA = ((PHOTOC(JE,JT,JP,1,1)-PHOTOC(JE,JT,JP,2,1)
     &                       )
     &                       +(PHOTOC(JE,JT,JP,1,2)-PHOTOC(JE,JT,JP,2,2)
     &                       ))/DBLE(IROTMAX)
C
C                    WRITE (36,99008) EBIND(JE),JT,JP,PHOTOA,SPINA,
C    &                                PHOTOD
                  END DO
               END DO
            END DO
         CASE (3)
C             other orthogonal photon polarisation and dichroism
C             or magnetic dichroism
C           WRITE (36,99010)
            DO JT = 1,NT
               DO JP = 1,NP
C
                  DO IS = 1,2
                     DO JS = 1,2
                        DO JE = 1,NEW
                           XI(JE) = EBIND(JE)
                           YI(JE) = PHOTOC(JE,JT,JP,IS,JS)
                        END DO
                        CALL MAKECONV(XI,YI,NEW,IGCONV,FTEMP,FWHM,YC)
                        DO JE = 1,NEW
                           PHOTOC(JE,JT,JP,IS,JS) = YC(JE)
                        END DO
                     END DO
                  END DO
C
                  DO JE = 1,NEW
                     PHOTO1 = (PHOTOC(JE,JT,JP,1,1)+PHOTOC(JE,JT,JP,2,1)
     &                        )/DBLE(2*IROTMAX)
C
                     PHOTO2 = (PHOTOC(JE,JT,JP,1,2)+PHOTOC(JE,JT,JP,2,2)
     &                        )/DBLE(2*IROTMAX)
C
                     PHOTOA = PHOTO1 + PHOTO2
                     SPINA = ((PHOTOC(JE,JT,JP,1,1)-PHOTOC(JE,JT,JP,2,1)
     &                       )
     &                       +(PHOTOC(JE,JT,JP,1,2)-PHOTOC(JE,JT,JP,2,2)
     &                       ))/DBLE(IROTMAX)
C
C                    WRITE (36,99007) EBIND(JE),JT,JP,PHOTO1,SPIN1,
C    &                                PHOTO2,SPIN2,PHOTOA,PHOTOD,SPINA
                  END DO
               END DO
            END DO
C
         CASE DEFAULT
            IF ( MCD.EQ.0 ) THEN
C             single photon polarisation
C              WRITE (36,99009)
               DO JT = 1,NT
                  DO JP = 1,NP
C
                     DO IS = 1,2
                        DO JE = 1,NEW
                           XI(JE) = EBIND(JE)
                           YI(JE) = PHOTOC(JE,JT,JP,IS,1)
                        END DO
                        CALL MAKECONV(XI,YI,NEW,IGCONV,FTEMP,FWHM,YC)
                        DO JE = 1,NEW
                           PHOTOC(JE,JT,JP,IS,1) = YC(JE)
                        END DO
                     END DO
C
                     DO JE = 1,NEW
                        PHOTOA = (PHOTOC(JE,JT,JP,1,1)
     &                           +PHOTOC(JE,JT,JP,2,1))/DBLE(2*IROTMAX)
                        SPINA = (PHOTOC(JE,JT,JP,1,1)
     &                          -PHOTOC(JE,JT,JP,2,1))/DBLE(IROTMAX)
C
C                       WRITE (36,99005) EBIND(JE),JT,JP,PHOTOA,SPINA
                     END DO
                  END DO
               END DO
            ELSE
C             mcd with single photon polarisation!
C              WRITE (36,99010)
               DO JT = 1,NT
                  DO JP = 1,NP
C
                     DO IS = 1,2
                        DO JS = 1,2
                           DO JE = 1,NEW
                              XI(JE) = EBIND(JE)
                              YI(JE) = PHOTOC(JE,JT,JP,IS,JS)
                           END DO
                           CALL MAKECONV(XI,YI,NEW,IGCONV,FTEMP,FWHM,YC)
                           DO JE = 1,NEW
                              PHOTOC(JE,JT,JP,IS,JS) = YC(JE)
                           END DO
                        END DO
                     END DO
C
                     DO JE = 1,NEW
                        PHOTO1 = (PHOTOC(JE,JT,JP,1,1)
     &                           +PHOTOC(JE,JT,JP,2,1))/DBLE(2*IROTMAX)
C
                        PHOTO2 = (PHOTOC(JE,JT,JP,1,2)
     &                           +PHOTOC(JE,JT,JP,2,2))/DBLE(2*IROTMAX)
C
                        PHOTOA = PHOTO1 + PHOTO2
                        SPINA = ((PHOTOC(JE,JT,JP,1,1)-PHOTOC(JE,JT,JP,2
     &                          ,1))
     &                          +(PHOTOC(JE,JT,JP,1,2)-PHOTOC(JE,JT,JP,
     &                          2,2)))/DBLE(IROTMAX)
C
C                       WRITE (36,99007) EBIND(JE),JT,JP,PHOTO1,SPIN1,
C    &                         PHOTO2,SPIN2,PHOTOA,PHOTOD,SPINA
                     END DO
                  END DO
               END DO
            END IF
C
         END SELECT
C        CLOSE (36)
         WRITE (*,99001) IBLOCH,OUTFILE
         WRITE (NOUT1,99001) IBLOCH,OUTFILE
C
      CASE (3)
         WRITE (*,99002) IBLOCH,OUTFILE
         WRITE (NOUT1,99002) IBLOCH,OUTFILE
C
      CASE (5)
         OPEN (38,FILE=OUTFILE,STATUS='unknown')
         WRITE (38,99006)
         DO IS = 1,2
            DO JE = 1,NEW
               XI(JE) = EPHOTO(JE)
               YI(JE) = XESINT(JE,IS)
            END DO
            CALL MAKECONV(XI,YI,NEW,IGCONV,FTEMP,FWHM,YC)
            DO JE = 1,NEW
               XESINT(JE,IS) = YC(JE)
            END DO
         END DO
         DO JE = 1,NEW
            PHOTOA = (XESINT(JE,1)+XESINT(JE,2))/DBLE(2*IROTMAX)
            SPINA = (XESINT(JE,1)-XESINT(JE,2))/DBLE(IROTMAX)
C
            WRITE (38,99005) EPHOTO(JE),XESINT(JE,1),XESINT(JE,2),
     &                       PHOTOA,SPINA
         END DO
         CLOSE (38)
         WRITE (*,99003) IBLOCH,OUTFILE
         WRITE (NOUT1,99003) IBLOCH,OUTFILE
C
      CASE (6)
         OPEN (38,FILE=OUTFILE,STATUS='unknown')
         WRITE (38,99007)
         DO IS = 1,2
            DO JE = 1,NEW
               XI(JE) = EPHOTO(JE)
               YI(JE) = XESINT(JE,IS)
            END DO
            CALL MAKECONV(XI,YI,NEW,IGCONV,FTEMP,FWHM,YC)
            DO JE = 1,NEW
               XESINT(JE,IS) = YC(JE)
            END DO
         END DO
C
         DO JE = 1,NEW
            PHOTOA = (XESINT(JE,1)+XESINT(JE,2))/DBLE(2*IROTMAX)
            SPINA = (XESINT(JE,1)-XESINT(JE,2))/DBLE(IROTMAX)
            WRITE (38,99005) EPHOTO(JE),XESINT(JE,1),XESINT(JE,2),
     &                       PHOTOA,SPINA
         END DO
         CLOSE (38)
         WRITE (*,99003) IBLOCH,OUTFILE
         WRITE (NOUT1,99003) IBLOCH,OUTFILE
C
      CASE DEFAULT
         WRITE (*,99004) IBLOCH
         WRITE (NOUT1,99004) IBLOCH
      END SELECT
C
      RETURN
C
99001 FORMAT (1x,'ibloch =',i2,2x,'output written to: ',a20)
99002 FORMAT (1x,'ibloch =',i2,2x,'aes output written to :',a20)
99003 FORMAT (1x,'ibloch =',i2,2x,'xes, xas output written to :',a20)
99004 FORMAT (1x,'ibloch =',i2,2x,'convolution not defined')
C
99005 FORMAT (1x,e14.7,4(2x,e14.7))
C
99006 FORMAT (3x,'eph',4x,'ircp',3x,'ilcp',3x,'i0',5x,'polari')
99007 FORMAT (3x,'eph',4x,'ircp',3x,'ilcp',3x,'i0',5x,'dichro')
C
      END
C*==sortbands.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE SORTBANDS(OUTFILE,FOUTFILE)
C     /****************************************************************/
C     # purpose       :  sort bands output file                        *
C     /****************************************************************/
C
      USE MOD_SPEC_OUTPUT,ONLY:NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER JMAX
      PARAMETER (JMAX=9000)
C
C Dummy arguments
C
      CHARACTER*20 FOUTFILE,OUTFILE
C
C Local variables
C
      REAL EB(JMAX),EBO(JMAX),ER(JMAX),ERO(JMAX),KR(JMAX),KRO(JMAX),X1,
     &     X2,X3
      INTEGER I,JN,NN
C
C*** End of declarations rewritten by SPAG
C
C      CHARACTER*20 DUMMY
C
C*** End of declarations rewritten by SPAG
C
      OPEN (36,FILE=OUTFILE,STATUS='unknown')
C      READ (36,*) DUMMY
      JN = 1
      READ (36,*) KR(JN),ER(JN),EB(JN)
      DO I = 1,JMAX
         READ (36,*) X1,X2,X3
C         IF ( X1.NE.KR(JN) .OR. X2.NE.ER(JN) ) THEN
         IF ( ABS(X1-KR(JN)).GT.1.0D-16 .OR. ABS(X2-ER(JN)).GT.1.0D-16 )
     &        THEN
            JN = JN + 1
            KR(JN) = X1
            ER(JN) = X2
            EB(JN) = X3
         END IF
         IF ( JN.EQ.JMAX ) THEN
            WRITE (*,99002)
            WRITE (*,99003) JMAX
            WRITE (NOUT1,99002)
            WRITE (NOUT1,99003) JMAX
            EXIT
         END IF
      END DO
      CLOSE (36)
C
C     sort for k-values
      NN = JN
      CALL PIKSORT3(NN,KR,ER,EB)
C
      JN = 1
      KRO(JN) = KR(JN)
      ERO(JN) = ER(JN)
      EBO(JN) = EB(JN)
      DO I = 1,NN
C         IF ( KRO(JN).NE.KR(I) .OR. ERO(JN).NE.ER(I) ) THEN
         IF ( ABS(KRO(JN)-KR(I)).GT.1.0D-16 .OR. ABS(ERO(JN)-ER(I))
     &        .GT.1.0D-16 ) THEN
            JN = JN + 1
            KRO(JN) = KR(I)
            ERO(JN) = ER(I)
            EBO(JN) = EB(I)
         END IF
      END DO
C
C     sort for e values
      NN = JN
      DO I = 1,JN
         KR(I) = KRO(I)
         ER(I) = ERO(I)
         EB(I) = EBO(I)
      END DO
      CALL PIKSORT3(NN,ER,KR,EB)
C
      JN = 1
      KRO(JN) = KR(JN)
      ERO(JN) = ER(JN)
      EBO(JN) = EB(JN)
      DO I = 1,NN
C         IF ( KRO(JN).NE.KR(I) .OR. ERO(JN).NE.ER(I) ) THEN
         IF ( ABS(KRO(JN)-KR(I)).GT.1.0D-16 .OR. ABS(ERO(JN)-ER(I))
     &        .GT.1.0D-16 ) THEN
            JN = JN + 1
            KRO(JN) = KR(I)
            ERO(JN) = ER(I)
            EBO(JN) = EB(I)
         END IF
      END DO
C
      NN = JN
C
C     OPEN (37,FILE=FOUTFILE,STATUS='unknown')
C     WRITE (37,99001)
C     DO JN = 1,NN
C        WRITE (37,99002) KRO(JN),ERO(JN),EBO(JN)
C     END DO
C     CLOSE (37)
C
      WRITE (*,99001) FOUTFILE
      WRITE (NOUT1,99001) FOUTFILE
C
      RETURN
C
99001 FORMAT (1x,'ibloch = 1',2x,'output written to: ',a20)
99002 FORMAT (1x,'warning: too many data in bands.dat')
99003 FORMAT (1x,'data truncated to n=',i5)
C
      END
C*==sortcbands.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE SORTCBANDS(OUTFILE,FOUTFILE)
C     /****************************************************************/
C     # purpose       :  sort complex bands output file                *
C     /****************************************************************/
C
      USE MOD_SPEC_OUTPUT,ONLY:NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER JMAX
      PARAMETER (JMAX=9000)
C
C Dummy arguments
C
      CHARACTER*20 FOUTFILE,OUTFILE
C
C Local variables
C
      REAL EB(JMAX),EBO(JMAX),ER(JMAX),ERO(JMAX),KI(JMAX),KIO(JMAX),
     &     KR(JMAX),KRO(JMAX),X11,X12,X2,X3
      INTEGER I,JN,NN
C
C*** End of declarations rewritten by SPAG
C
C      CHARACTER*20 DUMMY
C
C*** End of declarations rewritten by SPAG
C
      OPEN (36,FILE=OUTFILE,STATUS='unknown')
C      READ (36,*) DUMMY
      JN = 1
      READ (36,*) KR(JN),KI(JN),ER(JN),EB(JN)
      DO I = 1,JMAX
         READ (36,*) X11,X12,X2,X3
C         IF ( X11.NE.KR(JN) .OR. X2.NE.ER(JN) ) THEN
         IF ( ABS(X11-KR(JN)).GT.1.0D-16 .OR. ABS(X2-ER(JN))
     &        .GT.1.0D-16 ) THEN
            JN = JN + 1
            KR(JN) = X11
            KI(JN) = X12
            ER(JN) = X2
            EB(JN) = X3
         END IF
         IF ( JN.EQ.JMAX ) THEN
            WRITE (*,99002)
            WRITE (*,99003) JMAX
            WRITE (NOUT1,99002)
            WRITE (NOUT1,99003) JMAX
            EXIT
         END IF
      END DO
      CLOSE (36)
C
C     sort for k-values
      NN = JN
      CALL PIKSORT4(NN,KR,KI,ER,EB)
C
      JN = 1
      KRO(JN) = KR(JN)
      KIO(JN) = KR(JN)
      ERO(JN) = ER(JN)
      EBO(JN) = EB(JN)
      DO I = 1,NN
C         IF ( KRO(JN).NE.KR(I) .OR. ERO(JN).NE.ER(I) ) THEN
         IF ( ABS(KRO(JN)-KR(I)).GT.1.0D-16 .OR. ABS(ERO(JN)-ER(I))
     &        .GT.1.0D-16 ) THEN
            JN = JN + 1
            KRO(JN) = KR(I)
            KIO(JN) = KI(I)
            ERO(JN) = ER(I)
            EBO(JN) = EB(I)
         END IF
      END DO
C
C     sort for e values
      NN = JN
      DO I = 1,JN
         KR(I) = KRO(I)
         KI(I) = KIO(I)
         ER(I) = ERO(I)
         EB(I) = EBO(I)
      END DO
      CALL PIKSORT4(NN,ER,KR,KI,EB)
C
      JN = 1
      KRO(JN) = KR(JN)
      KIO(JN) = KI(JN)
      ERO(JN) = ER(JN)
      EBO(JN) = EB(JN)
      DO I = 1,NN
C         IF ( KRO(JN).NE.KR(I) .OR. ERO(JN).NE.ER(I) ) THEN
         IF ( ABS(KRO(JN)-KR(I)).GT.1.0D-16 .OR. ABS(ERO(JN)-ER(I))
     &        .GT.1.0D-16 ) THEN
            JN = JN + 1
            KRO(JN) = KR(I)
            KIO(JN) = KI(I)
            ERO(JN) = ER(I)
            EBO(JN) = EB(I)
         END IF
      END DO
C
      NN = JN
C
C     OPEN (37,FILE=FOUTFILE,STATUS='unknown')
C
C     DO JN = 1,NN
C        WRITE (37,99002) KRO(JN),KIO(JN),ERO(JN),EBO(JN)
C     END DO
C     CLOSE (37)
C
      WRITE (*,99001) FOUTFILE
      WRITE (NOUT1,99001) FOUTFILE
C
      RETURN
C
99001 FORMAT (1x,'ibloch = 1',2x,'output written to: ',a20)
99002 FORMAT (1x,'warning: too many data in bands.dat')
99003 FORMAT (1x,'data truncated to n=',i5)
C
      END
C*==makeconv.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE MAKECONV(XI,YI,NEW,IGCONV,FTEMP,FWHM,YOUT)
C     /****************************************************************/
C     # purpose      : driver for convolution of a data set
C     # parameters:
C         xi      =>  abscissa
C         yi      =>  function values
C         new     =>  number of function values
C         igconv  =>  kind of convolution
C         fwhm    |
C         ftemp   |=>  parameter for convolution
C         fwhm    |
C
C         yout    =>  convoluted yi
C
C     # note:  works only for equidistant data
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:MEW
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 FTEMP,FWHM
      INTEGER IGCONV,NEW
      REAL*8 XI(MEW),YI(MEW),YOUT(MEW)
C
C Local variables
C
      INTEGER I,IOCC
C
C*** End of declarations rewritten by SPAG
C
C                 iocc => kind of fermi-dirac distribution
C                         0: unoccupied;  1: occupied
C
      SELECT CASE (IGCONV)
      CASE (1)
         CALL GAUSCONV(XI,YI,NEW,FWHM,YOUT)
      CASE (2)
         IOCC = 1
         CALL FERMICONV(XI,YI,NEW,FTEMP,IOCC,YOUT)
      CASE (3)
         IOCC = 0
         CALL FERMICONV(XI,YI,NEW,FTEMP,IOCC,YOUT)
      CASE (4)
         IOCC = 1
         CALL FERMIGAUSCONV(XI,YI,NEW,FTEMP,IOCC,FWHM,YOUT)
      CASE (5)
         IOCC = 0
         CALL FERMIGAUSCONV(XI,YI,NEW,FTEMP,IOCC,FWHM,YOUT)
      CASE DEFAULT
         DO I = 1,NEW
            YOUT(I) = YI(I)
         END DO
      END SELECT
C
      END
C*==gausconv.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE GAUSCONV(XI,YI,N,FWHM,YOUT)
C     /****************************************************************/
C     # purpose      : convolution of a data set with
C                      gaussian
C     # parameters:
C         xi      =>  abscissa
C         yi      =>  function values
C         n       =>  number of function values
C         fwhm    =>  full width at half maximum for gaussian
C         yout    =>  yi convoluted with gaussian
C
C     # note:  works only for equidistant data
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:MEW,PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 FWHM
      INTEGER N
      REAL*8 XI(MEW),YI(MEW),YOUT(MEW)
C
C Local variables
C
      REAL*8 A,ARG,CONV(N),CSUM,DH,EM,GAUSS,H,WPI,X0
      INTEGER DK,I,IMAX,K,KMAX,KMIN
C
C*** End of declarations rewritten by SPAG
C
      H = ABS(XI(1)-XI(2))
      WPI = SQRT(PI)
      A = 2.D0*SQRT(LOG(2.D0))/FWHM
C
      IMAX = IDINT(15.*FWHM/H)
      IF ( IMAX.LE.0 ) IMAX = 1
C
C     convolution of complex dos with gaussian
C     ========================================
      DO I = 1,N
         X0 = XI(I)
         KMIN = MAX(1,I-IMAX)
         KMAX = MIN(N,I+IMAX)
         DO K = KMIN,KMAX
            EM = A*(XI(K)-X0)
            ARG = EM**2
            IF ( ARG.GT.50.D0 ) THEN
               GAUSS = EXP(-50.D0)
            ELSE
               GAUSS = (A/WPI)*EXP(-ARG)
            END IF
            CONV(K) = YI(K)*GAUSS
         END DO
C
C        integrate just by summation (may be improved later)
C        ===================================================
         CSUM = 0.D0
         DK = MIN(1,ABS(KMIN-KMAX))
         DH = H/DK
         DO K = KMIN,KMAX
            CSUM = CSUM + CONV(K)
         END DO
         YOUT(I) = CSUM*DH
      END DO
C
      END
C*==fermigausconv.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE FERMIGAUSCONV(XI,YI,N,FTEMP,IOCC,FWHM,YOUT)
C     /****************************************************************/
C     # purpose      : convolution of a data set with
C                      fermi-dirac distribution and gaussian
C     # parameters:
C         xi      =>  abscissa
C         yi      =>  function values
C         n       =>  number of function values
C         ftemp   =>  temperature in k
C         iocc    =>  kind of fermi distribution
C                     =0  unoccupied (e.g. ipe)
C                     =1  occupied (e.g. ups)
C         fwhm    =>  full width at half maximum for gaussian
C         yout    =>  yi multiplied by fermi-dirac
C                        and convoluted with gaussian
C
C     # note:  works only for equidistant data
C              fermi energy is assumed to be at zero
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:MEW,PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 FTEMP,FWHM
      INTEGER IOCC,N
      REAL*8 XI(MEW),YI(MEW),YOUT(MEW)
C
C Local variables
C
      REAL*8 A,ARG,CONV(:),CSUM,DH,EM,GAUSS,H,WPI,X0,YTEMP(:)
      INTEGER DK,I,IMAX,K,KMAX,KMIN
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE CONV,YTEMP
      ALLOCATE (CONV(N),YTEMP(N))
C
C*** End of declarations rewritten by SPAG
C
C      kbolz=0.862d-4
C      ktemp = ftemp * kbolz
C
C     multiplication with fermi-function:
C     ===================================
      CALL FERMICONV(XI,YI,N,FTEMP,IOCC,YTEMP)
C
C     convolution with gaussian:
C     ==========================
C
      H = ABS(XI(1)-XI(2))
      WPI = SQRT(PI)
      A = 2.D0*SQRT(LOG(2.D0))/FWHM
C
      IMAX = IDINT(15.*FWHM/H)
      IF ( IMAX.LE.0 ) IMAX = 1
C
      DO I = 1,N
         X0 = XI(I)
         KMIN = MAX(1,I-IMAX)
         KMAX = MIN(N,I+IMAX)
         DO K = KMIN,KMAX
            EM = A*(XI(K)-X0)
            ARG = EM**2
            IF ( ARG.GT.50.D0 ) THEN
               GAUSS = EXP(-50.D0)
            ELSE
               GAUSS = (A/WPI)*EXP(-ARG)
            END IF
            CONV(K) = YTEMP(K)*GAUSS
         END DO
C
C        integrate just by summation (may be improved later)
C        ===================================================
         CSUM = 0.D0
         DK = MIN(1,ABS(KMIN-KMAX))
         DH = H/DK
         DO K = KMIN,KMAX
            CSUM = CSUM + CONV(K)
         END DO
         YOUT(I) = CSUM*DH
      END DO
C
      END
C*==fermiconv.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE FERMICONV(XI,YI,N,FTEMP,IOCC,YOUT)
C     /****************************************************************/
C     # purpose      : multiplication of a data set with
C                      fermi-dirac distribution
C     # parameters:
C         xi      =>  abscissa
C         yi      =>  function values
C         n       =>  number of function values
C         ftemp   =>  temperature in k
C         iocc    =>  kind of fermi distribution
C                     =0  unoccupied (e.g. ipe)
C                     =1  occupied (e.g. ups)
C         yout    =>  yi multiplied by fermi-dirac
C
C     # note: fermi energy is assumed to be at zero
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:MEW
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 FTEMP
      INTEGER IOCC,N
      REAL*8 XI(MEW),YI(MEW),YOUT(MEW)
C
C Local variables
C
      INTEGER I
      REAL*8 KBOLZ,KTEMP,XF,YF
C
C*** End of declarations rewritten by SPAG
C
      KBOLZ = 0.862D-4
      KTEMP = FTEMP*KBOLZ
C
C     multiplication with fermi-function; ef=0:
C     =========================================
C     &: occupied   :  diracfermi = 1 /(1 + exp(e/kt))
C     0: unoccupied :  diracfermi = 1 /(1 + exp(-e/kt))
      DO I = 1,N
         XF = XI(I)/KTEMP
         YF = YI(I)
         IF ( IOCC.EQ.0 ) THEN
C         ipe case (unoccupied dos)
            IF ( XF.GT.50.0D0 ) THEN
               YOUT(I) = YF
            ELSE IF ( XF.LT.-50.0D0 ) THEN
               YOUT(I) = 0.D0
            ELSE
               YOUT(I) = YF/(1.D0+EXP(-XF))
            END IF
C         ups case (occupied dos)
         ELSE IF ( XF.GT.50.0D0 ) THEN
            YOUT(I) = 0.D0
         ELSE IF ( XF.LT.-50.0D0 ) THEN
            YOUT(I) = YF
         ELSE
            YOUT(I) = YF/(1.D0+EXP(XF))
         END IF
      END DO
C
      END
C*==piksort3.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE PIKSORT3(N,ARR,ARR2,ARR3)
C     /****************************************************************/
C     # purpose      : sort 3 vectors such that the first has          *
C                      ascending values                                *
C     /****************************************************************/
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N
      REAL ARR(N),ARR2(N),ARR3(N)
C
C Local variables
C
      REAL A,B,C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO J = 2,N
         A = ARR(J)
         B = ARR2(J)
         C = ARR3(J)
         DO I = J - 1,1, - 1
            IF ( ARR(I).LE.A ) GOTO 50
            ARR(I+1) = ARR(I)
            ARR2(I+1) = ARR2(I)
            ARR3(I+1) = ARR3(I)
         END DO
         I = 0
 50      CONTINUE
         ARR(I+1) = A
         ARR2(I+1) = B
         ARR3(I+1) = C
      END DO
C
      END
C*==piksort4.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE PIKSORT4(N,ARR,ARR2,ARR3,ARR4)
C     /****************************************************************/
C     # purpose      : sort 3 vectors such that the first has          *
C                      ascending values                                *
C     /****************************************************************/
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N
      REAL ARR(N),ARR2(N),ARR3(N),ARR4(N)
C
C Local variables
C
      REAL A,B,C,D
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO J = 2,N
         A = ARR(J)
         B = ARR2(J)
         C = ARR3(J)
         D = ARR4(J)
         DO I = J - 1,1, - 1
            IF ( ARR(I).LE.A ) GOTO 50
            ARR(I+1) = ARR(I)
            ARR2(I+1) = ARR2(I)
            ARR3(I+1) = ARR3(I)
            ARR4(I+1) = ARR4(I)
         END DO
         I = 0
 50      CONTINUE
         ARR(I+1) = A
         ARR2(I+1) = B
         ARR3(I+1) = C
         ARR4(I+1) = D
      END DO
C
      END
