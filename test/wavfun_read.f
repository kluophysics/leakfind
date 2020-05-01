C*==wavfun_write_rel.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      SUBROUTINE WAVFUN_WRITE_REL(IFIL,IT,KIRR,ZG,ZF,JG,JF,IRTOP,NCPLWF,
     &                            IKMCPLWF)
C   ********************************************************************
C   *                                                                  *
C   *  WRITE  the full vector of wave functions solutions for type IT  *
C   *                                                                  *
C   *      ZG(r,K,ISOL)  ISOL = 1,..,NKM_T  lin. indep. solutions      *
C   *                      K  = 1,..,NCPLWF coupled partial waves      *
C   *                                                                  *
C   *  if KIRR=1  read regular AND irregular wave fuinctions           *
C   *                                                                  *
C   *  the wave functions are written by:                              *
C   *     - <FPSSITE>  for FULLPOT calculations                        *
C   *     - <SSITE>    in case of ASA                                  *
C   *                                                                  *
C   *  fully relativistic version                                      *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX
      USE MOD_TYPES,ONLY:NT,IMT,NCPLWFMAX,NKM_T
      USE MOD_RMESH,ONLY:NRMAX,FULLPOT,JRWS,JRCRI
      IMPLICIT NONE
C*--WAVFUN_WRITE_REL25
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='WAVFUN_WRITE_REL')
C
C Dummy arguments
C
      INTEGER IFIL,IRTOP,IT,KIRR
      INTEGER IKMCPLWF(NCPLWFMAX,NKMMAX),NCPLWF(NKMMAX)
      COMPLEX*16 JF(NRMAX,NCPLWFMAX,NKMMAX),JG(NRMAX,NCPLWFMAX,NKMMAX),
     &           ZF(NRMAX,NCPLWFMAX,NKMMAX),ZG(NRMAX,NCPLWFMAX,NKMMAX)
C
C Local variables
C
      INTEGER IM,IR,ISOL,K
C
C*** End of declarations rewritten by SPAG
C
      IM = IMT(IT)
C
C========================================================================
C                    check consistency of parameters
C========================================================================
      IF ( FULLPOT ) THEN
         IF ( IRTOP.NE.JRCRI(IM) )
     &         CALL STOP_MESSAGE(ROUTINE,'for FP: IRTOP .NE. JRCRI(IM)')
      ELSE
         IF ( NCPLWFMAX.NE.2 )
     &         CALL STOP_MESSAGE(ROUTINE,'for ASA: NCPLWFMAX .NE. 2')
         IF ( IRTOP.NE.JRWS(IM) )
     &         CALL STOP_MESSAGE(ROUTINE,'for ASA: IRTOP .NE. JRWS(IM)')
      END IF
C========================================================================
C
C
C========================================================================
C                       read wave functions
C========================================================================
C
      DO ISOL = 1,NKM_T(IT)
C
         WRITE (IFIL,REC=ISOL+(IT-1)*NKM,ERR=100) IT,'REG',ISOL,IRTOP,
     &          NCPLWF(ISOL),(IKMCPLWF(K,ISOL),K=1,NCPLWF(ISOL)),
     &          ((ZG(IR,K,ISOL),IR=1,IRTOP),K=1,NCPLWF(ISOL)),
     &          ((ZF(IR,K,ISOL),IR=1,IRTOP),K=1,NCPLWF(ISOL))
C
C
         IF ( KIRR.EQ.1 ) WRITE (IFIL,REC=ISOL+(IT-1+NT)*NKM,ERR=100)
     &                           IT,'IRR',ISOL,IRTOP,
     &                           ((JG(IR,K,ISOL),IR=1,IRTOP),K=1,
     &                           NCPLWF(ISOL)),
     &                           ((JF(IR,K,ISOL),IR=1,IRTOP),K=1,
     &                           NCPLWF(ISOL))
C
      END DO
C========================================================================
C
      RETURN
C
 100  CONTINUE
      CALL STOP_MESSAGE(ROUTINE,'TROUBLE reading WF - ERROR occured')
      END
C*==wavfun_read_rel.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      SUBROUTINE WAVFUN_READ_REL(IFIL,IT,KIRR,ZG,ZF,JG,JF,IRTOP,NCPLWF,
     &                           IKMCPLWF)
C   ********************************************************************
C   *                                                                  *
C   *  reread the full vector of wave functions solutions for type IT  *
C   *                                                                  *
C   *      ZG(r,K,ISOL)  ISOL = 1,..,NKM_T  lin. indep. solutions      *
C   *                      K  = 1,..,NCPLWF coupled partial waves      *
C   *                                                                  *
C   *  if KIRR=1  read regular AND irregular wave fuinctions           *
C   *                                                                  *
C   *  the wave functions are written by:                              *
C   *     - <FPSSITE>  for FULLPOT calculations                        *
C   *     - <SSITE>    in case of ASA                                  *
C   *                                                                  *
C   *  fully relativistic version                                      *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IFILCORWF,IFILGFWF,IFILLDAU
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX
      USE MOD_TYPES,ONLY:NT,IMT,NCPLWFMAX,NKM_T
      USE MOD_RMESH,ONLY:NRMAX,FULLPOT,JRWS,JRCRI
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--WAVFUN_READ_REL131
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='WAVFUN_READ_REL')
C
C Dummy arguments
C
      INTEGER IFIL,IRTOP,IT,KIRR
      INTEGER IKMCPLWF(NCPLWFMAX,NKMMAX),NCPLWF(NKMMAX)
      COMPLEX*16 JF(NRMAX,NCPLWFMAX,NKMMAX),JG(NRMAX,NCPLWFMAX,NKMMAX),
     &           ZF(NRMAX,NCPLWFMAX,NKMMAX),ZG(NRMAX,NCPLWFMAX,NKMMAX)
C
C Local variables
C
      INTEGER IM,IR,IRTOPP,ISOL,ISOLP,ITP,K
      CHARACTER*3 STRP
C
C*** End of declarations rewritten by SPAG
C
      IM = IMT(IT)
C
C========================================================================
C                    check consistency of parameters
C========================================================================
      IF ( FULLPOT ) THEN
         IF ( IRTOP.NE.JRCRI(IM) )
     &         CALL STOP_MESSAGE(ROUTINE,'for FP: IRTOP .NE. JRCRI(IM)')
      ELSE
         IF ( NCPLWFMAX.NE.2 )
     &         CALL STOP_MESSAGE(ROUTINE,'for ASA: NCPLWFMAX .NE. 2')
         IF ( IRTOP.NE.JRWS(IM) )
     &         CALL STOP_MESSAGE(ROUTINE,'for ASA: IRTOP .NE. JRWS(IM)')
      END IF
C========================================================================
C
      NCPLWF(:) = 0
C
C========================================================================
C                       read wave functions
C========================================================================
C
      DO ISOL = 1,NKM_T(IT)
C
         READ (IFIL,REC=ISOL+(IT-1)*NKM,ERR=100) ITP,STRP,ISOLP,IRTOPP,
     &         NCPLWF(ISOL),(IKMCPLWF(K,ISOL),K=1,NCPLWF(ISOL)),
     &         ((ZG(IR,K,ISOL),IR=1,IRTOP),K=1,NCPLWF(ISOL)),
     &         ((ZF(IR,K,ISOL),IR=1,IRTOP),K=1,NCPLWF(ISOL))
C
         IF ( (STRP.NE.'REG' .AND. STRP.NE.'COR') .OR. IT.NE.ITP .OR. 
     &        IRTOP.NE.IRTOPP ) THEN
            WRITE (6,*) 'IT STRP ISOL',IT,ITP,'REG',STRP,ISOL,ISOLP
            CALL STOP_MESSAGE(ROUTINE,'TROUBLE reading REG WF')
         END IF
C
         IF ( KIRR.EQ.1 ) THEN
            IF ( (IFIL.EQ.IFILCORWF) .OR. (IFIL.EQ.IFILGFWF) .OR. 
     &           (IFIL.EQ.IFILLDAU) ) THEN
C
               JG(:,:,ISOL) = C0
               JF(:,:,ISOL) = C0
C
            ELSE
C
               READ (IFIL,REC=ISOL+(IT-1+NT)*NKM,ERR=100) ITP,STRP,
     &               ISOLP,IRTOPP,
     &               ((JG(IR,K,ISOL),IR=1,IRTOP),K=1,NCPLWF(ISOL)),
     &               ((JF(IR,K,ISOL),IR=1,IRTOP),K=1,NCPLWF(ISOL))
C
               IF ( STRP.NE.'IRR' .OR. IT.NE.ITP .OR. IRTOP.NE.IRTOPP )
     &              THEN
                  WRITE (6,*) 'IT STRP ISOL',IT,ITP,'IRR',STRP,ISOL,
     &                        ISOLP
                  CALL STOP_MESSAGE(ROUTINE,'TROUBLE reading IRR WF')
               END IF
C
            END IF
         END IF
C
      END DO
C========================================================================
C
      RETURN
C
 100  CONTINUE
      CALL STOP_MESSAGE(ROUTINE,'TROUBLE reading WF - ERROR occured')
      END
C*==wavfun_write_sra.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      SUBROUTINE WAVFUN_WRITE_SRA(IFIL,IT,KIRR,ZG,JG,IRTOP,NCPLWF,
     &                            IKMCPLWF)
C   ********************************************************************
C   *                                                                  *
C   *  WRITE  the full vector of wave functions solutions for type IT  *
C   *                                                                  *
C   *                                                                  *
C   *      ZG(r,K,ISOL)  ISOL = 1,..,NKM_T  lin. indep. solutions      *
C   *                      K  = 1,..,NCPLWF coupled partial waves      *
C   *                                                                  *
C   *  if KIRR=1  read regular AND irregular wave fuinctions           *
C   *                                                                  *
C   *  the wave functions are written by:                              *
C   *     - <FPNRSSITE>  for FULLPOT calculations                      *
C   *     - <NRSSITE>    for ASA     calculations                      *
C   *                                                                  *
C   *  scalar relativistic version                                     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,NL,NSPIN
      USE MOD_TYPES,ONLY:NT,NKM_T,NCPLWFMAX,IMT
      USE MOD_RMESH,ONLY:NRMAX,JRWS,FULLPOT,JRCRI
      IMPLICIT NONE
C*--WAVFUN_WRITE_SRA260
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='WAVFUN_WRITE_SRA')
C
C Dummy arguments
C
      INTEGER IFIL,IRTOP,IT,KIRR
      INTEGER IKMCPLWF(NCPLWFMAX,NKMMAX),NCPLWF(NKMMAX)
      COMPLEX*16 JG(NRMAX,NCPLWFMAX,NKMMAX),ZG(NRMAX,NCPLWFMAX,NKMMAX)
C
C Local variables
C
      INTEGER I,IKM,IL,IM,IR,IS,ISOL,K,L
C
C*** End of declarations rewritten by SPAG
C
      IM = IMT(IT)
C
C========================================================================
C                    check consistency of parameters
C========================================================================
      IF ( FULLPOT ) THEN
         IF ( IRTOP.NE.JRCRI(IM) )
     &         CALL STOP_MESSAGE(ROUTINE,'for FP: IRTOP .NE. JRCRI(IM)')
      ELSE
         IF ( IRTOP.NE.JRWS(IM) )
     &         CALL STOP_MESSAGE(ROUTINE,'for ASA: IRTOP .NE. JRWS(IM)')
      END IF
C========================================================================
C
      IF ( FULLPOT ) THEN
C========================================================================
C                               FULLPOT
C========================================================================
C
C --------------------------------------- write wavefunctions for type IT
C
         DO ISOL = 1,NKM_T(IT)
C
            WRITE (IFIL,REC=ISOL+(IT-1)*NKM,ERR=100) IT,'REG',ISOL,
     &             NCPLWF(ISOL),(IKMCPLWF(K,ISOL),K=1,NCPLWF(ISOL)),
     &             ((ZG(IR,IKM,ISOL),IR=1,IRTOP),IKM=1,NCPLWF(ISOL))
C
            IF ( KIRR.EQ.1 ) WRITE (IFIL,REC=ISOL+(IT-1+NT)*NKM,ERR=100)
     &                              IT,'IRR',ISOL,
     &                              ((JG(IR,IKM,ISOL),IR=1,IRTOP),IKM=1,
     &                              NCPLWF(ISOL))
C
         END DO
C
C========================================================================
C                                 ASA
C========================================================================
      ELSE
C
         DO IL = 1,NL
            L = IL - 1
C
C-----------------------------------------------------------------------
C                                                  REGULAR wave function
            WRITE (IFIL,REC=IL+(IT-1)*NL) IT,L,NL,NSPIN,'REG',
     &             ((ZG(I,1,NL*(IS-1)+IL),I=1,IRTOP),IS=1,NSPIN)
C
C-----------------------------------------------------------------------
C                                                IRREGULAR wave function
            IF ( KIRR.EQ.1 ) WRITE (IFIL,REC=IL+(IT-1+NT)*NL) IT,L,NL,
     &                              NSPIN,'IRR',
     &                              ((JG(I,1,NL*(IS-1)+IL),I=1,IRTOP),
     &                              IS=1,NSPIN)
C
         END DO
C
      END IF
C========================================================================
C
      RETURN
C
 100  CONTINUE
      CALL STOP_MESSAGE(ROUTINE,'TROUBLE reading WF - ERROR occured')
      END
C*==wavfun_read_sra.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      SUBROUTINE WAVFUN_READ_SRA(IFIL,IT,KIRR,ZG,JG,IRTOP,NCPLWF,
     &                           IKMCPLWF)
C   ********************************************************************
C   *                                                                  *
C   *  reread the full vector of wave functions solutions for type IT  *
C   *                                                                  *
C   *                                                                  *
C   *      ZG(r,K,ISOL)  ISOL = 1,..,NKM_T  lin. indep. solutions      *
C   *                      K  = 1,..,NCPLWF coupled partial waves      *
C   *                                                                  *
C   *  if KIRR=1  read regular AND irregular wave fuinctions           *
C   *                                                                  *
C   *  the wave functions are written by:                              *
C   *     - <FPNRSSITE>  for FULLPOT calculations                      *
C   *                                                                  *
C   *  scalar relativistic version                                     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,NL,NSPIN
      USE MOD_TYPES,ONLY:NT,NKM_T,NCPLWFMAX,IMT
      USE MOD_RMESH,ONLY:NRMAX,JRWS,FULLPOT,JRCRI
      IMPLICIT NONE
C*--WAVFUN_READ_SRA383
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='WAVFUN_READ_SRA')
C
C Dummy arguments
C
      INTEGER IFIL,IRTOP,IT,KIRR
      INTEGER IKMCPLWF(NCPLWFMAX,NKMMAX),NCPLWF(NKMMAX)
      COMPLEX*16 JG(NRMAX,NCPLWFMAX,NKMMAX),ZG(NRMAX,NCPLWFMAX,NKMMAX)
C
C Local variables
C
      INTEGER I,IKM,IL,IM,IR,IS,ISOL,ISOLP,ITP,K,L,LP,NLP,NSPINP
      CHARACTER*3 STRP
C
C*** End of declarations rewritten by SPAG
C
      IM = IMT(IT)
C
C========================================================================
C                    check consistency of parameters
C========================================================================
      IF ( FULLPOT ) THEN
C
         NCPLWF(:) = 0
C
         IF ( IRTOP.NE.JRCRI(IM) )
     &         CALL STOP_MESSAGE(ROUTINE,'for FP: IRTOP .NE. JRCRI(IM)')
      ELSE
         IF ( NCPLWFMAX.NE.2 )
     &         CALL STOP_MESSAGE(ROUTINE,'for ASA: NCPLWFMAX .NE. 2')
         IF ( IRTOP.NE.JRWS(IM) )
     &         CALL STOP_MESSAGE(ROUTINE,'for ASA: IRTOP .NE. JRWS(IM)')
      END IF
C========================================================================
C
      IF ( FULLPOT ) THEN
C========================================================================
C                               FULLPOT
C========================================================================
C
C ------------------------------------ read in wavefunctions for type IT
C
         DO ISOL = 1,NKM_T(IT)
C
            READ (IFIL,REC=ISOL+(IT-1)*NKM,ERR=100) ITP,STRP,ISOLP,
     &            NCPLWF(ISOL),(IKMCPLWF(K,ISOL),K=1,NCPLWF(ISOL)),
     &            ((ZG(IR,IKM,ISOL),IR=1,IRTOP),IKM=1,NCPLWF(ISOL))
C
            IF ( STRP.NE.'REG' .OR. IT.NE.ITP .OR. ISOL.NE.ISOLP ) THEN
               WRITE (6,*) 'IT STRP ISOL',IT,ITP,'REG',STRP,ISOL,ISOLP
               CALL STOP_MESSAGE(ROUTINE,'TROUBLE reading REG WF')
            END IF
C
            IF ( KIRR.EQ.1 ) THEN
               READ (IFIL,REC=ISOL+(IT-1+NT)*NKM,ERR=100) ITP,STRP,
     &               ISOLP,
     &               ((JG(IR,IKM,ISOL),IR=1,IRTOP),IKM=1,NCPLWF(ISOL))
C
               IF ( STRP.NE.'IRR' .OR. IT.NE.ITP .OR. ISOL.NE.ISOLP )
     &              THEN
                  WRITE (6,*) 'IT STRP ISOL',IT,ITP,'IRR',STRP,ISOL,
     &                        ISOLP
                  CALL STOP_MESSAGE(ROUTINE,'TROUBLE reading IRR WF')
               END IF
            END IF
C
         END DO
C
C========================================================================
C                                 ASA
C========================================================================
      ELSE
C
         DO IL = 1,NL
            L = IL - 1
C
C-----------------------------------------------------------------------
C                                                  REGULAR wave function
            READ (IFIL,REC=IL+(IT-1)*NL) ITP,LP,NLP,NSPINP,STRP,
     &            ((ZG(I,1,NL*(IS-1)+IL),I=1,IRTOP),IS=1,NSPIN)
C
            IF ( STRP.NE.'NRR' .OR. IT.NE.ITP .OR. LP.NE.L .OR. 
     &           NLP.NE.NL .OR. NSPINP.NE.NSPIN ) THEN
               WRITE (6,*) 'IT STRP L',IT,ITP,'REG',STRP,L,LP
               CALL STOP_MESSAGE(ROUTINE,'TROUBLE reading REG WF')
            END IF
C
C-----------------------------------------------------------------------
C                                                IRREGULAR wave function
            IF ( KIRR.EQ.1 ) THEN
               READ (IFIL,REC=IL+(IT-1+NT)*NL) ITP,LP,NLP,NSPINP,STRP,
     &               ((JG(I,1,NL*(IS-1)+IL),I=1,IRTOP),IS=1,NSPIN)
C
               IF ( STRP.NE.'NRI' .OR. IT.NE.ITP .OR. LP.NE.L .OR. 
     &              NLP.NE.NL .OR. NSPINP.NE.NSPIN ) THEN
                  WRITE (6,*) 'IT STRP L',IT,ITP,'REG',STRP,L,LP
                  CALL STOP_MESSAGE(ROUTINE,'TROUBLE reading IRR WF')
               END IF
            END IF
C
         END DO
C
      END IF
C========================================================================
C
      RETURN
C
 100  CONTINUE
      CALL STOP_MESSAGE(ROUTINE,'TROUBLE reading WF - ERROR occured')
      END
