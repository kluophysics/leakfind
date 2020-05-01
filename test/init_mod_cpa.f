C*==init_mod_cpa.f    processed by SPAG 6.70Rc at 13:38 on 14 Dec 2016
      SUBROUTINE INIT_MOD_CPA(NCPA)
C   ********************************************************************
C   *                                                                  *
C   *  initialize all variables used for   CPA   calculations          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:PROGNAME
      USE MOD_CPA,ONLY:CPAFIL,CPAMIX,CPATOL,CPALVL,CPASTART,ICPAALG,
     &    ITCPAMAX,NLCPAWRCFG,USENLCPA,WRCPA,ALPHASRO
      USE MOD_SITES,ONLY:NQ,IQAT
      USE MOD_TYPES,ONLY:NT
      USE MOD_ANGMOM,ONLY:NL,NKMQ,NKQ,NLINQ
      USE MOD_FILES,ONLY:DATSET,LDATSET,FOUND_SECTION,FOUND_INTEGER
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NCPA
C
C Local variables
C
      CHARACTER*10 CTXTALG(4),CTXTSTRT(2),STR10
      INTEGER IT
C
C*** End of declarations rewritten by SPAG
C
      DATA CTXTALG/' MILLS - m',' NESBET   ',' BROYDEN  ',' mixing   '/
      DATA CTXTSTRT/' ATA      ',' AKAI     '/
C
      ICPAALG = 1
      CPASTART = 1
      ITCPAMAX = 20
      CPATOL = 0.00001D0
C
      IF ( NCPA.NE.0 ) THEN
C
         CALL INPUT_FIND_SECTION('CPA',0)
C
         IF ( FOUND_SECTION ) THEN
            CALL SECTION_SET_INTEGER('NITER',ITCPAMAX,9999,0)
            CALL SECTION_SET_REAL('TOL',CPATOL,9999D0,0)
            CALL SECTION_SET_STRING('CPAALG',STR10,'9999',0)
            IF ( STR10(1:6).EQ.'NESBET' ) ICPAALG = 2
            IF ( STR10(1:6).EQ.'MIXING' ) ICPAALG = 4
            CALL SECTION_FIND_KEYWORD('NLCPA',USENLCPA)
            IF ( USENLCPA ) THEN
               CALL SECTION_SET_INTEGER('CPALVL',CPALVL,9999,0)
               CPALVL = MAX(1,CPALVL)
               IF ( .NOT.FOUND_INTEGER ) CPALVL = 2
               CALL SECTION_SET_REAL('MIX',CPAMIX,9999D0,0)
               CALL SECTION_SET_REAL('SRO',ALPHASRO,9999D0,0)
               CALL SECTION_FIND_KEYWORD('WRCFG',NLCPAWRCFG)
            END IF
         END IF
C
C
         IF ( ICPAALG.LT.1 .OR. ICPAALG.GT.4 ) ICPAALG = 1
         IF ( CPASTART.LT.1 .OR. CPASTART.GT.2 ) CPASTART = 1
C
         WRITE (6,99001) CTXTALG(ICPAALG),CTXTSTRT(CPASTART)
         WRITE (6,99002) ITCPAMAX,CPATOL
         IF ( USENLCPA ) THEN
            WRITE (6,99003) CPALVL,CPAMIX,ALPHASRO
            IF ( ICPAALG.NE.1 .AND. ICPAALG.NE.4 ) THEN
               WRITE (6,*) 'CPA-ALG not allowed for NLCPA'
               STOP
            END IF
            IF ( PROGNAME(4:6).EQ.'GEN' ) WRITE (6,99004) NLCPAWRCFG
         END IF
C
         CALL INPUT_FIND_SECTION('CONTROL',0)
C
         IF ( FOUND_SECTION ) THEN
            CALL SECTION_FIND_KEYWORD('WRCPA',WRCPA)
            IF ( WRCPA ) THEN
               CPAFIL = DATSET(1:LDATSET)//'.cpa'
               OPEN (UNIT=12,FILE=CPAFIL)
               WRITE (6,99005) CPAFIL
C
               WRITE (12,'(''SET  '',A)') DATSET
               WRITE (12,99007) 'NQ   ',NQ
               WRITE (12,99007) 'NT   ',NT
               WRITE (12,99007) 'NL   ',NL
               WRITE (12,99007) 'NKT  ',(NKQ(IQAT(1,IT)),IT=1,NT)
               WRITE (12,99007) 'NKMT ',(NKMQ(IQAT(1,IT)),IT=1,NT)
               WRITE (12,99007) 'NLINT',(NLINQ(IQAT(1,IT)),IT=1,NT)
C
            END IF
         END IF
C
         WRITE (6,*) ' '
C
      ELSE
C
         WRITE (6,99006)
C
      END IF
C
99001 FORMAT (/,1X,79('c'),//,10X,'parameters for CPA-iteration:',/,10X,
     &        'algorithm                       ',A10,/,10X,
     &        'start iteration by              ',A10,/)
99002 FORMAT (10X,'max. number of CPA iterations   ',I10,/,10X,
     &        'CPA - tolerance           ',F16.12)
99003 FORMAT (/,10X,'the NLCPA will be used ',/,10X,
     &        'NLCPA-level (cluster size)      ',I10,/,10X,
     &        'mixing parameter                ',F10.5,/,10X,
     &        'SRO - parameter ALPHA           ',F10.5)
99004 FORMAT (10X,'output for configurations       ',L10)
99005 FORMAT (10X,'CPA - cycle written to file     ',A)
99006 FORMAT (//,10X,28('*'),/,10X,'no CPA-iteration required !!',/,10X,
     &        28('*'),/)
99007 FORMAT (A5,15I5,:,(5X,15I5))
      END
