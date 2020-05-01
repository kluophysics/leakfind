C*==mecheck.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MECHECK(TSST,MSST,IPRINT,SSST,MEZZ,MEZJ,NCSTMAX)
C   ********************************************************************
C   *                                                                  *
C   *  subroutine to check the calculation of    MATRIX ELEMENTS       *
C   *                                                                  *
C   *         MZAZB  =   < Z^+(E_a) | H_lam | Z(E_b) >                 *
C   *                                                                  *
C   *         MZBZA  =   < Z^+(E_b) | H_lam | Z(E_a) >                 *
C   *                                                                  *
C   *         MIRR_2 =   < Z^+(E_b) | H_lam * IZAZB | J(E_a) >         *
C   *                  + < Z^+(E_b) | H_lam * IJAZB | Z(E_a) >         *
C   *                                                                  *
C   *         MIRR_3 =   < J^+(E_b) | H_lam * IZAZB | J(E_a) >         *
C   *                  + < Z^+(E_b) | H_lam * IZAJB | Z(E_a) >         *
C   *                                                                  *
C   *         MIRR_4 =   < J^+(E_b) | H_lam * IZAZB | J(E_a) >         *
C   *                  + < Z^+(E_b) | H_lam * IZAJB | Z(E_a) >         *
C   *                                                                  *
C   *    polarisation lam:        ipol = 1,2,3  ==  (-),(z),(+)        *
C   *                                                                  *
C   * ---------------------------------------------------------------- *
C   *                                                                  *
C   *  Checks:                                                         *
C   *                                                                  *
C   *    1  compare MEs among each other                               *
C   *    2  check interal consistency by overwriting Js with Zs        *
C   *         and compare IRR MEs with products of REG MEs             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:ORBPOL
      USE MOD_FILES,ONLY:LDATSET,LRECREAL8,DATSET,IOTMP,FOUND_SECTION,
     &    FOUND_STRING,FOUND_REAL_ARRAY,N_FOUND
      USE MOD_TYPES,ONLY:NTMAX,NCPLWFMAX,IMT,CTL,BT,LCXRAY,NCXRAY
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX,NKM,NLM,NPOL
      USE MOD_RMESH,ONLY:FULLPOT,NRMAX,JRWS
      IMPLICIT NONE
C*--MECHECK39
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER K_CALC_ME
      PARAMETER (K_CALC_ME=3)
C
C Dummy arguments
C
      INTEGER IPRINT,NCSTMAX
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),SSST(NKMMAX,NKMMAX,NTMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      REAL*8 ADEL,BCOR(NTMAX),BCORS(NTMAX),C,EAUX(2),ECOR(NCSTMAX),
     &       FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX),MAXDEV(3),
     &       SUMDEV(3),SZCOR(NCSTMAX),T
      LOGICAL CALCINT,GETIRRSOL,ME_CC_BRA_RWF,NONMAG
      CHARACTER*1 CHPOL(NPOL)
      CHARACTER*2 CL
      COMPLEX*16 ERYDF,ERYDI,MEXME(:,:,:,:),MIRR_2(:,:,:,:,:),
     &           MIRR_3(:,:,:,:,:),MIRR_4(:,:,:,:,:),MZAZB(:,:,:,:),
     &           MZBZA(:,:,:,:),P
      CHARACTER*80 FILNAM
      INTEGER I,IA_ERR,ICST,IFILF,IFILI,IKMCOR(NCSTMAX,2),IKMF1,IKMF2,
     &        IKMI1,IKMI2,IM,IMEFORM,IMEFORM1,IOL,IPOL,IT,ITEST,ITSEL,
     &        ITXRAY,IWME,IXX,IZERO(NCSTMAX),J,JPOL,KAPCOR(NCSTMAX),LFN,
     &        M,MM05COR(NCSTMAX),N,NC,NCST,NKPCOR(NCSTMAX),NMEFORM,NSUM,
     &        W
      CHARACTER*4 INITIALSTATE
      CHARACTER*3 MEFORM,MEFORMTAB(3)
      CHARACTER*15 T15
C
C*** End of declarations rewritten by SPAG
C
      DATA CHPOL/'-','z','+'/
      DATA MEFORMTAB/'NAB','ADA','GRV'/,CALCINT/.TRUE./
C
      ALLOCATABLE MEXME
      ALLOCATABLE MZAZB,MZBZA,MIRR_2,MIRR_3,MIRR_4
C
      CALL ME_CHEBY(TSST,MSST,SSST,MEZZ,MEZJ)
C
      IF ( FULLPOT ) THEN
         NC = NCPLWFMAX
      ELSE
         NC = 2
      END IF
      IF ( NC.LE.0 ) STOP '<MECHECK>  array sizes ???'
C
      WRITE (6,99001)
C
      GETIRRSOL = .TRUE.
      NMEFORM = 3
      IXX = 1
C
      M = NKMMAX
      N = NMEFORM
C
      ALLOCATE (MZAZB(M,M,3,N),MZBZA(M,M,3,N))
      ALLOCATE (MIRR_2(M,M,3,3,N),MIRR_3(M,M,3,3,N),MIRR_4(M,M,3,3,N))
      ALLOCATE (MEXME(NKMMAX,NKMMAX,3,3),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:mecheck -> MIRR'
C
      NONMAG = .FALSE.
      ITSEL = 1
      NCXRAY(1) = 2
      LCXRAY(1) = 1
      INITIALSTATE = 'BAND'
C      INITIALSTATE = 'CORE'
      ERYDI = 0D0
      ERYDF = 0D0
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_FIND_KEYWORD('NONMAG',NONMAG)
         CALL SECTION_SET_INTEGER('ITSEL',ITSEL,9999,0)
         CALL SECTION_SET_INTEGER('NMEFORM',NMEFORM,9999,0)
         CALL SECTION_SET_INTEGER('PRINT',IPRINT,9999,0)
         CALL SECTION_SET_STRING('INITIALSTATE',INITIALSTATE,'9999',0)
         IF ( INITIALSTATE.EQ.'CORE' ) THEN
            CALL SECTION_SET_STRING('CL',CL,'9999',0)
            IF ( FOUND_STRING ) THEN
               NCXRAY(1) = ICHAR(CL(1:1)) - ICHAR('1') + 1
               IF ( CL(2:2).EQ.'S' ) LCXRAY(1) = 0
               IF ( CL(2:2).EQ.'P' ) LCXRAY(1) = 1
               IF ( CL(2:2).EQ.'D' ) LCXRAY(1) = 2
               IF ( CL(2:2).EQ.'F' ) LCXRAY(1) = 3
            END IF
         ELSE IF ( INITIALSTATE.EQ.'BAND' ) THEN
            CALL SECTION_SET_REAL_ARRAY('ERYDI',EAUX,N_FOUND,2,0,9999D0,
     &                                  0)
            IF ( FOUND_REAL_ARRAY ) ERYDI = DCMPLX(EAUX(1),EAUX(2))
         ELSE
            WRITE (6,*) 'INITIALSTATE = ',INITIALSTATE
            WRITE (6,*) 'should be   CORE   or   BAND '
            STOP ' in <MECHECK>'
         END IF
         CALL SECTION_SET_REAL_ARRAY('ERYDF',EAUX,N_FOUND,2,0,9999D0,0)
         IF ( FOUND_REAL_ARRAY ) ERYDF = DCMPLX(EAUX(1),EAUX(2))
C
      END IF
C
      IF ( FULLPOT ) THEN
         WRITE (6,99012)
         NMEFORM = 2
      END IF
      NMEFORM = 2
C
      IT = ITSEL
      IM = IMT(IT)
      C = CTL(IT,1)
C
      IMEFORM1 = 1
C
C=======================================================================
C
      IFILI = 87
      IFILF = 88
C
      IF ( FULLPOT ) THEN
         IOL = (3+NCPLWFMAX+1+NCPLWFMAX*2*NRMAX*2)*LRECREAL8
      ELSE
         IOL = (5+2*(2+2*NRMAX*2))*LRECREAL8
      END IF
C
      OPEN (UNIT=IFILI,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=IOL)
      OPEN (UNIT=IFILF,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=IOL)
C
C-----------------------------------------------------------------------
C
      WRITE (6,99003) 'NPOL      ',NPOL,(CHPOL(I),I=1,NPOL)
      WRITE (6,99002) 'NLM       ',NLM
      WRITE (6,99002) 'NKM       ',NKM
      WRITE (6,99006) ITSEL
C
C ======================================================================
C                     suppress magnetic field - if requested
C ======================================================================
C
C      NONMAG = .TRUE.
      IF ( NONMAG ) THEN
         WRITE (6,99007)
         DO I = 1,JRWS(IM)
            BT(I,IT) = 0D0
         END DO
      END IF
C ======================================================================
C                          initial state
C ======================================================================
C
C
C ======================================================================
C                        initial CORE state
C ======================================================================
      IF ( INITIALSTATE.EQ.'CORE' ) THEN
C
         ITXRAY = IT
         NCXRAY(IT) = NCXRAY(1)
         LCXRAY(IT) = LCXRAY(1)
         NCST = 4*LCXRAY(IT) + 2
C
         WRITE (6,99004) INITIALSTATE,NCXRAY(IT),LCXRAY(IT)
C
         CALL CORE(IPRINT,GCOR,FCOR,ECOR,SZCOR,KAPCOR,MM05COR,NKPCOR,
     &             IKMCOR,IZERO,ITXRAY,BCOR,BCORS,NCSTMAX)
C
         IKMI1 = IKMCOR(1,1)
         IKMI2 = IKMCOR(NCST,1)
C
         ERYDI = DCMPLX(ECOR(1))
C
         DO ICST = 1,NCST
            WRITE (6,*) 'icst, ECOR',ICST,ECOR(ICST)
            WRITE (6,*) MM05COR(ICST),NKPCOR(ICST),KAPCOR(ICST),
     &                  IKMCOR(ICST,1),IKMCOR(ICST,2)
         END DO
C
         WRITE (6,99005) 'initial',ERYDI,IKMI1,IKMI2
C ======================================================================
C                        initial BAND state
C ======================================================================
      ELSE
C
         IF ( ABS(ERYDI).LT.1D-6 ) ERYDI = (0.5D0,0.01D0)
         IKMI1 = 1
         IKMI2 = NKM
C
         WRITE (6,99004) INITIALSTATE
         WRITE (6,99005) 'initial',ERYDI,IKMI1,IKMI2
C
         IF ( FULLPOT ) THEN
C
            CALL FPSSITE(1,1,IFILI,GETIRRSOL,ERYDI,P,IPRINT,TSST,MSST,
     &                   SSST,MEZZ,MEZJ,ORBPOL)
C
         ELSE
C
            CALL SSITE(1,1,IFILI,CALCINT,GETIRRSOL,ERYDI,P,IPRINT,NKM,
     &                 TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
         END IF
C
      END IF
C ======================================================================
C                         final BAND state
C ======================================================================
      IF ( ABS(ERYDF).LT.1D-6 ) ERYDF = (0.8D0,0.01D0)
      IKMF1 = 1
      IKMF2 = NKM
      WRITE (6,99005) 'final  ',ERYDF,IKMF1,IKMF2
C
      IF ( FULLPOT ) THEN
C
         CALL FPSSITE(1,1,IFILF,GETIRRSOL,ERYDF,P,IPRINT,TSST,MSST,SSST,
     &                MEZZ,MEZJ,ORBPOL)
C
      ELSE
C
         CALL SSITE(1,1,IFILF,CALCINT,GETIRRSOL,ERYDF,P,IPRINT,NKM,TSST,
     &              MSST,SSST,MEZZ,MEZJ,ORBPOL)
      END IF
C
C
C ======================================================================
C     Loop over tests       1  compare MEs
C                           2  overwrite J with Z and check IRR MEs
C ======================================================================
C
      LOOP_TEST:DO ITEST = 1,2
C
C ======================================================================
C         replace    J   by   Z  to check irregular subroutines
C ======================================================================
         IF ( ITEST.EQ.2 ) THEN
C
            CALL WF_J_TO_Z(IFILI,IKMI1,IKMI2,IT)
            CALL WF_J_TO_Z(IFILF,IKMF1,IKMF2,IT)
C
         END IF
C ======================================================================
C                            matrix elements
C ======================================================================
C
C        ME_CC_BRA_RWF = .FALSE.
         ME_CC_BRA_RWF = .TRUE.
C
         IT = ITSEL
         C = CTL(IT,1)
C
         LOOP_MEFORM_A:DO IMEFORM = IMEFORM1,NMEFORM
C
            MEFORM = MEFORMTAB(IMEFORM)
C
            IWME = IOTMP
C
            IF ( ITEST.EQ.1 ) THEN
               WRITE (6,99008) ITEST,MEFORM,
     &                         'calculate the matrix elements'
               FILNAM = DATSET(1:LDATSET)//'_check_'//MEFORM//'.dat'
               LFN = LDATSET + 7 + 3 + 4
            ELSE
               WRITE (6,99008) ITEST,MEFORM,
     &                         'calculate the matrix elements for Z = J'
               FILNAM = DATSET(1:LDATSET)//'_check_'//MEFORM//
     &                  '_ZtoJ.dat'
               LFN = LDATSET + 7 + 3 + 4 + 5
            END IF
C
            IF ( IWME.NE.6 ) OPEN (IWME,FILE=FILNAM(1:LFN))
C
            IWME = 0
C
C ================================================================== ADA
            IF ( MEFORM.EQ.'ADA' ) THEN
C
               CALL ME_ALF_ALF(1,NKM,IFILF,ERYDF,1,NKM,IFILI,ERYDI,
     &                         ME_CC_BRA_RWF,IT,MZAZB(1,1,1,IMEFORM),
     &                         MZBZA(1,1,1,IMEFORM),
     &                         MIRR_2(1,1,1,1,IMEFORM),
     &                         MIRR_3(1,1,1,1,IMEFORM),
     &                         MIRR_4(1,1,1,1,IMEFORM),C,K_CALC_ME)
C
C ================================================================== NAB
C
            ELSE IF ( MEFORM.EQ.'NAB' ) THEN
C
               CALL AME_INIT('NAB       ',IWME)
C
               CALL ME_NAB_NAB(IKMF1,IKMF2,IFILF,ERYDF,IKMI1,IKMI2,
     &                         IFILI,ERYDI,ME_CC_BRA_RWF,IT,
     &                         MZAZB(1,1,1,IMEFORM),MZBZA(1,1,1,IMEFORM)
     &                         ,MIRR_2(1,1,1,1,IMEFORM),
     &                         MIRR_3(1,1,1,1,IMEFORM),
     &                         MIRR_4(1,1,1,1,IMEFORM),C,K_CALC_ME)
C
C ================================================================== GRV
C
            ELSE IF ( MEFORM.EQ.'GRV' ) THEN
C
               CALL AME_INIT('GRV       ',IWME)
C
               CALL ME_GRV_GRV(IFILF,ERYDF,IFILI,ERYDI,ME_CC_BRA_RWF,IT,
     &                         MZAZB(1,1,1,IMEFORM),MZBZA(1,1,1,IMEFORM)
     &                         ,MIRR_2(1,1,1,1,IMEFORM),
     &                         MIRR_3(1,1,1,1,IMEFORM),
     &                         MIRR_4(1,1,1,1,IMEFORM),C)
C
            END IF
C
            W = IOTMP
            N = NKM
            M = NKMMAX
            T = 1D-8
            I = IMEFORM
            T15 = MEFORM//'   MZAZB '
            CALL CMATSTRUCT(T15//'(-)',MZAZB(1,1,1,I),N,M,3,3,0,T,W)
            CALL CMATSTRUCT(T15//'(z)',MZAZB(1,1,2,I),N,M,3,3,0,T,W)
            CALL CMATSTRUCT(T15//'(+)',MZAZB(1,1,3,I),N,M,3,3,0,T,W)
            T15 = MEFORM//'   MZBZA '
            CALL CMATSTRUCT(T15//'(-)',MZBZA(1,1,1,I),N,M,3,3,0,T,W)
            CALL CMATSTRUCT(T15//'(z)',MZBZA(1,1,2,I),N,M,3,3,0,T,W)
            CALL CMATSTRUCT(T15//'(+)',MZBZA(1,1,3,I),N,M,3,3,0,T,W)
C
            DO IPOL = 1,3
               DO JPOL = 1,3
                  WRITE (IWME,*) 'MIRR2 IPOL JPOL',IPOL,JPOL
                  CALL CMATSTRUCT('MIRR2  '//MEFORM,
     &                            MIRR_2(1,1,IPOL,JPOL,I),N,M,3,3,0,T,W)
               END DO
            END DO
C
            DO IPOL = 1,3
               DO JPOL = 1,3
                  WRITE (IWME,*) 'MIRR3 IPOL JPOL',IPOL,JPOL
                  CALL CMATSTRUCT('MIRR3  '//MEFORM,
     &                            MIRR_3(1,1,IPOL,JPOL,I),N,M,3,3,0,T,W)
               END DO
            END DO
C
            IF ( IWME.NE.6 ) THEN
               WRITE (6,99011) MEFORM,FILNAM(1:LFN)
               CLOSE (IWME)
            END IF
C
         END DO LOOP_MEFORM_A
C
C ======================================================================
C                  compare results for different ME-forms
C                      using IMEFORM1 as reference
C ======================================================================
C
         IF ( IMEFORM1.EQ.1 .AND. ITEST.EQ.1 ) THEN
C
            WRITE (6,99008) ITEST,'ALL',
     &                  'comparing the results for the various ME-forms'
C
            DO IMEFORM = 1,NMEFORM
               SUMDEV(IMEFORM) = 0D0
               MAXDEV(IMEFORM) = 0D0
            END DO
            NSUM = 0
C
            DO IPOL = 1,NPOL
C
               DO J = IKMI1,IKMI2
                  DO I = IKMF1,IKMF2
C
                     IF ( ABS(MZAZB(I,J,IPOL,1)).GT.1D-6 ) THEN
                        WRITE (6,99010) CHPOL(IPOL),I,J,MEFORMTAB(1),
     &                                  MZAZB(I,J,IPOL,1)
C
                        NSUM = NSUM + 1
                        DO IMEFORM = 2,NMEFORM
                           ADEL = ABS(1D0-MZAZB(I,J,IPOL,IMEFORM)
     &                            /MZAZB(I,J,IPOL,1))
                           SUMDEV(IMEFORM) = SUMDEV(IMEFORM) + ADEL
                           MAXDEV(IMEFORM) = MAX(MAXDEV(IMEFORM),ADEL)
                           WRITE (6,99010) CHPOL(IPOL),I,J,
     &                            MEFORMTAB(IMEFORM),
     &                            MZAZB(I,J,IPOL,IMEFORM),ADEL
                        END DO
                        WRITE (6,*) ' '
C
                     END IF
C
                  END DO
               END DO
C
            END DO
C
            DO IMEFORM = 2,NMEFORM
               WRITE (6,*) 'IMEFORM ',IMEFORM
               SUMDEV(IMEFORM) = SUMDEV(IMEFORM)/DBLE(NSUM)
               WRITE (6,99009) MEFORMTAB(IMEFORM),NSUM,NSUM,
     &                         SUMDEV(IMEFORM),MAXDEV(IMEFORM)
            END DO
         END IF
C
C ======================================================================
C                    check irregular term
C ======================================================================
C
         LOOP_MEFORM_B:DO IMEFORM = 1,NMEFORM
C
            IF ( ITEST.NE.2 ) CYCLE LOOP_MEFORM_B
C
            WRITE (6,99008) ITEST,MEFORMTAB(IMEFORM),
     &                      'checking the irregular terms for  Z = J'
C
C ======================================================================
C  TEST 4    MIRR_2_mn = MZBZA_m*MZAZB_n (overwriting J with Z)
C ======================================================================
C
            DO IPOL = 1,3
               DO JPOL = 1,3
C
                  MEXME(1:NKM,1:NKM,IPOL,JPOL)
     &               = MATMUL(MZBZA(1:NKM,1:NKM,IPOL,IMEFORM),
     &               MZAZB(1:NKM,1:NKM,JPOL,IMEFORM))
C
               END DO
            END DO
C
            WRITE (6,*) 'IMEFORM',IMEFORM,MEFORMTAB(IMEFORM)
            CALL CMPMAT9(MEXME,MIRR_2(1,1,1,1,IMEFORM),'TEST4',IXX,IT)
C
C ======================================================================
C  TEST 5    MIRR_3_nm = (MZAZB_n*MZBZA_m) (overwriting J with Z)
C ======================================================================
C
            DO IPOL = 1,3
               DO JPOL = 1,3
C
                  MEXME(1:NKM,1:NKM,IPOL,JPOL)
     &               = MATMUL(MZAZB(1:NKM,1:NKM,JPOL,IMEFORM),
     &               MZBZA(1:NKM,1:NKM,IPOL,IMEFORM))
C
               END DO
            END DO
C
            WRITE (6,*) 'IMEFORM',IMEFORM,MEFORMTAB(IMEFORM)
            CALL CMPMAT9(MEXME,MIRR_3(1,1,1,1,IMEFORM),'TEST5',IXX,IT)
C
         END DO LOOP_MEFORM_B
C
C ======================================================================
C
      END DO LOOP_TEST
C
      STOP
C
99001 FORMAT (' ',//,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '* *     *  ******   ****   *     *  ******   ****   *     *  *'
     &  ,/,10X,
     &  '* **   **  *       *    *  *     *  *       *    *  *    *   *'
     &  ,/,10X,
     &  '* * * * *  *       *       *     *  *       *       *  *     *'
     &  ,/,10X,
     &  '* *  *  *  *****   *       *******  *****   *       * *      *'
     &  ,/,10X,
     &  '* *     *  *       *       *     *  *       *       *  *     *'
     &  ,/,10X,
     &  '* *     *  *       *    *  *     *  *       *    *  *    *   *'
     &  ,/,10X,
     &  '* *     *  ******   ****   *     *  ******   ****   *     *  *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99002 FORMAT (10X,A10,I10)
99003 FORMAT (10X,A10,I10,5(:,'  (',A,')'))
99004 FORMAT (/,10X,'initial state type: ',A,/,:,/,10X,
     &        'quantum numbers:     n = ',I1,'  l = ',I1,/)
99005 FORMAT (/,10X,A,' state    energy  E = ',2F12.6,//,33X,'IKM =',I9,
     &        '   ...   ',I3,/)
99006 FORMAT (//,10X,'selected atom type   IT =',I3,/)
99007 FORMAT (/,10X,'spin-dependent part B of potential suppressed',/)
99008 FORMAT (//,2(1X,79('#'),/),20X,'<MECHECK>',5X,'ITEST =',I2,5X,
     &        'MEFORM = ',A,/,2(1X,79('#'),/),/,10X,A,//)
99009 FORMAT (/,10X,A,'  number of non-0 terms:       ',I12,/,:,13X,
     &        '  number of deviations:        ',I12,/,13X,
     &        '  average relative deviation:  ',F12.8,/,:,13X,
     &        '  maximum relative deviation:  ',F12.8,/)
99010 FORMAT (10X,'(',A,')',2I3,2X,A,3F15.10)
99011 FORMAT (/,10x,'matrix elements for   MEFORM = ',A,/,10x,
     &        'written to file   ',A,/)
99012 FORMAT (//,10x,32('<'),' Warning ',32('>'),//,10x,'MEFORM = GRV',
     &        ' not working in FULLPOT',/,10x,'will not be tested',//,
     &        10x,36('<'),36('>'),//)
      END
