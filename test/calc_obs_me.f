C*==calc_obs_me.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE CALC_OBS_ME(IT,ZGLB,ZFLB,JGLB,JFLB,NCPLWF_LB,
     &                       IKMCPLWF_LB,ZGRA,ZFRA,JGRA,JFRA,NCPLWF_RA,
     &                       IKMCPLWF_RA,MEZZ,MEZJ,MZBZA,MZBJA,NOBS,KEY,
     &                       CHECK_OBS_ME)
C   ********************************************************************
C   *                                                                  *
C   *    read wave function and the calculate matrix elements          *
C   *                                                                  *
C   *             MZBZA  =   < Z^+(E_b) | H_lam | Z(E_a) >             *
C   *             MZBJA  =   < Z^+(E_b) | H_lam | J(E_a) >             *
C   *                      + < J^+(E_b) | H_lam | Z(E_a) >             *
C   *                                                                  *
C   *     with    H_lam =  1, -> sigma, -> l, etc.                     *
C   *                                                                  *
C   *     Z (J) regular (irregular) solutions to Dirac equation        *
C   *                                                                  *
C   *    polarisation lam:        ipol = 1,2,3  ==  (-),(0),(+)        *
C   *                                                                  *
C   *    NOBS: number of observables to deal with                      *
C   *                                                                  *
C   *    KEY = 'S>C':  convert polarisation: SPHERICAL (-),(0),(+)     *
C   *                                     to CARTESIAN (x),(y),(z)     *
C   *                                                                  *
C   *    CHECK_OBS_ME = .TRUE.: the matrix elments will be checked     *
C   *          against the standard ones MEZZ, MEZJ calculated         *
C   *          in  <SSITE>  or  <FPSSITE>                              *
C   *                                                                  *
C   *    NOTE: in the most general case ZGLB, ... are LHS CB solutions *
C   *          ZGLB, ... and  ZGRA, ... may have different couplings   *
C   *                                                                  *
C   *    NOTE: the weight  R2DRDI_W_RADINT = R2DRDI * W_RADINT         *
C   *          will be multiplied to a LOCAL COPY of the radial wave   *
C   *          functions ZGLB, ... for radial integration              *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,L_IKM,AME_G,AME_F,AME_RLM,IMKM_IKM,
     &    NLM_AME_RLM_EXT,K_AME,K_AME_RLM,NMEMAX
      USE MOD_RMESH,ONLY:NRMAX,JRWS,JRMT,JRCRI,FULLPOT,NSF,LMISF,FLMSF,
     &    KLMSF,R2DRDI_W_RADINT
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_TYPES,ONLY:IMT,NCPLWFMAX,NTMAX
      IMPLICIT NONE
C*--CALC_OBS_ME45
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CALC_OBS_ME')
      REAL*8 TOL_CMATCMP,THRESH_CMATCMP
      PARAMETER (TOL_CMATCMP=1D-12,THRESH_CMATCMP=1D-12)
C
C Dummy arguments
C
      LOGICAL CHECK_OBS_ME
      INTEGER IT,NOBS
      CHARACTER*3 KEY
      INTEGER IKMCPLWF_LB(NCPLWFMAX,NKMMAX),
     &        IKMCPLWF_RA(NCPLWFMAX,NKMMAX),NCPLWF_LB(NKMMAX),
     &        NCPLWF_RA(NKMMAX)
      COMPLEX*16 JFLB(NRMAX,NCPLWFMAX,NKMMAX),
     &           JFRA(NRMAX,NCPLWFMAX,NKMMAX),
     &           JGLB(NRMAX,NCPLWFMAX,NKMMAX),
     &           JGRA(NRMAX,NCPLWFMAX,NKMMAX),
     &           MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MZBJA(NKMMAX,NKMMAX,3,NOBS),MZBZA(NKMMAX,NKMMAX,3,NOBS)
     &           ,ZFLB(NRMAX,NCPLWFMAX,NKMMAX),
     &           ZFRA(NRMAX,NCPLWFMAX,NKMMAX),
     &           ZGLB(NRMAX,NCPLWFMAX,NKMMAX),
     &           ZGRA(NRMAX,NCPLWFMAX,NKMMAX)
C
C Local variables
C
      COMPLEX*16 A_FF,A_GG,CADD,S_A_FF,W_JFLB(:),W_JGLB(:),W_ZFLB(:),
     &           W_ZGLB(:),ZFJF,ZFZF,ZGJG,ZGZG
      REAL*8 FLM
      INTEGER IA,IB,IKMA,IKMB,IM,IMKMA,IMKMB,IOBS,IPOL,IR,IRSF,IRTOP,
     &        IRTOP_SPHERE,ISF,LAMA,LAMB,LM,LMSF,NOBS_INIT
      LOGICAL INITIALIZE,SAME_CMATCMP
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE W_JFLB,W_JGLB,W_ZFLB,W_ZGLB
      DATA INITIALIZE/.TRUE./,NOBS_INIT/0/
C
      ALLOCATE (W_JFLB(NRMAX),W_JGLB(NRMAX),W_ZFLB(NRMAX),W_ZGLB(NRMAX))
C
      IM = IMT(IT)
C
      IF ( FULLPOT ) THEN
         IRTOP_SPHERE = JRMT(IM)
         IRTOP = JRCRI(IM)
      ELSE
         IRTOP_SPHERE = JRWS(IM)
         IRTOP = JRWS(IM)
      END IF
C
      IF ( NOBS.GT.NOBS_INIT ) INITIALIZE = .TRUE.
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
C
      IF ( INITIALIZE ) THEN
C
         IF ( ALLOCATED(K_AME) ) DEALLOCATE (K_AME,K_AME_RLM)
         ALLOCATE (K_AME(NKMMAX,NKMMAX))
         ALLOCATE (K_AME_RLM(NKMMAX,NKMMAX,NLM_AME_RLM_EXT))
         K_AME(:,:) = .FALSE.
C
C-----------------------------------------------------------------------
C   flags whether a matrix element for any OBS is non-0
C               K_AME:         <LAM|OBS|LAM'>
C               K_AME_RLM:     <LAM|OBS*RLM|LAM'>
C   for <LAM|OBS|LAM'> diagonal in l, i.e. l{LAM} = l{LAM'}
C   one can restrict LM via l{LM} <= 2*l{LAM} to LM <= LMFP
C-----------------------------------------------------------------------
C
         DO IKMB = 1,NKM
            IMKMB = IMKM_IKM(IKMB)
            DO IKMA = 1,NKM
               IMKMA = IMKM_IKM(IKMA)
C
C-------- include only terms diagonal in l as done so far in case of ASA
               IF ( .NOT.FULLPOT .AND. (L_IKM(IKMA).NE.L_IKM(IKMB)) )
     &              CYCLE
C
               DO IOBS = 1,NOBS
                  DO IPOL = 1,3
C
                     A_GG = AME_G(IKMB,IKMA,IPOL,IOBS)
                     A_FF = AME_F(IKMB,IKMA,IPOL,IOBS)
C
                     IF ( ABS(A_GG).GT.1D-6 .OR. ABS(A_FF).GT.1D-6 )
     &                    K_AME(IKMB,IKMA) = .TRUE.
C
                     IF ( .NOT.FULLPOT ) CYCLE
C
                     DO LM = 1,NLM_AME_RLM_EXT
                        A_GG = AME_RLM(IKMB,IKMA,IPOL,LM,IOBS)
                        A_FF = AME_RLM(IMKMB,IMKMA,IPOL,LM,IOBS)
                        IF ( ABS(A_GG).GT.1D-6 .OR. ABS(A_FF).GT.1D-6 )
     &                       K_AME_RLM(IKMB,IKMA,LM) = .TRUE.
                     END DO
C
                  END DO
               END DO
C
            END DO
         END DO
C
         NOBS_INIT = NOBS
C
         INITIALIZE = .FALSE.
C
      END IF
C=======================================================================
C
      IF ( NOBS.GT.NOBS_INIT )
     &      CALL STOP_MESSAGE(ROUTINE,'NOBS .GT. NOBS_INIT')
C
      MZBZA(:,:,:,:) = C0
      MZBJA(:,:,:,:) = C0
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C                                                     energy B -- LAMB
      DO LAMB = 1,NKM
C
C-----------------------------------------------------------------------
         LOOP_IB:DO IB = 1,NCPLWF_LB(LAMB)
            IKMB = IKMCPLWF_LB(IB,LAMB)
            IMKMB = IMKM_IKM(IKMB)
C
            DO IR = 1,IRTOP
               W_ZGLB(IR) = R2DRDI_W_RADINT(IR,IM)*ZGLB(IR,IB,LAMB)
               W_ZFLB(IR) = R2DRDI_W_RADINT(IR,IM)*ZFLB(IR,IB,LAMB)
               W_JGLB(IR) = R2DRDI_W_RADINT(IR,IM)*JGLB(IR,IB,LAMB)
               W_JFLB(IR) = R2DRDI_W_RADINT(IR,IM)*JFLB(IR,IB,LAMB)
            END DO
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                     energy A -- LAMA
            DO LAMA = 1,NKM
               LOOP_IA:DO IA = 1,NCPLWF_RA(LAMA)
                  IKMA = IKMCPLWF_RA(IA,LAMA)
                  IMKMA = IMKM_IKM(IKMA)
C
C-----------------------------------------------------------------------
C                    regime of muffin tin or ASA sphere
C-----------------------------------------------------------------------
C
                  IF ( K_AME(IKMB,IKMA) ) THEN
C
                     ZGZG = 0D0
                     ZFZF = 0D0
                     ZGJG = 0D0
                     ZFJF = 0D0
C
                     IF ( LAMB.EQ.LAMA ) THEN
C
                        DO IR = 1,IRTOP_SPHERE
                           ZGZG = ZGZG + W_ZGLB(IR)*ZGRA(IR,IA,LAMA)
                           ZFZF = ZFZF + W_ZFLB(IR)*ZFRA(IR,IA,LAMA)
C
C---------------------------------------------------- take average ZJ+JZ
                           ZGJG = ZGJG + W_ZGLB(IR)*JGRA(IR,IA,LAMA)
                           ZGJG = ZGJG + W_JGLB(IR)*ZGRA(IR,IA,LAMA)
C
                           ZFJF = ZFJF + W_ZFLB(IR)*JFRA(IR,IA,LAMA)
                           ZFJF = ZFJF + W_JFLB(IR)*ZFRA(IR,IA,LAMA)
                        END DO
C
                     ELSE
C
                        DO IR = 1,IRTOP_SPHERE
                           ZGZG = ZGZG + W_ZGLB(IR)*ZGRA(IR,IA,LAMA)
                           ZFZF = ZFZF + W_ZFLB(IR)*ZFRA(IR,IA,LAMA)
                        END DO
C
                     END IF
C
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
                     DO IOBS = 1,NOBS
                        DO IPOL = 1,3
C
                           A_GG = AME_G(IKMB,IKMA,IPOL,IOBS)
                           A_FF = AME_F(IKMB,IKMA,IPOL,IOBS)
                           IF ( IOBS.EQ.1 ) THEN
                              S_A_FF = +A_FF
                           ELSE
                              S_A_FF = -A_FF
                           END IF
C
                           CADD = A_GG*ZGZG + S_A_FF*ZFZF
C
                           MZBZA(LAMB,LAMA,IPOL,IOBS)
     &                        = MZBZA(LAMB,LAMA,IPOL,IOBS) + CADD
C
                           IF ( LAMB.EQ.LAMA ) THEN
C
C---------------------------------------------------- take average ZJ+JZ
                              CADD = (A_GG*ZGJG+S_A_FF*ZFJF)*0.5D0
C
                              MZBJA(LAMB,LAMA,IPOL,IOBS)
     &                           = MZBJA(LAMB,LAMA,IPOL,IOBS) + CADD
                           END IF
C
                        END DO
                     END DO
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
C
                  END IF
C
C-----------------------------------------------------------------------
C                interstitial regime in case of  FULLPOT
C-----------------------------------------------------------------------
C
                  LOOP_ISF:DO ISF = 1,NSF(IM)
C
                     IF ( .NOT.FULLPOT ) CYCLE LOOP_ISF
C
                     LMSF = LMISF(ISF,IM)
C
                     IF ( LMSF.GT.NLM_AME_RLM_EXT ) CYCLE LOOP_ISF
                     IF ( KLMSF(LMSF,IM).NE.1 ) CYCLE LOOP_ISF
                     IF ( .NOT.K_AME_RLM(IKMB,IKMA,LMSF) )
     &                    CYCLE LOOP_ISF
C
                     ZGZG = 0D0
                     ZFZF = 0D0
                     ZGJG = 0D0
                     ZFJF = 0D0
C
                     IF ( LAMB.EQ.LAMA ) THEN
C
                        DO IR = IRTOP_SPHERE + 1,JRCRI(IM)
                           IRSF = IR - IRTOP_SPHERE
                           FLM = FLMSF(IRSF,ISF,IM)
C
                           ZGZG = ZGZG + W_ZGLB(IR)*FLM*ZGRA(IR,IA,LAMA)
                           ZFZF = ZFZF + W_ZFLB(IR)*FLM*ZFRA(IR,IA,LAMA)
C
C---------------------------------------------------- take average ZJ+JZ
                           ZGJG = ZGJG + W_ZGLB(IR)*FLM*JGRA(IR,IA,LAMA)
C                          ZGJG = ZGJG + W_JGLB(IR)  NO
C   &                             *FLM*ZGRA(IR,IA,LAMA)
C
                           ZFJF = ZFJF + W_ZFLB(IR)*FLM*JFRA(IR,IA,LAMA)
C                          ZFJF = ZFJF + W_JFLB(IR)  NO
C    &                              *FLM*ZFRA(IR,IA,LAMA)
                        END DO
C
                     ELSE
C
                        DO IR = IRTOP_SPHERE + 1,JRCRI(IM)
                           IRSF = IR - IRTOP_SPHERE
                           FLM = FLMSF(IRSF,ISF,IM)
C
                           ZGZG = ZGZG + W_ZGLB(IR)*FLM*ZGRA(IR,IA,LAMA)
                           ZFZF = ZFZF + W_ZFLB(IR)*FLM*ZFRA(IR,IA,LAMA)
                        END DO
C
                     END IF
C
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
                     DO IOBS = 1,NOBS
                        DO IPOL = 1,3
C
                           A_GG = AME_RLM(IKMB,IKMA,IPOL,LMSF,IOBS)
                           A_FF = AME_RLM(IMKMB,IMKMA,IPOL,LMSF,IOBS)
                           IF ( IOBS.EQ.1 ) THEN
                              S_A_FF = +A_FF
                           ELSE
                              S_A_FF = -A_FF
                           END IF
C
                           CADD = A_GG*ZGZG + S_A_FF*ZFZF
C
                           MZBZA(LAMB,LAMA,IPOL,IOBS)
     &                        = MZBZA(LAMB,LAMA,IPOL,IOBS) + CADD
C
                           IF ( LAMB.EQ.LAMA ) THEN
C
C---------------------------------------------------- take average ZJ+JZ
C                             CADD = (A_GG*ZGJG+S_A_FF*ZFJF)*0.5D0  NO
                              CADD = (A_GG*ZGJG+S_A_FF*ZFJF)
C
                              MZBJA(LAMB,LAMA,IPOL,IOBS)
     &                           = MZBJA(LAMB,LAMA,IPOL,IOBS) + CADD
                           END IF
C
                        END DO
                     END DO
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
C
                  END DO LOOP_ISF
C-----------------------------------------------------------------------
C
               END DO LOOP_IA
C-----------------------------------------------------------------------
            END DO
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
         END DO LOOP_IB
C-----------------------------------------------------------------------
      END DO
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
C
C-----------------------------------------------------------------------
C  convert polarisation: SPHERICAL (-),(0),(+) to CARTESIAN (x),(y),(z)
C-----------------------------------------------------------------------
C
      IF ( KEY.EQ.'S>C' ) THEN
         DO IOBS = 2,NOBS
            CALL CMAT_CONVERT_POLAR(MZBZA(1,1,1,IOBS),'S>C')
            CALL CMAT_CONVERT_POLAR(MZBJA(1,1,1,IOBS),'S>C')
         END DO
      END IF
C
C=======================================================================
      IF ( .NOT.CHECK_OBS_ME ) RETURN
C=======================================================================
C
C???????????????????????????????????????????????????????????????????????
C                 compare with the standard matrix elements
C???????????????????????????????????????????????????????????????????????
C
      CALL CMATCMP(NKM,NKMMAX,3,'MEZZ    CHR',MEZZ(1,1,IT,1),
     &             'MZBZA   CHR',MZBZA(1,1,2,1),THRESH_CMATCMP,
     &             TOL_CMATCMP,SAME_CMATCMP)
C
      CALL CMATCMP(NKM,NKMMAX,3,'MEZJ    CHR',MEZJ(1,1,IT,1),
     &             'MZBJA   CHR',MZBJA(1,1,2,1),THRESH_CMATCMP,
     &             TOL_CMATCMP,SAME_CMATCMP)
C
      IF ( NOBS.EQ.1 ) RETURN
C
      CALL CMATCMP(NKM,NKMMAX,3,'MEZZ    SPN',MEZZ(1,1,IT,2),
     &             'MZBZA   SPN',MZBZA(1,1,2,2),THRESH_CMATCMP,
     &             TOL_CMATCMP,SAME_CMATCMP)
C
      CALL CMATCMP(NKM,NKMMAX,3,'MEZJ    SPN',MEZJ(1,1,IT,2),
     &             'MZBJA   SPN',MZBJA(1,1,2,2),THRESH_CMATCMP,
     &             TOL_CMATCMP,SAME_CMATCMP)
C
      IF ( NOBS.EQ.2 ) RETURN
C
      CALL CMATCMP(NKM,NKMMAX,3,'MEZZ    ORB',MEZZ(1,1,IT,3),
     &             'MZBZA   ORB',MZBZA(1,1,2,3),THRESH_CMATCMP,
     &             TOL_CMATCMP,SAME_CMATCMP)
C
      CALL CMATCMP(NKM,NKMMAX,3,'MEZJ    ORB',MEZJ(1,1,IT,3),
     &             'MZBJA   ORB',MZBJA(1,1,2,3),THRESH_CMATCMP,
     &             TOL_CMATCMP,SAME_CMATCMP)
C
      END
C*==calc_obs_me_selrule.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE CALC_OBS_ME_SELRULE(KEY,NOBS)
C   ********************************************************************
C   *                                                                  *
C   *   flags whether a matrix element for any OBS is non-0            *
C   *               K_AME:         <LAM|OBS|LAM'>                      *
C   *               K_AME_RLM:     <LAM|OBS*RLM|LAM'>                  *
C   *   for <LAM|OBS|LAM'> diagonal in l, i.e. l{LAM} = l{LAM'}        *
C   *   one can restrict LM via l{LM} <= 2*l{LAM}                      *
C   *   to LM <= NLM_AME_RLM_EXTLMFP                                   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,L_IKM,AME_G,AME_F,AME_RLM,IMKM_IKM,
     &    NLM_AME_RLM_EXT,K_AME,K_AME_RLM,K_AME_Z,K_AME_RLM_Z
      USE MOD_RMESH,ONLY:FULLPOT
      IMPLICIT NONE
C*--CALC_OBS_ME_SELRULE433
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CALC_OBS_ME_SELRULE')
C
C Dummy arguments
C
      CHARACTER*1 KEY
      INTEGER NOBS
C
C Local variables
C
      COMPLEX*16 A_FF,A_GG
      INTEGER IKMA,IKMB,IMKMA,IMKMB,IOBS,IPOL,LM,NOBS_INIT_A,NOBS_INIT_Z
C
C*** End of declarations rewritten by SPAG
C
      DATA NOBS_INIT_A/0/,NOBS_INIT_Z/0/
C
      IF ( KEY.EQ.'A' ) THEN
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C     deal with ALL components of vector observables up to NOBS
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
         IF ( NOBS.LE.NOBS_INIT_A ) RETURN
C
         NOBS_INIT_A = NOBS
C
         IF ( .NOT.ALLOCATED(K_AME) ) ALLOCATE (K_AME(NKMMAX,NKMMAX))
         K_AME(:,:) = .FALSE.
C
         IF ( FULLPOT ) THEN
            IF ( .NOT.ALLOCATED(K_AME_RLM) )
     &           ALLOCATE (K_AME_RLM(NKMMAX,NKMMAX,NLM_AME_RLM_EXT))
            K_AME_RLM(:,:,:) = .FALSE.
         END IF
C
         DO IKMB = 1,NKM
            IMKMB = IMKM_IKM(IKMB)
            DO IKMA = 1,NKM
               IMKMA = IMKM_IKM(IKMA)
C
C-------- include only terms diagonal in l as done so far in case of ASA
               IF ( .NOT.FULLPOT .AND. (L_IKM(IKMA).NE.L_IKM(IKMB)) )
     &              CYCLE
C
               DO IOBS = 1,NOBS
                  DO IPOL = 1,3
C
                     A_GG = AME_G(IKMB,IKMA,IPOL,IOBS)
                     A_FF = AME_F(IKMB,IKMA,IPOL,IOBS)
C
                     IF ( ABS(A_GG).GT.1D-6 .OR. ABS(A_FF).GT.1D-6 )
     &                    K_AME(IKMB,IKMA) = .TRUE.
C
                     IF ( .NOT.FULLPOT ) CYCLE
C
                     DO LM = 1,NLM_AME_RLM_EXT
                        A_GG = AME_RLM(IKMB,IKMA,IPOL,LM,IOBS)
                        A_FF = AME_RLM(IMKMB,IMKMA,IPOL,LM,IOBS)
                        IF ( ABS(A_GG).GT.1D-6 .OR. ABS(A_FF).GT.1D-6 )
     &                       K_AME_RLM(IKMB,IKMA,LM) = .TRUE.
                     END DO
C
                  END DO
               END DO
C
            END DO
         END DO
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
      ELSE IF ( KEY.EQ.'Z' ) THEN
CZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
C     deal only with Z component of vector observables up to NOBS
CZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
C
         IF ( NOBS.LE.NOBS_INIT_Z ) RETURN
C
         NOBS_INIT_Z = NOBS
C
         IF ( .NOT.ALLOCATED(K_AME_Z) )
     &        ALLOCATE (K_AME_Z(NKMMAX,NKMMAX))
         K_AME_Z(:,:) = .FALSE.
C
         IF ( FULLPOT ) THEN
            IF ( .NOT.ALLOCATED(K_AME_RLM_Z) )
     &           ALLOCATE (K_AME_RLM_Z(NKMMAX,NKMMAX,NLM_AME_RLM_EXT))
            K_AME_RLM_Z(:,:,:) = .FALSE.
         END IF
C
         IPOL = 2
C
         DO IKMB = 1,NKM
            IMKMB = IMKM_IKM(IKMB)
            DO IKMA = 1,NKM
               IMKMA = IMKM_IKM(IKMA)
C
C-------- include only terms diagonal in l as done so far in case of ASA
               IF ( .NOT.FULLPOT .AND. (L_IKM(IKMA).NE.L_IKM(IKMB)) )
     &              CYCLE
C
               DO IOBS = 1,NOBS
C
                  A_GG = AME_G(IKMB,IKMA,IPOL,IOBS)
                  A_FF = AME_F(IKMB,IKMA,IPOL,IOBS)
C
                  IF ( ABS(A_GG).GT.1D-6 .OR. ABS(A_FF).GT.1D-6 )
     &                 K_AME_Z(IKMB,IKMA) = .TRUE.
C
                  IF ( .NOT.FULLPOT ) CYCLE
C
                  DO LM = 1,NLM_AME_RLM_EXT
                     A_GG = AME_RLM(IKMB,IKMA,IPOL,LM,IOBS)
                     A_FF = AME_RLM(IMKMB,IMKMA,IPOL,LM,IOBS)
                     IF ( ABS(A_GG).GT.1D-6 .OR. ABS(A_FF).GT.1D-6 )
     &                    K_AME_RLM_Z(IKMB,IKMA,LM) = .TRUE.
                  END DO
C
               END DO
C
            END DO
         END DO
CZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
      ELSE
         CALL STOP_MESSAGE(ROUTINE,'KEY not found ')
      END IF
C
      END
