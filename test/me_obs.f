C*==me_obs_rel.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      SUBROUTINE ME_OBS_REL(IFILB,IFILA,TYPE_B_EQ_A,DZZ,DZJ,DZZBA,NOBS)
C   ********************************************************************
C   *                                                                  *
C   *   general routine to calculate matriy elements of the type       *
C   *                                                                  *
C   *    DZZBA = <Z(B)|Z(A)>   and   <Z(B)|sigma_i|Z(A)>  i=x,y,z      *
C   *                                                                  *
C   *   for TYPE_B_EQ_A:                                               *
C   *                                                                  *
C   *    .T.:   A and B are the same atom types                        *
C   *           DZZ and DZJ are returned                               *
C   *           IFILA and IFILB may differ                             *
C   *                                                                  *
C   *    .F.:   A and B may be different atom types (for NT>1)         *
C   *           DZZBA and DZJ are returned                             *
C   *           IFILA and IFILB may differ                             *
C   *                                                                  *
C   *   NOTE: the atoms type A and B have to be on the same site IQ    *
C   *                                                                  *
C   *   CHECK_OBS_ME always FALSE ! no CHECK for  MEZZ and MEZJ        *
C   *                                                                  *
C   ********************************************************************
      USE MOD_RMESH,ONLY:NRMAX,JRWS,JRCRI,FULLPOT
      USE MOD_ANGMOM,ONLY:NKMMAX,NMVECMAX,NCPLWF,IDOS,ISMT
      USE MOD_TYPES,ONLY:NTMAX,IMT,ITBOT,ITTOP,NCPLWFMAX,IKMCPLWF
      USE MOD_SITES,ONLY:IQAT
      USE MOD_CONSTANTS,ONLY:CI,SQRT_2
      USE MOD_CALCMODE,ONLY:LHS_SOL_EQ_RHS_SOL
      IMPLICIT NONE
C*--ME_OBS_REL31
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ME_OBS_REL')
      LOGICAL CHECK_OBS_ME
      PARAMETER (CHECK_OBS_ME=.FALSE.)
C
C Dummy arguments
C
      INTEGER IFILA,IFILB,NOBS
      LOGICAL TYPE_B_EQ_A
      COMPLEX*16 DZJ(NKMMAX,NKMMAX,NTMAX,NOBS),
     &           DZZ(NKMMAX,NKMMAX,NTMAX,NOBS),
     &           DZZBA(NKMMAX,NKMMAX,NTMAX,NTMAX,NOBS)
C
C Local variables
C
      INTEGER IA_ERR,IFIL_LHSB,IFIL_RHSA,IM,IQA,IQB,IRTOP,ITA,ITB,NMVEC
      COMPLEX*16 JFLB(:,:,:),JFRA(:,:,:),JGLB(:,:,:),JGRA(:,:,:),
     &           MEZJ(:,:,:,:),MEZZ(:,:,:,:),MZBJA_PO(:,:,:,:),
     &           MZBZA_PO(:,:,:,:),ZFLB(:,:,:),ZFRA(:,:,:),ZGLB(:,:,:),
     &           ZGRA(:,:,:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE ZGRA,ZFRA,ZGLB,ZFLB,JGRA,JFRA,JGLB,JFLB
      ALLOCATABLE MZBJA_PO,MZBZA_PO,MEZZ,MEZJ
C
      ALLOCATE (MEZZ(1,1,1,1),MEZJ(1,1,1,1))
      ALLOCATE (ZGRA(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZFRA(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JGRA(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JFRA(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZGRA')
C
      IF ( .NOT.TYPE_B_EQ_A .OR. .NOT.LHS_SOL_EQ_RHS_SOL .OR. 
     &     (IFILA.NE.IFILB) ) THEN
         ALLOCATE (ZGLB(NRMAX,NCPLWFMAX,NKMMAX))
         ALLOCATE (ZFLB(NRMAX,NCPLWFMAX,NKMMAX))
         ALLOCATE (JGLB(NRMAX,NCPLWFMAX,NKMMAX))
         ALLOCATE (JFLB(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZGLB')
      END IF
C
      ALLOCATE (MZBJA_PO(NKMMAX,NKMMAX,3,NMVECMAX))
      ALLOCATE (MZBZA_PO(NKMMAX,NKMMAX,3,NMVECMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MZBZA_PO')
C
C=======================================================================
C
C         set the chanel number for the LHS wave functions
C
      IFIL_RHSA = IFILA
C
      CALL SET_IFIL_LHS(IFILB,IFIL_LHSB)
C
C=======================================================================
C
      IF ( NOBS.EQ.1 ) THEN
         NMVEC = 1
      ELSE IF ( NOBS.EQ.4 ) THEN
         NMVEC = 2
      ELSE IF ( NOBS.EQ.7 ) THEN
         NMVEC = 3
      ELSE
         CALL STOP_MESSAGE(ROUTINE,'NOBS > 7')
      END IF
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      LOOP_ITA:DO ITA = ITBOT,ITTOP
C
         IQA = IQAT(1,ITA)
C
         IM = IMT(ITA)
         IF ( FULLPOT ) THEN
            IRTOP = JRCRI(IM)
         ELSE
            IRTOP = JRWS(IM)
         END IF
C
C ----------------------------------------- read in wavefunctions for RA
C
         CALL WAVFUN_READ_REL(IFIL_RHSA,ITA,1,ZGRA,ZFRA,JGRA,JFRA,IRTOP,
     &                        NCPLWF,IKMCPLWF)
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
         LOOP_ITB:DO ITB = ITBOT,ITTOP
C
            IQB = IQAT(1,ITB)
C
            IF ( IQA.NE.IQB .OR. (TYPE_B_EQ_A .AND. ITA.NE.ITB) )
     &           CYCLE LOOP_ITB
C
C-----------------------------------------------------------------------
C          calculate matrix elements  MZBZA_PO and MZBJA_PO
C-----------------------------------------------------------------------
C
            IF ( ITA.EQ.ITB .AND. LHS_SOL_EQ_RHS_SOL .AND. 
     &           (IFILA.EQ.IFILB) ) THEN
C
               CALL CALC_OBS_ME(ITB,ZGRA,ZFRA,JGRA,JFRA,NCPLWF,IKMCPLWF,
     &                          ZGRA,ZFRA,JGRA,JFRA,NCPLWF,IKMCPLWF,
     &                          MEZZ,MEZJ,MZBZA_PO,MZBJA_PO,NMVEC,'XXX',
     &                          CHECK_OBS_ME)
C
            ELSE
C ----------------------------------------- read in wavefunctions for LB
C
               CALL WAVFUN_READ_REL(IFIL_LHSB,ITB,1,ZGLB,ZFLB,JGLB,JFLB,
     &                              IRTOP,NCPLWF,IKMCPLWF)
C
               CALL CALC_OBS_ME(ITB,ZGLB,ZFLB,JGLB,JFLB,NCPLWF,IKMCPLWF,
     &                          ZGRA,ZFRA,JGRA,JFRA,NCPLWF,IKMCPLWF,
     &                          MEZZ,MEZJ,MZBZA_PO,MZBJA_PO,NMVEC,'XXX',
     &                          CHECK_OBS_ME)
C
            END IF
C
C=======================================================================
C                       store overlap matrix element
C    convert to cartesian polarisation and store spin matrix element
C=======================================================================
C
C-----------------------------------------------------------------------
C                    atom type A and B are the same
C-----------------------------------------------------------------------
            IF ( TYPE_B_EQ_A ) THEN
C
               DZZ(:,:,ITA,IDOS) = MZBZA_PO(:,:,1,IDOS)
C
               IF ( NOBS.GE.4 ) THEN
                  DZZ(:,:,ITA,2) = (MZBZA_PO(:,:,1,ISMT)-MZBZA_PO(:,:,3,
     &                             ISMT))/SQRT_2
                  DZZ(:,:,ITA,3) = (MZBZA_PO(:,:,1,ISMT)+MZBZA_PO(:,:,3,
     &                             ISMT))*CI/SQRT_2
                  DZZ(:,:,ITA,4) = MZBZA_PO(:,:,2,ISMT)
               END IF
C
            ELSE
C
C-----------------------------------------------------------------------
C                  atom type A and B may be different
C-----------------------------------------------------------------------
C
               DZZBA(:,:,ITB,ITA,IDOS) = MZBZA_PO(:,:,1,IDOS)
C
               IF ( NOBS.GE.4 ) THEN
                  DZZBA(:,:,ITB,ITA,2)
     &               = (MZBZA_PO(:,:,1,ISMT)-MZBZA_PO(:,:,3,ISMT))
     &               /SQRT_2
                  DZZBA(:,:,ITB,ITA,3)
     &               = (MZBZA_PO(:,:,1,ISMT)+MZBZA_PO(:,:,3,ISMT))
     &               *CI/SQRT_2
                  DZZBA(:,:,ITB,ITA,4) = MZBZA_PO(:,:,2,ISMT)
               END IF
C
            END IF
C
C-----------------------------------------------------------------------
C                   irregular term: type diagonal
C-----------------------------------------------------------------------
            IF ( ITA.EQ.ITB ) THEN
C
               DZJ(:,:,ITA,IDOS) = MZBJA_PO(:,:,1,IDOS)
C
               IF ( NOBS.GE.4 ) THEN
                  DZJ(:,:,ITA,2) = (MZBJA_PO(:,:,1,ISMT)-MZBJA_PO(:,:,3,
     &                             ISMT))/SQRT_2
                  DZJ(:,:,ITA,3) = (MZBJA_PO(:,:,1,ISMT)+MZBJA_PO(:,:,3,
     &                             ISMT))*CI/SQRT_2
                  DZJ(:,:,ITA,4) = MZBJA_PO(:,:,2,ISMT)
               END IF
C
            END IF
C
         END DO LOOP_ITB
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      END DO LOOP_ITA
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      END
C*==me_obs_sra.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      SUBROUTINE ME_OBS_SRA(IFIL,MEZZ,MEZJ,DZZ,DZJ)
C   ********************************************************************
C   *                                                                  *
C   *   calculate the overlap matrix elements                          *
C   *                                                                  *
C   *   DZZ = <Z(B)|Z(A)>    and    DZJ = <Z(B)|J(A)>                  *
C   *                                                                  *
C   *   CHECK: the results have to agree with the                      *
C   *          standard matrix elements MEZZ and MEZJ                  *
C   *                                                                  *
C   *                      -- SRA VERSION --                           *
C   *                                                                  *
C   ********************************************************************
      USE MOD_TYPES,ONLY:ITBOT,ITTOP,NTMAX,NT,IMT,NCPLWFMAX,IKMCPLWF
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,NMEMAX,NCPLWF
      USE MOD_RMESH,ONLY:JRCRI,NRMAX
      USE MOD_FILES,ONLY:IPRINT
      IMPLICIT NONE
C*--ME_OBS_SRA250
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ME_OBS_SRA')
      REAL*8 TOL
      PARAMETER (TOL=1D-12)
      LOGICAL CHECK
      PARAMETER (CHECK=.FALSE.)
C
C Dummy arguments
C
      INTEGER IFIL
      COMPLEX*16 DZJ(NKMMAX,NKMMAX,NTMAX),DZZ(NKMMAX,NKMMAX,NTMAX),
     &           MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX)
C
C Local variables
C
      INTEGER I,I1,I2,IA_ERR,IKM,IKMP,IM,IRTOP,IT,ITP,J,K,NKAME,NKMAME
      COMPLEX*16 JGB(:,:,:),ZGA(:,:,:),ZGB(:,:,:)
      COMPLEX*16 DZZTMP(:,:,:,:)
      CHARACTER*3 STRP
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE JGB,ZGA,ZGB,DZZTMP
C
C------------------------------------------- NKAME should be > NK=2*NL-1
      NKAME = 2*NINT(SQRT(DBLE(NKM/2))) + 1
      NKMAME = 2*((NKAME+1)/2)**2
      IF ( NKM.GT.NKMAME ) CALL STOP_MESSAGE(ROUTINE,'NKM > NKAME')
C
      ALLOCATE (JGB(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZGA(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (DZZTMP(NKMMAX,NKMMAX,NTMAX,NTMAX))
      ALLOCATE (ZGB(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZFB')
C
C-----------------------------------------------------------------------
C
      IF ( CHECK ) WRITE (6,99002) ROUTINE(1:LEN_TRIM(ROUTINE))
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = ITBOT,ITTOP
         IM = IMT(IT)
         IRTOP = JRCRI(IM)
C
C ----------------------------------------- read in wavefunctions for IT
C
         DO IKM = 1,NKM
C
            READ (IFIL,REC=IKM+(IT-1)*NKM) ITP,STRP,IKMP,
     &            NCPLWF(IKM),
     &            (IKMCPLWF(J,IKM),J=1,NCPLWF(IKM)),
     &            ((ZGA(I,K,IKM),I=1,IRTOP),K=1,NCPLWFMAX)

C
            IF ( STRP.NE.'REG' .OR. IT.NE.ITP ) THEN
               WRITE (6,*) '##### TROUBLE reading WF:',ITP,STRP,IKMP
               CALL STOP_MESSAGE(ROUTINE,'TROUBLE reading WF')
            END IF
C
         END DO
C
C ----------------------------------------- read in wavefunctions for IT
         DO IKM = 1,NKM
C
            READ (IFIL,REC=IKM+(IT-1)*NKM) ITP,STRP,IKMP,NCPLWF(IKM),
     &            (IKMCPLWF(J,IKM),J=1,NCPLWF(IKM)),
     &            ((ZGB(I,K,IKM),I=1,IRTOP),K=1,NCPLWFMAX)
C
            IF ( STRP.NE.'REG' .OR. IT.NE.ITP ) THEN
               WRITE (6,*) '##### TROUBLE reading WF:',ITP,STRP,IKMP
               CALL STOP_MESSAGE(ROUTINE,'TROUBLE reading WF')
            END IF
C
            READ (IFIL,REC=IKM+(IT-1+NT)*NKM) ITP,STRP,IKMP,
     &            ((JGB(I,K,IKM),I=1,IRTOP),K=1,NCPLWF(IKM))
C
            IF ( STRP.NE.'IRR' .OR. IT.NE.ITP ) THEN
               WRITE (6,*) '##### TROUBLE reading WF:',ITP,STRP,IKMP
               CALL STOP_MESSAGE(ROUTINE,'TROUBLE reading WF')
            END IF
C
         END DO
C
         CALL ME_OBS_AUX_SRA(ZGA,ZGB,JGB,DZZTMP,DZJ,IT,IT)
		 DZZ(:,:,IT)=DZZTMP(:,:,IT,IT)
C
C-----------------------------------------------------------------------
C
C cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C                     check consistency of data sets
C cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
         IF ( CHECK .OR. IPRINT.GE.2 ) THEN
            DO I1 = 1,NKM
               DO I2 = 1,NKM
                  IF ( ABS(MEZZ(I1,I2,IT,1)-DZZ(I1,I2,IT)).GT.TOL )
     &                 WRITE (6,99001) ROUTINE(1:LEN_TRIM(ROUTINE)),IT,
     &                                 I1,I2,'DZZ',MEZZ(I1,I2,IT,1),
     &                                 DZZ(I1,I2,IT)
                  IF ( ABS(MEZJ(I1,I2,IT,1)-DZJ(I1,I2,IT)).GT.TOL )
     &                 WRITE (6,99001) ROUTINE(1:LEN_TRIM(ROUTINE)),IT,
     &                                 I1,I2,'DZJ',MEZJ(I1,I2,IT,1),
     &                                 DZJ(I1,I2,IT)
               END DO
            END DO
         END IF
C cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      DEALLOCATE (JGB,ZGA,ZGB)
C
99001 FORMAT (' WARNING from <',A,'> IT=',I2,' (I1,I2)=',2I3,1X,A,1X,
     &        2E14.7,/,(51X,2E14.7))
99002 FORMAT (/,10X,'<',A,'>  called with CHECK=.TRUE. ',/)
      END
C*==me_obs_aux_sra.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      SUBROUTINE ME_OBS_AUX_SRA(ZGA,ZGB,JGB,DZZ,DZJ,ITA,ITB)
C   ********************************************************************
C   *                                                                  *
C   *  calculate all matrix elements                                   *
C   *               MEZZ =  < Z_LAM | A | Z_LAM' >                     *
C   *               MEZJ =  < Z_LAM | A | J_LAM' >                     *
C   *  with A = 1, sig_z             (index IME)                       *
C   *                                                                  *
C   *  this is a modified version of <FPNRMATELM> that calculates also *
C   *  the matrix elements for the different types in bra and ket      *
C   *  it also includes the spin loop                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKMMAX,NCPLWF,NL,NLM,L_LM,M_LM,NSPIN
      USE MOD_RMESH,ONLY:NRMAX,NSF,LMISF,NPAN,JRCUT,WINTLM,
     &    R2DRDI_W_RADINT,NM
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_TYPES,ONLY:IMT,NTMAX,NCPLWFMAX,IKMCPLWF
      IMPLICIT NONE
C*--ME_OBS_AUX_SRA404
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ME_OBS_AUX_SRA')
C
C Dummy arguments
C
      INTEGER ITA,ITB
      COMPLEX*16 DZJ(NKMMAX,NKMMAX,NTMAX),DZZ(NKMMAX,NKMMAX,NTMAX,NTMAX)
     &           ,JGB(NRMAX,NCPLWFMAX,NKMMAX),
     &           ZGA(NRMAX,NCPLWFMAX,NKMMAX),ZGB(NRMAX,NCPLWFMAX,NKMMAX)
C
C Local variables
C
      REAL*8 G123,W,ZZGNT(:,:,:,:)
      REAL*8 GAUNT_RYLM
      INTEGER I1,I2,IA_ERR,ILM1,ILM2,IM,IS,ISF,J,JLMS1,JLMS2,JOFF,JSF,
     &        L1,L2,LM,LM1,LM2,LMAX,LMSOFF,LSF,M1,M2,MSF,NSFZZMAX
      LOGICAL INITIALIZE
      COMPLEX*16 SIRR,SREG,ZGJG,ZGZG
      SAVE NSFZZMAX,ZZGNT
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
C
      ALLOCATABLE ZZGNT
C
C=======================================================================
C                           INITIALIZE
C=======================================================================
      IF ( INITIALIZE ) THEN
C
C-----------------------------------------------------------------------
C             Gaunt coefficients for REAL spherical harmonics
C-----------------------------------------------------------------------
C
         LMAX = NL - 1
C
         NSFZZMAX = 0
         DO IM = 1,NM
            DO ISF = 1,NSF(IM)
               LM = LMISF(ISF,IM)
               LSF = L_LM(LM)
               IF ( LSF.LE.2*LMAX ) NSFZZMAX = MAX(NSFZZMAX,ISF)
            END DO
         END DO
C
         ALLOCATE (ZZGNT(NLM,NLM,NSFZZMAX,NM),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZZGNT')
         ZZGNT(:,:,:,:) = 0D0
C
         DO IM = 1,NM
            DO ISF = 1,MIN(NSFZZMAX,NSF(IM))
               LM = LMISF(ISF,IM)
               LSF = L_LM(LM)
               MSF = M_LM(LM)
               DO LM2 = 1,NLM
                  L2 = L_LM(LM2)
                  M2 = M_LM(LM2)
                  DO LM1 = 1,NLM
                     L1 = L_LM(LM1)
                     M1 = M_LM(LM1)
                     ZZGNT(LM1,LM2,ISF,IM)
     &                  = GAUNT_RYLM(L1,M1,L2,M2,LSF,MSF)
                  END DO
               END DO
            END DO
         END DO
C
         INITIALIZE = .FALSE.
C
      END IF
C=======================================================================
C
      IM = IMT(ITB)
C
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      DO IS = 1,NSPIN
C
         LMSOFF = NLM*(IS-1)
C
         DO JLMS2 = LMSOFF + 1,LMSOFF + NLM
            DO I2 = 1,NCPLWF(JLMS2)
               ILM2 = IKMCPLWF(I2,JLMS2)
C
               DO JLMS1 = LMSOFF + 1,LMSOFF + NLM
                  DO I1 = 1,NCPLWF(JLMS1)
                     ILM1 = IKMCPLWF(I1,JLMS1)
C
C-----------------------------------------------------------------------
C               integrals within muffin-tin regime
C-----------------------------------------------------------------------
C
                     IF ( ILM1.EQ.ILM2 ) THEN
                        SREG = C0
                        IF ( JLMS1.EQ.JLMS2 ) THEN
                           SIRR = C0
                           DO J = 1,JRCUT(1,IM)
                              W = R2DRDI_W_RADINT(J,IM)
                              ZGZG = ZGA(J,I2,JLMS2)*ZGB(J,I1,JLMS1)
                              ZGJG = ZGA(J,I2,JLMS2)*JGB(J,I1,JLMS1)
                              SREG = SREG + W*ZGZG
                              SIRR = SIRR + W*ZGJG
                           END DO
                           DZJ(JLMS1,JLMS2,ITA) = DZJ(JLMS1,JLMS2,ITA)
     &                        + SIRR
                        ELSE
                           DO J = 1,JRCUT(1,IM)
                              W = R2DRDI_W_RADINT(J,IM)
                              ZGZG = ZGA(J,I2,JLMS2)*ZGB(J,I1,JLMS1)
                              SREG = SREG + W*ZGZG
                           END DO
                        END IF
                        DZZ(JLMS1,JLMS2,ITA,ITB)
     &                     = DZZ(JLMS1,JLMS2,ITA,ITB) + SREG
                     END IF
C
C-----------------------------------------------------------------------
C            integrals within interstitial regime
C-----------------------------------------------------------------------
C
                     DO ISF = 1,MIN(NSFZZMAX,NSF(IM))
                        G123 = ZZGNT(ILM1,ILM2,ISF,IM)
                        IF ( ABS(G123).GT.1D-8 ) THEN
C
                           JOFF = JRCUT(1,IM)
                           SREG = C0
                           IF ( JLMS1.EQ.JLMS2 ) THEN
                              SIRR = C0
                              DO J = JRCUT(1,IM) + 1,JRCUT(NPAN(IM),IM)
                                 JSF = J - JOFF
                                 W = WINTLM(JSF,ISF,IM)
                                 ZGZG = ZGA(J,I2,JLMS2)*ZGB(J,I1,JLMS1)
                                 ZGJG = ZGA(J,I2,JLMS2)*JGB(J,I1,JLMS1)
                                 SREG = SREG + W*ZGZG
                                 SIRR = SIRR + W*ZGJG
                              END DO
                              DZJ(JLMS1,JLMS2,ITA)
     &                           = DZJ(JLMS1,JLMS2,ITA) + SIRR*G123
                           ELSE
                              DO J = JRCUT(1,IM) + 1,JRCUT(NPAN(IM),IM)
                                 JSF = J - JOFF
                                 W = WINTLM(JSF,ISF,IM)
                                 ZGZG = ZGA(J,I2,JLMS2)*ZGB(J,I1,JLMS1)
                                 SREG = SREG + W*ZGZG
                              END DO
                           END IF
                           DZZ(JLMS1,JLMS2,ITA,ITB)
     &                        = DZZ(JLMS1,JLMS2,ITA,ITB) + SREG*G123
                        END IF
                     END DO
C
                  END DO
               END DO
C
            END DO
         END DO
C
      END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
      END
C*==blochint.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      SUBROUTINE BLOCHINT(MEZZ,MEZJ,DZZAB,DZJAB,NBSFOP)
C   ********************************************************************
C   *                                                                  *
C   *   calculate the overlap and spin matrix elements                 *
C   *                                                                  *
C   *        <Z(B)|Z(A)>   and   <Z(B)|sigma_i|Z(A)>  i=x,y,z          *
C   *                                                                  *
C   *   used for the Bloch spectral function  A(k,E)                   *
C   *                                                                  *
C   *   NOTE: the atoms type A and B have to be on the same site IQ    *
C   *                                                                  *
C   *   CHECK_OBS_ME: for A = B the results have to agree with the     *
C   *                 standard matrix elements MEZZ and MEZJ           *
C   *                                                                  *
C   ********************************************************************
      USE MOD_RMESH,ONLY:NRMAX,JRWS,JRCRI,FULLPOT
      USE MOD_FILES,ONLY:IFILCBWF
      USE MOD_ANGMOM,ONLY:NKMMAX,NMVECMAX,NCPLWF,IDOS,ISMT
      USE MOD_TYPES,ONLY:NTMAX,IMT,ITBOT,ITTOP,NCPLWFMAX,IKMCPLWF
      USE MOD_SITES,ONLY:IQAT
      USE MOD_CONSTANTS,ONLY:CI,SQRT_2
      USE MOD_CALCMODE,ONLY:LHS_SOL_EQ_RHS_SOL
      IMPLICIT NONE
C*--BLOCHINT610
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='BLOCHINT')
      LOGICAL CHECK_OBS_ME
      PARAMETER (CHECK_OBS_ME=.FALSE.)
C
C Dummy arguments
C
      INTEGER NBSFOP
      COMPLEX*16 DZJAB(NKMMAX,NKMMAX,NTMAX,NBSFOP),
     &           DZZAB(NKMMAX,NKMMAX,NTMAX,NTMAX,NBSFOP),
     &           MEZJ(NKMMAX,NKMMAX,NTMAX,NMVECMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMVECMAX)
C
C Local variables
C
      INTEGER IA_ERR,IFIL_LHSB,IFIL_RHSA,IM,IQA,IQB,IRTOP,ITA,ITB,NMVEC
      COMPLEX*16 JFLB(:,:,:),JFRA(:,:,:),JGLB(:,:,:),JGRA(:,:,:),
     &           MZBJA_PO(:,:,:,:),MZBZA_PO(:,:,:,:),ZFLB(:,:,:),
     &           ZFRA(:,:,:),ZGLB(:,:,:),ZGRA(:,:,:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE ZGRA,ZFRA,ZGLB,ZFLB,JGRA,JFRA,JGLB,JFLB
      ALLOCATABLE MZBJA_PO,MZBZA_PO
C
      ALLOCATE (ZGRA(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZFRA(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JGRA(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JFRA(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZGRA')
C
      ALLOCATE (ZGLB(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZFLB(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JGLB(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JFLB(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZGLB')
C
      ALLOCATE (MZBJA_PO(NKMMAX,NKMMAX,3,NMVECMAX))
      ALLOCATE (MZBZA_PO(NKMMAX,NKMMAX,3,NMVECMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MZBZA_PO')
C
C=======================================================================
C
C         set the chanel number for the LHS wave functions
C
      IFIL_RHSA = IFILCBWF
C
      CALL SET_IFIL_LHS(IFILCBWF,IFIL_LHSB)
C
C=======================================================================
C
      NMVEC = 2
C
      DZZAB(:,:,:,:,:) = 0D0
      DZJAB(:,:,:,:) = 0D0
C
      MZBZA_PO(:,:,:,:) = 0D0
      MZBJA_PO(:,:,:,:) = 0D0
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      LOOP_ITA:DO ITA = ITBOT,ITTOP
C
         IQA = IQAT(1,ITA)
C
         IM = IMT(ITA)
         IF ( FULLPOT ) THEN
            IRTOP = JRCRI(IM)
         ELSE
            IRTOP = JRWS(IM)
         END IF
C
C ----------------------------------------- read in wavefunctions for RA
C
         CALL WAVFUN_READ_REL(IFIL_RHSA,ITA,1,ZGRA,ZFRA,JGRA,JFRA,IRTOP,
     &                        NCPLWF,IKMCPLWF)
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
         LOOP_ITB:DO ITB = ITBOT,ITTOP
C
            IQB = IQAT(1,ITB)
C
            IF ( IQA.NE.IQB ) CYCLE LOOP_ITB
C
C-----------------------------------------------------------------------
C          calculate matrix elements  MZBZA_PO and MZBJA_PO
C-----------------------------------------------------------------------
C
            IF ( ITA.EQ.ITB .AND. LHS_SOL_EQ_RHS_SOL ) THEN
C
               CALL CALC_OBS_ME(ITB,ZGRA,ZFRA,JGRA,JFRA,NCPLWF,IKMCPLWF,
     &                          ZGRA,ZFRA,JGRA,JFRA,NCPLWF,IKMCPLWF,
     &                          MEZZ,MEZJ,MZBZA_PO,MZBJA_PO,NMVEC,'XXX',
     &                          CHECK_OBS_ME)
C
            ELSE
C ----------------------------------------- read in wavefunctions for LB
C
               CALL WAVFUN_READ_REL(IFIL_LHSB,ITB,1,ZGLB,ZFLB,JGLB,JFLB,
     &                              IRTOP,NCPLWF,IKMCPLWF)
C
               CALL CALC_OBS_ME(ITB,ZGLB,ZFLB,JGLB,JFLB,NCPLWF,IKMCPLWF,
     &                          ZGRA,ZFRA,JGRA,JFRA,NCPLWF,IKMCPLWF,
     &                          MEZZ,MEZJ,MZBZA_PO,MZBJA_PO,NMVEC,'XXX',
     &                          CHECK_OBS_ME)
C
            END IF
C
C-----------------------------------------------------------------------
C                       store overlap matrix element
C-----------------------------------------------------------------------
C
            DZZAB(:,:,ITA,ITB,IDOS) = MZBZA_PO(:,:,1,IDOS)
C
            IF ( ITA.EQ.ITB ) DZJAB(:,:,ITA,IDOS) = MZBJA_PO(:,:,1,IDOS)
C
C-----------------------------------------------------------------------
C    convert to cartesian polarisation and store spin matrix element
C-----------------------------------------------------------------------
C
            DZZAB(:,:,ITA,ITB,2) = (MZBZA_PO(:,:,1,ISMT)-MZBZA_PO(:,:,3,
     &                             ISMT))/SQRT_2
            DZZAB(:,:,ITA,ITB,3) = (MZBZA_PO(:,:,1,ISMT)+MZBZA_PO(:,:,3,
     &                             ISMT))*CI/SQRT_2
            DZZAB(:,:,ITA,ITB,4) = MZBZA_PO(:,:,2,ISMT)
C
            IF ( ITA.EQ.ITB ) THEN
               DZJAB(:,:,ITA,2) = (MZBJA_PO(:,:,1,ISMT)-MZBJA_PO(:,:,3,
     &                            ISMT))/SQRT_2
               DZJAB(:,:,ITA,3) = (MZBJA_PO(:,:,1,ISMT)+MZBJA_PO(:,:,3,
     &                            ISMT))*CI/SQRT_2
               DZJAB(:,:,ITA,4) = MZBJA_PO(:,:,2,ISMT)
C
            END IF
C
         END DO LOOP_ITB
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      END DO LOOP_ITA
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      END
