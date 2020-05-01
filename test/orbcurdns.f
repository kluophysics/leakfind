C*==orbcurdns.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE ORBCURDNS(EFCORRECT,RENORMALIZE,SHFTEF,SCLNOS,WEINP,
     &                     MSST,TAUT)
C   ********************************************************************
C   *                                                                  *
C   *            ROUTINE TO CALCULATE THE RELATIVISTIC                 *
C   *          PARAMAGNETIC CURRENT DENSITY CONTRIBUTION               *
C   *            FROM THE VALENCE ELECTRONS IN ATOMIC UNITS            *
C   *            IN FULL POTENTIAL FORMALISM                           *
C   *                                                                  *
C   *  JP = -1/PI * Im * (-i) * Int    TRACE GRAD(>) G(R,R,E)          *
C   *                         E=0..EF      - GRAD(<) G(R,R,E) dE       *
C   *                                                                  *
C   *         GRAD(< (>)) means gradient acting to the left (right)    *
C   *                                                                  *
C   *  adopted from T Huhne                                            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:NRMAX,R,JRWS,FULLPOT,JRCRI
      USE MOD_ANGMOM,ONLY:NL,NKM,NKMMAX,NLABIMAX,NCPLWF,L_IKM,IMKM_IKM,
     &    LB_IKM,MUEM05_IKM
      USE MOD_TYPES,ONLY:NTMAX,IMT,NCPLWFMAX,IKMCPLWF,ITBOT,ITTOP,
     &    JORB_LMCT
      USE MOD_FILES,ONLY:IFILCBWF,IOTMP,DATSET0,LDATSET0
      USE MOD_ENERGY,ONLY:NETAB
      USE MOD_CONSTANTS,ONLY:C0,CI,PI,C1,SQRT_2
      IMPLICIT NONE
C*--ORBCURDNS29
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*4 CURRSEL
      PARAMETER (CURRSEL='dwn')
      INTEGER NLAFPMAX
      PARAMETER (NLAFPMAX=8)
C
C Dummy arguments
C
      LOGICAL EFCORRECT,RENORMALIZE
      REAL*8 SCLNOS,SHFTEF
      COMPLEX*16 WEINP
      COMPLEX*16 MSST(NKMMAX,NKMMAX,NTMAX),TAUT(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 AF11,AF12,AF21,AF22,AG11,AG12,AG21,AG22,
     &           AMEA1(:,:,:,:,:,:),AMEA2(:,:,:,:,:,:),CTEST,CTRFO(3,3),
     &           CURRC1(:),CURRC2(:),CURRCH1(:),CURRCH2(:),
     &           JF(NRMAX,NCPLWFMAX,NKMMAX),JFP(:),
     &           JG(NRMAX,NCPLWFMAX,NKMMAX),JGP(:),TAU,TMAT(:,:),WE,
     &           ZF(NRMAX,NCPLWFMAX,NKMMAX),ZFJF1(:),ZFJF2(:),ZFP(:),
     &           ZFZF1(:),ZFZF2(:),ZG(NRMAX,NCPLWFMAX,NKMMAX),ZGJG1(:),
     &           ZGJG2(:),ZGP(:),ZGZG1(:),ZGZG2(:)
      REAL*8 CURR(:,:,:,:),JDNST(NRMAX,NLABIMAX,-1:+1,NTMAX),MJA,MJB
      CHARACTER*80 FILNAM
      INTEGER IA,IB,IC,ICURR,ICURRLM(:,:),ICURRLM1,ICURRLM2,IKMA,IKMB,
     &        IL0,IM,IMKMA,IMKMB,INDCURRLM(:,:),IPOL1,IPOL2,IR,IRTOP,
     &        ISPIN,IT,L0,L01MAX,L01MIN,L02MAX,L02MIN,LA,LAMA,LAMB,LB,
     &        LBARA,LBARB,LM0,LM01,LM02,LM0MAX(NTMAX),LMCURRI(:,:),LMR,
     &        LMR_MAX,LR_MAX,ML0,MPOL(3),MPOL1,NCURRLM(:),NCURRMAX
      LOGICAL INITIALIZE
      SAVE AMEA1,AMEA2,CURR,ICURRLM,INDCURRLM,LMCURRI,LMR_MAX,LR_MAX,
     &     NCURRLM,NCURRMAX
C
C*** End of declarations rewritten by SPAG
C
C select current contribution for spin component (up,dwn,tot)
C contributions are not calculated altogether at one time
C to save storage place
C
C
      ALLOCATABLE AMEA1,AMEA2,TMAT
      ALLOCATABLE CURRC1,CURRC2,CURRCH1,CURRCH2
      ALLOCATABLE JFP,ZFJF1,ZFJF2,ZFP,ZFZF1,ZFZF2
      ALLOCATABLE JGP,ZGJG1,ZGJG2,ZGP,ZGZG1,ZGZG2
      ALLOCATABLE CURR,ICURRLM,INDCURRLM,LMCURRI,NCURRLM
C
      DATA INITIALIZE/.TRUE./
C
      DATA MPOL/ + 1, - 1,0/
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
      IF ( INITIALIZE ) THEN
C
         LR_MAX = 2*(NL-1) + 1
         LMR_MAX = (LR_MAX+1)**2
         NCURRMAX = LMR_MAX
C
         ALLOCATE (CURR(NRMAX,NCURRMAX,NTMAX,3))
         ALLOCATE (ICURRLM(LMR_MAX,NTMAX),INDCURRLM(LMR_MAX,NTMAX))
         ALLOCATE (LMCURRI(NCURRMAX,NTMAX),NCURRLM(NTMAX))
C
         CURR(:,:,:,:) = 0D0
C
C         ALLOCATE (JORB_LMCT(NRMAX,LMR_MAX,3,NTMAX))
C         JORB_LMCT(:,:,:,:) = 0D0
C
         DO IT = 1,NTMAX
            NCURRLM(IT) = 0
            DO IB = 1,LMR_MAX
               INDCURRLM(IB,IT) = 0
               ICURRLM(IB,IT) = 0
            END DO
            DO IB = 1,NCURRMAX
               LMCURRI(IB,IT) = 0
            END DO
         END DO
C
C calculate angular parts of derivatives
C
         ALLOCATE (AMEA1(NKMMAX,NKMMAX,NLAFPMAX,3,2,2))
         ALLOCATE (AMEA2(NKMMAX,NKMMAX,NLAFPMAX,3,2,2))
C
         CALL AMECURRFP(AMEA1,AMEA2,NLAFPMAX)
C
         INITIALIZE = .FALSE.
C
      END IF
C=======================================================================
C
      ALLOCATE (TMAT(NKMMAX,NKMMAX))
      ALLOCATE (CURRC1(NRMAX),CURRC2(NRMAX),CURRCH1(NRMAX))
      ALLOCATE (CURRCH2(NRMAX),JFP(NRMAX),JGP(NRMAX),ZFJF1(NRMAX))
      ALLOCATE (ZFJF2(NRMAX),ZFP(NRMAX),ZFZF1(NRMAX),ZFZF2(NRMAX))
      ALLOCATE (ZGJG1(NRMAX),ZGJG2(NRMAX))
      ALLOCATE (ZGP(NRMAX),ZGZG1(NRMAX),ZGZG2(NRMAX))
C
C transformation gradient from spherical basis to
C rectangular basis (end of calculation)
C
      CALL CINIT(3*3,CTRFO)
C
      CTRFO(1,1) = 1/SQRT_2
      CTRFO(2,1) = 1/SQRT_2
C
      CTRFO(1,2) = +1/(SQRT_2*CI)
      CTRFO(2,2) = -1/(SQRT_2*CI)
C
      CTRFO(3,3) = C1
C
C select spin - resolved current density - component
C
      IF ( CURRSEL.EQ.'up' ) THEN
         ISPIN = 2
      ELSE
         ISPIN = 1
      END IF
C
C ======================================================================
C           renormalize current density JDNST
C           either after E-loop over semi core band
C           or after valence band E-loop with proper
C           charge determined by Lloyd formula
C ======================================================================
C
      IF ( RENORMALIZE ) THEN
C
         DO IT = ITBOT,ITTOP
            IM = IMT(IT)
            IRTOP = JRCRI(IM)
C
            JDNST(1:IRTOP,1:NLABIMAX,-1:+1,IT)
     &         = JDNST(1:IRTOP,1:NLABIMAX,-1:+1,IT)*SCLNOS
         END DO
C
      END IF
C
C ======================================================================
C
      IF ( EFCORRECT ) THEN
         WE = DCMPLX(SHFTEF,0.0D0)
      ELSE
         WE = WEINP
      END IF
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      DO IT = ITBOT,ITTOP
C
C--------------------- site diagonal multiple scattering matrix G or TAU
C
         CALL TAUGFCONV(MSST(1,1,IT),TAUT(1,1,IT),TMAT)
C
         IM = IMT(IT)
         IRTOP = JRCRI(IM)
C
         CALL WAVFUN_READ_REL(IFILCBWF,IT,1,ZG,ZF,JG,JF,IRTOP,NCPLWF,
     &                        IKMCPLWF)
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
         LOOP_LAMB:DO LAMB = 1,NKM
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            LOOP_LAMA:DO LAMA = 1,NKM
C
               IF ( (LAMB.EQ.LAMA) .OR. 
     &              (ABS(TMAT(LAMB,LAMA)).GT.1.0D-12) ) THEN
C
                  TAU = TMAT(LAMB,LAMA)
C
                  LM0MAX(IT) = 0
C
C-------------------------------------------------------------------- IB
                  LOOP_IB:DO IB = 1,NCPLWF(LAMB)
                     IKMB = IKMCPLWF(IB,LAMB)
                     IMKMB = IMKM_IKM(IKMB)
C
                     LB = L_IKM(IKMB)
                     LBARB = LB_IKM(IKMB)
                     MJB = MUEM05_IKM(IKMB) + 0.5D0
C
C-------------------------------------------------------------------- IA
                     LOOP_IA:DO IA = 1,NCPLWF(LAMA)
C
                        IKMA = IKMCPLWF(IA,LAMA)
                        IMKMA = IMKM_IKM(IKMA)
C
                        LA = L_IKM(IKMA)
                        LBARA = LB_IKM(IKMA)
                        MJA = MUEM05_IKM(IKMA) + 0.5D0
C
C calculate radial parts of derivatives
C
                        CALL CDIFFER(IM,ZG(1,IA,LAMA),ZGP)
C
                        CALL CDIFFER(IM,ZF(1,IA,LAMA),ZFP)
C
                        CALL CDIFFER(IM,JG(1,IA,LAMA),JGP)
C
                        CALL CDIFFER(IM,JF(1,IA,LAMA),JFP)
C
C
                        DO IR = 1,IRTOP
C
                           ZGZG1(IR) = ZG(IR,IB,LAMB)
     &                                 *(ZGP(IR)-LA/R(IR,IM)
     &                                 *ZG(IR,IA,LAMA))
C
                           ZGZG2(IR) = ZG(IR,IB,LAMB)
     &                                 *(ZGP(IR)+(LA+1)/R(IR,IM)
     &                                 *ZG(IR,IA,LAMA))
C
                           ZFZF1(IR) = -ZF(IR,IB,LAMB)
     &                                 *(ZFP(IR)-LBARA/R(IR,IM)
     &                                 *ZF(IR,IA,LAMA))
C
                           ZFZF2(IR) = -ZF(IR,IB,LAMB)
     &                                 *(ZFP(IR)+(LBARA+1)/R(IR,IM)
     &                                 *ZF(IR,IA,LAMA))
C
                           ZGJG1(IR) = 0.5D0*(ZG(IR,IB,LAMB)*(JGP(IR)-LA
     &                                 /R(IR,IM)*JG(IR,IA,LAMA))
     &                                 +JG(IR,IB,LAMB)
     &                                 *(ZGP(IR)-LA/R(IR,IM)
     &                                 *ZG(IR,IA,LAMA)))
C
                           ZGJG2(IR) = 0.5D0*(ZG(IR,IB,LAMB)*(JGP(IR)+(
     &                                 LA+1)/R(IR,IM)*JG(IR,IA,LAMA))
     &                                 +JG(IR,IB,LAMB)
     &                                 *(ZGP(IR)+(LA+1)/R(IR,IM)
     &                                 *ZG(IR,IA,LAMA)))
C
                           ZFJF1(IR) = -0.5D0*(ZF(IR,IB,LAMB)*(JFP(IR)-
     &                                 LBARA/R(IR,IM)*JF(IR,IA,LAMA))
     &                                 +JF(IR,IB,LAMB)
     &                                 *(ZFP(IR)-LBARA/R(IR,IM)
     &                                 *JF(IR,IA,LAMA)))
C
                           ZFJF2(IR) = -0.5D0*(ZF(IR,IB,LAMB)*(JFP(IR)+(
     &                                 LBARA+1)/R(IR,IM)*JF(IR,IA,LAMA))
     &                                 +JF(IR,IB,LAMB)
     &                                 *(ZFP(IR)+(LBARA+1)/R(IR,IM)
     &                                 *ZF(IR,IA,LAMA)))
C
                        END DO
C
C
                        DO IPOL1 = 1,3
C
                           MPOL1 = MPOL(IPOL1)
C
                           ML0 = INT(MJB-MJA-MPOL1)
C
                           L01MIN = MIN(ABS(LB-(LA+1)),
     &                              ABS(LBARB-(LBARA+1)))
C
                           L02MIN = MIN(ABS(LB-(LA-1)),
     &                              ABS(LBARB-(LBARA-1)))
C
                           L01MAX = MAX(ABS(LB+(LA+1)),
     &                              ABS(LBARB+(LBARA+1)))
C
                           L02MAX = MAX(ABS(LB+(LA-1)),
     &                              ABS(LBARB+(LBARA-1)))
C
                           DO L0 = MIN(L01MIN,L02MIN),MAX(L01MAX,L02MAX)
C
                              IF ( ABS(ML0).LE.L0 ) THEN
                                 IL0 = L0 + 1
                                 LM0 = L0**2 + L0 + 1 + ML0
                                 IF ( LM0.GT.LM0MAX(IT) ) LM0MAX = LM0
C
                                 LM01 = L0**2 + L0 + 1 - ABS(ML0)
                                 LM02 = L0**2 + L0 + 1 + ABS(ML0)
C
                                 CTEST = C0
C
                                 AG11 = AMEA1(IKMB,IKMA,IL0,IPOL1,1,
     &                                  ISPIN)
                                 AG21 = AMEA2(IKMB,IKMA,IL0,IPOL1,1,
     &                                  ISPIN)
                                 AG12 = AMEA1(IKMB,IKMA,IL0,IPOL1,2,
     &                                  ISPIN)
                                 AG22 = AMEA2(IKMB,IKMA,IL0,IPOL1,2,
     &                                  ISPIN)
C
                                 IF ( CURRSEL.EQ.'tot' ) THEN
                                    AG11 = AG11 + 
     &                                 AMEA1(IKMB,IKMA,IL0,IPOL1,1,
     &                                 -ISPIN+3)
                                    AG21 = AG21 + 
     &                                 AMEA2(IKMB,IKMA,IL0,IPOL1,1,
     &                                 -ISPIN+3)
                                    AG12 = AG12 + 
     &                                 AMEA1(IKMB,IKMA,IL0,IPOL1,2,
     &                                 -ISPIN+3)
                                    AG22 = AG22 + 
     &                                 AMEA2(IKMB,IKMA,IL0,IPOL1,2,
     &                                 -ISPIN+3)
                                 END IF
C
                                 IF ( (IMKMB.LE.NKMMAX) .AND. 
     &                                (IMKMA.LE.NKMMAX) ) THEN
C
                                    AF11 = AMEA1(IMKMB,IMKMA,IL0,IPOL1,
     &                                 1,-ISPIN+3)
                                    AF21 = AMEA2(IMKMB,IMKMA,IL0,IPOL1,
     &                                 1,-ISPIN+3)
                                    AF12 = AMEA1(IMKMB,IMKMA,IL0,IPOL1,
     &                                 2,-ISPIN+3)
                                    AF22 = AMEA2(IMKMB,IMKMA,IL0,IPOL1,
     &                                 2,-ISPIN+3)
C
                                    IF ( CURRSEL.EQ.'tot' ) THEN
C
                                       AF11 = AF11 + 
     &                                    AMEA1(IMKMB,IMKMA,IL0,IPOL1,1,
     &                                    -ISPIN+3)
                                       AF21 = AF21 + 
     &                                    AMEA2(IMKMB,IMKMA,IL0,IPOL1,1,
     &                                    -ISPIN+3)
                                       AF12 = AF12 + 
     &                                    AMEA1(IMKMB,IMKMA,IL0,IPOL1,2,
     &                                    -ISPIN+3)
                                       AF22 = AF22 + 
     &                                    AMEA2(IMKMB,IMKMA,IL0,IPOL1,2,
     &                                    -ISPIN+3)
C
                                    END IF
C
                                 ELSE
C
                                    AF11 = C0
                                    AF21 = C0
                                    AF12 = C0
                                    AF22 = C0
C
                                 END IF
C
                                 DO IR = 1,IRTOP
C
C "complex" current density, regular part
C  extension "H" means parts with gradient acting to the left
C
                                    CURRC1(IR)
     &                                 = TAU*WE*1/CI*(ZGZG1(IR)*AG11-
     &                                 ZGZG2(IR)*AG21+ZFZF1(IR)
     &                                 *AF11-ZFZF2(IR)*AF21)
C
                                    CURRCH1(IR)
     &                                 = TAU*WE*1/CI*(ZGZG1(IR)*DCONJG
     &                                 (AG11)-ZGZG2(IR)*DCONJG(AG21)
     &                                 +ZFZF1(IR)*DCONJG(AF11)-ZFZF2(IR)
     &                                 *DCONJG(AF21))
C
                                    CURRC2(IR)
     &                                 = TAU*WE*1/CI*(ZGZG1(IR)*AG12-
     &                                 ZGZG2(IR)*AG22+ZFZF1(IR)
     &                                 *AF12-ZFZF2(IR)*AF22)
C
                                    CURRCH2(IR)
     &                                 = TAU*WE*1/CI*(ZGZG1(IR)*DCONJG
     &                                 (AG12)-ZGZG2(IR)*DCONJG(AG22)
     &                                 +ZFZF1(IR)*DCONJG(AF12)-ZFZF2(IR)
     &                                 *DCONJG(AF22))
C
C irregular part
C
                                    IF ( LAMB.EQ.LAMA ) THEN
C
                                       CURRCH1(IR) = CURRCH1(IR)
     &                                    - WE*1/CI*(ZGJG1(IR)
     &                                    *DCONJG(AG11)-ZGJG2(IR)
     &                                    *DCONJG(AG21)+ZFJF1(IR)
     &                                    *DCONJG(AF11)-ZFJF2(IR)
     &                                    *DCONJG(AF21))
C
                                       CURRC1(IR) = CURRC1(IR)
     &                                    - WE*1/CI*(ZGJG1(IR)
     &                                    *AG11-ZGJG2(IR)*AG21+ZFJF1(IR)
     &                                    *AF11-ZFJF2(IR)*AF21)
C
C
                                       CURRC2(IR) = CURRC2(IR)
     &                                    - WE*1/CI*(ZGJG1(IR)
     &                                    *AG12-ZGJG2(IR)*AG22+ZFJF1(IR)
     &                                    *AF12-ZFJF2(IR)*AF22)
C
                                       CURRCH2(IR) = CURRCH2(IR)
     &                                    - WE*1/CI*(ZGJG1(IR)
     &                                    *DCONJG(AG12)-ZGJG2(IR)
     &                                    *DCONJG(AG22)+ZFJF1(IR)
     &                                    *DCONJG(AF12)-ZFJF2(IR)
     &                                    *DCONJG(AF22))
C
                                    END IF
C
                                    CTEST = CTEST + CURRC1(IR)
     &                                 + CURRCH1(IR)
C
                                    IF ( LM01.NE.LM02 ) CTEST = CTEST + 
     &                                 CURRC2(IR) + CURRCH2(IR)
C
                                 END DO
C
                                 IF ( ABS(CTEST).GT.1.0D-10 ) THEN
C
                                    IF ( INDCURRLM(LM0,IT).EQ.0 ) THEN
                                       INDCURRLM(LM01,IT) = 1
                                       INDCURRLM(LM02,IT) = 1
                                       IF ( LM02.NE.LM01 ) THEN
                                         ICURRLM(LM01,IT) = NCURRLM(IT)
     &                                      + 1
                                         ICURRLM(LM02,IT) = NCURRLM(IT)
     &                                      + 2
                                         LMCURRI(ICURRLM(LM01,IT),IT)
     &                                      = LM01
                                         LMCURRI(ICURRLM(LM02,IT),IT)
     &                                      = LM02
                                         NCURRLM(IT) = NCURRLM(IT) + 2
                                       ELSE
                                         ICURRLM(LM01,IT) = NCURRLM(IT)
     &                                      + 1
                                         NCURRLM(IT) = NCURRLM(IT) + 1
                                         LMCURRI(ICURRLM(LM01,IT),IT)
     &                                      = LM01
                                       END IF
                                    END IF
C
                                    ICURRLM1 = ICURRLM(LM01,IT)
                                    ICURRLM2 = ICURRLM(LM02,IT)
C
C transform derivatives from angular- to rectangular basis:
C d+,d-,d0 -> dx,dy,dz and take -1/pi*imaginary part of
C formerly calculated complex terms to get correct current density
C
                                    DO IPOL2 = 1,3
                                       DO IR = 1,IRTOP
C
                                         CURR(IR,ICURRLM1,IT,IPOL2)
     &                                      = CURR(IR,ICURRLM1,IT,IPOL2)
     &                                      - 
     &                                      1/PI*DIMAG(CTRFO(IPOL1,IPOL2
     &                                      )*CURRC1(IR)
     &                                      -DCONJG(CTRFO(IPOL1,IPOL2))
     &                                      *CURRCH1(IR))
C
                                         IF ( ICURRLM1.NE.ICURRLM2 )
     &                                      CURR(IR,ICURRLM2,IT,IPOL2)
     &                                      = CURR(IR,ICURRLM2,IT,IPOL2)
     &                                      - 
     &                                      1/PI*DIMAG(CTRFO(IPOL1,IPOL2
     &                                      )*CURRC2(IR)
     &                                      -DCONJG(CTRFO(IPOL1,IPOL2))
     &                                      *CURRCH2(IR))
C
                                       END DO
                                    END DO
                                 END IF
C
                              END IF
C     l0
                           END DO
C     ipol1
                        END DO
C
                     END DO LOOP_IA
C-----------------------------------------------------------------------
                  END DO LOOP_IB
C-----------------------------------------------------------------------
C
               END IF
C
            END DO LOOP_LAMA
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
         END DO LOOP_LAMB
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
C-----------------------------------------------------------------------
      IF ( .NOT.EFCORRECT ) RETURN
C-----------------------------------------------------------------------
C
C JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
C   *      JDNST(IR,ILR,M.IT) = j^c_{lm}(r)                            *
C   *                                                                  *
C   * For rotational symmetry one has                                  *
C   *                                                                  *
C   * c=x  m=-1     l=1,3,5,...   ILR=(l+1)/2                          *
C   * c=y  m=+1                                                        *
C JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
C
C
      ALLOCATE (JORB_LMCT(NRMAX,LMR_MAX,3,NTMAX))
      JORB_LMCT(:,:,:,:) = 0D0
C
      FILNAM = DATSET0(1:LDATSET0)//'_JORB.dat'
C
      CALL WRHEAD(IOTMP,FILNAM,'J_orb',NETAB(1))
C
      DO IT = ITBOT,ITTOP
         IM = IMT(IT)
         IF ( FULLPOT ) THEN
            IRTOP = JRCRI(IM)
         ELSE
            IRTOP = JRWS(IM)
         END IF
         DO ICURR = 1,NCURRLM(IT)
            LMR = LMCURRI(ICURR,IT)
            DO IC = 1,3
               DO IR = 1,IRTOP
                  JORB_LMCT(IR,LMR,IC,IT) = CURR(IR,ICURR,IT,IC)
               END DO
            END DO
         END DO
C
         WRITE (IOTMP,'(''TYPE IT'',/,4I10)') IT,IRTOP,LR_MAX,LMR_MAX
         DO LMR = 1,LMR_MAX
            DO IR = 1,IRTOP
               WRITE (IOTMP,'(3E25.10)')
     &                (JORB_LMCT(IR,LMR,IC,IT),IC=1,3)
            END DO
         END DO
C
      END DO
C
      CLOSE (IOTMP)
C
      END
C*==amecurrfp.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE AMECURRFP(AMEA1,AMEA2,NLAFPMAX)
C   ********************************************************************
C   *                                                                  *
C   *   calculate the angular matrix elements needed                   *
C   *                                                                  *
C   *   for the          current - density - calculation               *
C   *                                                                  *
C   *                       ->   ->                                    *
C   *          A1 = < LAM| NAB . e(POL) | LAM'>                        *
C   *          A2 = < LAM| NAB . e(POL) | LAM'>                        *
C   *                                                                  *
C   *   ipol= 1,2,3  ==  (+),(-),(z)                                   *
C   *                                                                  *
C   *  adopted from T Huhne                                            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_ANGMOM,ONLY:NKMMAX,CGC,NLMAX
      IMPLICIT NONE
C*--AMECURRFP602
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NLAFPMAX
      COMPLEX*16 AMEA1(NKMMAX,NKMMAX,NLAFPMAX,3,2,2),
     &           AMEA2(NKMMAX,NKMMAX,NLAFPMAX,3,2,2)
C
C Local variables
C
      COMPLEX*16 AME1C,AME2C,YR(2)
      REAL*8 CGC1,GAUNT_CYLM
      INTEGER IKM1,IKM2,ILA,IPOL1,J1P05,J2P05,K1,K2,KAP1,KAP2,L1,L2,LA,
     &        M1A,M1B,M2A,M2B,MA,MM,MPOL(3),MUE1M05,MUE2M05,NK
      REAL*8 SGN(-1:+1)
C
C*** End of declarations rewritten by SPAG
C
      DATA SGN/ + 1.0D0, + 1.0D0, - 1.0D0/
      DATA MPOL/ + 1, - 1,0/
C
      NK = 2*NLMAX - 1
C
      AMEA1(:,:,:,:,:,:) = C0
      AMEA2(:,:,:,:,:,:) = C0
C
      IKM1 = 0
      DO K1 = 1,NK
         L1 = K1/2
         IF ( MOD(K1,2).EQ.0 ) THEN
            KAP1 = L1
         ELSE
            KAP1 = -L1 - 1
         END IF
         J1P05 = IABS(KAP1)
C
         DO MUE1M05 = -J1P05,J1P05 - 1
            IKM1 = IKM1 + 1
C
            IKM2 = 0
            DO K2 = 1,NK
               L2 = K2/2
               IF ( MOD(K2,2).EQ.0 ) THEN
                  KAP2 = L2
               ELSE
                  KAP2 = -L2 - 1
               END IF
               J2P05 = IABS(KAP2)
C
               DO MUE2M05 = -J2P05,J2P05 - 1
                  IKM2 = IKM2 + 1
C
C-----------------------------------------------------------------------
C
                  DO IPOL1 = 1,3
C
                     MM = MPOL(IPOL1)
C
                     M1A = MUE1M05 + 1
                     M2A = MUE2M05 + 1 + MM
C
                     M1B = MUE1M05
                     M2B = MUE2M05 + MM
C
                     DO LA = 0,(2*NLMAX-1)
C
                        ILA = LA + 1
                        MA = MUE1M05 - MUE2M05 - MM
C
                        IF ( ABS(MA).LE.ABS(LA) ) THEN
C
                           AME1C = SGN(MM)*SQRT(DBLE(L2+1)/DBLE(2*L2+3))
     &                             *DCMPLX
     &                             ((CGC(IKM1,1)*CGC(IKM2,1)*CGC1(L2,
     &                             L2+1,(MUE2M05+1),MM)
     &                             *GAUNT_CYLM(L1,M1A,LA,MA,(L2+1),M2A))
     &                             )
C
                           AME2C = SGN(MM)*SQRT(DBLE(L2+1)/DBLE(2*L2+3))
     &                             *DCMPLX
     &                             ((CGC(IKM1,2)*CGC(IKM2,2)*CGC1(L2,
     &                             L2+1,(MUE2M05+0),MM)
     &                             *GAUNT_CYLM(L1,M1B,LA,MA,(L2+1),M2B))
     &                             )
C
                           CALL YCR(AME1C,YR,MA)
C
                           AMEA1(IKM1,IKM2,ILA,IPOL1,1,1) = YR(1)
                           AMEA1(IKM1,IKM2,ILA,IPOL1,2,1) = YR(2)
C
                           CALL YCR(AME2C,YR,MA)
C
                           AMEA1(IKM1,IKM2,ILA,IPOL1,1,2) = YR(1)
                           AMEA1(IKM1,IKM2,ILA,IPOL1,2,2) = YR(2)
C
C
                           IF ( L2.NE.0 ) THEN
C
                              AME1C = SGN(MM)
     &                                *SQRT(DBLE(L2)/DBLE(2*L2-1))
     &                                *DCMPLX
     &                                ((CGC(IKM1,1)*CGC(IKM2,1)*CGC1(L2,
     &                                L2-1,(MUE2M05+1),MM)
     &                                *GAUNT_CYLM(L1,M1A,LA,MA,(L2-1),
     &                                M2A)))
C
                              AME2C = SGN(MM)
     &                                *SQRT(DBLE(L2)/DBLE(2*L2-1))
     &                                *DCMPLX
     &                                ((+CGC(IKM1,2)*CGC(IKM2,2)*CGC1
     &                                (L2,L2-1,(MUE2M05+0),MM)
     &                                *GAUNT_CYLM(L1,M1B,LA,MA,(L2-1),
     &                                M2B)))
C
                              CALL YCR(AME1C,YR,MA)
C
                              AMEA2(IKM1,IKM2,ILA,IPOL1,1,1) = YR(1)
                              AMEA2(IKM1,IKM2,ILA,IPOL1,2,1) = YR(2)
C
                              CALL YCR(AME2C,YR,MA)
C
                              AMEA2(IKM1,IKM2,ILA,IPOL1,1,2) = YR(1)
                              AMEA2(IKM1,IKM2,ILA,IPOL1,2,2) = YR(2)
C
                           END IF
C
                        END IF
                     END DO
                  END DO
C
C-----------------------------------------------------------------------
C
               END DO
            END DO
         END DO
      END DO
C
      END
C*==ycr.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE YCR(YC,YR,M)
C   ********************************************************************
C   *                                                                  *
C   *   routine to perform transformation from                         *
C   *   complex spherical harmonic representation                      *
C   *   to real spherical harmonic representation                      *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_CONSTANTS,ONLY:C0,CI,C1,SQRT_2
      IMPLICIT NONE
C*--YCR766
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER M
      COMPLEX*16 YC
      COMPLEX*16 YR(2)
C
C Local variables
C
      COMPLEX*16 FAC,FACI
      INTEGER MABS
      REAL*8 RFAC
C
C*** End of declarations rewritten by SPAG
C
      RFAC = 1.0D0
C
      FACI = CI/SQRT_2*RFAC
      FAC = C1/SQRT_2*RFAC
C
      MABS = ABS(M)
      IF ( MABS.EQ.0 ) THEN
         YR(1) = YC*RFAC
         YR(2) = C0
      ELSE IF ( MABS.EQ.M ) THEN
         YR(1) = FACI*DCMPLX((-1)**(MABS+1))*YC
         YR(2) = FAC*DCMPLX((-1)**MABS)*YC
      ELSE IF ( MABS.EQ.-M ) THEN
         YR(1) = FACI*YC
         YR(2) = FAC*YC
      END IF
C
      END
