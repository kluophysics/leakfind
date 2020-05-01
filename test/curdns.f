C*==curdns.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CURDNS(EFCORRECT,RENORMALIZE,SHFTEF,SCLNOS,JDNST,WEINP,
     &                  TAUT)
C   ********************************************************************
C   *                                                                  *
C   * To avoid problems with complex expansion coefficients the vector *
C   * field ->j(r) is expanded w.r.t. REAL spherical harmonics         *
C   *                                                                  *
C   * 1st step: the expansion coefficients w.r.t. to COMPLEX spherical *
C   * harmonics are calculated for spherical coordiantes M:            *
C   * j^{M}_{lm}(r)                                                    *
C   *                                                                  *
C   * 2nd step:  the expansion coefficients are transferred to refer   *
C   * to REAL spherical harmonics and cartesian coordinates            *
C   *                                                                  *
C   *      JDNST(IR,ILR,M.IT) = j^c_{lm}(r)                            *
C   *                                                                  *
C   * For rotational symmetry one has                                  *
C   *                                                                  *
C   * c=x  m=-1     l=1,3,5,...   ILR=(l+1)/2                          *
C   * c=y  m=+1                                                        *
C   *                                                                  *
C   * For the following use is made of the relation between expansion  *
C   * coefficients w.r.t REAL and COMPLEX spherical harmonics.         *
C   *                                                                  *
C   * j^{+}_{l,+1} = -j^{-}_{l,-1} = -i*j^{x}_{l,-1} = i*j^{y}_{l,+1}  *
C   *                                                                  *
C   *                                                                  *
C   * NOTE: concerning the contributions of the irregular              *
C   *       wavefunctions the average of JFZG and ZFJG is taken        *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ENERGY,ONLY:NETAB
      USE MOD_CALCMODE,ONLY:LHS_SOL_EQ_RHS_SOL
      USE MOD_FILES,ONLY:IFILCBWF,LDATSET0,DATSET0,IOTMP
      USE MOD_RMESH,ONLY:JRWS,JRCRI,RWS,NRMAX,R2DRDI_W_RADINT,R,DRDI,
     &    FULLPOT
      USE MOD_TYPES,ONLY:IMT,NTMAX,ITBOT,ITTOP,IKMCPLWF,NCPLWFMAX,
     &    JALF_LMCT
      USE MOD_ANGMOM,ONLY:NL,NKM,NLABIMAX,NKMMAX,AMEBI1,AMEBI2,NCPLWF
      USE MOD_CONSTANTS,ONLY:PI,A0_CGS,M0_CGS,C_CGS,HBAR_CGS,CI,C0,
     &    SQRT_2,C_AU
      IMPLICIT NONE
C*--CURDNS45
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CURDNS')
C
C Dummy arguments
C
      LOGICAL EFCORRECT,RENORMALIZE
      REAL*8 SCLNOS,SHFTEF
      COMPLEX*16 WEINP
      REAL*8 JDNST(NRMAX,NLABIMAX,-1:+1,NTMAX)
      COMPLEX*16 TAUT(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      REAL*8 AFG,AGF,APHI1(:,:),APHIAUX(:),APHIL(:,:),DAPHI,DRA,
     &       DTETADEG,JPHI,MAGMOMT(:),MUEA,MUEJ,PRE,RA,RHAT(3),RINT(:),
     &       RJ,RSUM,TETA,UNITS,WGAU(:),XGAU(:),Y_JGRID(:,:)
      COMPLEX*16 APHI2(:,:),DA(:,:,:),DB(:,:,:,:),DC(:,:,:,:),JFL(:,:,:)
     &           ,JFR(:,:,:),JFZG(:),JGL(:,:,:),JGR(:,:,:),JGZF(:),WE,
     &           WJFZG,WJGZF,WM,WP,WX,WY,WZFJG,WZFZG,WZGJF,WZGZF,WZJ,
     &           WZZ,ZFJG(:),ZFL(:,:,:),ZFR(:,:,:),ZFZG(:),ZGJF(:),
     &           ZGL(:,:,:),ZGR(:,:,:),ZGZF(:)
      LOGICAL CALCVECPOT
      CHARACTER*80 FILNAM
      INTEGER IA,IA_ERR,IB,IC,IFIL_LHSB,IFIL_RHSA,IKMA,IKMB,ILR,IM,
     &        IMUEA,IMUEJ,IR,IRA,IRTOP,IT,ITETA,LAMA,LAMB,LMA,LMAX,LMJ,
     &        LMR,LMR_M1,LMR_MAX,LMR_P1,LR,LR_MAX,M,NLMR,NLR,NMUEA,
     &        NMUEJ,NRA,NTETA
      CHARACTER*60 OUTPUT
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE ZFR,ZGR,JFR,JGR,ZFL,ZGL,JFL,JGL
      ALLOCATABLE DA,DB,DC,MAGMOMT,ZGZF,ZFZG,ZGJF,ZFJG
      ALLOCATABLE APHI1,APHI2,XGAU,WGAU,RINT,Y_JGRID,APHIAUX,APHIL
      ALLOCATABLE JFZG,JGZF
C
      ALLOCATE (DA(NRMAX,NLABIMAX,-1:+1))
      ALLOCATE (DB(NRMAX,NLABIMAX,-1:+1,3))
      ALLOCATE (DC(NRMAX,NLABIMAX,-1:+1,3),MAGMOMT(NTMAX))
      ALLOCATE (JGR(NRMAX,NCPLWFMAX,NKM),JFR(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (ZGR(NRMAX,NCPLWFMAX,NKM),ZFR(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (JGL(NRMAX,NCPLWFMAX,NKM),JFL(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (ZGL(NRMAX,NCPLWFMAX,NKM),ZFL(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (ZGZF(NRMAX),ZFZG(NRMAX))
      ALLOCATE (ZGJF(NRMAX),ZFJG(NRMAX))
      ALLOCATE (JFZG(NRMAX),JGZF(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: JFZG')
C
      CALCVECPOT = .FALSE.
C
      IFIL_RHSA = IFILCBWF
      IF ( LHS_SOL_EQ_RHS_SOL ) THEN
         IFIL_LHSB = IFILCBWF
      ELSE
         IFIL_LHSB = IFILCBWF + 1
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
            IRTOP = JRWS(IM)
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
      WZJ = -CI*WE/PI
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = ITBOT,ITTOP
         IM = IMT(IT)
         IF ( FULLPOT ) THEN
            IRTOP = JRCRI(IM)
         ELSE
            IRTOP = JRWS(IM)
         END IF
C
         LMAX = NL - 1
         DA(1:NRMAX,1:NLABIMAX,-1:+1) = C0
C
C ======================================================================
C     1st step: calculate current density profile
C               w.r.t  sigma_(+,0,-) and COMPLEX spherical harmonics
C ======================================================================
C
C ------------------------------------------ read in wavefunctions for B
C
         CALL WAVFUN_READ_REL(IFIL_LHSB,IT,1,ZGL,ZFL,JGL,JFL,IRTOP,
     &                        NCPLWF,IKMCPLWF)
C
C ------------------------------------------ read in wavefunctions for A
C
         CALL WAVFUN_READ_REL(IFIL_RHSA,IT,1,ZGR,ZFR,JGR,JFR,IRTOP,
     &                        NCPLWF,IKMCPLWF)
C
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
         LOOP_LAMB:DO LAMB = 1,NKM
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            LOOP_LAMA:DO LAMA = 1,NKM
C
               WZZ = -CI*WE*TAUT(LAMA,LAMB,IT)/PI
C
C-------------------------------------------------------------------- IB
               LOOP_IB:DO IB = 1,NCPLWF(LAMB)
                  IKMB = IKMCPLWF(IB,LAMB)
C
C-------------------------------------------------------------------- IA
                  LOOP_IA:DO IA = 1,NCPLWF(LAMA)
                     IKMA = IKMCPLWF(IA,LAMA)
C
C---------------------------------------------------- set up radial part
                     DO IR = 1,IRTOP
                        ZGZF(IR) = ZGL(IR,IB,LAMB)*ZFR(IR,IA,LAMA)
                        ZFZG(IR) = ZFL(IR,IB,LAMB)*ZGR(IR,IA,LAMA)
C
                        JGZF(IR) = JGL(IR,IB,LAMB)*ZFR(IR,IA,LAMA)
                        JFZG(IR) = JFL(IR,IB,LAMB)*ZGR(IR,IA,LAMA)
C
                        ZGJF(IR) = ZGL(IR,IB,LAMB)*JFR(IR,IA,LAMA)
                        ZFJG(IR) = ZFL(IR,IB,LAMB)*JGR(IR,IA,LAMA)
C
                     END DO
C
C-------------------------------------------------------------------- LR
                     LOOP_LR:DO LR = 1,(2*LMAX+1),2
                        ILR = (LR+1)/2
C
C--------------------------------------------------------------------- M
                        LOOP_M:DO M = -1, + 1
C
                           AGF = AMEBI1(IKMB,IKMA,ILR,-M)
                           AFG = AMEBI2(IKMB,IKMA,ILR,-M)
C
C------------------------------------------------- check selection rules
C
                           IF ( ABS(AGF).GT.1D-8 .OR. ABS(AFG).GT.1D-8 )
     &                          THEN
C
                              WZGZF = WZZ*AGF
                              WZFZG = WZZ*AFG
C
                              DO IR = 1,IRTOP
                                 DA(IR,ILR,M) = DA(IR,ILR,M)
     &                              + WZGZF*ZGZF(IR) - WZFZG*ZFZG(IR)
                              END DO
C
C----------------- no irregular contributions to the backscattering term
C
                              IF ( LAMB.EQ.LAMA ) THEN
C
                                 WZGJF = WZJ*AGF
                                 WZFJG = WZJ*AFG
C
                                 WJGZF = WZJ*AGF
                                 WJFZG = WZJ*AFG
C
                                 DO IR = 1,IRTOP
                                    DA(IR,ILR,M) = DA(IR,ILR,M)
     &                                 - 0.5D0*(WJGZF*JGZF(IR)
     &                                 -WJFZG*JFZG(IR))
C
                                    DA(IR,ILR,M) = DA(IR,ILR,M)
     &                                 - 0.5D0*(WZGJF*ZGJF(IR)
     &                                 -WZFJG*ZFJG(IR))
                                 END DO
C
                              END IF
C
                           END IF
C
                        END DO LOOP_M
C-----------------------------------------------------------------------
                     END DO LOOP_LR
C-----------------------------------------------------------------------
                  END DO LOOP_IA
C-----------------------------------------------------------------------
               END DO LOOP_IB
C-----------------------------------------------------------------------
            END DO LOOP_LAMA
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
         END DO LOOP_LAMB
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
C
C ======================================================================
C     2nd step: convert current density profile
C               w.r.t  sigma_(x,y,z) and REAL spherical harmonics
C ======================================================================
C
         DB(1:NRMAX,1:NLABIMAX,-1:+1,1:3) = C0
         WX = -1D0/SQRT_2
         WY = CI/SQRT_2
         DO LR = 1,(2*LMAX+1),2
            ILR = (LR+1)/2
            DO IR = 1,IRTOP
               DB(IR,ILR,-1,1) = +WX*(DA(IR,ILR,+1))
               DB(IR,ILR,+1,1) = -WX*(DA(IR,ILR,-1))
C
               DB(IR,ILR,-1,2) = +WY*DA(IR,ILR,+1)
               DB(IR,ILR,+1,2) = +WY*DA(IR,ILR,-1)
            END DO
         END DO
C
C
         DC(1:NRMAX,1:NLABIMAX,-1:+1,1:3) = C0
         WP = CI/SQRT_2
         WM = -1D0/SQRT_2
         DO LR = 1,(2*LMAX+1),2
            ILR = (LR+1)/2
            DO IC = 1,2
               DO IR = 1,IRTOP
                  DC(IR,ILR,-1,IC) = WP*(DB(IR,ILR,-1,IC)+DB(IR,ILR,+1,
     &                               IC))
                  DC(IR,ILR,+1,IC) = WM*(DB(IR,ILR,+1,IC)-DB(IR,ILR,-1,
     &                               IC))
               END DO
            END DO
         END DO
C
         DO LR = 1,(2*LMAX+1),2
            ILR = (LR+1)/2
            DO IR = 1,IRTOP
C
C----------------------------------------------------------- j = c*alpha
C
               JDNST(IR,ILR,+1,IT) = JDNST(IR,ILR,+1,IT)
     &                               + DIMAG(DC(IR,ILR,+1,2))*C_AU
               JDNST(IR,ILR,-1,IT) = JDNST(IR,ILR,-1,IT)
     &                               + DIMAG(DC(IR,ILR,-1,1))*C_AU
            END DO
         END DO
C ======================================================================
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
      LR_MAX = 2*LMAX + 1
      LMR_MAX = (LR_MAX+1)**2
C
      ALLOCATE (JALF_LMCT(NRMAX,LMR_MAX,3,NTMAX))
      JALF_LMCT(:,:,:,:) = 0D0
C
      FILNAM = DATSET0(1:LDATSET0)//'_JALF.dat'
C
      CALL WRHEAD(IOTMP,FILNAM,'J_alpha',NETAB(1))
C
      DO IT = ITBOT,ITTOP
         IM = IMT(IT)
         IF ( FULLPOT ) THEN
            IRTOP = JRCRI(IM)
         ELSE
            IRTOP = JRWS(IM)
         END IF
C
         DO LR = 1,(2*LMAX+1),2
            ILR = (LR+1)/2
            LMR_M1 = LR*LR + LR + (-1) + 1
            LMR_P1 = LR*LR + LR + (+1) + 1
            DO IR = 1,IRTOP
               JALF_LMCT(IR,LMR_M1,1,IT) = JDNST(IR,ILR,-1,IT)
               JALF_LMCT(IR,LMR_P1,2,IT) = JDNST(IR,ILR,+1,IT)
            END DO
         END DO
C
         WRITE (IOTMP,'(''TYPE IT'',/,4I10)') IT,IRTOP,LR_MAX,LMR_MAX
         DO LMR = 1,LMR_MAX
            DO IR = 1,IRTOP
               WRITE (IOTMP,'(3E25.10)') JALF_LMCT(IR,LMR,1:3,IT)
            END DO
         END DO
C
      END DO
C
      CLOSE (IOTMP)
C JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
C
C ======================================================================
C     calculate total magnetic moment in Bohr magneton using relations
C     j^{+}_{l,+1} = -j^{-}_{l,-1} = -i*j^{x}_{l,-1} = i*j^{y}_{l,+1}
C ======================================================================
C
C-------------------- Prefactor dealing with spatial integral over j_phi
      PRE = SQRT((8D0*PI)/3D0)/SQRT(2D0)
C
      LR = 1
      ILR = (LR+1)/2
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = ITBOT,ITTOP
         IM = IMT(IT)
         IRTOP = JRWS(IM)
C
         MAGMOMT(IT) = 0D0
         DO IR = 1,IRTOP
            JPHI = JDNST(IR,ILR,+1,IT) - JDNST(IR,ILR,-1,IT)
            MAGMOMT(IT) = MAGMOMT(IT) + R2DRDI_W_RADINT(IR,IM)
     &                    *JPHI*R(IR,IM)
         END DO
C
         MAGMOMT(IT) = -PRE*MAGMOMT(IT)
         UNITS = (M0_CGS*C_CGS/HBAR_CGS)*A0_CGS/C_AU
         WRITE (6,*) 'UNITS',UNITS
         UNITS = 0.5D0
         WRITE (6,*) 'UNITS',UNITS
         MAGMOMT(IT) = UNITS*MAGMOMT(IT)
C
         WRITE (6,*) 'MAGNETIC MOMENT FROM CURRENT DENSITY FOR IT=',IT,
     &               MAGMOMT(IT)
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO LR = 1,3
         OUTPUT = 'curdns'//CHAR(ICHAR('1')-1+LR)
         OPEN (60,FILE=OUTPUT)
         WRITE (6,*) ' writing current density to  ',OUTPUT
         DO IR = 1,IRTOP
            WRITE (60,'(a,2i5,f10.7,3e15.6,3x,e15.6)') ' A ',LR,IR,
     &             R(IR,1),(JDNST(IR,LR,M,1),M=-1,+1)
         END DO
         CLOSE (60)
      END DO
C
      IF ( .NOT.CALCVECPOT ) RETURN
C ======================================================================
C ======================================================================
C                 calculate vector potential
C ======================================================================
C ======================================================================
C
      DTETADEG = 15D0
      NTETA = NINT(180D0/DTETADEG)
      NRA = 721
      ALLOCATE (APHI1(NRA,NTETA))
      ALLOCATE (APHI2(NRA,NTETA))
C
      NMUEJ = 12
C
      ALLOCATE (XGAU(NMUEJ),WGAU(NMUEJ),RINT(NRMAX))
C
      CALL GAULEG(-1D0,+1D0,XGAU,WGAU,NMUEJ)
C
      NLR = (2*LMAX+1) + 1
      NLMR = NLR*NLR
C
      ALLOCATE (Y_JGRID(NLMR,NMUEJ),STAT=IA_ERR)
C
      DO IMUEJ = 1,NMUEJ
C
         MUEJ = XGAU(IMUEJ)
C
         RHAT(1) = SQRT(1D0-MUEJ*MUEJ)
         RHAT(2) = 0D0
         RHAT(3) = MUEJ
C
         CALL CALC_RHPLM(RHAT(1),RHAT(2),RHAT(3),Y_JGRID(1,IMUEJ),NLR-1,
     &                   NLMR)
      END DO
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = ITBOT,ITTOP
         IM = IMT(IT)
         IRTOP = JRWS(IM)
C
C ======================================================================
C                find expansion coefficients for APHI
C ======================================================================
C
         DRA = RWS(IM)/DBLE(NRA)
C
         NMUEA = NMUEJ
         ALLOCATE (APHIAUX(NMUEA),APHIL(NRMAX,NLABIMAX))
         APHIL(1:NRMAX,1:NLABIMAX) = 0D0
C
         DAPHI = 0D0
C
         DO IRA = 1,IRTOP
C
            RA = R(IRA,IM)
C
            DO IMUEA = 1,NMUEA
C
               MUEA = XGAU(IMUEA)
C
C-----------------------------------------------------------------------
C                integrate over current density distribution
C-----------------------------------------------------------------------
               DO IR = 1,IRTOP
                  RJ = R(IR,IM)
C
                  RSUM = 0D0
                  DO IMUEJ = 1,NMUEJ
C
                     MUEJ = XGAU(IMUEJ)
C
                     JPHI = 0D0
                     DO LR = 1,(2*LMAX+1),2
                        LMJ = LR*LR + LR + (+1) + 1
                        ILR = (LR+1)/2
                        JPHI = JPHI - JDNST(IR,ILR,+1,IT)
     &                         *Y_JGRID(LMJ,IMUEJ)
                     END DO
C
                     CALL APHIRING(RJ,MUEJ,JPHI,RA,MUEA,DAPHI)
C
                     RSUM = RSUM + WGAU(IMUEJ)*DAPHI
C
                  END DO
C
                  RINT(IR) = RSUM*R(IR,IM)*DRDI(IR,IM)
C
               END DO
C
               CALL RRADINT(IM,RINT,APHIAUX(IMUEA))
C-----------------------------------------------------------------------
C
            END DO
C
            DO LR = 1,(2*LMAX+1),2
               LMA = LR*LR + LR + (+1) + 1
               ILR = (LR+1)/2
               RSUM = 0D0
               DO IMUEA = 1,NMUEA
                  RSUM = RSUM + WGAU(IMUEA)*APHIAUX(IMUEA)
     &                   *Y_JGRID(LMA,IMUEA)*PI
               END DO
               APHIL(IRA,ILR) = RSUM/C_AU*SQRT(2.0D0)
            END DO
            WRITE (990,*) R(IRA,IM),APHIL(IRA,1),APHIL(IRA,2),
     &                    APHIL(IRA,3)
C
         END DO
C
C ======================================================================
C              CHECK: compare exact A_phi with dipole approx.
C ======================================================================
C
         DRA = RWS(IM)/DBLE(NRA)
C
         APHI1(1:NRA,1:NTETA) = 0D0
         DAPHI = 0D0
C
         DO ITETA = 1,NTETA
C
            TETA = ITETA*DTETADEG*(PI/180D0)
            TETA = PI/2.0D0
            MUEA = COS(TETA)
C
            DO IRA = 1,NRA
C
               RA = IRA*DRA
C
C-----------------------------------------------------------------------
C                integrate over current density distribution
C-----------------------------------------------------------------------
               DO IR = 1,IRTOP
                  RJ = R(IR,IM)
C
                  RSUM = 0D0
                  DO IMUEJ = 1,NMUEJ
C
                     MUEJ = XGAU(IMUEJ)
C
                     JPHI = 0D0
                     DO LR = 1,(2*LMAX+1),2
                        LMJ = LR*LR + LR + (+1) + 1
                        ILR = (LR+1)/2
                        JPHI = JPHI - JDNST(IR,ILR,+1,IT)
     &                         *Y_JGRID(LMJ,IMUEJ)
                     END DO
C
                     CALL APHIRING(RJ,MUEJ,JPHI,RA,MUEA,DAPHI)
C
                     RSUM = RSUM + WGAU(IMUEJ)*DAPHI
C
                  END DO
C
                  RINT(IR) = RSUM*R(IR,IM)*DRDI(IR,IM)
C
               END DO
C
               CALL RRADINT(IM,RINT,APHI1(IRA,ITETA))
C-----------------------------------------------------------------------
C
               WRITE (998,*) RA,APHI1(IRA,ITETA)/C_AU*SQRT(2.0D0)
               WRITE (999,*) RA,MAGMOMT(IT)*SIN(TETA)/(RA**2)
     &                       /C_AU*SQRT(2.0D0)
            END DO
C
         END DO
C
         WRITE (6,*) 'vecpot current:',APHI1(721,1)/C_AU*SQRT(2.0D0)
         WRITE (6,*) 'vecpot approxi:',MAGMOMT(IT)/(1.0D0*RWS(IM))
     &               **2/C_AU*SQRT(2.0D0)
C
C ======================================================================
C              CHECK: compare exact A_phi with expansion
C ======================================================================
C
         DRA = RWS(IM)/DBLE(NRA)
C
         APHI1(1:NRA,1:NTETA) = 0D0
         APHI2(1:NRA,1:NTETA) = 0D0
         DAPHI = 0D0
C
         DO IMUEA = 1,NMUEA
C
            MUEA = XGAU(IMUEA)
C
            DO IRA = 1,IRTOP
C
               RA = R(IRA,IM)
C
C-----------------------------------------------------------------------
C                integrate over current density distribution
C-----------------------------------------------------------------------
               DO IR = 1,IRTOP
                  RJ = R(IR,IM)
C
                  RSUM = 0D0
                  DO IMUEJ = 1,NMUEJ
C
                     MUEJ = XGAU(IMUEJ)
C
                     JPHI = 0D0
                     DO LR = 1,(2*LMAX+1),2
                        LMJ = LR*LR + LR + (+1) + 1
                        ILR = (LR+1)/2
                        JPHI = JPHI - JDNST(IR,ILR,+1,IT)
     &                         *Y_JGRID(LMJ,IMUEJ)
                     END DO
C
                     CALL APHIRING(RJ,MUEJ,JPHI,RA,MUEA,DAPHI)
C
                     RSUM = RSUM + WGAU(IMUEJ)*DAPHI
C
                  END DO
C
                  RINT(IR) = RSUM*R(IR,IM)*DRDI(IR,IM)
C
               END DO
C
               CALL RRADINT(IM,RINT,APHI1(IRA,IMUEA))
C-----------------------------------------------------------------------
C
               RSUM = 0D0
               DO LR = 1,(2*LMAX+1),2
                  LMJ = LR*LR + LR + (+1) + 1
                  ILR = (LR+1)/2
                  RSUM = RSUM + APHIL(IRA,ILR)*Y_JGRID(LMJ,IMUEA)
               END DO
               APHI2(IRA,IMUEA) = RSUM
C
               WRITE (995,*) RA,APHI1(IRA,IMUEA)/C_AU*SQRT(2.0D0),
     &                       REAL(APHI2(IRA,IMUEA))
            END DO
C
         END DO
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      END
C*==aphiring.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE APHIRING(RJ,MUEJ,JPHI,RA,MUEA,DAPHI)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   *                                                                  *
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C*--APHIRING670
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 DAPHI,JPHI,MUEA,MUEJ,RA,RJ
C
C Local variables
C
      REAL*8 ARS,AUX,EEI,FAC,KEI,KSQ,RHO,RPA,SUM1,SUM2,TETPA,XA,XJ,XPA,
     &       ZA,ZJ,ZPA
      INTEGER KK
C
C*** End of declarations rewritten by SPAG
C
C---position of current loop in cell-centred coordinates
C
      ZJ = RJ*MUEJ
      XJ = SQRT(RJ**2-ZJ**2)
C
C---position of vector potential in cell-centred coordinates
C
      ZA = RA*MUEA
      XA = SQRT(RA**2-ZA**2)
C
C---position of vector potential w.r.t. current loop
C
      ZPA = ZA - ZJ
      XPA = XA
      RPA = SQRT(XPA**2+ZPA**2)
C
      TETPA = ASIN(XPA/RPA)
C
C---calculation of vector potential at RPA,TETPA
C
      KEI = 0.0D0
      EEI = 0.0D0
      RHO = RPA*SIN(TETPA)
      ARS = 4.0*XJ*RHO
C
      AUX = (XJ+RHO)**2 + ZPA**2
      KSQ = ARS/AUX
C
      IF ( RHO.GT.0.0 ) THEN
         SUM1 = 1.0D0
         SUM2 = 1.0D0
         FAC = KSQ*(1.0/2.0)**2
         DO KK = 1,200
            SUM1 = SUM1 + FAC
            SUM2 = SUM2 - FAC/(2.0*KK-1.0)
            FAC = FAC*KSQ*((2.0*KK+1.0)**2/(2.0*KK+2.0)**2)
            IF ( KK.EQ.200 .AND. ABS(FAC).GT.1.D-8 ) WRITE (6,*)
     &            'not converged',ABS(FAC),KSQ,SUM1,SUM2
         END DO
C
         KEI = SUM1*PI/2.0
         EEI = SUM2*PI/2.0
C
         DAPHI = 4*JPHI*XJ*((2D0-KSQ)*KEI-2*EEI)/(KSQ*SQRT(AUX))
C
      ELSE
         DAPHI = 0.0D0
      END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      END
C*==etotavec.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
C=======================================================================
      SUBROUTINE ETOTAVEC(ABITNEW)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   *                                                                  *
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:JRWS,NRMAX,R2DRDI_W_RADINT
      USE MOD_TYPES,ONLY:IMT,NTMAX,ITBOT,ITTOP,NT,JDNST,NAT
      USE MOD_ANGMOM,ONLY:NL,NLABIMAX
      USE MOD_CONSTANTS,ONLY:C_AU
      IMPLICIT NONE
C*--ETOTAVEC763
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ABITNEW(NRMAX,NLABIMAX,-1:+1,NTMAX)
C
C Local variables
C
      REAL*8 EJA,EJAT(NT)
      INTEGER ILA,IM,IR,IRTOP,IT,LA,LMAX,M
C
C*** End of declarations rewritten by SPAG
C
      EJA = 0.0D0
      EJAT(1:NT) = 0.0D0
      LMAX = NL - 1
C
      DO IT = ITBOT,ITTOP
         IM = IMT(IT)
         IRTOP = JRWS(IM)
C
         DO LA = 1,(2*LMAX+1),2
            ILA = (LA+1)/2
            DO M = -1, + 1
C
               DO IR = 1,IRTOP
                  EJAT(IT) = EJAT(IT) + R2DRDI_W_RADINT(IR,IM)
     &                       *SQRT(2.0D0)*JDNST(IR,ILA,-M,IT)
     &                       *ABITNEW(IR,ILA,M,IT)
               END DO
C
            END DO
C
         END DO
C
         EJAT(IT) = 0.5D0*EJAT(IT)/C_AU
         EJA = EJA + EJAT(IT)*NAT(IT)
C
      END DO
      WRITE (*,'(a,F20.10)') 'total EJA  =',EJA
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      END
