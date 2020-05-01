C*==me_nab_nab_sra.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE ME_NAB_NAB_SRA(IFILA,IFILB,IT,MZAZB,MZBZA,MIRR_2,
     &                          MIRR_3,MIRR_4)
C   ********************************************************************
C   *                                                                  *
C   *    read wave function and the calculate matrix elements          *
C   *                                                                  *
C   *             MREG = < Z^+_f | H_lam | Z_i >                       *
C   *             MBAR = < Z^x_i | H_lam | Z_f >                       *
C   *             MIRR = < Z^+_f | H_lam | G^irr_i | H^+_lam' | Z_f >  *
C   *                                                                  *
C   *     with    H_lam = mc alfa_lam = mc alfa * A_lam                *
C   *                                                                  *
C   *                G^irr_L =   Z_L(r) J_L(r') Theta(r'-r)            *
C   *                          + J_L(r) Z_L(r') Theta(r-r')            *
C   *                                                                  *
C   *    according to ELECTRIC DIPOLE selection rules                  *
C   *                                                                  *
C   *    polarisation lam:         ipol= 1,2,3  ==  (+),(-),(z)        *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:IMT
      USE MOD_RMESH,ONLY:NRMAX,JRCRI,FULLPOT,JRWS,JRCUT,NPANMAX,NMMAX,
     &    NPAN,NM,JRNSMIN
      USE MOD_CONSTANTS,ONLY:C0,CI,SQRT_2
      USE MOD_ANGMOM,ONLY:NL,NLM,NKM,NKMMAX
      USE MOD_STR,ONLY:LLMAX
      IMPLICIT NONE
C*--ME_NAB_NAB_SRA30
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IFILA,IFILB,IT
      COMPLEX*16 MIRR_2(NKMMAX,NKMMAX,3,3),MIRR_3(NKMMAX,NKMMAX,3,3),
     &           MIRR_4(NKMMAX,NKMMAX,3,3),MZAZB(NKMMAX,NKMMAX,3),
     &           MZBZA(NKMMAX,NKMMAX,3)
C
C Local variables
C
      LOGICAL CHECK,K_GWFLM(:,:,:),K_WFLM(:,:)
      COMPLEX*16 CRWK1(:),CRWK2(:),CX,CY,FIRR(:),FJAJB(:),FJAZB(:),
     &           FREG(:),FZAJB(:),FZAZB(:),FZBZA(:),GJA(:,:,:,:),
     &           GJB(:,:,:,:),GZA(:,:,:,:),GZB(:,:,:,:),IJAJB(:),
     &           IJAZB(:),IZAJB(:),IZAZB(:),IZBZA(:),JA(:,:,:),JB(:,:,:)
     &           ,MEAUX,ZA(:,:,:),ZB(:,:,:)
      INTEGER IA_ERR,IFLAG1,IFLAG2,IKM,IM,IPOL,IPOLREV,IR,IRTOP,JPOL,
     &        JPOLREV,LAMA1,LAMA2,LAMB3,LAMB4,LM,NLM_GWF,NL_GWF,NPOL
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE FZAZB,FZAJB,FJAZB,FJAJB,FZBZA
      ALLOCATABLE IZBZA,IZAZB,IZAJB,IJAZB,IJAJB
      ALLOCATABLE FREG,FIRR
      ALLOCATABLE CRWK1,CRWK2
      ALLOCATABLE JA,JB,ZA,ZB,K_WFLM,GZA,GZB,GJB,GJA,K_GWFLM
C
      CHECK = .FALSE.
C
      ALLOCATE (CRWK1(NRMAX),CRWK2(NRMAX))
      ALLOCATE (K_WFLM(NLM,NLM))
C
C
      NL_GWF = NL + 1
      NLM_GWF = NL_GWF**2
      ALLOCATE (K_GWFLM(NLM_GWF,3,NKM))
      ALLOCATE (GZA(NRMAX,NLM_GWF,3,NKM))
      ALLOCATE (GZB(NRMAX,NLM_GWF,3,NKM))
      ALLOCATE (GJA(NRMAX,NLM_GWF,3,NKM))
      ALLOCATE (GJB(NRMAX,NLM_GWF,3,NKM))
      ALLOCATE (FZAZB(NRMAX),FZAJB(NRMAX),FJAZB(NRMAX),FJAJB(NRMAX))
      ALLOCATE (IZAZB(NRMAX),IZAJB(NRMAX),IJAZB(NRMAX),IJAJB(NRMAX))
      ALLOCATE (IZBZA(NRMAX),FZBZA(NRMAX),FREG(NRMAX),FIRR(NRMAX))
      ALLOCATE (ZA(NRMAX,NKM,NKM),JA(NRMAX,NKM,NKM))
      ALLOCATE (JB(NRMAX,NKM,NKM),ZB(NRMAX,NKM,NKM),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:meadairr->ZAF'
C
C=======================================================================
C   supply dummy mesh parameters used in full potential calculations
C=======================================================================
C
      IF ( .NOT.FULLPOT ) THEN
         DEALLOCATE (NPAN,JRCUT,JRCRI)
         NPANMAX = 1
         ALLOCATE (NPAN(NMMAX),JRCUT(0:NPANMAX,NMMAX),JRCRI(NMMAX))
         NPAN(1:NM) = 1
         JRCUT(0,1:NM) = 0
         JRCUT(1,1:NM) = JRWS(1:NM)
         JRCRI(1:NM) = JRWS(1:NM)
         LLMAX = 1
         JRNSMIN = 1
      END IF
C
      NPOL = 3
      IM = IMT(IT)
      IRTOP = JRCRI(IM)
      LLMAX = 1
      JRNSMIN = 1
C
      MZAZB = C0
      MZBZA = C0
C=======================================================================
C
C-----------------------------------------------------------------------
C            read wave functions for energies A and B
C-----------------------------------------------------------------------
C
      CALL READ_WFMAT_SRA(IT,IFILA,IRTOP,GZA,ZA,JA,K_WFLM)
C
      CALL READ_WFMAT_SRA(IT,IFILB,IRTOP,GZB,ZB,JB,K_WFLM)
C
C-----------------------------------------------------------------------
C                  calculate gradient of wave functions
C-----------------------------------------------------------------------
      IF ( CHECK ) THEN
         DO IKM = 1,NKM
            JA(1:IRTOP,1:NKM,IKM) = ZA(1:IRTOP,1:NKM,IKM)
            JB(1:IRTOP,1:NKM,IKM) = ZB(1:IRTOP,1:NKM,IKM)
         END DO
      END IF
C
      DO IKM = 1,NKM
C
         CALL GRAD_WF_SRA(ZA(1,1,IKM),K_WFLM(1,IKM),GZA(1,1,1,IKM),
     &                    K_GWFLM(1,1,IKM),NLM_GWF,IM)
C
         CALL GRAD_WF_SRA(JA(1,1,IKM),K_WFLM(1,IKM),GJA(1,1,1,IKM),
     &                    K_GWFLM(1,1,IKM),NLM_GWF,IM)
C
         CALL GRAD_WF_SRA(ZB(1,1,IKM),K_WFLM(1,IKM),GZB(1,1,1,IKM),
     &                    K_GWFLM(1,1,IKM),NLM_GWF,IM)
C
         CALL GRAD_WF_SRA(JB(1,1,IKM),K_WFLM(1,IKM),GJB(1,1,1,IKM),
     &                    K_GWFLM(1,1,IKM),NLM_GWF,IM)
C
      END DO
C-----------------------------------------------------------------------
C                change to circulaer polarisation
C                index 3:  ipol= 1,2,3  ==  (+),(-),(z)
C                 M(X) =   [  M(+) + M(-) ] / SQRT(2)
C                 M(Y) = I*[ -M(+) + M(-) ] / SQRT(2)
C-----------------------------------------------------------------------
C
      CX = 1D0/SQRT_2
      CY = CI/SQRT_2
      DO IKM = 1,NKM
C
         DO LM = 1,NLM
C
            IF ( K_GWFLM(LM,1,IKM) .OR. K_GWFLM(LM,2,IKM) ) THEN
C
               K_GWFLM(LM,1,IKM) = .TRUE.
               K_GWFLM(LM,2,IKM) = .TRUE.
C
               CRWK1(1:IRTOP) = GZA(1:IRTOP,LM,1,IKM)
               CRWK2(1:IRTOP) = GZA(1:IRTOP,LM,2,IKM)
               GZA(1:IRTOP,LM,1,IKM) = CX*CRWK1(1:IRTOP)
     &                                 + CY*CRWK2(1:IRTOP)
               GZA(1:IRTOP,LM,2,IKM) = CX*CRWK1(1:IRTOP)
     &                                 - CY*CRWK2(1:IRTOP)
C
               CRWK1(1:IRTOP) = GJA(1:IRTOP,LM,1,IKM)
               CRWK2(1:IRTOP) = GJA(1:IRTOP,LM,2,IKM)
               GJA(1:IRTOP,LM,1,IKM) = CX*CRWK1(1:IRTOP)
     &                                 + CY*CRWK2(1:IRTOP)
               GJA(1:IRTOP,LM,2,IKM) = CX*CRWK1(1:IRTOP)
     &                                 - CY*CRWK2(1:IRTOP)
C
               CRWK1(1:IRTOP) = GZB(1:IRTOP,LM,1,IKM)
               CRWK2(1:IRTOP) = GZB(1:IRTOP,LM,2,IKM)
               GZB(1:IRTOP,LM,1,IKM) = CX*CRWK1(1:IRTOP)
     &                                 + CY*CRWK2(1:IRTOP)
               GZB(1:IRTOP,LM,2,IKM) = CX*CRWK1(1:IRTOP)
     &                                 - CY*CRWK2(1:IRTOP)
C
               CRWK1(1:IRTOP) = GJB(1:IRTOP,LM,1,IKM)
               CRWK2(1:IRTOP) = GJB(1:IRTOP,LM,2,IKM)
               GJB(1:IRTOP,LM,1,IKM) = CX*CRWK1(1:IRTOP)
     &                                 + CY*CRWK2(1:IRTOP)
               GJB(1:IRTOP,LM,2,IKM) = CX*CRWK1(1:IRTOP)
     &                                 - CY*CRWK2(1:IRTOP)
C
            END IF
C
         END DO
C
      END DO
C
      GZA = 0.5D0*GZA/CI
      GZB = 0.5D0*GZB/CI
      GJA = 0.5D0*GJA/CI
      GJB = 0.5D0*GJB/CI
C
      CRWK1(:) = C0
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C                                                       energy B -- LAM3
      DO LAMB3 = 1,NKM
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                       energy A -- LAM2
         DO LAMA2 = 1,NKM
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP JPOL
            DO JPOL = 1,NPOL
C
               IF ( JPOL.LT.3 ) THEN
                  JPOLREV = 3 - JPOL
               ELSE
                  JPOLREV = JPOL
               END IF
C
               DO IR = 1,IRTOP
                  FZAZB(IR) = C0
                  FZAJB(IR) = C0
                  FJAZB(IR) = C0
                  FJAJB(IR) = C0
C
                  FZBZA(IR) = C0
               END DO
C
               IFLAG1 = 0
C
               CALL MENAB_SRA_INT(FZAZB,ZA,GZA,ZB,GZB,LAMA2,LAMB3,JPOL,
     &                            K_WFLM,K_GWFLM,IFLAG1,IM,CRWK1,
     &                            .FALSE.,NLM_GWF)
C
               CALL MENAB_SRA_INT(FZAJB,ZA,GZA,JB,GJB,LAMA2,LAMB3,JPOL,
     &                            K_WFLM,K_GWFLM,IFLAG1,IM,CRWK1,
     &                            .FALSE.,NLM_GWF)
C
               CALL MENAB_SRA_INT(FJAZB,JA,GJA,ZB,GZB,LAMA2,LAMB3,JPOL,
     &                            K_WFLM,K_GWFLM,IFLAG1,IM,CRWK1,
     &                            .FALSE.,NLM_GWF)
C
               CALL MENAB_SRA_INT(FJAJB,JA,GJA,JB,GJB,LAMA2,LAMB3,JPOL,
     &                            K_WFLM,K_GWFLM,IFLAG1,IM,CRWK1,
     &                            .FALSE.,NLM_GWF)
C
               CALL MENAB_SRA_INT(FZBZA,ZB,GZB,ZA,GZA,LAMB3,LAMA2,
     &                            JPOLREV,K_WFLM,K_GWFLM,IFLAG1,IM,
     &                            CRWK1,.FALSE.,NLM_GWF)
C
C#######################################################################
C
               IF ( IFLAG1.EQ.1 ) THEN
C
                  CALL CRADINT_R(IM,FZAZB,IZAZB)
C
                  CALL CRADINT_INW_R(IM,FZAJB,IZAJB)
C
                  CALL CRADINT_INW_R(IM,FJAZB,IJAZB)
C
                  CALL CRADINT_INW_R(IM,FJAJB,IJAJB)
C
                  CALL CRADINT_R(IM,FZBZA,IZBZA)
C
C***********************************************************************
C***********************************************************************
C     TERM 1    MZBZA(4,1) * TAU(1,2,a) * MZAZB(2,3) * TAU(3,4,b)
C***********************************************************************
C***********************************************************************
C
                  MZAZB(LAMA2,LAMB3,JPOL) = MZAZB(LAMA2,LAMB3,JPOL)
     &               + IZAZB(IRTOP)
C
                  MZBZA(LAMB3,LAMA2,JPOLREV)
     &               = MZBZA(LAMB3,LAMA2,JPOLREV) + IZBZA(IRTOP)
C
C
C***********************************************************************
C***********************************************************************
C     TERM 2    - TAU(3,4,b) * MIRR_2(3,4)                 LAMA1 = LAMA2
C***********************************************************************
C***********************************************************************
C
                  LAMA1 = LAMA2
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C                                                       energy B -- LAM4
                  DO LAMB4 = 1,NKM
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
C
                     DO IPOL = 1,NPOL
C
                        IFLAG2 = 0
C
                        DO IR = 1,IRTOP
                           FREG(IR) = C0
                           FIRR(IR) = C0
                        END DO
C
                        CALL MENAB_SRA_INT(FREG,ZB,GZB,JA,GJA,LAMB4,
     &                     LAMA1,IPOL,K_WFLM,K_GWFLM,IFLAG2,IM,IZAZB,
     &                     .TRUE.,NLM_GWF)
C
                        CALL MENAB_SRA_INT(FIRR,ZB,GZB,ZA,GZA,LAMB4,
     &                     LAMA1,IPOL,K_WFLM,K_GWFLM,IFLAG2,IM,IJAZB,
     &                     .TRUE.,NLM_GWF)
C
                        IF ( IFLAG2.EQ.1 ) THEN
C
                           CRWK1(1:IRTOP) = FREG(1:IRTOP)
     &                        + FIRR(1:IRTOP)
C
                           CALL CRADINT(IM,CRWK1,MEAUX)
C
                           MIRR_2(LAMB4,LAMB3,IPOL,JPOL)
     &                        = MIRR_2(LAMB4,LAMB3,IPOL,JPOL) + MEAUX
C
                        END IF
C
C***********************************************************************
                     END DO
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
                  END DO
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB LAMB4
C
C
C***********************************************************************
C***********************************************************************
C     TERM 3    - TAU(1,2,a) * MIRR_3(1,2)                 LAMB3 = LAMB4
C***********************************************************************
C***********************************************************************
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                       energy A -- LAM1
                  DO LAMA1 = 1,NKM
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
C
                     DO IPOL = 1,NPOL
                        IF ( IPOL.LT.3 ) THEN
                           IPOLREV = 3 - IPOL
                        ELSE
                           IPOLREV = IPOL
                        END IF
C
                        DO IR = 1,IRTOP
                           FREG(IR) = C0
                           FIRR(IR) = C0
                        END DO
C
                        IFLAG2 = 0
C
                        CALL MENAB_SRA_INT(FREG,JB,GJB,ZA,GZA,LAMB3,
     &                     LAMA1,IPOL,K_WFLM,K_GWFLM,IFLAG2,IM,IZAZB,
     &                     .TRUE.,NLM_GWF)
C
                        CALL MENAB_SRA_INT(FIRR,ZB,GZB,ZA,GZA,LAMB3,
     &                     LAMA1,IPOL,K_WFLM,K_GWFLM,IFLAG2,IM,IZAJB,
     &                     .TRUE.,NLM_GWF)
C
                        IF ( IFLAG2.EQ.1 ) THEN
C
                           CRWK1(1:IRTOP) = FREG(1:IRTOP)
     &                        + FIRR(1:IRTOP)
C
                           CALL CRADINT(IM,CRWK1,MEAUX)
C
                           MIRR_3(LAMA1,LAMA2,IPOLREV,JPOLREV)
     &                        = MIRR_3(LAMA1,LAMA2,IPOLREV,JPOLREV)
     &                        + MEAUX
C
                        END IF
C
C***********************************************************************
                     END DO
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
                  END DO
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA LAMA1
C
C
C***********************************************************************
C***********************************************************************
C     TERM 4    MIRR_4(1,3)           LAMA1 = LAMA2  and   LAMB3 = LAMB4
C***********************************************************************
C***********************************************************************
C
                  LAMA1 = LAMA2
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
C
                  DO IPOL = 1,NPOL
                     IF ( IPOL.LT.3 ) THEN
                        IPOLREV = 3 - IPOL
                     ELSE
                        IPOLREV = IPOL
                     END IF
C
                     DO IR = 1,IRTOP
                        FREG(IR) = C0
                        FIRR(IR) = C0
                     END DO
C
                     IFLAG2 = 0
C
                     CALL MENAB_SRA_INT(FREG,JB,GJB,JA,GJA,LAMB3,LAMA1,
     &                                  IPOL,K_WFLM,K_GWFLM,IFLAG2,IM,
     &                                  IZAZB,.TRUE.,NLM_GWF)
C
                     CALL MENAB_SRA_INT(FIRR,ZB,GZB,ZA,GZA,LAMB3,LAMA1,
     &                                  IPOL,K_WFLM,K_GWFLM,IFLAG2,IM,
     &                                  IJAJB,.TRUE.,NLM_GWF)
C
                     IF ( IFLAG2.EQ.1 ) THEN
C
                        CRWK1(1:IRTOP) = FREG(1:IRTOP) + FIRR(1:IRTOP)
C
                        CALL CRADINT(IM,CRWK1,MEAUX)
C
                        MIRR_4(LAMA1,LAMB3,IPOLREV,JPOLREV) = MEAUX
C
                     END IF
C
                  END DO
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
C
C***********************************************************************
C
               END IF
C################################################################ IFLAG1
            END DO
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP JPOL
         END DO
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA LAMA2
      END DO
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
C
      END
C*==read_wfmat_sra.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE READ_WFMAT_SRA(IT,IFIL,IRTOP,WF,WR,WI,K_WFLM)
C   ********************************************************************
C   *                                                                  *
C   *    read regular and irregular wave functions   WR in WI          *
C   *    in compressed form and  expand to the full matrix form        *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:NCPLWFMAX,NT,NLT
      USE MOD_RMESH,ONLY:NRMAX,FULLPOT
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_ANGMOM,ONLY:NLM,NKM,NSPIN,NL
      IMPLICIT NONE
C*--READ_WFMAT_SRA460
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IFIL,IRTOP,IT
      LOGICAL K_WFLM(NLM,NLM)
      COMPLEX*16 WF(NRMAX,NLM),WI(NRMAX,NLM,NKM),WR(NRMAX,NLM,NKM)
C
C Local variables
C
      INTEGER CNLP,IKM,IKMCPLWF(NKM),IKMP,IL,IR,IS,ITP,J,JKM,K,L,
     &        LMOFFSET,LP,NCPLWF,NLML,NSPINP
      CHARACTER*3 STRP
C
C*** End of declarations rewritten by SPAG
C
      WR(:,:,:) = C0
      WI(:,:,:) = C0
      K_WFLM(:,:) = .FALSE.
C
C=======================================================================
C                        FULLPOT
C=======================================================================
      IF ( FULLPOT ) THEN
C
         DO IKM = 1,NKM
C
C-----------------------------------------------------------------------
C                   read regular wave functions
C-----------------------------------------------------------------------
C
            READ (IFIL,REC=IKM+(IT-1)*NKM) ITP,STRP,IKMP,
     &            (IKMCPLWF(J),J=1,NCPLWFMAX),NCPLWF,
     &            ((WF(IR,K),IR=1,IRTOP),K=1,NCPLWF)
C
            IF ( STRP.NE.'REG' .OR. IT.NE.ITP ) THEN
               WRITE (6,*) '##### TROUBLE reading WF:',ITP,STRP,IKMP
               STOP 'in <MENABIRR_SRA>'
            END IF
C
C-----------------------------------------------------------------------
C                   expand to matrix form
C-----------------------------------------------------------------------
C
            DO K = 1,NCPLWF
               JKM = IKMCPLWF(K)
               IF ( JKM.LE.NLM ) K_WFLM(JKM,IKM) = .TRUE.
               WR(1:IRTOP,JKM,IKM) = WF(1:IRTOP,K)
            END DO
C
C-----------------------------------------------------------------------
C                   read irregular wave functions
C-----------------------------------------------------------------------
C
            READ (IFIL,REC=IKM+(IT-1+NT)*NKM) ITP,STRP,IKMP,
     &            ((WF(IR,K),IR=1,IRTOP),K=1,NCPLWF)
C
C-----------------------------------------------------------------------
C                   expand to matrix form
C-----------------------------------------------------------------------
C
            DO K = 1,NCPLWF
               JKM = IKMCPLWF(K)
               WI(1:IRTOP,JKM,IKM) = WF(1:IRTOP,K)
            END DO
C
         END DO
C
C=======================================================================
C                            ASA
C=======================================================================
      ELSE
C
         LMOFFSET = 0
         DO IL = 1,NLT(IT)
            L = IL - 1
            NLML = 2*L + 1
C
C-----------------------------------------------------------------------
C                                                  REGULAR wave function
            READ (IFIL,REC=IL+(IT-1)*NL) ITP,LP,CNLP,NSPINP,STRP,
     &            ((WF(IR,IS),IR=1,IRTOP),IS=1,NSPIN)
C
            IF ( ITP.NE.IT .OR. LP.NE.L .OR. CNLP.NE.NL .OR. 
     &           NSPINP.NE.NSPIN .OR. STRP.NE.'NRR' ) STOP 'READ REG WF'
C
            DO K = 1,NLML
               IKM = K + LMOFFSET
               WR(1:IRTOP,IKM,IKM) = WF(1:IRTOP,1)
            END DO
C
C-----------------------------------------------------------------------
C                                                IRREGULAR wave function
            READ (IFIL,REC=IL+(IT-1+NT)*NL) ITP,LP,CNLP,NSPINP,STRP,
     &            ((WF(IR,IS),IR=1,IRTOP),IS=1,NSPIN)
C
            IF ( ITP.NE.IT .OR. LP.NE.L .OR. CNLP.NE.NL .OR. 
     &           NSPINP.NE.NSPIN .OR. STRP.NE.'NRI' ) STOP 'READ IRR WF'
C
            DO K = 1,NLML
               IKM = K + LMOFFSET
               WI(1:IRTOP,IKM,IKM) = WF(1:IRTOP,1)
            END DO
C
            LMOFFSET = LMOFFSET + NLML
C
         END DO
C
         DO IKM = 1,NKM
            K_WFLM(IKM,IKM) = .TRUE.
         END DO
C
      END IF
C
      END
C*==grad_wf_sra.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE GRAD_WF_SRA(FLM,K_FLM,GFLM,K_GWFLM,NLM_GWF,IM)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the gradient of    SUM_LM  FLM(IR,LM) * Y(LM)         *
C   *                                                                  *
C   *  grad SUM_LM FLM(IR,LM) * Y(LM)                                  *
C   *                                                                  *
C   *     = SUM_LM' GFLM(IR,LM',IC) * Y(LM') * ->E(IC)                 *
C   *                                                                  *
C   *                   with IC = 1,2,3 == x,y,z                       *
C   *                   and REAL spherical harmonics Y(LM)             *
C   *                                                                  *
C   *  making use of gradient formula for COMPLEX spherical harmonics  *
C   *                                                                  *
C   *  K_FLM(LM) indicates non-0 terms  FLM(IR,LM)                     *
C   *                                                                  *
C   *  NOTE: the conversion factors are initialized in first call      *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_ANGMOM,ONLY:L_LM,NL,NLM
      USE MOD_RMESH,ONLY:NRMAX,R,NPAN,JRCUT,DRDI,JRCRI
      USE MOD_CONSTANTS,ONLY:C0,C1,CI,SQRT_2
      IMPLICIT NONE
C*--GRAD_WF_SRA614
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IM,NLM_GWF
      COMPLEX*16 FLM(NRMAX,NLM),GFLM(NRMAX,NLM_GWF,3)
      LOGICAL K_FLM(NLM),K_GWFLM(NLM_GWF,3)
C
C Local variables
C
      REAL*8 CGC_RACAH
      COMPLEX*16 DFDR(:),F,FNAB,FNRC,FOVR,FRMIN(:),FRPLS(:),RC_GWF(:,:),
     &           SMIN,SNAB(-1:+1,3),SPLS,TGMIN(:,:),TGPLS(:,:)
      REAL*8 FAC_LMIN1,FAC_LPLS1,SGMIN(:,:,:),SGPLS(:,:,:)
      INTEGER IC,IPAN,IR,IR1,IRTOP,L,L2,L3,LM,LM1,LM2,LM3,M,M1,M1LIM,
     &        M1STEP,M2,M3,M3LIM,M3STEP,MM,MMLIM,NL_GWF,NR
      LOGICAL INITIALIZE
      SAVE NL_GWF,SGMIN,SGPLS
      EXTERNAL ZAXPY
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE RC_GWF,SGPLS,SGMIN,TGPLS,TGMIN
      ALLOCATABLE FRPLS,FRMIN,DFDR
C
      DATA INITIALIZE/.TRUE./
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
      IF ( INITIALIZE ) THEN
C
C ----------------------------------------------------------------------
C     allow to deal with full wave functions up to NL
C ----------------------------------------------------------------------
C
         NL_GWF = NL + 1
         IF ( NLM_GWF.NE.NL_GWF**2 ) STOP '<GRAD_WF_SRA>  --> NLM_GWF'
C
C-----------------------------------------------------------------------
C   transformation matrix  RC_GWF = RC:  complex -> real spher. harm.
C-----------------------------------------------------------------------
C
         ALLOCATE (RC_GWF(NLM_GWF,NLM_GWF))
C
         CALL CALC_RC(NL_GWF-1,RC_GWF,NLM_GWF,1)
C
C-----------------------------------------------------------------------
C            transformation matrices  SGPLS and SGMIN
C-----------------------------------------------------------------------
C
         SNAB(:,:) = C0
         SNAB(-1,1) = +C1/SQRT_2
         SNAB(+1,1) = -C1/SQRT_2
         SNAB(-1,2) = +CI/SQRT_2
         SNAB(+1,2) = +CI/SQRT_2
         SNAB(0,3) = C1
C
         ALLOCATE (SGPLS(NLM_GWF,NLM_GWF,3),SGMIN(NLM_GWF,NLM_GWF,3))
         ALLOCATE (TGPLS(NLM_GWF,NLM_GWF),TGMIN(NLM_GWF,NLM_GWF))
         SGPLS(:,:,:) = 0D0
         SGMIN(:,:,:) = 0D0
C
         LM = 0
         DO L = 0,NL_GWF - 1
            FAC_LPLS1 = SQRT(DBLE(L+1)/DBLE(2*L+3))
            FAC_LMIN1 = SQRT(DBLE(L)/DBLE(2*L-1))
C
            DO M = -L, + L
               LM = LM + 1
C
               DO IC = 1,3
C
                  TGPLS(:,:) = 0D0
                  TGMIN(:,:) = 0D0
C
                  IF ( IC.LE.2 ) THEN
                     MMLIM = 1
                  ELSE
                     MMLIM = 0
                  END IF
C
                  DO MM = -MMLIM, + MMLIM,2
C
                     FNAB = SNAB(MM,IC)
C
                     M1LIM = ABS(M)
                     M1STEP = MAX(1,2*M1LIM)
                     DO M1 = -M1LIM, + M1LIM,M1STEP
                        LM1 = L*(L+1) + M1 + 1
                        FNRC = FNAB*DCONJG(RC_GWF(LM,LM1))
C
C-------------------------------------------------------------- L2 = L+1
                        L2 = L + 1
                        M2 = M1 + MM
                        LM2 = L2*(L2+1) + M2 + 1
C
                        IF ( L2.LE.NL_GWF-1 .AND. ABS(M2).LE.L2 ) THEN
C
                           F = FNRC*FAC_LPLS1*CGC_RACAH(DBLE(L),1D0,
     &                         DBLE(L2),DBLE(M1),DBLE(MM),DBLE(M2))
C
                           L3 = L2
                           M3LIM = ABS(M2)
                           M3STEP = MAX(1,2*M3LIM)
                           DO M3 = -M3LIM, + M3LIM,M3STEP
                              LM3 = L3*(L3+1) + M3 + 1
                              TGPLS(LM,LM3) = TGPLS(LM,LM3)
     &                           + F*RC_GWF(LM3,LM2)
                           END DO
                        END IF
C
C-------------------------------------------------------------- L2 = L-1
                        L2 = L - 1
                        M2 = M1 + MM
                        LM2 = L2*(L2+1) + M2 + 1
C
                        IF ( L2.GE.0 .AND. ABS(M2).LE.L2 ) THEN
C
                           F = FNRC*FAC_LMIN1*CGC_RACAH(DBLE(L),1D0,
     &                         DBLE(L2),DBLE(M1),DBLE(MM),DBLE(M2))
C
                           L3 = L2
                           M3LIM = ABS(M2)
                           M3STEP = MAX(1,2*M3LIM)
                           DO M3 = -M3LIM, + M3LIM,M3STEP
                              LM3 = L3*(L3+1) + M3 + 1
                              TGMIN(LM,LM3) = TGMIN(LM,LM3)
     &                           + F*RC_GWF(LM3,LM2)
                           END DO
                        END IF
C-----------------------------------------------------------------------
C
                     END DO
C
                  END DO
C
C-----------------------------------------------------------------------
                  DO LM3 = 1,NLM_GWF
C
                     SGMIN(LM,LM3,IC) = DREAL(TGMIN(LM,LM3))
                     SGPLS(LM,LM3,IC) = DREAL(TGPLS(LM,LM3))
C
                     IF ( ABS(DIMAG(TGMIN(LM,LM3))).GT.1D-10 )
     &                    WRITE (6,*) 'TGMIN ',LM,LM3,TGMIN(LM,LM3)
                     IF ( ABS(DIMAG(TGPLS(LM,LM3))).GT.1D-10 )
     &                    WRITE (6,*) 'TGPLS ',LM,LM3,TGPLS(LM,LM3)
                  END DO
C-----------------------------------------------------------------------
C
               END DO
C
            END DO
         END DO
C
         DEALLOCATE (RC_GWF)
C
         IF ( IPRINT.GE.1 ) THEN
            M = NLM_GWF
            CALL RMATSTRUCT('SGMIN1',SGMIN(1,1,1),M,M,1,1,0,1D-8,6)
            CALL RMATSTRUCT('SGMIN2',SGMIN(1,1,2),M,M,1,1,0,1D-8,6)
            CALL RMATSTRUCT('SGMIN3',SGMIN(1,1,3),M,M,1,1,0,1D-8,6)
C
            CALL RMATSTRUCT('SGPLS1',SGPLS(1,1,1),M,M,1,1,0,1D-8,6)
            CALL RMATSTRUCT('SGPLS2',SGPLS(1,1,2),M,M,1,1,0,1D-8,6)
            CALL RMATSTRUCT('SGPLS3',SGPLS(1,1,3),M,M,1,1,0,1D-8,6)
         END IF
C
         INITIALIZE = .FALSE.
C
      END IF
C=======================================================================
C
      ALLOCATE (FRPLS(NRMAX),FRMIN(NRMAX),DFDR(NRMAX))
C
      IRTOP = JRCRI(IM)
      GFLM(:,:,:) = 0D0
C
      DO LM = 1,NLM
C
         IF ( .NOT.K_FLM(LM) ) CYCLE
C
         DO IPAN = 1,NPAN(IM)
            IR1 = JRCUT(IPAN-1,IM) + 1
            NR = JRCUT(IPAN,IM) - JRCUT(IPAN-1,IM)
C
            CALL CALC_DCFDR(FLM(IR1,LM),DFDR(IR1),DRDI(IR1,IM),NR)
C
         END DO
C
         L = L_LM(LM)
C
         DO IR = 1,IRTOP
            FOVR = FLM(IR,LM)/R(IR,IM)
            FRMIN(IR) = -(DFDR(IR)+(L+1)*FOVR)
            FRPLS(IR) = DFDR(IR) - L*FOVR
         END DO
C
C-----------------------------------------------------------------------
         DO IC = 1,3
            DO LM3 = 1,NLM_GWF
C
               SMIN = DCMPLX(SGMIN(LM,LM3,IC))
               IF ( ABS(SMIN).GT.1D-8 ) THEN
                  K_GWFLM(LM3,IC) = .TRUE.
                  CALL ZAXPY(IRTOP,SMIN,FRMIN,1,GFLM(1,LM3,IC),1)
               END IF
C
               SPLS = DCMPLX(SGPLS(LM,LM3,IC))
               IF ( ABS(SPLS).GT.1D-8 ) THEN
                  K_GWFLM(LM3,IC) = .TRUE.
                  CALL ZAXPY(IRTOP,SPLS,FRPLS,1,GFLM(1,LM3,IC),1)
               END IF
C
            END DO
         END DO
C-----------------------------------------------------------------------
C
      END DO
C
      END
C*==menab_sra_int.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MENAB_SRA_INT(FZAZB,ZA,GZA,ZB,GZB,LAMA,LAMB,IPOL,
     &                         K_WFLM,K_GWFLM,IFLAG,IM,WFUN,K_WFUN,
     &                         NLM_GWF)
C   ********************************************************************
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:NRMAX,JRCUT,NPAN,NM,NSF,LMISF,FLMSF,JRCRI,
     &    R2DRDI,FULLPOT
      USE MOD_ANGMOM,ONLY:L_LM,M_LM,NLM,NL,NKM
      IMPLICIT NONE
C*--MENAB_SRA_INT860
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IFLAG,IM,IPOL,LAMA,LAMB,NLM_GWF
      LOGICAL K_WFUN
      COMPLEX*16 FZAZB(NRMAX),GZA(NRMAX,NLM_GWF,3,NKM),
     &           GZB(NRMAX,NLM_GWF,3,NKM),WFUN(NRMAX),ZA(NRMAX,NLM,NKM),
     &           ZB(NRMAX,NLM,NKM)
      LOGICAL K_GWFLM(NLM_GWF,3,NKM),K_WFLM(NLM,NLM)
C
C Local variables
C
      REAL*8 G123,ZZGNT(:,:,:,:)
      REAL*8 GAUNT_RYLM
      INTEGER IA_ERR,ICALL,ILMA,ILMB,IM_INIT,IR,IROFF,IRSF,ISF,L1,L2,LM,
     &        LM1,LM2,LMAX,LSF,M1,M2,MSF,NSFZZMAX
      COMPLEX*16 W,WR(:),ZAZB
      SAVE NSFZZMAX,ZZGNT
C
C*** End of declarations rewritten by SPAG
C
      DATA ICALL/0/
C
      ALLOCATABLE ZZGNT,WR
C
C=======================================================================
C                           INITIALIZE
C=======================================================================
      IF ( ICALL.EQ.0 .AND. FULLPOT ) THEN
C
         ICALL = 1
C
C-----------------------------------------------------------------------
C             Gaunt coefficients for REAL spherical harmonics
C-----------------------------------------------------------------------
C
         LMAX = NL - 1
C
         NSFZZMAX = 0
         DO IM_INIT = 1,NM
            DO ISF = 1,NSF(IM_INIT)
               LM = LMISF(ISF,IM_INIT)
               LSF = L_LM(LM)
               IF ( LSF.LE.2*LMAX ) NSFZZMAX = MAX(NSFZZMAX,ISF)
            END DO
         END DO
C
         ALLOCATE (ZZGNT(NLM,NLM,NSFZZMAX,NM),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP 'alloc:FPNRSSITE-> ZZGNT'
C
         DO IM_INIT = 1,NM
            DO ISF = 1,MIN(NSFZZMAX,NSF(IM_INIT))
               LM = LMISF(ISF,IM_INIT)
               LSF = L_LM(LM)
               MSF = M_LM(LM)
               DO LM2 = 1,NLM
                  L2 = L_LM(LM2)
                  M2 = M_LM(LM2)
                  DO LM1 = 1,NLM
                     L1 = L_LM(LM1)
                     M1 = M_LM(LM1)
                     ZZGNT(LM1,LM2,ISF,IM_INIT)
     &                  = GAUNT_RYLM(L1,M1,L2,M2,LSF,MSF)
                  END DO
               END DO
            END DO
         END DO
C
      END IF
C=======================================================================
C
C
C=======================================================================
C                weight for radial integrals
C=======================================================================
C
      ALLOCATE (WR(NRMAX))
C
      IF ( K_WFUN ) THEN
         DO IR = 1,JRCRI(IM)
            WR(IR) = WFUN(IR)*R2DRDI(IR,IM)
         END DO
      ELSE
         DO IR = 1,JRCRI(IM)
            WR(IR) = R2DRDI(IR,IM)
         END DO
      END IF
C
C=======================================================================
C                   set up radial integrands
C=======================================================================
C
      DO ILMA = 1,NLM
C
         DO ILMB = 1,NLM
C
C-----------------------------------------------------------------------
C     integrals within first panel
C     ASA: full integration range      FP:  muffin-tin regime
C-----------------------------------------------------------------------
C
            IF ( ILMA.EQ.ILMB .AND. K_WFLM(ILMA,LAMA) .AND. 
     &           K_GWFLM(ILMB,IPOL,LAMB) ) THEN
               DO IR = 1,JRCUT(1,IM)
                  ZAZB = ZA(IR,ILMA,LAMA)*GZB(IR,ILMB,IPOL,LAMB)
                  FZAZB(IR) = FZAZB(IR) + WR(IR)*ZAZB
               END DO
               IFLAG = 1
            END IF
C
            IF ( ILMA.EQ.ILMB .AND. K_WFLM(ILMB,LAMB) .AND. 
     &           K_GWFLM(ILMA,IPOL,LAMA) ) THEN
               DO IR = 1,JRCUT(1,IM)
                  ZAZB = GZA(IR,ILMA,IPOL,LAMA)*ZB(IR,ILMB,LAMB)
                  FZAZB(IR) = FZAZB(IR) - WR(IR)*ZAZB
               END DO
               IFLAG = 1
            END IF
C
C-----------------------------------------------------------------------
C     integrals within interstitial regime  --  only FP
C-----------------------------------------------------------------------
C
            IROFF = JRCUT(1,IM)
C
            DO ISF = 1,MIN(NSFZZMAX,NSF(IM))
               G123 = ZZGNT(ILMB,ILMA,ISF,IM)
               IF ( ABS(G123).GT.1D-8 ) THEN
C
                  IF ( K_WFLM(ILMA,LAMA) .AND. K_GWFLM(ILMB,IPOL,LAMB) )
     &                 THEN
                     DO IR = JRCUT(1,IM) + 1,JRCUT(NPAN(IM),IM)
                        IRSF = IR - IROFF
                        W = WR(IR)*FLMSF(IRSF,ISF,IM)*G123
                        ZAZB = ZA(IR,ILMA,LAMA)*GZB(IR,ILMB,IPOL,LAMB)
                        FZAZB(IR) = FZAZB(IR) + W*ZAZB
                     END DO
                     IFLAG = 1
                  END IF
C
                  IF ( K_WFLM(ILMB,LAMB) .AND. K_GWFLM(ILMA,IPOL,LAMA) )
     &                 THEN
                     DO IR = JRCUT(1,IM) + 1,JRCUT(NPAN(IM),IM)
                        IRSF = IR - IROFF
                        W = WR(IR)*FLMSF(IRSF,ISF,IM)*G123
                        ZAZB = GZA(IR,ILMA,IPOL,LAMA)*ZB(IR,ILMB,LAMB)
                        FZAZB(IR) = FZAZB(IR) - W*ZAZB
                     END DO
                     IFLAG = 1
                  END IF
C
               END IF
            END DO
C
         END DO
      END DO
C
      END
