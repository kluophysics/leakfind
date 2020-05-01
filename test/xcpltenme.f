C*==xcpltenme.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      SUBROUTINE XCPLTENME(DELMSST)
C   ********************************************************************
C   *                                                                  *
C   *  matrix elements for the exchange tensor  J_ij                   *
C   *                                                                  *
C   *        < Z(LAM) | sigma_alpha  B | Z(LAM') >   alpha=x,y,z       *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_ANGMOM,ONLY:NXM,NKMMAX,IMKM_IKM,NCPLWF_LB,NCPLWF_RA,NKM,
     &    AME_G,A_SIG_RLM,AG_RGNT,ISMT,NLM_AME_RLM_EXT,NKM_EXT,NPOL
      USE MOD_TYPES,ONLY:NT,IMT,BT,BNST,NKM_T,IKMCPLWF_LB,IKMCPLWF_RA,
     &    NCPLWFMAX,ZGLB,ZFLB,ZGRA,ZFRA,NLMFP,NLMFPT,KLMFP,ITBOT,ITTOP,Z
      USE MOD_RMESH,ONLY:NRMAX,JRWS,FULLPOT,JRCRI,JRNS1,R2DRDI_W_RADINT
      USE MOD_FILES,ONLY:IPRINT,IFILCBWF
      USE MOD_CONSTANTS,ONLY:CI,C0,SQRT_2,Y00
      IMPLICIT NONE
C*--XCPLTENME20
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 DELMSST(NXM,NXM,3,NT)
C
C Local variables
C
      REAL*8 BWLM(:,:)
      INTEGER IA,IB,IFIL_LHSB,IFIL_RHSA,IFIL_RHSB,IKMA,IKMB,IM,IMKMA,
     &        IMKMB,IPOL,IR,IRBOT,IRTOP,IT,J1,J2,LAMA,LAMB,LM,N
      COMPLEX*16 MESPH(3),ZFZF(:,:),ZGZG(:,:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE BWLM,ZFZF,ZGZG
C
C    AME polarisation lamda:      ipol = 1, 2, 3  ==  (-), (0), (+)
C    ME  polarisation lamda:      ipol = 1, 2, 3  ==  (x), (y), (z)
C
      IF ( IREL.LE.2 ) RETURN
C
      ALLOCATE (ZGRA(NRMAX,NCPLWFMAX,NKM),ZFRA(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (ZFLB(NRMAX,NCPLWFMAX,NKM),ZGLB(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (BWLM(NRMAX,NLMFP))
      ALLOCATE (ZFZF(NCPLWFMAX,NCPLWFMAX),ZGZG(NCPLWFMAX,NCPLWFMAX))
C
C=======================================================================
C
C         set the chanel number for the LHS wave functions
C
      IFIL_RHSA = IFILCBWF
      IFIL_RHSB = IFILCBWF
C
      CALL SET_IFIL_LHS(IFIL_RHSB,IFIL_LHSB)
C
C=======================================================================
      IF ( .NOT.ALLOCATED(A_SIG_RLM) ) THEN
C
         ALLOCATE (A_SIG_RLM(NKM_EXT,NKM_EXT,3,NLM_AME_RLM_EXT))
C
         N = NKM_EXT
C
         DO IPOL = 1,3
C
            DO LM = 1,NLM_AME_RLM_EXT
               A_SIG_RLM(1:N,1:N,IPOL,LM)
     &            = MATMUL(AME_G(1:N,1:N,IPOL,ISMT),AG_RGNT(1:N,1:N,LM))
            END DO
C
         END DO
C
      END IF
C
      DELMSST(:,:,:,:) = C0
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = ITBOT,ITTOP
C
         IF ( Z(IT).EQ.0 ) CYCLE
C
         N = NKM_T(IT)
         IM = IMT(IT)
         IF ( FULLPOT ) THEN
            IRTOP = JRCRI(IM)
         ELSE
            IRTOP = JRWS(IM)
            NLMFPT(IT) = 1
            KLMFP(1,IT) = 1
         END IF
C
C------------------------- multiply B with radial and integration weight
C--- in case of FULLPOT: BNST is already convoluted with shape functions
C
         BWLM(:,:) = 0D0
C
         DO IR = 1,IRTOP
            BWLM(IR,1) = R2DRDI_W_RADINT(IR,IM)*BT(IR,IT)/Y00
         END DO
C
         DO LM = 2,NLMFPT(IT)
            IF ( KLMFP(LM,IT).NE.0 ) THEN
               DO IR = JRNS1(IM),IRTOP
                  BWLM(IR,LM) = R2DRDI_W_RADINT(IR,IM)*BNST(IR,LM,IT)
               END DO
            END IF
         END DO
C
C ----------------------------------------- read in wavefunctions for LB
C
         CALL WAVFUN_READ_REL(IFIL_LHSB,IT,0,ZGLB,ZFLB,ZGLB,ZFLB,IRTOP,
     &                        NCPLWF_LB,IKMCPLWF_LB)
C
C ----------------------------------------- read in wavefunctions for RA
C
         CALL WAVFUN_READ_REL(IFIL_RHSA,IT,0,ZGRA,ZFRA,ZGRA,ZFRA,IRTOP,
     &                        NCPLWF_RA,IKMCPLWF_RA)
C
C=======================================================================
C                      calculate matrix elements
C=======================================================================
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
         LOOP_LAMA:DO LAMA = 1,N
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            LOOP_LAMB:DO LAMB = 1,N
C
               MESPH(:) = 0D0
C
               LOOP_LM:DO LM = 1,NLMFPT(IT)
C
                  DO IA = 1,NCPLWF_RA(LAMA)
                     J2 = IKMCPLWF_RA(IA,LAMA)
                     DO IB = 1,NCPLWF_LB(LAMB)
                        J1 = IKMCPLWF_LB(IB,LAMB)
                        DO IPOL = 1,NPOL
                           IF ( ABS(A_SIG_RLM(J1,J2,IPOL,LM)).GT.1D-8 )
     &                          GOTO 5
                        END DO
                     END DO
                  END DO
C
                  CYCLE LOOP_LM
C ----------------------------------------------------------------------
C
C ---------------------------------- non-0 angular matrix elements found
C ------------------------------------- calculate radial matrix elements
C
 5                CONTINUE
                  IF ( LM.EQ.1 ) THEN
                     IRBOT = 1
                  ELSE
                     IRBOT = JRNS1(IM)
                  END IF
C
                  ZGZG(:,:) = 0D0
                  ZFZF(:,:) = 0D0
C
                  DO IA = 1,NCPLWF_RA(LAMA)
                     DO IB = 1,NCPLWF_LB(LAMB)
                        DO IR = IRBOT,IRTOP
                           ZGZG(IB,IA) = ZGZG(IB,IA) + ZGLB(IR,IB,LAMB)
     &                        *ZGRA(IR,IA,LAMA)*BWLM(IR,LM)
                           ZFZF(IB,IA) = ZFZF(IB,IA) + ZFLB(IR,IB,LAMB)
     &                        *ZFRA(IR,IA,LAMA)*BWLM(IR,LM)
                        END DO
                     END DO
                  END DO
C
C -------------------------------------- calculate total matrix elements
C
                  DO IA = 1,NCPLWF_RA(LAMA)
                     IKMA = IKMCPLWF_RA(IA,LAMA)
                     IMKMA = IMKM_IKM(IKMA)
C
                     DO IB = 1,NCPLWF_LB(LAMB)
                        IKMB = IKMCPLWF_LB(IB,LAMB)
                        IMKMB = IMKM_IKM(IKMB)
C
                        DO IPOL = 1,NPOL
                           MESPH(IPOL) = MESPH(IPOL)
     &                        + A_SIG_RLM(IKMB,IKMA,IPOL,LM)*ZGZG(IB,IA)
                        END DO
C
                        IF ( (IMKMB.LE.N) .AND. (IMKMA.LE.N) ) THEN
                           DO IPOL = 1,NPOL
                              MESPH(IPOL) = MESPH(IPOL)
     &                           - A_SIG_RLM(IMKMB,IMKMA,IPOL,LM)
     &                           *ZFZF(IB,IA)
                           END DO
C
                        END IF
C
                     END DO
                  END DO
C
               END DO LOOP_LM
C
C------------------------ store and convert spin polarisation to x, y, z
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C                            TO BE CHECKED
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
               DELMSST(LAMB,LAMA,1,IT) = (-MESPH(3)+MESPH(1))/SQRT_2
C
               DELMSST(LAMB,LAMA,2,IT) = (MESPH(3)+MESPH(1))*CI/SQRT_2
C
               DELMSST(LAMB,LAMA,3,IT) = MESPH(2)
C
            END DO LOOP_LAMB
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
         END DO LOOP_LAMA
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
C------------------------------------------------------------------------
C
         IF ( IPRINT.GE.3 ) THEN
            CALL CMATSTRUCT('DELM X  ',DELMSST(1,1,1,IT),NKMMAX,NKMMAX,
     &                      3,3,0,1D-8,6)
            CALL CMATSTRUCT('DELM Y  ',DELMSST(1,1,2,IT),NKMMAX,NKMMAX,
     &                      3,3,0,1D-8,6)
            CALL CMATSTRUCT('DELM Z  ',DELMSST(1,1,3,IT),NKMMAX,NKMMAX,
     &                      3,3,0,1D-8,6)
         END IF
C
C-----------------------------------------------------------------------
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
C=======================================================================
      DEALLOCATE (ZGRA,ZFRA,ZFLB,ZGLB)
C
      END
