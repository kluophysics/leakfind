C*==linresp_vertex.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_VERTEX(TAUQA,TAUQB,TSSTA,TSSTB,UMAT_VTA,
     &                          UMAT_VTB,MSSQA,MSSQB)
C   ********************************************************************
C   *                                                                  *
C   *   calculate the vertex-corrected CHI matrix                      *
C   *                                                                  *
C   *           chi(vertex) = [ 1 - chi * w ]^-1 * chi                 *
C   *                                                                  *
C   *   NOTE: the matrix  CHIZ  is overwritten                         *
C   *                                                                  *
C   ********************************************************************
      USE MOD_LINRESP,ONLY:CHIZ,NZ12
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,WKM1,IPIVKM
      USE MOD_TYPES,ONLY:NTMAX
      USE MOD_SITES,ONLY:NQMAX,NQ_CHI,NOQ,ITOQ,IQBOT_CHI,IQTOP_CHI
      USE MOD_THERMAL,ONLY:NVIBRA,NFLUCT,FMAT_FT,X_VFT,NVTMAX
      IMPLICIT NONE
C*--LINRESP_VERTEX19
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='LINRESP_VERTEX')
C
C Dummy arguments
C
      COMPLEX*16 MSSQA(NKMMAX,NKMMAX,NQMAX),MSSQB(NKMMAX,NKMMAX,NQMAX),
     &           TAUQA(NKMMAX,NKMMAX,NQMAX),TAUQB(NKMMAX,NKMMAX,NQMAX),
     &           TSSTA(NKMMAX,NKMMAX,NTMAX),TSSTB(NKMMAX,NKMMAX,NTMAX),
     &           UMAT_VTA(NKMMAX,NKMMAX,NVTMAX),
     &           UMAT_VTB(NKMMAX,NKMMAX,NVTMAX)
C
C Local variables
C
      COMPLEX*16 AMAT(:,:),BMAT(:,:),CXAXB,DMAT_VFTA(:,:),DMAT_VFTB(:,:)
     &           ,DMSSA(:,:),DMSSB(:,:),MSS_VFTA(:,:),MSS_VFTB(:,:),
     &           TSS_FTA(:,:),TSS_FTB(:,:),TSS_VFTA(:,:),TSS_VFTB(:,:),
     &           XMAT_VFTA(:,:),XMAT_VFTB(:,:,:)
      INTEGER I,I1,I2,IA_ERR,IFLUCT,IFT,INFO,IO,IPIV(:),IQ,IT,IVFT,
     &        IVIBRA,IVT,J,K1,K2,L1,L2,L3,L4,LWORK,M,N,NDIMCHI,NKMSQ,Z2
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE AMAT,BMAT,DMSSA,DMSSB,XMAT_VFTA,XMAT_VFTB,IPIV
      ALLOCATABLE TSS_FTA,TSS_VFTA,MSS_VFTA,DMAT_VFTA
      ALLOCATABLE TSS_FTB,TSS_VFTB,MSS_VFTB,DMAT_VFTB
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (TSS_VFTA(NKMMAX,NKMMAX),TSS_FTA(NKMMAX,NKMMAX))
      ALLOCATE (TSS_VFTB(NKMMAX,NKMMAX),TSS_FTB(NKMMAX,NKMMAX))
      ALLOCATE (MSS_VFTA(NKMMAX,NKMMAX),DMAT_VFTA(NKMMAX,NKMMAX))
      ALLOCATE (MSS_VFTB(NKMMAX,NKMMAX),DMAT_VFTB(NKMMAX,NKMMAX))
      ALLOCATE (DMSSA(NKMMAX,NKMMAX),DMSSB(NKMMAX,NKMMAX))
      ALLOCATE (XMAT_VFTA(NKMMAX,NKMMAX))
      ALLOCATE (XMAT_VFTB(NKMMAX,NKMMAX,2))
C
      NDIMCHI = NKMMAX*NKMMAX*NQ_CHI
      ALLOCATE (IPIV(NDIMCHI))
      ALLOCATE (AMAT(NDIMCHI,NDIMCHI),BMAT(NDIMCHI,NDIMCHI),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate AMAT')
C
      NKMSQ = NKM*NKM
C
C-----------------------------------------------------------------------
C NOTE: indexing of CHI in <SIGKLOOPS>, <LINRESP_VERTEX> and <SIG1>
C       as well as of  w  in <LINRESP_VERTEX>
C       have to be consistent with loop sequence (L1,L4,L2,L3)
C-----------------------------------------------------------------------
C
C      IF ( CHECK ) THEN
C        CALL  CMATSTR('MSS_VFTA 1',MSS_VFTA(1,1,1),NN,MM,3,3,0,1D-8,6)
C        CALL  CMATSTR('MSS_VFTA 2',MSS_VFTA(1,1,2),NN,MM,3,3,0,1D-8,6)
C        CALL  CMATSTR('MSS_VFTB 1',MSS_VFTB(1,1,1),NN,MM,3,3,0,1D-8,6)
C        CALL  CMATSTR('MSS_VFTB 2',MSS_VFTB(1,1,2),NN,MM,3,3,0,1D-8,6)
C        CALL  CMATSTR('MSSQA 1',MSSQA(1,1,1),NN,MM,3,3,0,1D-8,6)
C        CALL  CMATSTR('MSSQB 1',MSSQB(1,1,1),NN,MM,3,3,0,1D-8,6)
C
C CALL  CMATSTR('DMAT_VFTA 1  ',DMAT_VFTA(1,1,1),NN,MM,3,3,0,1D-8,6)
C CALL  CMATSTR('XMATTA 1     ',XMAT_VFTA(1,1,1),NN,MM,3,3,0,1D-8,6)
C CALL  CMATSTR('DMAT_VFTB 1  ',DMAT_VFTB(1,1,1),NN,MM,3,3,0,1D-8,6)
C CALL  CMATSTR('XMATTB 1     ',XMAT_VFTB(1,1,1,1),NN,MM,3,3,0,1D-8,6)
C CALL  CMATSTR('DMAT_VFTA 2  ',DMAT_VFTA(1,1,2),NN,MM,3,3,0,1D-8,6)
C CALL  CMATSTR('XMATTA 2     ',XMAT_VFTA(1,1,2),NN,MM,3,3,0,1D-8,6)
C CALL  CMATSTR('DMAT_VFTB 2  ',DMAT_VFTB(1,1,2),NN,MM,3,3,0,1D-8,6)
C CALL  CMATSTR('XMATTB 2     ',XMAT_VFTB(1,1,2,1),NN,MM,3,3,0,1D-8,6)
C CALL  CMATSTR('DMAT_VFTA 1 *',DMAT_VFTA(1,1,1),NN,MM,3,3,0,1D-8,6)
C CALL  CMATSTR('XMATTA 1 *   ',XMAT_VFTA(1,1,1),NN,MM,3,3,0,1D-8,6)
C CALL  CMATSTR('DMAT_VFTB 1 *',DMAT_VFTB(1,1,1),NN,MM,3,3,0,1D-8,6)
C CALL  CMATSTR('XMATTB 1 *   ',XMAT_VFTB(1,1,1,2),NN,MM,3,3,0,1D-8,6)
C CALL  CMATSTR('DMAT_VFTA 2 *',DMAT_VFTA(1,1,2),NN,MM,3,3,0,1D-8,6)
C CALL  CMATSTR('XMATTA 2 *   ',XMAT_VFTA(1,1,2),NN,MM,3,3,0,1D-8,6)
C CALL  CMATSTR('DMAT_VFTB 2 *',DMAT_VFTB(1,1,2),NN,MM,3,3,0,1D-8,6)
C CALL  CMATSTR('XMATTB 2 *   ',XMAT_VFTB(1,1,2,2),NN,MM,3,3,0,1D-8,6)
C      END IF
C
      N = NKM
      M = NKMMAX
C
      DO Z2 = 1,NZ12
C
C=======================================================================
C    set up effective interaction matrix  w  and store in AMAT
C
C                 w = sum(IT) conc(IT) * x(IT) * x(IT)
C=======================================================================
C
         CALL CINIT(NDIMCHI*NDIMCHI,AMAT)
C
C================================================================= IQ ==
         LOOP_IQ:DO IQ = IQBOT_CHI,IQTOP_CHI
C
C================================================================= IO ==
            LOOP_IO:DO IO = 1,NOQ(IQ)
C
               IT = ITOQ(IO,IQ)
C
C-----------------------------------------------------------------------
C         thermal lattice vibrations and/or spin fluctuations
C-----------------------------------------------------------------------
C
C============================================================= IFLUCT ==
C                                                 perform local rotation
               IFT = (IT-1)*NFLUCT
               LOOP_IFLUCT:DO IFLUCT = 1,NFLUCT
                  IFT = IFT + 1
C
C----------- keep intermediate t-matrix TSS_FT: t_ft = F_f * t_t * F_f^+
C
                  CALL CMATTRANS(N,M,WKM1,FMAT_FT(1,1,IFT),.TRUE.,
     &                           TSSTA(1,1,IT),'SAS+',TSS_FTA)
C
                  CALL CMATTRANS(N,M,WKM1,FMAT_FT(1,1,IFT),.TRUE.,
     &                           TSSTB(1,1,IT),'SAS+',TSS_FTB)
C
C============================================================= IVIBRA ==
C                                             perform local displacement
                  IVT = (IT-1)*NVIBRA
                  LOOP_IVIBRA:DO IVIBRA = 1,NVIBRA
                     IVT = IVT + 1
C
                     IVFT = (IT-1)*NVIBRA*NFLUCT + (IVIBRA-1)
     &                      *NFLUCT + IFLUCT
C
C-------------------------------- TSS_VFT: t_vft = U_v * t_ft * U_v^(-1)
C
                     CALL CMAT_U_TRANS(UMAT_VTA(1,1,IVT),TSS_FTA,'UAUT',
     &                                 TSS_VFTA)
C
                     CALL CMAT_U_TRANS(UMAT_VTB(1,1,IVT),TSS_FTB,'UAUT',
     &                                 TSS_VFTB)
C
C-------------------------------------------- MSS_VFT: m_vft = 1 / t_vft
C
                     CALL CMATINV3(N,M,IPIVKM,TSS_VFTA,WKM1,MSS_VFTA)
C
                     CALL CMATINV3(N,M,IPIVKM,TSS_VFTB,WKM1,MSS_VFTB)
C
C---------------- get projection matrices DMATTG, DTILTG in global frame
C
                     CALL GETDMAT(TAUQA(1,1,IQ),DMAT_VFTA,WKM1,DMSSA,N,
     &                            MSSQA(1,1,IQ),MSS_VFTA,M)
C
                     CALL GETDMAT(TAUQB(1,1,IQ),DMAT_VFTB,WKM1,DMSSB,N,
     &                            MSSQB(1,1,IQ),MSS_VFTB,M)
C
C-----------------------------------------------------------------------
C    set up auxilary x-matrix  XMAT  for each atom type  IT
C
C                x(IT) = [ m(IT) - m(CPA) ] * D(IT)
C-----------------------------------------------------------------------
C
                     CALL CMATMUL(N,M,DMSSA,DMAT_VFTA,XMAT_VFTA)
                     CALL CMATMUL(N,M,DMSSB,DMAT_VFTB,XMAT_VFTB)
C
                     DO I2 = 1,NKM
                        DO I1 = 1,NKM
                           XMAT_VFTB(I1,I2,2)
     &                        = DCONJG(XMAT_VFTB(I2,I1,1))
                        END DO
                     END DO
C
C-----------------------------------------------------------------------
C
                     K1 = (IQ-IQBOT_CHI)*NKMSQ
                     DO L1 = 1,NKM
                        DO L4 = 1,NKM
                           K1 = K1 + 1
C
                           K2 = (IQ-IQBOT_CHI)*NKMSQ
                           DO L2 = 1,NKM
                              DO L3 = 1,NKM
                                 K2 = K2 + 1
C
                                 CXAXB = X_VFT(IVFT)*XMAT_VFTA(L1,L2)
     &                              *XMAT_VFTB(L3,L4,Z2)
C
                                 AMAT(K1,K2) = AMAT(K1,K2) + CXAXB
C
                              END DO
                           END DO
                        END DO
                     END DO
C
                  END DO LOOP_IVIBRA
C============================================================= IVIBRA ==
C
               END DO LOOP_IFLUCT
C============================================================= IFLUCT ==
C
            END DO LOOP_IO
C================================================================= IO ==
C
         END DO LOOP_IQ
C================================================================= IQ ==
C
C-----------------------------------------------------------------------
C
C        CALL ZGEMM('N','N',NDIMCHI,NDIMCHI,NDIMCHI,C1,CHIZ(1,1,Z2),
C                   NDIMCHI,AMAT,NDIMCHI,C0,BMAT,NDIMCHI)
C
         CALL CMATMUL(NDIMCHI,NDIMCHI,CHIZ(1,1,Z2),AMAT,BMAT)
C
         DO I = 1,NDIMCHI
            BMAT(I,I) = -1D0 + BMAT(I,I)
         END DO
C
C        CALL CINVLU(BMAT,AMAT,NDIMCHI,NDIMCHI)
C
         LWORK = NDIMCHI*NDIMCHI
C
         CALL ZGETF2(NDIMCHI,NDIMCHI,BMAT,NDIMCHI,IPIV,INFO)
         CALL ZGETRI(NDIMCHI,BMAT,NDIMCHI,IPIV,AMAT,LWORK,INFO)
C
         CALL CMATMUL(NDIMCHI,NDIMCHI,BMAT,CHIZ(1,1,Z2),AMAT)
C
         DO J = 1,NDIMCHI
            DO I = 1,NDIMCHI
               CHIZ(I,J,Z2) = -AMAT(I,J)
            END DO
         END DO
C
      END DO
C
      DEALLOCATE (AMAT,BMAT,DMSSA,DMSSB,XMAT_VFTA,XMAT_VFTB)
C
      END
