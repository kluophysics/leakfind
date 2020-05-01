C*==negfkloop_stt.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE NEGFKLOOP_STT(GLESQ_TR,NKKR_TB,DELTA,LMAT2,IK,IKTOP,WK,
     &                         WKSUM,KVEC,TRANS,STT)
C   ********************************************************************
C   *                                                                  *
C   *  Calculate the Spin Transfer Torque from lesser Green's function *
C   *  of the transporting electrons G^<_tr and the matrix elements of *
C   *  the exchange splitting Delta_xc^nrel:                           *
C   *                                                                  *
C   *    - first rotate GLESQX_TR to the local FOR                     *
C   *    - then calculate the magnetisation m_tr(k_II,E_F) as          *
C   *      Trace_spin sigma_x/y/z G^<_tr(k_II,E_F)                     *
C   *    - then evaluate the vector product Delta_xc x m_tr = Tq_loc   *
C   *    - integrate over k_II and rotate back to global FOR           *
C   *                                                                  *
C   *  The result should be the layer-resolved STT for a given THETA   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,NKMQ,NLMMAX,NXM,WKM1
      USE MOD_CALCMODE,ONLY:IREL,KMROT
      USE MOD_SITES,ONLY:DROTQ,IQBOT_TB,IQTOP_TB,NQTB
      USE MOD_TYPES,ONLY:NTMAX
      IMPLICIT NONE
C*--NEGFKLOOP_STT25
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IK,IKTOP,NKKR_TB
      LOGICAL STT
      REAL*8 WK,WKSUM
      COMPLEX*16 DELTA(NLMMAX,NLMMAX,NTMAX),GLESQ_TR(NXM,NXM,NQTB),
     &           LMAT2(NXM,NXM),TRANS(3)
      REAL*8 KVEC(3)
C
C Local variables
C
      COMPLEX*16 DROTQ2(:,:,:),GLESQ_TR_LOC(:,:,:)
      INTEGER IQ,IQTB,M,N
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DROTQ2,GLESQ_TR_LOC
C
      ALLOCATE (DROTQ2(NXM,NXM,NQTB))
      ALLOCATE (GLESQ_TR_LOC(NXM,NXM,NQTB))
C
      IF ( STT .AND. IK.EQ.1 ) WRITE (6,*) 
     &                       '<NEGFKLOOP_STT>:   UNDER CONSTRUCTION !!!'
C
      DO IQ = IQBOT_TB,IQTOP_TB
C
         IQTB = IQ - IQBOT_TB + 1
C
C rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
C
C   ********************************************************************
C   *                                                                  *
C   *  rotate the "transporting electrons''s" lesser Green's function  *
C   *  G^<_tr(k_II,E_F) (GLESQ_TR) from the global to the local FOR    *
C   *                                                                  *
C   *  Since GLESQ_TR is in the l,m_l,m_s-representation and the       *
C   *  rotation matrices DROTQ are in kappa,mue, the latter have to    *
C   *  be transformed first (-> BASTRANS !? or rather CHANGEREP)       *
C   *                                                                  *
C   ********************************************************************
C
C
         IF ( KMROT.NE.0 ) THEN
C
            M = NKMMAX
C
            IF ( IREL.NE.2 ) THEN
               N = NKMQ(IQ)
            ELSE
               N = NKM
            END IF
C
            CALL CHANGEREP(NKM,NKMMAX,DROTQ(1,1,IQ),'REL>RLM',
     &                     DROTQ2(1,1,IQ))
C
            CALL ROTATE(GLESQ_TR(1,1,IQTB),'G->L',WKM1,N,DROTQ2(1,1,IQ),
     &                  M)
            GLESQ_TR_LOC(1:N,1:N,IQTB) = WKM1(1:N,1:N)
C
         END IF
C
C rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
C
         IF ( IK.EQ.IKTOP ) THEN
C
            WRITE (6,*) 'IQ = ',IQ
            WRITE (6,*) 'IQTB = ',IQTB
C
            WRITE (6,*) 'G^<_tr(k,E) at IKTOP in IREL2 and global FOR'
            CALL CMATSTRUCT('GLESQ_TR',GLESQ_TR(1,1,IQTB),18,18,2,2,1,
     &                      1.D-8,6)
C
            WRITE (6,*) 'G^<_tr(k,E) at IKTOP in IREL2 and local FOR'
            CALL CMATSTRUCT('GLESQ_TR_LOC',GLESQ_TR_LOC(1,1,IQTB),18,18,
     &                      2,2,1,1.D-8,6)
C
         END IF
C
      END DO
C
      END
