C*==ssite_rotate.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE SSITE_ROTATE(TSST,MSST,SSST)
C   ********************************************************************
C   *                                                                  *
C   *  rotate ALL type dependent ssite-matrices from the               *
C   *  LOCAL to the GLOBAL frame                                       *
C   *                                                                  *
C   *  NO ROTATION in case or thermal fluctuations                     *
C   *  this is done later on in combination with fluctuations          *
C   *                                                                  *
C   *     starting with version 7.4.0 MEZZ and MEZJ are excluded       *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:IREL,KMROT
      USE MOD_THERMAL,ONLY:NFLUCT
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,WKM1,NKMQ
      USE MOD_SITES,ONLY:IQAT,DROTQ,MAGROT_Q
      USE MOD_TYPES,ONLY:ITBOT,ITTOP,NTMAX
      IMPLICIT NONE
C*--SSITE_ROTATE21
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SSITE_ROTATE')
C
C Dummy arguments
C
      COMPLEX*16 MSST(NKMMAX,NKMMAX,NTMAX),SSST(NKMMAX,NKMMAX,NTMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      INTEGER IQ,IT,M,N
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
      IF ( KMROT.EQ.0 .OR. NFLUCT.GT.1 ) RETURN
C
      M = NKMMAX
C
      LOOP_IT:DO IT = ITBOT,ITTOP
C
         IQ = IQAT(1,IT)
C
         IF ( .NOT.MAGROT_Q(IQ) ) CYCLE LOOP_IT
C
         IF ( IREL.NE.2 ) THEN
            N = NKMQ(IQ)
         ELSE
            N = NKM
         END IF
C
C------------------------------------------------------------------- TSST
         CALL ROTATE(TSST(1,1,IT),'L->G',WKM1,N,DROTQ(1,1,IQ),M)
         TSST(1:N,1:N,IT) = WKM1(1:N,1:N)
C
C------------------------------------------------------------------- MSST
         CALL ROTATE(MSST(1,1,IT),'L->G',WKM1,N,DROTQ(1,1,IQ),M)
         MSST(1:N,1:N,IT) = WKM1(1:N,1:N)
C
C------------------------------------------------------------------- SSST
         CALL ROTATE(SSST(1,1,IT),'L->G',WKM1,N,DROTQ(1,1,IQ),M)
         SSST(1:N,1:N,IT) = WKM1(1:N,1:N)
C
      END DO LOOP_IT
C
      END
C*==me_rotate.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE ME_ROTATE(IROT,NROT,DROT,MROT,M_LOC,M_GLO)
C   ********************************************************************
C   *                                                                  *
C   *  rotate the vector matrix elements  M_LOC --> M_GLO              *
C   *  from a LOCAL to the GLOBAL frame of reference                   *
C   *  using the rotation matrices  DROT  and MROT                     *
C   *                                                                  *
C   *  NOTE:  the polarisation refers to CARTESIAN coordinates         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM
      IMPLICIT NONE
C*--ME_ROTATE103
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IROT,NROT
      COMPLEX*16 DROT(NKMMAX,NKMMAX,NROT),M_GLO(NKMMAX,NKMMAX,3),
     &           M_LOC(NKMMAX,NKMMAX,3)
      REAL*8 MROT(3,3,NROT)
C
C Local variables
C
      INTEGER IPOL,JPOL,M,N
      COMPLEX*16 WRK(:,:,:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE WRK
C
      ALLOCATE (WRK(NKMMAX,NKMMAX,3))
C
      M = NKMMAX
      N = NKM
C
C-----------------------------------------------------------------------
C               transform w.r.t. angular momentum indices
C-----------------------------------------------------------------------
C
      DO IPOL = 1,3
         CALL ROTATE(M_LOC(1,1,IPOL),'L->G',WRK(1,1,IPOL),N,
     &               DROT(1,1,IROT),M)
      END DO
C
C-----------------------------------------------------------------------
C                     transform polarisation
C-----------------------------------------------------------------------
C
      M_GLO(:,:,:) = C0
C
      DO IPOL = 1,3
         DO JPOL = 1,3
            M_GLO(1:N,1:N,IPOL) = M_GLO(1:N,1:N,IPOL)
     &                            + MROT(JPOL,IPOL,IROT)
     &                            *WRK(1:N,1:N,JPOL)
         END DO
      END DO
C
      END
C*==cmat_u_trans.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE CMAT_U_TRANS(U,A,KEY,B)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the matrix-matrix transformation for displacement     *
C   *                                                                  *
C   *     B = U      * A * U^(-1)    for  KEY = 'UAUT'                 *
C   *     B = U^(-1) * A * U         for  KEY = 'UTAU'                 *
C   *                                                                  *
C   *                                                                  *
C   *   A       matrix to be transformed  - unchanged on exit          *
C   *   B       result of transformation                               *
C   *   U       transformation matrix                                  *
C   *                                                                  *
C   *   A,B,U   complex  SQUARE  NKMAX x NKMMAX - matrices             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,WKM1,RREL
      USE MOD_CALCMODE,ONLY:IREL
      IMPLICIT NONE
C*--CMAT_U_TRANS184
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CME_TRANS_U')
C
C Dummy arguments
C
      CHARACTER*4 KEY
      COMPLEX*16 A(NKMMAX,NKMMAX),B(NKMMAX,NKMMAX),U(NKMMAX,NKMMAX)
C
C Local variables
C
      LOGICAL INITIALIZE
      COMPLEX*16 UT(:,:),XM(:,:),XMD(:,:)
      SAVE XM,XMD
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
C
      ALLOCATABLE XM,XMD
      ALLOCATABLE UT
C
      IF ( NKM.NE.NKMMAX ) CALL STOP_MESSAGE(ROUTINE,'NKM <> NKMMAX')
C
      IF ( INITIALIZE ) THEN
         ALLOCATE (XM(NKMMAX,NKMMAX),XMD(NKMMAX,NKMMAX))
         XM = MATMUL(TRANSPOSE(RREL),RREL)
         XMD = TRANSPOSE(DCONJG(XM))
         INITIALIZE = .FALSE.
      END IF
C
      ALLOCATE (UT(NKMMAX,NKMMAX))
C
      IF ( IREL.EQ.3 ) THEN
C
C=======================================================================
C    construct U(-s) = U^(-1) = X^+  U^T  X
C    for relativistic case, kappa-mu representation
C
C    with X = Y^T Y
C    where Y is matrix for tranfomation from LMS to KM representation
C
C=======================================================================
C
         WKM1 = TRANSPOSE(U)
         UT = MATMUL(XMD,MATMUL(WKM1,XM))
C
      ELSE
C
C=======================================================================
C    U(-s) = U^(-1) =  U^T
C    for lms - representation
C=======================================================================
C
         UT = TRANSPOSE(U)
C
      END IF
C
C=======================================================================
C                    perform transformation
C=======================================================================
C
C---------------------------------------------------  B = U * A * U^(-1)
      IF ( KEY.EQ.'UAUT' ) THEN
C
         B = MATMUL(U,MATMUL(A,UT))
C
C---------------------------------------------------  B = U^(-1) * A * U
      ELSE IF ( KEY.EQ.'UTAU' ) THEN
C
         B = MATMUL(UT,MATMUL(A,U))
C
      ELSE
C
         CALL STOP_MESSAGE(ROUTINE,' KEY  not allowed  ')
C
      END IF
C-----------------------------------------------------------------------
C
      END
C*==me_rotate_reg.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE ME_ROTATE_REG(IROT,NROT,DROT,MROT,MAB_LOC,MBA_LOC,
     &                         MAB_GLO,MBA_GLO)
C   ********************************************************************
C   *                                                                  *
C   *  rotate the regular matrix elements  MAB  and  MBA               *
C   *  from a LOCAL to the GLOBAL frame of reference                   *
C   *  using the rotation matrices  DROT  and MROT                     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM
      IMPLICIT NONE
C*--ME_ROTATE_REG297
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IROT,NROT
      COMPLEX*16 DROT(NKMMAX,NKMMAX,NROT),MAB_GLO(NKMMAX,NKMMAX,3),
     &           MAB_LOC(NKMMAX,NKMMAX,3),MBA_GLO(NKMMAX,NKMMAX,3),
     &           MBA_LOC(NKMMAX,NKMMAX,3)
      REAL*8 MROT(3,3,NROT)
C
C Local variables
C
      INTEGER IPOL,JPOL,M,N
      COMPLEX*16 WRK(:,:,:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE WRK
C
      ALLOCATE (WRK(NKMMAX,NKMMAX,3))
C
      M = NKMMAX
      N = NKM
C
C------------------------------------------------------------------- MAB
C
      DO IPOL = 1,3
         CALL ROTATE(MAB_LOC(1,1,IPOL),'L->G',WRK(1,1,IPOL),N,
     &               DROT(1,1,IROT),M)
      END DO
C
      MAB_GLO(:,:,:) = C0
C
      DO IPOL = 1,3
         DO JPOL = 1,3
            MAB_GLO(1:N,1:N,IPOL) = MAB_GLO(1:N,1:N,IPOL)
     &                              + MROT(JPOL,IPOL,IROT)
     &                              *WRK(1:N,1:N,JPOL)
         END DO
      END DO
C
C------------------------------------------------------------------- MBA
C
      DO IPOL = 1,3
         CALL ROTATE(MBA_LOC(1,1,IPOL),'L->G',WRK(1,1,IPOL),N,
     &               DROT(1,1,IROT),M)
      END DO
C
      MBA_GLO(:,:,:) = C0
C
      DO IPOL = 1,3
         DO JPOL = 1,3
            MBA_GLO(1:N,1:N,IPOL) = MBA_GLO(1:N,1:N,IPOL)
     &                              + MROT(JPOL,IPOL,IROT)
     &                              *WRK(1:N,1:N,JPOL)
         END DO
      END DO
C
      END
C*==me_rotate_irr.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE ME_ROTATE_IRR(IROT,NROT,DROT,MROT,MIRR2_LOC,MIRR3_LOC,
     &                         MIRR4_LOC,MIRR2_GLO,MIRR3_GLO,MIRR4_GLO)
C   ********************************************************************
C   *                                                                  *
C   *  rotate the irrregular matrix elements  MIRR                     *
C   *  from a LOCAL to the GLOBAL frame of reference                   *
C   *  using the rotation matrices  DROT  and MROT                     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM
      IMPLICIT NONE
C*--ME_ROTATE_IRR383
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IROT,NROT
      COMPLEX*16 DROT(NKMMAX,NKMMAX,NROT),MIRR2_GLO(NKMMAX,NKMMAX,3,3),
     &           MIRR2_LOC(NKMMAX,NKMMAX,3,3),
     &           MIRR3_GLO(NKMMAX,NKMMAX,3,3),
     &           MIRR3_LOC(NKMMAX,NKMMAX,3,3),
     &           MIRR4_GLO(NKMMAX,NKMMAX,3,3),
     &           MIRR4_LOC(NKMMAX,NKMMAX,3,3)
      REAL*8 MROT(3,3,NROT)
C
C Local variables
C
      INTEGER IPOL1,IPOL2,JPOL1,JPOL2,M,N
      COMPLEX*16 WRK(:,:,:,:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE WRK
C
      ALLOCATE (WRK(NKMMAX,NKMMAX,3,3))
C
      M = NKMMAX
      N = NKM
C
C-----------------------------------------------------------------------
C               transform w.r.t. angular momentum indices          MIRR2
C
C               map:   MIRR2      -->       D MIRR2 D^+
C
C-----------------------------------------------------------------------
C
      DO IPOL1 = 1,3
         DO IPOL2 = 1,3
            CALL ROTATE(MIRR2_LOC(1,1,IPOL1,IPOL2),'L->G',
     &                  WRK(1,1,IPOL1,IPOL2),N,DROT(1,1,IROT),M)
         END DO
      END DO
C
C-----------------------------------------------------------------------
C                     transform polarisation
C-----------------------------------------------------------------------
C
      MIRR2_GLO(:,:,:,:) = C0
C
      DO IPOL1 = 1,3
         DO IPOL2 = 1,3
C
            DO JPOL1 = 1,3
               DO JPOL2 = 1,3
                  MIRR2_GLO(1:N,1:N,IPOL1,IPOL2)
     &               = MIRR2_GLO(1:N,1:N,IPOL1,IPOL2)
     &               + MROT(JPOL1,IPOL1,IROT)*MROT(JPOL2,IPOL2,IROT)
     &               *WRK(1:N,1:N,JPOL1,JPOL2)
C
               END DO
            END DO
C
         END DO
      END DO
C-----------------------------------------------------------------------
C               transform w.r.t. angular momentum indices          MIRR3
C
C               map:   MIRR3      -->       D MIRR3 D^+
C
C-----------------------------------------------------------------------
C
      DO IPOL1 = 1,3
         DO IPOL2 = 1,3
            CALL ROTATE(MIRR3_LOC(1,1,IPOL1,IPOL2),'L->G',
     &                  WRK(1,1,IPOL1,IPOL2),N,DROT(1:M,1:M,IROT),M)
         END DO
      END DO
C
C-----------------------------------------------------------------------
C                     transform polarisation
C-----------------------------------------------------------------------
C
      MIRR3_GLO(:,:,:,:) = C0
C
      DO IPOL1 = 1,3
         DO IPOL2 = 1,3
C
            DO JPOL1 = 1,3
               DO JPOL2 = 1,3
                  MIRR3_GLO(1:N,1:N,IPOL1,IPOL2)
     &               = MIRR3_GLO(1:N,1:N,IPOL1,IPOL2)
     &               + MROT(JPOL1,IPOL1,IROT)*MROT(JPOL2,IPOL2,IROT)
     &               *WRK(1:N,1:N,JPOL1,JPOL2)
C
               END DO
            END DO
C
         END DO
      END DO
C
C-----------------------------------------------------------------------
C               transform w.r.t. angular momentum indices          MIRR4
C-----------------------------------------------------------------------
C
      DO IPOL1 = 1,3
         DO IPOL2 = 1,3
            CALL ROTATE(MIRR4_LOC(1,1,IPOL1,IPOL2),'L->G',
     &                  WRK(1,1,IPOL1,IPOL2),N,DROT(1,1,IROT),M)
         END DO
      END DO
C
C-----------------------------------------------------------------------
C                     transform polarisation
C-----------------------------------------------------------------------
C
      MIRR4_GLO(:,:,:,:) = C0
C
      DO IPOL1 = 1,3
         DO IPOL2 = 1,3
C
            DO JPOL1 = 1,3
               DO JPOL2 = 1,3
                  MIRR4_GLO(1:N,1:N,IPOL1,IPOL2)
     &               = MIRR4_GLO(1:N,1:N,IPOL1,IPOL2)
     &               + MROT(JPOL1,IPOL1,IROT)*MROT(JPOL2,IPOL2,IROT)
     &               *WRK(1:N,1:N,JPOL1,JPOL2)
               END DO
            END DO
C
         END DO
      END DO
C
      END
C*==me_utrans_reg.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE ME_UTRANS_REG(UMATA,UMATB,MAB_LOC,MBA_LOC,MAB_GLO,
     &                         MBA_GLO)
C   ********************************************************************
C   *                                                                  *
C   *  rotate the regular matrix elements  MAB  and  MBA               *
C   *  from a LOCAL to the GLOBAL frame of reference                   *
C   *  using the rotation matrices  UMATA  and  UMATB                  *
C   *                                                                  *
C   *  NOTE: the polarisation does not change unter the U-trans.       *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0,C1
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,WKM1,RREL
      USE MOD_CALCMODE,ONLY:IREL
      IMPLICIT NONE
C*--ME_UTRANS_REG544
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 MAB_GLO(NKMMAX,NKMMAX,3),MAB_LOC(NKMMAX,NKMMAX,3),
     &           MBA_GLO(NKMMAX,NKMMAX,3),MBA_LOC(NKMMAX,NKMMAX,3),
     &           UMATA(NKMMAX,NKMMAX),UMATB(NKMMAX,NKMMAX)
C
C Local variables
C
      LOGICAL INITIALIZE
      INTEGER IPOL,M,N
      COMPLEX*16 UMATINV(:,:),XM(:,:),XMD(:,:)
      SAVE XM,XMD
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
C
      ALLOCATABLE XM,XMD,UMATINV
C
      M = NKMMAX
      N = NKM
C
      ALLOCATE (UMATINV(M,M))
C
      IF ( INITIALIZE .AND. IREL.EQ.3 ) THEN
         ALLOCATE (XM(M,M),XMD(M,M))
         XM = MATMUL(TRANSPOSE(RREL),RREL)
         XMD = TRANSPOSE(DCONJG(XM))
         INITIALIZE = .FALSE.
      END IF
C
      DO IPOL = 1,3
C
C------------------------------------------------------------------- MAB
C
         IF ( IREL.EQ.3 ) THEN
            UMATINV = MATMUL(XMD,MATMUL(TRANSPOSE(UMATB),XM))
         ELSE
            UMATINV = TRANSPOSE(UMATB)
         END IF
C
         CALL ZGEMM('N','N',N,N,N,C1,UMATA,M,MAB_LOC(1,1,IPOL),M,C0,
     &              WKM1,M)
C
         CALL ZGEMM('N','N',N,N,N,C1,WKM1,M,UMATINV,M,C0,
     &              MAB_GLO(1,1,IPOL),M)
C
C------------------------------------------------------------------- MBA
C
         IF ( IREL.EQ.3 ) THEN
            UMATINV = MATMUL(XMD,MATMUL(TRANSPOSE(UMATA),XM))
         ELSE
            UMATINV = TRANSPOSE(UMATA)
         END IF
C
         CALL ZGEMM('N','N',N,N,N,C1,UMATB,M,MBA_LOC(1,1,IPOL),M,C0,
     &              WKM1,M)
C
         CALL ZGEMM('N','N',N,N,N,C1,WKM1,M,UMATINV,M,C0,
     &              MBA_GLO(1,1,IPOL),M)
C
      END DO
C
      END
C*==me_utrans_irr.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE ME_UTRANS_IRR(UMATA,UMATB,MIRR2_LOC,MIRR3_LOC,
     &                         MIRR4_LOC,MIRR2_GLO,MIRR3_GLO,MIRR4_GLO)
C   ********************************************************************
C   *                                                                  *
C   *  rotate the regular matrix elements  MAB  and  MBA               *
C   *  from a LOCAL to the GLOBAL frame of reference                   *
C   *  using the rotation matrices  UMATA  and  UMATB                  *
C   *                                                                  *
C   *  NOTE: the polarisation does not change under the U-trans.       *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKMMAX
      IMPLICIT NONE
C*--ME_UTRANS_IRR638
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 MIRR2_GLO(NKMMAX,NKMMAX,3,3),
     &           MIRR2_LOC(NKMMAX,NKMMAX,3,3),
     &           MIRR3_GLO(NKMMAX,NKMMAX,3,3),
     &           MIRR3_LOC(NKMMAX,NKMMAX,3,3),
     &           MIRR4_GLO(NKMMAX,NKMMAX,3,3),
     &           MIRR4_LOC(NKMMAX,NKMMAX,3,3),UMATA(NKMMAX,NKMMAX),
     &           UMATB(NKMMAX,NKMMAX)
C
C Local variables
C
      INTEGER IPOL,JPOL
C
C*** End of declarations rewritten by SPAG
C
C=======================================================================
      DO IPOL = 1,3
         DO JPOL = 1,3
C
C=======================================================================
C----------------------------------------------------------------- MIRR2
C
C      MIRR2 is used within S2 as   MIRR2(E_A,E_B) * TAU(E_B)
C      accordingly the transformation has to be
C
C               MIRR2  <-  U(E_B) * MIRR2(E_A,E_B) * U(E_B)^(-1)
C
            CALL CMAT_U_TRANS(UMATB,MIRR2_LOC(1,1,IPOL,JPOL),'UAUT',
     &                        MIRR2_GLO(1,1,IPOL,JPOL))
C
C----------------------------------------------------------------- MIRR3
C
C      MIRR3 is used within S3 as    TAU(E_A) * MIRR3(E_A,E_B)
C      accordingly the transformation has to be
C
C               MIRR3  <-  U(E_A) * MIRR3(E_A,E_B) * U(E_A)^(-1)
C
            CALL CMAT_U_TRANS(UMATA,MIRR3_LOC(1,1,IPOL,JPOL),'UAUT',
     &                        MIRR3_GLO(1,1,IPOL,JPOL))
C
C----------------------------------------------------------------- MIRR4
C
C     as the trace is taken as a next step no transformtion is applied
C
            MIRR4_GLO(:,:,IPOL,JPOL) = MIRR4_LOC(:,:,IPOL,JPOL)
C
C=======================================================================
C
         END DO
      END DO
C=======================================================================
C
      END
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
