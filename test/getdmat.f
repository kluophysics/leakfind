C*==getdmat.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE GETDMAT(TAUQ,DMATT,DTILT,DM,N,MSSQ,MSST,M)
C   ********************************************************************
C   *                                                                  *
C   *   calculate projection matrices   DMATT  and  DTILT              *
C   *   for preselected atom type  IT  on preselected site  IQ         *
C   *                                                                  *
C   *    DM = m(t)-m(c)                                                *
C   *                                                                  *
C   *    D~(t) = ( 1 + ( m(t) - m(c) ) * TAU )**(-1)                   *
C   *    D(t)  = ( 1 + TAU * ( m(t) - m(c) ) )**(-1)                   *
C   *                                                                  *
C   *   NOTE:   for IREL <> 2:     N = NKM   M = NKMMAX   N_Q = NKMQ   *
C   *           for IREL =  2:     N = NLM   M = NLMMAX   N_Q = NLMQ   *
C   *                           or N = NKM   M = NKMMAX   N_Q = NKMQ   *
C   *                           depending on calling routine           *
C   *                                                                  *
C   *           the working arrays  WKM1,IPIVKM  will always fit       *
C   *                                                                  *
C   * 14/08/2014 HE                                                    *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:WKM1,IPIVKM
      USE MOD_CONSTANTS,ONLY:C0,C1
      IMPLICIT NONE
C*--GETDMAT26
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER M,N
      COMPLEX*16 DM(M,M),DMATT(M,M),DTILT(M,M),MSSQ(M,M),MSST(M,M),
     &           TAUQ(M,M)
C
C Local variables
C
      INTEGER I,INFO,J
C
C*** End of declarations rewritten by SPAG
C
      DO J = 1,N
         DO I = 1,N
            DM(I,J) = MSST(I,J) - MSSQ(I,J)
         END DO
      END DO
C
C     -------------------------------------------
C                   ( m(t) - m(c) ) * TAU
C     -------------------------------------------
      CALL ZGEMM('N','N',N,N,N,C1,DM,M,TAUQ,M,C0,DTILT,M)
      CALL ZGEMM('N','N',N,N,N,C1,TAUQ,M,DM,M,C0,DMATT,M)
C
C     -------------------------------------------
C               1 + ( m(t) - m(c) ) * TAU
C     -------------------------------------------
      DO I = 1,N
         DTILT(I,I) = C1 + DTILT(I,I)
         DMATT(I,I) = C1 + DMATT(I,I)
      END DO
C
C     -------------------------------------------
C     D~(t) = ( 1 + ( m(t) - m(c) ) * TAU )**(-1)
C     -------------------------------------------
      CALL ZGETRF(N,N,DTILT,M,IPIVKM,INFO)
      CALL ZGETRI(N,DTILT,M,IPIVKM,WKM1,M*M,INFO)
C
      CALL ZGETRF(N,N,DMATT,M,IPIVKM,INFO)
      CALL ZGETRI(N,DMATT,M,IPIVKM,WKM1,M*M,INFO)
C
      END
C*==get_dtilt.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE GET_DTILT(MSST,MSSQ,TAUQ,DTILT)
C   ********************************************************************
C   *                                                                  *
C   *   calculate projection matrix     DTILT                          *
C   *   for preselected atom type  IT  on preselected site  IQ         *
C   *                                                                  *
C   *    D~(t) = ( 1 + ( m(t) - m(c) ) * TAU )**(-1)                   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,WKM1,WKM2,IPIVKM
      USE MOD_CONSTANTS,ONLY:C1,C0
      IMPLICIT NONE
C*--GET_DTILT97
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 DTILT(NKMMAX,NKMMAX),MSSQ(NKMMAX,NKMMAX),
     &           MSST(NKMMAX,NKMMAX),TAUQ(NKMMAX,NKMMAX)
C
C Local variables
C
      INTEGER I,M,N
C
C*** End of declarations rewritten by SPAG
C
      N = NKM
      M = NKMMAX
C
      WKM1(1:N,1:N) = MSST(1:N,1:N) - MSSQ(1:N,1:N)
C
C     -------------------------------------------
C                   ( m(t) - m(c) ) * TAU
C     -------------------------------------------
C
      CALL ZGEMM('N','N',N,N,N,C1,WKM1,M,TAUQ,M,C0,DTILT,M)
C
C     -------------------------------------------
C               1 + ( m(t) - m(c) ) * TAU
C     -------------------------------------------
C
      DO I = 1,N
         DTILT(I,I) = C1 + DTILT(I,I)
      END DO
C
C     -------------------------------------------
C     D~(t) = ( 1 + ( m(t) - m(c) ) * TAU )**(-1)
C     -------------------------------------------
C
      CALL CMATINV2(N,M,IPIVKM,WKM2,DTILT)
C
      END
C*==get_dmatt.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE GET_DMATT(MSST,MSSQ,TAUQ,DMATT)
C   ********************************************************************
C   *                                                                  *
C   *   calculate projection matrix     DMATT                          *
C   *   for preselected atom type  IT  on preselected site  IQ         *
C   *                                                                  *
C   *    D(t)  = ( 1 + TAU * ( m(t) - m(c) ) )**(-1)                   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,WKM1,WKM2,IPIVKM
      USE MOD_CONSTANTS,ONLY:C1,C0
      IMPLICIT NONE
C*--GET_DMATT163
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 DMATT(NKMMAX,NKMMAX),MSSQ(NKMMAX,NKMMAX),
     &           MSST(NKMMAX,NKMMAX),TAUQ(NKMMAX,NKMMAX)
C
C Local variables
C
      INTEGER I,M,N
C
C*** End of declarations rewritten by SPAG
C
      N = NKM
      M = NKMMAX
C
      WKM1(1:N,1:N) = MSST(1:N,1:N) - MSSQ(1:N,1:N)
C
C     -------------------------------------------
C                   TAU * ( m(t) - m(c) )
C     -------------------------------------------
C
      CALL ZGEMM('N','N',N,N,N,C1,TAUQ,M,WKM1,M,C0,DMATT,M)
C
C     -------------------------------------------
C               1 + TAU * ( m(t) - m(c) )
C     -------------------------------------------
C
      DO I = 1,N
         DMATT(I,I) = C1 + DMATT(I,I)
      END DO
C
C     -------------------------------------------
C     D(t) = ( 1 + TAU * ( m(t) - m(c) ) )**(-1)
C     -------------------------------------------
C
      CALL CMATINV2(N,M,IPIVKM,WKM2,DMATT)
C
      END
C*==get_taut.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE GET_TAUT(MSST,MSSQ,TAUQ,TAUT)
C   ********************************************************************
C   *                                                                  *
C   *   calculate procected TAU matrix  TAUT                           *
C   *   for preselected atom type  IT  on preselected site  IQ         *
C   *                                                                  *
C   *    TAU(t) =  TAU * D~(t)                                         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,WKM1,WKM2,WKM3,IPIVKM
      USE MOD_CONSTANTS,ONLY:C1,C0
      IMPLICIT NONE
C*--GET_TAUT229
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX),MSST(NKMMAX,NKMMAX),
     &           TAUQ(NKMMAX,NKMMAX),TAUT(NKMMAX,NKMMAX)
C
C Local variables
C
      INTEGER I,M,N
C
C*** End of declarations rewritten by SPAG
C
      N = NKM
      M = NKMMAX
C
      WKM1(1:N,1:N) = MSST(1:N,1:N) - MSSQ(1:N,1:N)
C
C     -------------------------------------------
C                   ( m(t) - m(c) ) * TAU
C     -------------------------------------------
C
      CALL ZGEMM('N','N',N,N,N,C1,WKM1,M,TAUQ,M,C0,WKM3,M)
C
C     -------------------------------------------
C               1 + ( m(t) - m(c) ) * TAU
C     -------------------------------------------
C
      DO I = 1,N
         WKM3(I,I) = C1 + WKM3(I,I)
      END DO
C
C     -------------------------------------------
C     D~(t) = ( 1 + ( m(t) - m(c) ) * TAU )**(-1)
C     -------------------------------------------
C
      CALL CMATINV2(N,M,IPIVKM,WKM2,WKM3)
C
C     -------------------------------------------
C              TAU(t) = TAU * D~(t)
C     -------------------------------------------
C
      CALL ZGEMM('N','N',N,N,N,C1,TAUQ,M,WKM3,M,C0,TAUT,M)
C
      END
C*==get_matz.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE GET_MATZ(MAT,MATZ)
C   ********************************************************************
C   *                                                                  *
C   *     perform transformation for z -> z*                           *
C   *                                                                  *
C   *                     M(z*) = L M(z)^+ L                           *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,WKM1,LMAT3
      USE MOD_CONSTANTS,ONLY:C1,C0
      IMPLICIT NONE
C*--GET_MATZ301
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 MAT(NKMMAX,NKMMAX),MATZ(NKMMAX,NKMMAX,2)
C
C Local variables
C
      INTEGER M,N
C
C*** End of declarations rewritten by SPAG
C
      N = NKM
      M = NKMMAX
C
C-------------------------------------------------------------------- Z*
C
      CALL ZGEMM('N','C',N,N,N,C1,LMAT3,M,MAT,M,C0,WKM1,M)
      CALL ZGEMM('N','N',N,N,N,C1,WKM1,M,LMAT3,M,C0,MATZ(1,1,1),M)
C
C--------------------------------------------------------------------- Z
C
      MATZ(:,:,2) = MAT(:,:)
C
      END
C*==get_matz_umat.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE GET_MATZ_UMAT(MAT,MATZ)
C   ********************************************************************
C   *                                                                  *
C   *  perform transformation for z -> z*                              *
C   *                                                                  *
C   *                     U(z*) =  U^+(z)                              *
C   *                                                                  *
C   * DK, 2014-05                                                      *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKMMAX
      IMPLICIT NONE
C*--GET_MATZ_UMAT352
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 MAT(NKMMAX,NKMMAX),MATZ(NKMMAX,NKMMAX,2)
C
C*** End of declarations rewritten by SPAG
C
      MATZ(:,:,1) = TRANSPOSE(DCONJG(MAT(:,:)))
      MATZ(:,:,2) = MAT(:,:)
C
      END
C
