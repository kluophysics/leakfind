C*==init_mod_str.f    processed by SPAG 6.70Rc at 16:42 on 30 Nov 2016
      SUBROUTINE INIT_MOD_STR(NGRLMAX,NRDLMAX,NQQP_STRLIM,LRGNT123_ARG,
     &                        LRGNT12_ARG,LLARR_ARG,NLLMMMAX_ARG,
     &                        LLMAX_ARG,MMLLMAX_ARG,J13MAX_ARG,
     &                        J22MAX_ARG,J22MAX_UPPER_LIMIT,
     &                        NQQP_STR_ARG,NGRL_ARG,NSDL_ARG,NRDL_ARG,
     &                        ALPHA0_ARG,GMAXSQ_ARG,ETA_ARG,CIPWL_ARG,
     &                        RGNT_ARG,NRGNT_ARG,IRGNT_ARG,SRREL_ARG,
     &                        NRREL_ARG,IRREL_ARG,HP_ARG,G123MAX_ARG,
c    &                        R123MAX_ARG,NIJQ_ARG,IJQ_ARG,QQPX_ARG,
     &                        R123MAX_ARG,QQPX_ARG,
c    &                        QQPY_ARG,QQPZ_ARG,SMAX_ARG,R1_ARG,R2_ARG,
     &                        QQPY_ARG,QQPZ_ARG,R1_ARG,R2_ARG,
c    &                        R3_ARG,G1_ARG,G2_ARG,G3_ARG,EXPGNQ_ARG,
     &                        R3_ARG,G1_ARG,G2_ARG,G3_ARG)
c    &                        INDR_ARG,QQMLRS_ARG,GGJLRS_ARG)
C   ********************************************************************
C   *                                                                  *
C   *  called from <STRINIT> to transfer primary structure constants   *
C   *  data, indicated by "_ARG" in the argument list,                 *
C   *  to module  MOD_STR                                              *
C   *                                                                  *
C   *  initialize in addition:                                         *
C   *  D1TERM3, IILERS, DLLMMKE, CHP, SHP, WK, CILMAT                  *
C   *  PWEX2K1, PWEX2K2, PWEX2K3, PWEXIK1, PWEXIK2, PWEXIK3            *
C   *                                                                  *
C   *  NOTE: some big local integers expected, therefore: INTEGER*8    *
C   *                                                                  *
C   ********************************************************************
C
C*** Start of declarations rewritten by SPAG
C
      USE MOD_ANGMOM,ONLY:NLMAX,NKMMAX
      USE MOD_STR,ONLY:LRGNT123,LRGNT12,LLARR,NLLMMMAX,LLMAX,MMLLMAX,
     &    J13MAX,J22MAX,NQQP_STR,NGRL,NSDL,NRDL,ALPHA0,GMAXSQ,ETA,RGNT,
c    &    NRGNT,IRGNT,SRREL,NRREL,IRREL,HP,G123MAX,R123MAX,NIJQ,IJQ,
     &    NRGNT,IRGNT,SRREL,NRREL,IRREL,HP,G123MAX,R123MAX,
c    &    QQPX,QQPY,QQPZ,SMAX,R1,R2,R3,G1,G2,G3,EXPGNQ,INDR,QQMLRS,
     &    QQPX,QQPY,QQPZ,R1,R2,R3,G1,G2,G3,
c    &    GGJLRS,D1TERM3,IILERS,DLLMMKE,CHP,SHP,WK,CILMAT,PWEX2K1,
     &    D1TERM3,IILERS,DLLMMKE,CHP,SHP,WK,CILMAT,PWEX2K1,
     &    PWEX2K2,PWEX2K3,PWEXIK1,PWEXIK2,PWEXIK3,PWP,M1PWL,
     &    USE_NEW_BBDD_VERSION,NQQP_STR_CC,NQQP_STR_RED
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALPHA0_ARG,ETA_ARG,GMAXSQ_ARG
      INTEGER G123MAX_ARG,J13MAX_ARG,J22MAX_ARG,J22MAX_UPPER_LIMIT,
     &        LLARR_ARG,LLMAX_ARG,LRGNT123_ARG,LRGNT12_ARG,MMLLMAX_ARG,
     &        NGRLMAX,NGRL_ARG,NLLMMMAX_ARG,NQQP_STRLIM,NQQP_STR_ARG,
     &        NRDLMAX,NRDL_ARG,NSDL_ARG,R123MAX_ARG
      COMPLEX*16 CIPWL_ARG((2*NLMAX)**2),
c    &           EXPGNQ_ARG(NGRL_ARG,NQQP_STR_ARG),
c    &           QQMLRS_ARG(NLLMMMAX_ARG,NSDL_ARG,NQQP_STR_ARG),
     &           SRREL_ARG(2,2,NKMMAX)
      INTEGER G1_ARG(NGRLMAX),G2_ARG(NGRLMAX),G3_ARG(NGRLMAX),
c    &        IJQ_ARG(NQQP_STRLIM,NQQP_STRLIM),
c    &        INDR_ARG(NSDL_ARG,NQQP_STR_ARG),IRGNT_ARG(LRGNT123_ARG),
     &        IRGNT_ARG(LRGNT123_ARG),
c    &        IRREL_ARG(2,2,NKMMAX),NIJQ_ARG(NQQP_STRLIM),
     &        IRREL_ARG(2,2,NKMMAX),
     &        NRGNT_ARG(LRGNT12_ARG),NRREL_ARG(2,NKMMAX),R1_ARG(NRDLMAX)
c    &        ,R2_ARG(NRDLMAX),R3_ARG(NRDLMAX),SMAX_ARG(NQQP_STRLIM)
     &        ,R2_ARG(NRDLMAX),R3_ARG(NRDLMAX)
c     REAL*8 GGJLRS_ARG(-J22MAX_UPPER_LIMIT:LLARR_ARG,NSDL_ARG,
      REAL*8 HP_ARG(NLLMMMAX_ARG),QQPX_ARG(NQQP_STRLIM),
     &       QQPY_ARG(NQQP_STRLIM),QQPZ_ARG(NQQP_STRLIM),
     &       RGNT_ARG(LRGNT123_ARG)

c     INTEGER NIJQMAX
C
C Local variables
C
      INTEGER*8 IA_ERR,LL,LM1,LM2,MM,MMLL,N,N_INTEGER,N_REAL
      REAL*8 M1PW
C
C*** End of declarations rewritten by SPAG
C
      N_INTEGER = 0
      N_REAL = 0
C
C-----------------------------------------------------------------------
C     deallocate arrays in module MOD_STR if the structure
C     constants are re-initialized
C-----------------------------------------------------------------------
C
      IF ( ALLOCATED(IRGNT) ) THEN
         DEALLOCATE (IRGNT,NRGNT,RGNT,IRREL,SRREL,NRREL,HP,M1PWL)
c        DEALLOCATE (NIJQ,IJQ,SMAX,QQPX,QQPY,QQPZ,R1,R2,R3,G1,G2,G3)
         DEALLOCATE (QQPX,QQPY,QQPZ,R1,R2,R3,G1,G2,G3)
c        DEALLOCATE (EXPGNQ,INDR,QQMLRS,GGJLRS,D1TERM3,DLLMMKE)
         DEALLOCATE (D1TERM3,DLLMMKE)
         DEALLOCATE (PWEX2K1,PWEX2K2,PWEX2K3,CHP,SHP,WK,CILMAT)
         DEALLOCATE (PWEXIK1,PWEXIK2,PWEXIK3,PWP,IILERS)
      END IF
C
C  =====================================================================
C                  copy array dimensions and allocate arrays
C  =====================================================================
C
      LRGNT123 = LRGNT123_ARG
      LRGNT12 = LRGNT12_ARG
      LLARR = LLARR_ARG
      NLLMMMAX = NLLMMMAX_ARG
      LLMAX = LLMAX_ARG
      MMLLMAX = MMLLMAX_ARG
      J13MAX = J13MAX_ARG
      J22MAX = J22MAX_ARG
      NQQP_STR = NQQP_STR_ARG
      NGRL = NGRL_ARG
      NSDL = NSDL_ARG
      NRDL = NRDL_ARG
C
C  =====================================================================
C            allocate arrays of module and copy arguments
C  =====================================================================
C
      ALPHA0 = ALPHA0_ARG
      GMAXSQ = GMAXSQ_ARG
      ETA = ETA_ARG
C
      ALLOCATE (IRGNT(LRGNT123),NRGNT(LRGNT12),RGNT(LRGNT123))
      RGNT(1:LRGNT123) = RGNT_ARG(1:LRGNT123)
      NRGNT(1:LRGNT12) = NRGNT_ARG(1:LRGNT12)
      IRGNT(1:LRGNT123) = IRGNT_ARG(1:LRGNT123)
      N_INTEGER = N_INTEGER + LRGNT12 + LRGNT123
      N_REAL = N_REAL + LRGNT123
C
      ALLOCATE (IRREL(2,2,NKMMAX),SRREL(2,2,NKMMAX),NRREL(2,NKMMAX))
      IRREL(1:2,1:2,1:NKMMAX) = IRREL_ARG(1:2,1:2,1:NKMMAX)
      SRREL(1:2,1:2,1:NKMMAX) = SRREL_ARG(1:2,1:2,1:NKMMAX)
      NRREL(1:2,1:NKMMAX) = NRREL_ARG(1:2,1:NKMMAX)
      N_INTEGER = N_INTEGER + 2*2*NKMMAX + 2*2*NKMMAX + 2*NKMMAX
C
      ALLOCATE (HP(NLLMMMAX))
      HP(1:NLLMMMAX) = HP_ARG(1:NLLMMMAX)
      N_REAL = N_REAL + NLLMMMAX
C
      G123MAX = G123MAX_ARG
      R123MAX = R123MAX_ARG
C
c     NIJQMAX = MAXVAL(NIJQ_ARG, NQQP_STR)
c     ALLOCATE (NIJQ(NQQP_STR),IJQ(NQQP_STR,NQQP_STR))
c     ALLOCATE (NIJQ(NQQP_STR),IJQ(NIJQMAX,NQQP_STR))
c     NIJQ(1:NQQP_STR) = NIJQ_ARG(1:NQQP_STR)
c     IJQ(1:NQQP_STR,1:NQQP_STR) = IJQ_ARG(1:NQQP_STR,1:NQQP_STR)
c     IJQ(1:NIJQMAX,1:NQQP_STR) = IJQ_ARG(1:NIJQMAX,1:NQQP_STR)
c     N_INTEGER = N_INTEGER + NQQP_STR + NQQP_STR*NQQP_STR
c     N_INTEGER = N_INTEGER + NQQP_STR + NQQP_STR*NIJQMAX
C
c     ALLOCATE (SMAX(NQQP_STR))
      ALLOCATE (QQPX(NQQP_STR),QQPY(NQQP_STR),QQPZ(NQQP_STR))
      QQPX(1:NQQP_STR) = QQPX_ARG(1:NQQP_STR)
      QQPY(1:NQQP_STR) = QQPY_ARG(1:NQQP_STR)
      QQPZ(1:NQQP_STR) = QQPZ_ARG(1:NQQP_STR)
c     SMAX(1:NQQP_STR) = SMAX_ARG(1:NQQP_STR)
c     N_INTEGER = N_INTEGER + 4*NQQP_STR
      N_INTEGER = N_INTEGER + 3*NQQP_STR
C
      ALLOCATE (R1(NRDL),R2(NRDL),R3(NRDL))
      R1(1:NRDL) = R1_ARG(1:NRDL)
      R2(1:NRDL) = R2_ARG(1:NRDL)
      R3(1:NRDL) = R3_ARG(1:NRDL)
      N_REAL = N_REAL + 3*NRDL
C
c     ALLOCATE (G1(NGRL),G2(NGRL),G3(NGRL),EXPGNQ(NGRL,NQQP_STR))
      ALLOCATE (G1(NGRL),G2(NGRL),G3(NGRL))
      G1(1:NGRL) = G1_ARG(1:NGRL)
      G2(1:NGRL) = G2_ARG(1:NGRL)
      G3(1:NGRL) = G3_ARG(1:NGRL)
c     EXPGNQ(1:NGRL,1:NQQP_STR) = EXPGNQ_ARG(1:NGRL,1:NQQP_STR)
c     N_REAL = N_REAL + 3*NGRL + NGRL*NQQP_STR
      N_REAL = N_REAL + 3*NGRL 
C
c     ALLOCATE (INDR(NSDL,NQQP_STR))
c     INDR(1:NSDL,1:NQQP_STR) = INDR_ARG(1:NSDL,1:NQQP_STR)
c     N_INTEGER = N_INTEGER + NSDL*NQQP_STR
C
C--------------------------------------- restrict site pairs to JQ <= IQ
C------------------- for QQMLRS, GGJLRS and IILERS dealt with in part CC
      IF ( USE_NEW_BBDD_VERSION ) THEN
         NQQP_STR_CC = NQQP_STR_RED
      ELSE
         NQQP_STR_CC = NQQP_STR
      END IF
C
c     ALLOCATE (QQMLRS(NLLMMMAX,NSDL,NQQP_STR_CC),STAT=IA_ERR)
c     IF ( IA_ERR.NE.0 ) STOP 'alloc:INIT_MOD_STR->QQMLRS'
c     QQMLRS(1:NLLMMMAX,1:NSDL,1:NQQP_STR_CC)
c    &   = QQMLRS_ARG(1:NLLMMMAX,1:NSDL,1:NQQP_STR_CC)
c     N_REAL = N_REAL + 2*NLLMMMAX*NSDL*NQQP_STR_CC
C
c     ALLOCATE (GGJLRS(-J22MAX:LLARR,NSDL,NQQP_STR_CC))
C
c     GGJLRS(-J22MAX:LLARR,1:NSDL,1:NQQP_STR_CC)
c    &   = GGJLRS_ARG(-J22MAX:LLARR,1:NSDL,1:NQQP_STR_CC)
c     N_REAL = N_REAL + 2*(1+J22MAX)*(1+LLARR)*NSDL*NQQP_STR_CC
C
C  =====================================================================
C            allocate and initialize additional arrays
C  =====================================================================
C
      ALLOCATE (D1TERM3(0:LLARR),IILERS(0:LLARR,NSDL,NQQP_STR_CC))
      ALLOCATE (DLLMMKE(NLLMMMAX,NQQP_STR),CHP(0:LLARR),SHP(0:LLARR))
      ALLOCATE (WK(NKMMAX,NKMMAX),CILMAT(NLMAX**2,NLMAX**2),STAT=IA_ERR)
      N_REAL = N_REAL + 2*((1+LLARR)+(1+LLARR)*NSDL*NQQP_STR_CC)
      N_REAL = N_REAL + 2*(NLLMMMAX*NQQP_STR+(1+LLARR)+(1+LLARR))
      N_REAL = N_REAL + 2*(NKMMAX*NKMMAX+(NLMAX**2)*(NLMAX**2))
      IF ( IA_ERR.NE.0 ) STOP 'alloc:INIT_MOD_STR->CILMAT'
C
      DO LM1 = 1,NLMAX**2
         DO LM2 = 1,NLMAX**2
            CILMAT(LM1,LM2) = CIPWL_ARG(LM1)/CIPWL_ARG(LM2)
         END DO
      END DO
C
      ALLOCATE (M1PWL(NLLMMMAX))
C
      MMLL = 0
      M1PW = -1D0
      DO LL = 0,LLMAX
         M1PW = -M1PW
         DO MM = -LL, + LL
            MMLL = MMLL + 1
            M1PWL(MMLL) = M1PW
         END DO
      END DO
C
      N = G123MAX
      ALLOCATE (PWEX2K1(-N:N),PWEX2K2(-N:N),PWEX2K3(-N:N))
      N = R123MAX
      ALLOCATE (PWEXIK1(-N:N),PWEXIK2(-N:N),PWEXIK3(-N:N),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:STRINIT->PWEXIK2'
      ALLOCATE (PWP(NLLMMMAX))
      N_REAL = N_REAL + 2*6*(2*N+1)
C
      WRITE (6,99001) N_INTEGER,N_REAL,NINT(N_INTEGER*2/(1024D0**2)),
     &                N_REAL*8/(1024**2)
C
99001 FORMAT (/,10X,'<INIT_MOD_STR>: section 1 of module MOD_STR',
     &        ' initialized',/,10X,'allocated  ',I10,'  integers',I14,
     &        '  reals (incl. complex)',/,10X,'corresp. to',I10,
     &        '  Mbytes  ',I14,'  Mbytes',/)
      END
