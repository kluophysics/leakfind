C*==strinit.f    processed by SPAG 6.70Rc at 07:55 on 26 Apr 2017
      SUBROUTINE STRINIT(ETOP)
C   ********************************************************************
C   *                                                                  *
C   *               KKR - structure constant routines                  *
C   *                                                                  *
C   *                       ===============                            *
C   *                          <STRINIT>                               *
C   *                       ===============                            *
C   *                                                                  *
C   *  call once to initialize all                                     *
C   *  <STRINIT> calls <VECGEN>, <GAUNT>, <INITS>, <STRAA>             *
C   *  supply:                                                         *
C   *  ARGUMENTS       ETA, RMAX, GMAX, QBAS                           *
C   *                                                                  *
C   *              NT,NQHOST, NL,NLM,NLMQ, NK,NKM,NLIN, NKKR           *
C   *              ^   ^                                               *
C   *                  ALAT                                            *
C   *                  ^                                               *
C   *  <STRAA>  calculates all energy- and k-indep. terms of D[L,M]    *
C   *  <INITS>  creates coefficients to transform real non-relat.      *
C   *           G[L,L'] to relat. G[Lam,Lam']                          *
C   *                                                                  *
C   *                                                                  *
C   *                       ===============                            *
C   *                           <STRCC>                                *
C   *                       ===============                            *
C   *                                                                  *
C   *  call once per energy-value                                      *
C   *  <STRC> calls NO other routines                                  *
C   *  supply:                                                         *
C   *                  ERYD,P                                          *
C   *                  ^                                               *
C   *  ERYD is converted to  EDU  used in <STRCC> and <STRBBDD>        *
C   *                                                                  *
C   *                       ===============                            *
C   *                          <STRSET>                                *
C   *                       ===============                            *
C   *  supply:                                                         *
C   *  arguments: DLLMMKE, KX,KY,KZ                                    *
C   *  KX,KY,KZ   in multiples of 2*pi/a  i.e.  in  d.u.'s             *
C   *                                                                  *
C   *  call once per k-point                                           *
C   *  <STRBBDD> is called to get the  D[L,M]'s                        *
C   *  the G[L,L']'s are then set up in a.u.                           *
C   *  conversion is managed by multiplying the gaunts with  2*pi/a    *
C   *  the factor  4*pi  is also included in the gaunts                *
C   *                                                                  *
C   *                       ===============                            *
C   *                          <STRBBDD>                               *
C   *                       ===============                            *
C   *                                                                  *
C   *  call once per k-point                                           *
C   *  supply:                                                         *
C   *  arguments: DLLMMKE, KX,KY,KZ                                    *
C   *  KX,KY,KZ   in multiples of 2*pi/a  i.e.  in  d.u.'s             *
C   *  D[L,M]     is created in  d.u.'s                                *
C   *                                                                  *
C   *  NOTE:   all str-routines work internally in d.u.'s              *
C   *                                                                  *
C   ********************************************************************
C   *                                                                  *
C   *  INITIALIZE CALCULATION OF STRUCTURE CONSTANTS                   *
C   *                                                                  *
C   ********************************************************************
C   *                                                                  *
C   *   ETA    :EWALD PARAMETER                                        *
C   *   RMAX   :RADIUS OF CONVERGENCE SPHERE IN REAL SPACE             *
C   *   GMAX   :RADIUS OF CONVERGENCE SPHERE IN RECIPROCAL SPACE       *
C   *   ABAS   :primitive vectors in real space  (UNITS OF A)          *
C   *   BBAS   :primitive vectors in reciprocal space  (in 2*PI/A)     *
C   *   NL     :NUMBER OF L-QUANTUM NUMBERS,LMAX=NL-1                  *
C   *   NQHOST :NUMBER OF ATOMS IN THE PRIMITIVE CELL                  *
C   *   QBAS   :basis vectors in real space      (UNITS OF A)          *
C   *                                                                  *
C   ********************************************************************
C   *                                                                  *
C   *   NPT    :NUMBER OF K-POINTS IN IRREDUCIBLE ZONE                 *
C   *   KX                                                             *
C   *   KY     :K-VECTOR     (UNITS OF 2*PI/A)   RECTANGULAR COORD.    *
C   *   KZ                                                             *
C   *                                                                  *
C   ********************************************************************
C   *                                                                  *
C   * NOTE: to avoid conflicts the variable names in the argument list *
C   *       have 'P' added at the beginning                            *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_ANGMOM,ONLY:NL,NLMAX,NKMMAX
      USE MOD_FILES,ONLY:IPRINT,FOUND_SECTION
      USE MOD_LATTICE,ONLY:ABAS,BBAS,ALAT,VOLUC
      USE MOD_SITES,ONLY:QBAS,NQHOST,NQMAX
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_CONSTANTS,ONLY:PI
      USE MOD_STR,ONLY:J13MAX_UPPER_LIMIT,J22MAX_UPPER_LIMIT
      !KL:  let strvecgen.f do the init to IJQ and NIJQ
      !USE MOD_STR,ONLY: IJQ, NIJQ
      USE MOD_STR,ONLY: USE_NEW_BBDD_VERSION
      USE MOD_STR,ONLY: GGJLRS, QQMLRS, INDR, SMAX, EXPGNQ
      USE MOD_STR,ONLY: NQQP_STR_RED, NQQP_STR_CC
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='STRINIT')
      REAL*8 SMALL,ALIM
      PARAMETER (SMALL=1D-15,ALIM=225D0)
C
C Dummy arguments
C
      REAL*8 ETOP
C
C Local variables
C
c     REAL*8 ALPHA0,B,ETA,ETA_INI,GA,GGJLRS(:,:,:),GMAX,GMAXSQ,GMAX_INI,
      REAL*8 ALPHA0,B,ETA,ETA_INI,GA,GMAX,GMAXSQ,GMAX_INI,
     &       HP(:),PRETA,QQPX(:),QQPY(:),QQPZ(:),RA,RGNT(:),RMAX,
     &       RMAX_INI
c     COMPLEX*16 CIPWL(:),EXPGNQ(:,:),QQMLRS(:,:,:),SRREL(:,:,:)
      COMPLEX*16 CIPWL(:),SRREL(:,:,:)
c     INTEGER G1(:),G123MAX,G2(:),G3(:),I,IA_ERR,IG123,IJQ(:,:),
      INTEGER G1(:),G123MAX,G2(:),G3(:),I,IA_ERR,IG123,
c    &        INDR(:,:),IQ,IRGNT(:),IRREL(:,:,:),J,J13MAX,J22MAX,LLARR,
     &        IQ,IRGNT(:),IRREL(:,:,:),J,J13MAX,J22MAX,LLARR,
c    &        LLMAX,LMAX,LRGNT12,LRGNT123,MMLLMAX,NGRL,NGRLMAX,NIJQ(:),
     &        LLMAX,LMAX,LRGNT12,LRGNT123,MMLLMAX,NGRL,NGRLMAX,
     &        NLLMMMAX,NQQP_STR,NQQP_STRMAX,NQ_STR,NRDL,NRDLMAX,NRGNT(:)
     &        ,NRGNT123TAB(20),NRREL(:,:),NSDL,NUMGH,NUMRH,R1(:),
c    &        R123MAX,R2(:),R3(:),SMAX(:),SMAXMAX,SMAXMIN
     &        R123MAX,R2(:),R3(:),SMAXMAX,SMAXMIN
C
C*** End of declarations rewritten by SPAG
C
      DATA NRGNT123TAB/1,15,96,388,1181,2917,6342,12452,22525,38289,
     &     61912,95914,143531,208371,294744,407644,552931,736829,966544,
     &     1250346/,IA_ERR/0/
C
c     ALLOCATABLE RGNT,IRGNT,NRGNT,CIPWL,IJQ,NIJQ,GGJLRS,EXPGNQ
      ALLOCATABLE RGNT,IRGNT,NRGNT,CIPWL
c     ALLOCATABLE SRREL,NRREL,IRREL,QQPX,QQPY,QQPZ,QQMLRS,INDR,SMAX
      ALLOCATABLE SRREL,NRREL,IRREL,QQPX,QQPY,QQPZ
      ALLOCATABLE G1,G2,G3,R1,R2,R3,HP
C
C-----------------------------------------------------------------------
C    get convergence parameter ETA, RMAX, GMAX from input if set
C-----------------------------------------------------------------------
C
      ETA_INI = 0D0
      RMAX_INI = 0D0
      GMAX_INI = 0D0
C
      CALL INPUT_FIND_SECTION('STRCONST',0)
C
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_SET_REAL('ETA',ETA_INI,9999D0,0)
         CALL SECTION_SET_REAL('RMAX',RMAX_INI,9999D0,0)
         CALL SECTION_SET_REAL('GMAX',GMAX_INI,9999D0,0)
      END IF
C
C-----------------------------------------------------------------------
C
      NQ_STR = NQHOST
C
      LLARR = 2*(NLMAX-1)
      NLLMMMAX = (LLARR+1)**2
C
C   ********************************************************************
C
      WRITE (6,99002)
      ETA = ETA_INI
      RMAX = RMAX_INI
      GMAX = GMAX_INI
C
C -------------------------------- set ETA according to AKAI's algorithm
C
      IF ( ETA.LT.1D-3 ) THEN
         B = -LOG(SMALL)
         PRETA = 1D0/VOLUC**(2D0/3D0)/PI
         PRETA = PRETA*0.75D0
         ETA = MAX(PRETA,ABS(ETOP)/ALIM)
         GMAX = SQRT(B*ETA)
         RMAX = SQRT(B/ETA)/PI
         PRETA = 5D-1*SQRT(PRETA)
         ETA = 5D-1*SQRT(ETA)
      END IF
C
      WRITE (6,99004)
      DO J = 1,3
         WRITE (6,99003) (ABAS(I,J),I=1,3)
      END DO
      WRITE (6,'(/,10X,''basis vectors in units of a''/)')
      DO IQ = 1,NQ_STR
         WRITE (6,99003) (QBAS(I,IQ),I=1,3)
      END DO
C
C------------------- primitive vectors     BBAS      of reciprocal space
C
      WRITE (6,99001)
      DO J = 1,3
         WRITE (6,99003) (BBAS(I,J),I=1,3)
      END DO
C
C -------------------------- find set of inequivelent lattice pairs Q-QP
C ------------------------------- and range for R- and K-lattice vectors
C     NQQP_STR: inequivalent combination of lattice sites
C               in structure constants matrix
C
      NQQP_STRMAX = NQ_STR*(NQ_STR-1) + 1
C
c     ALLOCATE (IJQ(NQQP_STRMAX,NQQP_STRMAX),NIJQ(NQQP_STRMAX))
      ALLOCATE (QQPX(NQQP_STRMAX))
      ALLOCATE (QQPY(NQQP_STRMAX),QQPZ(NQQP_STRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: QQPZ')
C
c     NIJQ(:) = 0
c     IJQ(:,:) = 0
C
c     CALL STRQQPLIM(NQ_STR,IPRINT,RMAX,GMAX,NQQP_STR,NIJQ,IJQ,NUMRH,
      CALL STRQQPLIM(NQ_STR,IPRINT,RMAX,GMAX,NQQP_STR,NUMRH,
     &               NUMGH,RA,GA,ABAS,BBAS,QBAS,QQPX,QQPY,QQPZ,NQMAX,
     &               NQQP_STRMAX)
C
C ---------------------- generate the  SMAX  and  NGRL  shortest vectors
C ----------------------------------------- of real and reciprocal space
C
      NRDLMAX = (2*NUMRH+1)**3
      ALLOCATE (R1(NRDLMAX),R2(NRDLMAX),R3(NRDLMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: R1')
      NGRLMAX = (2*NUMGH+1)**3
      ALLOCATE (G1(NGRLMAX),G2(NGRLMAX),G3(NGRLMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: G1')
      ALLOCATE (SMAX(NQQP_STRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: SMAX')
C
      DO I = 1,NRDLMAX
         R1(I) = 0
         R2(I) = 0
         R3(I) = 0
      END DO
      DO I = 1,NGRLMAX
         G1(I) = 0
         G2(I) = 0
         G3(I) = 0
      END DO
C
      CALL STRVECGEN(IPRINT,ALAT,RMAX,GMAX,GMAXSQ,NRDL,SMAX,NGRL,
     &               NQQP_STR,NUMGH,NUMRH,GA,RA,ABAS,BBAS,QQPX,QQPY,
     &               QQPZ,G1,G2,G3,R1,R2,R3,NGRLMAX,NRDLMAX,NQQP_STRMAX)
C
C ----------------------------------------- calculate Gaunt coefficients
C
      LRGNT123 = NRGNT123TAB(NL)
      LRGNT12 = (NL**2*(NL**2+1)/2)
C
      ALLOCATE (IRGNT(LRGNT123),NRGNT(LRGNT12),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: IRGNT')
      ALLOCATE (RGNT(LRGNT123),CIPWL((2*NLMAX)**2),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: RGNT')
C
      LMAX = NL - 1
      LLMAX = 2*LMAX
      MMLLMAX = (LLMAX+1)**2
C
      CALL STRGAUNT(LMAX,ALAT,RGNT,NRGNT,IRGNT,CIPWL,NLMAX,IG123,
     &              LRGNT12,LRGNT123)
C
C ------------------------------------ calculate transformation matrices
C
      ALLOCATE (IRREL(2,2,NKMMAX),NRREL(2,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: IRREL')
      ALLOCATE (SRREL(2,2,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: SRREL')
      CALL CINIT(2*2*NKMMAX,SRREL)
      IRREL(:,:,:) = 0
      NRREL(:,:) = 0
C
      IF ( IREL.GE.3 ) CALL STRSMAT(SRREL,NRREL,IRREL)
C
C -------------------------------------------- calculate STR-CONST terms
C
      SMAXMIN = 100000
      SMAXMAX = 0
      DO I = 1,NQQP_STR
         SMAXMIN = MIN(SMAX(I),SMAXMIN)
         SMAXMAX = MAX(SMAX(I),SMAXMAX)
      END DO
      NSDL = SMAXMAX
C
c     KL: use RED version
      IF ( USE_NEW_BBDD_VERSION ) THEN
         NQQP_STR_CC = NQQP_STR_RED
      ELSE
         NQQP_STR_CC = NQQP_STR
      END IF

c     ALLOCATE (GGJLRS(-J22MAX_UPPER_LIMIT:LLARR,NSDL,NQQP_STR),
c    &          STAT=IA_ERR)
c     IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: GGJLRS')
c     ALLOCATE (EXPGNQ(NGRL,NQQP_STR),STAT=IA_ERR)
c     IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: EXPGNQ')
c     ALLOCATE (INDR(NSDL,NQQP_STR),STAT=IA_ERR)
c     IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: INDR')
c     ALLOCATE (QQMLRS(NLLMMMAX,NSDL,NQQP_STR),STAT=IA_ERR)
c     IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: QQMLRS')
c     ALLOCATE (HP(NLLMMMAX),STAT=IA_ERR)
c     IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: HP')
c     CALL RINIT(NLLMMMAX,HP)

      ALLOCATE (GGJLRS(-J22MAX_UPPER_LIMIT:LLARR,NSDL,NQQP_STR_CC),
     &          STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: GGJLRS')
      ALLOCATE (EXPGNQ(NGRL,NQQP_STR_CC),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: EXPGNQ')
      ALLOCATE (INDR(NSDL,NQQP_STR_CC),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: INDR')
      ALLOCATE (QQMLRS(NLLMMMAX,NSDL,NQQP_STR_CC),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: QQMLRS')
      ALLOCATE (HP(NLLMMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: HP')
      CALL RINIT(NLLMMMAX,HP)

C
c     CALL STRAA(IPRINT,RMAX,NRDL,SMAX,NGRL,NQQP_STR,LLMAX,ETA,ALPHA0,
      CALL STRAA(IPRINT,RMAX,NRDL,SMAX,NGRL,NQQP_STR_CC,LLMAX,ETA,
     &         ALPHA0,ABAS,BBAS,QQPX,QQPY,QQPZ,EXPGNQ,QQMLRS,G1,G2,G3,
     &           G123MAX,R1,R2,R3,R123MAX,NSDL,INDR,HP,GGJLRS,NGRLMAX,
c    &           NRDLMAX,J13MAX,J22MAX,NQQP_STRMAX,LLARR,NLLMMMAX)
     &           NRDLMAX,J13MAX,J22MAX,NQQP_STR_CC,LLARR,NLLMMMAX)
C
      WRITE (6,99005) NGRL,NGRLMAX,NRDL,NRDLMAX,SMAXMAX,NSDL,SMAXMIN,NL,
     &                NLMAX,NQ_STR,NQMAX,NQQP_STR,NQQP_STRMAX,R123MAX,
     &                G123MAX,J13MAX,J13MAX_UPPER_LIMIT,J22MAX,
     &                J22MAX_UPPER_LIMIT,LRGNT12,LRGNT123,ETA
C
C=======================================================================
C        pass primary structure constants data to module  MOD_STR
C=======================================================================
C
      NRDL = NRDLMAX
C
      CALL INIT_MOD_STR(NGRLMAX,NRDLMAX,NQQP_STRMAX,LRGNT123,LRGNT12,
     &                  LLARR,NLLMMMAX,LLMAX,MMLLMAX,J13MAX,J22MAX,
     &                  J22MAX_UPPER_LIMIT,NQQP_STR,NGRL,NSDL,NRDL,
     &                  ALPHA0,GMAXSQ,ETA,CIPWL,RGNT,NRGNT,IRGNT,SRREL,
c    &                  NRREL,IRREL,HP,G123MAX,R123MAX,NIJQ,IJQ,QQPX,
     &                  NRREL,IRREL,HP,G123MAX,R123MAX,QQPX,
c    &                  QQPY,QQPZ,SMAX,R1,R2,R3,G1,G2,G3,EXPGNQ,INDR,
     &                  QQPY,QQPZ,R1,R2,R3,G1,G2,G3)
c    &                  QQMLRS,GGJLRS)
C
c     DEALLOCATE (RGNT,IRGNT,NRGNT,CIPWL,IJQ,NIJQ,GGJLRS,EXPGNQ)
c     DEALLOCATE (RGNT,IRGNT,NRGNT,CIPWL,GGJLRS,EXPGNQ)
      DEALLOCATE (RGNT,IRGNT,NRGNT,CIPWL)
c     DEALLOCATE (SRREL,NRREL,IRREL,QQPX,QQPY,QQPZ,QQMLRS,INDR,SMAX)
      DEALLOCATE (SRREL,NRREL,IRREL,QQPX,QQPY,QQPZ)
      DEALLOCATE (G1,G2,G3,R1,R2,R3,HP)
C
      RETURN
C=======================================================================
99001 FORMAT (/,10X,'primitive vectors of reciprocal space',/,10X,
     &        '       in units of 2*pi/a',/)
99002 FORMAT (/,1X,79('*'),/,35X,'<STRINIT>',/,1X,79('*'),//,10X,
     &        'parameters for calculation of structure constants:')
99003 FORMAT (12X,'(',F10.5,',',F10.5,',',F10.5,' )')
99004 FORMAT (/,10X,'primitive vectors for Bravais lattice',/)
99005 FORMAT (/,1X,79('*'),/,12X,
     &        'array sizes for calculation of structure constants',/,1X,
     &        79('*'),/,10X,32X,'used',14X,'available',/,10X,
     &        'G-vectors table       NGRL   ',I7,8X,'NGRLMAX ',I7,/,10X,
     &        'R-vectors table       NRDL   ',I7,8X,'NRDLMAX ',I7,/,10X,
     &        'R-vectors (max.)      SMAX   ',I7,8X,'NSDL    ',I7,/,10X,
     &        'R-vectors (min.)      SMAX   ',I7,/,10X,
     &        'l-max + 1             NL     ',I7,8X,'NLMAX   ',I7,/,10X,
     &        'sites                 NQ_STR ',I7,8X,'NQMAX   ',I7,/,10X,
     &        'q-q''-blocks           NQQP_STR',I6,8X,'LIMIT   ',I7,/,
     &        10X,'auxilary array        R123MAX',I7,/,10X,
     &        'auxilary array        G123MAX',I7,/,10X,
     &        'Davis eq. (13)        J13MAX ',I7,8X,'LIMIT   ',I7,/,10X,
     &        'Davis eq. (22)        J22MAX ',I7,8X,'LIMIT   ',I7,/,10X,
     &        'Gaunts                LRGNT12',I7,/,10X,
     &        'Gaunts                LRGNT123',I6,/,10X,
     &        'Ewald parameter       ETA    ',F7.2,/,1X,79('*'),/)
C    &     10X,'lattice sum  evaluated by continued fractions',/,
      END
