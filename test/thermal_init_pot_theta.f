C*==thermal_init_pot_theta.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE THERMAL_INIT_POT_THETA
C   ********************************************************************
C   *                                                                  *
C   *   write out the starting SCF-potential                           *
C   *                                                                  *
C   *   extend the list of types according to   NTET_FLUCT             *
C   *                                                                  *
C   *   for each atom type   NTET_FLUCT  orientations                  *
C   *   of the magnetic moment are allowed                             *
C   *                                                                  *
C   ********************************************************************
      USE MOD_SCF,ONLY:ITRSCF,SCFVXC,SCFTOL,SCFMIX,SCFALG,RMSAVV,RMSAVB,
     &    SCFSTATUS,SCF_THETA_DEPENDENT_POT
      USE MOD_KSPACE,ONLY:IBZINT,NKTAB
      USE MOD_ENERGY,ONLY:EFERMI,NETAB,EWORK
      USE MOD_RMESH,ONLY:DX,NPAN,RMESHTYPE,BRMSH,ARMSH,JRCUT,JRCRI,
     &    JRNS1,SPHERCELL,FULLPOT,RWS,RMT,R,NM,JRWS,JRMT,BASSCALE,
     &    CARTESIAN
      USE MOD_FILES,ONLY:IOTMP,TITLE,SYSTEM,POTFIL,LTITLE,LSYSTEM,
     &    LPOTFIL,LINFO,INFO,IDUMMY,POTFMTOUT,PACKAGE
      USE MOD_LATTICE,ONLY:SYSTEM_DIMENSION,SYSTEM_TYPE,ABAS,ABAS_L,
     &    ABAS_R,BRAVAIS,ALAT
      USE MOD_CALCMODE,ONLY:PROGNAME,NONMAG,SEMICORE,ORBPOL,KMROT,IREL,
     &    EXTFIELD,LLOYD,BREITINT,BLCOUPL
      USE MOD_LATTICE,ONLY:TXTBRAVAIS
      USE MOD_SITES,ONLY:NQ_R,NQ_L,QBAS,QMVEC,QMTET,QMPHI,NQ,NOQ,ITOQ,
     &    VLMMAD_HOST,CMNTQ,NLMVMAD,NLMQMAD
      USE MOD_TYPES,ONLY:QEL,ABIT,AOPT,KLMFP,BNST,VNST,NFPT,NLMFPT,VMTZ,
     &    Z,VT,TXT_T,NVALT,NCORT,RHO2NS,RHOSPN,RHOCHR,LOPT,NT,
     &    NSEMCORSHLT,IMT,CONC,BT,BEXT,OBS_T
      USE MOD_ANGMOM,ONLY:NLABIMAX,NSPIN,IDOS,ISMT,IOMT,IHFF
      USE MOD_TB,ONLY:NREF,RMTREF,VREF,IREFQ
      USE MOD_THERMAL,ONLY:NTET_FLUCT
      IMPLICIT NONE
C*--THERMAL_INIT_POT_THETA36
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='THERMAL_INIT_POT_THETA')
C
C Local variables
C
      REAL*8 DTET,MGAM,MPHI,MTET,X_IT_TET(:)
      CHARACTER*40 FMT02,FMT03,FMT04,FMT05,FMT06,FMT07,FMT08,FMT12
      INTEGER IC,ILA,IM,IO,IO_TET,IPAN,IQ,IR,IRSHFT,IRTOP,IT,ITET,
     &        IT_TET,IT_TET_OQ(:,:),J,LM,MS,N,NT_TET
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE X_IT_TET,IT_TET_OQ
C
C*** End of declarations rewritten by SPAG
C
C-----------------------------------------------------------------------
C                    set appropriate output format
C-----------------------------------------------------------------------
C
      IF ( POTFMTOUT.NE.7 ) POTFMTOUT = 8
C
      IF ( POTFMTOUT.GE.7 ) THEN
C >> CONC
         FMT02 = '(I10,3I10,10(I6,F16.12))'
C >> QBAS,VREF,RMTREF,VLMMAD_HOST,QMTET,QMPHI
         FMT03 = '(I10,1P,3E22.14)'
C >> VLMMAD_HOST
         FMT04 = '(6X,2I4,1P,3E22.14)'
C >> R1,DX,RMT,RWS
         FMT05 = '(1X,I4,1P,2E22.14,2(I5,1P,E22.14))'
C >> NVALT
         FMT06 = '(1X,I4,5X,A,I10,I10,I10,9X,I4)'
C >> V,B,...
         FMT07 = '(1P,5E22.14)'
C >> QMVEC
         FMT08 = '(A10,3(1P,E22.14))'
C >> IM,JRNS1,JRCRI,NPAN,JRCUT
         FMT12 = '(I5,3I6,20I5)'
      END IF
C
C ======================================================================
      KMROT = 2
C ======================================================================
      NT_TET = NT*NTET_FLUCT
      N = NTET_FLUCT
C
      ALLOCATE (X_IT_TET(NT_TET))
      ALLOCATE (IT_TET_OQ(NT_TET,NQ))
C
      DO IT_TET = 1,NT_TET
         IT = (IT_TET-1)/NTET_FLUCT + 1
         X_IT_TET(IT_TET) = CONC(IT)/DBLE(NTET_FLUCT)
      END DO
C
      DO IQ = 1,NQ
         IO_TET = 0
         DO IO = 1,NOQ(IQ)
            IT = ITOQ(IO,IQ)
            IT_TET = (IT-1)*NTET_FLUCT
            DO ITET = 1,NTET_FLUCT
               IT_TET = IT_TET + 1
               IO_TET = IO_TET + 1
               IT_TET_OQ(IO_TET,IQ) = IT_TET
            END DO
         END DO
      END DO
C ======================================================================
C
C-----------------------------------------------------------------------
C
      CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,POTFIL(1:LPOTFIL-4)
     &                     //'_THETA.pot')
C
C ======================================================================
      WRITE (IOTMP,99020) 'HEADER    ','SPR-KKR dataset created by ',
     &                    PROGNAME
      WRITE (IOTMP,99019)
      WRITE (IOTMP,99015) 'TITLE     ',TITLE(1:LTITLE)
      WRITE (IOTMP,99014) 'SYSTEM    ',SYSTEM(1:LSYSTEM)
      WRITE (IOTMP,99014) 'PACKAGE   ',PACKAGE
      WRITE (IOTMP,'(A,I2,A)') 'FORMAT    ',POTFMTOUT,' (28.03.2008)'
C ======================================================================
      WRITE (IOTMP,99001) 'GLOBAL SYSTEM PARAMETER'
      WRITE (IOTMP,99017) 'NQ        ',NQ
      WRITE (IOTMP,99017) 'NT        ',NT_TET
      WRITE (IOTMP,99017) 'NM        ',NM
      WRITE (IOTMP,99017) 'IREL      ',IREL
      WRITE (IOTMP,99017) 'NSPIN     ',NSPIN
C ======================================================================
      WRITE (IOTMP,99001) 'SCF-INFO  '
      WRITE (IOTMP,99015) 'INFO      ',INFO(1:LINFO)
      WRITE (IOTMP,99015) 'SCFSTATUS ',SCFSTATUS
      WRITE (IOTMP,99016) 'FULLPOT   ',FULLPOT
      IF ( FULLPOT ) WRITE (IOTMP,99016) 'SPHERCELL ',SPHERCELL
C
      WRITE (IOTMP,99016) 'BREITINT  ',BREITINT
      WRITE (IOTMP,99016) 'NONMAG    ',NONMAG
      WRITE (IOTMP,99014) 'ORBPOL    ',ORBPOL
      WRITE (IOTMP,99016) 'EXTFIELD  ',EXTFIELD
      WRITE (IOTMP,99016) 'BLCOUPL   ',BLCOUPL
      WRITE (IOTMP,FMT=FMT08) 'BEXT      ',BEXT
      WRITE (IOTMP,99016) 'SEMICORE  ',SEMICORE
      WRITE (IOTMP,99016) 'LLOYD     ',LLOYD
      WRITE (IOTMP,99017) 'NE        ',NETAB
      WRITE (IOTMP,99017) 'IBZINT    ',IBZINT
      WRITE (IOTMP,99017) 'NKTAB     ',NKTAB
      IF ( POTFMTOUT.GE.8 ) WRITE (IOTMP,99016) 'TETDEPPOT ',
     &                             SCF_THETA_DEPENDENT_POT
      WRITE (IOTMP,99014) 'XC-POT    ',SCFVXC
      WRITE (IOTMP,99014) 'SCF-ALG   ',SCFALG
      WRITE (IOTMP,99017) 'SCF-ITER  ',ITRSCF
      WRITE (IOTMP,FMT=FMT08) 'SCF-MIX   ',SCFMIX
      WRITE (IOTMP,FMT=FMT08) 'SCF-TOL   ',SCFTOL
      WRITE (IOTMP,FMT=FMT08) 'RMSAVV    ',RMSAVV
      WRITE (IOTMP,FMT=FMT08) 'RMSAVB    ',RMSAVB
      WRITE (IOTMP,FMT=FMT08) 'EF        ',EFERMI
      WRITE (IOTMP,FMT=FMT08) 'VMTZ      ',VMTZ
      IF ( SYSTEM_TYPE(1:3).EQ.'LIV' ) WRITE (IOTMP,FMT=FMT08)
     &      'EWORK     ',EWORK
      WRITE (IOTMP,99001) 'LATTICE   '
      WRITE (IOTMP,99015) 'SYSDIM    ',SYSTEM_DIMENSION
      WRITE (IOTMP,99015) 'SYSTYPE   ',SYSTEM_TYPE
      WRITE (IOTMP,99018) 'BRAVAIS   ',BRAVAIS,TXTBRAVAIS(BRAVAIS)
      WRITE (IOTMP,FMT=FMT08) 'ALAT      ',ALAT
      DO J = 1,3
         WRITE (IOTMP,FMT=FMT08) 'A('//CHAR(ICHAR('1')+J-1)//')      ',
     &                           (ABAS(IC,J),IC=1,3)
      END DO
      IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' ) THEN
         WRITE (IOTMP,99017) 'NQ_L      ',NQ_L
         WRITE (IOTMP,FMT=FMT08) 'A_L(3)    ',ABAS_L(1:3,3)
         WRITE (IOTMP,99017) 'NQ_R      ',NQ_R
         WRITE (IOTMP,FMT=FMT08) 'A_R(3)    ',ABAS_R(1:3,3)
      END IF
C ======================================================================
      WRITE (IOTMP,99001) 'SITES     '
      WRITE (IOTMP,99016) 'CARTESIAN ',CARTESIAN
      WRITE (IOTMP,FMT=FMT08) 'BASSCALE  ',(BASSCALE(IO),IO=1,3)
      WRITE (IOTMP,99002)
      DO IQ = 1,NQ
         WRITE (IOTMP,FMT=FMT03) IQ,(QBAS(IC,IQ),IC=1,3)
      END DO
C ======================================================================
      WRITE (IOTMP,99001) 'OCCUPATION'
      WRITE (IOTMP,99003)
      DO IQ = 1,NQ
         DO IO = 1,NOQ(IQ)
            IDUMMY = IMT(ITOQ(IO,IQ))
         END DO
C
         WRITE (IOTMP,FMT=FMT02) IQ,IREFQ(IQ),IDUMMY,NOQ(IQ)*N,
     &                           (IT_TET_OQ(IO_TET,IQ),
     &                           X_IT_TET(IT_TET_OQ(IO_TET,IQ)),
     &                           IO_TET=1,NOQ(IQ)*N)
      END DO
C ======================================================================
      WRITE (IOTMP,99001) 'REFERENCE SYSTEM FOR TIGHT BINDING MODE'
      WRITE (IOTMP,99017) 'NREF      ',NREF
      WRITE (IOTMP,99009)
      DO IQ = 1,NREF
         WRITE (IOTMP,FMT=FMT03) IQ,VREF(IQ),RMTREF(IQ)
      END DO
C ======================================================================
      WRITE (IOTMP,99001) 'HOST MADELUNG POTENTIAL'
      WRITE (IOTMP,99010) 'VLMMAD'
      WRITE (IOTMP,99017) 'NLMTOP-POT',NLMVMAD
      DO IQ = 1,NQ
         DO LM = 1,NLMVMAD
            WRITE (IOTMP,FMT=FMT04) IQ,LM,VLMMAD_HOST(LM,IQ)
         END DO
      END DO
C ======================================================================
      WRITE (IOTMP,99001) 'CHARGE MOMENTS'
      WRITE (IOTMP,99010) 'CMNTQ '
      WRITE (IOTMP,99017) 'NLMTOP-CHR',NLMQMAD
      DO IQ = 1,NQ
         DO LM = 1,NLMVMAD
            WRITE (IOTMP,FMT=FMT04) IQ,LM,CMNTQ(LM,IQ)
         END DO
      END DO
C ======================================================================
      WRITE (IOTMP,99001) 'MAGNETISATION DIRECTION'
      WRITE (IOTMP,99017) 'KMROT     ',KMROT
      WRITE (IOTMP,FMT=FMT08) 'QMVEC     ',(QMVEC(IC),IC=1,3)
      WRITE (IOTMP,99004)
      DO IQ = 1,NQ
         WRITE (IOTMP,FMT=FMT03) IQ,QMTET(IQ),QMPHI(IQ)
      END DO
      WRITE (IOTMP,99005)
      DTET = 180.D0/(NTET_FLUCT-1)
      MPHI = 0D0
      MGAM = 0D0
      DO IT_TET = 1,NT_TET
         IT = (IT_TET-1)/NTET_FLUCT + 1
         ITET = IT_TET - (IT-1)*NTET_FLUCT
         MTET = (ITET-1)*DTET
         WRITE (IOTMP,FMT=FMT03) IT_TET,MTET,MPHI,MGAM
      END DO
C======================================================================
      WRITE (IOTMP,99001) 'MESH INFORMATION'
      WRITE (IOTMP,99014) 'MESH-TYPE ',RMESHTYPE
      IRSHFT = 0
      IF ( RMESHTYPE.EQ.'EXPONENTIAL ' ) THEN
         WRITE (IOTMP,99006)
         DO IM = 1,NM
            WRITE (IOTMP,FMT=FMT05) IM,R(1,IM),DX(IM),JRMT(IM),RMT(IM),
     &                              JRWS(IM),RWS(IM)
         END DO
      ELSE IF ( RMESHTYPE.EQ.'JUELICH     ' ) THEN
         IF ( PACKAGE(1:7).EQ.'SPR-KKR' ) IRSHFT = 1
         WRITE (IOTMP,99007)
         DO IM = 1,NM
            WRITE (IOTMP,FMT=FMT05) IM,ARMSH(IM),BRMSH(IM),JRMT(IM)
     &                              + IRSHFT,RMT(IM),JRWS(IM) + IRSHFT,
     &                              RWS(IM)
         END DO
      ELSE
         WRITE (6,99015) 'MESH-TYPE ',RMESHTYPE
         STOP 'in <POTWR>  mesh-type  not allowed'
      END IF
C ----------------------------------------------------------------------
      IF ( FULLPOT ) THEN
         WRITE (IOTMP,99008)
         DO IM = 1,NM
            WRITE (IOTMP,FMT=FMT12) IM,JRNS1(IM) + IRSHFT,JRCRI(IM)
     &                              + IRSHFT,NPAN(IM),
     &                              (JRCUT(IPAN,IM)+IRSHFT,IPAN=1,
     &                              NPAN(IM))
         END DO
      END IF
C ----------------------------------------------------------------------
C
C ======================================================================
      WRITE (IOTMP,99001) 'TYPES'
      WRITE (IOTMP,99011)
      DO IT_TET = 1,NT_TET
         IT = (IT_TET-1)/NTET_FLUCT + 1
         WRITE (IOTMP,FMT=FMT06) IT_TET,TXT_T(IT),Z(IT),NCORT(IT),
     &                           NVALT(IT),NSEMCORSHLT(IT)
      END DO
C=======================================================================
      WRITE (IOTMP,99001) 'POTENTIAL'
      DO IT_TET = 1,NT_TET
         IT = (IT_TET-1)/NTET_FLUCT + 1
         IM = IMT(IT)
         IF ( FULLPOT ) THEN
            IRTOP = JRCRI(IM)
         ELSE
            IRTOP = JRWS(IM)
         END IF
         WRITE (IOTMP,99017) 'TYPE      ',IT_TET
         WRITE (IOTMP,FMT=FMT07) (0D0,IR=1,IRSHFT),
     &                           (VT(IR,IT),IR=1,IRTOP)
         WRITE (IOTMP,FMT=FMT07) (0D0,IR=1,IRSHFT),
     &                           (BT(IR,IT),IR=1,IRTOP)
C ----------------------------------------------------------------------
         IF ( FULLPOT ) THEN
            WRITE (IOTMP,99017) 'NFP       ',NFPT(IT)
            WRITE (IOTMP,99017) 'NLMFPT    ',NLMFPT(IT)
            DO LM = 2,NLMFPT(IT)
               IF ( KLMFP(LM,IT).NE.0 ) THEN
                  WRITE (IOTMP,99017) 'LM        ',LM
                  WRITE (IOTMP,FMT=FMT07)
     &                   (VNST(IR,LM,IT),IR=JRNS1(IM),IRTOP)
                  WRITE (IOTMP,FMT=FMT07)
     &                   (BNST(IR,LM,IT),IR=JRNS1(IM),IRTOP)
               END IF
            END DO
         END IF
C ----------------------------------------------------------------------
         WRITE (IOTMP,99012)
      END DO
C ======================================================================
      WRITE (IOTMP,99001) 'CHARGE'
      DO IT_TET = 1,NT_TET
         IT = (IT_TET-1)/NTET_FLUCT + 1
         WRITE (IOTMP,99017) 'TYPE      ',IT_TET
         IM = IMT(IT)
         IF ( FULLPOT ) THEN
            IRTOP = JRCRI(IM)
         ELSE
            IRTOP = JRWS(IM)
         END IF
         WRITE (IOTMP,FMT=FMT07) (RHOCHR(IR,IT),IR=1,IRTOP)
         WRITE (IOTMP,FMT=FMT07) (RHOSPN(IR,IT),IR=1,IRTOP)
C ----------------------------------------------------------------------
         IF ( FULLPOT ) THEN
            WRITE (IOTMP,99013)
            WRITE (IOTMP,99017) 'NFP       ',NFPT(IT)
            WRITE (IOTMP,99017) 'NLMFPT    ',NLMFPT(IT)
            DO LM = 1,NLMFPT(IT)
               IF ( KLMFP(LM,IT).NE.0 ) THEN
                  WRITE (IOTMP,99017) 'LM        ',LM
                  WRITE (IOTMP,FMT=FMT07)
     &                   (RHO2NS(IR,LM,IT,1),IR=1,IRTOP)
                  WRITE (IOTMP,FMT=FMT07)
     &                   (RHO2NS(IR,LM,IT,2),IR=1,IRTOP)
               END IF
            END DO
         END IF
C ----------------------------------------------------------------------
         WRITE (IOTMP,99012)
      END DO
C ======================================================================
      WRITE (IOTMP,99001) 'MOMENTS        QEL  NOS  SMT  OMT  HFF'
      DO IT_TET = 1,NT_TET
         IT = (IT_TET-1)/NTET_FLUCT + 1
         WRITE (IOTMP,99017) 'TYPE      ',IT_TET
         WRITE (IOTMP,FMT=FMT07) QEL(IT),OBS_T(0,IDOS,IT),
     &                           OBS_T(0,ISMT,IT),OBS_T(0,IOMT,IT),
     &                           OBS_T(0,IHFF,IT)
C
         WRITE (IOTMP,99012)
      END DO
C ======================================================================
      IF ( ORBPOL(1:6).EQ.'BROOKS' ) THEN
         WRITE (IOTMP,99001) 'VECTOR-POTENTIAL'
         DO IT_TET = 1,NT_TET
            IT = (IT_TET-1)/NTET_FLUCT + 1
            WRITE (IOTMP,99017) 'TYPE      ',IT_TET
            WRITE (IOTMP,99017) 'LOPT      ',LOPT(IT)
C ----------------------------------------------------------------------
            IF ( LOPT(IT).GT.0 ) THEN
               IM = IMT(IT)
               DO MS = 1,2
                  WRITE (IOTMP,FMT=FMT07) (0D0,IR=1,IRSHFT),
     &                   (AOPT(IR,MS,IT),IR=1,JRWS(IM))
               END DO
            END IF
            WRITE (IOTMP,99012)
         END DO
C ======================================================================
      ELSE IF ( BREITINT ) THEN
         WRITE (IOTMP,99001) 'VECPOT-BREITINT '
         DO IT_TET = 1,NT_TET
            IT = (IT_TET-1)/NTET_FLUCT + 1
            WRITE (IOTMP,99017) 'TYPE      ',IT_TET
C ----------------------------------------------------------------------
            IM = IMT(IT)
            DO MS = -1, + 1
               DO ILA = 1,NLABIMAX
                  WRITE (IOTMP,FMT=FMT07)
     &                   (ABIT(IR,ILA,MS,IT),IR=1,JRWS(IM))
               END DO
            END DO
            WRITE (IOTMP,99012)
         END DO
C ======================================================================
      ELSE IF ( ORBPOL(1:4).EQ.'DMFT' ) THEN
         WRITE (IOTMP,99001) 'SELF-ENERGY'
         DO IT_TET = 1,NT_TET
            IT = (IT_TET-1)/NTET_FLUCT + 1
            WRITE (IOTMP,99017) 'TYPE      ',IT_TET
            WRITE (IOTMP,99017) 'LOPT      ',LOPT(IT)
C ----------------------------------------------------------------------
            IF ( LOPT(IT).GT.0 ) THEN
            END IF
            WRITE (IOTMP,99012)
         END DO
C ======================================================================
      END IF
C ======================================================================
C HE 16/02/14     taken out
C ======================================================================
Cc      WRITE (IOTMP,99001) 'CORE STATES'
Cc      DO IT = 1,NT
Cc         WRITE (IOTMP,99018) 'TYPE      ',IT_TET
Cc         WRITE (IOTMP,99018) 'NCORE     ',NCORE(IT)
Cc         IF ( NCORE(IT).GT.0 ) THEN
Cc            WRITE (IOTMP,99014)
Cc            DO IS = 1,NSPIN
Cc               IQ = (IT-1)*NSPIN + IS
Cc               DO IM = 1,NCORE(IT)
Cc                  WRITE (IOTMP,FMT=FMT09) LCORE(IM,IT),ECORE(IM,IQ)
Cc               END DO
Cc            END DO
Cc         END IF
Cc         WRITE (IOTMP,99012)
Cc      END DO
C ======================================================================
C
      CLOSE (IOTMP)
C
99001 FORMAT (79('*'),/,A)
99002 FORMAT (8X,'IQ        QBAS(X)               QBAS(Y)',
     &        '               QBAS(Z)')
99003 FORMAT (8X,'IQ','     IREFQ','       IMQ','       NOQ','  ITOQ',
     &        '  CONC')
99004 FORMAT (8X,'IQ       MTET_Q',16X,'MPHI_Q',16X,'MGAM_Q')
99005 FORMAT (8X,'IT       MTET_T',16X,'MPHI_T',16X,'MGAM_T')
99006 FORMAT (3X,'IM       R(1)                  DX              JRMT',
     &        '       RMT             JRWS       RWS')
99007 FORMAT (3X,'  IM','      A       ','       B      ',' JRMT',
     &        '     RMT      ',' JRWS','     RWS      ')
99008 FORMAT (3X,'IM',' JRNS1',' JRCRI','  NPAN','  JRCUT')
99009 FORMAT (6X,'IREF       VREF                  RMTREF')
99010 FORMAT (6X,'  IQ      ',A)
99011 FORMAT ('   IT','     TXT_T','            ZT','     NCORT',
     &        '     NVALT','    NSEMCORSHLT')
99012 FORMAT (79('='))
99013 FORMAT (79('-'))
99014 FORMAT (A10,A)
99015 FORMAT (A10,'''',A,'''')
99016 FORMAT (A10,L1)
99017 FORMAT (A10,6I10)
99018 FORMAT (A10,I10,5X,A)
99019 FORMAT (79('*'))
99020 FORMAT (79('*'),/,A,'''',A,A,'''')
      END
