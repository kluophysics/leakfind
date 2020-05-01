C*==transpho_rslab.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE TRANSPHO_RSLAB(MEZJ,MEZZ,MSSQ,SSST,TAUQ,TAUT,PHASK)
C
C   ********************************************************************
C   *                                                                  *
C   *  main driver  for       KKRSPEC                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,ML,MQD,RSTEP,MLS,LL,XMAXE,NFULLPOT,
     &    MLQNAT,MEW,LH,MHDIM2,MLZP,HARTRE,PI,CZERO,NTPHOMAX
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX
      USE MOD_TYPES,ONLY:NTMAX,LTXT_T,TXT_T,LCXRAY,NCXRAY
      USE MOD_SITES,ONLY:NQMAX,NQ,NOQ,ITOQ
      USE MOD_ENERGY,ONLY:NEMAX
      USE MOD_SPEC_MESH,ONLY:RADMESH,RADMESH3,RADMESHM,RADMESH3M
      USE MOD_SPEC_COE1,ONLY:LSPOO,IIOO,POOO,POSOO,IGZOO,IROOO,IAOO
      USE MOD_SPEC_COE2,ONLY:IAOOS,LSPEE,IIEE,PEEE,PESEE,IGZEE,IROEE,
     &    IAEE
      USE MOD_SPEC_COE3,ONLY:IAEES,LSPOE,IIOE,POOE,PESOE,IGZOE,IROOE,
     &    IAOE
      USE MOD_SPEC_COE4,ONLY:IAOES,LSPEO,IIEO,PEEO,POSEO,IGZEO,IROEO,
     &    IAEO
      USE MOD_SPEC_COE5,ONLY:IAEOS,COO,CEE,COE,CEO
      USE MOD_SPEC_COMPOT,ONLY:ZZ,A,B,VV
      USE MOD_SPEC_GEOM,ONLY:SEP,AR1,AR2,RAR
      USE MOD_SPEC_LMKMS,ONLY:LME,LMO,KME,KMO,SW,TAUW,LKAE,LKAO
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      USE MOD_SPEC_POTLM,ONLY:MESH,EB,BB,ALPHA
      USE MOD_SPEC_RINDC,ONLY:KAP,MUD,NRL,NRM,NRS,NR,REL
      USE MOD_FILES,ONLY:IFILSPECSTR,IFILSPECBAR,IFILSPECPOT,
     &    IFILSPECINP,IFILSPECOU3,DATSET,LDATSET,IPRINT,IFILSPECOU1,
     &    IFILSPECOU2,FOUND_SECTION
      USE MOD_LATTICE,ONLY:SYSTEM_DIMENSION,SYSTEM_TYPE,SUB_SYSTEM
      USE MOD_CALCMODE,ONLY:KKRMODE
      USE MOD_MPI,ONLY:MPI,MPI_ID
      USE MOD_CONSTANTS,ONLY:C0,C1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER MLR,MLT,MLA,MLB,MLRNAT,MLTNAT,MAKM
      PARAMETER (MLR=(ML+1)*ML/2,MLT=ML*(ML-1)/2,MLA=2*ML-1,MLB=MLA*MLA,
     &           MLRNAT=MLR*NATLM,MLTNAT=MLT*NATLM,
     &           MAKM=4*MLQNAT*MLQNAT*LAYSM)
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SPEC_RSLAB')
C
C Dummy arguments
C
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),PHASK(NEMAX),
     &           SSST(NKMMAX,NKMMAX,NTMAX),TAUQ(NKMMAX,NKMMAX,NQMAX),
     &           TAUT(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 AA(3),AMAT1(:,:,:,:),AMAT2(:,:,:),AWR(:,:,:,:),BK(:,:),
     &           BULKX(:,:,:),CBX,PHOTOCSD(:,:,:,:,:),PSI2G(:,:,:),Q1,
     &           Q2,Q3,Q4,RMATO(:,:,:,:),RMATS(:,:,:,:),ROT(:,:),
     &           ROTM1(:,:),TSEH(:,:,:),TSEL(:,:,:),TSEOH(:,:,:),
     &           TSEOL(:,:,:),TSOEH(:,:,:),TSOEL(:,:,:),TSOH(:,:,:),
     &           TSOL(:,:,:),VLM(:,:,:,:,:),ZMAT0(:,:,:,:)
      REAL*8 AD1(2),AD2(2),ADD(3),ADU(3),ALQ,ASYM,BARABD,BARABU,BPAR(3),
     &       BRUSEP(3),CELM(:),CLIGHT,DELQ,DELTAP,DELTAT,DTEMP,EBIND(:),
     &       EEND,EFERM,EFEV,EMAX,EMESH,EMIN,EPSX,ESTART,ESTEP,EZ,
     &       FAC1(:),FAC2(:),FIQ,FIXTHQ,FROMSIGMA,FTEMP,FWHM,LOEF,LOIHM,
     &       LOWHM,MASS,OMEV,OMHAR,PHI,PHOTOC(:,:,:,:,:),PKEMAX,PKEMIN,
     &       PMAX,PMAXA,PMIN,PMINA,POL0(3),POL0L(3),POL0_INITIAL(3),
     &       POS(:,:,:),PSPIN(3),RD,ROATLA(:),TEMP,THETA,THQ,TMAX,TMAXA,
     &       TMIN,TMINA,TOLZ,TOSIGMA,VBH,VBL,VBSTP,VIH,VIL,VPR,VPREV,
     &       WFEV,YLM(MLB),Z(:,:,:),ZAD,ZAU,ZBD,ZBU,ZBX,ZMESH(6)
      INTEGER AMAT1V(:,:,:,:),AMAT1X(:,:,:,:),AMAT1Y(:,:,:,:),
     &        AMAT2V(:,:,:),AMAT2X(:,:,:),AMAT2Y(:,:,:),ASIG,CIV,CORE,
     &        GANZ,I,IAN,IAT(:,:),IBAR,IBLOCH,IBXY,ICIRC,IDCALC,IDORA,
     &        IDREH,IEDLAYER(:),IEODLAYER(:),IE_ICST,IFM,IFSP,IGCONV,
     &        IGPRO(:),IGVAL(:),ILATT,INPVER,IODLAYER(:),IOEDLAYER(:),
     &        IPARA,IQSURF,IREL,IROT,IROTMAX,IRUMP(:),ISEQ(:),ISTACK,
     &        ISTR(2),IT,ITXRAY,IUSEDOS,IXAS,IZ,J,K,L,LANZ1,LANZ2,LAYB,
     &        LAYER1,LAYP,LAYS,LFIL(5),LOUT,LSTR,LU(5),MAXL,MCD,
     &        MILLER_INDICES(3),NAT,NATL(:),NCLM,NE(:),NEV(:),NEW,NIN,
     &        NO(:),NOD(:),NP,NPOL,NREL,NROT,NSPIN,NSPLEE,NT,NTOT(:),
     &        PKSCAN
      CHARACTER*80 DATAFILE,INFILE(5),OUTFILE
      CHARACTER*10 IOUT,TASK
      LOGICAL POL0_VEC_TILT,UDT,USE_CRY_PRIM_VECS
      INTEGER PTYPE(:,:,:),SPOL,STRVER,TYP,USEEULER,XMATOFF
      REAL*8 ZPARD(3),ZPARU(3),ZSD,ZSU
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE IEDLAYER,IODLAYER,IEODLAYER,IOEDLAYER,BK,NE,NO
      ALLOCATABLE IAT,AWR,NEV,POS,ROT,NOD,FAC1,FAC2,CELM,NATL
      ALLOCATABLE TSEH,TSEL,TSOH,TSOL,ISEQ,AMAT1
      ALLOCATABLE AMAT2,NTOT,PSI2G,ROTM1,EBIND,IGVAL,IGPRO
      ALLOCATABLE RMATO,BULKX,RMATS,TSEOH,TSEOL,TSOEH,TSOEL,IRUMP
      ALLOCATABLE PTYPE,AMAT1V,AMAT2V,AMAT1X,AMAT1Y,AMAT2X,AMAT2Y
      ALLOCATABLE ROATLA,PHOTOC,PHOTOCSD
      ALLOCATABLE vlm,z,ZMAT0
C
      USE_CRY_PRIM_VECS = .TRUE.
      MILLER_INDICES(1) = 0
      MILLER_INDICES(2) = 0
      MILLER_INDICES(3) = 1
C
      IQSURF = 1
C
      TASK(1:10) = '          '
C     VERSION OF INPUT FILE: 0=OLD INPUT, 1=NEW INPUT
      INPVER = 1
      INFILE(1) = 'in_structur.inp'
      LU(1) = IFILSPECSTR
      LFIL(1) = 15
      INFILE(3) = 'in_bar.inp'
      LU(3) = IFILSPECBAR
      LFIL(3) = 10
      INFILE(4) = 'in_pot.inp'
      LU(4) = IFILSPECPOT
      LFIL(4) = 10
      INFILE(2) = 'input.inp'
      LU(2) = IFILSPECINP
      LFIL(2) = 9
C
C--------------------------------------------------------
C     HOW WE SHELL DEAL WITH THE STRUCTURE FILE:
C     STRVER=0 code will look for in_structure.inp
C              or if set then STRFIL
C              This will work in all SPR-KKR modi e.g.:
C              3D: Standard, TB
C              2D: LIV,VIV and LIR
C
C     STRVER=1 code will try to create own structure file
C              this can work only in 2D:LIV modi
C----------------------------------------------------------
      STRVER = 0
C
      NOUT1 = IFILSPECOU3
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_FIND_KEYWORD('ARPES',UDT)
         IF ( UDT ) TASK(1:10) = 'ARPES     '
         CALL SECTION_FIND_KEYWORD('AIPES',UDT)
         IF ( UDT ) TASK(1:10) = 'AIPES     '
         CALL SECTION_FIND_KEYWORD('SPLEED',UDT)
         IF ( UDT ) TASK(1:10) = 'SPLEED    '
         CALL SECTION_FIND_KEYWORD('BAND',UDT)
         IF ( UDT ) TASK(1:10) = 'BAND      '
         CALL SECTION_FIND_KEYWORD('CLPES',UDT)
         IF ( UDT ) TASK(1:10) = 'CLPES     '
C
         CALL SECTION_SET_INTEGER('INPVER',INPVER,9999,0)
C
         CALL SECTION_SET_INTEGER('STRVER',STRVER,9999,0)
C
         CALL SECTION_SET_STRING('INPFIL',INFILE(2),'9999',0)
         CALL SECTION_SET_STRING('STRFIL',INFILE(1),'9999',0)
         CALL SECTION_SET_STRING('BARFIL',INFILE(3),'9999',0)
C
         IF ( SYSTEM_DIMENSION(1:2).EQ.'3D' ) THEN
            CALL SECTION_FIND_KEYWORD('CRYS_VECS',UDT)
            IF ( UDT ) USE_CRY_PRIM_VECS = .FALSE.
            CALL SECTION_SET_INTEGER('IQ_AT_SURF',IQSURF,9999,0)
            IF ( IQSURF.GT.NQ .OR. IQSURF.LT.1 )
     &            STOP '<SPEC_RSLAB> IQ_AT_SURF .GT. NQ'
            NIN = 3
            CALL SECTION_SET_INTEGER_ARRAY('MILLER_HKL',MILLER_INDICES,
     &         NIN,3,1,9999,0)
            CALL SECTION_SET_REAL('DEL_Z_RUMPLED',TOLZ,99999.0D0,0)
         END IF
C
      END IF
      IF ( TASK(1:2).EQ.'  ' ) STOP 
     &                 '<SPEC_RSLAB> TASK NOT FOUND OR INVALID FOR SPEC'
C
C---------------------------------------------------------------------
C    Select task for the spec calculations
C---------------------------------------------------------------------
C
C     ibloch  = 0  band structure calc.
C                  calculate core-level if core = 1
C             = 1  spleed calculation,
C             = 2  xps (corelevel) calc.
C             = 3  auger calc.
C             = 4  ups (valence-band) calc. (arups)
C             = 5  xes calc.
C             = 6  xas calc.
C             = 7  angular integrated ups
C             = 8  secondary emission spectrum
C             = 9  2ppe => two photon photoemission (work in progress)
C            >= 10 special purpose (see rslab)
C
      NSPLEE = 1
      IBLOCH = 99
      IF ( TASK(1:5).EQ.'ARPES' ) IBLOCH = 4
      IF ( TASK(1:5).EQ.'AIPES' ) IBLOCH = 2
      IF ( TASK(1:6).EQ.'SPLEED' ) THEN
         IBLOCH = 1
         NSPLEE = 0
         IF ( MPI ) CALL STOP_MESSAGE(ROUTINE,'MPI and SPLEED not yet')
      END IF
      IF ( TASK(1:4).EQ.'BAND' ) IBLOCH = 0
      IF ( TASK(1:5).EQ.'CLPES' ) IBLOCH = 3
      IF ( IBLOCH.EQ.99 ) THEN
         WRITE (6,*) '<SPEC_RSLAB>: NOT YET INPLEMENTED TASK=',TASK
         STOP
      END IF
      IF ( TASK(1:5).EQ.'ARPES' ) THEN
         WRITE (6,99004) 'Angle resolved photoemission'
      ELSE IF ( TASK(1:5).EQ.'AIPES' ) THEN
         WRITE (6,99004) 'Angle integrated photoemission (XPS)'
      ELSE IF ( TASK(1:6).EQ.'SPLEED' ) THEN
         WRITE (6,99004) 'Spin polarised LEED'
      ELSE IF ( TASK(1:4).EQ.'BAND' ) THEN
         WRITE (6,99004) 'Band structure calculations (using spec)'
      ELSE IF ( TASK(1:5).EQ.'CLPES' ) THEN
         WRITE (6,99004) 
     &                 'Angle resolved core level photoemission and XPD'
      END IF
C
      IF ( STRVER.EQ.0 ) THEN
         WRITE (6,99007)
         WRITE (6,99008) 
     &          'Old format for structure file of surface will be used.'
         WRITE (6,99008) 'Potentials, scattering t- and tau- matrices'
         WRITE (6,99010) 'are calculated by given KKR MODE:'
C
         WRITE (6,99010) 'SYSDIM    = ',SYSTEM_DIMENSION
         WRITE (6,99010) 'SYSTYPE   = ',SYSTEM_TYPE
         WRITE (6,99010) 'SUBSYSTEM = ',SUB_SYSTEM
         WRITE (6,99010) 'KKRMODE   = ',KKRMODE
C
         WRITE (6,99008) 'Structure file used:'
         WRITE (6,99010) INFILE(1)
C
         WRITE (6,99008) 
     &           'Please check occupation table of given layers bellow:'
C
      ELSE
         IF ( SYSTEM_DIMENSION(1:2).EQ.'3D' .OR. 
     &        (SYSTEM_DIMENSION(1:2).EQ.'2D' .AND. SYSTEM_TYPE(1:3)
     &        .EQ.'LIV') ) THEN
            WRITE (6,99007)
            WRITE (6,99008) 
     &                 '------------ EXPERIMENTAL FEATURE--------------'
            WRITE (6,99010) 
     &                   'automatic creation of structure file for spec'
C
            IF ( SYSTEM_DIMENSION(1:2).EQ.'3D' ) THEN
               WRITE (6,99008) '3D bulk input used to create surface'
               IF ( USE_CRY_PRIM_VECS ) THEN
                  WRITE (6,99009) 
     &                      'Miller indices with respect to cryst. cell'
               ELSE
                  WRITE (6,99009) 
     &                       'Miller indices with respect to prim. cell'
               END IF
               WRITE (6,99009) 'Miller indices (hkl)',MILLER_INDICES(1),
     &                         MILLER_INDICES(2),MILLER_INDICES(3)
C
               WRITE (6,99009) 'Surface will be terminated with IQ:',
     &                         IQSURF
C
               WRITE (6,99010) 'occupied by: ',
     &                         (TXT_T(ITOQ(K,IQSURF))(1:LTXT_T
     &                         (ITOQ(K,IQSURF))),K=1,NOQ(IQSURF))
C
            ELSE
               WRITE (6,99008) '2D input used to create surface'
C
            END IF
         ELSE
            WRITE (6,99007)
            WRITE (6,99008) 
     &              '<SPEC_RSLAB>: Automatic creation of structure file'
            WRITE (6,99010) 
     &              'works only for KKRMODE=TB-KKR with SYSTEM_TYPE=LIV'
            WRITE (6,99010)
            WRITE (6,99007)
            STOP
         END IF
         WRITE (6,99008) 
     &                 '-----------------------------------------------'
C
      END IF
C---------------------------------------------------------------------
C     Open additional output files
C
      OUTFILE = DATSET(1:LDATSET)//'_SPEC.out'
      LOUT = LDATSET + 9
      IF ( MPI_ID.EQ.0 ) THEN
C
         OPEN (UNIT=NOUT1,FILE=OUTFILE(1:LOUT),STATUS='replace')
C
C
      ELSE IF ( IPRINT.GE.1 ) THEN
C
         IOUT = '0000'
         CALL STRING_ADD_N(IOUT,MPI_ID)
         LSTR = LEN_TRIM(IOUT)
         OUTFILE = 'SCRATCH.out'//IOUT((LSTR-3):LSTR)
C
         OPEN (UNIT=NOUT1,FILE=OUTFILE(1:LOUT+4),STATUS='replace')
      ELSE
         OPEN (UNIT=NOUT1,FILE='/dev/null')
C
      END IF
C
      CALL OPENFILES(INPVER,INFILE,LU,LFIL)
C---------------------------------------------------------------------
C
C     Create and open structure file for surface
C
      IF ( STRVER.EQ.1 ) THEN
         IF ( MPI_ID.EQ.0 ) THEN
            OPEN (UNIT=LU(1),FILE=INFILE(1)(1:LFIL(1)))
            CALL SPEC_INPSTRU(INFILE(1),IFILSPECSTR,USE_CRY_PRIM_VECS,
     &                        MILLER_INDICES,IQSURF,TOLZ)
         END IF
         IF ( MPI ) CALL DRV_MPI_BARRIER
         IF ( MPI .AND. MPI_ID.NE.0 )
     &        OPEN (UNIT=LU(1),FILE=INFILE(1)(1:LFIL(1)))
      ELSE
         INQUIRE (FILE=INFILE(1)(1:LFIL(1)),EXIST=UDT)
         IF ( .NOT.UDT .AND. STRVER.EQ.0 ) THEN
            WRITE (6,99011) INFILE(1)(1:LFIL(1))
            WRITE (NOUT1,99011) INFILE(1)(1:LFIL(1))
            STOP
         ELSE
            OPEN (UNIT=LU(1),FILE=INFILE(1)(1:LFIL(1)))
         END IF
C
      END IF
C
      CALL SPEC_CONVERTPOT()
C
      ALLOCATE (IEDLAYER(0:LAYSM),IODLAYER(0:LAYSM))
      ALLOCATE (IEODLAYER(0:LAYSM),IOEDLAYER(0:LAYSM),BK(LH,2))
      ALLOCATE (NE(LAYSM),NO(LAYSM))
      ALLOCATE (IAT(NATLM,LAYSM),AWR(MLS,LL,NATLM,2))
      ALLOCATE (NEV(LAYSM),POS(3,NATLM,LAYSM),ROT(MQD,MQD))
      ALLOCATE (NOD(LAYSM),FAC1(MLA),FAC2(MLB))
      ALLOCATE (CELM(MLZP),NATL(LAYSM),TSEH(MLQNAT,MLQNAT,LAYSM))
      ALLOCATE (TSEL(MLQNAT,MLQNAT,LAYSM),TSOH(MLQNAT,MLQNAT,LAYSM))
      ALLOCATE (TSOL(MLQNAT,MLQNAT,LAYSM))
      ALLOCATE (ISEQ(0:LAYSM))
      ALLOCATE (AMAT1(XMAXE,4,NFULLPOT,2))
      ALLOCATE (AMAT2(XMAXE,4,NFULLPOT),NTOT(LAYSM))
      ALLOCATE (PSI2G(MHDIM2,LH,4),ROTM1(MQD,MQD),EBIND(MEW))
      ALLOCATE (IGVAL(LH),IGPRO(LH),RMATO(MQD,MQD,LAYSM,NATLM))
      ALLOCATE (BULKX(LH,LH,LAYSM),RMATS(MQD,MQD,LAYSM,NATLM))
      ALLOCATE (TSEOH(MLQNAT,MLQNAT,LAYSM))
      ALLOCATE (TSEOL(MLQNAT,MLQNAT,LAYSM))
      ALLOCATE (TSOEH(MLQNAT,MLQNAT,LAYSM))
      ALLOCATE (TSOEL(MLQNAT,MLQNAT,LAYSM),IRUMP(LAYSM))
      ALLOCATE (PTYPE(NATLM,LAYSM,3),AMAT1V(XMAXE,4,NFULLPOT,2))
      ALLOCATE (AMAT2V(XMAXE,4,NFULLPOT),AMAT1X(MQD+1,4,NFULLPOT,2))
      ALLOCATE (AMAT1Y(MQD+1,4,NFULLPOT,2),AMAT2X(MQD+1,4,NFULLPOT))
      ALLOCATE (AMAT2Y(MQD+1,4,NFULLPOT),ROATLA(LL))
      ALLOCATE (VLM(LAYSM,NATLM,MQD,NRMAX,NTPHOMAX))
      ALLOCATE (ZMAT0(MQD,MQD,LAYSM,NATLM))
C
      ALLOCATE (Z(NATLM,LAYSM,NTMAX))
      TSEL(:,:,:) = C0
      TSEOL(:,:,:) = C0
      TSOEL(:,:,:) = C0
      TSOL(:,:,:) = C0
      TSEH(:,:,:) = C0
      TSEOH(:,:,:) = C0
      TSOEH(:,:,:) = C0
      TSOH(:,:,:) = C0
C
C
C
      ALLOCATE (RADMESHM(NRMAX,NATLM,LAYSM))
      ALLOCATE (RADMESH3M(NRMAX,NATLM,LAYSM))
      ALLOCATE (RADMESH(RSTEP,NATLM,LAYSM))
      ALLOCATE (RADMESH3(RSTEP,NATLM,LAYSM))
      RADMESHM(1:NRMAX,1:NATLM,1:LAYSM) = 0.0D0
      RADMESH3M(1:NRMAX,1:NATLM,1:LAYSM) = 0.0D0
      RADMESH(1:RSTEP,1:NATLM,1:LAYSM) = 0.0D0
      RADMESH3(1:RSTEP,1:NATLM,1:LAYSM) = 0.0D0
C
      ALLOCATE (A(4,2),ALPHA(LAYSM,NATLM),AR1(2),AR2(2),B(4,2),BB(3))
      ALLOCATE (CEE(MAKM),CEO(MAKM),COE(MAKM))
      ALLOCATE (COO(MAKM),EB(3),KME(MLQNAT,2))
      ALLOCATE (KMO(MLQNAT,2),MESH(LAYSM,NATLM,RSTEP),NRS(MQD),RAR(4))
      ALLOCATE (SEP(3,0:LAYSM),SW(2),TAUW(2),VV(5,2,2),ZZ(7,2))
      ALLOCATE (IAEE(MAKM),IAEES(MAKM),IAEO(MAKM))
      ALLOCATE (IAEOS(MAKM),IAOE(MAKM))
      ALLOCATE (IAOES(MAKM),IAOO(MAKM),IAOOS(MAKM))
      ALLOCATE (IGZEE(MAKM),IGZEO(MAKM))
      ALLOCATE (IGZOE(MAKM),IGZOO(MAKM),IIEE(MAKM))
      ALLOCATE (IIEO(MAKM),IIOE(MAKM))
      ALLOCATE (IIOO(MAKM),IROEE(MAKM),IROEO(MAKM),IROOE(MAKM))
      ALLOCATE (IROOO(MAKM),KAP(MQD),LKAE(MLQNAT),LKAO(MLQNAT))
      ALLOCATE (LME(MLRNAT,2),LMO(MLTNAT,2),LSPEE(MAKM),LSPEO(MAKM))
      LME(1:MLRNAT,1:2) = 0
      LMO(1:MLTNAT,1:2) = 0
      ALLOCATE (LSPOE(MAKM),LSPOO(MAKM),MUD(MQD),NR(MQD),NRL(MQD))
      ALLOCATE (NRM(MQD),PEEE(MAKM),PEEO(MAKM),PESEE(MAKM),PESOE(MAKM))
      ALLOCATE (POOE(MAKM),POOO(MAKM),POSEO(MAKM))
      ALLOCATE (POSOO(MAKM),REL(MQD,2))
C
      PSPIN(1:3) = 0.0D0
      PSPIN(3) = 1.0D0
C
C     read parameter for barrier potential
C     ====================================
      CALL POTBAR(IBAR,ZPARU,ZPARD,BPAR,EPSX,BARABU,BARABD,INPVER)
C
C
C     read structure-file
C
C
      NROT = 0
      IF ( STRVER.EQ.0 .OR. STRVER.EQ.1 )
     &     CALL INPSTRU(NATL,POS,IRUMP,ISEQ,LAYS,LAYB,IAT,BRUSEP,1.0D-6,
     &     AD1,AD2,ILATT,CIV,NROT,ISTACK,BARABU,BARABD)
C
C     read input-data-file
C
      IF ( INPVER.EQ.0 ) THEN
         CALL INPUT(EMESH,ESTEP,EFEV,WFEV,VIH,VIL,OMEV,EFERM,VPREV,THQ,
     &              FIQ,ALQ,NOUT1,IP,NREL,NPOL,NSPLEE,NSPIN,SPOL,LAYER1,
     &              GANZ,LANZ1,LANZ2,TMIN,TMAX,PMIN,PMAX,NT,NP,TYP,ISTR,
     &              TEMP,DTEMP,MASS,POL0,POL0L,IDORA,IDREH,ICIRC,IFSP,
     &              Q1,Q2,Q3,Q4,IBLOCH,CORE,ESTART,EEND,ASYM,XMATOFF,
     &              PKSCAN,PKEMIN,PKEMAX,ASIG,FROMSIGMA,TOSIGMA,EMIN,
     &              EMAX,VBH,VBL,VBSTP,NATL,LAYS,DELQ,MCD,USEEULER,
     &              IPARA,IXAS,LOWHM,LOIHM,LOEF,IUSEDOS,IDCALC,IGCONV,
     &              FWHM,FTEMP)
      ELSE
         CALL SPEC_INPUT(EMESH,ESTEP,EFEV,VIH,VIL,OMEV,EFERM,VPREV,THQ,
     &                   FIQ,ALQ,NOUT1,IP,NREL,NPOL,NSPLEE,NSPIN,SPOL,
     &                   LAYER1,GANZ,LANZ1,LANZ2,TMIN,TMAX,PMIN,PMAX,NT,
     &                   NP,TYP,ISTR,TEMP,DTEMP,MASS,POL0,POL0L,IDORA,
     &                   IDREH,ICIRC,IFSP,Q1,Q2,Q3,Q4,IBLOCH,CORE,
     &                   ESTART,EEND,ASYM,XMATOFF,PKSCAN,PKEMIN,PKEMAX,
     &                   ASIG,FROMSIGMA,TOSIGMA,EMIN,EMAX,VBH,VBL,VBSTP,
     &                   NATL,LAYS,DELQ,MCD,USEEULER,IPARA,IXAS,LOWHM,
     &                   LOIHM,LOEF,IUSEDOS,IDCALC,IGCONV,FWHM,FTEMP,
     &                   FIXTHQ,ITXRAY,LCXRAY,NCXRAY,IE_ICST,PSPIN,
     &                   POL0_VEC_TILT,POL0_INITIAL)
C
      END IF
C
      CALL RIND(ML,MQD,KAP,MUD,NRL,NRM,NRS,NR,REL)
      RD = PI/180.D0
C
C     convert to atomic units
C
      EMIN = EMIN/HARTRE
      EMAX = EMAX/HARTRE
      ESTART = ESTART/HARTRE
      EEND = EEND/HARTRE
      VBSTP = VBSTP/HARTRE
      VBH = VBH/HARTRE
      VBL = VBL/HARTRE
      VPR = VPREV/HARTRE
      OMHAR = OMEV/HARTRE
C
      NEW = IDINT(EMESH)
C      NEW = EMESH
C
C     prepare potential barrier
      CALL POTIN(VPR,VIH,1,ZPARU,ZPARD,ZAU,ZSU,ZBU,ZAD,ZSD,ZBD)
      CALL POTIN(VPR,VIL,2,ZPARU,ZPARD,ZAU,ZSU,ZBU,ZAD,ZSD,ZBD)
C
C     print the potential barrier to a file
      IF ( IP.GT.1 ) THEN
         OPEN (30,FILE='potbaru.dat')
         WRITE (30,99003)
C
         EZ = 0.D0
         DO IZ = 0,100
            ZBX = (-8.+IZ*12./100.)/0.52917706
            CALL CBPOT(EZ,ZBX,1,CBX,1)
            WRITE (30,99002) ZBX*0.52917706,CBX*13.605804
         END DO
C
         OPEN (31,FILE='potbard.dat')
         WRITE (31,99003)
C
         EZ = 0.D0
         DO IZ = 0,100
            ZBX = (-8.+IZ*12./100.)/0.52917706
            CALL CBPOT(EZ,ZBX,1,CBX,2)
            WRITE (31,99002) ZBX*0.52917706,CBX*13.605804
         END DO
C
      END IF
      CLOSE (30)
      CLOSE (31)
C
      ALLOCATE (PHOTOC(NEW,NT,NP,4,2))
      ALLOCATE (PHOTOCSD(NEW,NT,NP,4,2))
C
      PHOTOC(1:NEW,1:NT,1:NP,1:4,1:2) = 0.0D0
      PHOTOCSD(1:NEW,1:NT,1:NP,1:4,1:2) = CZERO
C
C     avoid symmetry calculations if bandstructure is calculated
C
      IF ( IBLOCH.EQ.0 .OR. NROT.LT.1 ) THEN
         IROTMAX = 1
      ELSE
         IROTMAX = CIV
      END IF
C
C     start of loop that respects the symmetry of the surface
C
      DO IROT = 1,IROTMAX
C        rotate structure if necessary
         IF ( IROT.GT.1 ) CALL ROTSEP(ILATT,CIV,IROT,ISTACK,LAYS,LAYB,
     &                                ISEQ,SEP,BRUSEP)
C
         LAYP = LAYS + (LAYER1-1)*(LAYS-LAYB+1)
         IF ( (LAYP+1).GT.LL ) THEN
            WRITE (NOUT1,99001)
            STOP
         END IF
         CALL INPCNT(ISEQ,LAYER1,GANZ,LAYS,LAYB,LAYP)
C
         CALL CELMG(CELM,YLM,FAC2,FAC1,2,NCLM)
C
         IF ( NREL.NE.1 ) THEN
            CLIGHT = 1.D+15
            IREL = 1
         ELSE
            CLIGHT = 137.036*1.0D0
            IREL = 2
         END IF
         IF ( CLIGHT.GT.137.036D0 ) WRITE (*,*)
     &         'speed of light increased'
C
         DO IT = 1,LAYS
            NAT = NATL(IT)
            DO IAN = 1,NAT
               NO(IT) = MLT
               NOD(IT) = NATL(IT)*NO(IT)
               NE(IT) = MLR
               NEV(IT) = NATL(IT)*NE(IT)
               NTOT(IT) = NOD(IT) + NEV(IT)
            END DO
         END DO
         MAXL = ML - 1
C
         DO I = 1,MLQNAT
            DO J = 1,MLQNAT
               DO K = 1,LAYSM
                  TSEL(I,J,K) = CZERO
                  TSEOL(I,J,K) = CZERO
                  TSOEL(I,J,K) = CZERO
                  TSOL(I,J,K) = CZERO
                  TSEH(I,J,K) = CZERO
                  TSEOH(I,J,K) = CZERO
                  TSOEH(I,J,K) = CZERO
                  TSOH(I,J,K) = CZERO
               END DO
            END DO
         END DO
C
         DO I = 1,MQD
            DO J = 1,MQD
               DO K = 1,LAYSM
                  DO L = 1,NATLM
                     ZMAT0(I,J,K,L) = CZERO
                  END DO
               END DO
            END DO
         END DO
C
C        nmesh = radial mesh points for the bulk muffin-tin potentials
C        is fixed to rstep=721 and rstepp5=726 in parms.h
C
         IFM = 0
         IBXY = 0
C
C        determine radial meshs for all atoms:
C        =====================================
         CALL DETRADMESHM(RADMESHM,RADMESH3M,NATL,LAYS,Z)
C
         CALL POTFULLM2(LAYS,NATL,Z,IFM,VLM)
C
C        ifm is for later use
C        IFM = 1
C
         CALL CHECKPOT(LAYS,NATL,PTYPE)
C
C        set useeuler and ipara for selection rules
C        call selecteuler(ibloch, eb, bb, npara, ipara, mcd, useeuler,
C    *                     rot, rotm1, ibxy)
C
C     setup rotation matrices to identity
         ROT(:,:) = C0
         ROTM1(:,:) = C0
         DO I = 1,MQD
            ROT(I,I) = C1
            ROTM1(I,I) = C1
         END DO
C
C         IF ( BB(1).NE.0.D0 .OR. BB(2).NE.0.D0 ) IBXY = 1
         IF ( ABS(BB(1)).GT.1.0D-16 .OR. ABS(BB(2)).GT.1.0D-16 )
     &        IBXY = 1
C
C        set parameter for angular scans
C        ===============================
         TMINA = TMIN*RD
         TMAXA = TMAX*RD
         PMINA = PMIN*RD
         PMAXA = PMAX*RD
         DELTAT = TMAXA - TMINA
         IF ( NT.GT.1 ) DELTAT = DELTAT/DBLE(NT-1)
         DELTAP = PMAXA - PMINA
         IF ( NP.GT.1 ) DELTAP = DELTAP/DBLE(NP-1)
         THETA = TMINA
         PHI = PMINA
C
C
C
C---------------------------------------------------------------------
C     Open datafiles for results
C
C
         WRITE (6,99005) OUTFILE
C
         IF ( MPI_ID.EQ.0 ) THEN
C     for each photon polarisation (linear or circular dichroism)
C     open additional file
            IF ( NPOL.EQ.1 ) THEN
C
               DATAFILE = DATSET(1:LDATSET)//'_data.spc'
               LOUT = LDATSET + 10
C
               CALL WRHEAD(IFILSPECOU1,DATAFILE,TASK,NEW)
C     OPEN (UNIT=IFILSPECOU1,FILE=DATAFILE(1:LOUT),
C     &               STATUS='replace')
               WRITE (IFILSPECOU1,99012) 'NT        ',NT
               WRITE (IFILSPECOU1,99012) 'NP        ',NP
               WRITE (IFILSPECOU1,99013) '#Decription of the columns'
               IF ( IBLOCH.EQ.1 ) THEN
                  WRITE (IFILSPECOU1,99013) '#1.-Energy, 2.-Theta'
                  WRITE (IFILSPECOU1,99013) '#3.-RINTG(G), 4.ASYMM(G)'
                  WRITE (IFILSPECOU1,99013) '#5.-P(1,G), 6.P(2,G)'
                  WRITE (IFILSPECOU1,99013) '#7.-P(3,G)'
C
               ELSE IF ( NEW.GT.1 ) THEN
                  WRITE (IFILSPECOU1,99013) '#1.-Theta, 2.-Energy'
                  WRITE (IFILSPECOU1,99013)
     &                    '#3.-Total Intensity, 4.-Spin up'
                  WRITE (IFILSPECOU1,99013)
     &                    '#5.-Spin down, 6.-Spin pol.'
                  WRITE (IFILSPECOU1,99013)
     &                    '#7.-k_paralel (pi/A), 8.-determinant'
                  WRITE (IFILSPECOU1,99013) '#######################'
               ELSE IF ( NT.GT.1 .AND. NP.GT.1 ) THEN
                  WRITE (IFILSPECOU1,99013)
     &                    '#1.-k_x (pi/A), 2.-k_y (pi/A)'
                  WRITE (IFILSPECOU1,99013)
     &                    '#3.-Total Intensity, 4.-Spin up'
                  WRITE (IFILSPECOU1,99013)
     &                    '#5.-Spin down, 6.-Spin pol.'
                  WRITE (IFILSPECOU1,99013) '#######################'
               END IF
               WRITE (6,99006) DATAFILE
            ELSE
               DATAFILE = DATSET(1:LDATSET)//'_data_pol1.spc'
               LOUT = LDATSET + 14
               WRITE (6,99006) DATAFILE
C
               CALL WRHEAD(IFILSPECOU1,DATAFILE,TASK,NEW)
               WRITE (IFILSPECOU1,99012) 'NT        ',NT
               WRITE (IFILSPECOU1,99012) 'NP        ',NP
               WRITE (IFILSPECOU1,99013) '#Decription of the columns'
               IF ( IBLOCH.EQ.1 ) THEN
                  WRITE (IFILSPECOU1,99013) '#1.-Energy, 2.-Theta'
                  WRITE (IFILSPECOU1,99013) '#3.-RINTG(G), 4.ASYMM(G)'
                  WRITE (IFILSPECOU1,99013) '#5.-P(1,G), 6.P(2,G)'
                  WRITE (IFILSPECOU1,99013) '#7.-P(3,G)'
               ELSE IF ( NEW.GT.1 ) THEN
                  WRITE (IFILSPECOU1,99013) '#1.-Theta, 2.-Energy'
                  WRITE (IFILSPECOU1,99013)
     &                    '#3.-Total Intensity, 4.-Spin up'
                  WRITE (IFILSPECOU1,99013)
     &                    '#5.-Spin down, 6.-Spin pol.'
                  WRITE (IFILSPECOU1,99013)
     &                    '#7.-k_paralel (pi/A), 8.-determinant'
                  WRITE (IFILSPECOU1,99013) '#######################'
               ELSE IF ( NT.GT.1 .AND. NP.GT.1 ) THEN
                  WRITE (IFILSPECOU1,99013)
     &                    '#1.-k_x (pi/A), 2.-k_y (pi/A)'
                  WRITE (IFILSPECOU1,99013)
     &                    '#3.-Total Intensity, 4.-Spin up'
                  WRITE (IFILSPECOU1,99013)
     &                    '#5.-Spin down, 6.-Spin pol.'
                  WRITE (IFILSPECOU1,99013) '#######################'
               END IF
C
C
C     OPEN (UNIT=IFILSPECOU1,FILE=DATAFILE(1:LOUT),
C     &               STATUS='replace')
               DATAFILE = DATSET(1:LDATSET)//'_data_pol2.spc'
               LOUT = LDATSET + 14
C
               CALL WRHEAD(IFILSPECOU2,DATAFILE,TASK,NEW)
               WRITE (IFILSPECOU2,99012) 'NT        ',NT
               WRITE (IFILSPECOU2,99012) 'NP        ',NP
               WRITE (IFILSPECOU2,99013) '#Decription of the columns'
               IF ( IBLOCH.EQ.1 ) THEN
                  WRITE (IFILSPECOU2,99013) '#1.-Energy, 2.-Theta'
                  WRITE (IFILSPECOU2,99013) '#3.-RINTG(G), 4.ASYMM(G)'
                  WRITE (IFILSPECOU2,99013) '#5.-P(1,G), 6.P(2,G)'
                  WRITE (IFILSPECOU2,99013) '#7.-P(3,G)'
               ELSE IF ( NEW.GT.1 ) THEN
                  WRITE (IFILSPECOU2,99013) '#1.-Theta, 2.-Energy'
                  WRITE (IFILSPECOU2,99013)
     &                    '#3.-Total Intensity, 4.-Spin up'
                  WRITE (IFILSPECOU2,99013)
     &                    '#5.-Spin down, 6.-Spin pol.'
                  WRITE (IFILSPECOU2,99013)
     &                    '#7.-k_paralel (pi/A), 8.-determinant'
                  WRITE (IFILSPECOU1,99013) '#######################'
               ELSE IF ( NT.GT.1 .AND. NP.GT.1 ) THEN
                  WRITE (IFILSPECOU2,99013)
     &                    '#1.-k_x (pi/A), 2.-k_y (pi/A)'
                  WRITE (IFILSPECOU2,99013)
     &                    '#3.-Total Intensity, 4.-Spin up'
                  WRITE (IFILSPECOU2,99013)
     &                    '#5.-Spin down, 6.-Spin pol.'
                  WRITE (IFILSPECOU2,99013) '#######################'
               END IF
C
C     OPEN (UNIT=IFILSPECOU2,FILE=OUTFILE(1:LOUT),
C     &               STATUS='replace')
               WRITE (6,99006) DATAFILE
            END IF
         END IF
C
C
         CALL SPEC_INIT_TRANSPHO(NKMMAX,ML)
C
C        select subroutine for calculation
         SELECT CASE (IBLOCH)
C
C IBLOCH = 0: bandstructure calculation
C IBLOCH = 1: spleed calculation
C IBLOCH = 2: VB-XPS photoemission calculation
C IBLOCH = 4: valence band photoemission calculation
C
         CASE (0,1,2,4)
C
            CALL SPEC_UPSRUN(IBLOCH,NOUT1,IP,NREL,IREL,NSPIN,SPOL,
     &                       NSPLEE,EFEV,VPR,VIH,VIL,OMEV,OMHAR,THQ,FIQ,
     &                       ALQ,DELQ,ICIRC,IDREH,MCD,NATL,LAYS,LAYB,
     &                       LAYP,LAYER1,GANZ,LANZ1,LANZ2,POL0,POL0L,
     &                       TYP,ISTR,Q1,Q2,Q3,Q4,NEW,PKSCAN,PKEMIN,
     &                       PKEMAX,ESTEP,THETA,PHI,NT,NP,DELTAT,DELTAP,
     &                       TMINA,PMINA,IFM,POS,IRUMP,BRUSEP,NCLM,CELM,
     &                       CLIGHT,MAXL,NOD,NEV,NTOT,USEEULER,IBXY,Z,
     &                       ZPARU,ZPARD,EPSX,EBIND,PHOTOCSD,PHOTOC,
     &                       ROATLA,RMATO,RMATS,AA,AWR,ZMAT0,AMAT1X,
     &                       AMAT1Y,AMAT1V,AMAT1,AMAT2X,AMAT2Y,AMAT2V,
     &                       AMAT2,PSI2G,BULKX,IGPRO,IGVAL,TSEL,TSOL,
     &                       TSEOL,TSOEL,TSEH,TSOH,TSEOH,TSOEH,IEDLAYER,
     &                       IODLAYER,IEODLAYER,IOEDLAYER,ZMESH,ADU,ADD,
     &                       BK,IFSP,EFERM,RADMESHM,VLM,MEZJ,MEZZ,MSSQ,
     &                       SSST,TAUQ,TAUT,PHASK,FIXTHQ,PSPIN,
     &                       POL0_VEC_TILT,POL0_INITIAL)
C
C        core level photomeission and CL XPD
         CASE (3)
C
            CALL SPEC_CLPES(IBLOCH,NOUT1,IP,NREL,IREL,NSPIN,SPOL,NSPLEE,
     &                      EFEV,VPR,VIH,VIL,OMHAR,THQ,FIQ,ALQ,DELQ,
     &                      NPOL,ICIRC,IDREH,NATL,LAYS,LAYB,LAYP,LAYER1,
     &                      GANZ,LANZ1,POL0,POL0L,TYP,ISTR,Q1,Q2,Q3,Q4,
     &                      NEW,ESTEP,THETA,PHI,NT,NP,DELTAT,DELTAP,
     &                      TMINA,PMINA,IFM,POS,IRUMP,BRUSEP,NCLM,CELM,
     &                      CLIGHT,MAXL,NOD,NEV,NTOT,USEEULER,IBXY,ROT,
     &                      ROTM1,Z,ZPARU,ZPARD,EPSX,EBIND,PHOTOC,
     &                      ROATLA,AA,AWR,ZMAT0,AMAT1X,AMAT1Y,AMAT1V,
     &                      AMAT1,AMAT2X,AMAT2Y,AMAT2V,AMAT2,PSI2G,
     &                      BULKX,IGPRO,IGVAL,TSEL,TSOL,TSEOL,TSOEL,
     &                      TSEH,TSOH,TSEOH,TSOEH,IEDLAYER,IODLAYER,
     &                      IEODLAYER,IOEDLAYER,ZMESH,ADU,ADD,BK,IFSP,
     &                      VLM,MEZJ,MEZZ,MSSQ,SSST,TAUQ,TAUT,PHASK,
     &                      FIXTHQ,ITXRAY,LCXRAY,NCXRAY,ASIG,FROMSIGMA,
     &                      TOSIGMA,EMIN,EMAX,ASYM,IE_ICST,
     &                      POL0_VEC_TILT,POL0_INITIAL)
C
         CASE DEFAULT
            WRITE (6,*) 'no valid ibloch flag, ibloch=',IBLOCH
            WRITE (NOUT1,*) 'no valid ibloch flag, ibloch=',IBLOCH
            STOP
         END SELECT
C     end of rotation for civ symmetry
      END DO
C
      RETURN
C
99001 FORMAT (1x,'layp+1.gt.ll')
99002 FORMAT (3(1x,e11.4))
99003 FORMAT (1x,'zv',10x,'vre',9x,'vim')
99004 FORMAT (/,79('*'),/,24X,A,/,79('*'),//)
99005 FORMAT (/,10X,'*************************',/,10X,
     &        'calculating      SPECTRUM',/,10X,
     &        '*************************',/,10X,
     &        'for more datails see also file:',/,10X,A,/)
99006 FORMAT (/,10X,'*************************',/,10X,
     &        'Results writen in files: ',/,10X,
     &        '*************************',/,10X,' ',/,10X,A,/)
99007 FORMAT (/,1X,79('-'),//)
99008 FORMAT (//,10X,A)
99009 FORMAT (10X,A,5I10)
99010 FORMAT (10X,A,5A)
99011 FORMAT (2x,'** error file:',2x,a15,1x,'does not exist')
99012 FORMAT (A10,I10)
99013 FORMAT (10A)
C
      END
