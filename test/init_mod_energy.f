C*==init_mod_energy.f    processed by SPAG 6.70Rc at 08:49 on  8 Mar 2017
      SUBROUTINE INIT_MOD_ENERGY(IPRINT,TASK,ITEST,NCPA,INITELOOP,
     &                           SCFSTATUS,RDTAU,RDTAUMQ,NOWRDOS,NETAU)
C   ********************************************************************
C   *                                                                  *
C   *  initialize all variables used for   ENERGY   contours           *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKMMAX
      USE MOD_TYPES,ONLY:NTMAX
      USE MOD_CALCMODE,ONLY:PROGNAME
      USE MOD_ENERGY,ONLY:NETAB,EMIN,EMAX,EFERMI,EIMAG,SEARCHEF,IGRID,
     &    ETAB,WETAB,SPLITSS,IMEMIN,IMEMAX,EILOW,NEPATH,NEMAX,PHASA,
     &    PHASK,PHAST,NEPOL,NEFD1,NEFD2,NEFD3,emin_contn,emax_contn,
     &    necontn,lactive_contn,lepath_contn,ime_contn,nepol_contn,
     &    around_fermi,temp_around_fermi
      USE MOD_CONSTANTS,ONLY:RY_EV
      USE MOD_FILES,ONLY:FOUND_SECTION,FOUND_REAL
      USE MOD_KSPACE,ONLY:NGFEPMAX,GFEP,NGFEP
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='INIT_MOD_ENERGY')
C
C Dummy arguments
C
      LOGICAL INITELOOP,NOWRDOS,RDTAU,RDTAUMQ
      INTEGER IPRINT,ITEST,NCPA,NETAU
      CHARACTER*10 SCFSTATUS,TASK
C
C Local variables
C
      REAL*8 EMAXEV,EMINDEFAULT,EMINEV,ERANGE
      INTEGER IA_ERR,IE,IP,ISEARCHEF,N
      CHARACTER*40 PATH(0:11)
      LOGICAL SIGMA_TENSOR,UDT1
      real*8 kb
      data kb /0.6333659D-5/
C
C*** End of declarations rewritten by SPAG
C
      DATA EMINDEFAULT/ - 0.2D0/
      DATA PATH/'ONLY REAL ENERGIES   Gauss-Legendre    ',
     &     'ONLY REAL ENERGIES   Trapez rule       ',
     &     'RECTANGULAR COMPLEX PATH               ',
     &     'STRAIGHT CPLX. PATH PARA. TO REAL AXIS ',
     &     'RECT. CPLX. GRID, RETURN ON LOG SCALE  ',
     &     'ARC IN COMPLEX PLANE                   ',
     &     'STANDARD X-RAY E-MESH                  ',
     &     'X-RAY E-MESH for integration           ',
     &     'ARC IN COMPLEX PLANE - LOGARITHMIC     ',
     &     'ARC IN COMPLEX PLANE - AKAI            ',
     &     'ELLIPSE IN COMPLEX PLANE - DK          ',
     &     'FERMI-DIRAC PATH WITH POLES            '/
C
      EILOW = 999999D0
C-----------------------------------------------------------------------
C
      SEARCHEF = .FALSE.
      ISEARCHEF = 0
      EMIN = EMINDEFAULT
C
      IF ( TASK(1:4).EQ.'SCF ' ) THEN
         IGRID(1) = 5
         NETAB(1) = 30
      END IF
C
      IF ( TASK(1:7).EQ.'COMPTON' ) THEN
         SPLITSS = .TRUE.
         IGRID(1) = 5
         IGRID(2) = 5
         NETAB(1) = 30
         NETAB(2) = 50
      END IF
C-------------------- allow tasks that need ELOOP to be done in parallel
      IF ( TASK(1:2).EQ.'TL' ) THEN
         IGRID(1) = 5
         NETAB(1) = 30
      END IF
C
      IF ( TASK(1:3).EQ.'XAS' ) THEN
         IGRID(1) = 6
         NETAB(1) = 180
         EMAX = 4D0
      END IF
C
      IF ( TASK(1:3).EQ.'CHI' .OR. TASK(1:3).EQ.'OPM' ) THEN
         IGRID(1) = 5
         NETAB(1) = 30
         EMIN = EMINDEFAULT
         EMAX = EFERMI
         EIMAG = 0.01D0
      END IF
C
      IF ( TASK(1:5).EQ.'SIGMA' .OR. TASK(1:6).EQ.'STONER' ) THEN
         IGRID(1) = 0
         NETAB(1) = 1
         EMIN = EFERMI
         EMAX = EFERMI
         EIMAG = 0.0D0
C
         IF ( TASK(1:5).EQ.'SIGMA' ) THEN
            CALL INPUT_FIND_SECTION('TASK',0)
C
            CALL SECTION_FIND_KEYWORD('TENSOR',SIGMA_TENSOR)
            IF ( SIGMA_TENSOR ) THEN
               NETAB(1) = 30
               IGRID(1) = 0
               EMIN = EMINDEFAULT
            END IF
         END IF
C
      END IF
C
      IF ( TASK(1:5).EQ.'CLXPS' ) THEN
         IGRID(1) = 3
         NETAB(1) = 18
C                      !Possible maximum if initial core is f
         EIMAG = 0.01D0
      END IF
C
C=======================================================================
C                overwrite defaults for energy mesh parameters
C=======================================================================
C
      CALL INPUT_FIND_SECTION('ENERGY',0)
C
      IF ( FOUND_SECTION ) THEN
         IF ( .NOT.SPLITSS ) THEN
            CALL SECTION_SET_INTEGER('GRID',IGRID(1),9999,0)
            CALL SECTION_SET_INTEGER('NE',NETAB(1),9999,0)
         ELSE
            CALL SECTION_SET_INTEGER_ARRAY('GRID',IGRID,N,2,2,9999,0)
            CALL SECTION_SET_INTEGER_ARRAY('NE',NETAB,N,2,2,9999,0)
            IF ( IGRID(2).EQ.0 ) IGRID(2) = 0
            IF ( NETAB(2).EQ.0 ) NETAB(2) = 200
         END IF
C
         IF ( IGRID(1).LT.0 .OR. IGRID(1).GT.11 )
     &        CALL STOP_MESSAGE(ROUTINE,'IGRID(1) <0 OR > 9')
         IF ( IGRID(2).LT.0 .OR. IGRID(2).GT.9 )
     &        CALL STOP_MESSAGE(ROUTINE,'IGRID(2) <0 OR > 9')
C
         CALL SECTION_SET_REAL('EMIN',EMIN,9999D0,0)
         CALL SECTION_SET_REAL('EMAX',EMAX,9999D0,0)
         CALL SECTION_SET_REAL('EMINEV',EMINEV,9999D0,0)
         UDT1 = FOUND_REAL
         CALL SECTION_SET_REAL('EMAXEV',EMAXEV,9999D0,0)
         IF ( UDT1 .AND. FOUND_REAL .AND. EFERMI.GT.-9900D0 ) THEN
            EMIN = EFERMI + EMINEV/RY_EV
            EMAX = EFERMI + EMAXEV/RY_EV
         END IF
c modified by XJQ:
         CALL SECTION_SET_REAL('AROUND_FERMI',around_fermi,9999d0,0)
         CALL SECTION_SET_REAL('TEMP_AROUND_FERMI',temp_around_fermi,
     &                         9999d0,0)
         if(temp_around_fermi > 1e-8) then
           around_fermi = around_fermi * kb * temp_around_fermi
         endif
         if(around_fermi > 1e-8) then
           emin = efermi - around_fermi
           emax = efermi + around_fermi
         endif
c end-mod-xjq
         CALL SECTION_SET_REAL('IME',EIMAG,9999D0,0)
         IF ( ISEARCHEF.EQ.0 ) CALL SECTION_FIND_KEYWORD('SEARCHEF',
     &        SEARCHEF)
         IF (IGRID(1).EQ.11) THEN
            CALL SECTION_SET_INTEGER('NEPOL',NEPOL,5,0)
            CALL SECTION_SET_INTEGER('NEFD1',NEFD1,3,0)
            CALL SECTION_SET_INTEGER('NEFD2',NEFD2,20,0)
            CALL SECTION_SET_INTEGER('NEFD3',NEFD3,10,0)
            IF ( NEFD3 .GT. 16 ) THEN
               STOP 'GAUSS-FD FOR NEFD3 > 16 IS NOT AVAILABLE NOW!'
            ENDIF
            CALL SECTION_FIND_KEYWORD('ACTIVE_CONTN',lactive_contn)
            if(lactive_contn) then
              lepath_contn = .true.
            else
              CALL SECTION_FIND_KEYWORD('EPATH_CONTN',lepath_contn)
            endif
            if(lepath_contn) then
              necontn = netab(1) -nefd1 -nefd2 -nefd3 -nepol
              CALL SECTION_SET_REAL('EMIN_CONTN',EMIN_contn,emin,0)
              CALL SECTION_SET_REAL('EMAX_CONTN',emax_contn,emax,0)
              CALL SECTION_SET_REAL('IME_CONTN',ime_contn,0.0d0,0)
              CALL SECTION_SET_INTEGER('NEPOL_CONTN',nepol_contn,0,0)
              if(necontn < 2 .or. necontn < nepol_contn) 
     &          stop 'too few necontn'
              if(necontn-nepol_contn==1) stop 'necontn-nepol_contn==1'
            else
              nefd2 = netab(1) -nefd1 -nefd3 -nepol
            endif
          END IF
      END IF
      IF ( EMAX.LT.-9900D0 ) EMAX = EFERMI
C
C---------------------------------------------------- conductivity SIGMA
      IF ( TASK(1:5).EQ.'SIGMA' ) THEN
         IF ( NETAB(1).GT.1 ) THEN
            IF ( ABS(EFERMI-EMIN).LT.1D-3 ) EMIN = EMINDEFAULT
         END IF
      END IF
C
C------------------------------- allocate variables depending only on NE
      NEMAX = MAX(NETAB(1),NETAB(2))
      IF ( TASK(1:5).EQ.'SIGMA' .OR. TASK(1:7).EQ.'COMPTON' )
     &     NEMAX = NETAB(1) + NETAB(2)
      IF ( TASK(1:5).EQ.'EKREL' ) NEMAX = 1
      IF ( NEMAX.EQ.0 ) STOP '<MAIN>    NEMAX = 0 !!!!!!!! '
      ALLOCATE (ETAB(NEMAX,2),WETAB(NEMAX,2),STAT=IA_ERR)
      ETAB = 999999D0
      WETAB = 999999D0
      IF ( IA_ERR.NE.0 ) STOP 'alloc: main -> ETAB'
C
C---------------------------------------------------- test Lloyd-formula
      IF ( ITEST.EQ.5 ) IGRID(1) = 3
C
      IF ( (IGRID(1).EQ.1) .OR. (IGRID(1).EQ.3) .OR. (IGRID(1).EQ.6) )
     &     SEARCHEF = .FALSE.
C
      IF ( TASK(1:4).NE.'NONE' .AND. TASK(1:4).NE.'DOS ' )
     &     SEARCHEF = .FALSE.
C
      IF ( FOUND_SECTION ) CALL SECTION_SET_REAL('EF',EFERMI,9999D0,0)
C
      IF ( SCFSTATUS(1:5).EQ.'START' ) THEN
         EFERMI = 1D0
         EMAX = EFERMI
      ELSE IF ( (EFERMI.LT.-9000.D0) .AND. (TASK(1:4).NE.'NONE') ) THEN
         STOP 'specify Fermi energy '
      END IF
C
      IF ( TASK(1:5).EQ.'XAS  ' .OR. TASK(1:5).EQ.'XMO  ' .OR. TASK(1:5)
     &     .EQ.'XRS  ' .OR. TASK(1:5).EQ.'APS  ' .OR. TASK(1:5)
     &     .EQ.'RELAX' .OR. TASK(1:5).EQ.'BIS  ' .OR. TASK(1:5)
     &     .EQ.'T1   ' ) EMIN = EFERMI
C
      IF ( TASK(1:5).EQ.'AES  ' .OR. TASK(1:5).EQ.'NRAES' .OR. TASK(1:5)
     &     .EQ.'XES  ' .OR. TASK(1:5).EQ.'T1   ' .OR. TASK(1:7)
     &     .EQ.'COMPTON' .OR. TASK(1:2).EQ.'TL' ) EMAX = EFERMI
C
      IF ( TASK(1:5).EQ.'VBPES' .OR. TASK(1:5).EQ.'ARPES' ) THEN
         IF ( EMAX.LT.-9900D0 ) EMAX = EFERMI
         CALL INPUT_FIND_SECTION('ENERGY',0)
         IF ( FOUND_SECTION ) THEN
            CALL SECTION_SET_REAL('ERANGE',ERANGE,9999D0,0)
            IF ( FOUND_REAL ) EMIN = EFERMI - ERANGE/RY_EV
         END IF
      END IF
C
      IF ( TASK(1:3).EQ.'T1 ' ) THEN
         EIMAG = 0.0D0
         CALL INPUT_FIND_SECTION('ENERGY',0)
         IF ( FOUND_SECTION ) CALL SECTION_SET_REAL('IME',EIMAG,9999D0,
     &        0)
         IGRID(1) = 3
         NETAB(1) = 1
         ETAB(1,1) = DCMPLX(EFERMI,EIMAG)
         EMIN = EFERMI
         EMAX = EFERMI
      END IF
C
      IF ( TASK(1:3).EQ.'BSF' .AND. EIMAG.LT.1D-6 .AND. NCPA.EQ.0 )
     &     CALL STOP_MESSAGE(ROUTINE,
     &       'for TASK=BSF: Im(E)<10^-6  not allowed for ordered system'
     &       )
C
C-------------------------------- check consistency of supplied TAU-file
C
      IF ( RDTAU .OR. RDTAUMQ ) THEN
         IF ( NETAU.NE.NETAB(1) ) CALL STOP_MESSAGE(ROUTINE,
     &      'the number of E-points in TAU-file differs from input file'
     &      )
      END IF
C
      IF ( TASK(1:3).EQ.'SCF' .OR. TASK(1:5).EQ.'APS  ' .OR. TASK(1:5)
     &     .EQ.'AES  ' .OR. TASK(1:5).EQ.'NRAES' .OR. TASK(1:5)
     &     .EQ.'RELAX' .OR. TASK(1:6).EQ.'SOCPAR' .OR. TASK(1:6)
     &     .EQ.'PSHIFT' .OR. TASK(1:5).EQ.'FSCAT' .OR. TASK(1:7)
     &     .EQ.'MECHECK' .OR. TASK(1:6).EQ.'MEPLOT' .OR. TASK(1:6)
     &     .EQ.'WFPLOT' ) NOWRDOS = .TRUE.
C
C ======================================================================
      IF ( INITELOOP .OR. TASK(1:1).EQ.'X' .OR. TASK(1:6)
     &     .EQ.'SOCPAR' .OR. TASK(1:6).EQ.'PSHIFT' .OR. TASK(1:5)
     &     .EQ.'FSCAT' .OR. TASK(1:6).EQ.'MEPLOT' .OR. TASK(1:6)
     &     .EQ.'WFPLOT' ) THEN
C
         IF ( TASK(1:6).EQ.'SOCPAR' .OR. TASK(1:6).EQ.'PSHIFT' .OR. 
     &        TASK(1:5).EQ.'FSCAT' ) THEN
            SPLITSS = .FALSE.
            EIMAG = 0D0
            EILOW = 0D0
            SEARCHEF = .FALSE.
         END IF
C
         WRITE (6,99001) PATH(IGRID(1))
C
         IF ( NETAB(1).GT.NEMAX .OR. NETAB(2).GT.NEMAX ) THEN
            WRITE (6,99002) PROGNAME,NETAB,NEMAX
            NETAB(1) = MIN(NEMAX,NETAB(1))
            NETAB(2) = MIN(NEMAX,NETAB(2))
         END IF
C
         IF ( TASK(1:3).EQ.'SCF' ) THEN
            IF ( (IGRID(1).EQ.2) .OR. (IGRID(1).EQ.4) .OR. 
     &           (IGRID(1).EQ.5) ) EMAX = EFERMI
         END IF
C
         IF ( SPLITSS ) EIMAG = 0D0
C
         WRITE (6,99003) 'energies  :  ',NETAB(1)
         WRITE (6,99004) 'E_min     :  ',EMIN
         IF ( SCFSTATUS(1:5).NE.'START' ) THEN
            WRITE (6,99004) 'E_max     :  ',EMAX
            IF ( IGRID(1).NE.5 .OR. IGRID(1).NE.8 ) WRITE (6,99004)
     &            'Im(E)     :  ',EIMAG
            IF ( SEARCHEF ) WRITE (6,99004) 
     &                         'E_Fermi will be searched starting with:'
            WRITE (6,99004) 'E_Fermi   :  ',EFERMI
         ELSE
            WRITE (6,99004) 'E_Fermi will be searched'
         END IF
         WRITE (6,*) ' '
C
         IF ( SPLITSS ) THEN
            IF ( IGRID(1).NE.5 .AND. IGRID(2).NE.3 ) THEN
               WRITE (6,99005) PROGNAME,IGRID
               STOP
            END IF
            WRITE (6,*) ' '
            WRITE (6,99004) 'SS-contributions will be split using:'
            WRITE (6,99003) 'energies  :  ',NETAB(2)
            WRITE (6,99004) 'Im(E)     :  ',EIMAG
            WRITE (6,99003) 'path      :  '//PATH(IGRID(2))
            NEPATH = 2
         ELSE
            NEPATH = 1
         END IF
C
         DO IP = 1,NEPATH
            CALL EPATH(IGRID(IP),EMIN,EMAX,EIMAG,NETAB(IP),ETAB(1,IP),
     &                 WETAB(1,IP),EILOW,IPRINT,NEMAX)
         END DO
C
C
         IMEMIN = 0.0D0
         IMEMAX = 0.0D0
         DO IE = 1,NETAB(1)
            IMEMIN = MIN(IMEMIN,DIMAG(ETAB(IE,1)))
            IMEMAX = MAX(IMEMAX,DIMAG(ETAB(IE,1)))
         END DO
C
C ----------------------------------------------------------------------
      END IF
C ======================================================================
C               variables connected with Lloyd formula
C ======================================================================
C
      ALLOCATE (PHAST(NKMMAX,NTMAX,NEMAX),STAT=IA_ERR)
      ALLOCATE (PHASA(NKMMAX,NTMAX,NEMAX),PHASK(NEMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: PHASA')
C
      ALLOCATE (GFEP(3,NGFEPMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: GFEP')
C
      NGFEP = 0
      GFEP(:,:) = 999999D0
C
99001 FORMAT (/,1X,79('e'),//,10X,'energy path:  ',A40)
99002 FORMAT (2(/,1X,79('#')),/,10X,'WARNING from   ',A,/,10X,'NETAB =',
     &        2I6,' too high ----  reduced to  NEMAX =',I6,
     &        2(/,1X,79('#'))/)
99003 FORMAT (10X,A,5I10)
99004 FORMAT (10X,A,5F13.6)
99005 FORMAT (2(/,1X,79('#')),/,10X,'MESSAGE from   ',A,/,10X,
     &        'for splitted energy path   (SPLITSS is set):',/,10X,
     &        'EPATH = { 5, 3}  expected in the input file ',/,10X,
     &        'EPATH = {',I2,',',I2,'}  supplied by user instead ',
     &        2(/,1X,79('#'))/)
C
      END
