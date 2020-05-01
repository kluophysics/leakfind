C*==tutorial.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TUTORIAL(RUNELOOP,TASK,DOSREP)
C   ********************************************************************
C   *                                                                  *
C   *  main subroutine program for       KKRGEN                        *
C   *  adapted as tutorial on   multiple scattering theory             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:MEZJ,MEZZ,MSST,SSST,TSST
      USE MOD_FILES,ONLY:PLOTPRS,IPRINT,FOUND_SECTION,FOUND_REAL,
     &    FOUND_INTEGER
      IMPLICIT NONE
C*--TUTORIAL14
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*3 DOSREP
      LOGICAL RUNELOOP
      CHARACTER*10 TASK
C
C Local variables
C
      INTEGER L,LMAX,M,NX
      CHARACTER*10 NORMALISATION,TUT_TASK,X_RANGE
      REAL*8 XMAX,XMIN
C
C*** End of declarations rewritten by SPAG
C
C=======================================================================
C=======================================================================
C
      CALL INPUT_FIND_SECTION('TASK',0)
      IF ( FOUND_SECTION ) CALL SECTION_SET_STRING('TUT-TASK',TUT_TASK,
     &     '9999',0)
C
      WRITE (6,99002) TUT_TASK
C
C=======================================================================
C
      IF ( TUT_TASK(1:3).EQ.'YLM' ) THEN
C
         L = 2
         M = 0
         CALL SECTION_SET_INTEGER('L',L,9999,0)
         IF ( .NOT.FOUND_INTEGER ) WRITE (6,99001) 'L'
         CALL SECTION_SET_INTEGER('M',M,9999,0)
         M = SIGN(MIN(L,ABS(M)),M)
C
         CALL TUTORIAL_YLM(L,M)
C
C=======================================================================
C
      ELSE IF ( TUT_TASK(1:10).EQ.'JLX-ASYMPT' ) THEN
C
         LMAX = 2
         NX = 300
         XMIN = 0.001D0
         XMAX = 5D0
C
         CALL SECTION_SET_INTEGER('LMAX',LMAX,9999,0)
         IF ( .NOT.FOUND_INTEGER ) WRITE (6,99001) 'L'
         LMAX = MIN(LMAX,4)
         CALL SECTION_SET_INTEGER('NX',M,9999,0)
         IF ( .NOT.FOUND_INTEGER ) WRITE (6,99001) 'NX'
         CALL SECTION_SET_REAL('XMIN',XMIN,9999D0,0)
         IF ( .NOT.FOUND_REAL ) WRITE (6,99001) 'XMIN'
         IF ( XMIN.LT.0.001D0 ) XMIN = 0.001D0
         CALL SECTION_SET_REAL('XMAX',XMAX,9999D0,0)
         IF ( .NOT.FOUND_REAL ) WRITE (6,99001) 'XMAX'
         IF ( XMAX.LT.XMIN ) XMAX = XMIN + 1D0
C
         CALL SECTION_SET_STRING('X-RANGE',X_RANGE,'SMALL',0)
         IF ( X_RANGE(1:5).EQ.'SMALL' ) XMIN = 0.001D0
C
         CALL TUTORIAL_JLX(LMAX,NX,XMIN,XMAX,X_RANGE)
C
C=======================================================================
C
      ELSE IF ( TUT_TASK(1:3).EQ.'JLX' ) THEN
C
         LMAX = 2
         NX = 300
         XMIN = 0.001D0
         XMAX = 5D0
C
         CALL SECTION_SET_INTEGER('LMAX',LMAX,9999,0)
         IF ( .NOT.FOUND_INTEGER ) WRITE (6,99001) 'L'
         LMAX = MIN(LMAX,20)
         CALL SECTION_SET_INTEGER('NX',M,9999,0)
         IF ( .NOT.FOUND_INTEGER ) WRITE (6,99001) 'NX'
         CALL SECTION_SET_REAL('XMIN',XMIN,9999D0,0)
         IF ( .NOT.FOUND_REAL ) WRITE (6,99001) 'XMIN'
         IF ( XMIN.LT.0.001D0 ) XMIN = 0.001D0
         CALL SECTION_SET_REAL('XMAX',XMAX,9999D0,0)
         IF ( .NOT.FOUND_REAL ) WRITE (6,99001) 'XMAX'
         IF ( XMAX.LT.XMIN ) XMAX = XMIN + 1D0
C
         X_RANGE = 'ALL'
C
         CALL TUTORIAL_JLX(LMAX,NX,XMIN,XMAX,X_RANGE)
C
C=======================================================================
C
      ELSE IF ( TUT_TASK(1:7).EQ.'CHR-POT' ) THEN
C
         PLOTPRS = .TRUE.
C
         CALL FPPLOTPRS
C
C=======================================================================
C
      ELSE IF ( TUT_TASK(1:10).EQ.'POT-DECOMP' ) THEN
C
         CALL TUTORIAL_POT_DECOMP
C
C=======================================================================
C
      ELSE IF ( TUT_TASK(1:6).EQ.'POT-XC' ) THEN
C
         CALL TUTORIAL_POT_XC
C
C=======================================================================
C
      ELSE IF ( TUT_TASK(1:6).EQ.'RWF-L' ) THEN
C
         CALL SECTION_SET_STRING('NORMALISATION',NORMALISATION,'1',0)
C
         CALL TUTORIAL_RWF('L',NORMALISATION,L,TSST,MSST,IPRINT,SSST,
     &                     MEZZ,MEZJ)
C
C=======================================================================
C
      ELSE IF ( TUT_TASK(1:6).EQ.'RWF-E' ) THEN
C
         CALL SECTION_SET_STRING('NORMALISATION',NORMALISATION,'1',0)
C
         M = 0
         CALL SECTION_SET_INTEGER('L',L,2,0)
C
         CALL TUTORIAL_RWF('E',NORMALISATION,L,TSST,MSST,IPRINT,SSST,
     &                     MEZZ,MEZJ)
C
C=======================================================================
C
      ELSE IF ( TUT_TASK(1:8).EQ.'CORE-ALL' ) THEN
C
         CALL TUTORIAL_CORE_ALL
C
C=======================================================================
C
      ELSE IF ( TUT_TASK(1:3).EQ.'GEN' ) THEN
C
         CALL GEN(RUNELOOP,TASK,DOSREP,'LOC')
C
      END IF
C
      STOP
C=======================================================================
C
99001 FORMAT (5x,'variable ',A,' not supplied -- default will be used')
99002 FORMAT (' '//,10X,63('*'),/,10X,'*',61X,'*',/,10X,
     & '* *******  *    *  *******   ****   *****   *    **    *      *'
     & ,/,10X,
     & '*    *     *    *     *     *    *  *    *  *   *  *   *      *'
     & ,/,10X,
     & '*    *     *    *     *     *    *  *    *  *  *    *  *      *'
     & ,/,10X,
     & '*    *     *    *     *     *    *  *****   *  ******  *      *'
     & ,/,10X,
     & '*    *     *    *     *     *    *  *  *    *  *    *  *      *'
     & ,/,10X,
     & '*    *     *    *     *     *    *  *   *   *  *    *  *      *'
     & ,/,10X,
     & '*    *      ****      *      ****   *    *  *  *    *  ****** *'
     & ,/,10X,'*',61X,'*',/,10X,63('*'),//,10X,'TUT-TASK = ',A,//)
      END
C*==tutorial_ylm.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TUTORIAL_YLM(L,M)
C   ********************************************************************
C   *                                                                  *
C   *  TUTORIAL  YLM   plot real sperical harmonics                    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:PI
      USE MOD_FILES,ONLY:DATSET,LDATSET
      IMPLICIT NONE
C*--TUTORIAL_YLM202
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER L,M
C
C Local variables
C
      REAL*8 CP,CT,PHI,RHAT(3),RYLM(:),SP,ST,TET,YVEC(3)
      CHARACTER*80 FILNAM
      INTEGER IPHI,ITET,L1,L2,LFILNAM,LM,NLMTUT,NLTUT,NPHI,NTET
      CHARACTER*5 STRL,STRM
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE RYLM
C
      NLTUT = L + 1
      NLMTUT = NLTUT*NLTUT
C
      ALLOCATE (RYLM(NLMTUT))
C
      OPEN (1,FILE='zzzzzz_ypls.dat')
      OPEN (2,FILE='zzzzzz_ymin.dat')
C
      NTET = 90
      NPHI = 180
C
      DO ITET = 1,NTET
C
         TET = PI*(ITET-1)/FLOAT(NTET-1)
         ST = SIN(TET)
         CT = COS(TET)
C
         DO IPHI = 1,NPHI
C
            PHI = 2*PI*(IPHI-1)/FLOAT(NPHI-1)
            SP = SIN(PHI)
            CP = COS(PHI)
C
            RHAT(1) = ST*CP
            RHAT(2) = ST*SP
            RHAT(3) = CT
C
            CALL CALC_RHPLM(RHAT(1),RHAT(2),RHAT(3),RYLM,NLTUT-1,NLMTUT)
C
            LM = L*L + L + M + 1
C
            YVEC(1:3) = RHAT(1:3)*ABS(RYLM(LM))
C
            IF ( RYLM(LM).GT.0 ) THEN
C
               WRITE (1,99001) YVEC
C
            ELSE
C
               WRITE (2,99001) YVEC
C
            END IF
C
         END DO
C
         WRITE (1,*) ' '
         WRITE (2,*) ' '
      END DO
      CLOSE (1)
      CLOSE (2)
C
      WRITE (STRL,'(I5)') L
      CALL STRING_TRIM_LEFT(STRL)
      L1 = LEN_TRIM(STRL)
      WRITE (STRM,'(I5)') M
      CALL STRING_TRIM_LEFT(STRM)
      L2 = LEN_TRIM(STRM)
C
      FILNAM = DATSET(1:LDATSET)//'.gnp'
      LFILNAM = LDATSET + 4
      OPEN (1,FILE=FILNAM(1:LFILNAM))
      WRITE (1,99002) STRL(1:L1),STRM(1:L2)
C
      WRITE (6,99003) STRL(1:L1),STRM(1:L2),FILNAM(1:LFILNAM)
C
99001 FORMAT (3F12.5)
99002 FORMAT ('set title "SPR-KKR tutorial: real spherical harmonics',
     &        ' for l = ',A,'  m = ',A,' "',/,'unset key',/,
     &        'set size square',/,'set ticslevel -0.5',/,'set size 1,1',
     &        /,
     &     '#----------------------------------------------------------'
     &     ,/,'set xlabel "x"',/,'set ylabel "y"',/,'set zlabel "z"',/,
     &     '#unset border ; unset xlabel ; unset ylabel ; unset zlabel',
     &     /,
     &     '                unset xlabel ; unset ylabel ; unset zlabel',
     &     /,'set xrange [-1.0:1.0] ',/,'set yrange [-1.0:1.0] ',/,
     &     'set zrange [-1.0:1.0] ',/,
     &     'set arrow 1 from -1.0, 0, 0 to 1.0, 0, 0',
     &     '  back nofilled linetype -1 linewidth 1.000',/,
     &     'set arrow 2 from 0, -1.0, 0 to 0, 1.0, 0',
     &     '  back nofilled linetype -1 linewidth 1.000',/,
     &     'set arrow 3 from 0, 0, -1.0 to 0, 0, 1.0',
     &     '  back nofilled linetype -1 linewidth 1.000',/,
     &     'set label 1 "x" at 1.0, 0.05, 0.05 centre norotate  ',/,
     &     'set label 2 "y" at 0.05, 1.0, 0.05 centre norotate  ',/,
     &     'set label 3 "z" at 0.05, 0.05, 1.0 centre norotate  ',/,
     &     'set parametric',/,'set view 70, 20, 1, 1',/,
     &     'set samples 51, 51',/,'set isosamples 20, 20',/,
     &     'set style data lines',/,'#set hidden3d',/,
     &     'splot  "zzzzzz_ypls.dat" with lines linestyle 1,',
     &     '       "zzzzzz_ymin.dat" with lines linestyle 3',/,
     &     'pause -1',/)
99003 FORMAT (/,10X,'spherical harmonics tabulated for  l = ',A,
     &        '   m = ',A,//,10X,'plot calling      gnuplot ',A,/)
      END
C*==tutorial_jlx.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TUTORIAL_JLX(LMAX,NX,XMIN,XMAX,X_RANGE)
C   ********************************************************************
C   *                                                                  *
C   *  TUTORIAL  JLX    plot spherical Bessel and Neumann functions    *
C   *                                                                  *
C   ********************************************************************
      USE MOD_FILES,ONLY:LSYSTEM,SYSTEM,LDATSET,DATSET,IOTMP
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C*--TUTORIAL_JLX337
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NCOLORTAB
      PARAMETER (NCOLORTAB=14)
C
C Dummy arguments
C
      INTEGER LMAX,NX
      REAL*8 XMAX,XMIN
      CHARACTER*10 X_RANGE
C
C Local variables
C
      COMPLEX*16 ARG
      COMPLEX*16 CJLZ,CNLZ
      INTEGER COLORTAB(NCOLORTAB),I,IC,ICOLOR_C(:),ISTYLE_C(:),IX,L,
     &        LFILNAM,LHEADER1,LHEADER2,LSYS,LYTXT1,LYTXT2,NC,TLM1,TLP1
      REAL*8 DBLFACT(-1:100),DX,JMAX,JMIN,JTAB(:,:),JTAB_ASYM(:,:),NMAX,
     &       NMIN,NTAB(:,:),NTAB_ASYM(:,:),PI_HALF,X,XPWL,XTAB(:)
      CHARACTER*80 FILNAM,HEADER1,HEADER2,SYS,YTXT1,YTXT2
      CHARACTER*20 LEG(:)
C
C*** End of declarations rewritten by SPAG
C
      DATA COLORTAB/1,2,4,3,6,7,8,9,10,11,12,13,14,15/
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE JTAB,NTAB,XTAB,LEG,JTAB_ASYM,NTAB_ASYM
      ALLOCATABLE ICOLOR_C,ISTYLE_C
C
      IF ( X_RANGE(1:3).EQ.'ALL' ) THEN
         WRITE (6,99002) LMAX,NX,XMIN,XMAX
      ELSE
         WRITE (6,99002) LMAX,NX,XMIN,XMAX,X_RANGE
      END IF
C
      PI_HALF = PI/2D0
C
      DBLFACT(-1:100) = 0
      DBLFACT(-1) = 1D0
      DBLFACT(1) = 1D0
      DO I = 3,100,2
         DBLFACT(I) = DBLFACT(I-2)*I
      END DO
C
      SYS = SYSTEM
      LSYS = LSYSTEM
      CALL XMGRSUBSCRIPTS(SYS,LSYS,80)
C
      NC = LMAX + 1
      ALLOCATE (JTAB(NX,0:LMAX),JTAB_ASYM(NX,0:LMAX),XTAB(NX))
      ALLOCATE (NTAB(NX,0:LMAX),NTAB_ASYM(NX,0:LMAX),LEG(100))
      ALLOCATE (ICOLOR_C(2*NC),ISTYLE_C(2*NC))
C
      DX = (XMAX-XMIN)/DBLE(NX-1)
      JMIN = 0D0
      JMAX = 0D0
      NMIN = 0D0
      NMAX = 0D0
C
      DO IX = 1,NX
C
         X = XMIN + DX*(IX-1)
         ARG = DCMPLX(X,0D0)
C
         XTAB(IX) = X
C
         DO L = 0,LMAX
C
            JTAB(IX,L) = DREAL(CJLZ(L,ARG))
            NTAB(IX,L) = DREAL(CNLZ(L,ARG))
C
            JMIN = MIN(JMIN,JTAB(IX,L))
            JMAX = MAX(JMAX,JTAB(IX,L))
            NMIN = MIN(NMIN,NTAB(IX,L))
            NMAX = MAX(NMAX,NTAB(IX,L))
         END DO
C
      END DO
      JMIN = -0.4D0
      JMAX = +1.1D0
      NMIN = -1.1D0
      NMAX = +0.4D0
C
      IF ( X_RANGE(1:5).EQ.'SMALL' ) THEN
C
         DO IX = 1,NX
C
            X = XTAB(IX)
            XPWL = 1D0
C
            TLP1 = +1
            TLM1 = -1
            DO L = 0,LMAX
C
               JTAB_ASYM(IX,L) = XPWL/DBLFACT(TLP1)
               NTAB_ASYM(IX,L) = -DBLFACT(TLM1)/(XPWL*X)
C
               TLP1 = TLP1 + 2
               TLM1 = TLM1 + 2
               XPWL = XPWL*X
            END DO
C
         END DO
C
      ELSE IF ( X_RANGE(1:5).EQ.'LARGE' ) THEN
C
C
         DO IX = 1,NX
C
            X = XTAB(IX)
C
            DO L = 0,LMAX
               JTAB_ASYM(IX,L) = +SIN(X-L*PI_HALF)/X
               NTAB_ASYM(IX,L) = -COS(X-L*PI_HALF)/X
            END DO
         END DO
C
      END IF
C
C=======================================================================
C                      set up xmgrace file
C=======================================================================
C
      YTXT1 = 'j!sl!N(x)'
      LYTXT1 = 9
      YTXT2 = 'n!sl!N(x)'
      LYTXT2 = 9
C
      HEADER1 = 'SPR-KKR tutorial'
      LHEADER1 = 16
C
      IF ( X_RANGE(1:5).EQ.'SMALL' ) THEN
C
         HEADER2 = 
     &      'asymptotic behavior of j!sl!N(x) and n!sl!N(x) for  x -> 0'
         LHEADER2 = 59
C
      ELSE IF ( X_RANGE(1:6).EQ.'LARGE' ) THEN
C
         HEADER2 = 'asymptotic behavior of j!sl!N(x) and n!sl!N(x)'//
     &             ' for  x -> infinity'
         LHEADER2 = 66
C
      ELSE
C
         HEADER2 = 
     &  'spherical Bessel j!sl!N(x) and von Neumann n!sl!N(x) functions'
         LHEADER2 = 63
C
      END IF
C
      CALL XMGRHEAD(DATSET,LDATSET,'sphfun',6,'XXXX',0,FILNAM,80,
     &              LFILNAM,IOTMP,2,XMIN,0,XMAX,0,JMIN,0,JMAX,0,NMIN,0,
     &              NMAX,0,'x',1,YTXT1,LYTXT1,YTXT2,LYTXT2,HEADER1,
     &              LHEADER1,HEADER2,LHEADER2,.FALSE.)
C
      LEG(1) = 's'
      LEG(2) = 'p'
      LEG(3) = 'd'
      DO IC = 4,NC
         LEG(IC) = CHAR(ICHAR('f')+IC-4)
      END DO
C
      CALL XMGRLEGEND(IOTMP,2,NC,NC,LEG,LEG)
C
      CALL XMGRCURVES(IOTMP,2,2*NC,2*NC,3,0,0)
C
      IF ( NC.LE.NCOLORTAB ) THEN
         ICOLOR_C(1:NC) = COLORTAB(1:NC)
      ELSE
         ICOLOR_C(1:NCOLORTAB) = COLORTAB(1:NCOLORTAB)
         ICOLOR_C((NCOLORTAB+1):NC) = 1
      END IF
      ISTYLE_C(1:NC) = 1
C
      CALL XMGRCOLOR(IOTMP,0,NC,ICOLOR_C)
      CALL XMGRCOLOR(IOTMP,1,NC,ICOLOR_C)
      CALL XMGRSTYLE(IOTMP,0,NC,ISTYLE_C)
      CALL XMGRSTYLE(IOTMP,1,NC,ISTYLE_C)
C
      DO L = 0,LMAX
         CALL XMGRTABLE(0,L,XTAB,JTAB(1,L),1D0,NX,IOTMP)
      END DO
C
      DO L = 0,LMAX
         CALL XMGRTABLE(1,L,XTAB,NTAB(1,L),1D0,NX,IOTMP)
      END DO
C
C=======================================================================
      IF ( X_RANGE(1:3).NE.'ALL' ) THEN
C
         IF ( X_RANGE(1:5).EQ.'SMALL' ) THEN
C
            CALL XMGRTEXT(IOTMP,'j\sl\N(x) --> +x\Sl\N / (2l+1)!!!!',34,
     &                    0.45D0,0.35D0,'VIEW ')
C
            CALL XMGRTEXT(IOTMP,'n\sl\N(x) --> -(2l-1)!!!! / x\Sl+1',34,
     &                    0.45D0,0.65D0,'VIEW ')
C
         ELSE
C
            CALL XMGRTEXT(IOTMP,'j\sl\N(x) --> + sin(x-l!xp!0/2) / x',
     &                    35,0.35D0,0.45D0,'VIEW ')
C
            CALL XMGRTEXT(IOTMP,'n\sl\N(x) --> - cos(x-l!xp!0/2) / x',
     &                    35,0.35D0,0.55D0,'VIEW ')
C
         END IF
C
         IF ( NC.LE.NCOLORTAB ) THEN
            ICOLOR_C((NC+1):(NC+NC)) = COLORTAB(1:NC)
         ELSE
            ICOLOR_C((NC+1):(NC+NCOLORTAB)) = COLORTAB(1:NCOLORTAB)
            ICOLOR_C((NC+NCOLORTAB+1):(NC+NC)) = 1
         END IF
         ISTYLE_C((NC+1):(NC+NC)) = 3
C
         CALL XMGRCOLOR(IOTMP,0,2*NC,ICOLOR_C)
         CALL XMGRCOLOR(IOTMP,1,2*NC,ICOLOR_C)
         CALL XMGRSTYLE(IOTMP,0,2*NC,ISTYLE_C)
         CALL XMGRSTYLE(IOTMP,1,2*NC,ISTYLE_C)
C
         DO L = 0,LMAX
            CALL XMGRTABLE(0,NC+L,XTAB,JTAB_ASYM(1,L),1D0,NX,IOTMP)
         END DO
C
         DO L = 0,LMAX
            CALL XMGRTABLE(1,NC+L,XTAB,NTAB_ASYM(1,L),1D0,NX,IOTMP)
         END DO
C
      END IF
C=======================================================================
C
      WRITE (6,*) '  '
      WRITE (6,99001) 'functions ',FILNAM(1:LFILNAM)
      CLOSE (IOTMP)
C
99001 FORMAT (/,10X,A,' written to the file ',A,/)
99002 FORMAT (/,10X,
     &        'plotting spherical Bessel and von Neumann functions for',
     &        //,10X,'l_max =',I2,5X,'NX =',I4,5X,'X = ',F8.4,' ... ',
     &        F8.4,/,:,/,10X,
     &        'asmyptotic behaviour of j_l(x) and n_l(x) for  ',A,' x',
     &        /)
      END
C*==tutorial_pot_decomp.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TUTORIAL_POT_DECOMP
C   ********************************************************************
C   *                                                                  *
C   *  TUTORIAL  POT.DECOMP    decompose potential                     *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CALCMODE,ONLY:NONMAG,TUTORIAL
      USE MOD_FILES,ONLY:IOTMP,LDATSET0,DATSET0,LSYSTEM,SYSTEM
      USE MOD_SCF,ONLY:SCFVXC
      USE MOD_RMESH,ONLY:R,JRWS,NRMAX
      USE MOD_TYPES,ONLY:ITBOT,ITTOP,IMT,RHOCHR,RHOSPN,LTXT_T,TXT_T,VT,
     &    BT,Z
      USE MOD_RMESH,ONLY:FULLPOT
      IMPLICIT NONE
C*--TUTORIAL_POT_DECOMP617
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NCOLORTAB,NC
      PARAMETER (NCOLORTAB=14,NC=5)
C
C Local variables
C
      REAL*8 BTAB(:,:),BXCMAX,BXCMIN,BXCP,EXC(:),RHO4PI(:,:),VAUX(:,:),
     &       VTAB(:,:),VXAUX(:,:),VXCAUX(:,:),VXCMAX,VXCMIN,VXCP,
     &       WEXC(:,:),XMAX,XMIN
      CHARACTER*80 BXCTXT,FILNAM,HEADER1,HEADER2,VXCTXT
      INTEGER COLORTAB(NCOLORTAB),IC,ICOLOR_C(:),IGRAPH,IM,IR,IRBOT,
     &        IRTOP,ISTYLE_C(:),IT,LBXCTXT,LFILNAM,LHEADER1,LHEADER2,LS,
     &        LTXT_IT,LVXCTXT,NGRAPH,NX
      CHARACTER*20 LEG(NC),STR20
C
C*** End of declarations rewritten by SPAG
C
      DATA COLORTAB/1,2,4,3,6,7,8,9,10,11,12,13,14,15/
      DATA LEG/'V!stot','V!snuc','V!sH','V!sxc','V!sx'/
C
      ALLOCATABLE RHO4PI,WEXC,EXC,VTAB,BTAB
      ALLOCATABLE VXCAUX,VXAUX,VAUX,ICOLOR_C,ISTYLE_C
C
      ALLOCATE (RHO4PI(NRMAX,2),WEXC(NRMAX,2),EXC(NRMAX))
      ALLOCATE (VXCAUX(NRMAX,2),VXAUX(NRMAX,2),VAUX(NRMAX,2))
      ALLOCATE (VTAB(NRMAX,NC),BTAB(NRMAX,NC))
      ALLOCATE (ICOLOR_C(2*NC),ISTYLE_C(2*NC))
C
      IF ( FULLPOT ) RETURN
C
C-----------------------------------------------------------------------
C                      loop over types IT
C-----------------------------------------------------------------------
      DO IT = ITBOT,ITTOP
C
         IM = IMT(IT)
         IRBOT = 1
         IRTOP = JRWS(IM)
C
         VXCMIN = 0D0
         VXCMAX = 0D0
         BXCMIN = 0D0
         BXCMAX = 0D0
C
         DO IR = 1,IRTOP
            VXCAUX(IR,1) = 0D0
            VXCAUX(IR,2) = 0D0
            VXAUX(IR,1) = 0D0
            VXAUX(IR,2) = 0D0
            VAUX(IR,1) = 0D0
            VAUX(IR,2) = 0D0
C
            RHO4PI(IR,1) = RHOCHR(IR,IT)
            RHO4PI(IR,2) = RHOSPN(IR,IT)
C
            VTAB(IR,1) = VT(IR,IT)
            BTAB(IR,1) = BT(IR,IT)
C
            VTAB(IR,2) = -2*Z(IT)/R(IR,IM)
            BTAB(IR,2) = 0D0
C
         END DO
C
         DO IR = 1,IRTOP
            WEXC(IR,1) = 0D0
         END DO
C
C -------------------------------------------------- get MJW X-potential
         CALL EXCMJW(RHO4PI,VXAUX,EXC,WEXC,IRTOP,NRMAX,0)
C
C ------------------------------------------------- get MJW XC-potential
         CALL EXCMJW(RHO4PI,VXCAUX,EXC,WEXC,IRTOP,NRMAX,1)
C
C
C -------------------------------------------- get grevious XC-potential
         IF ( SCFVXC(1:3).EQ.'VBH' ) THEN
C
            CALL EXCVBH(RHO4PI,VAUX,EXC,WEXC,IRTOP,NRMAX)
C
         ELSE IF ( SCFVXC(1:3).EQ.'MJW' ) THEN
C
            CALL EXCMJW(RHO4PI,VAUX,EXC,WEXC,IRTOP,NRMAX,1)
C
         ELSE IF ( SCFVXC(1:3).EQ.'VWN' ) THEN
C
            CALL EXCVWN(RHO4PI,VAUX,EXC,WEXC,IRTOP,NRMAX)
C
         ELSE IF ( SCFVXC(1:3).EQ.'PBE' ) THEN
C
            STOP ' SCFVXC = PBE not allowed'
C
         ELSE
            WRITE (6,99001) SCFVXC
            STOP
         END IF
C
C -------------------- convention VXCAUX(I,IS) IS=1,2 == up, down (AKAI)
C
         DO IR = 1,IRTOP
C
C -------------------------------- V_H = V_tot - V_nuc - V_xc (previous)
C
            VXCP = (VAUX(IR,1)+VAUX(IR,2))/2.0D0
            BXCP = (VAUX(IR,1)-VAUX(IR,2))/2.0D0
C
            VTAB(IR,3) = VT(IR,IT) - VTAB(IR,2) - VXCP
            BTAB(IR,3) = BT(IR,IT) - BTAB(IR,2) - BXCP
C
            VXCMAX = MAX(VXCMAX,VTAB(IR,3))
C
C ----------------------------------------------------- MJW XC-potential
C
            VTAB(IR,4) = (VXCAUX(IR,1)+VXCAUX(IR,2))/2.0D0
            BTAB(IR,4) = (VXCAUX(IR,1)-VXCAUX(IR,2))/2.0D0
C
C -----------------------------------------------------  MJW X-potential
C
            VTAB(IR,5) = (VXAUX(IR,1)+VXAUX(IR,2))/2.0D0
            BTAB(IR,5) = (VXAUX(IR,1)-VXAUX(IR,2))/2.0D0
C
            BXCMIN = MIN(BXCMIN,BT(IR,IT))
            BXCMAX = MAX(BXCMAX,BT(IR,IT))
C
         END DO
C
         IF ( ABS(BXCMIN).LT.1D-4 .AND. ABS(BXCMAX).LT.1D-4 )
     &        NONMAG = .TRUE.
C
         VXCMIN = -VXCMAX
C
C ----------------------------------------------------------------------
C                     plot potentials
C ----------------------------------------------------------------------
C
         STR20 = TXT_T(IT)(1:LTXT_T(IT))//'!N(r) (Ry)'
         LS = LTXT_T(IT) + 10
         VXCTXT = 'V!m{1}!S  !M{1}!s'//STR20(1:LS)
         LVXCTXT = 17 + LS
         IF ( NONMAG ) THEN
            NGRAPH = 1
            BXCTXT = ' '
            LBXCTXT = 1
         ELSE
            NGRAPH = 2
            STR20 = TXT_T(IT)(1:LTXT_T(IT))//'!N(r) (Ry)'
            LS = LTXT_T(IT) + 10
            BXCTXT = 'B!m{1}!S  !M{1}!s'//STR20(1:LS)
            LBXCTXT = 17 + LS
         END IF
C
         NX = IRTOP - IRBOT + 1
C
         XMIN = R(IRBOT,IM)
         XMAX = R(IRTOP,IM)
C
         HEADER2 = 'potential decomposition for '//TXT_T(IT)
     &             (1:LTXT_T(IT))
         LHEADER2 = 28 + LTXT_T(IT)
         IF ( TUTORIAL ) THEN
            HEADER1 = 'SPR-KKR tutorial'
            LHEADER1 = 16
            LTXT_IT = 0
         ELSE
            HEADER1 = 'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM)
            LHEADER1 = 25 + LSYSTEM
            LTXT_IT = LTXT_T(IT)
            HEADER2 = HEADER2(1:LHEADER2)//' in '//SYSTEM(1:LSYSTEM)
            LHEADER2 = LHEADER2 + 4 + LSYSTEM
         END IF
C
         CALL XMGRHEAD(DATSET0,LDATSET0,'Vdecomp',7,TXT_T(IT),LTXT_IT,
     &                 FILNAM,80,LFILNAM,IOTMP,NGRAPH,XMIN,1,XMAX,1,
     &                 VXCMIN,1,VXCMAX,1,BXCMIN,1,BXCMAX,1,
     &                 'radius r (a.u.)',15,VXCTXT,LVXCTXT,BXCTXT,
     &                 LBXCTXT,HEADER1,LHEADER1,HEADER2,LHEADER2,
     &                 .FALSE.)
C
C        CALL XMGRLEGEND(IOTMP,NGRAPH,NVXCTAB,NVXCTAB,LEG,LEG)
         CALL XMGRLEG1(IOTMP,0,NC,LEG,0.6D0,0.45D0)
C
         CALL XMGRCURVES(IOTMP,NGRAPH,NC,NC,3,1,0)
C
         ICOLOR_C(1:NC) = COLORTAB(1:NC)
         ISTYLE_C(1:NC) = 1
C
         DO IGRAPH = 1,NGRAPH
            CALL XMGRCOLOR(IOTMP,(IGRAPH-1),NC,ICOLOR_C)
            CALL XMGRSTYLE(IOTMP,(IGRAPH-1),NC,ISTYLE_C)
         END DO
C
         DO IC = 1,NC
            CALL XMGRTABLE(0,(IC-1),R(IRBOT,IM),VTAB(IRBOT,IC),1D0,NX,
     &                     IOTMP)
         END DO
C
         IF ( .NOT.NONMAG ) THEN
            ICOLOR_C((NC+1):(NC+NC)) = COLORTAB(1:NC)
            ISTYLE_C((NC+1):(NC+NC)) = 3
C
            DO IGRAPH = 1,NGRAPH
               CALL XMGRCOLOR(IOTMP,(IGRAPH-1),2*NC,ICOLOR_C)
               CALL XMGRSTYLE(IOTMP,(IGRAPH-1),2*NC,ISTYLE_C)
            END DO
C
            DO IC = 1,NC
               CALL XMGRTABLE(1,(IC-1),R(IRBOT,IM),BTAB(IRBOT,IC),1D0,
     &                        NX,IOTMP)
            END DO
         END IF
C
         WRITE (6,99002) 'spherical     potential V,   B    ',
     &                   FILNAM(1:LFILNAM)
         CLOSE (IOTMP)
C
      END DO
C
C=======================================================================
C
99001 FORMAT (//,1X,79('#'),/,10X,'ERROR in <TUTORIAL_POT_XC> ',/,10X,
     &        'SCFVXC = ',A,'  not known ',/,1X,79('#'),/)
99002 FORMAT (/,10X,A,' written to file ',A)
      END
C*==tutorial_pot_xc.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TUTORIAL_POT_XC
C   ********************************************************************
C   *                                                                  *
C   *  TUTORIAL  POT-XC    plot various spherical XC potentials        *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CALCMODE,ONLY:NONMAG,TUTORIAL
      USE MOD_FILES,ONLY:IOTMP,LDATSET0,DATSET0,LSYSTEM,SYSTEM
      USE MOD_SCF,ONLY:SCFVXC
      USE MOD_RMESH,ONLY:DRDI,R,JRWS,NRMAX
      USE MOD_TYPES,ONLY:ITBOT,ITTOP,IMT,RHOCHR,RHOSPN,LTXT_T,TXT_T
      USE MOD_RMESH,ONLY:FULLPOT
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C*--TUTORIAL_POT_XC870
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NCOLORTAB,NVXCTAB
      PARAMETER (NCOLORTAB=14,NVXCTAB=4)
C
C Local variables
C
      REAL*8 AGRDRHO(:),AGRDRHOD(:),AGRDRHOU(:),BXC(:,:),BXCMAX,BXCMIN,
     &       EXC(:),GDGAG(:),GDGAGD(:),GDGAGU(:),LAPRHO(:),LAPRHOD(:),
     &       LAPRHOU(:),RHO(:),RHO4PI(:,:),RHOD(:),RHOU(:),VAUX(:,:),
     &       VXC(:,:),VXCMAX,VXCMIN,WEXC(:,:),X,XMAX,XMIN
      CHARACTER*80 BXCTXT,FILNAM,HEADER1,HEADER2,VXCTXT
      INTEGER COLORTAB(NCOLORTAB),ICOLOR_C(:),IM,IR,IRBOT,IRTOP,
     &        ISTYLE_C(:),IT,IVXCTAB,LBXCTXT,LFILNAM,LHEADER1,LHEADER2,
     &        LS,LTXT_IT,LVXCTXT,NGRAPH,NX
      LOGICAL GGA
      CHARACTER*20 LEG(NVXCTAB),STR20
      CHARACTER*10 VXCTAB(NVXCTAB)
C
C*** End of declarations rewritten by SPAG
C
      DATA VXCTAB/'VBH','MJW','VWN','PBE'/
      DATA LEG/'VBH    (LDA)','MJW   (LDA)','VWN   (LDA)',
     &     'PBE    (GGA)'/
      DATA COLORTAB/1,2,4,3,6,7,8,9,10,11,12,13,14,15/
C
      ALLOCATABLE RHO4PI,WEXC,EXC,RHO,RHOD,RHOU
      ALLOCATABLE VAUX,VXC,BXC,ICOLOR_C,ISTYLE_C
      ALLOCATABLE AGRDRHO,AGRDRHOD,AGRDRHOU
      ALLOCATABLE LAPRHO,LAPRHOD,LAPRHOU
      ALLOCATABLE GDGAG,GDGAGD,GDGAGU
C
      ALLOCATE (RHO4PI(NRMAX,2),WEXC(NRMAX,2),EXC(NRMAX))
      ALLOCATE (RHO(NRMAX),RHOD(NRMAX),RHOU(NRMAX))
      ALLOCATE (VAUX(NRMAX,2))
      ALLOCATE (VXC(NRMAX,NVXCTAB),BXC(NRMAX,NVXCTAB))
      ALLOCATE (ICOLOR_C(2*NVXCTAB),ISTYLE_C(2*NVXCTAB))
C
      IF ( FULLPOT ) RETURN
C=======================================================================
C                      GGA - parametrisations
C=======================================================================
      IF ( SCFVXC(1:3).EQ.'PBE' ) THEN
C
         GGA = .TRUE.
C
         ALLOCATE (AGRDRHO(NRMAX),AGRDRHOD(NRMAX),AGRDRHOU(NRMAX))
         ALLOCATE (LAPRHO(NRMAX),LAPRHOD(NRMAX),LAPRHOU(NRMAX))
         ALLOCATE (GDGAG(NRMAX),GDGAGD(NRMAX),GDGAGU(NRMAX))
C
      ELSE
C
         GGA = .FALSE.
C
      END IF
C=======================================================================
C
C-----------------------------------------------------------------------
C                      loop over types IT
C-----------------------------------------------------------------------
      DO IT = ITBOT,ITTOP
C
         IM = IMT(IT)
         IRBOT = 1
         IRTOP = JRWS(IM)
C
         VXCMIN = 0D0
         VXCMAX = 0D0
         BXCMIN = 0D0
         BXCMAX = 0D0
C
         DO IVXCTAB = 1,NVXCTAB
C
            SCFVXC = VXCTAB(IVXCTAB)
C
            DO IR = 1,IRTOP
               VAUX(IR,1) = 0D0
               VAUX(IR,2) = 0D0
               RHO4PI(IR,1) = RHOCHR(IR,IT)
               RHO4PI(IR,2) = RHOSPN(IR,IT)
            END DO
C
C --------------------------------------------- skip or add XC-potential
            DO IR = 1,IRTOP
               WEXC(IR,1) = 0D0
            END DO
C
C= LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA
            IF ( SCFVXC.EQ.'NONE      ' ) THEN
               DO IR = 1,IRTOP
                  WEXC(IR,1) = 0D0
               END DO
C
            ELSE IF ( SCFVXC(1:3).EQ.'VBH' ) THEN
C
               CALL EXCVBH(RHO4PI,VAUX,EXC,WEXC,IRTOP,NRMAX)
C
            ELSE IF ( SCFVXC(1:3).EQ.'MJW' ) THEN
C
               CALL EXCMJW(RHO4PI,VAUX,EXC,WEXC,IRTOP,NRMAX,1)
C
            ELSE IF ( SCFVXC(1:3).EQ.'VWN' ) THEN
C
               CALL EXCVWN(RHO4PI,VAUX,EXC,WEXC,IRTOP,NRMAX)
C
C= LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA
            ELSE IF ( GGA ) THEN
C= GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA
C
C-------------------- get spin-resolved densities and remove factor 4 PI
C
               X = 4.0D0*PI
               DO IR = 1,IRTOP
                  RHO(IR) = RHOCHR(IR,IT)/X
                  RHOU(IR) = 0.5D0*(RHOCHR(IR,IT)+RHOSPN(IR,IT))/X
                  RHOD(IR) = 0.5D0*(RHOCHR(IR,IT)-RHOSPN(IR,IT))/X
               END DO
C
C---------- calculate  AGRD: n' = d rho/dr  and  LAP: n'' = d^2 rho/dr^2
C
               CALL CALCDFDR(RHO,AGRDRHO,DRDI(1,IM),IRTOP)
               CALL CALCDFDR(AGRDRHO,LAPRHO,DRDI(1,IM),IRTOP)
               CALL CALCDFDR(RHOU,AGRDRHOU,DRDI(1,IM),IRTOP)
               CALL CALCDFDR(AGRDRHOU,LAPRHOU,DRDI(1,IM),IRTOP)
               CALL CALCDFDR(RHOD,AGRDRHOD,DRDI(1,IM),IRTOP)
               CALL CALCDFDR(AGRDRHOD,LAPRHOD,DRDI(1,IM),IRTOP)
C
               DO IR = 1,IRTOP
C
C------------------------------------------------- check whether n' <= 0
C
                  IF ( AGRDRHO(IR).GT.0.0 ) WRITE (6,*)
     &                  ' WARNING IN <SCFNEWPOT> D RHO/DR > 0.0 FOR R=',
     &                 RHO(IR),' IR=',IR
                  IF ( AGRDRHOU(IR).GT.0.0 ) WRITE (6,*) 
     &               ' WARNING IN <SCFNEWPOT> D RHO(UP)/DR > 0.0 FOR R='
     &               ,RHOU(IR),' IR=',IR
                  IF ( AGRDRHOD(IR).GT.0.0 ) WRITE (6,*) 
     &               ' WARNING IN <SCFNEWPOT> D RHO(DN)/DR > 0.0 FOR R='
     &               ,RHOD(IR),' IR=',IR
C
C ---------------------- GDGAG: grad RHO . grad | grad RHO | = n' (-n" )
C
                  GDGAG(IR) = -AGRDRHO(IR)*LAPRHO(IR)
                  GDGAGU(IR) = -AGRDRHOU(IR)*LAPRHOU(IR)
                  GDGAGD(IR) = -AGRDRHOD(IR)*LAPRHOD(IR)
C
C -------------------------------------- LAP: grad^2 RHO = n" + 2 n' / r
C
                  LAPRHO(IR) = LAPRHO(IR) + 2D0*AGRDRHO(IR)/R(IR,IM)
                  LAPRHOU(IR) = LAPRHOU(IR) + 2D0*AGRDRHOU(IR)/R(IR,IM)
                  LAPRHOD(IR) = LAPRHOD(IR) + 2D0*AGRDRHOD(IR)/R(IR,IM)
C
C ------------------------------------------ AGRD: | grad RHO | = | n' |
C
                  AGRDRHO(IR) = ABS(AGRDRHO(IR))
                  AGRDRHOU(IR) = ABS(AGRDRHOU(IR))
                  AGRDRHOD(IR) = ABS(AGRDRHOD(IR))
C
               END DO
C
               WEXC(1:NRMAX,1:2) = 0D0
C
               IF ( SCFVXC(1:3).EQ.'PBE' ) THEN
C
                  CALL EXCPBE(VAUX,EXC,WEXC,IRTOP,NRMAX,RHO,AGRDRHO,
     &                        LAPRHO,GDGAG,RHOU,AGRDRHOU,LAPRHOU,GDGAGU,
     &                        RHOD,AGRDRHOD,LAPRHOD,GDGAGD)
C
               ELSE
                  STOP '<SCFNEWPOT>:  no available GGA parametrisation'
               END IF
C
               X = 4.0D0*PI
               WEXC(1:IRTOP,1:2) = X*WEXC(1:IRTOP,1:2)
C
C= GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA
            ELSE
               WRITE (6,99001) SCFVXC
               STOP
            END IF
C-----------------------------------------------------------------------
C
C ---------------------- convention VAUX(I,IS) IS=1,2 == up, down (AKAI)
C
            DO IR = 1,IRTOP
C
               VXC(IR,IVXCTAB) = (VAUX(IR,1)+VAUX(IR,2))/2.0D0
               BXC(IR,IVXCTAB) = (VAUX(IR,1)-VAUX(IR,2))/2.0D0
C
               VXCMIN = MIN(VXCMIN,VXC(IR,IVXCTAB))
               VXCMAX = MAX(VXCMAX,VXC(IR,IVXCTAB))
               BXCMIN = MIN(BXCMIN,BXC(IR,IVXCTAB))
               BXCMAX = MAX(BXCMAX,BXC(IR,IVXCTAB))
C
            END DO
C
         END DO
C
         IF ( ABS(BXCMIN).LT.1D-4 .AND. ABS(BXCMAX).LT.1D-4 )
     &        NONMAG = .TRUE.
C
         VXCMIN = -15D0
         VXCMAX = 0D0
C
C ----------------------------------------------------------------------
C                     plot VXC, BXC
C ----------------------------------------------------------------------
C
         STR20 = TXT_T(IT)(1:LTXT_T(IT))//'!N(r) (Ry)'
         LS = LTXT_T(IT) + 10
         VXCTXT = 'V!m{1}!Sxc!M{1}!s'//STR20(1:LS)
         LVXCTXT = 17 + LS
         IF ( NONMAG ) THEN
            NGRAPH = 1
            BXCTXT = ' '
            LBXCTXT = 1
         ELSE
            NGRAPH = 2
            STR20 = TXT_T(IT)(1:LTXT_T(IT))//'!N(r) (Ry)'
            LS = LTXT_T(IT) + 10
            BXCTXT = 'B!m{1}!Sxc!M{1}!s'//STR20(1:LS)
            LBXCTXT = 17 + LS
         END IF
C
         NX = IRTOP - IRBOT + 1
C
         XMIN = R(IRBOT,IM)
         XMAX = R(IRTOP,IM)
C
         HEADER2 = 'spherical xc-potential of '//TXT_T(IT)(1:LTXT_T(IT))
         LHEADER2 = 28 + LTXT_T(IT)
         IF ( TUTORIAL ) THEN
            HEADER1 = 'SPR-KKR tutorial'
            LHEADER1 = 16
            LTXT_IT = 0
         ELSE
            HEADER1 = 'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM)
            LHEADER1 = 25 + LSYSTEM
            LTXT_IT = LTXT_T(IT)
            HEADER2 = HEADER2(1:LHEADER2)//' in '//SYSTEM(1:LSYSTEM)
            LHEADER2 = LHEADER2 + 4 + LSYSTEM
         END IF
C
         CALL XMGRHEAD(DATSET0,LDATSET0,'Vxc',3,TXT_T(IT),LTXT_IT,
     &                 FILNAM,80,LFILNAM,IOTMP,NGRAPH,XMIN,1,XMAX,1,
     &                 VXCMIN,1,VXCMAX,1,BXCMIN,1,BXCMAX,1,
     &                 'radius r (a.u.)',15,VXCTXT,LVXCTXT,BXCTXT,
     &                 LBXCTXT,HEADER1,LHEADER1,HEADER2,LHEADER2,
     &                 .FALSE.)
C
         CALL XMGRLEG1(IOTMP,0,NVXCTAB,LEG,0.45D0,0.45D0)
C
         CALL XMGRCURVES(IOTMP,NGRAPH,NVXCTAB,NVXCTAB,3,1,0)
         ICOLOR_C(1:NVXCTAB) = COLORTAB(1:NVXCTAB)
         ISTYLE_C(1:NVXCTAB) = 1
C
         CALL XMGRCOLOR(IOTMP,0,NVXCTAB,ICOLOR_C)
         CALL XMGRSTYLE(IOTMP,0,NVXCTAB,ISTYLE_C)
C
         DO IVXCTAB = 1,NVXCTAB
            CALL XMGRTABLE(0,(IVXCTAB-1),R(IRBOT,IM),VXC(IRBOT,IVXCTAB),
     &                     1D0,NX,IOTMP)
         END DO
C
         IF ( .NOT.NONMAG ) THEN
            ICOLOR_C((NVXCTAB+1):(2*NVXCTAB)) = COLORTAB(1:NVXCTAB)
            ISTYLE_C((NVXCTAB+1):(2*NVXCTAB)) = 3
C
            CALL XMGRCOLOR(IOTMP,1,2*NVXCTAB,ICOLOR_C)
            CALL XMGRSTYLE(IOTMP,1,2*NVXCTAB,ISTYLE_C)
C
            DO IVXCTAB = 1,NVXCTAB
               CALL XMGRTABLE(1,(IVXCTAB-1),R(IRBOT,IM),
     &                        BXC(IRBOT,IVXCTAB),1D0,NX,IOTMP)
            END DO
C
         END IF
C
         WRITE (6,99002) 'spherical     potential V,   B    ',
     &                   FILNAM(1:LFILNAM)
         CLOSE (IOTMP)
C
      END DO
C
C=======================================================================
C
99001 FORMAT (//,1X,79('#'),/,10X,'ERROR in <TUTORIAL_POT_XC> ',/,10X,
     &        'SCFVXC = ',A,'  not known ',/,1X,79('#'),/)
99002 FORMAT (/,10X,A,' written to file ',A)
      END
C*==tutorial_core_all.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TUTORIAL_CORE_ALL
C   ********************************************************************
C   *                                                                  *
C   *  subroutine to plot ALL core wave functions                      *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKM
      USE MOD_TYPES,ONLY:NTMAX,NCPLWFMAX,IMT,LTXT_T,TXT_T,LCXRAY,NCXRAY,
     &    Z,NCORT,VT,BT
      USE MOD_FILES,ONLY:IOTMP,LRECREAL8,LDATSET0,DATSET0,IPRINT
      USE MOD_RMESH,ONLY:FULLPOT,NRMAX,R,JRWS,RWS,R2DRDI
      USE MOD_CONSTANTS,ONLY:RY_EV
      USE MOD_CALCMODE,ONLY:NONMAG
      IMPLICIT NONE
C*--TUTORIAL_CORE_ALL1192
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NLSHELLMAX,NCSTMAX,NNKCORMAX
      PARAMETER (NLSHELLMAX=15,NCSTMAX=14,NNKCORMAX=2*NLSHELLMAX)
C
C Local variables
C
      REAL*8 BCOR(NTMAX),BCORS(NTMAX),EBOT,ECOR(NCSTMAX),E_NKCOR(:),
     &       FCOR(:,:,:),F_NKCOR(:,:),GCOR(:,:,:),G_NKCOR(:,:),
     &       HX_NKCOR(:,:),HY_NKCOR(:,:),MJ,RINT(NRMAX),SZCOR(NCSTMAX),
     &       X,XNORM(2)
      CHARACTER*80 FILNAM,YTXT0,YTXT1
      CHARACTER*5 FUNTXTMJ
      INTEGER I,ICOLOR_C(:),ICST,IFIL,IKMCOR(NCSTMAX,2),
     &        IKMCPLWFCOR(NCPLWFMAX),ILSHELL,IM,INKCOR,IOL,IR0_NKCOR(:),
     &        IRTOP,ISF,ISG,ISH,IT,ITXRAY,IVBOT,IZERO(NCSTMAX),J,K,
     &        KAPCOR(NCSTMAX),L,LFILNAM,LQNTAB(NLSHELLMAX),LYTXT,MM05,
     &        MM05COR(NCSTMAX),N,NC,NCST,NGRAPH,NKPCOR(NCSTMAX),NLSHELL,
     &        NNKCOR,NQNTAB(NLSHELLMAX),NSOL
      CHARACTER*20 LEG(1),STR20,TXT_NKCOR(:)
      CHARACTER*3 TXTK(4)
      CHARACTER*1 TXT_L(0:3)
C
C*** End of declarations rewritten by SPAG
C
      DATA NQNTAB/1,2,2,3,3,3,4,4,4,5,5,4,5,6,6/
      DATA LQNTAB/0,0,1,0,1,2,0,1,2,0,1,3,2,0,1/
      DATA TXT_L/'s','p','d','f'/
      DATA TXTK/'1/2','3/2','5/2','7/2'/
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE ICOLOR_C,FCOR,GCOR
      ALLOCATABLE E_NKCOR,G_NKCOR,F_NKCOR,IR0_NKCOR
      ALLOCATABLE HY_NKCOR,HX_NKCOR,TXT_NKCOR
C
      IF ( FULLPOT ) THEN
         NC = NCPLWFMAX
         STOP 'in <WFPLOT>    FULLPOT-option not yet available'
      ELSE
         NC = 2
      END IF
      IF ( NC.LE.0 ) STOP '<WFPLOT>  array sizes ???'
C
      ALLOCATE (FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX))
      ALLOCATE (ICOLOR_C(3*NNKCORMAX))
      ALLOCATE (E_NKCOR(NNKCORMAX),IR0_NKCOR(NNKCORMAX))
      ALLOCATE (G_NKCOR(NRMAX,NNKCORMAX),F_NKCOR(NRMAX,NNKCORMAX))
      ALLOCATE (HY_NKCOR(2,NNKCORMAX),HX_NKCOR(2,NNKCORMAX))
      ALLOCATE (TXT_NKCOR(NNKCORMAX))
      HX_NKCOR(1:2,1:NNKCORMAX) = 0D0
      E_NKCOR(1:NNKCORMAX) = 0D0
      IR0_NKCOR(1:NNKCORMAX) = 0
      G_NKCOR(1:NRMAX,1:NNKCORMAX) = 0D0
      F_NKCOR(1:NRMAX,1:NNKCORMAX) = 0D0
      TXT_NKCOR(1:NNKCORMAX) = ' '
C
      BT(1:NRMAX,1:NTMAX) = 0D0
      NONMAG = .TRUE.
C
      WRITE (6,99004)
C
C=======================================================================
C
      IFIL = 87
C
      IF ( FULLPOT ) THEN
         IOL = (3+NCPLWFMAX+1+NCPLWFMAX*2*NRMAX*2)*LRECREAL8
      ELSE
         IOL = (5+2*(2+2*NRMAX*2))*LRECREAL8
      END IF
C
      OPEN (UNIT=IFIL,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=IOL)
C
C-----------------------------------------------------------------------
C
      IT = 1
      IM = IMT(IT)
      WRITE (6,99007) IT
C
C ======================================================================
C                             CORE state
C ======================================================================
C
      ITXRAY = IT
C
      NLSHELL = 0
      IF ( Z(IT).GT.2 ) NLSHELL = 1
      IF ( Z(IT).GT.10 ) NLSHELL = 3
      IF ( Z(IT).GT.18 ) NLSHELL = 5
      IF ( Z(IT).GT.30 ) NLSHELL = 6
      IF ( Z(IT).GT.36 ) NLSHELL = 8
      IF ( Z(IT).GT.48 ) NLSHELL = 9
      IF ( Z(IT).GT.54 ) NLSHELL = 11
      IF ( Z(IT).GT.70 ) NLSHELL = 12
      IF ( Z(IT).GT.80 ) NLSHELL = 13
      IF ( Z(IT).GT.86 ) NLSHELL = 15
C
      IF ( NCORT(IT).NE.0 ) THEN
         NLSHELL = 0
         N = 0
         DO ILSHELL = 1,NLSHELLMAX
            L = LQNTAB(ILSHELL)
            N = N + 2*(2*L+1)
            IF ( N.EQ.NCORT(IT) ) NLSHELL = ILSHELL
         END DO
         IF ( NLSHELL.EQ.0 ) THEN
            WRITE (6,*) 'NLSHELL not found for IT=',IT,' NCORT=',
     &                  NCORT(IT)
            STOP ' in <CORE>'
         END IF
      END IF
C
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C                   ---------------------------------------
C                   INITIALIZE QUANTUM NUMBERS  NQN  AND  L
C                   ---------------------------------------
      INKCOR = 0
      EBOT = 0D0
      DO ILSHELL = 1,NLSHELL
         NCXRAY(IT) = NQNTAB(ILSHELL)
         LCXRAY(IT) = LQNTAB(ILSHELL)
C
         MJ = 0.5D0
C
         MM05 = NINT(MJ-0.5D0)
C
         NCST = 4*LCXRAY(IT) + 2
C
         STR20 = FUNTXTMJ(MJ)
         WRITE (6,99005) NCXRAY(IT),LCXRAY(IT),STR20(1:5)
C
         IF ( FULLPOT ) THEN
C
            CALL CORE(IPRINT,GCOR,FCOR,ECOR,SZCOR,KAPCOR,MM05COR,NKPCOR,
     &                IKMCOR,IZERO,ITXRAY,BCOR,BCORS,NCSTMAX)
C
            DO ICST = 1,NCST
               IF ( ABS(MM05COR(ICST)+0.5D0).GT.LCXRAY(IT) ) THEN
                  NSOL = 1
               ELSE
                  NSOL = 2
               END IF
               DO J = 1,NCPLWFMAX
                  IKMCPLWFCOR(J) = 0
               END DO
               DO J = 1,NSOL
                  IKMCPLWFCOR(J) = IKMCOR(ICST,J)
               END DO
               WRITE (IFIL,REC=IKMCOR(ICST,1)+(IT-1)*NKM) IT,'COR',
     &                IKMCOR(ICST,1),JRWS(IM),
     &                (IKMCPLWFCOR(J),J=1,NCPLWFMAX),NSOL,
     &                ((DCMPLX(GCOR(I,K,ICST),0.0D0),I=1,JRWS(IM)),K=1,
     &                NSOL),
     &                ((DCMPLX(FCOR(I,K,ICST),0.0D0),I=1,JRWS(IM)),K=1,
     &                NSOL)
C
            END DO
C
         ELSE
C
            CALL CORE(IPRINT,GCOR,FCOR,ECOR,SZCOR,KAPCOR,MM05COR,NKPCOR,
     &                IKMCOR,IZERO,ITXRAY,BCOR,BCORS,NCSTMAX)
C
            WRITE (6,99001)
C
            DO ICST = 1,NCST
               IF ( MM05.EQ.MM05COR(ICST) ) THEN
C
                  DO K = 1,NKPCOR(ICST)
                     DO N = 1,JRWS(IM)
                        RINT(N) = R2DRDI(N,IM)
     &                            *(GCOR(N,K,ICST)**2+FCOR(N,K,ICST)**2)
                     END DO
                     CALL RRADINT(IM,RINT,XNORM(K))
                  END DO
C
                  WRITE (6,99002) ICST,NCXRAY(IT),LCXRAY(IT),
     &                            KAPCOR(ICST),(2*MM05COR(ICST)+1),
     &                            IKMCOR(ICST,1),XNORM(1),ECOR(ICST),
     &                            ECOR(ICST)*RY_EV,SZCOR(ICST),
     &                            IZERO(ICST)
                  IF ( NKPCOR(ICST).EQ.2 ) WRITE (6,99003)
     &                 IKMCOR(ICST,2),XNORM(2)
C
                  INKCOR = INKCOR + 1
                  NNKCOR = INKCOR
                  IR0_NKCOR(INKCOR) = IZERO(ICST)
                  E_NKCOR(INKCOR) = ECOR(ICST)
                  EBOT = MIN(EBOT,E_NKCOR(INKCOR))
                  HY_NKCOR(1,INKCOR) = ECOR(ICST)
                  HY_NKCOR(2,INKCOR) = ECOR(ICST)
                  HX_NKCOR(1,INKCOR) = 0D0
                  HX_NKCOR(2,INKCOR) = R(IZERO(ICST),IM)
                  TXT_NKCOR(INKCOR) = ' '
                  CALL STRING_ADD_N(TXT_NKCOR(INKCOR),NCXRAY(IT))
                  TXT_NKCOR(INKCOR)(2:7) = TXT_L(LCXRAY(IT))
     &               //'!s'//TXTK(ABS(KAPCOR(ICST)))
C
                  DO I = 1,JRWS(IM)
                     G_NKCOR(I,INKCOR) = GCOR(I,1,ICST)
                     F_NKCOR(I,INKCOR) = FCOR(I,1,ICST)
                  END DO
C
               END IF
C
            END DO
C
         END IF
C
      END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
      NGRAPH = 1
      IF ( ABS(MJ).GT.LCXRAY(IT) ) YTXT0 = 'g(r)    (j=l+1/2)'
      YTXT0 = 'E!sn!xk!0!N (Ry)    g!sn!xk!0!N(r) f!sn!xk!0!N(r)'
      YTXT1 = ' '
      LYTXT = LEN_TRIM(YTXT0)
C
      EBOT = (EBOT)*1.05D0
C
      IVBOT = 1
      DO I = 1,JRWS(IM)
         IF ( EBOT.LT.VT(I,IT) ) THEN
            IVBOT = I
            EXIT
         END IF
      END DO
C
      CALL XMGRHEAD(DATSET0,LDATSET0,'core-all',8,' ',0,FILNAM,80,
     &              LFILNAM,IOTMP,NGRAPH,0.0D0,1,RWS(IM),1,EBOT,0,0D0,1,
     &              0D0,0,0D0,0,'radius (a!s0!N)',15,YTXT0,LYTXT,YTXT1,
     &              LYTXT,'SPR-KKR tutorial',16,
     &              'core wave functions of '//TXT_T(IT)(1:LTXT_T(IT)),
     &              (23+LTXT_T(IT)),.FALSE.)
C
      LEG(1) = 'V(r)'
      CALL XMGRLEG1(IOTMP,0,1,LEG,0.6D0,0.45D0)
C
      CALL XMGRCURVES(IOTMP,1,3*NNKCOR+1,0,2,0,0)
      CALL XMGRCURVES(IOTMP,1,NNKCOR+1,0,1,0,0)
      CALL XMGRCURVES(IOTMP,1,1,0,3,0,0)
C
C V: green4   E: black  g: blue  f: red
      ICOLOR_C(1) = 15
      ICOLOR_C(2:1+NNKCOR) = 1
      ICOLOR_C(2+NNKCOR:2*NNKCOR+1) = 4
      ICOLOR_C(2+2*NNKCOR:3*NNKCOR+1) = 2
C
      CALL XMGRCOLOR(IOTMP,0,3*NNKCOR+1,ICOLOR_C)
C
      CALL XMGRTABLE(0,0,R(IVBOT,IM),VT(IVBOT,IT),1D0,(JRWS(IT)-IVBOT+1)
     &               ,IOTMP)
C
      DO INKCOR = 1,NNKCOR
C
         X = 0.2D0
C
         IRTOP = IR0_NKCOR(INKCOR)
         G_NKCOR(1:IRTOP,INKCOR) = E_NKCOR(INKCOR)
     &                             + X*G_NKCOR(1:IRTOP,INKCOR)
         F_NKCOR(1:IRTOP,INKCOR) = E_NKCOR(INKCOR)
     &                             + X*F_NKCOR(1:IRTOP,INKCOR)
C
         CALL XMGRTEXT(IOTMP,TXT_NKCOR(INKCOR),
     &                 LEN_TRIM(TXT_NKCOR(INKCOR)),HX_NKCOR(2,INKCOR)
     &                 +0.1D0,HY_NKCOR(2,INKCOR),'WORLD')
C
         ISH = INKCOR
         CALL XMGRTABLE(0,ISH,HX_NKCOR(1,INKCOR),HY_NKCOR(1,INKCOR),
     &                  1.0D0,2,IOTMP)
C
         ISG = NNKCOR + ISH
         CALL XMGRTABLE(0,ISG,R(1,IM),G_NKCOR(1,INKCOR),1.0D0,
     &                  IR0_NKCOR(INKCOR),IOTMP)
C
         ISF = NNKCOR + ISG
         CALL XMGRTABLE(0,ISF,R(1,IM),F_NKCOR(1,INKCOR),1.0D0,
     &                  IR0_NKCOR(INKCOR),IOTMP)
C
      END DO
C
      WRITE (6,99006) 'core',FILNAM(1:LFILNAM)
C
      STOP
99001 FORMAT (/,' ICST  N   L  KAP  MUE  IKM     NORM   ',
     &        '     E(Ry)       E(eV)    <SIGMA_z>  I0')
99002 FORMAT (5I4,'/2',I4,F12.6,F12.4,F12.3,F12.4,I5)
99003 FORMAT (22X,I4,F12.6)
99004 FORMAT (/,1X,79('*'),/,33X,'<TUTORIAL_CORE_ALL>',/,28X,
     &        'plot ALL core wave functions',/,1X,79('*'),/,/,10X,
     &        'NOTE: non-magnetic calculation for  B=0',/)
99005 FORMAT (/,10X,'core wave functions to be plotted for',
     &        ' quantum numbers:',//,22X,'n = ',I1,'  l = ',I1,
     &        '  mj = ',A5,/)
99006 FORMAT (/,10X,A,' wave functions written to file: ',A,/)
99007 FORMAT (//,10X,'selected atom type   IT =',I3,/)
      END
C*==tutorial_rwf.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TUTORIAL_RWF(KEY,NORMALISATION,L,TSST,MSST,IPRINT,SSST,
     &                        MEZZ,MEZJ)
C   ********************************************************************
C   *                                                                  *
C   *  subroutine to plot    wave functions                            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:ORBPOL,IREL
      USE MOD_ENERGY,ONLY:ETAB,NETAB,NE
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX,NSPIN,NL,NCPLWF
      USE MOD_TYPES,ONLY:NTMAX,NCPLWFMAX,IMT,LTXT_T,TXT_T,NT,IKMCPLWF
      USE MOD_FILES,ONLY:IOTMP,LRECREAL8,LSYSTEM,SYSTEM,LDATSET,DATSET,
     &    IFILCBWF,FOUND_SECTION
      USE MOD_RMESH,ONLY:FULLPOT,NRMAX,R,JRWS,RWS,R2DRDI_W_RADINT
      IMPLICIT NONE
C*--TUTORIAL_RWF1524
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER FMATCH,NCOLORTAB
      PARAMETER (FMATCH=4,NCOLORTAB=14)
C
C Dummy arguments
C
      INTEGER IPRINT,L
      CHARACTER*1 KEY
      CHARACTER*10 NORMALISATION
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),SSST(NKMMAX,NKMMAX,NTMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 ARG,CJL,ERYD,JG(:,:,:),P,ZG(:,:,:)
      COMPLEX*16 CJLZ
      INTEGER COLORTAB(NCOLORTAB),I,ICOLOR_C(:),IE,IG,IL,ILS,IM,IOL,IR,
     &        IRTOP,IS,ISPIN,IT,ITSEL,LFILNAM,LYTXT,NC,NGRAPH,NMATCH,
     &        NSG(0:1)
      REAL*8 DR,NORM,XMATCH(FMATCH*NRMAX),YMATCH(FMATCH*NRMAX,2),YMAX,
     &       YMIN,YY(:,:,:)
      CHARACTER*80 FILNAM,YTXT0,YTXT1
      CHARACTER*1 FUNTXTL
      LOGICAL GETIRRSOL
      CHARACTER*20 LEG(:,:)
      CHARACTER*40 STR40
C
C*** End of declarations rewritten by SPAG
C
      DATA COLORTAB/1,2,4,3,6,7,8,9,10,11,12,13,14,15/
C
      ALLOCATABLE ZG,JG,YY,LEG,ICOLOR_C
C
      IF ( KEY.EQ.'L' ) THEN
         NC = MAX(1,NL+1)
      ELSE
         NC = MAX(1,NETAB(1)+1)
      END IF
C
      ALLOCATE (ZG(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JG(NRMAX,NCPLWFMAX,NKMMAX))
C
      GETIRRSOL = .TRUE.
      IREL = 1
C
C=======================================================================
C
      IF ( FULLPOT ) THEN
         IOL = (3+NCPLWFMAX+1+NCPLWFMAX*2*NRMAX*2)*LRECREAL8
      ELSE
         IOL = (5+2*(2+2*NRMAX*2))*LRECREAL8
      END IF
C
      OPEN (UNIT=IFILCBWF,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=IOL)
C
C-----------------------------------------------------------------------
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
      IF ( FOUND_SECTION ) THEN
         ITSEL = 1
         CALL SECTION_SET_INTEGER('ITSEL',ITSEL,9999,1)
         IT = ITSEL
         IM = IMT(IT)
         IRTOP = JRWS(IM)
         WRITE (6,99002) ITSEL
C
      ELSE
         STOP 'in <WFPLOT>   TASK not found'
      END IF
C
C ======================================================================
C                             BAND state
C ======================================================================
C
CLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      IF ( KEY.EQ.'L' ) THEN
C
         NSG(0:1) = NL + 1
C
         ALLOCATE (YY(NRMAX,0:NL,0:1),LEG(0:NL,0:1),ICOLOR_C(NC))
C
         ERYD = ETAB(1,1)
         IF ( ABS(ERYD).LT.1D-7 ) ERYD = 1D-5
C
         YMIN = 0D0
         YMAX = 0D0
C
         CALL RUNSSITE(.FALSE.,1,1,IFILCBWF,GETIRRSOL,ERYD,P,IPRINT,
     &                 TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
         CALL WAVFUN_READ_SRA(IFILCBWF,IT,1,ZG,JG,IRTOP,NCPLWF,IKMCPLWF)
C
         LEG(0:NL,0:1) = ''
C
         DO IL = 1,NL
            L = IL - 1
C
            LEG(IL,0) = FUNTXTL(L)
C
            DO ISPIN = 1,NSPIN
C
               ILS = NL*(ISPIN-1) + IL
C
               DO IR = 1,IRTOP
                  YY(IR,IL,ISPIN-1) = DREAL(ZG(IR,1,ILS))
               END DO
C
               IF ( NORMALISATION(1:1).EQ.'1' ) THEN
                  NORM = 0D0
                  DO IR = 1,IRTOP
                     NORM = NORM + R2DRDI_W_RADINT(IR,IM)
     &                      *YY(IR,IL,ISPIN-1)**2
                  END DO
                  NORM = 1D0/NORM
                  YY(1:IRTOP,IL,ISPIN-1) = YY(1:IRTOP,IL,ISPIN-1)*NORM
               END IF
            END DO
C
            DO ISPIN = 1,NSPIN
               DO IR = 1,IRTOP
                  YMIN = MIN(YMIN,YY(IR,IL,ISPIN-1))
                  YMAX = MAX(YMAX,YY(IR,IL,ISPIN-1))
               END DO
            END DO
C
         END DO
C
         WRITE (STR40,FMT='(''E = ('',F6.2,'','',F6.2,'') Ry'')') ERYD
C
CLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      ELSE
CEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
         NE = NETAB(1)
         NSG(0:1) = NE
C
         ALLOCATE (YY(NRMAX,0:NE,0:1),LEG(0:NE,0:1),ICOLOR_C(NE+1))
C
         YMIN = 0D0
         YMAX = 0D0
C
         LEG(0:NE,0:1) = ''
C
         DO IE = 1,NE
C
            ERYD = ETAB(IE,1)
            IF ( ABS(ERYD).LT.1D-7 ) ERYD = 1D-5
C
            WRITE (LEG(IE,0),'(''E= '',f6.2,'' Ry'')') DREAL(ERYD)
C
            CALL RUNSSITE(.FALSE.,1,1,IFILCBWF,GETIRRSOL,ERYD,P,IPRINT,
     &                    TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
            CALL WAVFUN_READ_SRA(IFILCBWF,IT,1,ZG,JG,IRTOP,NCPLWF,
     &                           IKMCPLWF)
C
            DO ISPIN = 1,NSPIN
C
               ILS = NL*(ISPIN-1) + IL
C
               DO IR = 1,IRTOP
                  YY(IR,IE,ISPIN-1) = DREAL(ZG(IR,1,ILS))
               END DO
C
               IF ( NORMALISATION(1:1).EQ.'1' ) THEN
                  NORM = 0D0
                  DO IR = 1,IRTOP
                     NORM = NORM + R2DRDI_W_RADINT(IR,IM)
     &                      *YY(IR,IE,ISPIN-1)**2
                  END DO
                  NORM = 1D0/NORM
                  YY(1:IRTOP,IE,ISPIN-1) = YY(1:IRTOP,IE,ISPIN-1)*NORM
               END IF
            END DO
C
            DO ISPIN = 1,NSPIN
               DO IR = 1,IRTOP
                  YMIN = MIN(YMIN,YY(IR,IE,ISPIN-1))
                  YMAX = MAX(YMAX,YY(IR,IE,ISPIN-1))
               END DO
            END DO
C
         END DO
C
         WRITE (STR40,FMT='(''E = ('',F6.2,'','',F6.2,'') Ry'')') ERYD
C
      END IF
CEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
C
      IF ( NC.LE.NCOLORTAB ) THEN
         ICOLOR_C(1:NC) = COLORTAB(1:NC)
      ELSE
         ICOLOR_C(1:NCOLORTAB) = COLORTAB(1:NCOLORTAB)
         ICOLOR_C((NCOLORTAB+1):NC) = 1
      END IF
C
      NGRAPH = NSPIN
C
      IF ( NGRAPH.EQ.1 ) THEN
         YTXT0 = 'R!sl!N(r)  (a.u.)'
         LYTXT = 17
      ELSE
         YTXT0 = 'R!m{1}!sl!N!M{1}!S!UP!N(r)  (a.u.)'
         YTXT1 = 'R!m{1}!sl!N!M{1}!S!DN!N(r)  (a.u.)'
         LYTXT = 34
      END IF
C
      CALL XMGRHEAD(DATSET,LDATSET,'max',3,' ',0,FILNAM,80,LFILNAM,
     &              IOTMP,NGRAPH,0.0D0,1,RWS(IM),1,YMIN,0,YMAX,1,YMIN,0,
     &              YMAX,0,'radius (a!s0!N)',15,YTXT0,LYTXT,YTXT1,LYTXT,
     &              'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM),
     &              25+LSYSTEM,'valence band wave functions of '//
     &              TXT_T(IT)(1:LTXT_T(IT)),(31+LTXT_T(IT)),.FALSE.)
C
      CALL XMGRLEGEND(IOTMP,NGRAPH,NSG(0),NSG(1),LEG(0,0),LEG(0,1))
C
      CALL XMGRCURVES(IOTMP,NGRAPH,NSG(0),NSG(1),2,1,0)
C
      DO IG = 0,(NGRAPH-1)
         DO IS = 0,(NSG(IG)-1)
            CALL XMGRTABLE(IG,IS,R(1,IM),YY(1,IS,IG),1.0D0,IRTOP,IOTMP)
         END DO
C
         CALL XMGRCOLOR(IOTMP,IG,NC,ICOLOR_C)
      END DO
C
      IF ( KEY.EQ.'L' ) THEN
         CALL XMGRTEXT(IOTMP,STR40,40,0.40D0,0.25D0,'VIEW ')
      ELSE
         STR40 = 'wave function for l = '
         CALL STRING_ADD_N(STR40,L)
         CALL XMGRTEXT(IOTMP,STR40,40,0.40D0,0.25D0,'VIEW ')
      END IF
C
      IF ( NORMALISATION(1:1).EQ.'1' ) THEN
         CALL XMGRTEXT(IOTMP,'normalized to 1',15,0.40D0,0.20D0,'VIEW ')
      ELSE
         CALL XMGRTEXT(IOTMP,'normalized according to MST',27,0.40D0,
     &                 0.20D0,'VIEW ')
      END IF
C
      WRITE (6,99001) 'valence band',FILNAM(1:LFILNAM)
      CLOSE (IOTMP)
C
      IF ( NT.GT.0 ) STOP
C
C-------------------------------------------------------------- MATCHING
C                            only for 1 wave function !!!!!!!!
C-------------------------------------------------------------- MATCHING
      NGRAPH = 1
C
      DR = R(IRTOP,IM) - R(IRTOP-1,IM)
      DO I = 1,IRTOP
         XMATCH(I) = R(I,IM)
         YMATCH(I,1) = YY(I,0,0)
         ARG = P*XMATCH(I)
         YMATCH(I,2) = DREAL(CJLZ(L,ARG))
      END DO
C
      NMATCH = FMATCH*NRMAX
      DO I = IRTOP + 1,NMATCH
         XMATCH(I) = XMATCH(I-1) + DR
         ARG = P*XMATCH(I)
         CJL = CJLZ(L,ARG)
         YMATCH(I,2) = DREAL(CJL)
         IF ( XMATCH(I).GT.15D0 ) THEN
            NMATCH = I
            EXIT
         END IF
      END DO
C
C
      CALL XMGRHEAD(DATSET,LDATSET,'WF_match_l'//FUNTXTL(L),11,TXT_T(IT)
     &              ,LTXT_T(IT),FILNAM,80,LFILNAM,IOTMP,NGRAPH,0.0D0,1,
     &              XMATCH(NMATCH),1,YMIN,0,YMAX,1,YMIN,0,YMAX,0,
     &              'radius (a!s0!N)',15,YTXT0,LYTXT,YTXT1,LYTXT,
     &              'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM),
     &              25+LSYSTEM,
     &              'matching of valence band wave function for '//
     &              TXT_T(IT)(1:LTXT_T(IT)),(43+LTXT_T(IT)),.FALSE.)
C
Cccc        CALL XMGRLEGEND(IOTMP,NGRAPH,NSG(0),NSG(1),LEG(1,0),LEG(1,1))
C
Cccc        CALL XMGRCURVES(IOTMP,NGRAPH,NSG(0),NSG(1),2,1,0)
C
C
      DO IG = 0,(NGRAPH-1)
         DO IS = 0,1
            CALL XMGRTABLE(IG,IS,XMATCH,YMATCH(1,IS+1),1.0D0,NMATCH,
     &                     IOTMP)
         END DO
      END DO
C
      WRITE (6,99001) 'valence band matching ',FILNAM(1:LFILNAM)
      CLOSE (IOTMP)
C
C-------------------------------------------------------------- MATCHING
C
C ======================================================================
C
      STOP
99001 FORMAT (/,10X,A,' wave functions written to file: ',A,/)
99002 FORMAT (//,10X,'selected atom type   IT =',I3,/)
      END
