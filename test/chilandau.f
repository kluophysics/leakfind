C*==chilandau.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHILANDAU(IESORT,IPROCE)
C   ********************************************************************
C   *                                                                  *
C   *  main subroutine for susceptibility calculations                 *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_CPA,ONLY:CPATOL,NCPA
      USE MOD_KSPACE,ONLY:IBZINT
      USE MOD_ANGMOM,ONLY:NL,NKMMAX,NLMAX,TXT_L,NKM,MEZJ,MEZZ,TSSQ,MSSQ,
     &    MSST,SSST,TAUQ,TAUT,TSST
      USE MOD_SITES,ONLY:NQ,IQAT,ITOQ,NOQ
      USE MOD_TYPES,ONLY:NT,NTMAX,NLT,NAT,TXT_T,CONC
      USE MOD_CALCMODE,ONLY:ORBPOL
      USE MOD_ENERGY,ONLY:ETAB,IGRID,NETAB,EILOW,EIMAG,EMAX,EMIN,EFERMI,
     &    WETAB,NEMAX
      USE MOD_FILES,ONLY:IFILCBWF,LSYSTEM,SYSTEM,LDATSET,DATSET,IPRINT,
     &    WRTAU,RDTAU,CMDUC,FOUND_SECTION
      USE MOD_MPI,ONLY:MPI,NPROCS,MPI_ID
      USE MOD_CONSTANTS,ONLY:CHI_AU2CGS
      IMPLICIT NONE
C*--CHILANDAU24
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CHILANDAU')
      REAL*8 CU
      PARAMETER (CU=CHI_AU2CGS*1D+6)
C
C Dummy arguments
C
      INTEGER IESORT(NEMAX),IPROCE(NEMAX)
C
C Local variables
C
      COMPLEX*16 CHILAN,ERYD,P,PHASK(1)
      REAL*8 CHILANLT(0:NLMAX,NTMAX,2),CLURAD,CPACHNG,CPACHTAB(NEMAX),
     &       RHOCHRC(NRMAX,NTMAX),RHOSPNC(NRMAX,NTMAX),
     &       RWORK(0:NLMAX,NTMAX,2),TIME1,TIME2
      INTEGER CHIPRINT,I,ICPACONV,ICPAFLAG,IDIMS,IE,IECPAFAIL(NEMAX),
     &        IECURR,IFMT,IL,INC,IPROC,IQ,IT,ITCPA,ITERM,IWRI,JE,
     &        LTAUFIL,NCPAFAIL,NQ9,NT9
      CHARACTER*80 LINE,TAUFIL
      CHARACTER*40 PATH(0:7)
C
C*** End of declarations rewritten by SPAG
C
      DATA PATH/'ONLY REAL ENERGIES   Gauss-Legendre    ',
     &     'ONLY REAL ENERGIES   Trapez rule       ',
     &     'RECTANGULAR COMPLEX PATH               ',
     &     'STRAIGHT CPLX. PATH PARA. TO REAL AXIS ',
     &     'RECT. CPLX. GRID, RETURN ON LOG SCALE  ',
     &     'ARC IN COMPLEX PLANE                   ',
     &     'STANDARD X-RAY E-MESH                  ',
     &     'X-RAY E-MESH for integration           '/
C
C ======================================================================
C
C-----------------------------------------------------------------------
C
      IWRI = 6
C
C-----------------------------------------------------------------------
C
      WRITE (6,99005)
C
C=======================================================================
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
C=======================================================================
C     prepare calculation of site off-diagonal TAU for CHI_Landau
C=======================================================================
C
      IPRINT = 5
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
      CALL SECTION_SET_REAL('CLURAD',CLURAD,1.5D0,0)
C
      CALL INIT_MOD_TAUIJ_STAR(CLURAD)
      IF ( IBZINT.EQ.2 ) THEN
         CALL INIT_MOD_TAUIJ_KMESH
      ELSE
         WRITE (6,*) 'IBZINT should be 2  but not',IBZINT
         WRITE (6,*) 'for calculation of TAU(i,j)'
         CALL STOP_MESSAGE(ROUTINE,'IBZINT not allowed')
      END IF
C
C
      CALL RINIT((1+NLMAX)*NTMAX*2,CHILANLT)
C
      IPRINT = 0
C
C ======================================================================
C === orbital polarisation =============================================
C ======================================================================
C
      CHIPRINT = IPRINT
      CALL SECTION_SET_INTEGER('PRINT',CHIPRINT,9999,0)
C
      WRITE (6,99006) CHIPRINT
C
C ======================================================================
C
C ======================================================================
C *** FORMAT OF TAU x TAU-FILE:
C    1)  0=FORMATTED   -> 1 file for tau and tautau formatted
C or 2)  1=UNFORMATTED -> 2 files: one for tau      formatted
C                                  one for tautau   UNformatted
C -----------------------------------------------------------------
C
      TAUFIL = DATSET(1:LDATSET)//'.tau'
      LTAUFIL = LDATSET + 4
C
      FOUND_SECTION = .FALSE.
      IF ( RDTAU ) THEN
         OPEN (9,FILE=TAUFIL(1:LTAUFIL),STATUS='old',FORM='formatted')
         REWIND 9
      END IF
C
      WRITE (6,99009) PATH(IGRID(1))
      WRITE (6,99007) 'energies  :  ',NETAB(1)
      WRITE (6,99008) 'E_min     :  ',EMIN
      WRITE (6,99008) 'E_max     :  ',EMAX
      IF ( IGRID(1).NE.5 ) WRITE (6,99008) 'Im(E)     :  ',EIMAG
      WRITE (6,99008) 'E_Fermi   :  ',EFERMI
C
      CALL EPATH(IGRID(1),EMIN,EMAX,EIMAG,NETAB(1),ETAB(1,1),WETAB(1,1),
     &           EILOW,IPRINT,NEMAX)
C
      CALL RINIT(NRMAX*NTMAX,RHOCHRC)
      CALL RINIT(NRMAX*NTMAX,RHOSPNC)
C
C ======================================================================
C                      run ssite for E_F
C ======================================================================
      ERYD = DCMPLX(EFERMI,0.0D0)
C
      CALL SSITE(0,0,IFILCBWF,.TRUE.,.TRUE.,ERYD,P,IPRINT,NKM,TSST,MSST,
     &           SSST,MEZZ,MEZJ,ORBPOL)
C
C ======================================================================
C
      IF ( RDTAU ) THEN
         REWIND 9
C
         READ (9,'(//,10X,I3,5X,I3,6X,I2)',END=200,ERR=100) NT9,NQ9,IFMT
         WRITE (IWRI,*) ' IFMT=',IFMT
         IF ( NQ9.NE.NQ .OR. NT9.GT.NT ) THEN
            WRITE (IWRI,*) 'TAU-FILE inconsistent '
            WRITE (IWRI,*) ' NQ NQ9 ',NQ,NQ9,' NT NT9 ',NT,NT9
            CALL STOP_MESSAGE(ROUTINE,'reading the head of taufile!')
         END IF
         DO I = 1,(NT9+NQ9)
            READ (9,'(A80)',END=200,ERR=100) LINE
            WRITE (IWRI,'(A80)') LINE
         END DO
      END IF
C
      IF ( WRTAU .AND. (MPI_ID.EQ.0) ) THEN
         CALL INPUT_FIND_SECTION('ENERGY',0)
         OPEN (9,FILE=TAUFIL(1:LTAUFIL),FORM='formatted')
C
         REWIND 9
C
         WRITE (9,99011) CMDUC(3:LEN_TRIM(CMDUC)),NT,NQ
         WRITE (9,'(5X,'' IT ='',I3,'' IQ ='',I3,'': '',A)')
     &          (IT,IQAT(1,IT),TXT_T(IT),IT=1,NT)
         DO IQ = 1,NQ
            WRITE (9,'(''site'',I3,'' of '',A)') IQ,SYSTEM(1:LSYSTEM)
         END DO
C
      END IF
C
      NCPAFAIL = 0
C
      CALL CPU_TIME(TIME1)
C
      IF ( MPI ) CALL DRV_MPI_BARRIER
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      DO IE = 1,NETAB(1)
         IESORT(IE) = IE
         IPROCE(IE) = 0
      END DO
C
      IF ( NPROCS.GT.1 ) THEN
         IPROC = 0
         INC = 1
         DO JE = 1,NETAB(1) - 1
            IE = IESORT(JE)
            IPROC = IPROC + INC
            IF ( IPROC.EQ.NPROCS ) THEN
               IPROC = NPROCS - 1
               INC = -1
            ELSE IF ( IPROC.EQ.-1 ) THEN
               IPROC = 0
               INC = 1
            END IF
            IPROCE(IE) = IPROC
         END DO
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop    START
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
      WRITE (IWRI,99010) RDTAU,WRTAU
C
      DO IE = 1,NETAB(1)
         WRITE (6,99015) IE,ERYD
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         IF ( MPI_ID.EQ.IPROCE(IE) ) THEN
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
            IECURR = IE
C
            ERYD = ETAB(IE,1)
C
            ICPAFLAG = 0
            CPACHNG = 0.0D0
C
            CALL CINIT(NKMMAX*NKMMAX*NTMAX,TAUT)
C
C ===================================== solve SS - differential equation
C
            CALL SSITE(1,1,IFILCBWF,.TRUE.,.TRUE.,ERYD,P,IPRINT,NKM,
     &                 TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
            CALL MSSINIT(TSST,MSST,TSSQ,MSSQ)
C
C ======================================================================
C                                                  START:  CALCULATE TAU
            IF ( .NOT.RDTAU ) THEN
C
               CALL TAU_STD_KPOINTS(ICPAFLAG,CPACHNG,ERYD,P,IPRINT,
     &                              ITCPA,ICPACONV,CONC,NOQ,ITOQ,PHASK,
     &                              IECURR,NTMAX,TSST,MSST,TSSQ,MSSQ,
     &                              TAUQ)
C
               IF ( ICPAFLAG.NE.0 ) THEN
                  NCPAFAIL = NCPAFAIL + 1
                  CPACHTAB(NCPAFAIL) = CPACHNG
                  IECPAFAIL(NCPAFAIL) = IE
               END IF
C
               CALL PROJTAU(ICPAFLAG,CPACHNG,ERYD,MSST,MSSQ,TAUQ,TAUT)
C
            END IF
C                                                    END:  CALCULATE TAU
C=======================================================================
C                      calculate  CHI_Landau
C=======================================================================
C
C
C
            CALL CHILANCALC(IECURR,ERYD,P,CHILANLT,TSSQ,MSSQ,TSST,MSST,
     &                      TAUQ)
C
C=======================================================================
C
            IF ( IPRINT.GE.3 ) CALL DUMPTAU(IE,ERYD,IWRI,MSST,MSSQ,TAUT,
     &           TAUQ)
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      END DO
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop    END
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
      CALL CPU_TIME(TIME2)
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C ======================================================================
C                        print Landau susceptibility
C ======================================================================
      IF ( MPI ) THEN
C
         CALL DRV_MPI_BARRIER
C
         IDIMS = (1+NLMAX)*NTMAX*2
C
         CALL DRV_MPI_REDUCE_R(CHILANLT(0,1,2),RWORK(0,1,2),IDIMS)
C
         IF ( MPI_ID.EQ.0 ) THEN
C
            WRITE (6,99001) (TXT_L(IL),IL=1,NL)
            CHILAN = 0D0
            DO ITERM = 1,2
               DO IT = 1,NT
                  CHILANLT(0,IT,ITERM) = 0D0
C
                  DO IL = 1,NLT(IT)
                     CHILANLT(0,IT,ITERM) = CHILANLT(0,IT,ITERM)
     &                  + CHILANLT(IL,IT,ITERM)
                  END DO
C
                  CHILAN = CHILAN + NAT(IT)*CONC(IT)
     &                     *CHILANLT(0,IT,ITERM)
               END DO
C
               DO IT = 1,NT
                  WRITE (6,99002) ITERM,IT,TXT_T(IT),
     &                            (CHILANLT(IL,IT,ITERM)*CU,IL=0,NLT(IT)
     &                            )
C
               END DO
               IF ( NT.GT.1 ) WRITE (6,99003) CHILAN*CU
            END DO
            WRITE (6,99004)
C
         END IF
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      WRITE (IWRI,99012) TIME2 - TIME1
C
      IF ( NCPAFAIL.NE.0 ) THEN
         WRITE (IWRI,99013) CPATOL,NCPAFAIL,
     &                      (IECPAFAIL(I),DREAL(ETAB(IECPAFAIL(I),1)),
     &                      CPACHTAB(I),I=1,NCPAFAIL)
         WRITE (IWRI,'(1X,79(''*''),/)')
      ELSE IF ( NCPA.NE.0 ) THEN
         WRITE (IWRI,99014)
      END IF
C
C ======================================================================
C
      WRITE (6,*) '          <CHILANDAU> - run completed'
      STOP
C
 100  CONTINUE
      WRITE (IWRI,*) ' IECURR= ',IECURR,' (IE= ',IE,')'
      CALL STOP_MESSAGE(ROUTINE,'READING FILE 9')
 200  CONTINUE
      WRITE (IWRI,*) ' IECURR= ',IECURR,' (IE= ',IE,')'
      CALL STOP_MESSAGE(ROUTINE,'END OF FILE 9 REACHED')
C
99001 FORMAT (/,1X,79('='),//,10X,
     &        'Landau magnetic susceptibility  in [10^(-6) cm^3/mol]',
     &        //,1X,79('='),//,10X,'type',20X,'sum',18X,5(A,:,20X))
99002 FORMAT ('term',5X,I5,//,10X,I2,1X,A,6F20.10)
99003 FORMAT (10X,'total  ',6F11.4)
99004 FORMAT (/,1X,79('='),/)
99005 FORMAT (///,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*  **     *      *        **   *     * *****    **   *    *  *'
     &  ,/,10X,
     &  '*    *   *       *       *  *  **    * *    *  *  *  *    *  *'
     &  ,/,10X,
     &  '*     * *        *      *    * * *   * *    * *    * *    *  *'
     &  ,/,10X,
     &  '*      *     *** *      ****** *  *  * *    * ****** *    *  *'
     &  ,/,10X,
     &  '*     * *        *      *    * *   * * *    * *    * *    *  *'
     &  ,/,10X,
     &  '*    *   *       *      *    * *    ** *    * *    * *    *  *'
     &  ,/,10X,
     &  '*   *     **     ****** *    * *     * *****  *    *  ****   *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99006 FORMAT (1X,79('*'),/,24X,'settings for CHI - calculation ',/,1X,
     &        79('*'),//,5X,'print level            ',I3,/)
99007 FORMAT (10X,A,5I10)
99008 FORMAT (10X,A,5F10.6)
99009 FORMAT (/,5X,'energy path:  ',A)
99010 FORMAT (/,10X,'RDTAU     =',L2,'  WRTAU     =',L2)
99011 FORMAT ('ENERGY ',A,/,1X,79('*'),/,5X,' NT =',I3,' NQ =',I3,
     &        ' FMT = 3')
99012 FORMAT (/,' CPU - TIME  USED      ',F10.3,/)
99013 FORMAT (/,1X,79('*'),/,' tolerance for CPA-cycle:',F15.7,/,
     &        ' CPA not converged for',I3,' energies:',/,
     &        3(' E:',I3,F7.4,' C:',F8.5,:,2X))
99014 FORMAT (/,1X,79('*'),/,25X,'no problems with','  CPA-cycle !',/,
     &        1X,79('*'),/)
99015 FORMAT (/,' IE =',I3,' ERYD =',2F15.10,/)
      END
