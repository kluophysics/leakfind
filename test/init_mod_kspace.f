C*==init_mod_kspace.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
      SUBROUTINE INIT_MOD_KSPACE(INITELOOP,MOL,KMROT,ITEST,NKTABMAX,
     &                           NELMTMAX,SNTAUUVMAX)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:PROGNAME,TASK,KKRMODE
      USE MOD_LATTICE,ONLY:BRAVAIS,BOA,COA,BBAS,BBAS_L,BBAS_R,BBAS_I,
     &    SUB_SYSTEM,SYSTEM_DIMENSION
      USE MOD_KSPACE,ONLY:NKTAB,IBZINT,KTAB,WKTAB,KTET,IKCTET,NTETS,
     &    NKTET,NELMT,ITBZ,JTBZ,QTBZ,NTAUUV,WTAUUV,UTAUUV,VTAUUV,NGBAD,
     &    GBAD,NPTMIN,NPTMAX,NKPTS0,NZOOM,WKSUM,NKTABINP
      USE MOD_SYMMETRY,ONLY:NSYM,NWEDGE,IWEDGEROT,NSYMCRYSYS,MROTK,
     &    SYMACCEPTED,NO_SYMMETRY_LINRESP
      USE MOD_FILES,ONLY:IPRINT,IOTMP,DATSET,LDATSET,WRTAU,RDTAU,
     &    WRTAUMQ,RDTAUMQ,FOUND_SECTION,FOUND_INTEGER,FOUND_STRING
      USE MOD_MPI,ONLY:MPI_KLOOP,MPI_ID,NPROCS
      USE MOD_SIG,ONLY:NKTABSIG
      USE MOD_CPA,ONLY:USENLCPA
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='INIT_MOD_KSPACE')
C
C Dummy arguments
C
      LOGICAL INITELOOP,MOL
      INTEGER ITEST,KMROT,NELMTMAX,NKTABMAX,SNTAUUVMAX
C
C Local variables
C
      LOGICAL FOUND,INITIALIZE,LCHECKK,UDT
      INTEGER I,IA_ERR,IDIMS,IK,IKBOT,IKTOP,IPROC,IX,NFTET,NKTABMAXX,
     &        NKTABTMP,NKTET0
      REAL*8 KTABTMP(:,:),WKTABTMP(:)
      CHARACTER*10 STR10
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
C
      ALLOCATABLE KTABTMP,WKTABTMP
C
C      LCHECKK = .TRUE.
      LCHECKK = .FALSE.
C
      IF ( ALLOCATED(KTAB) ) DEALLOCATE (KTAB)
      IF ( ALLOCATED(WKTAB) ) DEALLOCATE (WKTAB)
      IF ( ALLOCATED(IKCTET) ) DEALLOCATE (IKCTET)
      IF ( ALLOCATED(KTET) ) DEALLOCATE (KTET)
      IF ( ALLOCATED(ITBZ) ) DEALLOCATE (ITBZ)
      IF ( ALLOCATED(JTBZ) ) DEALLOCATE (JTBZ)
      IF ( ALLOCATED(QTBZ) ) DEALLOCATE (QTBZ)
      IF ( ALLOCATED(NTAUUV) ) DEALLOCATE (NTAUUV)
      IF ( ALLOCATED(UTAUUV) ) DEALLOCATE (UTAUUV)
      IF ( ALLOCATED(VTAUUV) ) DEALLOCATE (VTAUUV)
      IF ( ALLOCATED(WTAUUV) ) DEALLOCATE (WTAUUV)
      IF ( ALLOCATED(GBAD) ) DEALLOCATE (GBAD)
C
C=======================================================================
C                        DUMMY allocation
C=======================================================================
      IF ( .NOT.INITELOOP ) THEN
         ALLOCATE (KTAB(3,1),WKTAB(1),IKCTET(1,1),KTET(1,1))
         NKTAB = 1
         ALLOCATE (ITBZ(1),JTBZ(1),QTBZ(1),NTAUUV(1))
         ALLOCATE (UTAUUV(1),VTAUUV(1),WTAUUV(1))
         NELMT = 1
         ALLOCATE (GBAD(3,1))
         ITBZ = 999999
         JTBZ = 999999
         QTBZ = 999999
         NGBAD = 999999
         GBAD = 999999D0
         RETURN
      END IF
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
      IF ( INITIALIZE ) THEN
C
         CALL INPUT_FIND_SECTION('TAU',0)
C
C --------------- allow to use the old section name KMESH instead of TAU
         IF ( .NOT.FOUND_SECTION ) CALL INPUT_FIND_SECTION('KMESH',0)
C
         IF ( FOUND_SECTION ) THEN
            IBZINT = -999999
            CALL SECTION_SET_STRING('BZINT',STR10,'9999',0)
C
            IF ( .NOT.FOUND_STRING ) THEN
C
               CALL SECTION_FIND_KEYWORD('CLUSTER',FOUND)
               IF ( FOUND ) IBZINT = 0
               CALL SECTION_FIND_KEYWORD('MOL',MOL)
               IF ( MOL ) IBZINT = 0
C
            ELSE IF ( STR10(1:7).EQ.'CLUSTER' ) THEN
C
               IBZINT = 0
C
            ELSE IF ( STR10(1:4).EQ.'WEYL' ) THEN
C
               IBZINT = 1
               CALL SECTION_SET_INTEGER('NKMIN',NPTMIN,9999,0)
               UDT = FOUND_INTEGER
               CALL SECTION_SET_INTEGER('NKMAX',NPTMAX,9999,0)
               IF ( .NOT.UDT .AND. .NOT.FOUND_INTEGER ) THEN
                  CALL SECTION_SET_INTEGER('NKTAB',NKTAB,9999,0)
                  IF ( FOUND_INTEGER ) THEN
                     NPTMIN = NINT(NKTAB*0.75)
                     NPTMAX = NKTAB
                  END IF
               END IF
C
            ELSE IF ( (STR10(1:6).EQ.'POINTS') .OR. (STR10(1:2).EQ.'SP')
     &                ) THEN
C
               IBZINT = 2
               CALL SECTION_SET_INTEGER('NKTAB',NKTAB,9999,0)
C
            ELSE IF ( STR10(1:3).EQ.'TET' .OR. STR10(1:4).EQ.'ZOOM' )
     &                THEN
C
               IF ( STR10(1:3).EQ.'TET' ) THEN
                  IBZINT = 3
               ELSE
                  IBZINT = 4
                  CALL SECTION_SET_INTEGER('NZOOM',NZOOM,9999,0)
               END IF
               NFTET = 0
               CALL SECTION_SET_INTEGER('NF',NFTET,9999,0)
               IF ( FOUND_INTEGER ) THEN
                  NKTET0 = 0
               ELSE
                  NKTET0 = 500
                  CALL SECTION_SET_INTEGER('NKTAB',NKTET0,9999,0)
               END IF
C
            ELSE
               WRITE (6,FMT='(/,5X,''BZINT = '',A,''   ??? '') ')
               CALL STOP_MESSAGE(ROUTINE,
     &                           'setting for BZINT not allowed')
            END IF
C
            IF ( IBZINT.EQ.-999999 )
     &            CALL STOP_MESSAGE(ROUTINE,'IBZINT not set')
C
         END IF
C
         IF ( (ITEST.EQ.2) .OR. (ITEST.EQ.3) ) IBZINT = -1
C
         IF ( IBZINT.NE.0 ) THEN
            WRITE (6,99001) NWEDGE,(IWEDGEROT(I),I=1,NWEDGE)
         ELSE
            WRTAU = .TRUE.
            RDTAU = .FALSE.
            WRTAUMQ = .FALSE.
            RDTAUMQ = .FALSE.
C        NOWRDOS = .TRUE.
            NKTAB = 0
            WRITE (6,99001)
            WRITE (6,99002)
         END IF
C------------------------------------------- save nktab value from input
         NKTABINP = NKTAB
C
         INITIALIZE = .FALSE.
C
C--------------- in sigma calc we will generate k-mesh for every E-point
         IF ( TASK.EQ.'SIGMA' .AND. .NOT.USENLCPA ) RETURN
C
      END IF
C=======================================================================
C                    END  INITIALIZE
C=======================================================================
C
C-----------------------------------------------------------------------
C                     WEYL POINT SAMPLING METHOD
C-----------------------------------------------------------------------
C
      IF ( IBZINT.EQ.1 ) THEN
         NKTAB = NPTMIN
         WRITE (6,99003) NPTMIN,NPTMAX
      END IF
C
C-----------------------------------------------------------------------
C                    REGULAR POINT SAMPLING METHOD
C-----------------------------------------------------------------------
C
      IF ( IBZINT.EQ.2 ) THEN
C
         WRITE (6,99004)
C
C=======================================================================
C                            3D system
C=======================================================================
C
         IF ( SYSTEM_DIMENSION(1:2).EQ.'3D' ) THEN
C
            NKPTS0 = NKTAB*NSYMCRYSYS
C
            IF ( TASK.EQ.'SIGMA' .AND. .NOT.USENLCPA )
     &           NKPTS0 = NKTABSIG*NSYMCRYSYS
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
            IF ( MPI_KLOOP .AND. 
c modified by XJQ: close the distribution for sigma calculation
c
c     &           ((TASK.EQ.'SIGMA' .OR. TASK.EQ.'PHONONS') .AND. 
     &           (TASK.EQ.'PHONONS' .AND. 
     &           KKRMODE(1:6).NE.'TB-KKR') ) THEN
c end-mod-xjq
C
C
C           + this part is only used for SIGMA_SPIN calcs at the moment
C             in case one energy point (E_F) is used
C           + for paralell run this ensures that only the root process
C             reduces the number of k-points to the IBZ, then the
C             irreducible k-points are  broadcasted to all processes
C             where every process gets only his share of points
C             (optimized for large number of k-points)
C           + this is of particular use when running a parallel
C             job on a system with shared memory
C
               IF ( MPI_ID.EQ.0 ) THEN
C
                  CALL KMESHS(IOTMP,NSYM,SYMACCEPTED,MROTK,IPRINT,
     &                        NKPTS0,BBAS,DATSET,LDATSET)
C
                  REWIND (IOTMP)
                  READ (IOTMP) NKTAB
C
                  NKTABMAXX = NKTAB
                  NKTABMAX = NKTAB
                  ALLOCATE (KTABTMP(3,NKTAB),WKTABTMP(NKTAB),
     &                      STAT=IA_ERR)
                  IF ( IA_ERR.NE.0 )
     &                  CALL STOP_MESSAGE(ROUTINE,'ALLOC: KTABTMP')
C
                  READ (IOTMP) ((KTABTMP(IX,IK),IX=1,3),WKTABTMP(IK),
     &                         IK=1,NKTAB)
C
                  WKSUM = SUM(WKTABTMP(1:NKTAB))
C
                  IF ( LCHECKK ) THEN
                     DO IK = 1,NKTAB
                        WRITE (6,'(I5,4E20.10)') IK,
     &                         (KTABTMP(I,IK),I=1,3),WKTABTMP(IK)
                     END DO
                  END IF
C
               END IF
C
C------divide k-array in pieces and send them from master to slave procs
C
               CALL DRV_MPI_BCAST_I(0,NKTABMAXX,1)
               DO IPROC = 1,NPROCS - 1
C
                  NKTABTMP = NKTABMAXX/NPROCS
                  IKBOT = IPROC*NKTABTMP + 1
                  IKTOP = (IPROC+1)*NKTABTMP
C
                  IF ( IPROC.EQ.NPROCS-1 ) THEN
                     NKTABTMP = NKTABMAXX - IPROC*NKTABTMP
                     IKTOP = NKTABMAXX
                  END IF
C
                  IF ( MPI_ID.EQ.0 ) WRITE (6,99007) NKTABTMP,IPROC
C
                  IF ( LCHECKK ) WRITE (6,*) IKBOT,IKTOP
C
                  IF ( MPI_ID.EQ.IPROC ) THEN
C
                     NKTABMAX = NKTABTMP
                     NKTAB = NKTABTMP
                     ALLOCATE (KTAB(3,NKTABTMP),WKTAB(NKTABTMP),
     &                         STAT=IA_ERR)
                     IF ( IA_ERR.NE.0 )
     &                     CALL STOP_MESSAGE(ROUTINE,'ALLOC: KTAB (a)')
                  END IF
C
                  IDIMS = 3*NKTABTMP
C
                  IF ( MPI_ID.EQ.IPROC ) THEN
                     CALL DRV_MPI_RECV_R(KTAB,IDIMS,IPROC,IPROC)
                     CALL DRV_MPI_RECV_R(WKTAB,NKTABTMP,IPROC,
     &                  10000*IPROC)
                     WRITE (6,99008) IPROC,NKTABTMP
                  ELSE IF ( MPI_ID.EQ.0 ) THEN
                     CALL DRV_MPI_RECV_R(KTABTMP(1,IKBOT),IDIMS,IPROC,
     &                  IPROC)
                     CALL DRV_MPI_RECV_R(WKTABTMP(IKBOT),NKTABTMP,IPROC,
     &                  10000*IPROC)
                  END IF
C
                  IF ( LCHECKK .AND. (MPI_ID.EQ.IPROC) ) THEN
                     NKTABMAX = NKTABTMP
                     DO IK = 1,NKTABMAX
                        WRITE (6,'(I5,4E20.10)') IK,(KTAB(I,IK),I=1,3),
     &                         WKTAB(IK)
                     END DO
                  END IF
C
               END DO
C
C------------------------------ master itself gets its array of k-points
               IF ( MPI_ID.EQ.0 ) THEN
C
                  NKTABTMP = NKTABMAXX/NPROCS
                  ALLOCATE (KTAB(3,NKTABTMP),WKTAB(NKTABTMP),
     &                      STAT=IA_ERR)
                  IF ( IA_ERR.NE.0 )
     &                  CALL STOP_MESSAGE(ROUTINE,'ALLOC: KTAB (b)')
C
                  KTAB(1:3,1:NKTABTMP) = KTABTMP(1:3,1:NKTABTMP)
                  WKTAB(1:NKTABTMP) = WKTABTMP(1:NKTABTMP)
C
                  DEALLOCATE (KTABTMP,WKTABTMP)
C
                  NKTABMAX = NKTABTMP
                  NKTAB = NKTABTMP
                  IF ( MPI_ID.EQ.0 ) WRITE (6,99007) NKTABTMP,0
C
                  CLOSE (IOTMP)
C
                  IF ( LCHECKK ) THEN
                     WRITE (6,*) 'K-points master'
                     DO IK = 1,NKTABTMP
                        WRITE (6,'(i5,4e20.10)') IK,(KTAB(I,IK),I=1,3),
     &                         WKTAB(IK)
                     END DO
                  END IF
               END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
            ELSE
C
               CALL KMESHS(IOTMP,NSYM,SYMACCEPTED,MROTK,IPRINT,NKPTS0,
     &                     BBAS,DATSET,LDATSET)
C
               REWIND (IOTMP)
               READ (IOTMP) NKTAB
               NKTABMAX = NKTAB
               ALLOCATE (KTAB(3,NKTABMAX),WKTAB(NKTABMAX),STAT=IA_ERR)
               IF ( IA_ERR.NE.0 )
     &               CALL STOP_MESSAGE(ROUTINE,'ALLOC: KTAB (c)')
               READ (IOTMP) ((KTAB(IX,IK),IX=1,3),WKTAB(IK),IK=1,NKTAB)
               CLOSE (IOTMP)
C
               WKSUM = SUM(WKTAB(1:NKTAB))
C
            END IF
C
C=======================================================================
C                            2D-layered system
C=======================================================================
C
         ELSE IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' ) THEN
C
            IF ( SUB_SYSTEM(1:6).EQ.'L-BULK' ) THEN
C
               NKPTS0 = NKTAB*NSYMCRYSYS
C
               CALL KMESHS(IOTMP,NSYM,SYMACCEPTED,MROTK,IPRINT,NKPTS0,
     &                     BBAS_L,DATSET,LDATSET)
C
               REWIND (IOTMP)
               READ (IOTMP) NKTAB
               NKTABMAX = NKTAB
               ALLOCATE (KTAB(3,NKTABMAX),WKTAB(NKTABMAX),STAT=IA_ERR)
C
               IF ( IA_ERR.NE.0 )
     &               CALL STOP_MESSAGE(ROUTINE,'ALLOC: KTAB (d)')
               READ (IOTMP) ((KTAB(IX,IK),IX=1,3),WKTAB(IK),IK=1,NKTAB)
               CLOSE (IOTMP)
C
            ELSE IF ( SUB_SYSTEM(1:6).EQ.'R-BULK' ) THEN
               NKPTS0 = NKTAB*NSYMCRYSYS
C
               CALL KMESHS(IOTMP,NSYM,SYMACCEPTED,MROTK,IPRINT,NKPTS0,
     &                     BBAS_R,DATSET,LDATSET)
C
               REWIND (IOTMP)
               READ (IOTMP) NKTAB
               NKTABMAX = NKTAB
               ALLOCATE (KTAB(3,NKTABMAX),WKTAB(NKTABMAX),STAT=IA_ERR)
C
               IF ( IA_ERR.NE.0 )
     &               CALL STOP_MESSAGE(ROUTINE,'ALLOC: KTAB (e)')
               READ (IOTMP) ((KTAB(IX,IK),IX=1,3),WKTAB(IK),IK=1,NKTAB)
               CLOSE (IOTMP)
C
            ELSE IF ( SUB_SYSTEM(1:6).EQ.'I-ZONE' ) THEN
C
               IF ( TASK.EQ.'SIGMA' ) THEN
                  NKPTS0 = NKTABSIG*NSYMCRYSYS
               ELSE
                  NKPTS0 = NKTAB*NSYMCRYSYS
               END IF
C
               CALL KMESHS(IOTMP,NSYM,SYMACCEPTED,MROTK,IPRINT,NKPTS0,
     &                     BBAS_I,DATSET,LDATSET)
C
               REWIND (IOTMP)
               READ (IOTMP) NKTAB
               NKTABMAX = NKTAB
               ALLOCATE (KTAB(3,NKTABMAX),WKTAB(NKTABMAX),STAT=IA_ERR)
C
               IF ( IA_ERR.NE.0 )
     &               CALL STOP_MESSAGE(ROUTINE,'ALLOC: KTAB (f)')
               READ (IOTMP) ((KTAB(IX,IK),IX=1,3),WKTAB(IK),IK=1,NKTAB)
               CLOSE (IOTMP)
            END IF
C
         END IF
      END IF
C
C-----------------------------------------------------------------------
C                        TETRAHEDRON METHOD
C-----------------------------------------------------------------------
C
      IF ( IBZINT.EQ.3 .OR. IBZINT.EQ.4 ) THEN
C
         CALL TETARRSIZ(BRAVAIS,NKTET0,NFTET,NKTET,NTETS,BBAS)
C
         ALLOCATE (IKCTET(4,NTETS),KTET(3,NKTET),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: IKCTET')
C
         CALL TETGEN(BRAVAIS,BOA,COA,NFTET,KTET,NKTET,IKCTET,NTETS)
C
         WRITE (6,99005) NFTET
         IF ( KMROT.NE.0 ) CALL STOP_MESSAGE(ROUTINE,
     &        'rotation of magnetisation not allowed for IBZINT=3')
C
         CALL STRGBAD(KTET,NKTET,MROTK,IWEDGEROT,NWEDGE)
C
      ELSE
C
         ALLOCATE (GBAD(3,1))
         NGBAD = 999999
         GBAD = 999999D0
      END IF
C
C=======================================================================
C                   REAL SPACE  type calculation
C=======================================================================
C
      IF ( IBZINT.LE.0 ) THEN
         IF ( ITEST.NE.2 ) KKRMODE = 'REAL-SPACE'
         ALLOCATE (KTAB(3,1),WKTAB(1))
         NKTAB = 1
      END IF
C
C=======================================================================
C
      IF ( IBZINT.NE.3 .AND. IBZINT.NE.4 ) THEN
         ALLOCATE (IKCTET(1,1),KTET(1,1))
         KTET = 999999D0
         IKCTET = 999999
      END IF
C
C=======================================================================
C                           form of   TAU
C=======================================================================
C
c modified by XJQ: probably not needed for sigma calculation for non-symmetry cluster
c it is safe to turn off symtau when NO_SYMMETRY_LINRESP
      if(NO_SYMMETRY_LINRESP) return
c end-mod-xjq
      IF ( IBZINT.NE.0 ) THEN
C
         IF ( PROGNAME(4:6).EQ.'CHI' .OR. PROGNAME(4:6).EQ.'OPM' .OR. 
     &        IBZINT.EQ.3 ) THEN
C
            CALL SYMTAU(ITEST)
C
            REWIND (IOTMP)
            READ (IOTMP) NELMTMAX,SNTAUUVMAX
            WRITE (6,99006) ' TAU-elements     NELMT: ',NELMT,NELMTMAX
         END IF
C
         NELMT = NELMTMAX
C
      END IF
C
      ALLOCATE (ITBZ(NELMTMAX),JTBZ(NELMTMAX))
      ALLOCATE (QTBZ(NELMTMAX),NTAUUV(NELMTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: QTBZ')
C
      ALLOCATE (UTAUUV(SNTAUUVMAX))
      ALLOCATE (VTAUUV(SNTAUUVMAX),WTAUUV(SNTAUUVMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: WTAUUV')
C
      IF ( NELMTMAX.GT.1 .OR. SNTAUUVMAX.GT.1 ) THEN
         READ (IOTMP) (ITBZ(I),JTBZ(I),NTAUUV(I),QTBZ(I),I=1,NELMTMAX)
         READ (IOTMP) (UTAUUV(I),VTAUUV(I),WTAUUV(I),I=1,SNTAUUVMAX)
         CLOSE (IOTMP)
      ELSE
         ITBZ = 999999
         JTBZ = 999999
         QTBZ = 999999
      END IF
C
C=======================================================================
C
99001 FORMAT (/,1X,79('t'),/,:,/,10X,
     &        'parameters for k-space integration:',/,10X,'NWEDGE:',I4,
     &        '  wedges:',25I3)
99002 FORMAT (10X,'A real-space cluster calculation will be performed')
99003 FORMAT (10X,'WEYL - k-mesh used for point sampling',/,10X,
     &        'number of k-points/wedge:',I7,' --',I7,
     &        '   according to Im(E)')
99004 FORMAT (10X,'regular k-mesh used for point sampling')
99005 FORMAT (10X,'tetrahedron method used with mesh ','NF: ',I7)
99006 FORMAT (10X,A,10I5)
99007 FORMAT (10X,'Master sending ',I13,' kpoints to process number',i5)
99008 FORMAT (/,10X,'Slave process number',i5,' -  obtained',i13,
     &        ' kpoints from master')
      END
