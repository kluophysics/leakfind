C*==potrd1.f    processed by SPAG 6.70Rc at 08:47 on  8 Mar 2017
      SUBROUTINE POTRD1(IPRINT,POTFMTINP,SCFSTATUS,NQ,NTLIM,NM)
C   ********************************************************************
C   *                                                                  *
C   *    reading the potential file   STEP 1                           *
C   *                                                                  *
C   *    get  NQ,  NM and NTLIM; i.e.  an upper limit for NT           *
C   *    to set preliminary array dimensions                           *
C   *                                                                  *
C   *    POTFMTINP gives potential file format                         *
C   *    <= 5:  OLD SPR-KKR format                                     *
C   *    >= 6:  NEW common SPR-KKR / SPR-TB-KKR format                 *
C   *                                                                  *
C   *  - first the routine tries to read using the old SPR-KKR format  *
C   *  - next  it tries to use the new SPR-KKR / SPR-TB-KKR format     *
C   *  - finally it tries to import data from other sources            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:POTFIL,LPOTFIL,RDUMMY,IDUMMY
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IPRINT,NM,NQ,NTLIM,POTFMTINP
      CHARACTER*10 SCFSTATUS
C
C Local variables
C
      LOGICAL CHECK
      INTEGER I,IO,IOPOT,IPOS,IQ,IQIN,IWRI,LOOP,NOQ
      CHARACTER*80 LINE,LINUC
      REAL*8 QVEC(3)
      CHARACTER*10 STR10
C
C*** End of declarations rewritten by SPAG
C
      IF ( IPRINT.GE.3 ) THEN
         CHECK = .TRUE.
         IWRI = 6
      ELSE
         CHECK = .FALSE.
         IWRI = 0
      END IF
      IF ( IPRINT.GT.0 ) WRITE (6,99002) POTFIL(1:LPOTFIL)
C
      IOPOT = 4
C
      LOOP = 0
      LOOP = LOOP + 1
C
      REWIND IOPOT
C
      I = 0
 100  CONTINUE
      READ (IOPOT,'(A)',END=200) LINUC
      I = I + 1
      IF ( LINUC(1:6).NE.'FORMAT' .AND. I.LT.20 ) GOTO 100
C
      IF ( LINUC(1:6).EQ.'FORMAT' ) THEN
C
         IPOS = INDEX(LINUC,'VERSION')
C
C-------------------------------------------------------  POTFMTINP.LE.5
         IF ( IPOS.GE.1 ) THEN
            LINUC = LINUC((IPOS+7):80)
            READ (LINUC,*) POTFMTINP
         ELSE
            REWIND IOPOT
            CALL READKWINT(IOPOT,'FORMAT    ',POTFMTINP,0,IWRI,1)
         END IF
         IF ( POTFMTINP.LT.3 ) STOP 'in <POTRD1>:  POTFMTINP < 3'
         IF ( POTFMTINP.GT.8 ) STOP 'in <POTRD1>:  POTFMTINP > 8'
      ELSE
         RETURN
      END IF
C
      REWIND IOPOT
C
      IF ( POTFMTINP.LE.5 ) THEN
C
C=======================================================================
C                          SPRKKR FORMAT 5
C=======================================================================
C
         READ (IOPOT,'(A)') LINE
         IF ( CHECK ) WRITE (6,*) LINE
C
         IF ( LINE(1:10).EQ.'SPRKKR    ' ) THEN
C
            IF ( LINE(11:19).EQ.'SCF-start' ) THEN
               SCFSTATUS = 'START'
            ELSE
               SCFSTATUS = 'ITERATING'
            END IF
C
            READ (IOPOT,99004) LINE
C
            IF ( LINE(1:7).NE.'VERSION' ) RETURN
C
            CALL READKWINT(IOPOT,'NQ        ',NQ,0,IWRI,1)
            CALL READKWINT(IOPOT,'NM        ',NM,0,IWRI,1)
C
            CALL POSFIL(IOPOT,'SITES     ',IPOS,1)
C
            NTLIM = 0
            DO IQ = 1,NQ
               READ (IOPOT,99005) IQIN,QVEC,NOQ,(IDUMMY,RDUMMY,IO=1,NOQ)
               IF ( IPRINT.GT.2 ) WRITE (6,99005) IQIN,QVEC,NOQ,
     &              (IDUMMY,RDUMMY,IO=1,NOQ)
               NTLIM = NTLIM + NOQ
            END DO
         END IF
C
      ELSE
C
C=======================================================================
C                          SPRKKR FORMAT 6
C=======================================================================
C
         READ (IOPOT,*)
         READ (IOPOT,99006) STR10
         IF ( STR10(1:6).EQ.'HEADER' ) THEN
C
            CALL READKWSTR(IOPOT,'SCFSTATUS ',SCFSTATUS,'          ',10,
     &                     IWRI,1)
            REWIND IOPOT
            CALL READKWINT(IOPOT,'NQ        ',NQ,0,IWRI,1)
            CALL READKWINT(IOPOT,'NM        ',NM,0,IWRI,1)
C
            CALL POSFIL(IOPOT,'OCCUPATION',IPOS,1)
            READ (IOPOT,99006) LINE
C
            NTLIM = 0
C
            DO IQ = 1,NQ
               READ (IOPOT,*) IDUMMY,IDUMMY,IDUMMY,NOQ
               NTLIM = NTLIM + NOQ
            END DO
C
C=======================================================================
C                   TRY TO IMPORT POTENTIAL FROM OTHER SOURCES
C=======================================================================
C
         ELSE IF ( LOOP.EQ.1 ) THEN
            CALL POTFIT
C=======================================================================
            STOP 'in <POTRD1>:   trouble reading FORMAT'
         ELSE
            WRITE (6,99001)
            STOP
         END IF
C
      END IF
C=======================================================================
C
      IF ( IPRINT.GT.0 ) WRITE (6,99003) NQ,NTLIM,NM
      RETURN
C
 200  CONTINUE
      CALL STOP_MESSAGE('POTRD1','string FORMAT not found in pot file')
C
99001 FORMAT (//,2(1X,79('*'),/),33X,'<POTRD1>',/,2(1X,79('*'),/),10X,
     &        'all attempts to read the potential file FAILED',/,
     &        2(1X,79('*'),/),/)
99002 FORMAT (//,1X,79('*'),/,33X,'<POTRD1>',/,1X,79('*'),//,10X,
     &        'get variables       NQ,    NTLIM,    NM'/,10X,
     &        'from potential file ',A,/,10X,
     &        'to fix array sizes  NQMAX, NMMAX and NTLIM (temporarily)'
     &        /)
99003 FORMAT (10X,'NQ =',I4,4X,'NT <=',I3,4X,'NM =',I4)
99004 FORMAT (10X,A)
99005 FORMAT (I4,4X,2(F14.10,1X),F14.10,8X,I3,10X,20(I3,F6.3))
99006 FORMAT (A)
      END
C*==potrd2.f    processed by SPAG 6.70Rc at 08:47 on  8 Mar 2017
      SUBROUTINE POTRD2(POTFMTINP,IREL,BRAVAIS,QBAS,SYSTEM,LSYSTEM,NQ,
     &                  NT0,NM,Z_TAUX,NA_TAUX,IM_QAUX,IQ_ATAUX,NOQ,
     &                  IT_OQAUX,RMTRED0,SCFSTATUS,NONMAG,BREITINT,
     &                  ORBPOL,EXTFIELD,KMROT,QMTET,QMPHI,QMVEC,FULLPOT,
     &                  SPHERCELL,NQMAX,NTLIM,NMMAX,SYSTEM_DIMENSION,
     &                  SYSTEM_TYPE,ALAT,ABAS,ABAS_L,ABAS_R,NQ_L,NQ_R)
C   ********************************************************************
C   *                                                                  *
C   *    reading the potential file   STEP 2                           *
C   *                                                                  *
C   *  - only for standard SPRKKR format                               *
C   *  - read all quantities relevant to fix the number of types  NT   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IPRINT,POTFIL,LPOTFIL,RDUMMY,IDUMMY,NDUMMY,
     &    S10DUMMY
      USE MOD_SITES,ONLY:QMPHI_POTFIL,QMTET_POTFIL
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALAT
      INTEGER BRAVAIS,IREL,KMROT,LSYSTEM,NM,NMMAX,NQ,NQMAX,NQ_L,NQ_R,
     &        NT0,NTLIM,POTFMTINP
      LOGICAL BREITINT,EXTFIELD,FULLPOT,NONMAG,SPHERCELL
      CHARACTER*10 ORBPOL,SCFSTATUS,SYSTEM_DIMENSION,SYSTEM_TYPE
      CHARACTER*80 SYSTEM
      REAL*8 ABAS(3,3),ABAS_L(3,3),ABAS_R(3,3),QBAS(3,NQMAX),
     &       QMPHI(NQMAX),QMTET(NQMAX),QMVEC(3),RMTRED0(NMMAX)
      INTEGER IM_QAUX(NQMAX),IQ_ATAUX(NQMAX,NTLIM),IT_OQAUX(NTLIM,NQMAX)
     &        ,NA_TAUX(NTLIM),NOQ(NQMAX),Z_TAUX(NTLIM)
C
C Local variables
C
      REAL*8 BT(:),QMGAM(:),RMT(:),RWS(:),SWS
      LOGICAL CHECK
      CHARACTER*80 HEADER,LINE,TITLE
      INTEGER IA,IC,IM,IM_TAUX,IO,IOPOT,IPAN,IPOS,IQ,IR,IRSHFT,IRTOP,IT,
     &        IWRI,J,JRCRI(:),JRMT(:),JRNS1(:),JRWS(:)
      CHARACTER*12 RMESHTYPE
      CHARACTER*4 STR4
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE JRNS1,JRCRI,BT,JRWS,JRMT,RWS,RMT,QMGAM
C
      ALLOCATE (JRWS(NM),JRMT(NM),RWS(NM),RMT(NM),JRNS1(NM),JRCRI(NM))
      ALLOCATE (QMGAM(NQMAX))
      QMGAM(:) = 0D0
C
      CALL WRSUBNAMHEAD('POTRD2',6,' ',0)
      IF ( IPRINT.GT.0 ) WRITE (6,99001) POTFIL(1:LPOTFIL)
C
      IF ( IPRINT.GE.3 ) THEN
         CHECK = .TRUE.
         IWRI = 6
      ELSE
         CHECK = .FALSE.
         IWRI = 0
      END IF
C
      REWIND 4
      IOPOT = 4
C
C=======================================================================
C
C                             SPRKKR FORMAT 5
C
C=======================================================================
C
      IF ( POTFMTINP.LE.5 ) THEN
C
         READ (4,'(A)') LINE
         IF ( LINE(1:6).NE.'SPRKKR' ) STOP 'in <POTRD2> --> SPRKKR'
C
         FULLPOT = (LINE(11:24).EQ.'full potential')
         IF ( LINE(11:19).EQ.'SCF-start' ) THEN
            SCFSTATUS = 'START'
         ELSE
            SCFSTATUS = 'ITERATING'
         END IF
         WRITE (6,99005) LINE(1:LEN_TRIM(LINE))
         READ (4,99005) TITLE
         IF ( TITLE(1:7).NE.'VERSION' ) THEN
            POTFMTINP = 0
         ELSE
            POTFMTINP = ICHAR(TITLE(10:10)) - ICHAR('1') + 1
            READ (4,99005) TITLE
         END IF
         IF ( POTFMTINP.LT.3 ) STOP 'in <POTRD2> --> POTFMTINP'
         READ (4,99005) SYSTEM
         LSYSTEM = LEN_TRIM(SYSTEM)
         WRITE (6,99006) SYSTEM(1:LSYSTEM)
         READ (4,99007) NM
         READ (4,99007) NQ
         READ (4,99007) NT0
         IF ( NM.GT.NMMAX ) STOP 'in <POTRD2> --> NM'
         IF ( NQ.GT.NQMAX ) STOP 'in <POTRD2> --> NQ'
         IF ( NT0.GT.NTLIM ) STOP 'in <POTRD2> --> NT0'
C
         READ (4,*)
         READ (4,*)
         READ (4,*)
         READ (4,*)
         READ (4,*)
         READ (4,*)
         READ (4,*)
         READ (4,*)
         READ (4,*)
         READ (4,*)
         READ (4,*)
         IF ( POTFMTINP.GE.4 ) READ (4,99005) STR4
         READ (4,99010) BREITINT
         READ (4,99005) ORBPOL
         IF ( ORBPOL.EQ.'          ' ) ORBPOL = 'NONE'
         READ (4,99010) EXTFIELD
         READ (4,*)
         READ (4,*)
         READ (4,99007) KMROT
         READ (4,99009) (QMVEC(IC),IC=1,3)
         READ (4,*)
         READ (4,99008) STR4,BRAVAIS
         IF ( STR4.EQ.'LAT ' ) STOP 'in <POTRD2> --> found LAT '
         READ (4,*)
         READ (4,*)
         READ (4,*)
         READ (4,*)
         DO J = 1,3
            READ (4,99009) (ABAS(IC,J),IC=1,3)
         END DO
C --------------------------------------- read in radial mesh parameters
         READ (4,99005) STR4
         READ (4,99005) STR4
         READ (4,99005) RMESHTYPE
         IRTOP = 0
         DO IM = 1,NM
            READ (4,*)
            READ (4,*)
            READ (4,*)
            READ (4,99007) JRWS(IM)
            IRTOP = MAX(IRTOP,JRWS(IM))
            READ (4,*)
            IF ( RMESHTYPE.EQ.'EXPONENTIAL ' ) THEN
               READ (4,*)
               READ (4,*)
               IF ( FULLPOT ) THEN
                  READ (4,*)
                  READ (4,*)
                  READ (4,99007) NDUMMY
                  READ (4,99007) (IDUMMY,IPAN=1,NDUMMY)
                  READ (4,99010) SPHERCELL
               END IF
            ELSE IF ( RMESHTYPE.NE.'JUELICH     ' ) THEN
               WRITE (6,*) 'RMESHTYPE = ',RMESHTYPE
               STOP ' check RMESHTYPE in POTFILE '
            END IF
         END DO
C
         READ (4,'(/)')
         DO IQ = 1,NQ
            READ (4,99011) IDUMMY,(QBAS(IC,IQ),IC=1,3),NOQ(IQ),
     &                     (IT_OQAUX(IO,IQ),RDUMMY,IO=1,NOQ(IQ))
         END DO
         DO IQ = 1,NQ
            READ (4,99012) IDUMMY,QMTET(IQ),QMPHI(IQ)
         END DO
C
C--------------------------- do a final test and adjust KMROT eventually
C
         CALL MROT_CHECK_SETTINGS(NQ,IREL,QMPHI,QMTET,QMGAM,KMROT,QMVEC,
     &                            NQMAX)
C
C---------------------------------------------- IT-dependent information
C
         ALLOCATE (BT(IRTOP))
         IT = 0
 50      CONTINUE
         READ (4,'(A)',END=100) LINE
         IF ( LINE(1:10).NE.'TYPE      ' ) GOTO 50
         IT = IT + 1
         READ (4,*)
         READ (4,99007) Z_TAUX(IT)
         READ (4,99007) NA_TAUX(IT)
         READ (4,*)
         READ (4,99007) IM_TAUX
         IM = IM_TAUX
         READ (4,99007) (IQ_ATAUX(IA,IT),IA=1,NA_TAUX(IT))
         DO IA = 1,NA_TAUX(IT)
            IM_QAUX(IQ_ATAUX(IA,IT)) = IM
         END DO
         IF ( POTFMTINP.GE.4 ) READ (4,*)
         IF ( SCFSTATUS(1:5).NE.'START' ) THEN
            IRTOP = JRWS(IM)
            READ (4,99013) (RDUMMY,IR=1,IRTOP)
            READ (4,99013) (BT(IR),IR=1,IRTOP)
            DO IR = 1,IRTOP
               IF ( ABS(BT(IR)).GT.1D-7 ) NONMAG = .FALSE.
            END DO
         END IF
         IF ( IT.LT.NT0 ) GOTO 50
C
         IF ( SCFSTATUS(1:5).EQ.'START' ) NONMAG = .FALSE.
C
C=======================================================================
C
C                             SPRKKR FORMAT 6
C
C=======================================================================
C
      ELSE
C
         LINE = '          '
         REWIND IOPOT
         CALL READKWSTR(IOPOT,'HEADER    ',HEADER,LINE,80,IWRI,1)
         CALL READKWSTR(IOPOT,'TITLE     ',TITLE,LINE,80,IWRI,1)
         CALL READKWSTR(IOPOT,'SYSTEM    ',SYSTEM,LINE,80,IWRI,1)
         LSYSTEM = LEN_TRIM(SYSTEM)
C
         CALL READKWINT(IOPOT,'FORMAT    ',POTFMTINP,0,IWRI,1)
         IF ( POTFMTINP.LT.3 ) STOP 'in <POTRD2>:  POTFMTINP < 3'
         IF ( POTFMTINP.GT.8 ) STOP 'in <POTRD2>:  POTFMTINP > 8'
C
C ======================================================================
C
         CALL READKWINT(IOPOT,'NQ        ',NQ,0,IWRI,1)
         CALL READKWINT(IOPOT,'NT        ',NT0,0,IWRI,1)
         CALL READKWINT(IOPOT,'NM        ',NM,0,IWRI,1)
         CALL READKWINT(IOPOT,'IREL      ',IREL,3,IWRI,0)
C
         IF ( NM.GT.NMMAX ) STOP 'in <POTRD2> --> NM'
         IF ( NQ.GT.NQMAX ) STOP 'in <POTRD2> --> NQ'
         IF ( NT0.GT.NTLIM ) STOP 'in <POTRD2> --> NT0'
C
C ======================================================================
C
         CALL READKWSTR(IOPOT,'SCFSTATUS ',SCFSTATUS,'ITERATING ',10,
     &                  IWRI,1)
C
         CALL READKWLOG(IOPOT,'FULLPOT   ',FULLPOT,.FALSE.,IWRI,0)
         CALL READKWLOG(IOPOT,'SPHERCELL ',SPHERCELL,.FALSE.,IWRI,0)
C
         CALL READKWLOG(IOPOT,'BREITINT  ',BREITINT,.FALSE.,IWRI,0)
         CALL READKWSTR(IOPOT,'ORBPOL    ',ORBPOL,'NONE      ',10,IWRI,
     &                  0)
         CALL READKWLOG(IOPOT,'NONMAG    ',NONMAG,.FALSE.,IWRI,0)
         CALL READKWLOG(IOPOT,'EXTFIELD  ',EXTFIELD,.FALSE.,IWRI,0)
C
C ======================================================================
C
         CALL READKWSTR(IOPOT,'SYSDIM    ',SYSTEM_DIMENSION,
     &                  '3D        ',10,IWRI,1)
         CALL READKWSTR(IOPOT,'SYSTYPE   ',SYSTEM_TYPE,'BULK      ',10,
     &                  IWRI,1)
         CALL READKWINT(IOPOT,'BRAVAIS   ',BRAVAIS,0,IWRI,0)
         CALL READKWREAL(IOPOT,'ALAT      ',ALAT,0D0,IWRI,1)
         CALL READKWRARR(IOPOT,'A(1)      ',ABAS(1,1),3,0D0,IWRI,1)
         CALL READKWRARR(IOPOT,'A(2)      ',ABAS(1,2),3,0D0,IWRI,1)
         CALL READKWRARR(IOPOT,'A(3)      ',ABAS(1,3),3,0D0,IWRI,1)
C
         IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' ) THEN
            ABAS_L(1:3,1:2) = ABAS(1:3,1:2)
            ABAS_R(1:3,1:2) = ABAS(1:3,1:2)
C
            CALL READKWINT(IOPOT,'NQ_L      ',NQ_L,0,IWRI,1)
            CALL READKWRARR(IOPOT,'A_L(3)    ',ABAS_L(1,3),3,0D0,IWRI,1)
C
            CALL READKWINT(IOPOT,'NQ_R      ',NQ_R,0,IWRI,1)
            CALL READKWRARR(IOPOT,'A_R(3)    ',ABAS_R(1,3),3,0D0,IWRI,1)
         END IF
C
C ======================================================================
         CALL POSFIL(IOPOT,'BASSCALE   ',IPOS,1)
         READ (IOPOT,'(A80)') LINE
         IF ( CHECK ) WRITE (6,'(A)') LINE(1:LEN_TRIM(LINE))
         DO IQ = 1,NQ
            READ (IOPOT,*) IDUMMY,(QBAS(IC,IQ),IC=1,3)
            IF ( CHECK ) WRITE (6,99002) IDUMMY,(QBAS(IC,IQ),IC=1,3)
         END DO
C ======================================================================
         CALL POSFIL(IOPOT,'OCCUPATION',IPOS,1)
         READ (IOPOT,'(A80)') LINE
         IF ( CHECK ) WRITE (6,'(A)') LINE(1:LEN_TRIM(LINE))
         DO IQ = 1,NQ
            READ (IOPOT,*) IDUMMY,IDUMMY,IM_QAUX(IQ),NOQ(IQ),
     &                     (IT_OQAUX(IO,IQ),RDUMMY,IO=1,NOQ(IQ))
            IF ( CHECK ) WRITE (6,99003) IQ,IQ,IM_QAUX(IQ),NOQ(IQ),
     &                          (IT_OQAUX(IO,IQ),RDUMMY,IO=1,NOQ(IQ))
         END DO
C ======================================================================
         CALL READKWINT(IOPOT,'KMROT     ',KMROT,0,IWRI,0)
         READ (IOPOT,*) S10DUMMY,(QMVEC(IC),IC=1,3)
         IF ( CHECK ) WRITE (6,99014) S10DUMMY,(QMVEC(IC),IC=1,3)
C
         READ (IOPOT,'(A80)') LINE
         IF ( CHECK ) WRITE (6,'(A)') LINE(1:LEN_TRIM(LINE))
         DO IQ = 1,NQ
            READ (IOPOT,*) IDUMMY,QMTET(IQ),QMPHI(IQ)
            IF ( CHECK ) WRITE (6,99004) IQ,QMTET(IQ),QMPHI(IQ)
         END DO
C
C--------------------------- do a final test and adjust KMROT eventually
C
         CALL MROT_CHECK_SETTINGS(NQ,IREL,QMPHI,QMTET,QMGAM,KMROT,QMVEC,
     &                            NQMAX)
C
C ======================================================================
C
         CALL POSFIL(IOPOT,'MESH INFOR',IPOS,1)
         CALL READKWSTR(IOPOT,'MESH-TYPE ',RMESHTYPE,'            ',12,
     &                  IWRI,1)
C ----------------------------------------------------------------------
         IRSHFT = 0
         IRTOP = 0
         IF ( RMESHTYPE.EQ.'EXPONENTIAL ' ) THEN
            READ (IOPOT,*)
            DO IM = 1,NM
               READ (IOPOT,*) IDUMMY,RDUMMY,RDUMMY,JRMT(IM),RMT(IM),
     &                        JRWS(IM),RWS(IM)
               IRTOP = MAX(IRTOP,JRWS(IM))
            END DO
         ELSE IF ( RMESHTYPE.EQ.'JUELICH     ' ) THEN
C ----------------------------------------------------------------------
            IRSHFT = 1
            READ (IOPOT,*)
            DO IM = 1,NM
               READ (IOPOT,*) IDUMMY,RDUMMY,RDUMMY,JRMT(IM),RMT(IM),
     &                        JRWS(IM),RWS(IM)
               JRWS(IM) = JRWS(IM) - IRSHFT
               IRTOP = MAX(IRTOP,JRWS(IM))
            END DO
         ELSE
            WRITE (6,*) 'MESH-TYPE read in:',RMESHTYPE
            STOP 'in <POTRD>  mesh-type  not allowed'
         END IF
C ----------------------------------------------------------------------
         IF ( FULLPOT ) THEN
            READ (IOPOT,*)
            DO IM = 1,NM
               READ (IOPOT,*) IDUMMY,JRNS1(IM),JRCRI(IM)
               JRNS1(IM) = JRNS1(IM) - IRSHFT
               JRCRI(IM) = JRCRI(IM) - IRSHFT
            END DO
         END IF
C
C ======================================================================
C
         CALL POSFIL(IOPOT,'TYPES     ',IPOS,1)
         READ (IOPOT,'(A)') LINE
         DO IT = 1,NT0
            READ (4,*) NDUMMY,S10DUMMY,Z_TAUX(IT)
         END DO
C
         DO IT = 1,NTLIM
            NA_TAUX(IT) = 0
         END DO
C
         DO IQ = 1,NQ
            DO IO = 1,NOQ(IQ)
               IT = IT_OQAUX(IO,IQ)
               NA_TAUX(IT) = NA_TAUX(IT) + 1
               IA = NA_TAUX(IT)
               IQ_ATAUX(IA,IT) = IQ
            END DO
         END DO
C
C=======================================================================
C
         CALL POTIO_SWSRMT(NQ,NM,IM_QAUX,SWS,JRMT,RMT,JRWS,RWS,ABAS,
     &                     QBAS,ALAT,RMESHTYPE,SCFSTATUS,NQMAX)
C
C=======================================================================
C
         RMTRED0(1:NM) = RMT(1:NM)/ALAT
C
C ======================================================================
C
      END IF
C
      WRITE (6,99015) SYSTEM_DIMENSION,SYSTEM_TYPE,SCFSTATUS
C
C-----------------------------------------------------------------------
C  keep Euler angles for magnetisation in case of new setting in input
C-----------------------------------------------------------------------
C
      QMPHI_POTFIL(1:NQ) = QMPHI(1:NQ)
      QMTET_POTFIL(1:NQ) = QMTET(1:NQ)
C
      RETURN
C ======================================================================
 100  CONTINUE
      STOP 'in <POTRD2> --> EOF reached reading types'
C
C ======================================================================
99001 FORMAT (10X,'get variables that specify the site occupation'/,10X,
     &        'from potential file ',A,/,10X,
     &        'to determine the symmetry of the system '/)
99002 FORMAT (I10,3(F16.10))
99003 FORMAT (I10,3I10,10(I6,F6.3))
99004 FORMAT (I10,2(F16.10))
99005 FORMAT (10X,A)
99006 FORMAT (/,10X,'for system:   ',A)
99007 FORMAT (10X,7I10,:,/,(10X,7I10))
99008 FORMAT (A4,6X,20I10)
99009 FORMAT (10X,5F20.10)
99010 FORMAT (10X,L3)
99011 FORMAT (I4,4X,2(F14.10,1X),F14.10,8X,I3,10X,20(I3,F6.3))
99012 FORMAT (I4,9X,F15.10,9X,F15.10)
99013 FORMAT (5E16.9)
99014 FORMAT (12X,A,3(F14.10,1X))
99015 FORMAT (/,10X,'system dimension: ',A,/,10X,'system type:      ',A,
     &        /,10X,'SCF status:       ',A,/)
      END
C*==potrdspr.f    processed by SPAG 6.70Rc at 08:47 on  8 Mar 2017
      SUBROUTINE POTRDSPR(POTFMTINP)
C   ********************************************************************
C   *                                                                  *
C   *    reading the potential file   STEP 3                           *
C   *                                                                  *
C   *  - only for standard SPRKKR format                               *
C   *  - ALL  array sizes have been set now properly                   *
C   *  - reread the  WHOLE  potfile skipping those variables that      *
C   *    might have been overwritten via the standard input file       *
C   *                                                                  *
C   *  don't overwrite the following variables:                        *
C   *                                                                  *
C   *  IREL, FULLPOT, BREITINT, ORBPOL, EXTFIELD                       *
C   *  KMROT, QMVEC, QMTET, QMPHI                                      *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:FINITE_NUCLEUS,NRSFMAX,NLMSFMAX,NPANMAX,NMMAX,
     &    NRMAX,SPHERCELL,FULLPOT,JRNSMIN,NRSFTOT,LMISF,NSF,KLMSF,ISFLM,
     &    FLMSF,NPAN,JRNS1,NRNS,JRCUT,JRCRI,JRMT,JRWS,RMT,RWS,R,BRMSH,
     &    ARMSH,DX,NM,RMESHTYPE
      USE MOD_SITES,ONLY:NQ,NQ_R,NQ_L,NQMAX,ITOQ,NOQ,IMQ,IQAT,QBAS,
     &    VLMMAD_HOST,CMNTQ
      USE MOD_ENERGY,ONLY:EFERMI
      USE MOD_LATTICE,ONLY:TXTBRAVAIS,ABAS_R,ABAS_L,ABAS,SYSTEM_TYPE,
     &    SYSTEM_DIMENSION,SWS,ALAT,BRAVAIS
      USE MOD_CONSTANTS,ONLY:PI
      USE MOD_TB,ONLY:NREF,RMTREF,VREF,IREFQ
      USE MOD_TYPES,ONLY:VMTZ,QEL,AOPT,NLSHELLMAX,NLMFPMAX,NTMAX,NLMFPT,
     &    KLMFP,BNST,VNST,LMIFP,NFPT,SEMCORSHLT,NSEMCORSHLT,NAT,IMT,
     &    RHOSPNC,RHOCHRC,RHOSPN,RHOCHR,CONC,BEXT,LOPT,Z,BT,VT,NT,TXT_T,
     &    OBS_T
      USE MOD_FILES,ONLY:IPRINT,LSYSTEM,SYSTEM
      USE MOD_CALCMODE,ONLY:BLCOUPL,ORBPOL,BREITINT,SEMICORE
      USE MOD_SCF,ONLY:SCFSTATUS,SCFVXC
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='POTRDSPR')
C
C Dummy arguments
C
      INTEGER POTFMTINP
C
C Local variables
C
      REAL*8 AUX,DELR,DRSF(:,:),QMPHI(:),QMTET(:),QMVEC(3),RMSAVB,
     &       RMSAVV,SCFMIX,SCFTOL,SCLM(:),XRSF(:,:)
      LOGICAL EXTFIELD,LLOYD
      INTEGER IA_STAT,IBZINT,IFLAG,ILOOP,IM,IO,IOPOT,IPAN,IQ,IR,IREL,IT,
     &        IT1,IT2,ITRSCF,J,KMROT,LINFO,LM,LTITLE,N,NCORT(:),NETAB(2)
     &        ,NKTAB,NRPAN(:,:),NVALT(:)
      CHARACTER*80 INFO,TITLE
      CHARACTER*10 SCFALG
C
C*** End of declarations rewritten by SPAG
C
      DATA IA_STAT/0/
C
      ALLOCATABLE QMPHI,QMTET
      ALLOCATABLE NVALT,NCORT
      ALLOCATABLE DRSF,SCLM,XRSF,NRPAN
C
      FINITE_NUCLEUS = .FALSE.
      TITLE = '                    '
      LSYSTEM = 0
C
      DO IM = 1,NMMAX
         RMT(IM) = 0.0D0
         JRMT(IM) = 0
      END DO
      DO IT = 1,NTMAX
         IMT(IT) = 1
         NAT(IT) = 1
         IQAT(1,IT) = IT
         DO N = 1,NRMAX
            RHOCHR(N,IT) = 0.0D0
            RHOSPN(N,IT) = 0.0D0
            RHOCHRC(N,IT) = 0.0D0
            RHOSPNC(N,IT) = 0.0D0
         END DO
      END DO
      IF ( ORBPOL(1:6).EQ.'BROOKS' .OR. BLCOUPL )
     &     CALL RINIT(NRMAX*2*NTMAX,AOPT)
      DO IT = 1,NTMAX
         TXT_T(IT) = '    '
      END DO
C
C
      WRITE (6,99002) POTFMTINP
C
      IOPOT = 4
C
C ======================================================================
C                       READ POTENTIAL FILE
C ======================================================================
C
      IF ( POTFMTINP.LE.5 ) THEN
C
         CALL POTRDFMT5(IOPOT,EFERMI,BRAVAIS,ALAT,ABAS,SYSTEM,LSYSTEM,
     &                  TXT_T,NQ,NT,VT,BT,AOPT,Z,BLCOUPL,BEXT,CONC,
     &                  RHOCHR,RHOSPN,RMESHTYPE,NM,DX,R,QBAS,RWS,RMT,
     &                  JRWS,JRMT,JRCRI,IMT,NAT,IQAT,NOQ,ITOQ,SEMICORE,
     &                  NSEMCORSHLT,SEMCORSHLT,SCFSTATUS,JRCUT,NRNS,
     &                  JRNS1,NPAN,NLMFPT,NFPT,LMIFP,VNST,BNST,JRNSMIN,
     &                  FULLPOT,KLMFP,NQMAX,NTMAX,NRMAX,NMMAX,NLMFPMAX,
     &                  NPANMAX,NLSHELLMAX)
C
      ELSE
C
         ALLOCATE (QMPHI(NQMAX),QMTET(NQMAX))
         ALLOCATE (NVALT(NTMAX),NCORT(NTMAX))
         NREF = NM
         IF ( IA_STAT.NE.0 ) STOP 'alloc: POTRDSPR -> IREFQ'
C
         CALL POTRD(POTFMTINP,IPRINT,ALAT,BEXT,BLCOUPL,BRAVAIS,BREITINT,
     &              LLOYD,BT,CONC,DX,EFERMI,EXTFIELD,IBZINT,IMT,INFO,
     &              IREL,ITOQ,ITRSCF,JRMT,JRWS,KMROT,LINFO,LSYSTEM,
     &              LTITLE,NETAB,NKTAB,NM,NOQ,NQ,NSEMCORSHLT,NT,ORBPOL,
     &              LOPT,QMPHI,QMTET,QMVEC,QBAS,R,RHOCHR,RHOSPN,RMSAVB,
     &              RMSAVV,RMT,RWS,SCFALG,SCFMIX,SCFTOL,SCFVXC,SEMICORE,
     &              NCORT,NVALT,SYSTEM,TITLE,TXT_T,VLMMAD_HOST,CMNTQ,VT,
     &              Z,QEL,OBS_T,VMTZ,FULLPOT,SPHERCELL,JRNS1,JRCRI,
     &              JRCUT,ARMSH,BRMSH,RMESHTYPE,IOPOT,JRNSMIN,NPAN,
     &              IREFQ,VREF,RMTREF,NREF,NLMFPT,NFPT,LMIFP,VNST,BNST,
     &              KLMFP,NMMAX,NQMAX,NRMAX,NTMAX,NLMFPMAX,NPANMAX,
     &              SYSTEM_DIMENSION,SYSTEM_TYPE,SCFSTATUS,ABAS,ABAS_L,
     &              ABAS_R,NQ_L,NQ_R)
C
      END IF
C
C ======================================================================
C                       FULL POTENTIAL
C ======================================================================
      IF ( FULLPOT ) THEN
C
         ALLOCATE (DRSF(NRSFMAX,NMMAX),SCLM(NMMAX))
         ALLOCATE (XRSF(NRSFMAX,NMMAX),NRPAN(NPANMAX,NMMAX))
C
C-------------------------------------------------- read SHAPE FUNCTIONS
C
         IF ( SPHERCELL ) THEN
C
            DO IM = 1,NM
               NSF(IM) = 1
               LMISF(1,IM) = 1
               NRPAN(1,IM) = JRCUT(1,IM)
               ISFLM(1,IM) = 1
               DO LM = 1,NLMSFMAX
                  KLMSF(LM,IM) = 0
               END DO
               KLMSF(1,IM) = 1
               DO IPAN = 2,NPAN(IM)
                  NRPAN(IPAN,IM) = JRCUT(IPAN,IM) - JRCUT(IPAN-1,IM)
               END DO
C
               SCLM(IM) = 1D0
               NRSFTOT(IM) = JRCUT(NPAN(IM),IM) - JRCUT(1,IM)
               DELR = (RWS(IM)-RMT(IM))/DBLE(NRSFTOT(IM)-1)
               DO J = 1,NRSFTOT(IM)
                  XRSF(J,IM) = (RMT(IM)+DELR*(J-1))/ALAT
                  DRSF(J,IM) = DELR/ALAT
                  FLMSF(J,1,IM) = SQRT(4.0D0*PI)
               END DO
            END DO
C
         ELSE
C
            CALL SFNREAD(NM,NRPAN,XRSF,DRSF,SCLM)
C
         END IF
C
      END IF
C
C ======================================================================
C                potential initialized -- do final checks
C ======================================================================
C
      DO IT = 1,NT
         NAT(IT) = 0
      END DO
      DO IQ = 1,NQ
         DO IO = 1,NOQ(IQ)
            IT = ITOQ(IO,IQ)
            IMQ(IQ) = IMT(IT)
            NAT(IT) = NAT(IT) + 1
            IQAT(NAT(IT),IT) = IQ
         END DO
      END DO
C
      AUX = 0.0D0
      DO IT = 1,NT
         AUX = AUX + CONC(IT)*NAT(IT)
      END DO
      IF ( ABS(AUX-NQ).GT.1D-8 )
     &      CALL STOP_MESSAGE(ROUTINE,'sum(IT) CONC(IT) * NAT(IT) <> NQ'
     &     )
C
C=======================================================================
C
      CALL POTIO_SWSRMT(NQ,NM,IMQ,SWS,JRMT,RMT,JRWS,RWS,ABAS,QBAS,ALAT,
     &                  RMESHTYPE,SCFSTATUS,NQMAX)
C
C=======================================================================
C
C---------------------------------------------------- SCFSTART completed
C
C ======================================================================
C                        set up radial mesh
C ======================================================================
C
      IF ( FULLPOT ) THEN
C
         CALL RMESHFP(NRPAN,XRSF,DRSF,SCLM)
C
      ELSE
C
         CALL RMESHASA
C
      END IF
C
C ======================================================================
C                       add nuclear potential
C ======================================================================
      IF ( RMESHTYPE.EQ.'JUELICH     ' ) THEN
         DO IT = 1,NT
            IM = IMT(IT)
            DO IR = 1,JRWS(IM)
               VT(IR,IT) = VT(IR,IT) - 2D0*Z(IT)/R(IR,IM)
            END DO
         END DO
      END IF
C
C --------------------------- check consistency of muffin-tin radius RMT
      IFLAG = 0
      DO IM = 1,NM
         IF ( (R(JRMT(IM),IM)-RMT(IM)).GT.1D-6 .OR. 
     &        (RMT(IM)-R(JRMT(IM)+1,IM)).GT.1D-6 ) THEN
            IFLAG = 1
            WRITE (*,99001)
            WRITE (*,99003) IM,R(JRMT(IM),IM),RMT(IM),R(JRMT(IM)+1,IM)
         END IF
      END DO
      IF ( IFLAG.NE.0 ) CALL STOP_MESSAGE(ROUTINE,
     &     'R(JRMT) <= RMT <= R(JRMT+1)   not fulfilled')
C
      IF ( FULLPOT ) THEN
         WRITE (6,99004) TXTBRAVAIS(BRAVAIS)(1:25),ALAT,NQ
      ELSE
         WRITE (6,99004) TXTBRAVAIS(BRAVAIS)(1:25),ALAT,NQ,SWS
      END IF
C
C
      IF ( SCFSTATUS(1:5).NE.'START' ) THEN
         IT1 = 1 - 3
         DO ILOOP = 1,(NT-1)/3 + 1
            IT1 = IT1 + 3
            IT2 = MIN(IT1+2,NT)
            WRITE (6,99006) IT1,IT2
            DO IR = 1,JRWS(1),30
               WRITE (6,'(I5,3F12.6,4(2X,2F12.6))') IR,R(IR,IMT(1)),
     &                (VT(IR,IT)*R(IR,IMT(IT)),BT(IR,IT)*R(IR,IMT(IT)),
     &                IT=IT1,IT2)
            END DO
            IR = JRWS(1)
            WRITE (6,'(I5,3F12.6,4(2X,2F12.6))') IR,R(IR,IMT(1)),
     &             (VT(IR,IT)*R(IR,IMT(IT)),BT(IR,IT)*R(IR,IMT(IT)),
     &             IT=IT1,IT2)
         END DO
C
         WRITE (6,*) ' '
C
Cc         IF ( DABS(Z(1)+VT(1,1)*R(1,IMT(1))/2.0D0).GT.0.5D0 )
Cc     &        FINITE_NUCLEUS = .TRUE.
      END IF
C
      IF ( FINITE_NUCLEUS ) WRITE (6,99005)
C
      IF ( ALLOCATED(QMPHI) ) THEN
         DEALLOCATE (QMPHI,QMTET)
         DEALLOCATE (NVALT,NCORT)
      END IF
      IF ( ALLOCATED(DRSF) ) DEALLOCATE (DRSF,SCLM,XRSF,NRPAN)
C
99001 FORMAT (/,' ##### TROUBLE in <POTRDSPR> ',51('#'))
99002 FORMAT (//,1X,79('*'),/,33X,'<POTRDSPR>',/,1X,79('*'),//,10X,
     &        'POTFMTINP =',I3,/)
99003 FORMAT (10X,' MESH        ',I15,/,10X,' R <= R(mt)  ',F15.10,/,
     &        10X,' R(mt)       ',F15.10,/,10X,' R >= R(mt)  ',F15.10)
99004 FORMAT (/,10X,'Bravais lattice:                ',A,/,10X,
     &        'lattice constant  ALAT          ',F12.5,/,10X,
     &        'number of sites   NQ            ',I12,/,:,10X,
     &        'average Wigner-Seitz radius     ',F12.5,//)
99005 FORMAT (//,10X,'>>> FINITE NUCLEUS CALCULATION <<<',//)
99006 FORMAT (/,10X,'potential data for IT =',I3,' .. ',I3,/,10X,
     &        '=================================',/)
      END
C*==potrdfmt5.f    processed by SPAG 6.70Rc at 08:47 on  8 Mar 2017
      SUBROUTINE POTRDFMT5(IOPOT,EFERMI,BRAVAIS,ALAT,ABAS,SYSTEM,
     &                     LSYSTEM,TXT_T,NQ,NT,VT,BT,AOPT,Z,BLCOUPL,
     &                     BEXT,CONC,RHOCHR,RHOSPN,RMESHTYPE,NM,DX,R,
     &                     QBAS,RWS,RMT,JRWS,JRMT,JRCRI,IMT,NAT,IQAT,
     &                     NOQ,ITOQ,SEMICORE,NSEMCORSHLT,SEMCORSHLT,
     &                     SCFSTATUS,JRCUT,NRNS,JRNS1,NPAN,NLMFPT,NFPT,
     &                     LMIFP,VNST,BNST,JRNSMIN,FULLPOT,KLMFP,NQMAX,
     &                     NTMAX,NRMAX,NMMAX,NLMFPMAX,NPANMAX,
     &                     NLSHELLMAX)
C   ********************************************************************
C   *                                                                  *
C   *    reading the potential file   STEP 3                           *
C   *                                                                  *
C   *  - only for standard SPRKKR format up to FORMAT VERSION 5        *
C   *  - ALL  array sizes have been set now properly                   *
C   *  - reread the  WHOLE  potfile skipping those variables that      *
C   *    might have been overwritten via the standard input file       *
C   *                                                                  *
C   *  don't overwrite the following variables:                        *
C   *                                                                  *
C   *  IREL, FULLPOT, BREITINT, ORBPOL, EXTFIELD                       *
C   *  KMROT, QMVEC, QMTET, QMPHI                                      *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IDUMMY,RDUMMY
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='POTRDFMT5')
C
C Dummy arguments
C
      REAL*8 ALAT,BEXT,EFERMI
      LOGICAL BLCOUPL,FULLPOT,SEMICORE
      INTEGER BRAVAIS,IOPOT,JRNSMIN,LSYSTEM,NLMFPMAX,NLSHELLMAX,NM,
     &        NMMAX,NPANMAX,NQ,NQMAX,NRMAX,NT,NTMAX
      CHARACTER*12 RMESHTYPE
      CHARACTER*10 SCFSTATUS
      CHARACTER*80 SYSTEM
      REAL*8 ABAS(3,3),AOPT(NRMAX,2,NTMAX),
     &       BNST(JRNSMIN:NRMAX,NLMFPMAX,NTMAX),BT(NRMAX,NTMAX),
     &       CONC(NTMAX),DX(NMMAX),QBAS(3,NQMAX),R(NRMAX,NMMAX),
     &       RHOCHR(NRMAX,NTMAX),RHOSPN(NRMAX,NTMAX),RMT(NMMAX),
     &       RWS(NMMAX),VNST(JRNSMIN:NRMAX,NLMFPMAX,NTMAX),
     &       VT(NRMAX,NTMAX)
      INTEGER IMT(NTMAX),IQAT(NQMAX,NTMAX),ITOQ(NTMAX,NQMAX),
     &        JRCRI(NMMAX),JRCUT(0:NPANMAX,NMMAX),JRMT(NMMAX),
     &        JRNS1(NMMAX),JRWS(NMMAX),KLMFP(NLMFPMAX,NTMAX),
     &        LMIFP(NLMFPMAX,NTMAX),NAT(NTMAX),NFPT(NTMAX),NLMFPT(NTMAX)
     &        ,NOQ(NQMAX),NPAN(NMMAX),NRNS(NMMAX),NSEMCORSHLT(NTMAX),
     &        Z(NTMAX)
      CHARACTER*2 SEMCORSHLT(NLSHELLMAX,NTMAX)
      CHARACTER*8 TXT_T(NTMAX)
C
C Local variables
C
      INTEGER IA,IC,IFP,IM,IO,IPAN,IPOS,IQ,IR,IRTOP,ISHL,IT,ITP,J,JQ,LM,
     &        MS,POTFMTINP
      CHARACTER*80 LINE,TITLE
      CHARACTER*10 ORBPOLINP
      CHARACTER*4 STR4
C
C*** End of declarations rewritten by SPAG
C
      WRITE (6,99003)
C
C **********************************************************************
C                   read in the standard potential set
C **********************************************************************
C
      REWIND IOPOT
      READ (IOPOT,'(A)') LINE
      IF ( LINE(1:6).NE.'SPRKKR' )
     &      CALL STOP_MESSAGE(ROUTINE,'LINE(1:6).NE.SPRKKR')
C
      FULLPOT = (LINE(11:24).EQ.'full potential')
      IF ( LINE(11:19).EQ.'SCF-start' ) THEN
         SCFSTATUS = 'START'
      ELSE
         SCFSTATUS = 'ITERATING'
      END IF
      WRITE (6,99004) LINE(1:LEN_TRIM(LINE))
      READ (IOPOT,99004) TITLE
      IF ( TITLE(1:7).NE.'VERSION' ) THEN
         POTFMTINP = 0
      ELSE
         POTFMTINP = ICHAR(TITLE(10:10)) - ICHAR('1') + 1
         READ (IOPOT,99004) TITLE
      END IF
      IF ( POTFMTINP.LT.3 ) CALL STOP_MESSAGE(ROUTINE,'POTFMTINP < 3')
      READ (IOPOT,99004) SYSTEM
      LSYSTEM = LEN_TRIM(SYSTEM)
      WRITE (6,99005) SYSTEM(1:LSYSTEM)
      READ (IOPOT,99006) NM
      READ (IOPOT,99006) NQ
      READ (IOPOT,99006) NT
      IF ( NM.GT.NMMAX ) CALL STOP_MESSAGE(ROUTINE,'NM > NMMAX')
      IF ( NQ.GT.NQMAX ) CALL STOP_MESSAGE(ROUTINE,'NQ > NQMAX')
      IF ( NT.GT.NTMAX ) CALL STOP_MESSAGE(ROUTINE,'NT > NTMAX')
C
CcccccREAD (iopot,99018) NS
      READ (IOPOT,*)
CcccccREAD (iopot,99018) NETAB
      READ (IOPOT,*)
C>>>>>READ (iopot,99017) IREL
      READ (IOPOT,*)
CcccccREAD (iopot,99018) IBZINT
      READ (IOPOT,*)
CcccccREAD (iopot,99018) NPTBZ
      READ (IOPOT,*)
      READ (IOPOT,*)
CcccccREAD (iopot,99015) SCFVXC
      READ (IOPOT,*)
CcccccREAD (iopot,99015) SCFALG
      READ (IOPOT,*)
CcccccREAD (iopot,99018) ITRSCF
      READ (IOPOT,*)
CcccccREAD (iopot,99020) SCFMIX
      READ (IOPOT,*)
CcccccREAD (iopot,99020) SCFTOL
      READ (IOPOT,*)
CcccccREAD (iopot,99021) SEMICORE
      IF ( POTFMTINP.GE.4 ) THEN
         READ (IOPOT,99001) SEMICORE
      ELSE
         SEMICORE = .FALSE.
      END IF
C>>>>>READ (iopot,99020) BREITINT
      READ (IOPOT,*)
      READ (IOPOT,99004) ORBPOLINP
      IF ( ORBPOLINP.EQ.'          ' ) ORBPOLINP = 'NONE'
C     READ (iopot,99020) EXTFIELD
      READ (IOPOT,*)
      READ (IOPOT,99009) BLCOUPL
      READ (IOPOT,99008) BEXT
C     READ (iopot,99017) KMROT
      READ (IOPOT,99006) IDUMMY
C     READ (iopot,99019) (QMVEC(I),I=1,3)
      READ (IOPOT,*)
      READ (IOPOT,99008) EFERMI
      READ (IOPOT,99007) STR4,BRAVAIS
      IF ( STR4.EQ.'LAT ' ) CALL STOP_MESSAGE(ROUTINE,'found LAT')
      READ (IOPOT,99006) IDUMMY
      READ (IOPOT,99008) ALAT
      READ (IOPOT,*)
      READ (IOPOT,*)
      DO J = 1,3
         READ (IOPOT,99008) (ABAS(IC,J),IC=1,3)
      END DO
C --------------------------------------- read in radial mesh parameters
      READ (IOPOT,*)
      READ (IOPOT,*)
      READ (IOPOT,99004) RMESHTYPE
      DO IM = 1,NM
         READ (IOPOT,99006) IDUMMY
         READ (IOPOT,99008) RWS(IM)
         READ (IOPOT,99008) RMT(IM)
         READ (IOPOT,99006) JRWS(IM)
         READ (IOPOT,99006) JRMT(IM)
         IF ( RMESHTYPE.EQ.'EXPONENTIAL ' ) THEN
            READ (IOPOT,99008) R(1,IM)
            READ (IOPOT,99008) DX(IM)
            IF ( FULLPOT ) THEN
               READ (IOPOT,99006) JRNS1(IM)
               READ (IOPOT,99006) JRCRI(IM)
               READ (IOPOT,99006) NPAN(IM)
               READ (IOPOT,99006) (JRCUT(IPAN,IM),IPAN=1,NPAN(IM))
               JRCUT(0,IM) = 0
               READ (IOPOT,99004) STR4
               NRNS(IM) = JRCUT(NPAN(IM),IM) - JRNS1(IM)
            END IF
         ELSE
            WRITE (6,*) 'RMESHTYPE = ',RMESHTYPE
            CALL STOP_MESSAGE(ROUTINE,'check RMESHTYPE in POTFILE')
         END IF
      END DO
C
      READ (IOPOT,'(/)')
      DO IQ = 1,NQ
         READ (IOPOT,99012) IDUMMY,(QBAS(IC,IQ),IC=1,3),NOQ(IQ),
     &                      (ITOQ(IO,IQ),CONC(ITOQ(IO,IQ)),IO=1,NOQ(IQ))
      END DO
      DO IQ = 1,NQ
         READ (IOPOT,99013) IDUMMY
      END DO
C
      CALL EXTEND_TXT_T
C
C ------------------------------------------------------ MADELUNG MATRIX
      READ (IOPOT,'(/)')
      DO IQ = 1,NQ
         READ (IOPOT,99010) (RDUMMY,JQ=1,NQ)
      END DO
C
C -------------------------------------- read in atomic type information
      DO IT = 1,NT
         IM = IMT(IT)
         IF ( FULLPOT ) THEN
            IRTOP = JRCRI(IM)
         ELSE
            IRTOP = JRWS(IM)
         END IF
         READ (IOPOT,'(/,10X,I10,/,A)') ITP,TXT_T(IT)
         READ (IOPOT,99006) Z(IT)
         READ (IOPOT,99006) NAT(IT)
         READ (IOPOT,99008) CONC(IT)
         READ (IOPOT,99006) IMT(IT)
         READ (IOPOT,99006) (IQAT(IA,IT),IA=1,NAT(IT))
         IF ( POTFMTINP.GE.4 ) READ (IOPOT,99002) NSEMCORSHLT(IT),
     &                               (SEMCORSHLT(ISHL,IT),ISHL=1,
     &                               NSEMCORSHLT(IT))
         IM = IMT(IT)
         IF ( SCFSTATUS(1:5).EQ.'START' ) THEN
            DO IR = 1,JRWS(IM)
               VT(IR,IT) = 0D0
               BT(IR,IT) = 0D0
            END DO
         ELSE
C
            READ (IOPOT,99010) (VT(IR,IT),IR=1,IRTOP)
            READ (IOPOT,99010) (BT(IR,IT),IR=1,IRTOP)
C
            IF ( FULLPOT ) THEN
               NLMFPT(IT) = 1
               DO LM = 1,NLMFPMAX
                  KLMFP(LM,IT) = 0
               END DO
               KLMFP(1,IT) = 1
               LMIFP(1,IT) = 1
               READ (IOPOT,99006) NFPT(IT)
               DO IFP = 2,NFPT(IT)
                  READ (IOPOT,99006) LM
                  KLMFP(LM,IT) = 1
                  LMIFP(IFP,IT) = LM
                  NLMFPT(IT) = MAX(LM,NLMFPT(IT))
                  READ (IOPOT,99010) (VNST(IR,LM,IT),IR=JRNS1(IM),IRTOP)
                  READ (IOPOT,99010) (BNST(IR,LM,IT),IR=JRNS1(IM),IRTOP)
               END DO
            END IF
C
         END IF
      END DO
C
      IF ( ORBPOLINP(1:4).NE.'NONE' .AND. ORBPOLINP(1:4).NE.'DMFT' )
     &     CALL POSFIL(IOPOT,'V-ORB     ',IPOS,1)
      IF ( ORBPOLINP(1:6).EQ.'BROOKS' ) THEN
         DO IT = 1,NT
            IM = IMT(IT)
            READ (IOPOT,'(/,10X,I10,/,A)') ITP,TXT_T(IT)
            DO MS = 1,2
               READ (IOPOT,99010) (AOPT(IR,MS,IT),IR=1,JRWS(IM))
            END DO
            WRITE (6,*) '  vector potential read for IT=',ITP,'   ',
     &                  TXT_T(IT)
         END DO
      END IF
C
      CALL POSFIL(IOPOT,'RHO    SCF',IPOS,0)
C
      IF ( IPOS.EQ.1 ) THEN
         DO IT = 1,NT
            IM = IMT(IT)
            READ (IOPOT,'(/,10X,I10,/,A)') ITP,TXT_T(IT)
            READ (IOPOT,99010) (RHOCHR(IR,IT),IR=1,JRWS(IM))
            READ (IOPOT,99010) (RHOSPN(IR,IT),IR=1,JRWS(IM))
            WRITE (6,99011) ITP,TXT_T(IT)
         END DO
      END IF
C
      RETURN
C ======================================================================
99001 FORMAT (10X,L3)
99002 FORMAT (10X,I5,10(2X,A2))
99003 FORMAT (//,1X,79('*'),/,33X,'<POTRDFMT5>',/,1X,79('*'),/)
99004 FORMAT (10X,A)
99005 FORMAT (/,10X,'for system:   ',A)
99006 FORMAT (10X,7I10,:,/,(10X,7I10))
99007 FORMAT (A4,6X,20I10)
99008 FORMAT (10X,5F20.10)
99009 FORMAT (10X,L3)
99010 FORMAT (5E16.9)
99011 FORMAT (10X,'charge and spin density read for IT=',I2,3X,A)
99012 FORMAT (I4,4X,2(F14.10,1X),F14.10,8X,I3,10X,20(I3,F6.3))
99013 FORMAT (I4,9X,F15.10,9X,F15.10)
      END
C*==potrd.f    processed by SPAG 6.70Rc at 08:47 on  8 Mar 2017
      SUBROUTINE POTRD(POTFMTINP,IPRINT,ALAT,BEXT,BLCOUPL,BRAVAIS,
     &                 BREITINT,LLOYD,BT,CONC,DX,EFERMI,EXTFIELD,IBZINT,
     &                 IMT,INFO,IREL,ITOQ,ITRSCF,JRMT,JRWS,KMROT,LINFO,
     &                 LSYSTEM,LTITLE,NETAB,NKTAB,NM,NOQ,NQ,NSEMCORSHLT,
     &                 NT,ORBPOL,LOPT,QMPHI,QMTET,QMVEC,QBAS,R,RHOCHR,
     &                 RHOSPN,RMSAVB,RMSAVV,RMT,RWS,SCFALG,SCFMIX,
     &                 SCFTOL,SCFVXC,SEMICORE,NCORT,NVALT,SYSTEM,TITLE,
     &                 TXT_T,VLMMAD_HOST,CMNTQ,VT,Z,QEL,OBS_T,VMTZ,
     &                 FULLPOT,SPHERCELL,JRNS1,JRCRI,JRCUT,ARMSH,BRMSH,
     &                 RMESHTYPE,IOPOT,JRNSMIN,NPAN,IREFQ,VREF,RMTREF,
     &                 NREF,NLMFPT,NFPT,LMIFP,VNST,BNST,KLMFP,NMMAX,
     &                 NQMAX,NRMAX,NTMAX,NLMFPMAX,NPANMAX,
     &                 SYSTEM_DIMENSION,SYSTEM_TYPE,SCFSTATUS,ABAS,
     &                 ABAS_L,ABAS_R,NQ_L,NQ_R)
C   ********************************************************************
C   *                                                                  *
C   *   read the SCF-potential for  POTFMT > 5                         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:BASSCALE,CARTESIAN
      USE MOD_SCF,ONLY:SCF_THETA_DEPENDENT_POT
      USE MOD_FILES,ONLY:POTFIL,LPOTFIL,IDUMMY,RDUMMY,PACKAGE
      USE MOD_TYPES,ONLY:ABIT,AOPT,RHO2NS,MTET_T,MPHI_T,MGAM_T,X_CHEM_T
      USE MOD_ANGMOM,ONLY:NOBSMAX,IDOS,ISMT,IOMT,IHFF,NLABIMAX
      USE MOD_ENERGY,ONLY:EWORK
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='POTRD')
C
C Dummy arguments
C
      REAL*8 ALAT,BEXT,EFERMI,RMSAVB,RMSAVV,SCFMIX,SCFTOL,VMTZ
      LOGICAL BLCOUPL,BREITINT,EXTFIELD,FULLPOT,LLOYD,SEMICORE,SPHERCELL
      INTEGER BRAVAIS,IBZINT,IOPOT,IPRINT,IREL,ITRSCF,JRNSMIN,KMROT,
     &        LINFO,LSYSTEM,LTITLE,NKTAB,NLMFPMAX,NM,NMMAX,NPANMAX,NQ,
     &        NQMAX,NQ_L,NQ_R,NREF,NRMAX,NT,NTMAX,POTFMTINP
      CHARACTER*80 INFO,SYSTEM,TITLE
      CHARACTER*10 ORBPOL,SCFALG,SCFSTATUS,SCFVXC,SYSTEM_DIMENSION,
     &             SYSTEM_TYPE
      CHARACTER*12 RMESHTYPE
      REAL*8 ABAS(3,3),ABAS_L(3,3),ABAS_R(3,3),ARMSH(NMMAX),
     &       BNST(JRNSMIN:NRMAX,NLMFPMAX,NTMAX),BRMSH(NMMAX),
     &       BT(NRMAX,NTMAX),CMNTQ(NLMFPMAX,NQMAX),CONC(NTMAX),DX(NMMAX)
     &       ,OBS_T(0:3,NOBSMAX,NTMAX),QBAS(3,NQMAX),QEL(NTMAX),
     &       QMPHI(NQMAX),QMTET(NQMAX),QMVEC(3),R(NRMAX,NMMAX),
     &       RHOCHR(NRMAX,NTMAX),RHOSPN(NRMAX,NTMAX),RMT(NMMAX),
     &       RMTREF(NREF),RWS(NMMAX),VLMMAD_HOST(NLMFPMAX,NQMAX),
     &       VNST(JRNSMIN:NRMAX,NLMFPMAX,NTMAX),VREF(NREF),
     &       VT(NRMAX,NTMAX)
      INTEGER IMT(NTMAX),IREFQ(NQMAX),ITOQ(NTMAX,NQMAX),JRCRI(NMMAX),
     &        JRCUT(0:NPANMAX,NMMAX),JRMT(NMMAX),JRNS1(NMMAX),
     &        JRWS(NMMAX),KLMFP(NLMFPMAX,NTMAX),LMIFP(NLMFPMAX,NTMAX),
     &        LOPT(NTMAX),NCORT(NTMAX),NETAB(2),NFPT(NTMAX),
     &        NLMFPT(NTMAX),NOQ(NQMAX),NPAN(NMMAX),NSEMCORSHLT(NTMAX),
     &        NVALT(NTMAX),Z(NTMAX)
      CHARACTER*8 TXT_T(NTMAX)
C
C Local variables
C
      LOGICAL CHECK
      CHARACTER*40 FMT02,FMT03,FMT04,FMT05,FMT06,FMT07
      CHARACTER*80 HEADER,LINE
      INTEGER IC,IFLAG,IFP,ILA,IM,IO,IPAN,IPOS,IQ,IR,IRSHFT,IRTOP,IT,
     &        IWRI,LM,MS,NLMTOP,NSPIN
      CHARACTER*10 PACKAGEIN,STR10
      REAL*8 QMGAM(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE QMGAM
C
      ALLOCATE (QMGAM(NQMAX))
      QMGAM(:) = 0D0
C
      IF ( IPRINT.GE.3 ) THEN
         CHECK = .TRUE.
         IWRI = 6
      ELSE
         CHECK = .FALSE.
         IWRI = 0
      END IF
      LINE = '          '
C
C-----------------------------------------------------------------------
C                    set appropriate input format
C-----------------------------------------------------------------------
C
      IF ( POTFMTINP.EQ.6 ) THEN
C >> CONC
         FMT02 = '(I10,3I10,10(I6,1X,F5.3))'
C >> QBAS,VREF,RMTREF,VLMMAD_HOST,QMTET,QMPHI
         FMT03 = '(I10,3(F16.10))'
C >> VLMMAD_HOST
         FMT04 = '(6X,2I4,3(F16.10))'
C >> R1,DX,RMT,RWS
         FMT05 = '(1X,I4,2(F16.10),2(I5,F16.10))'
C >> NVALT
         FMT06 = '(1X,I4,5X,A,I10,I10,I10,9X,I4)'
C >> V,B,...
         FMT07 = '(1P,5E16.9)'
      ELSE
         FMT02 = '(I10,3I10,10(I6,F8.5))'
         FMT03 = '(I10,1P,3E22.14)'
         FMT04 = '(6X,2I4,1P,3E22.14)'
         FMT05 = '(1X,I4,1P,2E22.14,2(I5,1P,E22.14))'
         FMT06 = '(1X,I4,5X,A,I10,I10,I10,9X,I4)'
         FMT07 = '(1P,5E22.14)'
      END IF
C
C-----------------------------------------------------------------------
C
      REWIND IOPOT
      CALL READKWSTR(IOPOT,'HEADER    ',HEADER,LINE,80,IWRI,1)
      WRITE (6,*) HEADER
      LINE = '          '
      CALL READKWSTR(IOPOT,'TITLE     ',TITLE,LINE,80,IWRI,1)
      LTITLE = LEN_TRIM(TITLE)
      CALL READKWSTR(IOPOT,'SYSTEM    ',SYSTEM,LINE,80,IWRI,1)
      LSYSTEM = LEN_TRIM(SYSTEM)
      CALL READKWSTR(IOPOT,'PACKAGE   ',PACKAGEIN,'          ',10,IWRI,
     &               1)
C
      CALL READKWSTR(IOPOT,'INFO      ',INFO,LINE,80,IWRI,0)
      LINFO = LEN_TRIM(INFO)
C
C ======================================================================
C
      REWIND IOPOT
      CALL READKWINT(IOPOT,'NQ        ',NQ,0,IWRI,1)
      CALL READKWINT(IOPOT,'NT        ',NT,0,IWRI,1)
      CALL READKWINT(IOPOT,'NM        ',NM,0,IWRI,1)
      CALL READKWINT(IOPOT,'IREL      ',IREL,3,IWRI,0)
      CALL READKWINT(IOPOT,'NSPIN     ',NSPIN,1,IWRI,0)
C
      IF ( NM.GT.NMMAX ) CALL STOP_MESSAGE(ROUTINE,'NM > NMMAX')
      IF ( NQ.GT.NQMAX ) CALL STOP_MESSAGE(ROUTINE,'NQ > NQMAX')
      IF ( NT.GT.NTMAX ) CALL STOP_MESSAGE(ROUTINE,'NT > NTMAX')
C
C ======================================================================
C
      CALL READKWLOG(IOPOT,'FULLPOT   ',FULLPOT,.FALSE.,IWRI,0)
      CALL READKWLOG(IOPOT,'SPHERCELL ',SPHERCELL,.FALSE.,IWRI,0)
      CALL READKWLOG(IOPOT,'BREITINT  ',BREITINT,.FALSE.,IWRI,0)
      CALL READKWSTR(IOPOT,'ORBPOL    ',ORBPOL,'NONE      ',10,IWRI,0)
      CALL READKWLOG(IOPOT,'EXTFIELD  ',EXTFIELD,.FALSE.,IWRI,0)
      CALL READKWLOG(IOPOT,'BLCOUPL   ',BLCOUPL,.FALSE.,IWRI,0)
      CALL READKWREAL(IOPOT,'BEXT      ',BEXT,0D0,IWRI,0)
      CALL READKWLOG(IOPOT,'SEMICORE  ',SEMICORE,.FALSE.,IWRI,0)
      CALL READKWLOG(IOPOT,'LLOYD     ',LLOYD,.FALSE.,IWRI,0)
      CALL READKWINT(IOPOT,'NE        ',NETAB(1),0,IWRI,0)
      CALL READKWINT(IOPOT,'IBZINT    ',IBZINT,3,IWRI,0)
      CALL READKWINT(IOPOT,'NKTAB     ',NKTAB,100,IWRI,0)
      IF ( POTFMTINP.GE.8 )
     &     CALL READKWLOG(IOPOT,'TETDEPPOT ',SCF_THETA_DEPENDENT_POT,
     &     .FALSE.,IWRI,0)
      CALL READKWSTR(IOPOT,'XC-POT    ',SCFVXC,'VWN       ',10,IWRI,0)
      CALL READKWSTR(IOPOT,'SCF-ALG   ',SCFALG,'BROYDEN2  ',10,IWRI,0)
      CALL READKWINT(IOPOT,'SCF-ITER  ',ITRSCF,100,IWRI,0)
      CALL READKWREAL(IOPOT,'SCF-MIX   ',SCFMIX,0.2D0,IWRI,0)
      CALL READKWREAL(IOPOT,'SCF-TOL   ',SCFTOL,0.00001D0,IWRI,0)
      CALL READKWREAL(IOPOT,'RMSAVV    ',RMSAVV,999999.0D0,IWRI,0)
      CALL READKWREAL(IOPOT,'RMSAVB    ',RMSAVB,999999.0D0,IWRI,0)
      CALL READKWREAL(IOPOT,'EF        ',EFERMI,0.0D0,IWRI,0)
      CALL READKWREAL(IOPOT,'VMTZ      ',VMTZ,0.0D0,IWRI,0)
      CALL READKWREAL(IOPOT,'EWORK     ',EWORK,0.0D0,IWRI,0)
C
C ======================================================================
C
      CALL READKWSTR(IOPOT,'SYSDIM    ',SYSTEM_DIMENSION,'3D        ',
     &               10,IWRI,1)
      CALL READKWSTR(IOPOT,'SYSTYPE   ',SYSTEM_TYPE,'BULK      ',10,
     &               IWRI,1)
      CALL READKWINT(IOPOT,'BRAVAIS   ',BRAVAIS,0,IWRI,0)
      CALL READKWREAL(IOPOT,'ALAT      ',ALAT,0D0,IWRI,1)
      CALL READKWRARR(IOPOT,'A(1)      ',ABAS(1,1),3,0D0,IWRI,1)
      CALL READKWRARR(IOPOT,'A(2)      ',ABAS(1,2),3,0D0,IWRI,1)
      CALL READKWRARR(IOPOT,'A(3)      ',ABAS(1,3),3,0D0,IWRI,1)
C
      IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' ) THEN
         ABAS_L(1:3,1:2) = ABAS(1:3,1:2)
         ABAS_R(1:3,1:2) = ABAS(1:3,1:2)
C
         CALL READKWINT(IOPOT,'NQ_L      ',NQ_L,0,IWRI,1)
         CALL READKWRARR(IOPOT,'A_L(3)    ',ABAS_L(1,3),3,0D0,IWRI,1)
C
         CALL READKWINT(IOPOT,'NQ_R      ',NQ_R,0,IWRI,1)
         CALL READKWRARR(IOPOT,'A_R(3)    ',ABAS_R(1,3),3,0D0,IWRI,1)
      END IF
C
C ======================================================================
      CALL READKWLOG(IOPOT,'CARTESIAN ',CARTESIAN,.TRUE.,IWRI,0)
      CALL READKWRARR(IOPOT,'BASSCALE  ',BASSCALE,3,0D0,IWRI,1)
      READ (IOPOT,'(A80)') LINE
      IF ( CHECK ) WRITE (6,'(A)') LINE(1:LEN_TRIM(LINE))
      DO IQ = 1,NQ
         READ (IOPOT,*) IDUMMY,(QBAS(IC,IQ),IC=1,3)
         IF ( CHECK ) WRITE (6,FMT=FMT03) IDUMMY,(QBAS(IC,IQ),IC=1,3)
      END DO
C ======================================================================
      CALL POSFIL(IOPOT,'OCCUPATION',IPOS,1)
      READ (IOPOT,'(A80)') LINE
      IF ( CHECK ) WRITE (6,'(A)') LINE(1:LEN_TRIM(LINE))
      DO IQ = 1,NQ
         READ (IOPOT,*) IDUMMY,IREFQ(IQ),IDUMMY,NOQ(IQ),
     &                  (ITOQ(IO,IQ),X_CHEM_T(ITOQ(IO,IQ)),IO=1,NOQ(IQ))
C
         DO IO = 1,NOQ(IQ)
            IMT(ITOQ(IO,IQ)) = IDUMMY
         END DO
C
         IF ( CHECK ) WRITE (6,FMT=FMT02) IDUMMY,IREFQ(IQ),IDUMMY,
     &                       NOQ(IQ),
     &                       (ITOQ(IO,IQ),X_CHEM_T(ITOQ(IO,IQ)),IO=1,
     &                       NOQ(IQ))
      END DO
C ======================================================================
      CALL READKWINT(IOPOT,'NREF      ',NREF,1,IWRI,1)
      READ (IOPOT,'(A80)') LINE
      IF ( CHECK ) WRITE (6,'(A)') LINE(1:LEN_TRIM(LINE))
      DO IQ = 1,NREF
         READ (IOPOT,*) IDUMMY,VREF(IQ),RMTREF(IQ)
         IF ( CHECK ) WRITE (6,FMT=FMT03) IDUMMY,VREF(IQ),RMTREF(IQ)
      END DO
C ======================================================================
      IF ( SCFSTATUS(1:5).NE.'START' ) THEN
C
         CALL POSFIL(IOPOT,'HOST MADEL',IPOS,0)
         IFLAG = 1
         IF ( IPOS.EQ.1 ) THEN
            READ (IOPOT,'(A80)') LINE
            IF ( CHECK ) WRITE (6,'(A)') LINE(1:LEN_TRIM(LINE))
            CALL READKWINT(IOPOT,'NLMTOP-POT',NLMTOP,0,IWRI,0)
C           IF ( NLMTOP.GT.NLMFPMAX ) STOP '<POTRD>: nlmtop > NLMFPMAX'
C------ allow more LMs in pot file than requested for actual calculation
            DO IQ = 1,NQ
               DO LM = 1,MIN(NLMTOP,NLMFPMAX)
                  READ (IOPOT,*,ERR=50,END=50) IDUMMY,IDUMMY,
     &                  VLMMAD_HOST(LM,IQ)
                  IF ( CHECK ) WRITE (6,FMT=FMT04) IDUMMY,IDUMMY,
     &                                VLMMAD_HOST(LM,IQ)
               END DO
               DO LM = MIN(NLMTOP,NLMFPMAX) + 1,NLMTOP
                  READ (IOPOT,*,ERR=50,END=50) IDUMMY
               END DO
            END DO
            IFLAG = 0
         END IF
 50      CONTINUE
         IF ( IPOS.NE.1 .OR. IFLAG.EQ.1 ) THEN
            VLMMAD_HOST(1:NLMFPMAX,1:NQ) = 0D0
            WRITE (6,99001) POTFIL(1:LPOTFIL)
         END IF
C ----------------------------------------------------------------------
         CALL POSFIL(IOPOT,'CHARGE MOM',IPOS,0)
         IFLAG = 1
         IF ( IPOS.EQ.1 ) THEN
            READ (IOPOT,'(A80)') LINE
            IF ( CHECK ) WRITE (6,'(A)') LINE(1:LEN_TRIM(LINE))
            CALL READKWINT(IOPOT,'NLMTOP-CHR',NLMTOP,1,IWRI,1)
C           IF ( NLMTOP.GT.NLMFPMAX ) STOP '<POTRD>: nlmtop > NLMFPMAX'
C------ allow more LMs in pot file than requested for actual calculation
            DO IQ = 1,NQ
               DO LM = 1,MIN(NLMTOP,NLMFPMAX)
                  READ (IOPOT,*,ERR=100,END=100) IDUMMY,IDUMMY,
     &                  CMNTQ(LM,IQ)
                  IF ( CHECK ) WRITE (6,FMT=FMT04) IDUMMY,IDUMMY,
     &                                CMNTQ(LM,IQ)
               END DO
               DO LM = MIN(NLMTOP,NLMFPMAX) + 1,NLMTOP
                  READ (IOPOT,*,ERR=100,END=100) IDUMMY
               END DO
            END DO
            IFLAG = 0
         END IF
 100     CONTINUE
         IF ( IPOS.NE.1 .OR. IFLAG.EQ.1 ) THEN
            CMNTQ(1:NLMFPMAX,1:NQ) = 0D0
            WRITE (6,99002) POTFIL(1:LPOTFIL)
         END IF
C
      END IF
C ======================================================================
      CALL READKWINT(IOPOT,'KMROT     ',KMROT,0,IWRI,1)
      CALL READKWRARR(IOPOT,'QMVEC     ',QMVEC,3,0D0,IWRI,1)
      READ (IOPOT,'(A80)') LINE
      IF ( CHECK ) WRITE (6,'(A)') LINE(1:LEN_TRIM(LINE))
      DO IQ = 1,NQ
         READ (IOPOT,*) IDUMMY,QMTET(IQ),QMPHI(IQ)
         IF ( CHECK ) WRITE (6,FMT=FMT03) IDUMMY,QMTET(IQ),QMPHI(IQ)
      END DO
      IF ( POTFMTINP.GE.8 ) THEN
         READ (IOPOT,'(A80)') LINE
         IF ( CHECK ) WRITE (6,'(A)') LINE(1:LEN_TRIM(LINE))
         DO IT = 1,NT
            READ (IOPOT,*) IDUMMY,MTET_T(IT),MPHI_T(IT),MGAM_T(IT)
            IF ( CHECK ) WRITE (6,FMT=FMT03) IDUMMY,MTET_T(IT),
     &                          MPHI_T(IT),MGAM_T(IT)
         END DO
      END IF
C
C--------------------------- do a final test and adjust KMROT eventually
C
      CALL MROT_CHECK_SETTINGS(NQ,IREL,QMPHI,QMTET,QMGAM,KMROT,QMVEC,
     &                         NQMAX)
C
C ======================================================================
      CALL POSFIL(IOPOT,'MESH INFOR',IPOS,1)
      CALL READKWSTR(IOPOT,'MESH-TYPE ',RMESHTYPE,'            ',12,
     &               IWRI,1)
C ----------------------------------------------------------------------
      IRSHFT = 0
      IF ( RMESHTYPE.EQ.'EXPONENTIAL ' ) THEN
         READ (IOPOT,'(A80)') LINE
         IF ( CHECK ) WRITE (6,'(A)') LINE(1:LEN_TRIM(LINE))
         DO IM = 1,NM
            READ (IOPOT,*) IDUMMY,R(1,IM),DX(IM),JRMT(IM),RMT(IM),
     &                     JRWS(IM),RWS(IM)
            IF ( CHECK ) WRITE (6,FMT=FMT05) IDUMMY,R(1,IM),DX(IM),
     &                          JRMT(IM),RMT(IM),JRWS(IM),RWS(IM)
         END DO
      ELSE IF ( RMESHTYPE.EQ.'JUELICH     ' ) THEN
C ----------------------------------------------------------------------
         IF ( PACKAGE(1:7).EQ.'SPR-KKR' ) IRSHFT = 1
         READ (IOPOT,'(A80)') LINE
         IF ( CHECK ) WRITE (6,'(A)') LINE(1:LEN_TRIM(LINE))
         DO IM = 1,NM
            READ (IOPOT,*) IDUMMY,ARMSH(IM),BRMSH(IM),JRMT(IM),RMT(IM),
     &                     JRWS(IM),RWS(IM)
            IF ( PACKAGE(1:7).EQ.'SPR-KKR' ) JRMT(IM) = JRMT(IM) - 2
            JRWS(IM) = JRWS(IM) - IRSHFT
            IF ( CHECK ) WRITE (6,*) 'IDUMMY,ARMSH(IM),BRMSH(IM),',
     &                               'JRMT(IM),RMT(IM),JRWS(IM),RWS(IM)'
     &                               ,IDUMMY,ARMSH(IM),BRMSH(IM),
     &                               JRMT(IM),RMT(IM),JRWS(IM),RWS(IM)
         END DO
      ELSE
         READ (6,*) STR10,RMESHTYPE
         CALL STOP_MESSAGE(ROUTINE,'mesh-type  not allowed')
      END IF
C ----------------------------------------------------------------------
      IF ( FULLPOT ) THEN
         READ (IOPOT,*) STR10
         DO IM = 1,NM
            READ (IOPOT,*) IDUMMY,JRNS1(IM),JRCRI(IM),NPAN(IM),
     &                     (JRCUT(IPAN,IM),IPAN=1,NPAN(IM))
            JRNS1(IM) = JRNS1(IM) - IRSHFT
            JRCRI(IM) = JRCRI(IM) - IRSHFT
            DO IPAN = 1,NPAN(IM)
               JRCUT(IPAN,IM) = JRCUT(IPAN,IM) - IRSHFT
            END DO
         END DO
      END IF
C ----------------------------------------------------------------------
C
C ======================================================================
      CALL POSFIL(IOPOT,'TYPES     ',IPOS,1)
      READ (IOPOT,'(A80)') LINE
      IF ( CHECK ) WRITE (6,'(A)') LINE(1:LEN_TRIM(LINE))
      DO IT = 1,NT
         CONC(IT) = X_CHEM_T(IT)
         READ (IOPOT,*) IDUMMY,TXT_T(IT),Z(IT),NCORT(IT),NVALT(IT),
     &                  NSEMCORSHLT(IT)
         IF ( CHECK ) WRITE (6,FMT=FMT06) IDUMMY,TXT_T(IT),Z(IT),
     &                       NCORT(IT),NVALT(IT),NSEMCORSHLT(IT)
      END DO
C======================================================================-
C
      IF ( SCFSTATUS(1:5).EQ.'START' ) RETURN
C
C======================================================================-
C
      CALL POSFIL(IOPOT,'POTENTIAL ',IPOS,1)
      DO IT = 1,NT
         IM = IMT(IT)
         IF ( FULLPOT ) THEN
            IRTOP = JRCRI(IM)
         ELSE
            IRTOP = JRWS(IM)
         END IF
         READ (IOPOT,*) STR10,IDUMMY
         IF ( STR10.NE.'TYPE      ' .OR. IDUMMY.NE.IT )
     &        CALL STOP_MESSAGE(ROUTINE,'reading potential failed')
         READ (IOPOT,FMT=FMT07) (RDUMMY,IR=1,IRSHFT),
     &                          (VT(IR,IT),IR=1,IRTOP)
         READ (IOPOT,FMT=FMT07) (RDUMMY,IR=1,IRSHFT),
     &                          (BT(IR,IT),IR=1,IRTOP)
         IF ( CHECK ) WRITE (6,*) 'IT,IM,JTOP,IRSHFT',IT,IM,IRTOP,IRSHFT
C ----------------------------------------------------------------------
         IF ( FULLPOT ) THEN
            DO LM = 1,NLMFPMAX
               KLMFP(LM,IT) = 0
            END DO
            KLMFP(1,IT) = 1
            LMIFP(1,IT) = 1
            READ (IOPOT,*) STR10,NFPT(IT)
            READ (IOPOT,*) STR10,NLMFPT(IT)
C------ allow more LMs in pot file than requested for actual calculation
            IFP = 1
            IFLAG = 0
            DO WHILE ( IFP.LT.NFPT(IT) .AND. IFLAG.EQ.0 )
               READ (IOPOT,*) STR10,LM
               IF ( LM.LE.NLMFPMAX ) THEN
                  IFP = IFP + 1
                  KLMFP(LM,IT) = 1
                  LMIFP(IFP,IT) = LM
                  READ (IOPOT,FMT=FMT07)
     &                  (VNST(IR,LM,IT),IR=JRNS1(IM),IRTOP)
                  READ (IOPOT,FMT=FMT07)
     &                  (BNST(IR,LM,IT),IR=JRNS1(IM),IRTOP)
               ELSE
                  IFLAG = 1
               END IF
            END DO
            IF ( IFLAG.NE.0 ) THEN
               NFPT(IT) = IFP
               NLMFPT(IT) = NLMFPMAX
            END IF
         END IF
C ----------------------------------------------------------------------
         DO WHILE ( STR10.NE.'==========' )
            READ (IOPOT,*) STR10
         END DO
      END DO
C ======================================================================
      CALL POSFIL(IOPOT,'CHARGE    ',IPOS,1)
      DO IT = 1,NT
         READ (IOPOT,*) STR10,IDUMMY
         IF ( STR10.NE.'TYPE      ' .OR. IDUMMY.NE.IT )
     &        CALL STOP_MESSAGE(ROUTINE,'reading charge failed')
         IM = IMT(IT)
         IF ( FULLPOT ) THEN
            IRTOP = JRCRI(IM)
         ELSE
            IRTOP = JRWS(IM)
         END IF
         READ (IOPOT,FMT=FMT07) (RDUMMY,IR=1,IRSHFT),
     &                          (RHOCHR(IR,IT),IR=1,IRTOP)
         READ (IOPOT,FMT=FMT07) (RDUMMY,IR=1,IRSHFT),
     &                          (RHOSPN(IR,IT),IR=1,IRTOP)
         READ (IOPOT,*) STR10
C ----------------------------------------------------------------------
C       RHO2NS added to potfile Nov. 2009 - allow reading older potfiles
         IF ( FULLPOT .AND. STR10.NE.'==========' ) THEN
            READ (IOPOT,*) STR10,IDUMMY
            READ (IOPOT,*) STR10,IDUMMY
            DO LM = 1,NLMFPT(IT)
               IF ( KLMFP(LM,IT).NE.0 ) THEN
                  READ (IOPOT,*) STR10,IDUMMY
                  IF ( STR10.NE.'LM        ' .OR. IDUMMY.NE.LM )
     &                 CALL STOP_MESSAGE(ROUTINE,'reading charge  LM ')
                  READ (IOPOT,FMT=FMT07) (RHO2NS(IR,LM,IT,1),IR=1,IRTOP)
                  READ (IOPOT,FMT=FMT07) (RHO2NS(IR,LM,IT,2),IR=1,IRTOP)
               END IF
            END DO
            DO WHILE ( STR10.NE.'==========' )
               READ (IOPOT,*) STR10
            END DO
         END IF
C ----------------------------------------------------------------------
      END DO
C ======================================================================
      CALL POSFIL(IOPOT,'MOMENTS   ',IPOS,1)
      DO IT = 1,NT
         READ (IOPOT,*) STR10,IDUMMY
         READ (IOPOT,FMT=FMT07) QEL(IT),OBS_T(0,IDOS,IT),
     &                          OBS_T(0,ISMT,IT),OBS_T(0,IOMT,IT),
     &                          OBS_T(0,IHFF,IT)
C
         READ (IOPOT,*) STR10
      END DO
C ======================================================================
      IF ( ORBPOL(1:6).EQ.'BROOKS' ) THEN
C        CALL POSFIL(IOPOT,'VECTOR-POTENTIAL',IPOS,1)
         CALL POSFIL(IOPOT,'VECTOR-POT',IPOS,1)
         DO IT = 1,NT
            READ (IOPOT,*) STR10,IDUMMY
            READ (IOPOT,*) STR10,LOPT(IT)
C ----------------------------------------------------------------------
            IF ( LOPT(IT).GT.0 ) THEN
               IM = IMT(IT)
               DO MS = 1,2
                  READ (IOPOT,FMT=FMT07) (RDUMMY,IR=1,IRSHFT),
     &                  (AOPT(IR,MS,IT),IR=1,JRWS(IM))
               END DO
            END IF
            READ (IOPOT,*) STR10
         END DO
C ======================================================================
      ELSE IF ( BREITINT ) THEN
         CALL POSFIL(IOPOT,'VECPOT-BRE',IPOS,1)
         DO IT = 1,NT
            IM = IMT(IT)
            READ (IOPOT,*) STR10,IDUMMY
            DO MS = -1, + 1
               DO ILA = 1,NLABIMAX
                  READ (IOPOT,FMT=FMT07)
     &                  (ABIT(IR,ILA,MS,IT),IR=1,JRWS(IM))
               END DO
            END DO
            READ (IOPOT,*) STR10
         END DO
C ======================================================================
CC      ELSE IF ( ORBPOL(1:4).EQ.'DMFT' ) THEN
CCC        CALL POSFIL(IOPOT,'SELF-ENERGY',IPOS,1)
CC         CALL POSFIL(IOPOT,'SELF-ENERG',IPOS,1)
CC         DO IT = 1,NT
CC            READ (IOPOT,*) STR10,IDUMMY
CC            READ (IOPOT,*) STR10,LOPT(IT)
CCC ----------------------------------------------------------------------
CC            IF ( LOPT(IT).GT.0 ) THEN
CC            END IF
CC            READ (IOPOT,*) STR10
CC         END DO
C ======================================================================
      END IF
C ======================================================================
C HE 16/02/14     taken out
C ======================================================================
C     CALL POSFIL(IOPOT,'CORE STATES',IPOS,1)
Cc      CALL POSFIL(IOPOT,'CORE STATE',IPOS,0)
Cc      if( ipos .eq. 0 ) then
Cc      DO IT = 1,NT
Cc         READ (IOPOT,*) STR10,IDUMMY
Cc         READ (IOPOT,*) STR10,NCORE(IT)
Cc         IF ( NCORE(IT).GT.0 ) THEN
Cc            READ (IOPOT,*) STR10
Cc            DO IS = 1,NSPIN
Cc               IQ = (IT-1)*NSPIN + IS
Cc               DO IM = 1,NCORE(IT)
Cc                  READ (IOPOT,*) LCORE(IM,IT),ECORE(IM,IQ)
Cc               END DO
Cc            END DO
Cc         END IF
Cc         READ (IOPOT,*) STR10
Cc      END DO
Cc      end if
C======================================================================
C
      CLOSE (IOPOT)
C
99001 FORMAT (2(/,1X,79('W')),/,10X,'WARNING from <POTRD>',/,10X,
     &        'host Madelung potential not found in potfile',2X,A,/,10X,
     &        'VLMMAD_HOST set to 0 !',2(/,1X,79('W')))
99002 FORMAT (2(/,1X,79('W')),/,10X,'WARNING from <POTRD>',/,10X,
     &        'host charge moments not found in potfile',2X,A,/,10X,
     &        'CMNTQ set to 0 !',2(/,1X,79('W')))
      END
C*==potwrspr.f    processed by SPAG 6.70Rc at 08:47 on  8 Mar 2017
      SUBROUTINE POTWRSPR
C   ********************************************************************
C   *                                                                  *
C   *  Driver routine to < POTWR >                                     *
C   *                                                                  *
C   *  ALL output data are taken from modules                          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:WRPOT
      USE MOD_RMESH,ONLY:NM
      USE MOD_LATTICE,ONLY:SYSTEM_TYPE
      USE MOD_TB,ONLY:NREF
      USE MOD_MPI,ONLY:MPI_ID
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C*** End of declarations rewritten by SPAG
C
      IF ( .NOT.WRPOT ) RETURN
      IF ( MPI_ID.NE.0 ) RETURN
C
C=======================================================================
      IF ( SYSTEM_TYPE(1:16).EQ.'EMBEDDED-CLUSTER' ) THEN
C
         CALL CLUPOTWR
C
         RETURN
C
      END IF
C=======================================================================
C
C ----------------------------------------------------------------------
C  set tb-kkr variables: IREFQ, VREF, RMTREF
C
      NREF = NM
C
C ======================================================================
C
      CALL POTWR
C
C ======================================================================
C
      END
C*==potwr.f    processed by SPAG 6.70Rc at 08:47 on  8 Mar 2017
      SUBROUTINE POTWR
C   ********************************************************************
C   *                                                                  *
C   *   write out the SCF-potential for  POTFMT > 5                    *
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
     &    NSEMCORSHLT,IMT,X_CHEM_T,BT,BEXT,OBS_T,MTET_T,MPHI_T,MGAM_T
      USE MOD_ANGMOM,ONLY:NLABIMAX,NSPIN,IDOS,ISMT,IOMT,IHFF
      USE MOD_TB,ONLY:NREF,RMTREF,VREF,IREFQ
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='POTWR')
C
C Local variables
C
      CHARACTER*40 FMT02,FMT03,FMT04,FMT05,FMT06,FMT07,FMT08,FMT12
      INTEGER IC,ILA,IM,IO,IPAN,IQ,IR,IRSHFT,IRTOP,IT,J,LM,MS
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
         FMT02 = '(I10,3I10,10(I6,F15.12))'
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
C-----------------------------------------------------------------------
C
      CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,POTFIL(1:LPOTFIL))
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
      WRITE (IOTMP,99017) 'NT        ',NT
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
         WRITE (IOTMP,FMT=FMT02) IQ,IREFQ(IQ),IDUMMY,NOQ(IQ),
     &                           (ITOQ(IO,IQ),X_CHEM_T(ITOQ(IO,IQ)),
     &                           IO=1,NOQ(IQ))
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
      DO IT = 1,NT
c modified by XJQ: BUG bug !!!
c         WRITE (IOTMP,FMT=FMT03) IQ,MTET_T(IT),MPHI_T(IT),MGAM_T(IT)
         WRITE (IOTMP,FMT=FMT03) IT,MTET_T(IT),MPHI_T(IT),MGAM_T(IT)
c end-mod-xjq
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
         CALL STOP_MESSAGE(ROUTINE,'mesh-type  not allowed')
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
      DO IT = 1,NT
         WRITE (IOTMP,FMT=FMT06) IT,TXT_T(IT),Z(IT),NCORT(IT),NVALT(IT),
     &                           NSEMCORSHLT(IT)
      END DO
C=======================================================================
      WRITE (IOTMP,99001) 'POTENTIAL'
      DO IT = 1,NT
         IM = IMT(IT)
         IF ( FULLPOT ) THEN
            IRTOP = JRCRI(IM)
         ELSE
            IRTOP = JRWS(IM)
         END IF
         WRITE (IOTMP,99017) 'TYPE      ',IT
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
      DO IT = 1,NT
         WRITE (IOTMP,99017) 'TYPE      ',IT
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
      DO IT = 1,NT
         WRITE (IOTMP,99017) 'TYPE      ',IT
         WRITE (IOTMP,FMT=FMT07) QEL(IT),OBS_T(0,IDOS,IT),
     &                           OBS_T(0,ISMT,IT),OBS_T(0,IOMT,IT),
     &                           OBS_T(0,IHFF,IT)
C
         WRITE (IOTMP,99012)
      END DO
C ======================================================================
      IF ( ORBPOL(1:6).EQ.'BROOKS' ) THEN
         WRITE (IOTMP,99001) 'VECTOR-POTENTIAL'
         DO IT = 1,NT
            WRITE (IOTMP,99017) 'TYPE      ',IT
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
         DO IT = 1,NT
            WRITE (IOTMP,99017) 'TYPE      ',IT
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
         DO IT = 1,NT
            WRITE (IOTMP,99017) 'TYPE      ',IT
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
Cc         WRITE (IOTMP,99018) 'TYPE      ',IT
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
C*==potio_swsrmt.f    processed by SPAG 6.70Rc at 08:47 on  8 Mar 2017
      SUBROUTINE POTIO_SWSRMT(NQ,NM,IMQ,SWS,JRMT,RMT,JRWS,RWS,ABAS,QBAS,
     &                        ALAT,RMESHTYPE,SCFSTATUS,NQMAX)
C   ********************************************************************
C   *                                                                  *
C   *   fix the average Wigner-Seitz and mufin tin radii               *
C   *                                                                  *
C   *   called by <POTRD2>  and <POTRDSPR>                             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='POTIO_SWSRMT')
C
C Dummy arguments
C
      REAL*8 ALAT,SWS
      INTEGER NM,NQ,NQMAX
      CHARACTER*12 RMESHTYPE
      CHARACTER*10 SCFSTATUS
      REAL*8 ABAS(3,3),QBAS(3,NQMAX),RMT(NM),RWS(NM)
      INTEGER IMQ(NQMAX),JRMT(NM),JRWS(NM)
C
C Local variables
C
      REAL*8 AUX,DIJ,DQBAS(3),DRVEC(3),DX(:),R(:,:),VUC,VWS
      REAL*8 DNRM2
      INTEGER I1,I2,I3,IM,IQ,IRMTALG,ITOP,JM,JQ
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DX,R
C
      ALLOCATE (DX(NM),R(1,NM))
C
C dummy setting for R
      R(1,1:NM) = 1D-6
C
      CALL WRSUBNAMHEAD('POTIO_SWSRMT',12,' ',0)
      IF ( IPRINT.GT.0 ) WRITE (6,99001)
C
C=======================================================================
C
      SWS = 0.0D0
      DO IQ = 1,NQ
         SWS = SWS + RWS(IMQ(IQ))**3
      END DO
      SWS = (SWS/DBLE(NQ))**(1.0D0/3.0D0)
      VWS = SWS**3*4D0*PI/3D0
C
C ---------------------------------------------------- check radial mesh
      IF ( SCFSTATUS(1:5).EQ.'START' ) THEN
C --------------------------------------------------- Wigner-Seitz radii
C
         CALL RVECSPAT(ABAS(1,1),ABAS(1,2),ABAS(1,3),VUC,1)
         VUC = ABS(VUC)*ALAT**3
C
c modified by XJQ: in any case, dx is needed
c         IF ( ABS(1D0-NQ*VWS/VUC).GT.1D-14 ) THEN
c end-mod-xjq
            IF ( ABS(1D0-NQ*VWS/VUC).GT.1D-8 ) WRITE (6,99006) VUC,
     &           NQ*VWS,ABS(VUC-NQ*VWS)
C
            AUX = (VUC/(NQ*VWS))**(1D0/3D0)
            DO IM = 1,NM
               RWS(IM) = AUX*RWS(IM)
               DX(IM) = LOG(RWS(IM)/R(1,IM))/DBLE(JRWS(IM)-1)
            END DO
C
            SWS = 0.0D0
            DO IQ = 1,NQ
               SWS = SWS + RWS(IMQ(IQ))**3
            END DO
            SWS = (SWS/DBLE(NQ))**(1.0D0/3.0D0)
            VWS = SWS**3*4D0*PI/3D0
            IF ( ABS(1D0-NQ*VWS/VUC).GT.1D-8 ) THEN
               WRITE (6,99006) VUC,NQ*VWS,ABS(VUC-NQ*VWS)
               CALL STOP_MESSAGE(ROUTINE,'update of SWS failed')
            END IF
C
c         END IF
C
C ----------------------------------------------------- muffin-tin radii
C adjust  R_mt   if  1 <= JRMT <= 2  or  R_mt <= 0
C algorithm  JRMT = IRMTALG = 1:  scale R_ws
C                             2:  take midpoint
C
C         DO IM = 1,NM
C            IF ( JRMT(IM).GE.1 .AND. JRMT(IM).LE.2 ) IRMTALG = JRMT(IM)
C            IF ( RMT(IM).LT.1D-6 ) IRMTALG = 1
C         END DO
C
         IRMTALG = 2
C
         IF ( IRMTALG.NE.0 ) THEN
            DO IM = 1,NM
               RMT(IM) = RWS(IM)
            END DO
         END IF
C
         ITOP = 2
C
         DO I1 = -ITOP,ITOP
            DO I2 = -ITOP,ITOP
               DO I3 = -ITOP,ITOP
C
                  CALL RVECLCIB(I1,I2,I3,ABAS,DRVEC)
C
                  DO IQ = 1,NQ
                     IM = IMQ(IQ)
                     DO JQ = 1,NQ
                        JM = IMQ(JQ)
                        DQBAS(:) = DRVEC(:) + QBAS(:,IQ) - QBAS(:,JQ)
                        DIJ = DNRM2(3,DQBAS,1)*ALAT
                        IF ( DIJ.GT.1D-8 ) THEN
                           IF ( IRMTALG.EQ.1 ) THEN
                              AUX = DIJ/(RWS(IM)+RWS(JM))
                              RMT(IM) = MIN(RMT(IM),AUX*RWS(IM))
                              RMT(JM) = MIN(RMT(JM),AUX*RWS(JM))
                           ELSE IF ( IRMTALG.EQ.2 ) THEN
                              AUX = DIJ/2D0
                              RMT(IM) = MIN(RMT(IM),AUX)
                              RMT(JM) = MIN(RMT(JM),AUX)
                           END IF
                           IF ( DIJ-RMT(IM)-RMT(JM).LT.-1D-8 ) THEN
                              WRITE (6,*) '  IRMTALG  ',IRMTALG
                              WRITE (6,*) '  IM  ',IM,'  RMT:',RMT(IM)
                              WRITE (6,*) '  JM  ',JM,'  RMT:',RMT(JM)
                              WRITE (6,*) '  DIJ          ',DIJ
                              WRITE (6,*) '  RMT(I)+RMT(J)',RMT(IM)
     &                               + RMT(JM)
                              CALL STOP_MESSAGE(ROUTINE,
     &                           'mesh inconsistent')
                           END IF
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END DO
C
         WRITE (6,99002) IRMTALG
         IF ( IRMTALG.EQ.0 ) THEN
            WRITE (6,99003)
         ELSE IF ( IRMTALG.EQ.1 ) THEN
            WRITE (6,99004)
         ELSE IF ( IRMTALG.EQ.2 ) THEN
            WRITE (6,99005)
         ELSE
            CALL STOP_MESSAGE(ROUTINE,'IRMTALG > 2 ')
         END IF
C
         DO IM = 1,NM
C
            IF ( ABS(RMT(IM)-RWS(IM)).LT.1D-6 ) THEN
               RMT(IM) = 0.97D0*RWS(IM)
               WRITE (6,99007) IM
            END IF
C
            IF ( RMESHTYPE.EQ.'EXPONENTIAL ' ) THEN
               JRMT(IM) = INT(LOG(RMT(IM)/R(1,IM))/DX(IM)) + 1
            ELSE
               CALL STOP_MESSAGE(ROUTINE,
     &                           'JRMT to be done for RMESHTYPE')
            END IF
         END DO
C
      END IF
C
      RETURN
C
C ======================================================================
99001 FORMAT (10X,'fix the muffin tin radii according to control',
     &        ' variable IRMTALG',/)
99002 FORMAT (10X,'IRMTALG =',I3,'    read from potential file ',/)
99003 FORMAT (10X,'=> muffin tin radii taken from potential file ',/)
99004 FORMAT (10X,'=> muffin tin radii obtained scaling WS-radii',/)
99005 FORMAT (10X,'=> muffin tin radii = half interatomic distance',/)
99006 FORMAT (2(/,1X,79('W')),/,10X,'WARNING from <POTIO_SWSRMT>',/,10X,
     &        'radial mesh parameters in potfile  inconsistent '/,10X,
     &        'a1*(a2xa3) * a^3     ',F15.8,/,10X,
     &        'NQ * SWS^3 * 4*PI/3  ',F15.8,/,10X,
     &        'DEVIATION            ',F15.8,/,10X,
     &        'radial mesh parameters will be modified',2(/,1X,79('W')),
     &        /)
99007 FORMAT (10X,'INFO: for IM =',I3,'    R_mt = 0.97 * R_WS')
      END
C*==wrsubnamhead.f    processed by SPAG 6.70Rc at 08:47 on  8 Mar 2017
      SUBROUTINE WRSUBNAMHEAD(SUBNAME,LSUB,TXT,LTXT)
C   ********************************************************************
C   *                                                                  *
C   *   write subroutine name as a header  + a text line               *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER LSUB,LTXT
      CHARACTER*(*) SUBNAME,TXT
C
C Local variables
C
      CHARACTER*80 FMTHEAD
      INTEGER LL
C
C*** End of declarations rewritten by SPAG
C
      LL = (80-LSUB-2)/2
      WRITE (FMTHEAD,99001) LL
C
      WRITE (6,FMTHEAD) SUBNAME(1:LSUB)
C
      IF ( LTXT.GT.0 ) WRITE (6,99002) TXT(1:LTXT)
C
99001 FORMAT ('(//,1X,79(''*''),/,',I2,
     &        '('' ''),''<'',A ,''>'',/,1X,79(''*''),/)')
99002 FORMAT (/,10X,A,/)
C
      END
