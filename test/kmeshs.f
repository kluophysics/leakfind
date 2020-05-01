C*==kmeshs.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE KMESHS(IOTMP,NSYM,SYMACCEPTED,MROTK,IPRINT,NKPTS0,BBAS,
     &                  DATSET,LDATSET)
C   ********************************************************************
C   *                                                                  *
C   *  create a regular k-mesh by subdividing a parallelepiped         *
C   *  spanned by the reciprocal basis vectors  BBAS                   *
C   *  at the beginning around NKPTS0 k-points are created             *
C   *  this number is reduced making use of the symmetry operations    *
C   *                                                                  *
C   * 01/08/2000 HE                                                    *
C   ********************************************************************
C
      USE MOD_LATTICE,ONLY:SYSTEM_DIMENSION,SUB_SYSTEM
      USE MOD_SYMMETRY,ONLY:NSYMMAX
      IMPLICIT NONE
C*--KMESHS17
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='KMESHS')
      LOGICAL CHECK
      PARAMETER (CHECK=.FALSE.)
C
C Dummy arguments
C
      CHARACTER*80 DATSET
      INTEGER IOTMP,IPRINT,LDATSET,NKPTS0,NSYM
      REAL*8 BBAS(3,3),MROTK(3,3,NSYMMAX)
      LOGICAL SYMACCEPTED(NSYMMAX)
C
C Local variables
C
      REAL*8 BG(3,3),BGINV(3,3),BGL(3),BGMAT(3,3),BGP(3,3),BV(3),
     &       BZSTEP(3,3),CF(3),KTAB(:,:),PROBGL,RSUM,S,SCLFAC,SCLWGT,
     &       WKTAB(:),WMAX,XDIM
      REAL*8 DDOT
      INTEGER*8 I,I0,I1,I2,I3,IA_ERR,IK,IKREDUC,IKTAB(3),
     &          IKTAB_GEN_KREDUC(:),IND2(3),IROT,IS,ISYM_GEN_KREDUC(:),
     &          ITRY,IW,IWT,IX,J,LIN,LIN2,NBGL(3),NBGLMAX,
     &          NBGP(3,3,NSYMMAX),NDIM,NKTAB,NKTABMAX,NKTABMAXSYM,
     &          NKTAB_KREDUC,NLIN,NSYMACCEPTED,NVEC_BZSTEP_KREDUC(:,:)
      INTEGER*8 INT8
      LOGICAL KNOTDONE(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE KNOTDONE,KTAB,WKTAB
      ALLOCATABLE NVEC_BZSTEP_KREDUC,IKTAB_GEN_KREDUC,ISYM_GEN_KREDUC
C
      NKTABMAX = INT8(NKPTS0*1.2D0)
      IF ( NKTABMAX.LT.0 ) THEN
         WRITE (6,'("NKTABMAX",i50)') NKTABMAX
         CALL STOP_MESSAGE(ROUTINE,'Negative NKTABMAX - overflow??')
      END IF
C
      ITRY = 0
C
      NSYMACCEPTED = 0
      DO IROT = 1,NSYM
         IF ( SYMACCEPTED(IROT) ) NSYMACCEPTED = NSYMACCEPTED + 1
      END DO
C     ---------------- find a save estimate for the array size for IRRBZ
      NKTABMAXSYM = INT8((NKTABMAX/NSYMACCEPTED)*1.2D0)
C
 100  CONTINUE
      ITRY = ITRY + 1
C
C     ---------------- avoid problems for case of having small amount of
C     ---------------- kpoints when doing the estimate
      IF ( NKTABMAXSYM.LT.1000 ) NKTABMAXSYM = NKTABMAX
C
      IF ( CHECK ) THEN
         ALLOCATE (NVEC_BZSTEP_KREDUC(3,NKTABMAX))
         ALLOCATE (IKTAB_GEN_KREDUC(NKTABMAX),ISYM_GEN_KREDUC(NKTABMAX))
      END IF
      IF ( ALLOCATED(KTAB) ) DEALLOCATE (KTAB,WKTAB)
      ALLOCATE (KTAB(3,NKTABMAXSYM),WKTAB(NKTABMAXSYM))
      ALLOCATE (KNOTDONE(NKPTS0),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'KNOTDONE')
C
C=======================================================================
C                         system dimension
C=======================================================================
C
      IF ( SYSTEM_DIMENSION(1:2).EQ.'3D' ) THEN
         NDIM = 3
      ELSE IF ( SUB_SYSTEM(1:6).EQ.'L-BULK' ) THEN
         NDIM = 3
      ELSE IF ( SUB_SYSTEM(1:6).EQ.'R-BULK' ) THEN
         NDIM = 3
      ELSE IF ( SUB_SYSTEM(1:6).EQ.'I-ZONE' ) THEN
         NDIM = 2
      END IF
C
C=======================================================================
C                          construct the BZ-mesh
C=======================================================================
C                                BG    basis vectors in reciprocal space
C
      BG(1:3,1:3) = BBAS(1:3,1:3)
C
C--------------------------------------- for 2D-systems:  B(3) = (0,0,1)
      IF ( NDIM.EQ.2 ) THEN
         BG(1:2,3) = 0D0
         BG(3,3) = 1D0
      END IF
C
C--------------------------------------------------- fix the points grid
C     allow an individual grid only for basis vectors that are
C     orthogonal to all other basis vectors -
C     otherwise take max number of divisions
C
      PROBGL = 1D0
      DO I = 1,NDIM
         BGL(I) = SQRT(DDOT(3,BG(1,I),1,BG(1,I),1))
         PROBGL = PROBGL*BGL(I)
      END DO
C
      XDIM = 1D0/DBLE(NDIM)
      NBGLMAX = 0
      NBGL(3) = 1
      DO I = 1,NDIM
         NBGL(I) = MAX(1,NINT(BGL(I)*(DBLE(NKPTS0)/PROBGL)**XDIM))
         NBGLMAX = MAX(NBGLMAX,NBGL(I))
      END DO
C
      DO J = 1,NDIM
         DO I = 1,NDIM
            BGMAT(I,J) = DDOT(3,BG(1,I),1,BG(1,J),1)
         END DO
      END DO
Ci
      DO I = 1,NDIM
         DO J = 1,NDIM
            IF ( (I.NE.J) .AND. ABS(BGMAT(I,J)).GT.1D-8 ) NBGL(I)
     &           = NBGLMAX
         END DO
      END DO
C
 200  CONTINUE
      LIN = NBGL(1)*NBGL(2)*NBGL(3)
      IF ( LIN.GT.NKPTS0 ) THEN
         DO I = 1,3
            NBGL(I) = MAX(1,NBGL(I)-1)
         END DO
         IF ( IPRINT.GT.0 ) THEN
            WRITE (6,*) 'WARNING from <KMESHS>'
            WRITE (6,*) 'available array size exceeded '
            WRITE (6,*) 'NKPTS0 =',NKPTS0
            WRITE (6,*) 'NBGL reduced to ',NBGL
         END IF
         GOTO 200
      END IF
C
      LIN = 0
      DO I1 = 0,NBGL(1) - 1
         DO I2 = 0,NBGL(2) - 1
            DO I3 = 0,NBGL(3) - 1
               LIN = LIN + 1
               KNOTDONE(LIN) = .TRUE.
               LIN2 = I1*NBGL(3)*NBGL(2) + I2*NBGL(3) + I3 + 1
               IF ( LIN.NE.LIN2 ) WRITE (6,*) '###LIN-INDEX#######',I1,
     &              I2,I3,LIN,LIN2
            END DO
         END DO
      END DO
      NLIN = LIN
C
C-----------------------------------------------------------------------
C                                          steps along the basis vectors
      DO J = 1,3
         SCLFAC = 1D0/DBLE(NBGL(J))
         BZSTEP(1:3,J) = BG(1:3,J)*SCLFAC
      END DO
C
      DO J = 1,3
         DO I = 1,3
            BGMAT(I,J) = DDOT(3,BZSTEP(1,I),1,BZSTEP(1,J),1)
         END DO
      END DO
C
C---------------------- check for 2D-systems:  B(1)*B(3) = B(2)*B(3) = 0
      IF ( NDIM.EQ.2 ) THEN
         IF ( ABS(BGMAT(1,3)).GT.1D-6 .OR. ABS(BGMAT(2,3)).GT.1D-6 )
     &        THEN
            WRITE (*,*) ' ***************************************'
            WRITE (*,*) ' TEST  B(1)*B(3) = B(2)*B(3) = 0  failed'
            WRITE (*,'(3F10.5)') BGMAT
            CALL STOP_MESSAGE(ROUTINE,'check for 2D-systems')
         END IF
      END IF
C
      CALL RINVGJ(BGINV,BGMAT,3,3)
C
C-----------------------------------------------------------------------
C          rotate the 3 step vectors   BZSTEP  --  it should be possible
C          to express the rotated vectors  B(j) in terms of the old ones
C     B(i) = SUM(j) BZSTEP(j) * n(j,i)  with integer coefficients n(j,i)
C
      DO IROT = 1,NSYM
         IF ( SYMACCEPTED(IROT) ) THEN
            IF ( IPRINT.GT.0 ) WRITE (6,99002) IROT
C
            DO I = 1,3
               CALL DGEMV('N',3,3,1D0,MROTK(1,1,IROT),3,BZSTEP(1,I),1,
     &                    0D0,BGP(1,I),1)
C
               DO J = 1,3
                  BV(J) = DDOT(3,BZSTEP(1,J),1,BGP(1,I),1)
               END DO
               CALL DGEMV('N',3,3,1D0,BGINV,3,BV,1,0D0,CF,1)
C
               DO J = 1,3
                  IF ( ABS(NINT(CF(J))-CF(J)).GT.1D-8 ) THEN
                     WRITE (6,99010) I,J,CF(J)
                     CALL STOP_MESSAGE(ROUTINE,
     &                                 'ABS(NINT(CF(J))-CF(J))>1D-8')
                  END IF
                  NBGP(J,I,IROT) = NINT(CF(J))
               END DO
               IF ( IPRINT.GT.0 ) WRITE (6,99003) I,(BGP(J,I),J=1,3),BV,
     &              CF,(NBGP(J,I,IROT),J=1,3)
C
            END DO
C
         END IF
      END DO
C
C=======================================================================
C              scan the full k-mesh and remove equivalent points
C=======================================================================
C
      WMAX = 0D0
      NKTAB = 0
      NKTAB_KREDUC = 0
      LIN = 0
      DO I1 = 0,NBGL(1) - 1
         IKTAB(1) = I1
         DO I2 = 0,NBGL(2) - 1
            IKTAB(2) = I2
            DO I3 = 0,NBGL(3) - 1
               IKTAB(3) = I3
               LIN = LIN + 1
               IF ( KNOTDONE(LIN) ) THEN
                  NKTAB = NKTAB + 1
                  IWT = 0
C
C-----------------------------------------------------------------------
C scan symmetry operations and drop all points equivalent to present one
C
                  DO IROT = 1,NSYM
                     IF ( SYMACCEPTED(IROT) ) THEN
C
C               rotate k-vector  LIN   and transform into parallelepiped
C          ROT k = SUM(i) m(i) ROT bzstep(i)
C                = SUM(i) m(i) SUM(j) n(i,j) bzstep(j)
C                = SUM(j) [SUM(i) m(i) n(i,j)] bzstep(j)
C
                        DO J = 1,3
                           IS = 0
                           DO I = 1,3
                              IS = IS + IKTAB(I)*NBGP(J,I,IROT)
                           END DO
                           IS = MOD(IS,NBGL(J))
                           IF ( IS.LT.0 ) IS = IS + NBGL(J)
                           IND2(J) = IS
                        END DO
C
                        LIN2 = IND2(1)*NBGL(3)*NBGL(2) + IND2(2)*NBGL(3)
     &                         + IND2(3) + 1
C
                        IF ( KNOTDONE(LIN2) ) THEN
C------------------------------- drop equivalent k-point from to-do list
                           KNOTDONE(LIN2) = .FALSE.
C------------------------------------ increase weight of present k-point
                           IWT = IWT + 1
C------------------------ store information for the whole reducible list
                           IF ( CHECK ) THEN
                              NKTAB_KREDUC = NKTAB_KREDUC + 1
                              IKREDUC = NKTAB_KREDUC
                              IF ( IKREDUC.LT.NKTABMAX ) THEN
                                 NVEC_BZSTEP_KREDUC(1:3,IKREDUC)
     &                              = IND2(1:3)
                                 IKTAB_GEN_KREDUC(IKREDUC) = NKTAB
                                 ISYM_GEN_KREDUC(IKREDUC) = IROT
                              END IF
                           END IF
                        END IF
C
                     END IF
                  END DO
C-----------------------------------------------------------------------
C
                  IF ( NKTAB.LE.NKTABMAXSYM ) THEN
                     DO I = 1,3
                        RSUM = 0D0
                        DO J = 1,3
                           RSUM = RSUM + BZSTEP(I,J)*IKTAB(J)
                        END DO
                        KTAB(I,NKTAB) = RSUM
                     END DO
                     WKTAB(NKTAB) = DBLE(IWT)
                     WMAX = MAX(WKTAB(NKTAB),WMAX)
                  END IF
C
               END IF
            END DO
         END DO
      END DO
C
      WRITE (6,99004) NKPTS0,NLIN,NBGL,NKTAB,NKTABMAXSYM
      IF ( NKTAB.GT.NKTABMAXSYM ) THEN
         WRITE (6,99011)
         IF ( ITRY.EQ.1 ) THEN
            IF ( CHECK ) DEALLOCATE (NVEC_BZSTEP_KREDUC,IKTAB_GEN_KREDUC
     &                               ,ISYM_GEN_KREDUC)
            DEALLOCATE (KTAB,WKTAB,KNOTDONE)
            NKTABMAXSYM = NKTAB
            GOTO 100
         ELSE
            CALL STOP_MESSAGE(ROUTINE,'ITRY > 1')
         END IF
      END IF
C
      IF ( CHECK ) THEN
         IF ( NKTAB_KREDUC.GT.NKTABMAX ) CALL STOP_MESSAGE(ROUTINE,
     &        '<KMESHS>: array size exceeded -->>> increase  NKTABMAX')
         IF ( NKTAB_KREDUC.NE.NLIN )
     &         CALL STOP_MESSAGE(ROUTINE,'NKTAB_KREDUC.NE.NLIN')
      END IF
C
      RSUM = 0D0
      DO IK = 1,NKTAB
         KTAB(3,IK) = KTAB(3,IK) + 1D-9
         RSUM = RSUM + WKTAB(IK)
      END DO
C
      IF ( ABS(NLIN-RSUM).GT.1D-6 ) THEN
         WRITE (6,99001) RSUM,NLIN
         CALL STOP_MESSAGE(ROUTINE,'ABS(NLIN-RSUM).GT.1D-6')
      END IF
C
C-----------------------------------------------------------------------
C                  create input for the  rasmol  program
C-----------------------------------------------------------------------
      IF ( IPRINT.GT.3 ) THEN
C
         I0 = 16
         IW = 80
         S = 15D0
         OPEN (IW,FILE='rasmol_kmesh.pdb')
C
         WRITE (IW,FMT=99005) 'for '//DATSET(1:(LDATSET-1))
C
C
C------------------------------------------ specify corners of unit cell
C
         WRITE (IW,FMT=99008) 1,1,0D0,0D0,0D0,0D0
         DO I = 1,3
            WRITE (IW,FMT=99008) (I+1),(I+1),S*BBAS(1,I),S*BBAS(2,I),
     &                           S*BBAS(3,I),0D0
         END DO
C
         WRITE (IW,FMT=99008) 5,5,S*(BBAS(1,1)+BBAS(1,2)),
     &                        S*(BBAS(2,1)+BBAS(2,2)),
     &                        S*(BBAS(3,1)+BBAS(3,2)),0D0
         WRITE (IW,FMT=99008) 6,6,S*(BBAS(1,1)+BBAS(1,3)),
     &                        S*(BBAS(2,1)+BBAS(2,3)),
     &                        S*(BBAS(3,1)+BBAS(3,3)),0D0
         WRITE (IW,FMT=99008) 7,7,S*(BBAS(1,3)+BBAS(1,2)),
     &                        S*(BBAS(2,3)+BBAS(2,2)),
     &                        S*(BBAS(3,3)+BBAS(3,2)),0D0
         WRITE (IW,FMT=99008) 8,8,S*(BBAS(1,1)+BBAS(1,2)+BBAS(1,3)),
     &                        S*(BBAS(2,1)+BBAS(2,2)+BBAS(2,3)),
     &                        S*(BBAS(3,1)+BBAS(3,2)+BBAS(3,3)),0D0
C
C-------------------- specify corners of cube with edge length 2*PI/ALAT
C
         WRITE (IW,FMT=99008) 9,9,S*0D0,S*0D0,S*0D0,3D0
         WRITE (IW,FMT=99008) 10,10,S*1D0,S*0D0,S*0D0,3D0
         WRITE (IW,FMT=99008) 11,11,S*0D0,S*1D0,S*0D0,3D0
         WRITE (IW,FMT=99008) 12,12,S*0D0,S*0D0,S*1D0,3D0
         WRITE (IW,FMT=99008) 13,13,S*1D0,S*1D0,S*0D0,3D0
         WRITE (IW,FMT=99008) 14,14,S*1D0,S*0D0,S*1D0,3D0
         WRITE (IW,FMT=99008) 15,15,S*0D0,S*1D0,S*1D0,3D0
         WRITE (IW,FMT=99008) 16,16,S*1D0,S*1D0,S*1D0,3D0
C
         SCLWGT = 3D0/WMAX
C
         I = I0
         DO IK = 1,NKTAB
            I = I + 1
            WRITE (IW,FMT=99007) I,I,(S*KTAB(J,IK),J=1,3),
     &                           SCLWGT*WKTAB(IK)
         END DO
C
         WRITE (IW,FMT=99009)
         CLOSE (IW)
C
C--------------------------------------------------- write RASMOL script
C
         OPEN (IW,FILE='rasmol_kmesh.ras')
         WRITE (IW,*) 'load ''rasmol_kmesh.pdb'' '
         WRITE (IW,*) 'color temperature'
         WRITE (IW,*) 'select 1-16'
         WRITE (IW,*) 'cpk 0'
         WRITE (IW,*) 'select 17-',(NKTAB+16)
         WRITE (IW,*) 'cpk 100'
         WRITE (IW,*) 'select all '
         WRITE (IW,*) 'set axes on '
C
         CLOSE (IW)
C
         WRITE (6,99006)
C
      END IF
C
C=======================================================================
C                   write results for rereading by main program
C=======================================================================
C
      CALL OPEN_IOTMP_SCRATCH(ROUTINE,IOTMP)
      WRITE (IOTMP) NKTAB
      WRITE (IOTMP) ((KTAB(IX,IK),IX=1,3),WKTAB(IK),IK=1,NKTAB)
      IF ( CHECK ) THEN
         WRITE (IOTMP) BZSTEP(1:3,1:3)
         WRITE (IOTMP) NKTAB_KREDUC
         WRITE (IOTMP) (NVEC_BZSTEP_KREDUC(1:3,IKREDUC),IKTAB_GEN_KREDUC
     &                 (IKREDUC),ISYM_GEN_KREDUC(IKREDUC),IKREDUC=1,
     &                 NKTAB_KREDUC)
      END IF
C
C-----------------------------------------------------------------------
      DEALLOCATE (KNOTDONE,KTAB,WKTAB,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC')
C
      RETURN
99001 FORMAT (/,10X,'trouble in <KMESHS>:',/,10X,'SUM(WK)',F8.1,/,10X,
     &        'NLIN   ',I6,/)
99002 FORMAT (/,10X,'rotated BZSTEP  for IROT=',I3)
99003 FORMAT (10X,I3,3F7.3,2X,3F7.3,2X,3F7.3,2X,3I3)
99004 FORMAT (/,10X,'number of k-vectors to start  NKPTS0',I12,3X,
     &        'NLIN',I12,/,10X,'mesh parameters               NBGL',4X,
     &        3I5,/,10X,'number of k-vectors created   NKTAB ',I12,/,
     &        10X,'array size (NKTABMAX)       ',I20,/)
99005 FORMAT ('HEADER    k-mesh                 ',A,/,
     &        'SOURCE    SPRKKR - program       ',/,
     &        'AUTHOR    H. Ebert               ',/,
     &        'REMARK    None                   ')
99006 FORMAT (/,10X,'k-mesh stored in rasmol data-file ',
     &        ' rasmol_kmesh.pdb',/,10X,
     &        'view via:   rasmol  -script rasmol_kmesh.ras')
99007 FORMAT ('ATOM  ',I5,'          ',I5,'    ',3F8.3,'  0.00',F8.3)
99008 FORMAT ('HETATM',I5,'          ',I5,'    ',3F8.3,'  0.00',F8.3)
99009 FORMAT ('CONECT    1    2',/,'CONECT    1    3',/,
     &        'CONECT    1    4',/,'CONECT    5    2',/,
     &        'CONECT    5    3',/,'CONECT    5    8',/,
     &        'CONECT    6    2',/,'CONECT    6    4',/,
     &        'CONECT    6    8',/,'CONECT    7    3',/,
     &        'CONECT    7    4',/,'CONECT    7    8',/,
     &        'CONECT    9   10',/,'CONECT    9   11',/,
     &        'CONECT    9   12',/,'CONECT   13   10',/,
     &        'CONECT   13   11',/,'CONECT   13   16',/,
     &        'CONECT   14   10',/,'CONECT   14   12',/,
     &        'CONECT   14   16',/,'CONECT   15   11',/,
     &        'CONECT   15   12',/,'CONECT   15   16',/,'END')
99010 FORMAT (10X,'ERROR: in <KMESHS> CF should be integer: ',I3,I3,
     &        F10.5)
99011 FORMAT (/,10X,'<KMESHS>: array size NKTABMAXSYM exceeded ',/)
      END
