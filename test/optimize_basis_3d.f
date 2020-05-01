C*==optimize_basis_3d.f    processed by SPAG 6.70Rc at 13:55 on 27 Jan 2017
      SUBROUTINE OPTIMIZE_BASIS_3D(ABAS,QBAS,NQ,NQMAX,IFLAG,IOTMP)
C   ********************************************************************
C   *                                                                  *
C   *  optimize the list of basis vectors   ->q_i                      *
C   *  to have the largest interatomic distance |->q_i - ->q_j|        *
C   *  as small as possible                                            *
C   *                                                                  *
C   *  IFLAG = 0:   no changes done                                    *
C   *          1:   changes done                                       *
C   *         -1:   changes done - result questionable                 *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 TOL
      PARAMETER (TOL=1D-10)
C
C Dummy arguments
C
      INTEGER IFLAG,IOTMP,NQ,NQMAX
      REAL*8 ABAS(3,3),QBAS(3,NQMAX)
C
C Local variables
C
      LOGICAL CHANGES_MADE
      REAL*8 DNRM2
      REAL*8 DQIJ,DQIJMAX_NEW,DQIJMAX_OLD,DQIJMAX_Q,DQIPJMAX,
     &       QBAS_OLD(:,:),QIJVEC(3),QIPVEC(1:3)
      INTEGER I,I1,I123TOP,I1MIN,I2,I2MIN,I3,I3MIN,IQ,JQ
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE QBAS_OLD
C
      IFLAG = 0
C
      IF ( NQ.EQ.1 ) RETURN
C
      ALLOCATE (QBAS_OLD(3,NQ))
C
      QBAS_OLD(1:3,1:NQ) = QBAS(1:3,1:NQ)
C
C-----------------------------------------------------------------------
C          find max. distance for old set of basis vectors
C-----------------------------------------------------------------------
C
      DQIJMAX_OLD = 0D0
      DO IQ = 1,NQ
         DO JQ = IQ + 1,NQ
            QIJVEC(1:3) = QBAS(1:3,IQ) - QBAS(1:3,JQ)
            DQIJ = DNRM2(3,QIJVEC,1)
            DQIJMAX_OLD = MAX(DQIJ,DQIJMAX_OLD)
         END DO
      END DO
C
C-----------------------------------------------------------------------
C          find MINIMUM for max. distance for set of basis vectors
C-----------------------------------------------------------------------
C
      DQIJMAX_NEW = 0D0
      CHANGES_MADE = .FALSE.
      I123TOP = 2
C----------------------------------------------------- scan all sites IQ
      DO IQ = 2,NQ
C
         DQIJMAX_Q = 1D+20
C
C------------------------------------------- scan neighboring unit cells
         DO I1 = I123TOP, - I123TOP, - 1
            DO I2 = I123TOP, - I123TOP, - 1
               DO I3 = I123TOP, - I123TOP, - 1
C
                  QIPVEC(1:3) = QBAS(1:3,IQ) + I1*ABAS(1:3,1)
     &                          + I2*ABAS(1:3,2) + I3*ABAS(1:3,3)
C
C---------------------------------------- scan all sites JQ done (JQ<IQ)
                  DQIPJMAX = 0D0
                  DO JQ = 1,IQ - 1
                     QIJVEC(1:3) = QIPVEC(1:3) - QBAS(1:3,JQ)
                     DQIJ = DNRM2(3,QIJVEC,1)
                     DQIPJMAX = MAX(DQIJ,DQIPJMAX)
                  END DO
C
                  IF ( DQIPJMAX.LT.(DQIJMAX_Q-TOL) ) THEN
                     DQIJMAX_Q = DQIPJMAX
                     I1MIN = I1
                     I2MIN = I2
                     I3MIN = I3
                  END IF
C
               END DO
            END DO
         END DO
C------------------------------------------- scan neighboring unit cells
C
         IF ( I1MIN.NE.0 .OR. I2MIN.NE.0 .OR. I3MIN.NE.0 ) THEN
C
            QBAS(1:3,IQ) = QBAS(1:3,IQ) + I1MIN*ABAS(1:3,1)
     &                     + I2MIN*ABAS(1:3,2) + I3MIN*ABAS(1:3,3)
C
            CHANGES_MADE = .TRUE.
         END IF
C
         DQIJMAX_NEW = MAX(DQIJMAX_Q,DQIJMAX_NEW)
C
      END DO
C----------------------------------------------------- scan all sites IQ
C
      IF ( CHANGES_MADE ) THEN
C
         WRITE (6,99001)
         DO IQ = 1,NQ
            WRITE (6,99002) (QBAS_OLD(I,IQ),I=1,3),(QBAS(I,IQ),I=1,3)
         END DO
         WRITE (6,99003) DQIJMAX_OLD,DQIJMAX_NEW
         IFLAG = 1
C
         IF ( DQIJMAX_NEW.GT.(DQIJMAX_OLD+TOL) ) THEN
            WRITE (6,99005)
            IFLAG = -1
         END IF
C
         CALL OPTIMIZE_BASIS_RASMOL(NQ,ABAS,QBAS_OLD,QBAS,IOTMP,NQMAX)
C
      ELSE
C
         WRITE (6,99004) DQIJMAX_NEW
C
      END IF
C
99001 FORMAT (//,2(1X,79('*'),/),30X,'<OPTIMIZE_BASIS_3D>',/,
     &        2(1X,79('*'),/),//,15X,'OLD basis vectors',21X,
     &        'NEW basis vectors '/)
99002 FORMAT (2X,2(3X,'(',F10.5,',',F10.5,',',F10.5,' )'))
99003 FORMAT (/,10X,'largest distance between basis sites',/,10X,
     &        'originally:      ',F10.5,' a.u.',/,10X,
     &        'optimized:       ',F10.5,' a.u.',/)
99004 FORMAT (//,2(1X,79('*'),/),30X,'<OPTIMIZE_BASIS_3D>',/,
     &        2(1X,79('*'),/),//,10X,'no changes done on basis vectors',
     &        /,10X,'largest distance between basis sites:',F10.5,
     &        ' a.u.',/)
99005 FORMAT (//,2(1X,79('#'),/),10X,'TROUBLE finding optimized basis',
     &        /,2(1X,79('#'),/),//)
      END
C*==optimize_basis_rasmol.f    processed by SPAG 6.70Rc at 13:55 on 27 Jan 2017
      SUBROUTINE OPTIMIZE_BASIS_RASMOL(NQ,ABAS,QBAS_OLD,QBAS,IOTMP,
     &                                 NQMAX)
C   ********************************************************************
C   *                                                                  *
C   *    create data and script files for  rasmol                      *
C   *                                                                  *
C   *    - unit cell and polyhedra                                     *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_CONSTANTS,ONLY:A0_ANG
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='OPTIMIZE_BASIS_RASMOL')
C
C Dummy arguments
C
      INTEGER IOTMP,NQ,NQMAX
      REAL*8 ABAS(3,3),QBAS(3,NQMAX),QBAS_OLD(3,NQMAX)
C
C Local variables
C
      REAL*8 COLOR,S
      CHARACTER*80 FILPDB,FILRAS
      INTEGER I,IQ,IX,LFIL,NPOINTS_SITES,NPOINTS_TOT,NPOINTS_UC
C
C*** End of declarations rewritten by SPAG
C
C----- use larger scaling factor to avoid spurious lines drawn by rasmol
      S = 3*A0_ANG*8D0
      S = 6*S
C
C=======================================================================
C
      FILPDB = 'OPTIMIZE_BASIS'
      LFIL = 14
C
      FILPDB = FILPDB(1:LFIL)//'.pdb'
      LFIL = LFIL + 4
C
      CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILPDB(1:LFIL))
      WRITE (IOTMP,99003) 'unit cell with OLD and NEW site positions'
C
C------------------------------------------ specify corners of unit cell
C
      COLOR = 3D0
      WRITE (IOTMP,FMT=99002) 1,1,0D0,0D0,0D0,COLOR
      DO I = 1,3
         WRITE (IOTMP,FMT=99002) (I+1),(I+1),(ABAS(IX,I)*S,IX=1,3),COLOR
      END DO
C
      WRITE (IOTMP,FMT=99002) 5,5,S*(ABAS(1,1)+ABAS(1,2)),
     &                        S*(ABAS(2,1)+ABAS(2,2)),
     &                        S*(ABAS(3,1)+ABAS(3,2)),COLOR
      WRITE (IOTMP,FMT=99002) 6,6,S*(ABAS(1,1)+ABAS(1,3)),
     &                        S*(ABAS(2,1)+ABAS(2,3)),
     &                        S*(ABAS(3,1)+ABAS(3,3)),COLOR
      WRITE (IOTMP,FMT=99002) 7,7,S*(ABAS(1,3)+ABAS(1,2)),
     &                        S*(ABAS(2,3)+ABAS(2,2)),
     &                        S*(ABAS(3,3)+ABAS(3,2)),COLOR
      WRITE (IOTMP,FMT=99002) 8,8,S*(ABAS(1,1)+ABAS(1,2)+ABAS(1,3)),
     &                        S*(ABAS(2,1)+ABAS(2,2)+ABAS(2,3)),
     &                        S*(ABAS(3,1)+ABAS(3,2)+ABAS(3,3)),COLOR
C
      NPOINTS_UC = 8
      NPOINTS_SITES = 0
C
C------------------------------------------ specify OLD atomic positions
C
      COLOR = 1.0D0
      I = NPOINTS_UC + NPOINTS_SITES
      DO IQ = 1,NQ
         I = I + 1
         WRITE (IOTMP,99001) I,I,(QBAS_OLD(IX,IQ)*S,IX=1,3),COLOR
      END DO
      NPOINTS_SITES = NPOINTS_SITES + NQ
C
C------------------------------------------ specify NEW atomic positions
C
      COLOR = 2.0D0
      DO IQ = 1,NQ
         I = I + 1
         WRITE (IOTMP,99001) I,I,(QBAS(IX,IQ)*S,IX=1,3),COLOR
      END DO
      NPOINTS_SITES = NPOINTS_SITES + NQ
C
      NPOINTS_TOT = NPOINTS_UC + NPOINTS_SITES
C
      WRITE (IOTMP,FMT=99005)
      CLOSE (IOTMP)
C
C--------------------------------------------------- write RASMOL script
C
      FILRAS = FILPDB(1:(LFIL-4))//'.ras'
C
      CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILRAS(1:LFIL))
C
      WRITE (IOTMP,*) 'load '''//FILPDB(1:LFIL)//'''  '
      WRITE (IOTMP,*) 'set background white'
      WRITE (IOTMP,*) 'color temperature'
      WRITE (IOTMP,*) 'set axes on '
C
      WRITE (IOTMP,99006) 1,NPOINTS_UC
      WRITE (IOTMP,*) 'cpk  0'
      WRITE (IOTMP,99006) (NPOINTS_UC+1),(NPOINTS_UC+NPOINTS_SITES)
      WRITE (IOTMP,*) 'cpk 150'
      WRITE (IOTMP,99006) (NPOINTS_UC+NPOINTS_SITES+1),NPOINTS_TOT
      WRITE (IOTMP,*) 'cpk  0'
C
      CLOSE (IOTMP)
C
      WRITE (6,99004) FILPDB(1:LFIL),FILRAS(1:LFIL)
C
C=======================================================================
C
99001 FORMAT ('ATOM  ',I5,'          ',I5,'    ',3F8.3,'  0.00',F8.3)
99002 FORMAT ('HETATM',I5,'          ',I5,'    ',3F8.3,'  0.00',F8.3)
99003 FORMAT ('HEADER    ',A,/,'SOURCE    SPRKKR - program       ',/,
     &        'AUTHOR    H. Ebert               ',/,
     &        'REMARK    None                   ')
99004 FORMAT (/,5X,'unit cell data stored in rasmol data-file ',A,/,5X,
     &        'view via:   rasmol  -script ',A)
99005 FORMAT ('CONECT    1    2',/,'CONECT    1    3',/,
     &        'CONECT    1    4',/,'CONECT    5    2',/,
     &        'CONECT    5    3',/,'CONECT    5    8',/,
     &        'CONECT    6    2',/,'CONECT    6    4',/,
     &        'CONECT    6    8',/,'CONECT    7    3',/,
     &        'CONECT    7    4',/,'CONECT    7    8',/,'END')
99006 FORMAT ('SELECT ',I5,'-',I5)
      END
