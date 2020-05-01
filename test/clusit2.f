C*==clusit2.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CLUSIT2(IOTMP,NQCLU,RQCLU,IPRINT,NQCLUMAX)
C   ********************************************************************
C   *                                                                  *
C   * find:                                                            *
C   *     - equivalent  S-S'-blocks of G-matrix                        *
C   *     - S-S'-blocks related by inversion symmetry                  *
C   *     - different length of S-S' vectors                           *
C   *     - different directions of S-S' vectors                       *
C   *                                                                  *
C   * ---------------------------------------------------------------- *
C   *                                                                  *
C   * the subroutine determines the additional cluster parameters:     *
C   *                                                                  *
C   *   NSSP1,NIJSTABMAX,NSSP2A,NSSP2B,NSSPABS,NNSP4MAX,NSSPDIR        *
C   *   RSSP,DSSP,NIJSTAB,IQCLUTAB,JQCLUTAB                            *
C   *   ISSP2AB,ISSP2BB,NSSP4,ISSP4,ISDA4,JSDA4,ISSPDIR,ISSP5          *
C   *                                                                  *
C   * - a rough estimate is made first for the corresponding           *
C   *   array sizes and the storage is allocated                       *
C   * - ALL data are passed via the temporary file  IOTMP              *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--CLUSIT226
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CLUSIT2')
      REAL*8 TOL
      PARAMETER (TOL=1D-5)
      LOGICAL CHECK
      PARAMETER (CHECK=.FALSE.)
C
C Dummy arguments
C
      INTEGER IOTMP,IPRINT,NQCLU,NQCLUMAX
      REAL*8 RQCLU(3,NQCLUMAX)
C
C Local variables
C
      REAL*8 DIR(3),DSSP(:),RSSP(:,:),SD
      CHARACTER*1 GC(:,:)
      COMPLEX*16 GMAT(:,:)
      INTEGER I,I01,I02,I0A,I0B,I1,I2,I2A,I2B,I5,IA,IA_ERR,ID,IQCLU,
     &        IQCLUTAB(:,:),ISDA4(:,:),ISSP,ISSP2A(:),ISSP2AB(:),
     &        ISSP2BB(:),ISSP2C(:),ISSP4(:,:),ISSP5(:,:),ISSPDIR(:),J,
     &        J01,J02,J0A,J0B,JA,JABS,JQCLU,JQCLUTAB(:,:),JSDA4(:,:),
     &        JSSP,M,N,NIJSTAB(:),NIJSTABMAX,NN0,NN1,NN2,NN4,NN5,
     &        NNSP4MAX,NSSP1,NSSP2A,NSSP2B,NSSP4(:),NSSPABS,NSSPDIR
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE GMAT,GC
      ALLOCATABLE DSSP,RSSP,NIJSTAB,IQCLUTAB,JQCLUTAB
      ALLOCATABLE ISSP2A,ISSP2AB,ISSP2BB,ISSP2C
      ALLOCATABLE NSSP4,ISSP4,ISDA4,JSDA4
      ALLOCATABLE ISSPDIR,ISSP5
C
      ALLOCATE (GMAT(NQCLU,NQCLU),GC(NQCLU,NQCLU),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: GC')
C
C=======================================================================
C                  CHECK - CALCULATE ALL BLOCKS INDIVIDUALLY
C=======================================================================
      IF ( CHECK ) THEN
C
         NN0 = NQCLU
         NN1 = NQCLU*NQCLU
C
         ALLOCATE (DSSP(0:NN1),RSSP(3,0:NN1),NIJSTAB(NN1),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: DSSP')
         ALLOCATE (IQCLUTAB(NN0,NN1),JQCLUTAB(NN0,NN1),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: IQCLUTAB')
C
         NIJSTAB(:) = 0
         DSSP(:) = 0D0
         RSSP(:,:) = 0D0
C
         IQCLUTAB(:,:) = 9999
         JQCLUTAB(:,:) = 9999
C
C     NN2 = NN1/2 + 1
         NN2 = NQCLU*NQCLU
         ALLOCATE (ISSP2A(NN2),ISSP2AB(NN2))
         ALLOCATE (ISSP2BB(NN2),ISSP2C(NN1),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ISSP2BB')
         ISSP2A(:) = 9999
         ISSP2AB(:) = 9999
         ISSP2BB(:) = 9999
         ISSP2C(:) = 9999
C
         NN4 = NQCLU*NQCLU
         ALLOCATE (NSSP4(NN2),ISSP4(NN4,NN2))
         ALLOCATE (ISDA4(NN4,NN2),JSDA4(NN4,NN2),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ISDA4')
         NSSP4(:) = 9999
         ISSP4(:,:) = 9999
         ISDA4(:,:) = 9999
         JSDA4(:,:) = 9999
C
         NN5 = NN1
         ALLOCATE (ISSPDIR(NN5),ISSP5(NN5,NN4),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ISSPDIR')
         ISSPDIR(:) = 9999
         ISSP5(:,:) = 9999
C
         NSSP1 = 0
         DO IQCLU = 1,NQCLU
            DO JQCLU = 1,NQCLU
               DO ISSP = 1,NSSP1
                  SD = 0.0D0
                  DO I = 1,3
                     SD = SD + 
     &                    ABS(RSSP(I,ISSP)-(RQCLU(I,JQCLU)-RQCLU(I,IQCLU
     &                    )))
                  END DO
                  IF ( SD.LT.TOL ) THEN
                     NIJSTAB(ISSP) = NIJSTAB(ISSP) + 1
                     IQCLUTAB(NIJSTAB(ISSP),ISSP) = IQCLU
                     JQCLUTAB(NIJSTAB(ISSP),ISSP) = JQCLU
                     GOTO 20
                  END IF
               END DO
C
               IF ( (NSSP1+1).GT.NN1 ) THEN
                  WRITE (6,99001) 'NSSP1','NN1',NN1
                  CALL STOP_MESSAGE(ROUTINE,'NSSP1 NN1 ???')
               END IF
               NSSP1 = NSSP1 + 1
               NIJSTAB(NSSP1) = 1
               IQCLUTAB(1,NSSP1) = IQCLU
               JQCLUTAB(1,NSSP1) = JQCLU
               SD = 0.0D0
               DO I = 1,3
                  RSSP(I,NSSP1) = RQCLU(I,JQCLU) - RQCLU(I,IQCLU)
                  SD = SD + RSSP(I,NSSP1)**2
               END DO
               DSSP(NSSP1) = SQRT(SD)
C
 20         END DO
         END DO
C
C left from inversion symmetry
C
         NSSP2A = 0
         NSSP2B = 0
         DO ISSP = 1,NSSP1
            NSSP2A = NSSP2A + 1
            ISSP2A(NSSP2A) = ISSP
            ISSP2C(ISSP) = ISSP
         END DO
C
C left from |RSSP| check
C
         NSSPABS = 0
         DO I2A = 1,NSSP2A
            NSSPABS = NSSPABS + 1
            NSSP4(NSSPABS) = 1
            ISSP4(1,NSSPABS) = ISSP2A(I2A)
            ISDA4(1,NSSPABS) = IQCLUTAB(1,ISSP2A(I2A))
            JSDA4(1,NSSPABS) = JQCLUTAB(1,ISSP2A(I2A))
         END DO
C
C left from direction check
C
         NSSPDIR = 0
         DO IA = 1,NSSPABS
            DO ID = 1,NSSP4(IA)
               ISSP = ISSP4(ID,IA)
               IF ( DABS(DSSP(ISSP)).GE.1D-6 ) THEN
                  DO I = 1,3
                     DIR(I) = RSSP(I,ISSP)/DSSP(ISSP)
                  END DO
               ELSE
                  DO I = 1,3
                     DIR(I) = 0D0
                  END DO
               END IF
               NSSPDIR = NSSPDIR + 1
               ISSPDIR(NSSPDIR) = ISSP
               ISSP5(ID,IA) = NSSPDIR
            END DO
         END DO
C
C------------------------------------------------------- go to PRINT OUT
         GOTO 300
C
      END IF
C============================= CHECK - CALCULATE ALL BLOCKS INDIVIDUALLY
C
C=======================================================================
C                   find table of equivalent S-S' blocks
C=======================================================================
C
C11111111111111111111111111111111111111111111111111111111111111111111111
C  STEP 1:            set up table of pointers
C11111111111111111111111111111111111111111111111111111111111111111111111
C
      NN0 = NQCLU
      NN1 = NQCLU*NQCLU
C
      ALLOCATE (DSSP(0:NN1),RSSP(3,0:NN1),NIJSTAB(NN1),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: DSSP')
      ALLOCATE (IQCLUTAB(NN0,NN1),JQCLUTAB(NN0,NN1),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: IQCLUTAB')
C
      NIJSTAB(:) = 0
      DSSP(:) = 0D0
      RSSP(:,:) = 0D0
C
      IQCLUTAB(:,:) = 9999
      JQCLUTAB(:,:) = 9999
C
      NSSP1 = 0
      DO IQCLU = 1,NQCLU
         DO JQCLU = 1,NQCLU
            DO ISSP = 1,NSSP1
               SD = 0.0D0
               DO I = 1,3
                  SD = SD + 
     &                 ABS(RSSP(I,ISSP)-(RQCLU(I,JQCLU)-RQCLU(I,IQCLU)))
               END DO
               IF ( SD.LT.TOL ) THEN
                  NIJSTAB(ISSP) = NIJSTAB(ISSP) + 1
                  IQCLUTAB(NIJSTAB(ISSP),ISSP) = IQCLU
                  JQCLUTAB(NIJSTAB(ISSP),ISSP) = JQCLU
                  GOTO 50
               END IF
            END DO
C
            IF ( (NSSP1+1).GT.NN1 ) THEN
               WRITE (6,99001) 'NSSP1','NN1',NN1
               CALL STOP_MESSAGE(ROUTINE,'NSSP1 NN1 ???')
            END IF
            NSSP1 = NSSP1 + 1
            NIJSTAB(NSSP1) = 1
            IQCLUTAB(1,NSSP1) = IQCLU
            JQCLUTAB(1,NSSP1) = JQCLU
            SD = 0.0D0
            DO I = 1,3
               RSSP(I,NSSP1) = RQCLU(I,JQCLU) - RQCLU(I,IQCLU)
               SD = SD + RSSP(I,NSSP1)**2
            END DO
            DSSP(NSSP1) = SQRT(SD)
C
 50      END DO
      END DO
C
C
C22222222222222222222222222222222222222222222222222222222222222222222222
C  STEP 2:   check for inversion symmetry
C            table A:  representative set
C            table B:  related by pointer  ISSP2AB  to table A
C22222222222222222222222222222222222222222222222222222222222222222222222
C
      NN2 = (NSSP1+1)/2 + 1
C
      ALLOCATE (ISSP2A(NN2),ISSP2AB(NN2))
      ALLOCATE (ISSP2BB(NN2),ISSP2C(NSSP1),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ISSP2BB')
      ISSP2A(:) = 9999
      ISSP2AB(:) = 9999
      ISSP2BB(:) = 9999
      ISSP2C(:) = 9999
C
      NSSP2A = 0
      NSSP2B = 0
      DO ISSP = 1,NSSP1
         DO JA = 1,NSSP2A
            JSSP = ISSP2A(JA)
            SD = 0.0D0
            DO I = 1,3
               SD = SD + ABS(RSSP(I,ISSP)+RSSP(I,JSSP))
            END DO
            IF ( SD.LT.TOL ) THEN
               NSSP2B = NSSP2B + 1
               IF ( RSSP(1,JSSP).GT.RSSP(1,ISSP) ) THEN
                  ISSP2AB(NSSP2B) = JSSP
                  ISSP2BB(NSSP2B) = ISSP
                  ISSP2C(ISSP) = JSSP
               ELSE
                  ISSP2A(JA) = ISSP
                  ISSP2C(JSSP) = ISSP
                  ISSP2AB(NSSP2B) = ISSP
                  ISSP2BB(NSSP2B) = JSSP
                  ISSP2C(NSSP2B) = ISSP
               END IF
               GOTO 100
            END IF
         END DO
         NSSP2A = NSSP2A + 1
         IF ( NSSP2A.GT.NN2 ) THEN
            WRITE (6,99001) 'NSSP2A','NN2',NN2
            CALL STOP_MESSAGE(ROUTINE,'NSSP2A NN2 ???')
         END IF
         ISSP2A(NSSP2A) = ISSP
         ISSP2C(ISSP) = ISSP
 100  END DO
C
C
C33333333333333333333333333333333333333333333333333333333333333333333333
C  STEP 3:   check for subsets of S-S' blocks
C            with same |RSSP| within set 2A
C33333333333333333333333333333333333333333333333333333333333333333333333
C
      NN4 = NN2
C
      ALLOCATE (NSSP4(NN2),ISSP4(NN4,NN2))
      ALLOCATE (ISDA4(NN4,NN2),JSDA4(NN4,NN2),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ISDA4')
      NSSP4(:) = 9999
      ISSP4(:,:) = 9999
      ISDA4(:,:) = 9999
      JSDA4(:,:) = 9999
C
      NSSPABS = 0
      DO I2A = 1,NSSP2A
         SD = DSSP(ISSP2A(I2A))
         DO JABS = 1,NSSPABS
            IF ( ABS(SD-DSSP(ISSP4(1,JABS))).LT.TOL ) THEN
               IF ( JABS.GT.NN4 ) THEN
                  WRITE (6,99001) 'NSSP4','NN4',NN4
                  CALL STOP_MESSAGE(ROUTINE,'NSSP4 > NN4 ???')
               END IF
               NSSP4(JABS) = NSSP4(JABS) + 1
               ISSP4(NSSP4(JABS),JABS) = ISSP2A(I2A)
               ISDA4(NSSP4(JABS),JABS) = IQCLUTAB(1,ISSP2A(I2A))
               JSDA4(NSSP4(JABS),JABS) = JQCLUTAB(1,ISSP2A(I2A))
               GOTO 200
            END IF
         END DO
         NSSPABS = NSSPABS + 1
         IF ( NSSPABS.GT.NN2 ) THEN
            WRITE (6,99001) 'NSSPABS','NN2',NN2
            CALL STOP_MESSAGE(ROUTINE,'NSSPABS > NN2 ???')
         END IF
C
         NSSP4(NSSPABS) = 1
         ISSP4(1,NSSPABS) = ISSP2A(I2A)
         ISDA4(1,NSSPABS) = IQCLUTAB(1,ISSP2A(I2A))
         JSDA4(1,NSSPABS) = JQCLUTAB(1,ISSP2A(I2A))
 200  END DO
C
C
C44444444444444444444444444444444444444444444444444444444444444444444444
C  STEP 4:   check for the orientation of the vector     RSSP
C44444444444444444444444444444444444444444444444444444444444444444444444
C
      NN5 = NN2
C
      ALLOCATE (ISSPDIR(NN5),ISSP5(NN5,NN4),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ISSPDIR')
      ISSPDIR(:) = 9999
      ISSP5(:,:) = 9999
C
      NSSPDIR = 0
      DO IA = 1,NSSPABS
         DO ID = 1,NSSP4(IA)
            ISSP = ISSP4(ID,IA)
            IF ( DABS(DSSP(ISSP)).GE.1D-6 ) THEN
               DO I = 1,3
                  DIR(I) = RSSP(I,ISSP)/DSSP(ISSP)
               END DO
            ELSE
               DO I = 1,3
                  DIR(I) = 0D0
               END DO
            END IF
C
            DO I5 = 1,NSSPDIR
               JSSP = ISSPDIR(I5)
               SD = 0.0D0
               DO I = 1,3
                  IF ( DABS(DSSP(JSSP)).GE.1D-6 ) THEN
                     SD = SD + DABS(DIR(I)-RSSP(I,JSSP)/DSSP(JSSP))
                  ELSE
                     SD = SD + DABS(DIR(I))
                  END IF
               END DO
               IF ( SD.LT.TOL ) THEN
                  ISSP5(ID,IA) = I5
                  GOTO 250
               END IF
            END DO
            NSSPDIR = NSSPDIR + 1
C
            IF ( ID.GT.NN5 .OR. NSSPDIR.GT.NN5 ) THEN
               WRITE (6,99001) 'ID','NN5',NN5
               CALL STOP_MESSAGE(ROUTINE,'ID NN5 ???')
            END IF
            ISSPDIR(NSSPDIR) = ISSP
            ISSP5(ID,IA) = NSSPDIR
 250     END DO
      END DO
C
C
C=======================================================================
C           print results to standard output according to request
C=======================================================================
C
 300  CONTINUE
      IF ( IPRINT.GE.2 ) THEN
         WRITE (6,99002) NSSP1
         WRITE (6,99003)
         DO N = 1,NSSP1
            WRITE (6,99004) N,ISSP2C(N),NIJSTAB(N),(RSSP(I,N),I=1,3),
     &                      (IQCLUTAB(M,N),JQCLUTAB(M,N),M=1,NIJSTAB(N))
         END DO
      END IF
C
      IF ( IPRINT.GE.1 ) THEN
         WRITE (6,99005) NSSP2A,NSSP2B
         WRITE (6,99006) NSSPABS,NSSPDIR
      END IF
C
      IF ( IPRINT.GE.2 ) THEN
         WRITE (6,99007)
C
         DO IA = 1,NSSPABS
            WRITE (6,99008) DSSP(ISSP4(1,IA)),NSSP4(IA),
     &                      (ISSP4(ID,IA),ISSP5(ID,IA),ISDA4(ID,IA),
     &                      JSDA4(ID,IA),ID=1,NSSP4(IA))
         END DO
C---------------------------------------- check structure in <CLUTAUMAT>
C
         DO I = 1,NQCLU
            DO J = 1,NQCLU
               GC(I,J) = '-'
               GMAT(I,J) = C0
            END DO
         END DO
C
C------------------------------ set up all inequivalent G(S,S') - blocks
         DO IA = 2,NSSPABS
            DO ID = 1,NSSP4(IA)
               I = ISDA4(ID,IA)
               J = JSDA4(ID,IA)
               GMAT(I,J) = DCMPLX(DBLE(ID),DBLE(IA))
               GC(I,J) = 'A'
            END DO
         END DO
C
C------------------------------------------------ use inversion symmetry
         DO I2B = 1,NSSP2B
C
            I0A = IQCLUTAB(1,ISSP2AB(I2B))
            J0A = JQCLUTAB(1,ISSP2AB(I2B))
            I0B = IQCLUTAB(1,ISSP2BB(I2B))
            J0B = JQCLUTAB(1,ISSP2BB(I2B))
C
            GMAT(I0B,J0B) = -GMAT(I0A,J0A)
            GC(I0B,J0B) = 'B'
         END DO
C
C------------------------------------------------- copy identical blocks
         DO I1 = 2,NSSP1
            I01 = IQCLUTAB(1,I1)
            J01 = JQCLUTAB(1,I1)
C
            DO I2 = 2,NIJSTAB(I1)
               I02 = IQCLUTAB(I2,I1)
               J02 = JQCLUTAB(I2,I1)
C
               GMAT(I02,J02) = GMAT(I01,J01)
               GC(I02,J02) = 'C'
            END DO
C
         END DO
C
         WRITE (6,99009)
         DO I = 1,NQCLU
            WRITE (6,99010) I,(GC(I,J),J=1,NQCLU)
         END DO
C
         N = NQCLU
         CALL CMATSTRUCT('(ID,IA)-indices for G-matrix',GMAT,N,N,0,0,1,
     &                   1D-6,6)
C
      END IF
C
C=======================================================================
C  write results to temporary file for data transfer to calling routine
C=======================================================================
C
      NIJSTABMAX = 0
      DO I = 1,NSSP1
         NIJSTABMAX = MAX(NIJSTABMAX,NIJSTAB(I))
      END DO
      IF ( NIJSTABMAX.GT.NN0 ) THEN
         WRITE (6,99001) 'NIJSTABMAX','NN0',NN0
         CALL STOP_MESSAGE(ROUTINE,'NIJSTABMAX NN0 ???')
      END IF
C
      IF ( NSSP2B.GT.NSSP2A ) THEN
         WRITE (6,99001) 'NSSP2B','NSSP2A',NSSP2A
         CALL STOP_MESSAGE(ROUTINE,'NSSP2B NSSP2A ???')
      END IF
C
      NNSP4MAX = 0
      DO I = 1,NSSPABS
         NNSP4MAX = MAX(NNSP4MAX,NSSP4(I))
      END DO
C
      NN5 = NSSPDIR
C
      WRITE (IOTMP) NSSP1,NIJSTABMAX,NSSP2A,NSSP2B,NSSPABS,NNSP4MAX,
     &              NSSPDIR
      WRITE (IOTMP) ((RSSP(J,I),J=1,3),DSSP(I),NIJSTAB(I),I=1,NSSP1),
     &              ((IQCLUTAB(J,I),JQCLUTAB(J,I),J=1,NIJSTAB(I)),I=1,
     &              NSSP1),(ISSP2AB(I),ISSP2BB(I),I=1,NSSP2A),
     &              (NSSP4(I),I=1,NSSPABS),
     &              ((ISSP4(J,I),ISDA4(J,I),JSDA4(J,I),J=1,NSSP4(I)),
     &              I=1,NSSPABS),(ISSPDIR(I),I=1,NSSPDIR),
     &              ((ISSP5(J,I),J=1,NSSP4(I)),I=1,NSSPABS)
C
C ----------------------------------------------------------------------
C
      DEALLOCATE (GMAT,GC,DSSP,RSSP,NIJSTAB,IQCLUTAB,JQCLUTAB)
      DEALLOCATE (ISSP2A,ISSP2AB,ISSP2BB,ISSP2C)
      DEALLOCATE (NSSP4,ISSP4,ISDA4,JSDA4,ISSPDIR,ISSP5)
C
99001 FORMAT (/,10X,'STOP IN <CLUSIT2>   ',/,10X,A,'  >  ',A,' =',I5)
99002 FORMAT (/,10X,'number of inequivalent SS''-block-matrices',
     &        '  NSSP1 =',I3,/)
99003 FORMAT (' ISSP ISSP2 NIJS (  ->R[IS] - ->R[JS]  ) [IS,JS] ...')
99004 FORMAT (3I5,2X,'(',F7.3,',',F7.3,',',F7.3,' )',1X,
     &        5(:,'[',I2,',',I2,']'),:,/,(44X,5(:,'[',I2,',',I2,']')))
99005 FORMAT (/,10X,'number of SS''-blocks in subset A',I12,/,10X,
     &        'number of SS''-blocks in subset B',I12,/,10X,
     &        'related to subset A by inversion')
99006 FORMAT (/,10X,'number of blocks with different |RSSP|  NSSPABS =',
     &        I4,/,10X,'with NSSPDIR = ',I5,'  different orientations ')
99007 FORMAT (10X,'        (indicated in parenthesis)',//,10X,
     &        ' |RSSP|  NSSP4  blocks with same |RSSP| ')
99008 FORMAT (10X,F7.3,I6,2X,3(:,I5,'(',I3,')[',I2,',',I2,']'),:,/,
     &        (25X,3(:,I5,'(',I3,')[',I2,',',I2,']')))
99009 FORMAT (//,10X,'setting up of G-matrix',/)
99010 FORMAT (10X,I5,5X,200(1X,A1))
      END
