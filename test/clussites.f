C*==clussites.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE CLUSSITES(IOTMP,IPRINT,MOL,SYSTEM_DIMENSION,ABAS,
     &                     ABAS_L,ABAS_I,ABAS_R,ADAINV_L,ADAINV_I,
     &                     ADAINV_R,QBAS,CLURAD,IQCNTR,NQCLU,NQCLU_L,
     &                     NQCLU_I,NQCLU_R,NSHLCLU,NQ,NQ_L,NQ_R,NQMAX)
C   ********************************************************************
C   *                                                                  *
C   * find the configuration of a cluster around a central site IQCNTR *
C   * either specified by the cluster of radius  CLURAD                *
C   * or the number of atomic shells             NSHLCLU               *
C   *                                                                  *
C   * if MOL = .TRUE.:   do not scan neighboring units cells           *
C   *                    but use   ALL    NQ basis atoms that          *
C   *                    make the cluster                              *
C   *                                                                  *
C   * ---------------------------------------------------------------- *
C   *                                                                  *
C   * the subroutine determines the cluster parameters:                *
C   *                                                                  *
C   * RQCLU,DQCLU,IQ_IQCLU,NQSHLCLU                                    *
C   *                                                                  *
C   * - a rough estimate is made first for the corresponding           *
C   *   array sizes   NQCLUMAX  and  NSHLCLUMAX                        *
C   *   using these the storage is allocated                           *
C   * - the actual array sizes   NQCLU  and  NSHLCLU  are passed       *
C   *   to the calling subroutine via the argument list                *
C   * - the data are passed via the temporary file  IOTMP              *
C   *                                                                  *
C   *  04/03/08                                                        *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CLUSSITES')
C
C Dummy arguments
C
      REAL*8 CLURAD
      INTEGER IOTMP,IPRINT,IQCNTR,NQ,NQCLU,NQCLU_I,NQCLU_L,NQCLU_R,
     &        NQMAX,NQ_L,NQ_R,NSHLCLU
      LOGICAL MOL
      CHARACTER*10 SYSTEM_DIMENSION
      REAL*8 ABAS(3,3),ABAS_I(3,3),ABAS_L(3,3),ABAS_R(3,3),ADAINV_I(3,3)
     &       ,ADAINV_L(3,3),ADAINV_R(3,3),QBAS(3,NQMAX)
C
C Local variables
C
      REAL*8 BON(3,3),BR(3),CEXT,DC(3),DQ(3),DQCLU(:),DR,DSH(:),DTMP,
     &       D_C,D_I,D_L,D_R,RCOF(3),RCOR(3),RQCLU(:,:),RVEC_TMP(3),TOL,
     &       XTMP,YTMP,ZTMP
      REAL*8 DDOT,DNRM2
      INTEGER I,I1,I2,I3,IA_ERR,IBV3,IEXT,IFLAG,IQ,IQ1,IQ2,IQTMP,
     &        IQ_IQCLU(:),ISH,IX,IY,IZ,J,JJ,JQ,K,N5VEC_QCLU(:,:),
     &        N5VEC_TMP(5),NBR(3),NMAX(3),NMAX_I(3),NMAX_L(3),NMAX_R(3),
     &        NMIN(3),NMIN_I(3),NMIN_L(3),NMIN_R(3),NQCLUMAX,NQSHLCLU(:)
     &        ,NSH,NSHLCLUMAX
C
C*** End of declarations rewritten by SPAG
C
      DATA IA_ERR/0/
C
      ALLOCATABLE DSH,RQCLU,DQCLU,IQ_IQCLU,NQSHLCLU,N5VEC_QCLU
C
      IF ( IPRINT.GE.0 ) WRITE (6,99005)
C
      NQCLU_L = 999999
      NQCLU_I = 999999
      NQCLU_R = 999999
C
      IQ1 = 1
      IQ2 = NQ
C
      DO I = 1,3
         BR(I) = DNRM2(3,ABAS(1,I),1)
         NBR(I) = INT(CLURAD/BR(I)) + 2
      END DO
C
      IF ( MOL ) THEN
         NSHLCLU = 0
         CLURAD = 1D6
      END IF
C
C ----------------------------------------------------------------------
C       if NSHLCLU is given: find first corresponding radius CLURAD
C ----------------------------------------------------------------------
C
      TOL = 1.0D-4
      IF ( NSHLCLU.NE.0 ) THEN
C
         ALLOCATE (DSH(2*NSHLCLU),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: DSH')
C
         CLURAD = 3D0
 50      CONTINUE
         CLURAD = CLURAD*1.05D0
C
         DO I = 1,3
            NBR(I) = INT(CLURAD/BR(I)) + 2
         END DO
C
C------------------------------- first shell coincides with central atom
C
         NSH = 1
         DSH(1) = 0D0
C
         DO JQ = 1,NQ
C
            DO I = 1,3
               DQ(I) = QBAS(I,JQ) - QBAS(I,IQCNTR)
            END DO
            DO I1 = -NBR(1), + NBR(1)
               DO I2 = -NBR(2), + NBR(2)
                  DO I3 = -NBR(3), + NBR(3)
C
                     DO I = 1,3
                        DC(I) = I1*ABAS(I,1) + I2*ABAS(I,2)
     &                          + I3*ABAS(I,3) + DQ(I)
                     END DO
                     DR = DNRM2(3,DC,1)
C
                     IF ( DR.GT.TOL ) THEN
                        DO I = 1,NSH
                           IF ( ABS(DR-DSH(I)).LT.TOL ) GOTO 55
                           IF ( DR.LT.(DSH(I)-TOL) ) THEN
                              DO J = MIN(NSH,NSHLCLU-1),I, - 1
                                 DSH(J+1) = DSH(J)
                              END DO
                              DSH(I) = DR
                              IF ( NSH.LT.NSHLCLU ) NSH = NSH + 1
                              GOTO 55
                           END IF
                        END DO
                        IF ( NSH.LT.NSHLCLU ) THEN
                           NSH = NSH + 1
                           DSH(NSH) = DR
                        END IF
                     END IF
C
 55               END DO
               END DO
            END DO
C
         END DO
C
         IF ( NSH.LT.NSHLCLU ) GOTO 50
C
         CLURAD = DSH(NSHLCLU) + 1D-6
C
      END IF
C
C ----------------------------------------------------------------------
C                   allocate temporary arrays
C ----------------------------------------------------------------------
C
      IF ( MOL ) THEN
         NQCLUMAX = NQ
         NSHLCLUMAX = NQ
      ELSE
         NQCLUMAX = NQ
         DO I = 1,3
            NQCLUMAX = NQCLUMAX*(INT(2*CLURAD/BR(I))+1)
         END DO
         NSHLCLUMAX = NQCLUMAX
      END IF
C
      ALLOCATE (RQCLU(3,NQCLUMAX),DQCLU(NQCLUMAX))
      ALLOCATE (N5VEC_QCLU(5,NQCLUMAX))
      ALLOCATE (IQ_IQCLU(NQCLUMAX),NQSHLCLU(NSHLCLUMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: IQCLUS')
C
      IF ( .NOT.ALLOCATED(DSH) ) ALLOCATE (DSH(NSHLCLUMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: DSH')
      GOTO 200
C
 100  CONTINUE
      DEALLOCATE (RQCLU,DQCLU,N5VEC_QCLU,IQ_IQCLU,NQSHLCLU,DSH)
      NQCLUMAX = NINT(NQCLUMAX*1.2)
      NSHLCLUMAX = NINT(NSHLCLUMAX*1.2)
C
      ALLOCATE (RQCLU(3,NQCLUMAX),DQCLU(NQCLUMAX))
      ALLOCATE (N5VEC_QCLU(5,NQCLUMAX))
      ALLOCATE (IQ_IQCLU(NQCLUMAX),NQSHLCLU(NSHLCLUMAX),STAT=IA_ERR)
      ALLOCATE (DSH(NSHLCLUMAX),STAT=IA_ERR)
C
 200  CONTINUE
      IF ( .NOT.MOL ) THEN
         DO I = 1,3
            NBR(I) = INT(CLURAD/BR(I)) + 3
         END DO
      END IF
C
      NQCLU = 1
      N5VEC_QCLU(1:5,1:NQCLUMAX) = 0
      RQCLU(1:3,NQCLU) = 0.0D0
      DQCLU(NQCLU) = 0.0D0
      IQ_IQCLU(NQCLU) = IQCNTR
C
C=======================================================================
      IF ( SYSTEM_DIMENSION(1:2).EQ.'3D' ) THEN
C=======================================================================
C
         DO I = 1,3
            NMIN(I) = -NBR(I)
            NMAX(I) = NBR(I)
         END DO
         IF ( IPRINT.GE.2 ) THEN
            WRITE (6,99003) 'NMIN     = ',NMIN
            WRITE (6,99003) 'NMAX     = ',NMAX
         END IF
C
         IQ1 = 1
         IQ2 = NQ
         IBV3 = 3
         CALL CLUSSITES_ADD(IFLAG,CLURAD,MOL,NMIN,NMAX,ABAS,DQ,QBAS,IQ1,
     &                      IQ2,IQCNTR,RQCLU,DQCLU,N5VEC_QCLU,IBV3,
     &                      IQ_IQCLU,NQCLU,NQCLUMAX,NQMAX)
         IF ( IFLAG.EQ.1 ) GOTO 100
C
      ELSE
C=======================================================================
C          SYSTEM_DIMENSION(1:2) = '2D'
C=======================================================================
C
         DO I = 1,3
            NMIN_L(I) = 0
            NMAX_L(I) = 0
            NMIN_I(I) = 0
            NMAX_I(I) = 0
            NMIN_R(I) = 0
            NMAX_R(I) = 0
         END DO
C
         CALL RVECORTHO(ABAS,BON)
C
         CEXT = CLURAD*1.001D0
C
         DO IX = -1,1,2
            DO IY = -1,1,2
               DO IZ = -1,1,2
                  RCOR(1) = QBAS(1,IQCNTR) + IX*CEXT
                  RCOR(2) = QBAS(2,IQCNTR) + IY*CEXT
                  RCOR(3) = QBAS(3,IQCNTR) + IZ*CEXT
C
                  CALL RVECEXPAND(RCOR,ABAS_L,ADAINV_L,RCOF)
C
                  DO I = 1,3
                     IF ( RCOF(I).GT.0D0 ) THEN
                        IEXT = INT(RCOF(I)) + 2
                     ELSE
                        IEXT = INT(RCOF(I)) - 1
                     END IF
                     NMIN_L(I) = MIN(NMIN_L(I),IEXT)
                     NMAX_L(I) = MAX(NMAX_L(I),IEXT)
                  END DO
C
                  DO I = 1,3
                     RCOR(I) = RCOR(I) - ABAS_L(I,3)
                  END DO
C
                  CALL RVECEXPAND(RCOR,ABAS_I,ADAINV_I,RCOF)
C
                  DO I = 1,3
                     IF ( RCOF(I).GT.0D0 ) THEN
                        IEXT = INT(RCOF(I)) + 2
                     ELSE
                        IEXT = INT(RCOF(I)) - 1
                     END IF
                     NMIN_I(I) = MIN(NMIN_I(I),IEXT)
                     NMAX_I(I) = MAX(NMAX_I(I),IEXT)
                  END DO
C
                  DO I = 1,3
                     RCOR(I) = RCOR(I) - ABAS_I(I,3)
                  END DO
C
                  CALL RVECEXPAND(RCOR,ABAS_R,ADAINV_R,RCOF)
C
                  DO I = 1,3
                     IF ( RCOF(I).GT.0D0 ) THEN
                        IEXT = INT(RCOF(I)) + 2
                     ELSE
                        IEXT = INT(RCOF(I)) - 1
                     END IF
                     NMIN_R(I) = MIN(NMIN_R(I),IEXT)
                     NMAX_R(I) = MAX(NMAX_R(I),IEXT)
                  END DO
               END DO
            END DO
         END DO
C
         D_C = ABS(DDOT(3,QBAS(1,IQCNTR),1,BON(1,3),1))
         D_L = ABS(DDOT(3,ABAS_L(1,3),1,BON(1,3),1))
         D_I = ABS(DDOT(3,ABAS_I(1,3),1,BON(1,3),1))
         D_R = ABS(DDOT(3,ABAS_R(1,3),1,BON(1,3),1))
C
         NMIN_L(3) = -INT((CLURAD-D_C)/D_L) - 1
         NMAX_L(3) = 0
         NMIN_I(3) = 0
         NMAX_I(3) = 0
         NMIN_R(3) = 0
         NMAX_R(3) = INT((CLURAD-D_L-D_I+D_C)/D_R) + 1
C
         IF ( IPRINT.GE.2 ) THEN
            WRITE (6,99004) 'CLURAD   = ',CLURAD
            WRITE (6,99004) 'QCTR     = ',(QBAS(1,IQCNTR),I=1,3),D_C
            WRITE (6,99003) 'NMIN (L) = ',NMIN_L,D_L
            WRITE (6,99003) 'NMAX (L) = ',NMAX_L
            WRITE (6,99003) 'NMIN (I) = ',NMIN_I,D_I
            WRITE (6,99003) 'NMAX (I) = ',NMAX_I
            WRITE (6,99003) 'NMIN (R) = ',NMIN_R,D_R
            WRITE (6,99003) 'NMAX (R) = ',NMAX_R
         END IF
C
C ----------------------------------------------------------------------
C
C................................................................ L BULK
         IQ1 = 1
         IQ2 = NQ_L
         IBV3 = 4
         CALL CLUSSITES_ADD(IFLAG,CLURAD,MOL,NMIN_L,NMAX_L,ABAS_L,DQ,
     &                      QBAS,IQ1,IQ2,IQCNTR,RQCLU,DQCLU,N5VEC_QCLU,
     &                      IBV3,IQ_IQCLU,NQCLU,NQCLUMAX,NQMAX)
         IF ( IFLAG.EQ.1 ) GOTO 100
         NQCLU_L = NQCLU
C
C...................................................... INTERACTION ZONE
         IQ1 = NQ_L + 1
         IQ2 = NQ - NQ_R
         CALL CLUSSITES_ADD(IFLAG,CLURAD,MOL,NMIN_I,NMAX_I,ABAS_I,DQ,
     &                      QBAS,IQ1,IQ2,IQCNTR,RQCLU,DQCLU,N5VEC_QCLU,
     &                      IBV3,IQ_IQCLU,NQCLU,NQCLUMAX,NQMAX)
         IF ( IFLAG.EQ.1 ) GOTO 100
         NQCLU_I = NQCLU - NQCLU_L
C
C................................................................ R BULK
         IQ1 = NQ - NQ_R + 1
         IQ2 = NQ
         IBV3 = 5
         CALL CLUSSITES_ADD(IFLAG,CLURAD,MOL,NMIN_R,NMAX_R,ABAS_R,DQ,
     &                      QBAS,IQ1,IQ2,IQCNTR,RQCLU,DQCLU,N5VEC_QCLU,
     &                      IBV3,IQ_IQCLU,NQCLU,NQCLUMAX,NQMAX)
         IF ( IFLAG.EQ.1 ) GOTO 100
         NQCLU_R = NQCLU - NQCLU_L - NQCLU_I
C
      END IF
C=======================================================================
C
C
C ------------------------------  sort vectors in order of increasing D
C
      NSHLCLU = 1
      NQSHLCLU(1) = 1
      DSH(1) = 0D0
C
      DO I = 2,NQCLU
         K = I
         DTMP = DQCLU(I)
         ZTMP = RQCLU(3,I)
         YTMP = RQCLU(2,I)
         XTMP = RQCLU(1,I)
C
C -----------------------  look for shortest distance DTMP to the center
         DO J = I + 1,NQCLU
C
            IF ( DQCLU(J).LT.(DTMP-TOL) ) THEN
               K = J
               DTMP = DQCLU(K)
               ZTMP = RQCLU(3,K)
               YTMP = RQCLU(2,K)
               XTMP = RQCLU(1,K)
            ELSE IF ( ABS(DQCLU(J)-DTMP).LT.TOL ) THEN
               IF ( RQCLU(3,J).LT.(ZTMP-TOL) ) THEN
                  K = J
                  ZTMP = RQCLU(3,K)
                  YTMP = RQCLU(2,K)
                  XTMP = RQCLU(1,K)
               ELSE IF ( ABS(RQCLU(3,J)-ZTMP).LT.TOL ) THEN
                  IF ( RQCLU(2,J).LT.(YTMP-TOL) ) THEN
                     K = J
                     YTMP = RQCLU(2,K)
                     XTMP = RQCLU(1,K)
                  ELSE IF ( ABS(RQCLU(2,J)-YTMP).LT.TOL ) THEN
                     IF ( RQCLU(1,J).LT.(XTMP-TOL) ) THEN
                        K = J
                        XTMP = RQCLU(1,K)
                     END IF
                  END IF
               END IF
            END IF
C
         END DO
C
C --------------------------------------  interchange sites if necessary
         IF ( K.NE.I ) THEN
            DQCLU(K) = DQCLU(I)
            DQCLU(I) = DTMP
C
            IQTMP = IQ_IQCLU(K)
            IQ_IQCLU(K) = IQ_IQCLU(I)
            IQ_IQCLU(I) = IQTMP
C
            RVEC_TMP(1:3) = RQCLU(1:3,K)
            RQCLU(1:3,K) = RQCLU(1:3,I)
            RQCLU(1:3,I) = RVEC_TMP(1:3)
C
            N5VEC_TMP(1:5) = N5VEC_QCLU(1:5,K)
            N5VEC_QCLU(1:5,K) = N5VEC_QCLU(1:5,I)
            N5VEC_QCLU(1:5,I) = N5VEC_TMP(1:5)
         END IF
C
         IF ( I.GT.1 ) THEN
            IF ( ABS(DQCLU(I)-DQCLU(I-1)).GT.TOL ) THEN
               IF ( IPRINT.GT.1 ) WRITE (6,99001) NSHLCLU,
     &              NQSHLCLU(NSHLCLU),DQCLU(I-1)
               NSHLCLU = NSHLCLU + 1
               IF ( NSHLCLU.GT.NSHLCLUMAX ) THEN
                  WRITE (6,*) ' array dimension  NSHLCLUMAX exceeded',
     &                        NSHLCLUMAX
                  CALL STOP_MESSAGE(ROUTINE,'array dimension')
               END IF
               NQSHLCLU(NSHLCLU) = 1
               DSH(NSHLCLU) = DQCLU(I)
            ELSE
               NQSHLCLU(NSHLCLU) = NQSHLCLU(NSHLCLU) + 1
            END IF
         END IF
C
         IF ( IPRINT.GT.2 ) THEN
            IQ = IQ_IQCLU(I)
            WRITE (6,99002) I,IQ,(RQCLU(JJ,I),JJ=1,3)
         END IF
C
      END DO
C
      IF ( MOL ) THEN
         CLURAD = DQCLU(NQCLU) + 1D-6
         IF ( NQ.NE.NQCLU ) THEN
            WRITE (6,*) ' for MOL: NQ=',NQ,'  NQCLU=',NQCLU
            CALL STOP_MESSAGE(ROUTINE,'NQ.NE.NQCLU')
         END IF
      END IF
C
      IF ( IPRINT.GE.0 ) THEN
         IF ( IPRINT.GE.1 ) THEN
            DO ISH = 1,NSHLCLU
               WRITE (6,99006) ISH,NQSHLCLU(ISH),DSH(ISH)
            END DO
         END IF
         WRITE (6,99007) CLURAD,IQCNTR,NSHLCLU
         IF ( SYSTEM_DIMENSION(1:2).NE.'3D' ) WRITE (6,99009) NQCLU_L,
     &        NQCLU_I,NQCLU_R
         WRITE (6,99008) NQCLU
      END IF
C
C ----------------------------------------------------------------------
C  write results to temporary file for data transfer to calling routine
C ----------------------------------------------------------------------
C
      DO I = 1,NQCLU
         DQCLU(I) = DNRM2(3,RQCLU(1,I),1)
      END DO
C
      CALL OPEN_IOTMP_SCRATCH(ROUTINE,IOTMP)
      WRITE (IOTMP) ((RQCLU(J,I),J=1,3),DQCLU(I),IQ_IQCLU(I),I=1,NQCLU),
     &              (NQSHLCLU(I),I=1,NSHLCLU)
      WRITE (IOTMP) ((N5VEC_QCLU(J,I),J=1,5),I=1,NQCLU)
      REWIND IOTMP
C
C ----------------------------------------------------------------------
      IF ( ALLOCATED(DSH) ) DEALLOCATE (DSH,STAT=IA_ERR)
C
      DEALLOCATE (RQCLU,DQCLU,IQ_IQCLU,NQSHLCLU,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC')
C
C ----------------------------------------------------------------------
99001 FORMAT (10X,2I10,f10.6)
99002 FORMAT (5X,'site',I4,'  IQ=',I3,'   ->R=',3F7.3)
99003 FORMAT (10X,A,3I10,F10.3)
99004 FORMAT (10X,A,4F10.3)
99005 FORMAT (/,1X,79('*'),/,34X,'<CLUSSITES>',/,1X,79('*'),/)
99006 FORMAT (10X,'shell ',I3,' with ',I3,' sites ',' at distance ',
     &        F10.4)
99007 FORMAT (/,10X,'cluster of radius  ',F7.3,' a ',
     &        '  around central site  IQ =',I5,:,/,10X,
     &        'number of shells in cluster      NSHLCLU  = ',I5)
99008 FORMAT (10X,'total number of sites in cluster NQCLU    = ',I5,/,
     &        10X,'INCLUDING central site !',/)
99009 FORMAT (10X,'number of cluster sites  NQCLU  L | I | R = ',3I5)
      END
C*==clussites_add.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE CLUSSITES_ADD(IFLAG,CLURAD,MOL,NMIN,NMAX,ABAS,DQ,QBAS,
     &                         IQ1,IQ2,IQCNTR,RQCLU,DQCLU,N5VEC_QCLU,
     &                         IBV3,IQ_IQCLU,NQCLU,NQCLUMAX,NQMAX)
C   ********************************************************************
C   *                                                                  *
C   *   auxilary subroutine for <CLUSSITES> to cluster sites to table  *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 CLURAD
      INTEGER IBV3,IFLAG,IQ1,IQ2,IQCNTR,NQCLU,NQCLUMAX,NQMAX
      LOGICAL MOL
      REAL*8 ABAS(3,3),DQ(3),DQCLU(NQCLUMAX),QBAS(3,NQMAX),
     &       RQCLU(3,NQCLUMAX)
      INTEGER IQ_IQCLU(NQCLUMAX),N5VEC_QCLU(5,NQCLUMAX),NMAX(3),NMIN(3)
C
C Local variables
C
      REAL*8 CLURADSQ,DC(3),DRSQ
      REAL*8 DNRM2
      INTEGER I,IQ,L,M,N
C
C*** End of declarations rewritten by SPAG
C
      IFLAG = 0
      CLURADSQ = CLURAD*CLURAD
C
      DO IQ = IQ1,IQ2
C
         DO I = 1,3
            DQ(I) = QBAS(I,IQ) - QBAS(I,IQCNTR)
         END DO
C
         DO N = NMIN(3),NMAX(3)
            DO M = NMIN(2),NMAX(2)
               DO L = NMIN(1),NMAX(1)
C
                  DRSQ = 0D0
                  DO I = 1,3
                     DC(I) = L*ABAS(I,1) + M*ABAS(I,2) + N*ABAS(I,3)
     &                       + DQ(I)
                     DRSQ = DRSQ + DC(I)*DC(I)
                  END DO
C
                  IF ( DRSQ.GT.1D-3 .AND. (DRSQ.LE.CLURADSQ .OR. MOL) )
     &                 THEN
                     NQCLU = NQCLU + 1
                     IF ( NQCLU.LE.NQCLUMAX ) THEN
                        RQCLU(1:3,NQCLU) = DC(1:3)
                        DQCLU(NQCLU) = DNRM2(3,DC,1)
                        IQ_IQCLU(NQCLU) = IQ
                        N5VEC_QCLU(1,NQCLU) = L
                        N5VEC_QCLU(2,NQCLU) = M
                        N5VEC_QCLU(IBV3,NQCLU) = N
                     ELSE
                        WRITE (6,99001) NQCLUMAX
                        IFLAG = 1
                        RETURN
                     END IF
                  END IF
               END DO
            END DO
         END DO
C
      END DO
C
99001 FORMAT (/,5X,'WARNING: NQCLUMAX too small for  NQCLU=',I8,/,5X,
     &        'WARNING: NQCLUMAX will be increased')
      END
