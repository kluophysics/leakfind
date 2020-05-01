C*==symsites.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE SYMSITES(IQCNTR,MOL,IPRINT,NQ,QBAS0,ABAS,NONMAG,IREL,
     &                    KMROT,QMVEC,MROTR,NSYM,SYMSYMBL,ITBOT,ITTOP,
     &                    IQAT,NAT,QMPHI,QMTET,NTMAX,NQMAX,NSYMACCEPTED,
     &                    SYMACCEPTED,SYMUNITARY)
C   ********************************************************************
C   *                                                                  *
C   *                KMROT = 3    QMVEC = +/- QMVEC'   TET=90          *
C   *                KMROT = 4    QMVEC =     QMVEC'   TET<>90         *
C   *                implies IREL = 2                                  *
C   *                                                                  *
C   *   this subroutine is a modified version of SYMLATTICE            *
C   *   it finds out the allowed symmetry operations w.r.t one         *
C   *   of the sites  IQCNTR, i.e. its position is seen as the origin  *
C   *   this procedure gives the  LOCAL  site related point group      *
C   *   i.e. only operations of the type {R|0}                         *
C   *                                                                  *
C   *   NOTE: the possible symmetry point operations MROTR have        *
C   *         to be supplied as they should not change                 *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:PI
      USE MOD_SYMMETRY,ONLY:NSYMMAX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 TOL
      PARAMETER (TOL=1.0D-6)
C
C Dummy arguments
C
      INTEGER IPRINT,IQCNTR,IREL,ITBOT,ITTOP,KMROT,NQ,NQMAX,NSYM,
     &        NSYMACCEPTED,NTMAX
      LOGICAL MOL,NONMAG
      REAL*8 ABAS(3,3),MROTR(3,3,NSYMMAX),QBAS0(3,NQMAX),QMPHI(NQMAX),
     &       QMTET(NQMAX),QMVEC(3)
      INTEGER IQAT(NQMAX,NTMAX),NAT(NTMAX)
      LOGICAL SYMACCEPTED(NSYMMAX),SYMUNITARY(NSYMMAX)
      CHARACTER*4 SYMSYMBL(NSYMMAX)
C
C Local variables
C
      REAL*8 AUX,BR(3,3),BRINV(3,3),BRMAT(3,3),BRP(3),BV(3),CF(3),
     &       MAGMOM,MDOTMP,MQMVEC(3),MVECQ(:,:),MVECQP(:,:),QBAS(:,:),
     &       QJP(3),QMVECP(3),QSFT(3),QVEC(:,:),QVECP(:,:),STET
      REAL*8 DDOT,RMAT3X3DET
      INTEGER I,IA,IA_ERR,ICL,ICLQ0(:),IO,IQ,IQPSYMQ(:,:),IQREPCL0(:),
     &        IROT,ISYM,IT,IT_OQAUX(:,:),J,JCL,JQ,LQ,NCL,NO_QAUX(:),
     &        NQOK,NSYMH,SYMDET(NSYMMAX)
      LOGICAL LATSITE,SAMECL,SYMCRYSYS(NSYMMAX),SYMOPER
      LOGICAL RVEC_SAME
      CHARACTER*1 STRA,STRB
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE QVEC,ICLQ0,MVECQ,QVECP
      ALLOCATABLE IQPSYMQ,MVECQP,IQREPCL0,QBAS
      ALLOCATABLE IT_OQAUX,NO_QAUX
C
      ALLOCATE (IT_OQAUX(NTMAX,NQMAX),NO_QAUX(NQMAX))
      ALLOCATE (QVEC(3,NQMAX),ICLQ0(NQMAX))
      ALLOCATE (IQPSYMQ(NSYMMAX,NQMAX),MVECQ(3,NQMAX))
      ALLOCATE (QVECP(3,NQMAX),MVECQP(3,NQMAX))
      ALLOCATE (IQREPCL0(NQMAX),QBAS(3,NQMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:symsites.f->IQREPCL0'
C
      IF ( IPRINT.GE.0 ) WRITE (6,99001) IQCNTR
C
C-----------------------------------------------------------------------
C     fix temporarily the occupation information   NO_QAUX and  IT_OQAUX
C-----------------------------------------------------------------------
      DO IQ = 1,NQ
C
         NO_QAUX(IQ) = 0
         DO IT = ITBOT,ITTOP
            DO IA = 1,NAT(IT)
               IF ( IQAT(IA,IT).EQ.IQ ) THEN
                  NO_QAUX(IQ) = NO_QAUX(IQ) + 1
                  IT_OQAUX(NO_QAUX(IQ),IQ) = IT
               END IF
            END DO
         END DO
      END DO
C
      NSYMH = NSYM/2
C-----------------------------------------------------------------------
C                                     shift to temporary origin at QCNTR
      QSFT(1:3) = QBAS0(1:3,IQCNTR)
C
      DO IQ = 1,NQ
         QBAS(1:3,IQ) = QBAS0(1:3,IQ) - QSFT(1:3)
      END DO
C
C-----------------------------------------------------------------------
C                                      BR    basis vectors in real space
      BR(1:3,1:3) = ABAS(1:3,1:3)
C
      DO J = 1,3
         DO I = 1,3
            BRMAT(I,J) = DDOT(3,BR(1,I),1,BR(1,J),1)
         END DO
      END DO
C
      CALL RINVGJ(BRINV,BRMAT,3,3)
C
C-----------------------------------------------------------------------
C          rotate the 3 primitive lattice vectors and check
C          whether the rotated vector is a proper lattice vector
C          if not:  exclude operation
C-----------------------------------------------------------------------
C
      SYMCRYSYS(1:NSYM) = .TRUE.
C
      DO IROT = 1,NSYM
C
         DO I = 1,3
C
            CALL DGEMV('N',3,3,1D0,MROTR(1,1,IROT),3,BR(1,I),1,0D0,BRP,
     &                 1)
C
            DO J = 1,3
               BV(J) = DDOT(3,BR(1,J),1,BRP,1)
            END DO
            CALL DGEMV('N',3,3,1D0,BRINV,3,BV,1,0D0,CF,1)
C
            DO J = 1,3
C
               IF ( ABS(NINT(CF(J))-CF(J)).GT.1D-8 ) SYMCRYSYS(IROT)
     &              = .FALSE.
C
            END DO
C
         END DO
C
      END DO
C
C-----------------------------------------------------------------------
C             find class  CL  for every atomic site  Q
C                  the magnetisation is ignored
C-----------------------------------------------------------------------
      NCL = 1
      ICLQ0(1) = 1
      IQREPCL0(1) = 1
      DO JQ = 2,NQ
         JCL = 0
         DO ICL = 1,NCL
            IF ( JCL.EQ.0 ) THEN
               SAMECL = .FALSE.
               IQ = IQREPCL0(ICL)
               IF ( NO_QAUX(IQ).EQ.NO_QAUX(JQ) ) THEN
                  SAMECL = .TRUE.
                  DO IO = 1,NO_QAUX(IQ)
                     IF ( IT_OQAUX(IO,IQ).NE.IT_OQAUX(IO,JQ) )
     &                    SAMECL = .FALSE.
                  END DO
               END IF
               IF ( SAMECL ) JCL = ICL
            END IF
         END DO
         IF ( JCL.EQ.0 ) THEN
            NCL = NCL + 1
            IQREPCL0(NCL) = JQ
            JCL = NCL
         END IF
         ICLQ0(JQ) = JCL
      END DO
C
C-----------------------------------------------------------------------
C                find the allowed symmetry operations
C-----------------------------------------------------------------------
C
      IF ( NONMAG .OR. (IREL.LT.3) ) THEN
         MAGMOM = 0D0
      ELSE
         MAGMOM = 1D0
      END IF
C
      MQMVEC(1:3) = -QMVEC(1:3)
C
      DO IQ = 1,NQ
C
         QVEC(1:3,IQ) = QBAS(1:3,IQ)
C
         STET = SIN(QMTET(IQ)*PI/180D0)
         MVECQ(1,IQ) = MAGMOM*STET*COS(QMPHI(IQ)*PI/180D0)
         MVECQ(2,IQ) = MAGMOM*STET*SIN(QMPHI(IQ)*PI/180D0)
         MVECQ(3,IQ) = MAGMOM*COS(QMTET(IQ)*PI/180D0)
C
      END DO
C
C-----------------------------------------------------------------------
C       find equivalent position QVEC in parallelepiped spanned by BR''s
C
      IF ( .NOT.MOL ) THEN
         DO IQ = 1,NQ
C
            DO I = 1,3
               BV(I) = DDOT(3,BR(1,I),1,QVEC(1,IQ),1)
            END DO
C
            CALL DGEMV('N',3,3,1D0,BRINV,3,BV,1,0D0,CF,1)
C
            QVEC(1:3,IQ) = 0D0
            DO I = 1,3
               CF(I) = CF(I) - INT(CF(I)+1000D0) + 1000D0
               IF ( ABS(CF(I)-1D0).LT.TOL ) CF(I) = 0D0
               CALL DAXPY(3,CF(I),BR(1,I),1,QVEC(1,IQ),1)
            END DO
C
         END DO
      END IF
C
      NSYMACCEPTED = 0
      MDOTMP = 1D0
C
      DO ISYM = 1,NSYM
C
         SYMOPER = .FALSE.
C
         IF ( SYMCRYSYS(ISYM) ) THEN
C
            SYMDET(ISYM) = NINT(RMAT3X3DET(MROTR(1,1,ISYM)))
C
C-----------------------------------------------------------------------
C                   apply symmetry operation to spin spiral vector QMVEC
C
            IF ( KMROT.GE.3 ) THEN
               CALL DGEMV('N',3,3,1D0,MROTR(1,1,ISYM),3,QMVEC,1,0D0,
     &                    QMVECP,1)
               IF ( .NOT.RVEC_SAME(3,QMVEC,QMVECP,1D-7) ) THEN
                  IF ( KMROT.EQ.4 ) GOTO 50
                  IF ( .NOT.RVEC_SAME(3,MQMVEC,QMVECP,1D-7) ) GOTO 50
               END IF
            END IF
C
C-----------------------------------------------------------------------
C                     apply symmetry operation to every  QVEC  and  MVEC
            DO IQ = 1,NQ
C
               CALL DGEMV('N',3,3,1D0,MROTR(1,1,ISYM),3,QVEC(1,IQ),1,
     &                    0D0,QVECP(1,IQ),1)
C
               CALL DGEMV('N',3,3,1D0,MROTR(1,1,ISYM),3,MVECQ(1,IQ),1,
     &                    0D0,MVECQP(1,IQ),1)
C
               CALL DSCAL(3,DBLE(SYMDET(ISYM)),MVECQP(1,IQ),1)
C
               IF ( NONMAG .OR. (IREL.LT.3) ) THEN
                  MDOTMP = 1D0
               ELSE
                  MDOTMP = DDOT(3,MVECQ(1,IQ),1,MVECQP(1,IQ),1)
               END IF
C
C-----------------------------------------------------------------------
C                     reject symmetry operations that do not keep M-axis
               IF ( ABS(ABS(MDOTMP)-1D0).GT.TOL ) GOTO 50
C
            END DO
C
C-----------------------------------------------------------------------
C                                             in contrast to SYMLATTICE:
C                         don't allow for non-primitive translation QSFT
C
            MDOTMP = DDOT(3,MVECQ(1,1),1,MVECQP(1,1),1)
C
            NQOK = 0
            DO JQ = 1,NQ
               QJP(1:3) = QVECP(1:3,JQ)
C
C-------------------------- bring QJP in parallelepiped spanned by BR''s
C
               IF ( .NOT.MOL ) THEN
                  DO I = 1,3
                     BV(I) = DDOT(3,BR(1,I),1,QJP,1)
                  END DO
                  CALL DGEMV('N',3,3,1D0,BRINV,3,BV,1,0D0,CF,1)
                  QJP(1:3) = 0D0
                  DO I = 1,3
                     CF(I) = CF(I) - INT(CF(I)+1000D0) + 1000D0
                     IF ( ABS(CF(I)-1D0).LT.TOL ) CF(I) = 0D0
                     CALL DAXPY(3,CF(I),BR(1,I),1,QJP,1)
                  END DO
               END IF
C
               LATSITE = .FALSE.
               DO LQ = 1,NQ
                  IF ( .NOT.LATSITE ) THEN
                     IF ( ICLQ0(JQ).EQ.ICLQ0(LQ) ) THEN
                        AUX = DDOT(3,MVECQ(1,LQ),1,MVECQP(1,JQ),1)
                        IF ( NONMAG .OR. ABS(AUX-MDOTMP).LE.TOL ) THEN
                           IF ( RVEC_SAME(3,QVEC(1,LQ),QJP,1D-7) ) THEN
                              LATSITE = .TRUE.
                              NQOK = NQOK + 1
                              IQPSYMQ(ISYM,JQ) = LQ
                           END IF
                        END IF
                     END IF
                  END IF
               END DO
            END DO
            IF ( NQOK.EQ.NQ ) SYMOPER = .TRUE.
C
         END IF
C
C-----------------------------------------------------------------------
C
C
 50      CONTINUE
         SYMACCEPTED(ISYM) = SYMOPER
         IF ( MDOTMP.GT.-1D-6 ) THEN
            SYMUNITARY(ISYM) = .TRUE.
         ELSE
            SYMUNITARY(ISYM) = .FALSE.
         END IF
C
         IF ( SYMOPER ) NSYMACCEPTED = NSYMACCEPTED + 1
      END DO
C
      IF ( IPRINT.GE.0 ) THEN
         WRITE (6,99002) NSYM,NSYMACCEPTED
         WRITE (6,99003)
C
         DO ISYM = 1,NSYM
            IF ( SYMACCEPTED(ISYM) ) THEN
               IF ( ISYM.LE.NSYMH ) THEN
                  STRA = ' '
               ELSE
                  STRA = 'I'
               END IF
               IF ( SYMUNITARY(ISYM) ) THEN
                  STRB = ' '
               ELSE
                  STRB = 'T'
               END IF
C
               WRITE (6,99004) ISYM,STRA,SYMSYMBL(ISYM),STRB,
     &                         (IQPSYMQ(ISYM,IQ),IQ=1,MIN(5,NQ))
               IF ( NQ.GT.6 ) WRITE (6,99005) (IQPSYMQ(ISYM,IQ),IQ=6,NQ)
            END IF
         END DO
      END IF
C
      DEALLOCATE (QVEC,ICLQ0,MVECQ,QVECP)
      DEALLOCATE (MVECQP,IQREPCL0)
      RETURN
99001 FORMAT (//,1X,79('*'),/,34X,'<SYMSITES>',/,1X,79('*'),//,10X,
     &        'find out local point symmetry of site IQCNTR =',I4,/)
99002 FORMAT (/,10X,'crystal system allowed symmetry operations ',I8,/,
     &        10X,'accepted symmetry operations     ',I18,/)
99003 FORMAT (/,7X,'I',3X,'INV SYMBOL TIM',3X,'Q'' ')
99004 FORMAT (5X,I3,1X,3X,A,3X,A,3X,A,2X,(5I4))
99005 FORMAT ((59X,5I4))
      END
