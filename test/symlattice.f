C*==symlattice.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE SYMLATTICE(IWR,MOL,IPRINT,NWEDGE,IWEDGEROT,NQ,QBAS,
     &                      ABAS,NONMAG,IREL,KMROT,QMVEC,MROTR,MROTK,
     &                      NSYM,NSYMCRYSYS,SYMACCEPTED,SYMUNITARY,
     &                      SYMCRYSYS,SYMTVEC,SYMSYMBL,SYMEULANG,NOQ,
     &                      ITOQ,IT0,NT0,NT,NAT,IQAT,QMPHI,QMTET,
     &                      IQORGQP,IQPSYMQ,NSFTSYMQ,NSYMACCEPTED,
     &                      SYMDET,ISYMGENQ,IQREPQ,NQEQ,IQEQ,NMSYM,
     &                      IQREPMSYM,NTMAX,NQMAX)
C   ********************************************************************
C   *                                                                  *
C   *                KMROT = 3    QMVEC = +/- QMVEC'   TET=90          *
C   *                KMROT = 4    QMVEC =     QMVEC'   TET<>90         *
C   *                implies IREL = 2                                  *
C   *                                                                  *
C   * 16/06/05  check compatibility of symmetry operations             *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:PI,CI
      USE MOD_SYMMETRY,ONLY:NSYMMAX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SYMLATTICE')
      REAL*8 TOL
      PARAMETER (TOL=1.0D-6)
      INTEGER NLROT,NLMROT
      PARAMETER (NLROT=2,NLMROT=NLROT**2)
C
C Dummy arguments
C
      INTEGER IPRINT,IREL,IWR,KMROT,NMSYM,NQ,NQMAX,NSYM,NSYMACCEPTED,
     &        NSYMCRYSYS,NT,NT0,NTMAX,NWEDGE
      LOGICAL MOL,NONMAG
      REAL*8 ABAS(3,3),MROTK(3,3,NSYMMAX),MROTR(3,3,NSYMMAX),
     &       QBAS(3,NQMAX),QMPHI(NQMAX),QMTET(NQMAX),QMVEC(3),
     &       SYMEULANG(3,NSYMMAX),SYMTVEC(3,NSYMMAX)
      INTEGER IQAT(NQMAX,NTMAX),IQEQ(NQMAX,NQMAX),IQORGQP(NSYMMAX,NQMAX)
     &        ,IQPSYMQ(NSYMMAX,NQMAX),IQREPMSYM(NQMAX),IQREPQ(NQMAX),
     &        ISYMGENQ(NQMAX),IT0(NTMAX),ITOQ(NTMAX,NQMAX),
     &        IWEDGEROT(NSYMMAX),NAT(NTMAX),NOQ(NQMAX),NQEQ(NQMAX),
     &        NSFTSYMQ(3,NSYMMAX,NQMAX),SYMDET(NSYMMAX)
      LOGICAL SYMACCEPTED(NSYMMAX),SYMCRYSYS(NSYMMAX),
     &        SYMUNITARY(NSYMMAX)
      CHARACTER*4 SYMSYMBL(NSYMMAX)
C
C Local variables
C
      REAL*8 AUX,BR(3,3),BRINV(3,3),BRMAT(3,3),BRP(3),BV(3),CF(3),CF0(3)
     &       ,M3WK1(3,3),M3WK2(3,3),M4E(4,4),M4WK1(4,4),MAGMOM,MDOTMP,
     &       MQMVEC(3),MVECQ(:,:),MVECQP(:,:),QJP(3),QMVECP(3),QSFT(3),
     &       QVEC(:,:),QVECP(:,:),QVECPP(3),QVECTST(3),RMAT3(3,3),
     &       SGOP(:,:,:),SHIFT,STET,V(3),VP(3),W
      COMPLEX*16 CS,DSYM_CLM(:,:,:),USC(3,3),W3X3(3,3)
      REAL*8 DDOT,RMAT3X3DET
      LOGICAL FOUND,FOUND_INV,LATSITE,SAMECL,SYMOPER,WEDGEOK(NSYMMAX)
      INTEGER I,IA,IA_ERR,ICL,ICL1,ICLA,ICLOK,ICLQ(:),ICLQ0(:),IEQ,
     &        IFLAG,IINV,IMB,IO,IQ,IQAT0(:,:),IQP,IQREP,IQREPCL0(:),
     &        IQ_MBCL(:,:),IROT,IROTP,IRUN,ISYM,ISYMP,IT,J,JCL,JQ,JSYM,
     &        JT,K,KSYM,LQ,NAT0(:),NCL,NERR,NMBMAX,NMB_CL(:),NOQ0(:),
     &        NQOK,NSYMCRYSYS0,NSYMH,NTADD,NTP
      LOGICAL RVEC_SAME
      CHARACTER*1 STRA,STRB
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE ICLQ,IQ_MBCL,NMB_CL,QVEC,ICLQ0,IQAT0,MVECQ,QVECP,NOQ0
      ALLOCATABLE MVECQP,NAT0,IQREPCL0,DSYM_CLM,SGOP
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (SGOP(4,4,NSYMMAX))
      ALLOCATE (ICLQ(NQMAX),IQ_MBCL(NQMAX,NQMAX),NOQ0(NQMAX))
      ALLOCATE (NMB_CL(NQMAX),QVEC(3,NQMAX),ICLQ0(NQMAX))
      ALLOCATE (IQAT0(NQMAX,NTMAX),MVECQ(3,NQMAX))
      ALLOCATE (QVECP(3,NQMAX),MVECQP(3,NQMAX),NAT0(NTMAX))
      ALLOCATE (DSYM_CLM(NLMROT,NLMROT,NSYMMAX))
      ALLOCATE (IQREPCL0(NQMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'alloc: IQREPCL0')
C
      IF ( IPRINT.GE.0 ) THEN
         IF ( MOL ) THEN
            WRITE (6,99010) 'a cluster'
         ELSE
            WRITE (6,99010) 'the lattice'
         END IF
      END IF
C
C-----------------------------------------------------------------------
C  create transformation matrix   U  cartesian/spherical ccordinates
C-----------------------------------------------------------------------
C  RC,RCP  vectors in cartesian coordinates
C  RS,RSP  vectors in spherical coordinates
C         RS  = USC * RC                                 (4.40)
C         RSP = MS  * RS                                 (4.37)
C     MS(i,j) = D(j,i)                                   (4.42)
C     D  rotation matrix for complex spherical harmonics
C
      W = 1.0D0/SQRT(2.0D0)
C
C ordering of: m=-1,0,+1 >>> row 1 and 3 interchanged compared to (4.44)
      USC(1,1) = W
      USC(1,2) = -CI*W
      USC(1,3) = 0.0D0
      USC(2,1) = 0.0D0
      USC(2,2) = 0.0D0
      USC(2,3) = 1.0D0
      USC(3,1) = -W
      USC(3,2) = -CI*W
      USC(3,3) = 0.0D0
C-----------------------------------------------------------------------
C
      NSYMH = NSYM/2
C
C-----------------------------------------------------------------------
C create the rotation matrices DSYM_CLM for complex spherical harmonics
C-----------------------------------------------------------------------
C
      LOOP_CREATE_SYMMETRY_MATRICES:DO ISYM = 1,NSYMH
C
         CALL ROTMAT(NLROT,1,SYMEULANG(1,ISYM),SYMEULANG(2,ISYM),
     &               SYMEULANG(3,ISYM),DSYM_CLM(1,1,ISYM),NLMROT)
C
         IF ( IPRINT.GT.1 ) THEN
C
            IF ( ISYM.EQ.1 ) WRITE (6,99011) NSYM/2,NSYM
C
            WRITE (6,99012) ISYM,(SYMEULANG(I,ISYM),I=1,3)
C
            CALL CMATSTRUCT('4x4 rotation matrix D for l<=1 '//
     &                      '(spherical basis)',DSYM_CLM(1,1,ISYM),
     &                      NLMROT,NLMROT,1,1,0,1D-8,6)
         END IF
C
      END DO LOOP_CREATE_SYMMETRY_MATRICES
C
C-----------------------------------------------------------------------
C create the rotation matrix MROTR for vectors in cartesian coordinates
C NOTE:  U^+ D^T U gives the inverse of the real matrix  M
C        for that reason  the transposed matrix is stored as MROTR(J,I)
C-----------------------------------------------------------------------
      LOOP_CREATE_THE_SYMMETRY_MATRIX_MROTR:DO ISYM = 1,NSYMH
C
         DO I = 1,3
            DO J = 1,3
               CS = 0.0D0
               DO K = 1,3
                  CS = CS + DSYM_CLM(K+1,I+1,ISYM)*USC(K,J)
               END DO
               W3X3(I,J) = CS
            END DO
         END DO
C
         DO I = 1,3
            DO J = 1,3
               CS = 0.0D0
               DO K = 1,3
                  CS = CS + DCONJG(USC(K,I))*W3X3(K,J)
               END DO
               IF ( DIMAG(CS).GT.1D-8 ) WRITE (6,*) 'ISYM=',ISYM,
     &              ' MROT',I,J,CS,' ???????????'
C-------------------------------------------------------- see NOTE above
               MROTR(J,I,ISYM) = DREAL(CS)
            END DO
         END DO
C
C-------------------------------- compare with M3WK1 calculated directly
C
         CALL GETMROT(SYMEULANG(1,ISYM),SYMEULANG(2,ISYM),
     &                SYMEULANG(3,ISYM),M3WK1)
C
         IF ( .NOT.RVEC_SAME(9,MROTR(1,1,ISYM),M3WK1,1D-8) ) THEN
            WRITE (6,*) 'ERROR for ISYM = ',ISYM
            WRITE (6,'(/,3(3F10.6,/))') ((MROTR(I,J,ISYM),J=1,3),I=1,3)
            WRITE (6,'(/,3(3F10.6,/))') ((M3WK1(I,J),J=1,3),I=1,3)
            CALL STOP_MESSAGE(ROUTINE,'comparison with M3WK1')
         END IF
C
         IF ( IPRINT.GT.0 ) THEN
            IF ( ISYM.EQ.1 ) THEN
               WRITE (6,99006)
               V(1) = 1D0
               V(2) = 2D0
               V(3) = 3D0
            END IF
C
            CALL DGEMV('N',3,3,1D0,MROTR(1,1,ISYM),3,V,1,0D0,VP,1)
C
            CALL DGEMM('N','N',3,3,3,1D0,MROTR(1,1,ISYM),3,ABAS,3,0D0,
     &                 M3WK1,3)
C
            CALL DGEMM('T','N',3,3,3,1D0,ABAS,3,M3WK1,3,0D0,M3WK2,3)
C
            WRITE (6,99007) ISYM,(NINT(SYMEULANG(I,ISYM)),I=1,3),
     &                      SYMSYMBL(ISYM),
     &                      ((MROTR(I,J,ISYM),J=1,3),VP(I),
     &                      (M3WK2(I,J),J=1,3),I=1,3)
C
         END IF
C
      END DO LOOP_CREATE_THE_SYMMETRY_MATRIX_MROTR
C
C-----------------------------------------------------------------------
C                     create matrix for inversion
C-----------------------------------------------------------------------
      IINV = NSYMH + 1
C
      CALL DGEMM('N','N',3,3,3,-1D0,MROTR(1,1,1),3,MROTR(1,1,1),3,0D0,
     &           MROTR(1,1,IINV),3)
C
C-----------------------------------------------------------------------
C                         include inversion
C-----------------------------------------------------------------------
C
      DO ISYM = 2,NSYMH
         ISYMP = NSYMH + ISYM
C
         CALL DGEMM('N','N',3,3,3,1D0,MROTR(1,1,ISYM),3,MROTR(1,1,IINV),
     &              3,0D0,MROTR(1,1,ISYMP),3)
C
      END DO
C
C=======================================================================
C            set up of transformation matrices completed
C=======================================================================
C
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
C-----------------------------------------------------------------------
C
      IFLAG = 0
      NSYMCRYSYS0 = NSYMCRYSYS
      IRUN = 1
C
C
 100  CONTINUE
      LOOP_ROTATE_PRIMITIVE_LATTICE_VECTOR:DO IROT = 1,NSYM
C
         IF ( SYMCRYSYS(IROT) ) THEN
C
            DO I = 1,3
C
               CALL DGEMV('N',3,3,1D0,MROTR(1,1,IROT),3,BR(1,I),1,0D0,
     &                    BRP,1)
C
               DO J = 1,3
                  BV(J) = DDOT(3,BR(1,J),1,BRP,1)
               END DO
               CALL DGEMV('N',3,3,1D0,BRINV,3,BV,1,0D0,CF,1)
C
C-----------------------------------------------------------------------
               DO J = 1,3
                  IF ( ABS(NINT(CF(J))-CF(J)).GT.1D-8 ) THEN
                     IF ( IRUN.EQ.1 ) THEN
                        NSYMCRYSYS = NSYM
                        SYMCRYSYS(1:NSYM) = .TRUE.
                        IRUN = 2
                        GOTO 100
                     END IF
C
                     IFLAG = 1
                     SYMCRYSYS(IROT) = .FALSE.
                     NSYMCRYSYS = NSYMCRYSYS - 1
                     IROTP = IROT + NSYM/2
                     IF ( IROTP.LE.NSYM ) THEN
                        SYMCRYSYS(IROTP) = .FALSE.
                        NSYMCRYSYS = NSYMCRYSYS - 1
                     END IF
                     GOTO 200
                  END IF
               END DO
C-----------------------------------------------------------------------
C
            END DO
C
         END IF
C
 200  END DO LOOP_ROTATE_PRIMITIVE_LATTICE_VECTOR
C
      IF ( IFLAG.EQ.1 ) WRITE (6,99026) NSYMCRYSYS0,NSYMCRYSYS
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
               IF ( NOQ(IQ).EQ.NOQ(JQ) ) THEN
                  SAMECL = .TRUE.
                  DO IO = 1,NOQ(IQ)
                     IF ( ITOQ(IO,IQ).NE.ITOQ(IO,JQ) ) SAMECL = .FALSE.
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
C                find the accepted symmetry operations
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
                  IF ( KMROT.EQ.4 ) GOTO 250
                  IF ( .NOT.RVEC_SAME(3,MQMVEC,QMVECP,1D-7) ) GOTO 250
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
               IF ( ABS(ABS(MDOTMP)-1D0).GT.TOL ) GOTO 250
C
            END DO
C
C-----------------------------------------------------------------------
C                           find possible non-primitive translation QSFT
C                     select site 1 and scan all sites of the same class
            DO IQ = 1,NQ
C
               IF ( (ICLQ0(1).EQ.ICLQ0(IQ)) .AND. (.NOT.SYMOPER) ) THEN
C
                  MDOTMP = DDOT(3,MVECQ(1,IQ),1,MVECQP(1,1),1)
C
                  IF ( MOL ) THEN
                     QSFT(1:3) = 0D0
                  ELSE
                     QSFT(1:3) = QVEC(1:3,IQ) - QVECP(1:3,1)
                  END IF
C
                  NQOK = 0
                  DO JQ = 1,NQ
                     QJP(1:3) = QVECP(1:3,JQ) + QSFT(1:3)
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
                              IF ( ABS(AUX-MDOTMP).LE.TOL ) THEN
                                 IF ( RVEC_SAME(3,QVEC(1,LQ),QJP,1D-7) )
     &                                THEN
                                    LATSITE = .TRUE.
                                    NQOK = NQOK + 1
                                    IQPSYMQ(ISYM,JQ) = LQ
                                    IQORGQP(ISYM,LQ) = JQ
                                 END IF
                              END IF
                           END IF
                        END IF
                     END DO
                  END DO
                  IF ( NQOK.EQ.NQ ) SYMOPER = .TRUE.
C
               END IF
            END DO
         END IF
C
C-----------------------------------------------------------------------
C
 250     CONTINUE
         SYMACCEPTED(ISYM) = SYMOPER
         IF ( MDOTMP.GT.-1D-6 ) THEN
            SYMUNITARY(ISYM) = .TRUE.
         ELSE
            SYMUNITARY(ISYM) = .FALSE.
         END IF
C
         IF ( SYMOPER ) THEN
C
            NSYMACCEPTED = NSYMACCEPTED + 1
C
            IF ( MOL ) THEN
               SYMTVEC(1:3,ISYM) = 0D0
            ELSE
               DO I = 1,3
                  BV(I) = DDOT(3,BR(1,I),1,QSFT,1)
               END DO
               CALL DGEMV('N',3,3,1D0,BRINV,3,BV,1,0D0,CF,1)
               SYMTVEC(1:3,ISYM) = 0D0
               DO I = 1,3
                  CF(I) = CF(I) - INT(CF(I)+1000D0) + 1000D0
                  IF ( ABS(CF(I)-1D0).LT.TOL ) CF(I) = 0D0
                  CALL DAXPY(3,CF(I),BR(1,I),1,SYMTVEC(1,ISYM),1)
               END DO
            END IF
C
         END IF
      END DO
C
      IF ( IPRINT.GE.0 ) THEN
         WRITE (6,99013) NSYMCRYSYS,NSYMACCEPTED
         IF ( IPRINT.GT.0 ) WRITE (6,99016)
         WRITE (6,99017)
C
         DO ISYM = 1,NSYM
            IF ( SYMACCEPTED(ISYM) ) THEN
               IF ( ISYM.LE.NSYMH ) THEN
                  JSYM = ISYM
                  STRA = ' '
               ELSE
                  JSYM = ISYM - NSYMH
                  STRA = 'I'
               END IF
               IF ( SYMUNITARY(ISYM) ) THEN
                  STRB = ' '
               ELSE
                  STRB = 'T'
               END IF
C
               WRITE (6,99018) ISYM,(NINT(SYMEULANG(I,JSYM)),I=1,3),
     &                         STRA,SYMSYMBL(ISYM),STRB,
     &                         (SYMTVEC(I,ISYM),I=1,3),
     &                         (IQPSYMQ(ISYM,IQ),IQ=1,MIN(5,NQ))
               IF ( NQ.GE.6 ) WRITE (6,99019) (IQPSYMQ(ISYM,IQ),IQ=6,NQ)
               IF ( IPRINT.GT.0 ) THEN
                  WRITE (6,99020) (IQORGQP(ISYM,IQ),IQ=1,NQ)
                  CALL DGEMV('N',3,3,1D0,MROTR(1,1,ISYM),3,
     &                       SYMTVEC(1,ISYM),1,0D0,VP,1)
                  WRITE (6,99021) SYMDET(ISYM),
     &                            (SYMTVEC(I,ISYM)-VP(I),I=1,3)
                  IF ( .NOT.RVEC_SAME(3,SYMTVEC(1,ISYM),VP,1D-8) )
     &                 WRITE (*,*) '???????????????????????'
               END IF
            END IF
         END DO
      END IF
C
C=======================================================================
C            check whether symmetry operations form a group
C=======================================================================
C
      IFLAG = 0
      M4E(1:4,1:4) = 0D0
      DO I = 1,4
         M4E(I,I) = 1D0
      END DO
C
      SGOP(1:4,1:4,1:NSYMMAX) = 0D0
C
      DO ISYM = 1,NSYM
         IF ( SYMACCEPTED(ISYM) ) THEN
            SGOP(1:3,1:3,ISYM) = MROTR(1:3,1:3,ISYM)
            SGOP(1:3,4,ISYM) = SYMTVEC(1:3,ISYM)
            SGOP(4,4,ISYM) = 1D0
         END IF
      END DO
C
      DO ISYM = 1,NSYM
         IF ( SYMACCEPTED(ISYM) ) THEN
            FOUND_INV = .FALSE.
C
            DO JSYM = 1,NSYM
               IF ( SYMACCEPTED(JSYM) ) THEN
C
C-----------------------------------------------------------------------
                  M4WK1(1:4,1:4) = MATMUL(SGOP(1:4,1:4,ISYM),SGOP(1:4,1:
     &                             4,JSYM))
C
                  DO I = 1,3
                     BV(I) = DDOT(3,BR(1,I),1,M4WK1(1,4),1)
                  END DO
C
                  CALL DGEMV('N',3,3,1D0,BRINV,3,BV,1,0D0,CF,1)
C
                  M4WK1(1:3,4) = 0D0
                  DO I = 1,3
                     CF(I) = CF(I) - INT(CF(I)+1000D0) + 1000D0
                     IF ( ABS(CF(I)-1D0).LT.TOL ) CF(I) = 0D0
                     CALL DAXPY(3,CF(I),BR(1,I),1,M4WK1(1,4),1)
                  END DO
C-----------------------------------------------------------------------
C
                  FOUND = .FALSE.
                  DO KSYM = 1,NSYM
                     IF ( SYMACCEPTED(KSYM) ) THEN
                        IF ( RVEC_SAME(16,SGOP(1,1,KSYM),M4WK1,1D-8) )
     &                       THEN
                           FOUND = .TRUE.
                           EXIT
                        END IF
                     END IF
                  END DO
                  IF ( FOUND ) THEN
                     IF ( RVEC_SAME(16,M4E,M4WK1,1D-8) )
     &                    FOUND_INV = .TRUE.
                  ELSE
                     IFLAG = 1
                  END IF
C
               END IF
            END DO
            IF ( .NOT.FOUND_INV ) THEN
               WRITE (6,99001) ISYM
               IFLAG = 1
            END IF
C
         END IF
      END DO
      IF ( IFLAG.EQ.1 ) THEN
         CALL STOP_MESSAGE(ROUTINE,
     &                     'symmetry operations don''t form a group')
      ELSE
         IF ( IPRINT.GE.0 ) WRITE (6,99002) NSYMACCEPTED
      END IF
C
C-----------------------------------------------------------------------
C      ONCE MORE:  find class  CL  for every atomic site  Q
C                  the magnetisation is accounted for now
C-----------------------------------------------------------------------
C
      NERR = 0
      NCL = 0
      NMB_CL(1) = 0
      IQ_MBCL(1,1) = 0
C
      DO IQ = 1,NQ
         ICLOK = 0
         DO ICL = 1,NCL
            DO I = 1,NMB_CL(ICL)
               IF ( IQ_MBCL(I,ICL).EQ.IQ ) THEN
                  ICLOK = ICLOK + 1
                  ICLQ(IQ) = ICL
               END IF
            END DO
         END DO
C
         IF ( ICLOK.EQ.0 ) THEN
            NCL = NCL + 1
            NMB_CL(NCL) = 1
            IQ_MBCL(1,NCL) = IQ
            ICLQ(IQ) = NCL
         ELSE IF ( ICLOK.GT.1 ) THEN
            WRITE (6,*) 'IQ found in more than 1 classes ICL '
            NERR = NERR + 1
         END IF
C
         DO ISYM = 1,NSYM
            IF ( SYMACCEPTED(ISYM) ) THEN
               JQ = IQPSYMQ(ISYM,IQ)
               ICLOK = 0
               DO ICL = 1,NCL
                  DO I = 1,NMB_CL(ICL)
                     IF ( IQ_MBCL(I,ICL).EQ.JQ ) THEN
                        ICLOK = ICLOK + 1
                        ICLQ(JQ) = ICL
                     END IF
                  END DO
               END DO
               IF ( ICLOK.EQ.0 ) THEN
                  ICL = ICLQ(IQ)
                  NMB_CL(ICL) = NMB_CL(ICL) + 1
                  IQ_MBCL(NMB_CL(ICL),ICL) = JQ
                  ICLQ(JQ) = ICL
               ELSE IF ( ICLOK.EQ.1 ) THEN
                  IF ( ICLQ(JQ).NE.ICLQ(IQ) ) THEN
                     WRITE (6,*) 'class(IQ) <> class(JQ)',IQ,ICLQ(IQ),
     &                           JQ,ICLQ(JQ),ISYM
                     NERR = NERR + 1
                  END IF
               ELSE
                  WRITE (6,*) 'JQ found in more than 1 classes ICL'
                  NERR = NERR + 1
               END IF
C
            END IF
         END DO
C
      END DO
C
C-----------------------------------------------------------------------
C  general symmetry operation:   {D|->t} ->q = D ->q + ->t = R_q + ->q''
C     find  R_q  as NSFTSYMQ that ensures ->q'' to be within unit cell
C-----------------------------------------------------------------------
C
      DO IQ = 1,NQ
         DO ISYM = 1,NSYM
            IF ( .NOT.SYMACCEPTED(ISYM) ) CYCLE
C
            CALL DGEMV('N',3,3,1D0,MROTR(1,1,ISYM),3,QVEC(1,IQ),1,0D0,
     &                 QVECP(1,IQ),1)
C
            CALL DAXPY(3,1D0,SYMTVEC(1,ISYM),1,QVECP(1,IQ),1)
C
            DO I = 1,3
               BV(I) = DDOT(3,BR(1,I),1,QVECP(1,IQ),1)
            END DO
C
            CALL DGEMV('N',3,3,1D0,BRINV,3,BV,1,0D0,CF0,1)
            CF(1:3) = CF0(1:3)
C
            QVECPP(1:3) = 0D0
            QSFT(1:3) = 0D0
            DO I = 1,3
               CF(I) = CF(I) - INT(CF(I)+1000D0) + 1000D0
               IF ( ABS(CF(I)-1D0).LT.TOL ) CF(I) = 0D0
               SHIFT = CF0(I) - CF(I)
               NSFTSYMQ(I,ISYM,IQ) = NINT(SHIFT)
               CALL DAXPY(3,CF(I),BR(1,I),1,QVECPP,1)
               CALL DAXPY(3,SHIFT,BR(1,I),1,QSFT,1)
            END DO
C
            QVECTST(1:3) = QSFT(1:3) + QVECPP(1:3)
C
            IF ( .NOT.RVEC_SAME(3,QVECTST,QVECP(1,IQ),1D-8) ) THEN
               WRITE (6,99008) IQ,ISYM,QVECTST,(QVECP(I,IQ),I=1,3)
               NERR = NERR + 1
            END IF
C
            IF ( MOL ) QVECPP(1:3) = QVECP(1:3,IQ)
C
            DO JQ = 1,NQ
               IF ( RVEC_SAME(3,QVECPP,QVEC(1,JQ),1D-8) ) THEN
                  IF ( JQ.NE.IQPSYMQ(ISYM,IQ) .OR. IQORGQP(ISYM,JQ)
     &                 .NE.IQ ) THEN
                     WRITE (6,99009) IQ,ISYM,JQ,IQPSYMQ(ISYM,IQ),
     &                               IQORGQP(ISYM,JQ)
                     NERR = NERR + 1
                  END IF
               END IF
            END DO
C
         END DO
      END DO
C
      IF ( NERR.GT.0 ) CALL STOP_MESSAGE(ROUTINE,'NERR > 0 (A)')
C
C-----------------------------------------------------------------------
C              readjust number of inequivalent atom types  NT
C-----------------------------------------------------------------------
      NT0 = NT
C
      DO IT = 1,NT0
         IT0(IT) = IT
         NAT0(IT) = NAT(IT)
         DO IA = 1,NAT0(IT)
            IQAT0(IA,IT) = IQAT(IA,IT)
         END DO
      END DO
C
      DO IT = 1,NT0
C-----------------------------------------------------------------------
C             the first lattice site IQAT(1,IT) with symmetry class ICL1
C                                            remains occupied by type IT
C
         ICL1 = ICLQ(IQAT(1,IT))
         NAT(IT) = 1
         IQAT(1,IT) = IQAT0(1,IT)
         NTADD = 0
         DO IA = 2,NAT0(IT)
C-----------------------------------------------------------------------
C       check all other sites  IQAT(IA,IT) whether they are still in the
C          same symmetry class ICL1; i.e. occupied by the same atom type
C
            ICLA = ICLQ(IQAT(IA,IT))
C
CC            write(6,'(A,5I3,2X,'' IAQT'',3I3)') 'it,nat,ia,icl1,icl2',
CC     *                 it,NAT0(IT),ia,icl1,icla,IQAT(1,IT),IQAT(IA,IT)
C
            IF ( ICL1.EQ.ICLA ) THEN
               NAT(IT) = NAT(IT) + 1
               IQAT(NAT(IT),IT) = IQAT0(IA,IT)
            ELSE
C--------------------------------------- the symmetry class is different
               FOUND = .FALSE.
               DO JT = (NT+1),(NT+NTADD)
                  JCL = ICLQ(IQAT(1,JT))
C ---------------------------------- the symmetry class ICLA is the same
C------------------------------- as for a previously added new atom type
                  IF ( ICLA.EQ.JCL ) THEN
                     NAT(JT) = NAT(JT) + 1
                     IQAT(NAT(JT),JT) = IQAT0(IA,IT)
                     FOUND = .TRUE.
                  END IF
               END DO
C --------------------------------- the symmetry class ICLA is not found
C -------------------------------------- add a new atom type to the list
               IF ( .NOT.FOUND ) THEN
                  NTADD = NTADD + 1
                  NTP = NT + NTADD
                  IF ( NTP.GT.NTMAX ) THEN
                     WRITE (6,99022) IT,NT,NTP,NTMAX
                     CALL STOP_MESSAGE(ROUTINE,'trouble checking types')
                  END IF
                  NAT(NTP) = 1
                  IQAT(1,NTP) = IQAT0(IA,IT)
                  IT0(NTP) = IT
               END IF
            END IF
         END DO
         NT = NT + NTADD
      END DO
C
C=======================================================================
C                find rotations of BZ-wegdes
C=======================================================================
      IF ( .NOT.MOL ) THEN
C
         WEDGEOK(1:NSYM) = SYMACCEPTED(1:NSYM)
C
         DO IROT = 1,NSYM
            IF ( SYMCRYSYS(IROT) ) THEN
               IF ( SYMUNITARY(IROT) ) THEN
                  MROTK(1:3,1:3,IROT) = MROTR(1:3,1:3,IROT)
               ELSE
                  MROTK(1:3,1:3,IROT) = -MROTR(1:3,1:3,IROT)
               END IF
            END IF
         END DO
C
         NWEDGE = 1
         IWEDGEROT(1) = 1
         DO IROT = 1,NSYM
            IF ( SYMCRYSYS(IROT) ) THEN
               IF ( .NOT.WEDGEOK(IROT) ) THEN
                  DO ISYM = 1,NSYM
                     IF ( SYMACCEPTED(ISYM) ) THEN
C
                        CALL DGEMM('N','N',3,3,3,1D0,MROTK(1,1,ISYM),3,
     &                             MROTK(1,1,IROT),3,0D0,RMAT3,3)
C
                        DO JSYM = 1,NSYM
                           IF ( RVEC_SAME(9,RMAT3,MROTK(1,1,JSYM),1D-8)
     &                          ) WEDGEOK(JSYM) = .TRUE.
                        END DO
                     END IF
                  END DO
                  NWEDGE = NWEDGE + 1
                  IWEDGEROT(NWEDGE) = IROT
                  WEDGEOK(IROT) = .TRUE.
               END IF
            END IF
         END DO
C-----------------------------------------------------------------------
C                                                     check the matrices
         NERR = 0
         RMAT3(1:3,1:3) = 0D0
         RMAT3(1,1) = 1D0
         RMAT3(2,2) = 1D0
         RMAT3(3,3) = 1D0
C
         IF ( .NOT.RVEC_SAME(9,RMAT3,MROTK(1,1,1),1D-8) ) THEN
            NERR = NERR + 1
            WRITE (6,*) 
     &                '<SYMLATTICE>: 1st symmetry operation should be E'
         END IF
C
         DO IROT = 1,NSYM
            IF ( SYMCRYSYS(IROT) ) THEN
               IF ( .NOT.WEDGEOK(IROT) ) THEN
                  WRITE (6,*) '<SYMLATTICE>: BZ-wedge ',IROT,' not set'
                  NERR = NERR + 1
               END IF
            END IF
         END DO
C
         IF ( IPRINT.GE.0 ) WRITE (6,99023) NWEDGE,
     &                             (IWEDGEROT(I),SYMSYMBL(IWEDGEROT(I)),
     &                             I=1,NWEDGE)
C
         IF ( NSYMCRYSYS.NE.NWEDGE*NSYMACCEPTED ) THEN
            WRITE (6,*) 
     &               '<SYMLATTICE>: NSYMCRYSYS <> NWEDGE * NSYMACCEPTED'
            NERR = NERR + 1
         END IF
C
      END IF
C
C=======================================================================
C                   update or set up tables
C=======================================================================
C
C---------------------------------------------------- IQAT --> ITOQ, NOQ
C
      IF ( NT.GT.NTMAX ) THEN
         WRITE (6,99003) NT0,NT,NTMAX
         NERR = NERR + 1
      ELSE
         DO IQ = 1,NQ
            NOQ0(IQ) = NOQ(IQ)
            NOQ(IQ) = 0
         END DO
         DO IT = 1,NT
            DO IA = 1,NAT(IT)
               IQ = IQAT(IA,IT)
               NOQ(IQ) = NOQ(IQ) + 1
               ITOQ(NOQ(IQ),IQ) = IT
            END DO
         END DO
         DO IQ = 1,NQ
            IF ( NOQ0(IQ).NE.NOQ(IQ) ) THEN
               WRITE (6,99004) IQ,NOQ0(IQ),NOQ(IQ)
               NERR = NERR + 1
            END IF
         END DO
      END IF
C
C=======================================================================
C                         ISYMGENQ  IQREPQ
C=======================================================================
C
C  set up the symmetry connection between sites
C  for each site  IQP   IQ = IQREPQ(IQP) gives the representative site
C  with the symmetry operation  ISYMGENQ   transfering IQ to IQP
C  ISYM is not unique as there may be several operations with q' = S q
C  to find ISYM a seach is first done within the UNITARY operations
C
      DO IQ = 1,NQ
         ISYMGENQ(IQ) = 0
      END DO
C
      DO IQ = 1,NQ
C
C----------------------------------------------- scan UNITARY operations
         DO ISYM = 1,NSYM
            IF ( SYMACCEPTED(ISYM) .AND. SYMUNITARY(ISYM) ) THEN
               IQP = IQPSYMQ(ISYM,IQ)
               IF ( ISYMGENQ(IQP).EQ.0 ) THEN
                  ISYMGENQ(IQP) = ISYM
                  IQREPQ(IQP) = IQ
               END IF
            END IF
         END DO
C
C------------------------------------------ scan ANTI-UNITARY operations
         DO ISYM = 1,NSYM
            IF ( SYMACCEPTED(ISYM) .AND. .NOT.SYMUNITARY(ISYM) ) THEN
               IQP = IQPSYMQ(ISYM,IQ)
               IF ( ISYMGENQ(IQP).EQ.0 ) THEN
                  ISYMGENQ(IQP) = ISYM
                  IQREPQ(IQP) = IQ
               END IF
            END IF
         END DO
C
      END DO
C
C=======================================================================
C                         NQEQ  IQEQ
C=======================================================================
C
      DO IQ = 1,NQ
         NQEQ(IQ) = 0
         DO JQ = 1,NQ
            IF ( NOQ(IQ).EQ.NOQ(JQ) ) THEN
               IFLAG = 1
               DO IO = 1,NOQ(IQ)
                  IF ( ITOQ(IO,IQ).NE.ITOQ(IO,JQ) ) IFLAG = 0
               END DO
               IF ( IFLAG.EQ.1 ) THEN
                  NQEQ(IQ) = NQEQ(IQ) + 1
                  IQEQ(NQEQ(IQ),IQ) = JQ
               END IF
            END IF
         END DO
      END DO
C
      IF ( IPRINT.GT.0 ) WRITE (6,99014)
      DO IQ = 1,NQ
         IQREP = IQREPQ(IQ)
         IFLAG = 1
         DO IEQ = 1,NQEQ(IQ)
            JQ = IQEQ(IEQ,IQ)
            IF ( IQREP.EQ.JQ ) IFLAG = 0
         END DO
         IF ( IFLAG.EQ.1 ) THEN
            WRITE (6,99005) IQ,IQREP,(IQEQ(IEQ,IQ),IEQ=1,NQEQ(IQ))
            NERR = NERR + 1
         END IF
         IF ( IPRINT.GT.0 ) WRITE (6,99015) IQ,IQREPQ(IQ),ISYMGENQ(IQ),
     &                             SYMUNITARY(ISYMGENQ(IQ)),NQEQ(IQ),
     &                             (IQEQ(IEQ,IQ),IEQ=1,NQEQ(IQ))
      END DO
C
C------------------------------------------------------ NMSYM  IQREPMSYM
C
      NMSYM = 0
      DO IQ = 1,NQ
         IQREP = IQREPQ(IQ)
         IF ( IQREP.EQ.IQ ) THEN
            NMSYM = NMSYM + 1
            IQREPMSYM(NMSYM) = IQ
         END IF
      END DO
C
C------------------------------------------------------ transfer NCL etc
C
      IF ( IWR.GT.0 ) THEN
         NMBMAX = 0
         DO ICL = 1,NCL
            NMBMAX = MAX(NMBMAX,NMB_CL(ICL))
         END DO
         WRITE (IWR) NCL,NMBMAX
         DO ICL = 1,NCL
            WRITE (IWR) NMB_CL(ICL)
            WRITE (IWR) (IQ_MBCL(IMB,ICL),IMB=1,NMB_CL(ICL))
         END DO
      END IF
C
      IF ( NERR.GT.0 ) CALL STOP_MESSAGE(ROUTINE,'NERR > 0 (B)')
C
C=======================================================================
C              check symmetry connectivity among sites
C=======================================================================
C
      NERR = 0
      DO IQP = 1,NQ
C
         ISYM = ISYMGENQ(IQP)
C
         IQ = IQORGQP(ISYM,IQP)
C
         IF ( IQ.NE.IQREPQ(IQP) ) THEN
            NERR = NERR + 1
            WRITE (6,99024) IQP,ISYMGENQ(IQP),IQ,IQREPQ(IQP)
         END IF
C
         DO ISYM = 1,NSYM
            IF ( SYMACCEPTED(ISYM) ) THEN
C
               IQ = IQORGQP(ISYM,IQP)
C
               IF ( IQP.NE.IQPSYMQ(ISYM,IQ) ) THEN
                  NERR = NERR + 1
                  WRITE (6,99025) IQP,ISYM,IQ,IQPSYMQ(ISYM,IQ)
               END IF
C
            END IF
         END DO
C
      END DO
C
      IF ( NERR.GT.0 ) CALL STOP_MESSAGE(ROUTINE,
     &     'trouble checking the symmetry connectivity among sites')
C
      DEALLOCATE (ICLQ,IQ_MBCL,NMB_CL,QVEC,ICLQ0,IQAT0,MVECQ,QVECP,NOQ0)
      DEALLOCATE (MVECQP,NAT0,IQREPCL0,DSYM_CLM)
C
99001 FORMAT (/,' ##### TROUBLE in <SYMLATTICE> ',49('#'),:,/,10X,
     &        'ISYM =',I3,'  no inverse operation found ')
99002 FORMAT (/,10X,'the detected ',I3,
     &        ' symmetry operations form a group',/)
99003 FORMAT (/,' ##### TROUBLE in <SYMLATTICE> ',49('#'),:,/,10X,
     &        'NT0=',I3,3X,'NT=',I3,' > NTMAX=',I3,
     &        '   ITOQ cannot be updated')
99004 FORMAT (/,' ##### TROUBLE in <SYMLATTICE> ',49('#'),:,/,10X,
     &        'number of occupants changed for IQ=',I3,'  NOQ=',I3,
     &        ' --> ',I3)
99005 FORMAT (/,' ##### TROUBLE in <SYMLATTICE> ',49('#'),:,/,10X,
     &        'for IQ=',I3,'   IQREP=',I3,
     &        ' not in list of equivalent sites',(/,10X,10I3))
99006 FORMAT (/,10X,'test of spatial rotation matrices'//,10X,
     &        'list for all accepted symmetry operations ',//,10X,
     &       'ISYM       index I within list of all possible operations'
     &       ,/,10X,
     &       'SYMEULANG  Euler angles for ACTIVE TEMPORARY convention',
     &       /,10X,'           in degrees - printed as rounded integer',
     &       /,10X,
     &       'SYMSYMBL   symmetry operation (see Bradley Cracknell)',/,
     &       10X,'MROTR      3x3 rotation matrix for cartesian vectors',
     &       /,10X,'VP',9X,
     &       'result of  MROT * V  for vector ->V = (1,2,3)^T',/,10X,
     &       'MROTC      3x3 rotation matrix w.r.t.',
     &       ' crystallographic units',//,7X,'I',4X,
     &       'Euler angles  SYMBOL',14X,'MROTR',12X,'->VP',18X,'MROTC',
     &       /)
99007 FORMAT (5X,I3,1X,3I5,3X,A,2X,3F8.2,F10.2,5X,3F8.2,/,
     &        2(33X,3F8.2,F10.2,5X,3F8.2,/))
99008 FORMAT ('##### IQ=',I4,' ISYM=',I3,' QVECTST = ',3F10.6,/22X,
     &        ' QVECP   = ',3F10.6,/)
99009 FORMAT ('##### IQ=',I4,' ISYM=',I3,3X,' JQ=',I3,3X,
     &        'IQPSYMQ(ISYM,IQ)=',I3,3X,'IQORGQP(ISYM,JQ)=',I3,/)
99010 FORMAT (//,1X,79('*'),/,34X,'<SYMLATTICE>',/,1X,79('*'),//,10X,
     &        'find out symmetry properties of ',A,/)
99011 FORMAT (//,10X,
     &        '* use the Euler angles for ACTIVE TEMPORARY convention ',
     &        '(Rose) to fix ',/,10X,
     &        '  the rotation matrices D for all accepted symmetry ',
     &        'operations',/,10X,'  (',I2,'  out of ',I2,
     &        ' possible ones)',
     &        ' w.r.t spherical coordinates up to l=1',/,10X,
     &        '* transform D to cartesian coordinates ',
     &        'to get 3x3-matrix MROTR')
99012 FORMAT (//,'ISYM=',I2,'  (x y z)  >> (x'' y'' z'')  ',
     &        ' Euler-angles=',3F6.1)
99013 FORMAT (/,10X,'crystal system accepted symmetry operations',I8,/,
     &        10X,'accepted symmetry operations     ',I18,/)
99014 FORMAT (//,10X,'list for all sites  IQ ',//,10X,
     &        'IQREPQ     representative site related by ',
     &        'symmetry operation ISYMGENQ',/,10X,
     &     'ISYMGENQ   symmetry operation that transfers site IQREPQ to'
     &     ,' site IQ',/,10X,
     &     'UNITARY    T indicates symmetry operation ISYMGENQ',
     &     ' to be unitary ',/,10X,'NQEQ       number of equivalent,',
     &     ' i.e. symmetry related sites',/,10X,
     &     'IQEQ       list of NQEQ equivalent,',
     &     ' i.e. symmetry related sites',//,10X,
     &     'IQ  IQREPQ ISYMGENQ UNITARY NQEQ IQEQ')
99015 FORMAT (I12,I7,I6,L9,I6,1X,8I4,:,/,(41X,8I4))
99016 FORMAT (10X,'list for all allowed symmetry operations ',//,10X,
     &       'ISYM       index I within list of all possible operations'
     &       ,/,10X,
     &       'SYMEULANG  Euler angles for ACTIVE TEMPORARY convention',
     &       /,10X,'           in degrees - printed as rounded integer',
     &       /,10X,'INV        I indicates: inversion used to create',
     &       ' symmetry operation',/,10X,
     &       'SYMSYMBL   symmetry operation (see Bradley Cracknell)',/,
     &       10X,'TIM        T indicates: time reversal included in',
     &       ' symmetry operation',/,10X,
     &       'SYMTVEC    primitive translation part of',
     &       ' symmetry operation',/,10X,
     &       'IQPSYMQ    symmetry operation moves site Q to site Q'' '/,
     &       10X,'IQORGQP    site Q'' was originally at Qorg  -',
     &       '  inverse of IQPSYMQ',/,10X,
     &       'SYMDET     determinant of   RMOTR',/)
99017 FORMAT (/,7X,'I',4X,'Euler angles','  INV SYMBOL TIM',5X,
     &        'shift  ->DQ',5X,'Q'' ')
99018 FORMAT (5X,I3,1X,3I5,3X,A,3X,A,3X,A,2X,3F6.2,(5I4))
99019 FORMAT ((59X,5I4))
99020 FORMAT (49X,'IQORGQP   ',5I4,:,/,(59X,5I4))
99021 FORMAT (49X,'SYMDET    ',I4,/,49X,'->p - R->p',3F6.2,/)
99022 FORMAT (/,' ##### TROUBLE in <SYMLATTICE> ',49('#'),/,10X,
     &        'trouble checking the types for IT=',I2,/,10X,'NT=',I3,
     &        '   NTP=',I3,'   NTMAX=',I3,/,10X,
     &        'increase array size   NTMAX')
99023 FORMAT (/,10X,'BZ-wedges to be treated  NWEDGE ',I3,/,10X,
     &        'created by operations  ',/,(10X,6(I2,2X,A4,3X)))
99024 FORMAT (10X,'A:  IQP =',I4,'    ISYMGENQ =',I3,'    IQORGQP =',I4,
     &        '    IQREPQ  =',I4)
99025 FORMAT (10X,'B:  IQP =',I4,'    ISYM     =',I3,'    IQORGQP =',I4,
     &        '    IQPSYMQ =',I4)
99026 FORMAT (//,1X,79('I'),/,10X,
     &        'basis vectors not conform with max. number of',
     &        ' symmetry operations',/,10X,'the full set of ',
     &        'possible symmetry operations will be scanned:',/,10X,
     &        'NSYMCRYSYS0 =',I3,'    >>>    NSYMCRYSYS =',I3,/,1X,
     &        79('I'),/)
      END
C*==xxsym_get_nsftsymq.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE XXSYM_GET_NSFTSYMQ(IQ,ISYM,NERR,MOL,NQ,QVEC,QVECP,BR,
     &                              BRINV,MROTR,SYMTVEC,IQORGQP,IQPSYMQ,
     &                              NSFTSYMQ,NQMAX)
C   ********************************************************************
C   *                                                                  *
C   *  general symmetry operation:                                     *
C   *                                                                  *
C   *        {D|->t} ->q = D ->q + ->t = R_q + ->q''                   *
C   *                                                                  *
C   * find  R_q  as NSFTSYMQ that ensures ->q'' to be within unit cell *
C   *                                                                  *
C   * NOTE: the error flag NERR has to be initialized                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SYMMETRY,ONLY:NSYMMAX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SYM_GET_NSFTSYMQ')
      REAL*8 TOL
      PARAMETER (TOL=1.0D-6)
C
C Dummy arguments
C
      INTEGER IQ,ISYM,NERR,NQ,NQMAX
      LOGICAL MOL
      REAL*8 BR(3,3),BRINV(3,3),MROTR(3,3,NSYMMAX),QVEC(3,NQMAX),
     &       QVECP(3,NQMAX),SYMTVEC(3,NSYMMAX)
      INTEGER IQORGQP(NSYMMAX,NQMAX),IQPSYMQ(NSYMMAX,NQMAX),
     &        NSFTSYMQ(3,NSYMMAX,NQMAX)
C
C Local variables
C
      REAL*8 BV(3),CF(3),CF0(3),QSFT(3),QVECPP(3),QVECTST(3),SHIFT
      REAL*8 DDOT
      INTEGER I,JQ,N_JQ
      LOGICAL RVEC_SAME
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
      N_JQ = 0
C
      CALL DGEMV('N',3,3,1D0,MROTR(1,1,ISYM),3,QVEC(1,IQ),1,0D0,
     &           QVECP(1,IQ),1)
C
      CALL DAXPY(3,1D0,SYMTVEC(1,ISYM),1,QVECP(1,IQ),1)
C
      DO I = 1,3
         BV(I) = DDOT(3,BR(1,I),1,QVECP(1,IQ),1)
      END DO
C
      CALL DGEMV('N',3,3,1D0,BRINV,3,BV,1,0D0,CF0,1)
      CF(1:3) = CF0(1:3)
C
      QVECPP(1:3) = 0D0
      QSFT(1:3) = 0D0
      DO I = 1,3
         CF(I) = CF(I) - INT(CF(I)+1000D0) + 1000D0
         IF ( ABS(CF(I)-1D0).LT.TOL ) CF(I) = 0D0
         SHIFT = CF0(I) - CF(I)
         NSFTSYMQ(I,ISYM,IQ) = NINT(SHIFT)
         CALL DAXPY(3,CF(I),BR(1,I),1,QVECPP,1)
         CALL DAXPY(3,SHIFT,BR(1,I),1,QSFT,1)
      END DO
C
      QVECTST(1:3) = QSFT(1:3) + QVECPP(1:3)
C
      IF ( .NOT.RVEC_SAME(3,QVECTST,QVECP(1,IQ),1D-8) ) THEN
         WRITE (6,99001) IQ,ISYM,QVECTST,(QVECP(I,IQ),I=1,3)
         NERR = NERR + 1
      END IF
C
      IF ( MOL ) QVECPP(1:3) = QVECP(1:3,IQ)
C
      DO JQ = 1,NQ
         IF ( RVEC_SAME(3,QVECPP,QVEC(1,JQ),1D-8) ) THEN
            N_JQ = N_JQ + 1
            IF ( JQ.NE.IQPSYMQ(ISYM,IQ) .OR. IQORGQP(ISYM,JQ).NE.IQ )
     &           THEN
               WRITE (6,99002) IQ,ISYM,JQ,IQPSYMQ(ISYM,IQ),
     &                         IQORGQP(ISYM,JQ)
               NERR = NERR + 1
            END IF
         END IF
      END DO
C
      IF ( N_JQ.NE.1 ) CALL STOP_MESSAGE(ROUTINE,'N_JQ .NE. 1 ')
C
99001 FORMAT ('##### IQ=',I4,' ISYM=',I3,' QVECTST = ',3F10.6,/22X,
     &        ' QVECP   = ',3F10.6,/)
99002 FORMAT ('##### IQ=',I4,' ISYM=',I3,3X,' JQ=',I3,3X,
     &        'IQPSYMQ(ISYM,IQ)=',I3,3X,'IQORGQP(ISYM,JQ)=',I3,/)
      END
