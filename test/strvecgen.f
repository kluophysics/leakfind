C*==strqqplim.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
c     SUBROUTINE STRQQPLIM(NQ_STR,IPRINT,RMAX,GMAX,NQQP_STR,NIJQ,IJQ,
      SUBROUTINE STRQQPLIM(NQ_STR,IPRINT,RMAX,GMAX,NQQP_STR,
     &                     NUMRH,NUMGH,RA,GA,ABAS,BBAS,QBAS,QQPX,QQPY,
     &                     QQPZ,NQMAX,NQQP_STRMAX)
C   ********************************************************************
C   *                                                                  *
C   *           find set of inequivelent lattice pairs Q-QP            *
C   *             and range for R- and K-lattice vectors               *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_STR,ONLY:NQQP_STR_RED
      USE MOD_STR,ONLY: IJQ, NIJQ
      IMPLICIT NONE
C*--STRQQPLIM14
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='STRQQPLIM')
C
C Dummy arguments
C
      REAL*8 GA,GMAX,RA,RMAX
      INTEGER IPRINT,NQMAX,NQQP_STR,NQQP_STRMAX,NQ_STR,NUMGH,NUMRH
      REAL*8 ABAS(3,3),BBAS(3,3),QBAS(3,NQMAX),QQPX(NQQP_STRMAX),
     &       QQPY(NQQP_STRMAX),QQPZ(NQQP_STRMAX)
c     INTEGER IJQ(NQQP_STRMAX,NQQP_STRMAX),NIJQ(NQQP_STRMAX)
C
C Local variables
C
      REAL*8 DD(3),DGMIN,DK(3),DQ,DQIJ,DQIJMAXNEW,DQIJMAXORG,DQIJMAXTMP,
     &       DQIPJMAX,DRMIN,GMAX1,QIJVEC(3),QIPVEC(1:3),QQPX_TMP(:),
     &       QQPY_TMP(:),QQPZ_TMP(:),RMAX1
      REAL*8 DNRM2
      INTEGER I,I1,I123TOP,I1MIN,I2,I2MIN,I3,I3MIN,IFLAG,IJQ_TMP(:,:),
     &        IQ,IQQP,J,JQ,M,N,NIJQ_TMP(:),NOFF,NUMG,NUMR
      INTEGER NIJQMAX
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE QQPX_TMP,QQPY_TMP,QQPZ_TMP,IJQ_TMP,NIJQ_TMP
c     ALLOCATABLE IJQ,NIJQ
C
C-----------------------------------------------------------------------
C modify the list of basis vectors to minimize the interatomic distance
C-----------------------------------------------------------------------
C
      IF ( NQ_STR.GT.1000 ) THEN
         DQIJMAXORG = 0D0
         DO IQ = 1,NQ_STR
            DO JQ = 2,NQ_STR
               QIJVEC(1:3) = QBAS(1:3,IQ) - QBAS(1:3,JQ)
               DQIJ = DNRM2(3,QIJVEC,1)
               DQIJMAXORG = MAX(DQIJ,DQIJMAXORG)
            END DO
         END DO
C
         I123TOP = 2
         DO IQ = 2,NQ_STR
C
            IFLAG = 0
            DQIJMAXTMP = 1D+20
C
            DO I1 = I123TOP, - I123TOP, - 1
               DO I2 = I123TOP, - I123TOP, - 1
                  DO I3 = I123TOP, - I123TOP, - 1
C
                     QIPVEC(1:3) = QBAS(1:3,IQ) + I1*ABAS(1:3,1)
     &                             + I2*ABAS(1:3,2) + I3*ABAS(1:3,3)
C
                     DQIPJMAX = 0D0
                     DO JQ = 1,IQ - 1
                        QIJVEC(1:3) = QIPVEC(1:3) - QBAS(1:3,JQ)
                        DQIJ = DNRM2(3,QIJVEC,1)
                        DQIPJMAX = MAX(DQIJ,DQIPJMAX)
                     END DO
C
                     IF ( DQIPJMAX.LT.DQIJMAXTMP ) THEN
                        DQIJMAXTMP = DQIPJMAX
                        I1MIN = I1
                        I2MIN = I2
                        I3MIN = I3
                        IFLAG = 1
                     END IF
C
                  END DO
               END DO
            END DO
C
            IF ( IFLAG.EQ.1 ) QBAS(1:3,IQ) = QBAS(1:3,IQ)
     &           + I1MIN*ABAS(1:3,1) + I2MIN*ABAS(1:3,2)
     &           + I3MIN*ABAS(1:3,3)
         END DO
C
         DQIJMAXNEW = 0D0
         DO IQ = 1,NQ_STR
            DO JQ = 2,NQ_STR
               QIJVEC(1:3) = QBAS(1:3,IQ) - QBAS(1:3,JQ)
               DQIJ = DNRM2(3,QIJVEC,1)
               DQIJMAXNEW = MAX(DQIJ,DQIJMAXNEW)
            END DO
         END DO
C
         WRITE (6,'(/,10X,''basis vectors in units of a''/)')
         DO IQ = 1,NQ_STR
            WRITE (6,99001) (QBAS(I,IQ),I=1,3)
         END DO
         WRITE (*,*) ' dqijmaxorg ',DQIJMAXORG
         WRITE (*,*) ' dqijmaxnew ',DQIJMAXNEW
C
      END IF
C
C   calculate radii RA,GA of spheres holding all vectors used in lattice
C   sums. RMAX1 is longest basis vector. GMAX1 is 2* longest vector
C   in Brillouin zone. must be reconsidered in any new applications
C
C-----------------------------------------------------------------------
C      create table of vectors    ->QQP = ->q(IQ) - ->q(JQ)
C      the first entry is         ->q(IQ) - ->q(IQ) = ->0
C      the first half of the table lists all cases derived for IQ <= JQ
C      the second half lists all vectors derived from IQ > JQ
C      (IQ,JQ) and (JQ,IQ) have their direction just reversed
C-----------------------------------------------------------------------
      ALLOCATE (QQPX_TMP(NQQP_STRMAX),QQPY_TMP(NQQP_STRMAX))
      ALLOCATE (QQPZ_TMP(NQQP_STRMAX),NIJQ_TMP(NQQP_STRMAX))
      ALLOCATE (IJQ_TMP(NQQP_STRMAX,NQQP_STRMAX))
C
      RMAX1 = 0.0D0
      NQQP_STR = 0
      NIJQ_TMP(1) = 1
      QQPX_TMP(1) = 0D0
      QQPY_TMP(1) = 0D0
      QQPZ_TMP(1) = 0D0
C
C---------------------------------- deal with all pairs of sites IQ - JQ
      DO IQ = 1,NQ_STR
         DO JQ = 1,NQ_STR
            DO IQQP = 1,NQQP_STR
               IF ( (ABS(QQPX_TMP(IQQP)-(QBAS(1,IQ)-QBAS(1,JQ)))
     &              .LT.1D-5) .AND. 
     &              (ABS(QQPY_TMP(IQQP)-(QBAS(2,IQ)-QBAS(2,JQ)))
     &              .LT.1D-5) .AND. 
     &              (ABS(QQPZ_TMP(IQQP)-(QBAS(3,IQ)-QBAS(3,JQ)))
     &              .LT.1D-5) ) THEN
                  NIJQ_TMP(IQQP) = NIJQ_TMP(IQQP) + 1
c modified by XJQ: allow a supercell containing > 99 atoms
c                  IJQ_TMP(NIJQ_TMP(IQQP),IQQP) = 100*IQ + JQ
                  IJQ_TMP(NIJQ_TMP(IQQP),IQQP) = max(100,NQ_STR+1)*IQ 
     &                                           + JQ
c end-mod-xjq
                  GOTO 50
               END IF
            END DO
C
            IF ( (NQQP_STR+1).GT.NQQP_STRMAX ) THEN
               WRITE (6,99002) NQQP_STR
               CALL STOP_MESSAGE(ROUTINE,'NQQP_STR ???')
            END IF
            NQQP_STR = NQQP_STR + 1
            NIJQ_TMP(NQQP_STR) = 1
c modified by XJQ: allow a supercell containing > 99 atoms
c            IJQ_TMP(1,NQQP_STR) = 100*IQ + JQ
            IJQ_TMP(1,NQQP_STR) = max(100,NQ_STR+1)*IQ + JQ
c end-mod-xjq
            QQPX_TMP(NQQP_STR) = QBAS(1,IQ) - QBAS(1,JQ)
            QQPY_TMP(NQQP_STR) = QBAS(2,IQ) - QBAS(2,JQ)
            QQPZ_TMP(NQQP_STR) = QBAS(3,IQ) - QBAS(3,JQ)
            DQ = SQRT(QQPX_TMP(NQQP_STR)**2+QQPY_TMP(NQQP_STR)
     &           **2+QQPZ_TMP(NQQP_STR)**2)
            IF ( DQ.GE.RMAX1 ) RMAX1 = DQ
C
 50      END DO
      END DO
C
      NQQP_STR_RED = NQQP_STR/2 + 1
      IF ( 2*NQQP_STR_RED-1.NE.NQQP_STR )
     &      CALL STOP_MESSAGE(ROUTINE,'NQQP_STR ???')
C
C-------------------------------------------------------------- now sort
C
      NIJQMAX = MAXVAL(NIJQ_TMP, NQQP_STR)
      ALLOCATE (NIJQ(NQQP_STR),IJQ(NIJQMAX,NQQP_STR))



      N = 1
      IQQP = 1
      NIJQ(N) = NIJQ_TMP(IQQP)
      IJQ(1:NIJQ(N),N) = IJQ_TMP(1:NIJQ(IQQP),IQQP)
      QQPX(N) = QQPX_TMP(IQQP)
      QQPY(N) = QQPY_TMP(IQQP)
      QQPZ(N) = QQPZ_TMP(IQQP)
C
      NOFF = NQQP_STR/2
      IF ( NQQP_STR.NE.(1+2*NOFF) )
     &     CALL STOP_MESSAGE(ROUTINE,'NOFF ???')
C
      DO IQQP = 2,NQQP_STR
C
         DO J = 2,N
            IF ( (ABS(QQPX_TMP(IQQP)+QQPX(J)).LT.1D-5) .AND. 
     &           (ABS(QQPY_TMP(IQQP)+QQPY(J)).LT.1D-5) .AND. 
     &           (ABS(QQPZ_TMP(IQQP)+QQPZ(J)).LT.1D-5) ) THEN
               M = J + NOFF
               NIJQ(M) = NIJQ_TMP(IQQP)
               IJQ(1:NIJQ(M),M) = IJQ_TMP(1:NIJQ_TMP(IQQP),IQQP)
               QQPX(M) = QQPX_TMP(IQQP)
               QQPY(M) = QQPY_TMP(IQQP)
               QQPZ(M) = QQPZ_TMP(IQQP)
               GOTO 100
            END IF
         END DO
C
         N = N + 1
         NIJQ(N) = NIJQ_TMP(IQQP)
         IJQ(1:NIJQ(N),N) = IJQ_TMP(1:NIJQ_TMP(IQQP),IQQP)
         QQPX(N) = QQPX_TMP(IQQP)
         QQPY(N) = QQPY_TMP(IQQP)
         QQPZ(N) = QQPZ_TMP(IQQP)
C
 100  END DO
      IF ( N.NE.(1+NOFF) ) CALL STOP_MESSAGE(ROUTINE,'N.NE.(NOFF+1)')

      DEALLOCATE(IJQ_TMP, NIJQ_TMP)

C
      IFLAG = 0
      DO I = 2,1 + NOFF
         J = I + NOFF
         IF ( (ABS(QQPX(I)+QQPX(J)).GT.1D-5) .OR. 
     &        (ABS(QQPY(I)+QQPY(J)).GT.1D-5) .OR. 
     &        (ABS(QQPZ(I)+QQPZ(J)).GT.1D-5) .OR. (NIJQ(I).NE.NIJQ(J)) )
     &        IFLAG = 1
      END DO
C
      WRITE (6,99003) NQQP_STR,NQQP_STR_RED
C
      IF ( IFLAG.NE.0 .OR. IPRINT.GT.0 ) THEN
         WRITE (6,99004)
         DO N = 1,NQQP_STR
            WRITE (6,99005) NIJQ(N),QQPX(N),QQPY(N),QQPZ(N),
c     &                      (INT(IJQ(M,N)/100),
     &                      (INT(IJQ(M,N)/max(100,NQ_STR+1)),
c     &                      (IJQ(M,N)-100*INT(IJQ(M,N)/100)),M=1,
     &                      (IJQ(M,N)-max(100,NQ_STR+1)*INT(IJQ(M,N)
     &                      /max(100,NQ_STR+1))),M=1,
     &                      MIN(4,NIJQ(N)))
            IF ( NIJQ(N).GT.5 ) WRITE (6,99006)
c     &                                 (INT(IJQ(M,N)/100),(IJQ(M,N)
c     &                                 -100*INT(IJQ(M,N)/100)),M=5,
     &                                 (INT(IJQ(M,N)/max(100,NQ_STR+1)),
     &                                 (IJQ(M,N)-max(100,NQ_STR+1)
     &                                 *INT(IJQ(M,N)
     &                                 /max(100,NQ_STR+1))),
     &                                 M=5,NIJQ(N))
         END DO
      END IF
      IF ( IFLAG.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'IFLAG .NE. 0')
C
      DO I = 1,3
         DD(I) = DNRM2(3,ABAS(1,I),1)
         DK(I) = DNRM2(3,BBAS(1,I),1)
      END DO
      DRMIN = DMIN1(DD(1),DD(2),DD(3))
      DGMIN = DMIN1(DK(1),DK(2),DK(3))
      GMAX1 = DMAX1(DK(1),DK(2),DK(3))
C
      RA = RMAX + RMAX1*1.001D0
      GA = GMAX + GMAX1
      GA = GMAX
C
      WRITE (6,99008) RMAX,GMAX,RMAX1,GMAX1,RA,GA
C
      NUMR = 2*(INT(RA/DRMIN)+1) + 1
      NUMG = 2*(INT(GA/DGMIN)+1) + 1
      NUMRH = NUMR/2 + 1
      NUMGH = NUMG/2 + 1
      WRITE (6,99007) NUMR,NUMG
99001 FORMAT (12X,'(',F10.5,',',F10.5,',',F10.5,' )')
C
99002 FORMAT (/,10X,'STOP in <STRQQPLIM> ',/,'NQQP_STR =',I8,
     &        ' > NQQP_STRMAX')
99003 FORMAT (/,10X,'number of inequivalent QQP-block-matrices',
     &        '  NQQP_STR     =',I8,/,10X,
     &        'reduced number of inequiv. block-matrices',
     &        '  NQQP_STR_RED =',I8,/)
99004 FORMAT (10X,'NIJQ (   ->Q[IQ] -  ->Q[JQ]   ) [ IQ, JQ] ...')
99005 FORMAT (10X,I3,2X,'(',F7.3,',',F7.3,',',F7.3,' )',1X,
     &        4(:,'[',I3,',',I3,']'))
99006 FORMAT ((42X,4(:,'[',I3,',',I3,']')))
99007 FORMAT (10X,'NUMR  =',I10,10X,'NUMG  =',I10)
99008 FORMAT (//,10X,'RMAX  =',F10.4,10X,'GMAX  =',F10.4,/,10X,
     &        'RMAX1 =',F10.4,10X,'GMAX1 =',F10.4,/,10X,'RA    =',F10.4,
     &        10X,'GA    =',F10.4,/)
      END
C*==strvecgen.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE STRVECGEN(IPRINT,ALAT,RMAX,GMAX,GMAXSQ,NRDL,SMAX,NGRL,
     &                     NQQP_STR,NUMGH,NUMRH,GA,RA,ABAS,BBAS,QQPX,
     &                     QQPY,QQPZ,G1,G2,G3,R1,R2,R3,NGRLMAX,NRDLMAX,
     &                     NQQP_STRMAX)
C   ********************************************************************
C   *                                                                  *
C   *   GENERATE VECTORS OF DIRECT AND RECIPROCAL SPACE FROM           *
C   *   BASIC TRANSLATION VECTORS  ABAS                                *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C*--STRVECGEN302
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='STRVECGEN')
C
C Dummy arguments
C
      REAL*8 ALAT,GA,GMAX,GMAXSQ,RA,RMAX
      INTEGER IPRINT,NGRL,NGRLMAX,NQQP_STR,NQQP_STRMAX,NRDL,NRDLMAX,
     &        NUMGH,NUMRH
      REAL*8 ABAS(3,3),BBAS(3,3),QQPX(NQQP_STRMAX),QQPY(NQQP_STRMAX),
     &       QQPZ(NQQP_STRMAX)
      INTEGER G1(NGRLMAX),G2(NGRLMAX),G3(NGRLMAX),R1(NRDLMAX),
     &        R2(NRDLMAX),R3(NRDLMAX),SMAX(NQQP_STR)
C
C Local variables
C
      REAL*8 DG(:),DR(:),DX,EDUMAX,EDUMIN,GX,GY,GZ,KNX,KNY,KNZ,KSQ,P,
     &       REDU,RX,RY,RZ,SX,SY,SZ
      INTEGER I,I1,I2,I3,IA_ERR,IE,IETOP,IG,II,INSIDE,IQQP,IR,J,J1,J2,
     &        J3,JHIG,JLOW,K,KK,NG,NR,NSG(:),NSH,NSR(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DG,DR,NSG,NSR
C
      ALLOCATE (DG(NGRLMAX),DR(NRDLMAX))
      ALLOCATE (NSG(NGRLMAX),NSR(NRDLMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: NSG')
C
C ======================================================================
C                            REAL SPACE
C ======================================================================
C
C  accept a R-point ->R[S] in list if
C
C                  |->R[S] - (->Q[i] - ->Q[i'])| < rmax
C
C  i.e. the center of the convergence sphere of radius rmax  is at
C  (->Q[i] - ->Q[i']) and is therfore different for every [i,i'].
C  the  SMAX(IQQP)  R-points ->R[S] within the sphere are picked
C  out of the list using the pointer INDR(S,IQQP) S=1,..,SMAX(IQQP)
C  INDR(S,IQQP) is set up in <STRAA> and used in <STRBBDD>.
C
      IF ( IPRINT.GT.2 ) WRITE (6,99001)
C
      DO IQQP = 1,NQQP_STR
         SMAX(IQQP) = 0
      END DO
C
      NR = 0
      DO I1 = -NUMRH, + NUMRH
         DO I2 = -NUMRH, + NUMRH
            DO I3 = -NUMRH, + NUMRH
               RX = I1*ABAS(1,1) + I2*ABAS(1,2) + I3*ABAS(1,3)
               RY = I1*ABAS(2,1) + I2*ABAS(2,2) + I3*ABAS(2,3)
               RZ = I1*ABAS(3,1) + I2*ABAS(3,2) + I3*ABAS(3,3)
               INSIDE = 0
C
               DO IQQP = 1,NQQP_STR
                  SX = RX - QQPX(IQQP)
                  SY = RY - QQPY(IQQP)
                  SZ = RZ - QQPZ(IQQP)
                  DX = DSQRT(SX*SX+SY*SY+SZ*SZ)
                  IF ( DX.LE.RMAX ) THEN
                     INSIDE = 1
                     SMAX(IQQP) = SMAX(IQQP) + 1
                     IF ( SMAX(IQQP).GT.NRDLMAX ) GOTO 300
                  END IF
               END DO
C
               IF ( INSIDE.NE.0 ) THEN
                  NR = NR + 1
                  IF ( NR.GT.NRDLMAX ) GOTO 100
                  DR(NR) = DSQRT(RX*RX+RY*RY+RZ*RZ)
                  R1(NR) = I1
                  R2(NR) = I2
                  R3(NR) = I3
               END IF
C
            END DO
         END DO
      END DO
C
C -------------------------------- SORT VECTORS IN ORDER OF INCREASING D
C
      NSH = 1
      NSR(1) = 1
C
      DO II = 2,NR
         I = II - 1
         K = I
         P = DR(I)
C
         DO J = II,NR
            IF ( DR(J).LT.P ) THEN
               K = J
               P = DR(J)
            END IF
         END DO
C
         IF ( K.NE.I ) THEN
            DR(K) = DR(I)
            DR(I) = P
C
            IR = R1(I)
            R1(I) = R1(K)
            R1(K) = IR
            IR = R2(I)
            R2(I) = R2(K)
            R2(K) = IR
            IR = R3(I)
            R3(I) = R3(K)
            R3(K) = IR
         END IF
C
         IF ( I.GT.1 ) THEN
            IF ( ABS(DR(I)-DR(I-1)).GT.1.0D-6 ) THEN
               IF ( IPRINT.GT.2 ) WRITE (6,99008) NSH,NSR(NSH)
               NSH = NSH + 1
               NSR(NSH) = 1
            ELSE
               NSR(NSH) = NSR(NSH) + 1
            END IF
         END IF
C
         IF ( IPRINT.GT.2 ) WRITE (6,99009) I,R1(I),R2(I),R3(I),DR(I)
C
         IF ( I.EQ.(NR-1) ) THEN
            J = I + 1
            IF ( ABS(DR(J)-DR(J-1)).GT.1.0D-6 ) THEN
               IF ( IPRINT.GT.2 ) WRITE (6,99008) NSH,NSR(NSH)
               NSH = NSH + 1
               NSR(NSH) = 1
            ELSE
               NSR(NSH) = NSR(NSH) + 1
            END IF
            IF ( IPRINT.GT.2 ) WRITE (6,99009) J,R1(I),R2(I),R3(I),DR(I)
            IF ( IPRINT.GT.2 ) WRITE (6,99008) NSH,NSR(NSH)
         END IF
C
      END DO
C
C
C ======================================================================
C                           RECIPROCAL SPACE
C ======================================================================
C
C  accept a k-point ->G[n] in list if
C
C                     ( (->k+->G[n])**2 - EDU ) < GMAX**2
C
C  i.e. the center of the convergence sphere of radius  GMAX  is at +EDU
C  EDU can vary between  edumin .. edumax  and ->k is within the BZ
C  the set parameter correspond to E=0..3 Ry
C
      EDUMIN = 0.0D0/(2*PI/ALAT)**2
      EDUMAX = 3.0D0/(2*PI/ALAT)**2
      IETOP = 3
      KK = 1
C
      IF ( IPRINT.GT.2 ) WRITE (6,99002)
C
      GMAXSQ = GMAX**2
      NG = 0
C
      DO I1 = -NUMGH, + NUMGH
         DO I2 = -NUMGH, + NUMGH
            DO I3 = -NUMGH, + NUMGH
               GX = I1*BBAS(1,1) + I2*BBAS(1,2) + I3*BBAS(1,3)
               GY = I1*BBAS(2,1) + I2*BBAS(2,2) + I3*BBAS(2,3)
               GZ = I1*BBAS(3,1) + I2*BBAS(3,2) + I3*BBAS(3,3)
               INSIDE = 0
C
               DO J1 = -KK, + KK
                  DO J2 = -KK, + KK
                     DO J3 = -KK, + KK
                        KNX = GX + 0.5D0*(J1*BBAS(1,1)+J2*BBAS(1,2)
     &                        +J3*BBAS(1,3))
                        KNY = GY + 0.5D0*(J1*BBAS(2,1)+J2*BBAS(2,2)
     &                        +J3*BBAS(2,3))
                        KNZ = GZ + 0.5D0*(J1*BBAS(3,1)+J2*BBAS(3,2)
     &                        +J3*BBAS(3,3))
                        KSQ = KNX**2 + KNY**2 + KNZ**2
C
                        DO IE = 0,IETOP
                           IF ( IETOP.NE.0 ) THEN
                              REDU = EDUMIN + IE*(EDUMAX-EDUMIN)
     &                               /DBLE(IETOP)
                           ELSE
                              REDU = 0.0D0
                           END IF
                           IF ( (KSQ-REDU).LE.GMAXSQ ) THEN
                              INSIDE = 1
                              GOTO 10
                           END IF
                        END DO
                     END DO
                  END DO
               END DO
C
 10            CONTINUE
               IF ( INSIDE.NE.0 ) THEN
                  NG = NG + 1
                  IF ( NG.GT.NGRLMAX ) GOTO 200
                  DG(NG) = DSQRT(GX*GX+GY*GY+GZ*GZ)
                  G1(NG) = I1
                  G2(NG) = I2
                  G3(NG) = I3
               END IF
            END DO
         END DO
      END DO
C
C -------------------------------- SORT VECTORS IN ORDER OF INCREASING D
C
      NSH = 1
      NSG(1) = 1
C
      DO II = 2,NG
         I = II - 1
         K = I
         P = DG(I)
C
         DO J = II,NG
            IF ( DG(J).LT.P ) THEN
               K = J
               P = DG(J)
            END IF
         END DO
C
         IF ( K.NE.I ) THEN
            DG(K) = DG(I)
            DG(I) = P
C
            IG = G1(I)
            G1(I) = G1(K)
            G1(K) = IG
            IG = G2(I)
            G2(I) = G2(K)
            G2(K) = IG
            IG = G3(I)
            G3(I) = G3(K)
            G3(K) = IG
         END IF
C
         IF ( I.GT.1 ) THEN
            IF ( ABS(DG(I)-DG(I-1)).GT.1.0D-7 ) THEN
               IF ( IPRINT.GT.2 ) WRITE (6,99008) NSH,NSG(NSH)
               NSH = NSH + 1
               NSG(NSH) = 1
            ELSE
               NSG(NSH) = NSG(NSH) + 1
            END IF
         END IF
C
         IF ( IPRINT.GT.2 ) WRITE (6,99009) I,G1(I),G2(I),G3(I),DG(I)
C
         IF ( I.EQ.(NG-1) ) THEN
            J = I + 1
            IF ( ABS(DG(J)-DG(J-1)).GT.1.0D-7 ) THEN
               IF ( IPRINT.GT.2 ) WRITE (6,99008) NSH,NSG(NSH)
               NSH = NSH + 1
               NSG(NSH) = 1
            ELSE
               NSG(NSH) = NSG(NSH) + 1
            END IF
            IF ( IPRINT.GT.2 ) WRITE (6,99009) J,G1(I),G2(I),G3(I),DG(I)
            IF ( IPRINT.GT.2 ) WRITE (6,99008) NSH,NSG(NSH)
         END IF
C
      END DO
C
      NGRL = NG
      NRDL = NR
C
      WRITE (6,99006) NRDL,NGRL
      IF ( NQQP_STR.GT.1 ) THEN
         IF ( IPRINT.GT.0 ) THEN
            WRITE (6,99007) (I,SMAX(I),I=1,NQQP_STR)
         ELSE
            JHIG = 0
            JLOW = 1000000
            DO I = 1,NQQP_STR
               IF ( SMAX(I).GT.JHIG ) JHIG = SMAX(I)
               IF ( SMAX(I).LT.JLOW ) JLOW = SMAX(I)
            END DO
         END IF
      END IF
C
      DEALLOCATE (DG,DR,NSG,NSR,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC')
      RETURN
C
C-----------------------------------------------------------------------
 100  CONTINUE
      WRITE (6,*) NR,DX,RA,RMAX
      WRITE (6,99003) NRDLMAX,DX,RA,RMAX
      STOP
 200  CONTINUE
      WRITE (6,99004) NGRLMAX,DX,GA,GMAX
      STOP
 300  CONTINUE
      WRITE (6,*) NR,DX,RA,RMAX
      WRITE (6,99005) NRDLMAX,SMAX
      STOP
C-----------------------------------------------------------------------
99001 FORMAT (10X,'result from <STRVECGEN> for real space vectors',//,
     &        14X,'NO',5X,'SX',8X,'SY',8X,'SZ',8X,'D',/)
99002 FORMAT (10X,'result from <STRVECGEN> for reciprocal vectors',//,
     &        14X,'NO',5X,'KX',8X,'KY',8X,'KZ',8X,'D',/)
99003 FORMAT (' *** NR > NRDLMAX = ',I5,/,'   DX,RA,RMAX: ',3F10.5)
99004 FORMAT (' *** NG > NGRLMAX  = ',I5,/,'   DX,GA,GMAX: ',3F10.5)
99005 FORMAT (' *** SMAX > NRDLMAX =',I5,/,(10X,14I5))
99006 FORMAT (10X,'NRDL  =',I10,10X,'NGRL  =',I10)
99007 FORMAT (10X,'IQQP=',I4,': SMAX =',I4)
99008 FORMAT (' ',//,15X,'shell number',I5,' with',I5,' points',//)
99009 FORMAT (' ',10X,I5,3I10,F10.6)
      END
