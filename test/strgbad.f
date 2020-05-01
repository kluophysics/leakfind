C*==strgbad.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE STRGBAD(KTET,NKTET,MROTK,IWEDGEROT,NWEDGE)
C   ********************************************************************
C   *                                                                  *
C   *   GENERATE SET OF RECIP. VECTORS THAT MIGHT CAUSE A FREE         *
C   *   ELECTRON POLE   (K+GN)**2 - E = 0                              *
C   *                                                                  *
C   ********************************************************************
      USE MOD_ENERGY,ONLY:ETAB,NETAB,NEMAX
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_KSPACE,ONLY:NGBAD,NGBADMAX,GBAD
      USE MOD_LATTICE,ONLY:ALAT,BBAS
      USE MOD_STR,ONLY:NGRL,G1,G2,G3
      USE MOD_CONSTANTS,ONLY:PI
      USE MOD_SYMMETRY,ONLY:NSYMMAX
      IMPLICIT NONE
C*--STRGBAD17
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='STRGBAD')
C
C Dummy arguments
C
      INTEGER NKTET,NWEDGE
      INTEGER IWEDGEROT(*)
      REAL*8 KTET(3,NKTET),MROTK(3,3,NSYMMAX)
C
C Local variables
C
      REAL*8 GBADTMP(:,:),GX,GY,GZ,KNSQ,KVEC(3),REDU,RRE,RRK
      INTEGER I,IA_ERR,IE,IIE,IIK,IK,IROT,IWEDGE,IX,NG,NNE,NNK
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE GBADTMP
C
      NGBADMAX = NGRL
C
      ALLOCATE (GBADTMP(3,NGBADMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: GBAD')
C
      IF ( NEMAX.NE.UBOUND(ETAB,1) )
     &      CALL STOP_MESSAGE(ROUTINE,'NEMAX .NE. DIM(ETAB)')
cc      IF ( NWEDGE.GT.UBOUND(IWEDGEROT,1) )
cc     &      CALL STOP_MESSAGE(ROUTINE,'NWEDGE > DIM(IWEDGEROT)')
C
C----------------------check only NNK k-points and NNE E-values in table
C
      NNE = MIN(5,NETAB(1))
      NNK = MIN(50,NKTET)
      RRE = DBLE(NETAB(1))/DBLE(NNE)
      RRK = DBLE(NKTET)/DBLE(NNK)
C
      NG = 0
C
      DO I = 1,NGRL
         GX = G1(I)*BBAS(1,1) + G2(I)*BBAS(1,2) + G3(I)*BBAS(1,3)
         GY = G1(I)*BBAS(2,1) + G2(I)*BBAS(2,2) + G3(I)*BBAS(2,3)
         GZ = G1(I)*BBAS(3,1) + G2(I)*BBAS(3,2) + G3(I)*BBAS(3,3)
C
         DO IWEDGE = 1,NWEDGE
            IROT = IWEDGEROT(IWEDGE)
C
            DO IIK = 1,NNK
               IK = MIN(NINT(IIK*RRK),NKTET)
C
               CALL DGEMV('N',3,3,1D0,MROTK(1,1,IROT),3,KTET(1,IK),1,
     &                    0D0,KVEC,1)
C
               KNSQ = (GX+KVEC(1))**2 + (GY+KVEC(2))**2 + (GZ+KVEC(3))
     &                **2
C
               DO IIE = 1,NNE
                  IE = MIN(NINT(IIE*RRE),NETAB(1))
C
                  REDU = DBLE(ETAB(IE,1)/(2*PI/ALAT)**2)
C
                  IF ( ABS(KNSQ-REDU).LE.1D0 ) THEN
                     NG = NG + 1
                     IF ( NG.LE.NGBADMAX ) THEN
                        GBADTMP(1,NG) = GX
                        GBADTMP(2,NG) = GY
                        GBADTMP(3,NG) = GZ
                     END IF
                     GOTO 100
                  END IF
               END DO
            END DO
         END DO
 100  END DO
C
      NGBAD = NG
C
      IF ( NG.GT.NGBADMAX ) THEN
         WRITE (6,*) ' STOP in <STRGBAD>  NGBAD =',NGBAD
         WRITE (6,*) ' array size:        NGBADMAX=',NGBADMAX
         STOP
      END IF
C
C=======================================================================
C                          store results
C=======================================================================
      NGBADMAX = NGBAD
C
      ALLOCATE (GBAD(3,NGBADMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'ALLOC: INIT_MOD_KSPACE -> GBAD'
C
      GBAD(1:3,1:NGBADMAX) = GBADTMP(1:3,1:NGBADMAX)
C
      DEALLOCATE (GBADTMP)
C
C-----------------------------------------------------------------------
      WRITE (6,99001) 'BAD G-vectors    NGBAD: ',NGBAD,NGBADMAX
      WRITE (6,99001)
C
      IF ( IPRINT.LT.2 ) RETURN
C
      DO I = 1,NGBAD
         WRITE (6,'(I5,3F10.4)') I,(GBAD(IX,I),IX=1,3)
      END DO
C
99001 FORMAT (10X,A,10I5)
      END
