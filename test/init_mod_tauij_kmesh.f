C*==init_mod_tauij_kmesh.f    processed by SPAG 6.70Rc at 09:36 on 30 Mar 2017
      SUBROUTINE INIT_MOD_TAUIJ_KMESH
C   ********************************************************************
C   *                                                                  *
C   *   parameters for Brillouin zone integration of TAUIJ             *
C   *                                                                  *
C   *   KTABTAUIJ(3,i)     k-point of the full mesh i=1,NKTABTAUIJ     *
C   *   IK_IKTABTAUIJ(i)   index IK of the symmetry reduced mesh       *
C   *   ISYM_IKTABTAUIJ(i) symmetry operation R with k(i) = R k(IK)    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SYMMETRY,ONLY:SYMUNITARY,SYMACCEPTED,MROTK,NSYM
      USE MOD_LATTICE,ONLY:BDBINV,BBAS,BBAS_I,SYSTEM_DIMENSION
      USE MOD_TAUIJ,ONLY:NKTABTAUIJ,IK_IKTABTAUIJ,ISYM_IKTABTAUIJ,
     &    KTABTAUIJ,WKTABTAUIJ
      USE MOD_KSPACE,ONLY:NKTAB,KTAB,WKTAB
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='INIT_MOD_TAUIJ_KMESH')
C
C Local variables
C
      REAL*8 BBAS_LOC(1:3,1:3),BV(3),CF(3),KVECP(3),WK,WKSUM
      REAL*8 DDOT
      INTEGER DELCF,I,IFLAG,IK,IKTABTAUIJ,IKTABTAUIJ_1,ISYM,
     &        NKTABTAUIJ_0,NKTABTAUIJ_IK,NN
      LOGICAL FOUND
      LOGICAL RVEC_SAME
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
C=======================================================================
C                         system dimension
C=======================================================================
C
      IF ( SYSTEM_DIMENSION(1:2).EQ.'3D' ) THEN
         BBAS_LOC(1:3,1:3) = BBAS(1:3,1:3)
      ELSE IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' ) THEN
         BBAS_LOC(1:3,1:3) = BBAS_I(1:3,1:3)
      END IF
C
      IFLAG = 0
C
      NKTABTAUIJ = 0
      DO IK = 1,NKTAB
         NN = NINT(WKTAB(IK))
         IF ( ABS(WKTAB(IK)-NN).GT.1D-8 ) THEN
            IFLAG = 1
            WRITE (6,*) 'INIT_MOD_TAUIJ_KMESH:  WKTAB not INT: ',IK,
     &                  WKTAB(IK)
         END IF
         NKTABTAUIJ = NKTABTAUIJ + NN
      END DO
C
      NKTABTAUIJ = NKTABTAUIJ*2
      NKTABTAUIJ_0 = NKTABTAUIJ
C
      ALLOCATE (IK_IKTABTAUIJ(NKTABTAUIJ),ISYM_IKTABTAUIJ(NKTABTAUIJ))
      ALLOCATE (KTABTAUIJ(3,NKTABTAUIJ),WKTABTAUIJ(NKTABTAUIJ))
C
C
      NKTABTAUIJ = 0
      WKSUM = 0D0
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      DO IK = 1,NKTAB
         IKTABTAUIJ_1 = NKTABTAUIJ + 1
         WKSUM = WKSUM + WKTAB(IK)
C
         NKTABTAUIJ_IK = 0
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
         DO ISYM = 1,NSYM
            IF ( SYMACCEPTED(ISYM) ) THEN
C
C-NOTE:  MROTK includes the INVERSION for ANTI - unitary transformations
C
               CALL DGEMV('N',3,3,1D0,MROTK(1,1,ISYM),3,KTAB(1,IK),1,
     &                    0D0,KVECP,1)
C
C----------- shift KVECP back to parallelepiped spanned by basis vectors
C
               DO I = 1,3
                  BV(I) = DDOT(3,BBAS_LOC(1,I),1,KVECP,1)
               END DO
               CALL DGEMV('N',3,3,1D0,BDBINV,3,BV,1,0D0,CF,1)
C
               DO I = 1,3
                  DELCF = INT(CF(I))
                  IF ( CF(I).LT.0D0 ) DELCF = DELCF - 1
                  CF(I) = CF(I) - DELCF
                  IF ( ABS(1D0-CF(I)).LT.1D-8 ) CF(I) = 0D0
                  IF ( CF(I).LT.0D0 .OR. CF(I).GE.1D0 ) THEN
                     WRITE (6,*) '********* CF ',IK,I,CF(I),DELCF
                     IFLAG = 1
                  END IF
               END DO
C
               CALL RVECLCRB(CF(1),CF(2),CF(3),BBAS_LOC,KVECP)
C
               FOUND = .FALSE.
               DO IKTABTAUIJ = IKTABTAUIJ_1,NKTABTAUIJ
                  FOUND = RVEC_SAME(3,KVECP,KTABTAUIJ(1,IKTABTAUIJ),
     &                    1D-6)
                  IF ( FOUND ) THEN
                     IF ( SYMUNITARY(ISYM) ) ISYM_IKTABTAUIJ(IKTABTAUIJ)
     &                    = ISYM
                     EXIT
                  END IF
               END DO
               IF ( .NOT.FOUND ) THEN
                  NKTABTAUIJ = NKTABTAUIJ + 1
                  IF ( NKTABTAUIJ.LE.NKTABTAUIJ_0 ) THEN
                     KTABTAUIJ(1:3,NKTABTAUIJ) = KVECP(1:3)
                     IK_IKTABTAUIJ(NKTABTAUIJ) = IK
                     ISYM_IKTABTAUIJ(NKTABTAUIJ) = ISYM
                     NKTABTAUIJ_IK = NKTABTAUIJ_IK + 1
                  ELSE
                     IFLAG = 1
                     WRITE (6,*) ' NKTABTAUIJ = ',NKTABTAUIJ,
     &                           ' > NKTABTAUIJ_0 = ',NKTABTAUIJ_0
                  END IF
               END IF
C
            END IF
         END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
         IF ( NKTABTAUIJ_IK.NE.NINT(WKTAB(IK)) ) THEN
            WRITE (6,99001) IK,NKTABTAUIJ_IK,WKTAB(IK)
            IFLAG = 1
         END IF
C
      END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
      WK = 1D0/DBLE(NKTABTAUIJ)
      WKTABTAUIJ(1:NKTABTAUIJ) = WK
C
      WRITE (6,99002) NKTAB,NKTABTAUIJ,NINT(WKSUM)
      IF ( IFLAG.EQ.1 ) THEN
         WRITE (6,99003)
         STOP '<INIT_MOD_TAUIJ_KMESH>'
      END IF
C
99001 FORMAT (5X,'IK =',I8,':     NKTABTAUIJ_IK =',I8,
     &        ' <> WKTAB(IK) = ',F8.2)
99002 FORMAT (/,1X,79('*'),/,29X,'<INIT_MOD_TAUIJ_KMESH>',/,1X,79('*'),
     &        //,10X,'NKTAB =',I8,'     NKTABTAUIJ =',I8,/,30X,
     &        'WKSUM      =',I8,/)
99003 FORMAT (/,1X,79('#'),/,29X,'<INIT_MOD_TAUIJ_KMESH>',/,1X,79('#'),
     &        //,5X,'stop due to inconsistencies -- see output above ',
     &        /)
      END
