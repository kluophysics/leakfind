C*==vmuftin.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE VMUFTIN(VMTZ)
C   ********************************************************************
C   *                                                                  *
C   *       muffintinize potential                                     *
C   *                                                                  *
C   *       empty spheres are excluded form the AVERAGING              *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:VT,CONC,NAT,IMT,Z
      USE MOD_RMESH,ONLY:R,R2DRDI,RWS,RMT,JRWS,NRMAX
      USE MOD_TYPES,ONLY:ITBOT,ITTOP
      USE MOD_LATTICE,ONLY:SYSTEM_DIMENSION,SUB_SYSTEM,SYSTEM_TYPE
      IMPLICIT NONE
C*--VMUFTIN16
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='VMUFTIN')
C
C Dummy arguments
C
      REAL*8 VMTZ
C
C Local variables
C
      LOGICAL CALCVMTZ
      REAL*8 F1(:),F2(:),RELDIFF,RSQINT,S1(:),S1MT,S2(:),S2MT,SUM_RSQ,
     &       SUM_V_RSQ
      INTEGER I,IA_ERR,IM,IRTOP,IT
      REAL*8 YLAG
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE F1,F2,S1,S2
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (F1(NRMAX),F2(NRMAX),S1(NRMAX),S2(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: S')
C
      IF ( SYSTEM_DIMENSION(1:2).EQ.'3D' .OR. SYSTEM_TYPE(1:3)
     &     .EQ.'VIV' ) THEN
         CALCVMTZ = .TRUE.
      ELSE IF ( SUB_SYSTEM(1:6).EQ.'I-ZONE' ) THEN
         CALCVMTZ = .FALSE.
      ELSE IF ( SUB_SYSTEM(1:6).EQ.'L-BULK' ) THEN
         CALCVMTZ = .TRUE.
      ELSE
         CALL STOP_MESSAGE(ROUTINE,'R-BULK not yet implemented')
      END IF
C
      IF ( SYSTEM_TYPE(1:16).EQ.'EMBEDDED-CLUSTER' ) CALCVMTZ = .FALSE.
C
      IF ( CALCVMTZ ) THEN
C
         SUM_V_RSQ = 0.0D0
         SUM_RSQ = 0.0D0
C
         DO IT = ITBOT,ITTOP
C
C---------------------------------- EXCLUDE empty spheres from averaging
            IF ( Z(IT).EQ.0 ) CYCLE
C
            IM = IMT(IT)
C
            DO I = 1,JRWS(IM)
               F1(I) = R2DRDI(I,IM)
               F2(I) = VT(I,IT)*F1(I)
            END DO
C
            CALL RRADINT_R(IM,F1,S1)
            CALL RRADINT_R(IM,F2,S2)
C
            S1MT = YLAG(RMT(IM),R(1,IM),S1,0,3,JRWS(IM))
            S2MT = YLAG(RMT(IM),R(1,IM),S2,0,3,JRWS(IM))
C
            RSQINT = (RWS(IM)**3-RMT(IM)**3)/3.0D0
C
            SUM_V_RSQ = SUM_V_RSQ + (S2(JRWS(IM))-S2MT)*CONC(IT)*NAT(IT)
            SUM_RSQ = SUM_RSQ + RSQINT*CONC(IT)*NAT(IT)
C
            RELDIFF = 1D0 - (S1(JRWS(IM))-S1MT)/RSQINT
            IF ( ABS(RELDIFF).GT.1D-5 ) WRITE (6,99003) IT,IM,
     &           (S1(JRWS(IM))-S1MT),RSQINT,RELDIFF,(S2(JRWS(IM))-S2MT)
         END DO
C
         VMTZ = SUM_V_RSQ/SUM_RSQ
C
         WRITE (6,99001) VMTZ
C
      ELSE
C
         WRITE (6,99002) VMTZ
C
      END IF
C
C-----------------------------------------------------------------------
C     shift potential to muffin tin zero for ALL sites
C-----------------------------------------------------------------------
      DO IT = ITBOT,ITTOP
         IM = IMT(IT)
         IRTOP = JRWS(IM)
         VT(1:IRTOP,IT) = VT(1:IRTOP,IT) - VMTZ
         VT((IRTOP+1):NRMAX,IT) = 0.0D0
      END DO
C
      DEALLOCATE (F1,F2,S1,S2,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC')
C
99001 FORMAT (/,1X,79('*'),/,10X,'shift of muffin-tin zero  VMTZ ',
     &        F10.6,/,1X,79('*'),/)
99002 FORMAT (/,1X,79('*'),/,10X,'shift of muffin-tin zero  VMTZ ',
     &        F10.6,'     fixed by host'/,1X,79('*'),/)
99003 FORMAT (10X,'INFO from <VMUFTIN>:  IT =',I3,' IM =',I3,/,20X,
     &        'V_is (num)           ',F20.10,/,20X,
     &        'V_is (ana)           ',F20.10,/,20X,
     &        'rel. diff.           ',F20.10,/,20X,
     &        'Int r^2 V [ws] - [mt]',F20.10)
      END
