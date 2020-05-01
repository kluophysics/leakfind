C*==torque.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TORQUE
C   ********************************************************************
C   *                                                                  *
C   *  calculate the magnetic torque                                   *
C   *                                                                  *
C   *  see for example:                                                *
C   *                                                                  *
C   *  Eq. (A6) in  Staunton et al. PRB,74,144411,(2006)               *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_CONSTANTS,ONLY:CI,C0,PI,RY_EV,CPRE
      USE MOD_SITES,ONLY:NQMAX,NQ,IQBOT,IQTOP
      USE MOD_ENERGY,ONLY:NETAB,WETAB
      USE MOD_ANGMOM,ONLY:NKMQ,NKMMAX,NKM,WKM1,WKM2,IOMT,ISMT,WKM3,WKM4,
     &    TAUQ,MSSQ,AME_G
      USE MOD_FILES,ONLY:IPRINT,FOUND_SECTION,FOUND_REAL_ARRAY,N_FOUND
      IMPLICIT NONE
C*--TORQUE21
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='TORQUE')
C
C Local variables
C
      LOGICAL ANGLES_SUPPLIED
      CHARACTER*1 CHPOL(3),CHPOL1(3)
      COMPLEX*16 CMATTRC
      COMPLEX*16 ERYD,J(:,:,:),UJ(:,:,:)
      INTEGER IE,IPOL,IQ,M,N,NE
      REAL*8 QMPHIMAE(:),QMTETMAE(:),RSUM,T(:),UX,UY,UZ
      CHARACTER*20 STR20
      CHARACTER*4 TXTMVEC(3)
C
C*** End of declarations rewritten by SPAG
C
      DATA CHPOL/'x','y','z'/
      DATA CHPOL1/'-','z','+'/
      DATA TXTMVEC/'spin','orb ','J   '/
      DATA ERYD/(999999D0,999999D0)/
C
      ALLOCATABLE J,T,UJ,QMPHIMAE,QMTETMAE
C
      WRITE (6,99001) 'Calculation of magnetic torque'
C
      IF ( IQTOP.GT.NQ ) CALL STOP_MESSAGE(ROUTINE,'IQTOP > NQ')
C
      ALLOCATE (J(NKMMAX,NKMMAX,3))
      ALLOCATE (UJ(NKMMAX,NKMMAX,NQ))
      ALLOCATE (T(NQ))
      ALLOCATE (QMPHIMAE(NQMAX),QMTETMAE(NQMAX))
C
      TAUQ(:,:,:) = C0
      MSSQ(:,:,:) = C0
      UJ(:,:,:) = C0
      T(:) = 0.0D0
C
      QMPHIMAE(:) = 0.0D0
      QMTETMAE(:) = 0.0D0
C
      NE = NETAB(1)
C
C-----------------------------------  Read in site diagonal TAU and MSSQ
C
      CALL READTAU(9,ERYD,0,NE,TAUQ,0,NQ,.FALSE.,WKM1,WKM2,0,NQ,0,
     &             NKMMAX,NKMMAX,IPRINT)
C
      IF ( NE.EQ.0 ) CALL STOP_MESSAGE(ROUTINE,'NE = 0 in TAU-file')
      IF ( NE.NE.NETAB(1) )
     &      CALL STOP_MESSAGE(ROUTINE,'NE <> NETAB in TAU-file')
C
C-----------------------------------------------------------------------
C                     read rotation angles
C-----------------------------------------------------------------------
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
      ANGLES_SUPPLIED = .FALSE.
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_SET_REAL_ARRAY('THETAQ',QMTETMAE,N_FOUND,NQ,0,
     &                               9999D0,0)
         ANGLES_SUPPLIED = ANGLES_SUPPLIED .OR. FOUND_REAL_ARRAY
         CALL SECTION_SET_REAL_ARRAY('PHIQ',QMPHIMAE,N_FOUND,NQ,0,
     &                               9999D0,0)
         ANGLES_SUPPLIED = ANGLES_SUPPLIED .OR. FOUND_REAL_ARRAY
      END IF
      IF ( .NOT.ANGLES_SUPPLIED ) THEN
         WRITE (6,99006)
         QMTETMAE(IQBOT:IQTOP) = 45D0
         QMPHIMAE(IQBOT:IQTOP) = 0D0
      END IF
C
      QMTETMAE(IQBOT:IQTOP) = QMTETMAE(IQBOT:IQTOP)*PI/180D0
      QMPHIMAE(IQBOT:IQTOP) = QMPHIMAE(IQBOT:IQTOP)*PI/180D0
C
C-- calculate J(IPOL)  angular matrix elements of TOTAL angular momentum
C--------------------------------- spherical coordinates   (-), (0), (+)
C
      J(1:NKM,1:NKM,:) = AME_G(1:NKM,1:NKM,:,IOMT)
     &                   + 0.5D0*AME_G(1:NKM,1:NKM,:,ISMT)
C
      IF ( IPRINT.GT.0 ) THEN
         DO IPOL = 1,3
            STR20 = 'A  '//TXTMVEC(3)//'  ('//CHPOL1(IPOL)//')'
            CALL CMATSTRUCT(STR20,J(1,1,IPOL),NKM,NKMMAX,3,3,0,1D-8,6)
         END DO
      END IF
C
C-- convert polarisation: SPHERICAL (-),(0),(+) to CARTESIAN (x),(y),(z)
C
      CALL CMAT_CONVERT_POLAR(J,'S>C')
C
C-------------------------------------------------------- calculate J**2
C
      IF ( IPRINT.GT.0 ) THEN
C
         WKM1(1:NKM,1:NKM) = 0D0
         DO IPOL = 1,3
            CALL CMATMUL(NKM,NKMMAX,J(1,1,IPOL),J(1,1,IPOL),WKM2)
            WKM1(1:NKM,1:NKM) = WKM1(1:NKM,1:NKM) + WKM2(1:NKM,1:NKM)
         END DO
C
         DO IPOL = 1,3
            STR20 = 'A  '//TXTMVEC(3)//'  ('//CHPOL(IPOL)//')'
            CALL CMATSTRUCT(STR20,J(1,1,IPOL),NKM,NKMMAX,3,3,0,1D-8,6)
         END DO
         CALL CMATSTRUCT('J**2',WKM1(1,1),NKM,NKMMAX,3,3,0,1D-8,6)
      END IF
C
C--------------------------------------------------------- calculate u*J
C
      DO IQ = IQBOT,IQTOP
C
         UX = DSIN(QMTETMAE(IQ))*DCOS(QMPHIMAE(IQ))
         UY = DSIN(QMTETMAE(IQ))*DSIN(QMPHIMAE(IQ))
         UZ = DCOS(QMTETMAE(IQ))
C
         UJ(1:NKM,1:NKM,IQ) = UX*J(1:NKM,1:NKM,1) + UY*J(1:NKM,1:NKM,2)
     &                        + UZ*J(1:NKM,1:NKM,3)
C
         IF ( IPRINT.GT.0 ) THEN
            WRITE (6,*) 'U_x',UX
            WRITE (6,*) 'U_y',UY
            WRITE (6,*) 'U_z',UZ
            CALL CMATSTRUCT('U*J',UJ(1,1,IQ),NKM,NKMMAX,3,3,0,1D-8,6)
         END IF
C
      END DO
C
C-----------------------------------------------------  Calculate torque
C
      LOOP_IE:DO IE = 1,NE
C
         LOOP_IQ:DO IQ = IQBOT,IQTOP
C
            M = NKMMAX
            N = NKMQ(IQ)
C
            CALL READTAU(9,ERYD,IE,NE,WKM1,0,0,.FALSE.,TAUQ(1,1,IQ),
     &                   MSSQ(1,1,IQ),IQ,NQ,1,N,M,IPRINT)
C
            CALL CMATMUL(N,M,UJ(1,1,IQ),MSSQ(1,1,IQ),WKM2)
            CALL CMATMUL(N,M,MSSQ(1,1,IQ),UJ(1,1,IQ),WKM3)
C
            WKM4(1:NKM,1:NKM) = WKM2(1:NKM,1:NKM) - WKM3(1:NKM,1:NKM)
C
            CALL CMATMUL(N,M,TAUQ(1,1,IQ),WKM4(1,1),WKM2)
C
            T(IQ) = T(IQ) + DIMAG(CPRE*CI*WETAB(IE,1)*CMATTRC(N,M,WKM2))
C
         END DO LOOP_IQ
C
      END DO LOOP_IE
C
      WRITE (6,99002) 'Site resolved magnetic torque:'
      RSUM = 0.0D0
      DO IQ = IQBOT,IQTOP
         WRITE (6,99003) IQ,T(IQ)
         RSUM = RSUM + T(IQ)
      END DO
      WRITE (6,99002) 'Total:'
C
      WRITE (6,99004) RSUM,'Ryd'
      WRITE (6,99004) RSUM*1000.0D0*RY_EV,'meV'
      WRITE (6,99005)
C
99001 FORMAT (//,1X,79('#'),/,10X,A,/,1X,79('#'),/)
99002 FORMAT (//,10X,A,/)
99003 FORMAT (10X,'IQ =',I3,2X,F20.10)
99004 FORMAT (10X,F20.10,5X,A)
99005 FORMAT (//,1X,79('#'))
99006 FORMAT (/,10X,'no tilting angles supplied by input file',/,10X,
     &        'the default values will be used: TET=45 PHI=0',/)
      END
