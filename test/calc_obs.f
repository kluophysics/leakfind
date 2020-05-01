C*==calc_obs.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE CALC_OBS(IT,TMAT,MZBZA,MZBJA,WE,OBS_T,NOBS)
C   ********************************************************************
C   *                                                                  *
C   *   evaluate the energy dependent contribution when calculating    *
C   *   NOBS observables according to                                  *
C   *                                                                  *
C   *      <A> = -1/PI Im Sum_i WE_i Trace MZBZA_i TAU_i - MZBJA_i     *
C   *                                                                  *
C   *   with the matrix elements                                       *
C   *                                                                  *
C   *             MZBZA  =   < Z^+(E_b) | H_lam | Z(E_a) >             *
C   *             MZBJA  =   < Z^+(E_b) | H_lam | J(E_a) >             *
C   *                      + < J^+(E_b) | H_lam | Z(E_a) >             *
C   *                                                                  *
C   *     with    H_lam =  1, -> sigma, -> l, etc.                     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:IQAT
      USE MOD_CONSTANTS,ONLY:PI
      USE MOD_TYPES,ONLY:NTMAX
      USE MOD_ANGMOM,ONLY:NKMMAX,WKM1,WKM2,NKMQ,NOBSMAX,IDOS,NPOL
      IMPLICIT NONE
C*--CALC_OBS25
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IT,NOBS
      COMPLEX*16 WE
      COMPLEX*16 MZBJA(NKMMAX,NKMMAX,3,NOBS),MZBZA(NKMMAX,NKMMAX,3,NOBS)
     &           ,TMAT(NKMMAX,NKMMAX)
      REAL*8 OBS_T(0:3,NOBSMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 CMATTRC
      COMPLEX*16 DOBS(3)
      INTEGER IOBS,IPOL,IQ,M,N
C
C*** End of declarations rewritten by SPAG
C
C-----------------------------------------------------------------------
C       calculate differential expectation value   DOBS  and integrate
C-----------------------------------------------------------------------
C
      IQ = IQAT(1,IT)
      M = NKMMAX
      N = NKMQ(IQ)
C
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      DO IOBS = 1,NOBS
C
         DO IPOL = 1,NPOL
C
            WKM1(1:N,1:N) = MATMUL(MZBZA(1:N,1:N,IPOL,IOBS),TMAT(1:N,1:N
     &                      ))
C
            WKM2(1:N,1:N) = WKM1(1:N,1:N) - MZBJA(1:N,1:N,IPOL,IOBS)
C
            DOBS(IPOL) = -WE*CMATTRC(N,M,WKM2)/PI
C
         END DO
C
         OBS_T(1:3,IOBS,IT) = OBS_T(1:3,IOBS,IT) + DIMAG(DOBS(1:3))
C
      END DO
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
C
C-------------------------------------- store charge only in z-component
      OBS_T(1:2,IDOS,IT) = 0D0
C
C
      END
C*==cmat_convert_polar.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE CMAT_CONVERT_POLAR(A,KEY)
C   ********************************************************************
C   *                                                                  *
C   *  convert a complex NxN matrix A depending on the polarisation    *
C   *  from SPHERICAL to CARTESIAN coordinates or the other direction  *
C   *                                                                  *
C   *  KEY='S>C'   change polarisation from SPHERICAL to CARTESIAN     *
C   *  KEY='C>S'   change polarisation from to CARTESIAN SPHERICAL     *
C   *                                                                  *
C   *  indexing:  spherical   (-), (0), (+)                            *
C   *             cartesian   (x), (y), (z)                            *
C   *                                                                  *
C   *  the 3x3 transformation matrices  U_SC  and  U_CS                *
C   *  are set up in <INIT_MOD_ANGMOM>                                 *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,U_SC,U_CS
      IMPLICIT NONE
C*--CMAT_CONVERT_POLAR108
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CMAT_CONVERT_POLAR')
C
C Dummy arguments
C
      CHARACTER*3 KEY
      COMPLEX*16 A(NKMMAX,NKMMAX,3)
C
C Local variables
C
      INTEGER I,J
      COMPLEX*16 WC(3),WS(3)
C
C*** End of declarations rewritten by SPAG
C
C-- convert polarisation: SPHERICAL (-),(0),(+) to CARTESIAN (x),(y),(z)
C
      IF ( KEY.EQ.'S>C' ) THEN
C
         DO J = 1,NKM
            DO I = 1,NKM
C
               WC(1:3) = MATMUL(U_CS,A(I,J,1:3))
               A(I,J,1:3) = WC(1:3)
C
            END DO
         END DO
C
C-- convert polarisation: CARTESIAN (x),(y),(z) to SPHERICAL (-),(0),(+)
C
      ELSE IF ( KEY.EQ.'C>S' ) THEN
C
         DO J = 1,NKM
            DO I = 1,NKM
C
               WS(1:3) = MATMUL(U_SC,A(I,J,1:3))
               A(I,J,1:3) = WS(1:3)
C
            END DO
         END DO
C
      ELSE
C
         CALL STOP_MESSAGE(ROUTINE,'WRONG KEY')
C
      END IF
C
      END
C*==cmat_convert_polar2.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE CMAT_CONVERT_POLAR2(A,KEY)
C   ********************************************************************
C   *                                                                  *
C   *  convert a complex NxN matrix A                                  *
C   *  depending twofold on the polarisation                           *
C   *  from SPHERICAL to CARTESIAN coordinates or the other direction  *
C   *                                                                  *
C   *  KEY='S>C'   change polarisation from SPHERICAL to CARTESIAN     *
C   *  KEY='C>S'   change polarisation from to CARTESIAN SPHERICAL     *
C   *                                                                  *
C   *  indexing:  spherical   (-), (0), (+)                            *
C   *             cartesian   (x), (y), (z)                            *
C   *                                                                  *
C   *  the 3x3 transformation matrices  U_SC  and  U_CS                *
C   *  are set up in <INIT_MOD_ANGMOM>                                 *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,U_SC,U_CS
      IMPLICIT NONE
C*--CMAT_CONVERT_POLAR2197
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CMAT_CONVERT_POLAR2')
C
C Dummy arguments
C
      CHARACTER*3 KEY
      COMPLEX*16 A(NKMMAX,NKMMAX,3,3)
C
C Local variables
C
      COMPLEX*16 B(3,3),C(3,3),U_CS_T(3,3),U_SC_T(3,3)
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
C-- convert polarisation: SPHERICAL (-),(0),(+) to CARTESIAN (x),(y),(z)
C
      IF ( KEY.EQ.'S>C' ) THEN
C
         U_CS_T = TRANSPOSE(U_CS)
C
         DO J = 1,NKM
            DO I = 1,NKM
C
               B(1:3,1:3) = A(I,J,1:3,1:3)
               C(1:3,1:3) = MATMUL(U_CS,B)
               A(I,J,1:3,1:3) = MATMUL(C,U_CS_T)
C
            END DO
         END DO
C
C-- convert polarisation: CARTESIAN (x),(y),(z) to SPHERICAL (-),(0),(+)
C
      ELSE IF ( KEY.EQ.'C>S' ) THEN
C
         U_SC_T = TRANSPOSE(U_SC)
C
         DO J = 1,NKM
            DO I = 1,NKM
C
               B(1:3,1:3) = A(I,J,1:3,1:3)
               C(1:3,1:3) = MATMUL(U_SC,B)
               A(I,J,1:3,1:3) = MATMUL(C,U_SC_T)
C
            END DO
         END DO
C
      ELSE
C
         CALL STOP_MESSAGE(ROUTINE,'WRONG KEY')
C
      END IF
C
      END
C*==cvec_convert_coordinates.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE CVEC_CONVERT_COORDINATES(V,KEY)
C   ********************************************************************
C   *                                                                  *
C   *  convert a complex 3 dimensional vector V                        *
C   *  from SPHERICAL to CARTESIAN coordinates or the other direction  *
C   *                                                                  *
C   *  KEY='S>C'   change polarisation from SPHERICAL to CARTESIAN     *
C   *  KEY='C>S'   change polarisation from to CARTESIAN SPHERICAL     *
C   *                                                                  *
C   *  indexing:  spherical   (-), (0), (+)                            *
C   *             cartesian   (x), (y), (z)                            *
C   *                                                                  *
C   *  the 3x3 transformation matrices  U_SC  and  U_CS                *
C   *  are set up in <INIT_MOD_ANGMOM>                                 *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:U_SC,U_CS
      IMPLICIT NONE
C*--CVEC_CONVERT_COORDINATES291
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CVEC_CONVERT_COORDINATES')
C
C Dummy arguments
C
      CHARACTER*3 KEY
      COMPLEX*16 V(3)
C
C Local variables
C
      COMPLEX*16 WC(3),WS(3)
C
C*** End of declarations rewritten by SPAG
C
C--- convert coordinates: SPHERICAL (-),(0),(+) to CARTESIAN (x),(y),(z)
C
      IF ( KEY.EQ.'S>C' ) THEN
C
         WC = MATMUL(U_CS,V(1:3))
C
         V(1:3) = WC(1:3)
C
C--- convert coordinates: CARTESIAN (x),(y),(z) to SPHERICAL (-),(0),(+)
C
      ELSE IF ( KEY.EQ.'C>S' ) THEN
C
         WS(1:3) = MATMUL(U_SC,V(1:3))
C
         V(1:3) = WS(1:3)
C
      ELSE
C
         CALL STOP_MESSAGE(ROUTINE,'WRONG KEY')
C
      END IF
C
      END
