C*==vbpesmetr.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE VBPESMETR(MA,MB,BAR,TREATIRR,MAIRR,MBIRR,MPOL,N,M)
C   ********************************************************************
C   *                                                                  *
C   *     transform the matrix elements   <Z | H_lam | Z>  for         *
C   *     polarisation  lam = (+),(-),(z) to new polarisation lam'     *
C   *     according to the rotation matrix  MPOL                       *
C   *                                                                  *
C   *     MA = <Z^+ | H_lam | Z>  -->  MB  = <Z^+ | H_lam' | Z>        *
C   *                                                                  *
C   *     MAIRR =  <Z+ | H_lam | Z J | H^+_lam" | Z>                   *
C   *         -->  MBIRR =  <Z+ | H_lam' | Z J | H^+_lam' | Z>         *
C   *     only the lam'-diagonal case is considered                    *
C   *                                                                  *
C   *     TREATIRR = .f.  skip treatment of irregular matrix           *
C   *                     <Z^+ | H_lam | Z J^+ | H^+_lam' | Z>         *
C   *                                                                  *
C   *     BAR = .t.  assume matrix MA to be MBAR = <Z^x | H_lam | Z>   *
C   *                and deliver  <Z^x | H^+_lam' | Z>                 *
C   *                treat adjoint operator H^+ by reversing pol. lam  *
C   *                skip treatment of irregular matrix                *
C   *                                                                  *
C   *     M: matrix array size  N: true dimension                      *
C   *     NOTE: the matrices are dealt with up to index M              *
C   *           i.e. the FULL matrices MA, MAIRR should be initialized *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL BAR,TREATIRR
      INTEGER M,N
      COMPLEX*16 MA(M,M,3),MAIRR(M,M,M,3,3),MB(M,M,3),MBIRR(M,M,M,3),
     &           MPOL(3,3)
C
C Local variables
C
      INTEGER IPOL,JPOL,JPOLREV,KPOL,KPOLREV,NUM
      COMPLEX*16 W
C
C*** End of declarations rewritten by SPAG
C
      NUM = M*M*3
      CALL CINIT(NUM,MB)
C
      NUM = M*N
C
C-----------------------------------------------------------------------
C          MA = <Z^x | H_lam | Z>  -->  MB  = <Z^x | H^+_lam' | Z>
C-----------------------------------------------------------------------
C
      IF ( BAR ) THEN
C
         DO IPOL = 1,3
            DO JPOL = 1,3
               IF ( JPOL.EQ.3 ) THEN
                  JPOLREV = JPOL
               ELSE
                  JPOLREV = 3 - JPOL
               END IF
               W = DCONJG(MPOL(IPOL,JPOL))
               CALL ZAXPY(NUM,W,MA(1,1,JPOLREV),1,MB(1,1,IPOL),1)
            END DO
         END DO
C
         RETURN
C
      ELSE
C
C-----------------------------------------------------------------------
C          MA = <Z^+ | H_lam | Z>  -->  MB  = <Z^+ | H_lam' | Z>
C-----------------------------------------------------------------------
C
         DO IPOL = 1,3
            DO JPOL = 1,3
               W = MPOL(IPOL,JPOL)
               CALL ZAXPY(NUM,W,MA(1,1,JPOL),1,MB(1,1,IPOL),1)
            END DO
         END DO
C
      END IF
C
C-----------------------------------------------------------------------
C          MAIRR =  <Z^+ | H_lam  | Z J^+ | H^+_lam" | Z>
C     -->  MBIRR =  <Z^+ | H_lam' | Z J^+ | H^+_lam' | Z>
C-----------------------------------------------------------------------
C
      IF ( TREATIRR ) THEN
C
         NUM = M*M*M*3
         CALL CINIT(NUM,MBIRR)
C
         NUM = M*M*N
C
         DO IPOL = 1,3
            DO JPOL = 1,3
               DO KPOL = 1,3
                  IF ( KPOL.EQ.3 ) THEN
                     KPOLREV = KPOL
                  ELSE
                     KPOLREV = 3 - KPOL
                  END IF
C
                  W = MPOL(IPOL,JPOL)*DCONJG(MPOL(IPOL,KPOL))
                  CALL ZAXPY(NUM,W,MAIRR(1,1,1,JPOL,KPOLREV),1,
     &                       MBIRR(1,1,1,IPOL),1)
               END DO
            END DO
         END DO
C
      END IF
C
      END
C*==vbpesinttr_0.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE VBPESINTTR_0(ISS,TOT,POL,NSPIN)
C   ********************************************************************
C   *                                                                  *
C   *             get polarisation and total intensity                 *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_CONSTANTS,ONLY:C0,CI,C1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 ISS(2,2,3)
      REAL*8 NSPIN(3),POL(2,3),TOT(3)
C
C Local variables
C
      REAL*8 DDOT
      INTEGER I,IMS,IMSP,IPOL
      REAL*8 P(3),PN
      COMPLEX*16 PAULI(2,2,1:3),POLMAT(2,2,3),RHO(2,2)
C
C*** End of declarations rewritten by SPAG
C
C
      PAULI(1,1,1) = C0
      PAULI(2,1,1) = C1
      PAULI(1,2,1) = C1
      PAULI(2,2,1) = C0
C
      PAULI(1,1,2) = C0
      PAULI(2,1,2) = CI
      PAULI(1,2,2) = -CI
      PAULI(2,2,2) = C0
C
      PAULI(1,1,3) = C1
      PAULI(2,1,3) = C0
      PAULI(1,2,3) = C0
      PAULI(2,2,3) = -C1
C
C=======================================================================
      DO IPOL = 1,3
C
C--------------------------------------------------- spin density matrix
C
         DO IMS = 1,2
            DO IMSP = 1,2
               RHO(IMS,IMSP) = (ISS(IMS,IMSP,IPOL)-DCONJG(ISS(IMSP,IMS,
     &                         IPOL)))/(2.0D0*CI)
            END DO
         END DO
C
C------------------------------------------------------- total intensity
C
         TOT(IPOL) = DREAL(RHO(1,1)+RHO(2,2))
C
C--------------------------------------------------- polarisation vector
C
         DO I = 1,3
            CALL CMATMUL(2,2,PAULI(1,1,I),RHO,POLMAT(1,1,IPOL))
            IF ( ABS(TOT(IPOL)).GT.1D-12 ) THEN
               P(I) = DREAL(POLMAT(1,1,IPOL)+POLMAT(2,2,IPOL))/TOT(IPOL)
            ELSE
               P(I) = 0.0D0
            END IF
         END DO
C
C   JM: X-component of the polarisation vector is reversed
C       This we got by experimental calculations. it still has
C       to be understand
C
         P(1) = -1D0*P(1)
C
C------------------------------------------ spin-projected photo-current
C
         PN = DDOT(3,NSPIN,1,P(1),1)
C
         POL(1,IPOL) = (1D0-PN)*TOT(IPOL)/2D0
         POL(2,IPOL) = (1D0+PN)*TOT(IPOL)/2D0
C
C-----------------------------------------------------------------------
      END DO
C
      END
C*==polconv.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE POLCONV(TETPHOT,PHIPHOT,QVECPHOT,QHAT,KGEO,MGEO)
C   ********************************************************************
C   *                                                                  *
C   *  - fix the photon wave vector  QVECPHOT and normalize (QHAT)     *
C   *    KGEO = 3:   (TET,PHI) supplied in degree                      *
C   *    KGEO = 4:   QVECPHOT supplied -> normalize and fix  (TET,PHI) *
C   *                                                                  *
C   *  - get polarisation axis E1 and E2                               *
C   *    with  E1 x E2 = QHAT   and   E1 in the xy-plane               *
C   *                                                                  *
C   *  - get the polarisation conversion matrix  MGEO                  *
C   *    that allows to express the matrix elements for (+,-,z)        *
C   *    w.r.t. the local frame of reference to that fixed by          *
C   *    E1, E2, and QHAT                                              *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C1,C0,CI,PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='POLCONV')
      REAL*8 R05
      PARAMETER (R05=1D0/1.4142135623730951D0)
      COMPLEX*16 CR05,IR05,MIR05
      PARAMETER (CR05=C1*R05,IR05=CI*R05,MIR05=-IR05)
C
C Dummy arguments
C
      INTEGER KGEO
      REAL*8 PHIPHOT,TETPHOT
      COMPLEX*16 MGEO(3,3)
      REAL*8 QHAT(3),QVECPHOT(3)
C
C Local variables
C
      COMPLEX*16 BMAT(3,3),CMAT(3,3)
      REAL*8 DDOT,DNRM2
      REAL*8 E1(3),E2(3),EX(3),EY(3),EZ(3),MEX(3),MEY(3),MEZ(3),QXY,
     &       RNORM,TOL,V
      INTEGER I,IE,J
      LOGICAL RVEC_SAME
      COMPLEX*16 ZDOTU
C
C*** End of declarations rewritten by SPAG
C
      DATA EX/1D0,0D0,0D0/,EY/0D0,1D0,0D0/,EZ/0D0,0D0,1D0/
      DATA MEX/ - 1D0,0D0,0D0/,MEY/0D0, - 1D0,0D0/,MEZ/0D0,0D0, - 1D0/
      DATA BMAT/CR05,IR05,C0,CR05,MIR05,C0,C0,C0,C1/
C
      TOL = 1D-6
C
      CALL CMATSTRUCT('BMAT ',BMAT,3,3,0,0,0,1D-8,6)
      IF ( KGEO.EQ.3 ) THEN
C=======================================================================
C      KGEO = 3  input:  TETPHOT,PHIPHOT
C=======================================================================
C
         QVECPHOT(3) = -COS(TETPHOT*PI/180D0)
         QXY = SIN(TETPHOT*PI/180D0)
         QVECPHOT(1) = -COS(PHIPHOT*PI/180D0)*QXY
         QVECPHOT(2) = -SIN(PHIPHOT*PI/180D0)*QXY
C
         CALL DCOPY(3,QVECPHOT,1,QHAT,1)
C
      ELSE
C=======================================================================
C      KGEO = 4  input:  QPHOTTET
C=======================================================================
C
         CALL DCOPY(3,QVECPHOT,1,QHAT,1)
         RNORM = 1D0/DNRM2(3,QHAT,1)
         CALL DSCAL(3,RNORM,QHAT,1)
C
         TETPHOT = ACOS(-QHAT(3))
         QXY = SQRT(QHAT(1)*QHAT(1)+QHAT(2)*QHAT(2))
         IF ( ABS(QXY).LE.TOL ) THEN
            PHIPHOT = 0D0
         ELSE IF ( (-QHAT(2)).GE.0D0 ) THEN
            PHIPHOT = ACOS(-QHAT(1)/QXY)*180D0/PI
         ELSE
            PHIPHOT = (2D0*PI-ACOS(-QHAT(1)/QXY))*180D0/PI
         END IF
      END IF
C=======================================================================
C
      IF ( RVEC_SAME(3,QHAT,EX,TOL) ) THEN
         CALL DCOPY(3,MEY,1,E1,1)
C=======================================================================
      ELSE IF ( RVEC_SAME(3,QHAT,MEX,TOL) ) THEN
         CALL DCOPY(3,EY,1,E1,1)
C=======================================================================
      ELSE IF ( RVEC_SAME(3,QHAT,EY,TOL) ) THEN
         CALL DCOPY(3,EX,1,E1,1)
C=======================================================================
      ELSE IF ( RVEC_SAME(3,QHAT,MEY,TOL) ) THEN
         CALL DCOPY(3,EX,1,E1,1)
C=======================================================================
      ELSE IF ( RVEC_SAME(3,QHAT,EZ,TOL) ) THEN
         CALL DCOPY(3,EX,1,E1,1)
C=======================================================================
      ELSE IF ( RVEC_SAME(3,QHAT,MEZ,TOL) ) THEN
         CALL DCOPY(3,EX,1,E1,1)
C=======================================================================
      ELSE
         E1(1) = QHAT(2)/QXY
         E1(2) = -QHAT(1)/QXY
         E1(3) = 0D0
      END IF
C=======================================================================
C
      CALL RVECXPRO(QHAT,E1,E2)
C
C=======================================================================
C                                            check consistency and print
C
      CALL RVECSPRO(E1,E2,QHAT,V)
      IF ( V.LT.0D0 ) THEN
         WRITE (6,99002)
         CALL DSCAL(3,-1D0,E2,1)
      END IF
C
      WRITE (6,99001) QHAT,E1,E2
C
      IF ( ABS(1D0-DDOT(3,E1,1,E1,1)).GT.TOL )
     &     CALL STOP_MESSAGE(ROUTINE,'case A')
      IF ( ABS(1D0-DDOT(3,E2,1,E2,1)).GT.TOL )
     &     CALL STOP_MESSAGE(ROUTINE,'case B')
      IF ( ABS(DDOT(3,E1,1,E2,1)).GT.TOL )
     &      CALL STOP_MESSAGE(ROUTINE,'case C')
      IF ( ABS(DDOT(3,E1,1,QHAT,1)).GT.TOL )
     &      CALL STOP_MESSAGE(ROUTINE,'case D')
      IF ( ABS(DDOT(3,E2,1,QHAT,1)).GT.TOL )
     &      CALL STOP_MESSAGE(ROUTINE,'case E')
C
C=======================================================================
C                              set up the polarisation conversion matrix
C
      CALL CINIT(3*3,CMAT)
      DO I = 1,3
         DO J = 1,3
            CMAT(I,1) = CMAT(I,1) + DCONJG(BMAT(J,I))*E1(J)
            CMAT(I,2) = CMAT(I,2) + DCONJG(BMAT(J,I))*E2(J)
            CMAT(I,3) = CMAT(I,3) + DCONJG(BMAT(J,I))*QHAT(J)
         END DO
      END DO
C
      CALL CMATSTRUCT('CMAT ',CMAT,3,3,0,0,0,1D-8,6)
C
      MGEO(1,1) = R05*(CMAT(1,1)+CI*CMAT(1,2))
      MGEO(1,2) = R05*(CMAT(2,1)+CI*CMAT(2,2))
      MGEO(1,3) = R05*(CMAT(3,1)+CI*CMAT(3,2))
C
      MGEO(2,1) = R05*(CMAT(1,1)-CI*CMAT(1,2))
      MGEO(2,2) = R05*(CMAT(2,1)-CI*CMAT(2,2))
      MGEO(2,3) = R05*(CMAT(3,1)-CI*CMAT(3,2))
C
      MGEO(3,1) = CMAT(1,3)
      MGEO(3,2) = CMAT(2,3)
      MGEO(3,3) = CMAT(3,3)
C
      WRITE (6,*) ZDOTU(3,MGEO(1,1),1,MGEO(1,1),1)
      WRITE (6,*) ZDOTU(3,MGEO(1,2),1,MGEO(1,2),1)
      WRITE (6,*) ZDOTU(3,MGEO(1,3),1,MGEO(1,3),1)
C
      IE = 0
C      IF ( ABS(1D0-ZDOTU(3,MGEO(1,1),1,MGEO(1,1),1)).GT.TOL ) IE = IE+1
C      IF ( ABS(1D0-ZDOTU(3,MGEO(1,2),1,MGEO(1,2),1)).GT.TOL ) IE = IE+1
C      IF ( ABS(1D0-ZDOTU(3,MGEO(1,3),1,MGEO(1,3),1)).GT.TOL ) IE = IE+1
C      IF ( ABS(ZDOTU(3,MGEO(1,1),1,MGEO(1,2),1)).GT.TOL ) IE = IE+1
C      IF ( ABS(ZDOTU(3,MGEO(1,1),1,MGEO(1,3),1)).GT.TOL ) IE = IE+1
C      IF ( ABS(ZDOTU(3,MGEO(1,2),1,MGEO(1,3),1)).GT.TOL ) IE = IE+1
C
      IF ( IE.NE.0 ) THEN
         WRITE (6,*) 'IE=',IE
         CALL CMATSTRUCT('MGEO ',MGEO,3,3,0,0,0,1D-8,6)
         CALL STOP_MESSAGE(ROUTINE,'MGEO - matrix not unitary ')
      END IF
C
C=======================================================================
C
99001 FORMAT (/,10X,' q(phot)    ',3F10.5,/,10X,' E1         ',3F10.5,/,
     &        10X,' E2         ',3F10.5)
99002 FORMAT (/,10X,'WARNING from <POLCONV>  (E1xE2)*Q < 0',/)
      END
