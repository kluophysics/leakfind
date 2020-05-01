C*==fpefg.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPEFG(CMNTMTT,R2RHO,V)
C   ********************************************************************
C   *                                                                  *
C   *  calculates the electric field gradients from a given non -      *
C   *  spherical  charge density . the lattice summation is done       *
C   *  implicitly by using the coulomb potential at the sphere         *
C   *  boundary - suggested by m. weinert (private com. 1984)          *
C   *  the different components of the electric field gradients        *
C   *  stored in the following way :                                   *
C   *                                                                  *
C   *           EFG(..,1) : intra atomic contribution                  *
C   *           EFG(..,2) : extra atomic contribution                  *
C   *           EFG(..,3) : total contribution                         *
C   *                                                                  *
C   *                            b.drittler   sep. 1988                *
C   *------------------------------------------------------------------*
C   *  transformation of EFGLM (spherical coordinates)                 *
C   *                 to EFG   (cartesian coordinates)  see:           *
C   *  P. Herzig, Theoret. Chim. Acta (1985), p. 323-333, eqs. 40, 41  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IOTMP,IFILBUILDBOT,WRBUILDBOT
      USE MOD_RMESH,ONLY:NRMAX,R,JRMT,DRDI_W_RADINT
      USE MOD_TYPES,ONLY:TXT_T,NLMFPMAX,NLFP,NTMAX,IMT,LTXT_T,ITBOT,
     &    ITTOP
      USE MOD_CONSTANTS,ONLY:PI,E0_CGS,A0_CGS,C_SI
      IMPLICIT NONE
C*--FPEFG30
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='FPEFG')
      INTEGER LWORK
      PARAMETER (LWORK=12)
      REAL*8 EFG_AU2CGS,EFG_CGS2SI
      PARAMETER (EFG_AU2CGS=E0_CGS/A0_CGS**3,EFG_CGS2SI=C_SI*1D-2)
C
C Dummy arguments
C
      REAL*8 CMNTMTT(NLMFPMAX,NTMAX),R2RHO(NRMAX,NLMFPMAX,NTMAX,3),
     &       V(NRMAX,NLMFPMAX,NTMAX)
C
C Local variables
C
      REAL*8 AUX,EFG(3,3,3),EFGLM(-2:2,3),ETA,EVAL(3,3),EVEC(3,3),FAC,
     &       RMTIN,RNORM,RPWM3(:),WORK(LWORK)
      INTEGER I,IM,INFO,IR,IRMTIN,IT,LM,M,N
      CHARACTER*12 TEXT(3)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE RPWM3
C
      DATA TEXT/'intra atomic','extra atomic','total      '/
C
      ALLOCATE (RPWM3(NRMAX))
C
      IF ( NLFP.LT.3 ) THEN
         WRITE (6,99001)
         RETURN
      END IF
C
      FAC = SQRT(5.0D0/(4.0D0*PI))/2.0D0
C
      WRITE (6,99002)
C
      DO IT = ITBOT,ITTOP
C
         IF ( IT.GT.ITBOT ) WRITE (6,99003)
         WRITE (6,99004) IT,TXT_T(IT)(1:LTXT_T(IT))
C
         IM = IMT(IT)
         IRMTIN = JRMT(IM)
         RMTIN = R(IRMTIN,IM)
C
         DO IR = 2,IRMTIN
            RPWM3(IR) = 1D0/R(IR,IM)**3
         END DO
C
C-----------------------------------------------------------------------
C                        plot radial data
C-----------------------------------------------------------------------
         CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,'efg.dat')
C
         DO IR = 2,IRMTIN
            WRITE (IOTMP,'(10E16.5)') R(IR,IM),R2RHO(IR,7,IT,1),
     &                                R2RHO(IR,7,IT,1)*RPWM3(IR)
         END DO
         CLOSE (IOTMP)
C-----------------------------------------------------------------------
C
         DO M = -2,2
            LM = 6 + M + 1
C
C---> integrate with simpson subroutine
            AUX = 0D0
            DO IR = 2,IRMTIN
               AUX = AUX + R2RHO(IR,LM,IT,1)*RPWM3(IR)
     &               *DRDI_W_RADINT(IR,IM)
            END DO
C
            EFGLM(M,1) = (8.0D0*PI/5.0D0)*AUX
C
C---> use coulomb potential to determine extra atomic contribution
C
            EFGLM(M,2) = (V(IRMTIN,LM,IT)-(8.0D0*PI/5.0D0)*CMNTMTT(LM,IT
     &                   )/(RMTIN**3))/(RMTIN**2)
C
            EFGLM(M,3) = EFGLM(M,1) + EFGLM(M,2)
         END DO
C
C
         LOOP_I:DO I = 1,3
C-----------------------------------------------------------------------
C                convert the EFG to cartesian representation
C-----------------------------------------------------------------------
C---> xx
            EFG(1,1,I) = (EFGLM(2,I)*SQRT(3.0D0)-EFGLM(0,I))*FAC
C---> yy
            EFG(2,2,I) = -(EFGLM(2,I)*SQRT(3.0D0)+EFGLM(0,I))*FAC
C---> zz
            EFG(3,3,I) = EFGLM(0,I)*FAC*2.0D0
C---> xy
            EFG(1,2,I) = EFGLM(-2,I)*FAC*SQRT(3.0D0)
            EFG(2,1,I) = EFG(1,2,I)
C---> xz
            EFG(1,3,I) = EFGLM(1,I)*FAC*SQRT(3.0D0)
            EFG(3,1,I) = EFG(1,3,I)
C---> yz
            EFG(2,3,I) = EFGLM(-1,I)*FAC*SQRT(3.0D0)
            EFG(3,2,I) = EFG(2,3,I)
C
C-----------------------------------------------------------------------
C                     diagonalize the EFG tensor
C-----------------------------------------------------------------------
C
            EVEC(1:3,1:3) = EFG(1:3,1:3,I)
C
            CALL DSYEV('V','U',3,EVEC,3,EVAL(1,I),WORK,LWORK,INFO)
C
            DO N = 1,3
               RNORM = SQRT(EVEC(1,N)**2+EVEC(2,N)**2+EVEC(3,N)**2)
               DO M = 1,3
                  EVEC(M,N) = EVEC(M,N)/RNORM
               END DO
            END DO
C
C-----------------------------------------------------------------------
C                     convert to cgs units
C-----------------------------------------------------------------------
C
            EFG(1:3,1:3,I) = EFG(1:3,1:3,I)*EFG_AU2CGS
C
            EVAL(1:3,I) = EVAL(1:3,I)*EFG_AU2CGS
C
C-----------------------------------------------------------------------
C                               print out
C-----------------------------------------------------------------------
C
            WRITE (6,99005) TEXT(I)
            WRITE (6,99006) ((EFG(M,N,I),N=1,3),M=1,3)
C
            WRITE (6,99007)
            DO M = 1,3
               WRITE (6,99008) (EVEC(N,M),N=1,3),EVAL(M,I)/E0_CGS,
     &                         EVAL(M,I),EVAL(M,I)*EFG_CGS2SI
            END DO
C
            IF ( I.EQ.3 ) THEN
               IF ( ABS(EVAL(3,I)).GT.1.0D-7 ) THEN
                  ETA = ABS((EVAL(2,I)-EVAL(3,I))/EVAL(1,I))
                  WRITE (6,99009) ETA
               END IF
            END IF
         END DO LOOP_I
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
         IF ( WRBUILDBOT ) THEN
            WRITE (IFILBUILDBOT,99010) ROUTINE(1:LEN_TRIM(ROUTINE)),IT,
     &                                 ((EVEC(N,M),N=1,3),M=1,3)
            WRITE (IFILBUILDBOT,99011) ROUTINE(1:LEN_TRIM(ROUTINE)),IT,
     &                                 (((EFG(M,N,I),N=1,3),M=1,3),I=1,
     &                                 3)
         END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
      END DO
C
C ----------------------------------------------------------------------
99001 FORMAT (13X,'WARNING from <FPEFG>:',
     &        ' the charge density has to contain non spherical',
     &        ' contributions up to l=2 at least !')
99002 FORMAT (/,1X,79('*'),/,36X,'<FPEFG>',/,1X,79('*'),/)
99003 FORMAT (/,1X,79('-'))
99004 FORMAT (/,5X,'EFG for atom type IT=',I3,3X,A,/)
99005 FORMAT (/,5X,A,2X,'EFG tensor ',
     &        'w.r.t. global frame of reference',/)
99006 FORMAT (9X,1P,3E15.6,/,9X,1P,3E15.6,3X,'[V/cm^2]   (cgs)',/,9X,1P,
     &        3E15.6)
99007 FORMAT (/,14X,'principal axes',10X,
     &        'q [1/cm^3]   EFG [V/cm^2]    EFG [V/m^2]',/)
99008 FORMAT (8X,3F8.4,1X,1P,3E15.4)
99009 FORMAT (47X,'eta =',F8.5)
99010 FORMAT ('# BUILDBOT: ',A,': EFG  eigen vectors ','for IT =',I5,/,
     &        (1PE22.14))
99011 FORMAT ('# BUILDBOT: ',A,': EFG  tensors ','for IT =',I5,/,
     &        (1PE22.14))
C
      END
