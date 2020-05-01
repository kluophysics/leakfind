C*==core_hff_split.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CORE_HFF_SPLIT(IM,RNUC,IRTOP,KAP1,KAP2,NSOL,MJ,GC,FC,
     &                          NRC,SHF,S,NMEHFMAX,NKMMAX,R,DRDI,SDIA,
     &                          SMDIA,SOFF,SMOFF,QDIA,QMDIA,QOFF,QMOFF,
     &                          JLIM)
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Calculates matrix elements of several hyperfine interaction
C     connected quantities in the core.
C     All the related arrays have a counting index as
C     the last index of the array indicates the corresponding physical
C     property.
C     Index-list
C     1      electron-Spin-electron-Spin Hyperfine field
C     2      nuclear-spin-electron-orbit hyperfine field
C     3      electron-spin-nulceus-spin-contact hyperfine field
C     4      expectation value of (1/r)^3
C     5      Total Hyperfine Field (see Rose (1961))
C     called by core
C
C PARAMETER definitions
C     BOHR-MAGNETON       IN ERG/GAUSS
C     CONVERSION FACTOR FOR HYPERFINE FIELDS FROM A.U. TO GAUSS
C     ELECTRON CHARGE     IN ESU
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      USE MOD_RMESH,ONLY:FINITE_NUCLEUS
      USE MOD_CONSTANTS,ONLY:A0_CGS,MB_CGS
      IMPLICIT NONE
C*--CORE_HFF_SPLIT28
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 F1,F2
      PARAMETER (F1=1.0D0,F2=2.0D0*MB_CGS/(A0_CGS*A0_CGS*A0_CGS))
C
C Dummy arguments
C
      INTEGER IM,IRTOP,JLIM,KAP1,KAP2,NKMMAX,NMEHFMAX,NRC,NSOL,S
      REAL*8 MJ,RNUC
      REAL*8 DRDI(NRC),FC(2,2,NRC),GC(2,2,NRC),QDIA(NKMMAX),
     &       QMDIA(NKMMAX),QMOFF(NKMMAX),QOFF(NKMMAX),R(NRC),
     &       SDIA(NKMMAX),SHF(2,2,NMEHFMAX),SMDIA(NKMMAX),SMOFF(NKMMAX),
     &       SOFF(NKMMAX)
C
C Local variables
C
      REAL*8 AMEHF(2,2),CFF(2,2),CFG(2,2),CGF(2,2),CGG(2,2),CQF(2,2),
     &       CQG(2,2),CSF(2,2),CSG(2,2),DOVR(NRC),DROVRN(NRC),
     &       DROVRN1(NRC),F(NRC,2),FF(2,2),FF1(2,2),FF2(2,2),FG(2,2),
     &       FG1(2,2),FG2(2,2),G(NRC,2),GF(2,2),GF1(2,2),GF2(2,2),
     &       GG(2,2),GG1(2,2),GG2(2,2)
      INTEGER I,IKM1,IKM2,J,K,K1,K2,N
      INTEGER IKAPMUE
C
C*** End of declarations rewritten by SPAG
C
      IF ( KAP2.EQ.0 ) KAP2 = KAP1
C
      CALL RINIT(4,GG)
      CALL RINIT(4,FF)
      CALL RINIT(4,GG1)
      CALL RINIT(4,FF1)
      CALL RINIT(4,GG2)
      CALL RINIT(4,FF2)
      CALL RINIT(4,GF)
      CALL RINIT(4,FG)
      CALL RINIT(4,GF1)
      CALL RINIT(4,FG1)
      CALL RINIT(4,GF2)
      CALL RINIT(4,FG2)
      CALL RINIT(2*NRC,G)
      CALL RINIT(2*NRC,F)
C
      DO K = 1,2
         DO N = 1,IRTOP
            G(N,K) = GC(K,S,N)
            F(N,K) = FC(K,S,N)
         END DO
      END DO
C     prepare meshes for finite nucleus calculation
      DO I = 1,IRTOP
         DOVR(I) = DRDI(I)/R(I)
         IF ( FINITE_NUCLEUS ) THEN
            DROVRN1(I) = (R(I)/RNUC)**3*DRDI(I)
            DROVRN(I) = DROVRN1(I)/R(I)
         END IF
      END DO
      IKM1 = IKAPMUE(KAP1,NINT(MJ-0.5D0))
      IKM2 = IKAPMUE(KAP2,NINT(MJ-0.5D0))
C     angular hyperfine matrix elements   see e.g.  E.M.Rose
C     the factor  i  has been omitted
      CALL RINIT(4,AMEHF)
      AMEHF(1,1) = 4.0D0*KAP1*MJ/(4.0D0*KAP1*KAP1-1.0D0)
      IF ( NSOL.EQ.2 ) THEN
         AMEHF(2,2) = 4.0D0*KAP2*MJ/(4.0D0*KAP2*KAP2-1.0D0)
         AMEHF(1,2) = DSQRT(0.25D0-(MJ/DBLE(KAP1-KAP2))**2)
         AMEHF(2,1) = AMEHF(1,2)
      END IF
C     coefficients for the spin-dipolar matrix elements
      CALL RINIT(4,CSF)
      CALL RINIT(4,CSG)
      CSG(1,1) = SDIA(IKM1)
      CSF(1,1) = SMDIA(IKM1)
      IF ( NSOL.EQ.2 ) THEN
         CSG(2,2) = SDIA(IKM2)
         CSG(1,2) = SOFF(IKM1)
         CSG(2,1) = CSG(1,2)
         CSF(2,2) = SMDIA(IKM2)
         CSF(1,2) = SMOFF(IKM1)
         CSF(2,1) = SMOFF(IKM1)
      END IF
C     COEFFICIENTS FOR THE QUADRUPOLAR MATRIX ELEMENTS
      CQG(1,1) = QDIA(IKM1)
      CQG(2,2) = QDIA(IKM2)
      CQG(1,2) = QOFF(IKM1)
      CQG(2,1) = CQG(1,2)
      CALL RINIT(4,CQF)
      CQF(1,1) = QMDIA(IKM1)
      CQF(2,2) = QMDIA(IKM2)
      CQF(1,2) = QMOFF(IKM1)
      CQF(2,1) = CQF(1,2)
C     coefficients to calculate the spin-spin field
      CALL RINIT(4,CGG)
      CALL RINIT(4,CGF)
      CGG(1,1) = -MJ/(KAP1+0.5D0)
      CGF(1,1) = -MJ/(-KAP1+0.5D0)
      IF ( NSOL.EQ.2 ) THEN
         CGG(1,2) = -DSQRT(1.0D0-(MJ/(KAP1+0.5D0))**2)
         CGG(2,1) = CGG(1,2)
         CGG(2,2) = -MJ/(KAP2+0.5D0)
         CGF(2,2) = -MJ/(-KAP2+0.5D0)
C     CGF(1,2) = -DSQRT( 1.0D0 - (MJ/(- KAP1+0.5D0))**2 )
C     CGF(2,1) = CGF(1,2)
      END IF
C     coefficients to calculate the orbital field
      CALL RINIT(4,CFG)
      CALL RINIT(4,CFF)
      CFG(1,1) = MJ*(KAP1+1.0D0)/(KAP1+0.5D0)
      CFF(1,1) = MJ*(-KAP1+1.0D0)/(-KAP1+0.5D0)
      IF ( NSOL.EQ.2 ) THEN
         CFG(2,2) = MJ*(KAP2+1.0D0)/(KAP2+0.5D0)
         CFG(1,2) = 0.5D0*DSQRT(1.0D0-(MJ/(KAP1+0.5D0))**2)
         CFG(2,1) = CFG(1,2)
         CFF(2,2) = MJ*(-KAP2+1.0D0)/(-KAP2+0.5D0)
C     CFF(1,2) = 0.5D0 * DSQRT( 1.0D0 - (MJ/(- KAP1 + 0.5D0))**2 )
C     CFF(2,1) = CFF(1,2)
      END IF
C Calculates integrals from 0.0 to jtop
      CALL HFFINT(IM,GG,G,G,DOVR,R,0.0D0,NSOL,IRTOP,NRC)
      CALL HFFINT(IM,FF,F,F,DOVR,R,0.0D0,NSOL,IRTOP,NRC)
      CALL HFFINT(IM,GF,G,F,DRDI,R,0.0D0,NSOL,IRTOP,NRC)
      CALL HFFINT(IM,FG,F,G,DRDI,R,0.0D0,NSOL,IRTOP,NRC)
      CALL RSUMUPINT(SHF(1,1,5),F1,GG,CQG,F1,FF,CQF,NSOL)
      IF ( FINITE_NUCLEUS ) THEN
C     calculates integrals inside nucleus at RNUC in order to get
C     contribution outside the nucleus
         CALL HFFINT(IM,GG1,G,G,DOVR,R,RNUC,NSOL,JLIM,NRC)
         CALL HFFINT(IM,FF1,F,F,DOVR,R,RNUC,NSOL,JLIM,NRC)
         CALL HFFINT(IM,GF1,G,F,DRDI,R,RNUC,NSOL,JLIM,NRC)
         CALL HFFINT(IM,FG1,F,G,DRDI,R,RNUC,NSOL,JLIM,NRC)
C     calculates contribution from RNUC to jtop
         DO I = 1,NSOL
            DO J = 1,NSOL
               GG(I,J) = GG(I,J) - GG1(I,J)
               FF(I,J) = FF(I,J) - FF1(I,J)
               GF(I,J) = GF(I,J) - GF1(I,J)
               FG(I,J) = FG(I,J) - FG1(I,J)
            END DO
         END DO
      END IF                    !end of nucleus.eq.0
C     calculates B_sp which is zero inside the nucleus
      CALL RSUMUPINT(SHF(1,1,1),F1,GG,CSG,-F1,FF,CSF,NSOL)
C     calculates hyperfine integrals from 0.0 to RNUC which are added
C     external integrals
      IF ( FINITE_NUCLEUS ) THEN
         CALL HFFINT(IM,GG2,G,G,DROVRN,R,RNUC,NSOL,JLIM,NRC)
         CALL HFFINT(IM,FF2,F,F,DROVRN,R,RNUC,NSOL,JLIM,NRC)
         CALL HFFINT(IM,GF2,G,F,DROVRN1,R,RNUC,NSOL,JLIM,NRC)
         CALL HFFINT(IM,FG2,F,G,DROVRN1,R,RNUC,NSOL,JLIM,NRC)
         DO I = 1,NSOL
            DO J = 1,NSOL
               GG(I,J) = GG(I,J) + GG2(I,J)
               FF(I,J) = FF(I,J) + FF2(I,J)
               GF(I,J) = GF(I,J) + GF2(I,J)
               FG(I,J) = FG(I,J) + FG2(I,J)
            END DO
         END DO
      END IF
C
C     calculates B_nseo and B_tot
C
      CALL RSUMUPINT(SHF(1,1,2),F2,GG,CFG,-F2,FF,CFF,NSOL)
C
C     modifications for B_ssc which is zero outside the nucleus
C
      IF ( FINITE_NUCLEUS ) THEN
         DO I = 1,NSOL
            DO J = 1,NSOL
               GG(I,J) = GG2(I,J)
               FF(I,J) = FF2(I,J)
            END DO
         END DO
      END IF
C
C     calculates B_ssc
C
      CALL RSUMUPINT(SHF(1,1,3),F2,GG,CGG,-F2,FF,CGF,NSOL)
C
C     for testing purposes write in 4 the sum of 1,2,3
      DO K1 = 1,NSOL
         DO K2 = 1,NSOL
            SHF(K1,K2,4) = SHF(K1,K2,1) + SHF(K1,K2,2) + SHF(K1,K2,3)
         END DO
      END DO
C
      END
C*==rsumupint.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
C
      SUBROUTINE RSUMUPINT(RSUM,VG,G,WG,VF,F,WF,N)
      IMPLICIT NONE
C*--RSUMUPINT236
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N
      REAL*8 VF,VG
      REAL*8 F(N,N),G(N,N),RSUM(N,N),WF(N,N),WG(N,N)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO J = 1,N
         DO I = 1,N
            RSUM(I,J) = VG*G(I,J)*WG(I,J) + VF*F(I,J)*WF(I,J)
         END DO
      END DO
      END
C*==hffint.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
C
      SUBROUTINE HFFINT(IM,GG,GA,GB,DR,R,RNUC,NSOL,IRTOP,NRC)
C     Calculates Hyperfine integrals, extrapolates to zero and
C     intrapolates to exact nuclear radius RNUC
      IMPLICIT NONE
C*--HFFINT275
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IM,IRTOP,NRC,NSOL
      REAL*8 RNUC
      REAL*8 DR(NRC),GA(NRC,2),GB(NRC,2),GG(2,2),R(NRC)
C
C Local variables
C
      INTEGER I,K1,K2
      REAL*8 X(5),Y(5),YI(NRC),ZI(NRC)
      REAL*8 YLAG
C
C*** End of declarations rewritten by SPAG
C
      DO K1 = 1,NSOL
         DO K2 = 1,NSOL
C
            DO I = 1,IRTOP
               YI(I) = GA(I,K1)*GB(I,K2)*DR(I)
            END DO
C
            CALL RRADINT_R(IM,YI,ZI)
C
C     Intrapolation
            IF ( ABS(RNUC).GT.1D-16 ) THEN
               DO I = 1,5
                  X(I) = R(IRTOP-5+I)
                  Y(I) = ZI(IRTOP-5+I)
               END DO
               ZI(IRTOP) = YLAG(RNUC,X,Y,0,4,5)
            END IF
C     Extrapolation
            X(1) = 1.0D0
            X(2) = 6.0D0
            X(3) = 11.0D0
            Y(1) = ZI(IRTOP) - ZI(1)
            Y(2) = ZI(IRTOP) - ZI(5)
            Y(3) = ZI(IRTOP) - ZI(9)
            GG(K1,K2) = YLAG(0.0D0,X,Y,0,2,3)
         END DO
      END DO
C
      END
C
