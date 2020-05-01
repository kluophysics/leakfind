C*==me_ada_cor.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE ME_ADA_COR(NCST,NKPCOR,IKMCOR,GCOR,FCOR,IKMCB,IKMTF,
     &                      KCB,NSOLCB,GCB,FCB,IT,C,NPOL,ME,WGTR,
     &                      NCSTMAX,NPOLMAX)
C   ********************************************************************
C   *                                                                  *
C   * read wave function and calculate matrix elements for             *
C   *                                                                  *
C   *        <E_f, LAM_f | m c ->alfa * ->A_lam | E_i, LAM_i>          *
C   *                                                                  *
C   * according to ELECTRIC DIPOLE selection rules                     *
C   *                                                                  *
C   * using the         ALFA * A - FORM                                *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:C0,CI
      USE MOD_ANGMOM,ONLY:A1_ADA,A2_ADA,NKMMAX
      USE MOD_RMESH,ONLY:JRWS,NRMAX
      USE MOD_TYPES,ONLY:IMT
      IMPLICIT NONE
C*--ME_ADA_COR21
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 TOL
      PARAMETER (TOL=1D-10)
C
C Dummy arguments
C
      REAL*8 C
      INTEGER IKMTF,IT,KCB,NCST,NCSTMAX,NPOL,NPOLMAX,NSOLCB
      COMPLEX*16 FCB(NRMAX,2,NKMMAX),GCB(NRMAX,2,NKMMAX),
     &           ME(NCSTMAX,NKMMAX,NPOLMAX),WGTR(NRMAX)
      REAL*8 FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX)
      INTEGER IKMCB(2),IKMCOR(NCSTMAX,2),NKPCOR(NCSTMAX)
C
C Local variables
C
      COMPLEX*16 CINT1(NRMAX),CINT2(NRMAX),IMC,R1A1,R2A2,RME1(2,2),
     &           RME2(2,2)
      INTEGER I,ICST,IM,IPOL,IRTOP,JF,JI,KF,KI
C
C*** End of declarations rewritten by SPAG
C
      IMC = CI*0.5D0*C
C
      IM = IMT(IT)
      IRTOP = JRWS(IM)
C
      DO ICST = 1,NCST
         DO IPOL = 1,NPOL
            ME(ICST,IKMTF,IPOL) = C0
         END DO
      END DO
C
      DO ICST = 1,NCST
C
         DO KI = 1,NKPCOR(ICST)
            JI = IKMCOR(ICST,KI)
            DO KF = 1,NSOLCB
               JF = IKMCB(KF)
               DO IPOL = 1,NPOL
                  IF ( ABS(A1_ADA(JF,JI,IPOL)).GT.TOL ) GOTO 50
                  IF ( ABS(A2_ADA(JF,JI,IPOL)).GT.TOL ) GOTO 50
               END DO
            END DO
         END DO
C -------------------------------------- all angular matrix elements = 0
         CYCLE
C ---------------------------------- non-0 angular matrix elements found
C ------------------------------------- calculate radial matrix elements
 50      CONTINUE
         DO KI = 1,NKPCOR(ICST)
            DO KF = 1,NSOLCB
               DO I = 1,IRTOP
                  CINT1(I) = GCB(I,KF,KCB)*FCOR(I,KI,ICST)*WGTR(I)
                  CINT2(I) = FCB(I,KF,KCB)*GCOR(I,KI,ICST)*WGTR(I)
               END DO
C
               CALL CRADINT(IM,CINT1,RME1(KF,KI))
               CALL CRADINT(IM,CINT2,RME2(KF,KI))
C
            END DO
         END DO
C
C -------------------------------------- calculate total matrix elements
         DO KI = 1,NKPCOR(ICST)
            JI = IKMCOR(ICST,KI)
            DO KF = 1,NSOLCB
               JF = IKMCB(KF)
               DO IPOL = 1,NPOL
                  R1A1 = RME1(KF,KI)*A1_ADA(JF,JI,IPOL)
                  R2A2 = RME2(KF,KI)*A2_ADA(JF,JI,IPOL)
C
                  ME(ICST,IKMTF,IPOL) = ME(ICST,IKMTF,IPOL)
     &                                  + IMC*(R1A1-R2A2)
C
               END DO
            END DO
         END DO
C
      END DO
C
      END
