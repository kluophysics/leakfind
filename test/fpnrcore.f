C*==fpnrcore.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPNRCORE(EMIN,ECOR_LT,NLCORE,BCOR,BCORS)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_RMESH,ONLY:NRMAX,JRCUT,DRDI,R,R2DRDI
      USE MOD_TYPES,ONLY:NT,NTMAX,NLMFPMAX,RHO2NS,IMT,BT,VT,NCORT,ITBOT,
     &    ITTOP,Z
      USE MOD_ANGMOM,ONLY:NSPIN
      USE MOD_CONSTANTS,ONLY:SQRT_4PI
      IMPLICIT NONE
C*--FPNRCORE14
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 EMIN
      INTEGER NLCORE
      REAL*8 BCOR(NTMAX),BCORS(NTMAX),ECOR_LT(NLCORE,NTMAX)
C
C Local variables
C
      REAL*8 ANC,AUX,CORLVL(18),RHO(:),RORG(20,2),SPNWGT,V1(:),W,WK(:,:)
      INTEGER I,IM,IR,IRTOP,IS,ISR,IT,MATCH(18)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE V1,WK,RHO
C
      ALLOCATE (V1(NRMAX),RHO(NRMAX),WK(NRMAX,5))
C
      CALL RINIT(NRMAX*NLMFPMAX*NTMAX*3,RHO2NS)
C
      IF ( NLCORE.GT.4 ) STOP 'in <FPNRCORE>: NLCORE > 4'
      IF ( IREL.GT.2 ) STOP 'in <FPNRCORE>: IREL > 2'
C
      BCOR(1:NTMAX) = 0D0
      BCORS(1:NTMAX) = 0D0
      ECOR_LT(1:NLCORE,1:NTMAX) = 0.0D0
      IF ( IREL.LE.0 ) THEN
         ISR = 0
      ELSE
         ISR = 1
      END IF
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = ITBOT,ITTOP
C
         IM = IMT(IT)
         IRTOP = JRCUT(1,IM) + 1
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
         DO IS = 1,NSPIN
C
            IF ( IREL.EQ.2 ) THEN
               SPNWGT = NINT((IS-1.5D0)*2D0)
            ELSE
               SPNWGT = 0D0
            END IF
C
C-----------------------------------------------------------------------
C              set up spin-dependent spherical potential
C-----------------------------------------------------------------------
C
            DO IR = 1,IRTOP - 1
               V1(IR) = VT(IR,IT) + SPNWGT*BT(IR,IT)
            END DO
            V1(IRTOP) = 0D0
C
            MATCH(1:18) = 0
            CORLVL(1:18) = -10000.0D0
C
            CALL CSTATE(Z(IT),ECOR_LT(1,IT),NCORT(IT),V1,RHO,RORG(1,IS),
     &                  CORLVL,ANC,DRDI(1,IM),R(1,IM),MATCH,IRTOP,WK,
     &                  ISR,EMIN)
C
            DO IR = 1,JRCUT(1,IM)
               RHO2NS(IR,1,IT,1) = RHO2NS(IR,1,IT,1) + RHO(IR)
               RHO2NS(IR,1,IT,2) = RHO2NS(IR,1,IT,2) + SPNWGT*RHO(IR)
            END DO
C
            IF ( ABS(ANC-NCORT(IT)/2D0).GT.0.2D0 )
     &            STOP 'in <FPNRCORE>: ANC'
         END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
         IF ( NSPIN.EQ.1 ) THEN
            ECOR_LT(1:NLCORE,IT) = 2*ECOR_LT(1:NLCORE,IT)
            DO IR = 1,JRCUT(1,IM)
               RHO2NS(IR,1,IT,1) = 2*RHO2NS(IR,1,IT,1)
               RHO2NS(IR,1,IT,2) = 0D0
            END DO
         END IF
C
         DO IR = 1,IRTOP
            WK(IR,1) = RHO2NS(IR,1,IT,1)*R2DRDI(IR,IM)
            WK(IR,2) = RHO2NS(IR,1,IT,2)*R2DRDI(IR,IM)
            RHO2NS(IR,1,IT,1) = RHO2NS(IR,1,IT,1)/SQRT_4PI
            RHO2NS(IR,1,IT,2) = RHO2NS(IR,1,IT,2)/SQRT_4PI
         END DO
C
         CALL RRADINT(IM,WK(1,1),AUX)
         WRITE (6,99001) 'charge',IT,AUX
C
         CALL RRADINT(IM,WK(1,2),AUX)
         WRITE (6,99001) ' spin ',IT,AUX
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      IF ( NT.EQ.0 ) THEN
         OPEN (91,FILE='RHOCOR-NR.dat')
         IT = 1
         IM = IMT(IT)
         DO I = 1,JRCUT(1,IM)
C
C ------- account for factor SQRT_4PI between RHOCHR (RHOSPN) and RHO2NS
C
            W = R(I,IM)**2*SQRT_4PI
C
            WRITE (91,'(12E14.5)') R(I,IM),
     &                             (RHO2NS(I,1,IT,1)*W,RHO2NS(I,1,IT,2)
     &                             *W,IT=ITBOT,ITTOP)
         END DO
         CLOSE (91)
         WRITE (6,*) 'core charge density written to RHOCOR-NR.dat '
      END IF
C
      DEALLOCATE (RHO)
C
      BCORS(1:NT) = BCOR(1:NT)
C
99001 FORMAT (10X,'integrated core ',A,' density for atom type ',I2,':',
     &        F12.8)
C
C
      END
C*==cstate.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE CSTATE(Z,ECOR_LT,NCORT,V,RO,RORG,CORLVL,ANC,DR,XR,
     &                  MATCH,MESHR,WK,ISR,EBTM)
C----------------------------------------------------------------------
C     Calculate core state by matching boundary condition.
C     coded by H.Akai, 1983, Juelich
C----------------------------------------------------------------------
      IMPLICIT NONE
C*--CSTATE160
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ANC,EBTM
      INTEGER ISR,MESHR,NCORT,Z
      REAL*8 CORLVL(18),DR(MESHR),ECOR_LT(4),RO(MESHR),RORG(20),V(MESHR)
     &       ,WK(MESHR,5),XR(MESHR)
      INTEGER MATCH(18)
C
C Local variables
C
      REAL*8 CONFIG(18),DLT,E,EB,EMAX,EMIN,G1,G2,RIN,SMALL,TINT,TOL
      LOGICAL INIT
      INTEGER ISTOP,ITR,J,JJ,K,L(18),LCST,NN,NODE,NPQ(18)
C
C*** End of declarations rewritten by SPAG
C
      DATA SMALL/2D-2/
C                1s 2s 2p 3s 3p 3d 4s 4p 4d 5s 5p 4f 5d 6s 6p 5f 6d 7s
      DATA NPQ/1,2,2,3,3,3,4,4,4,5,5,4,5,6,6,5,6,7/,L/0,0,1,0,1,2,0,1,2,
     &     0,1,3,2,0,1,3,2,0/
C
      DATA ISTOP/10000/,TOL/1D-8/,EB/ - 20D0/
C
      NN = 0
      DO J = 1,18
         JJ = 2*L(J) + 1
         NN = NN + JJ
         IF ( NN.LE.NCORT/2 ) THEN
            CONFIG(J) = JJ
         ELSE
            CONFIG(J) = 0
         END IF
C
         CORLVL(J) = -Z*Z/(1.5D0*NPQ(J)*NPQ(J))
C
      END DO
C
      RO(1:MESHR) = 0D0
C
      ANC = 0D0
      TINT = 0D0
C
      DO J = 1,18
         IF ( ABS(CONFIG(J)).GT.0D0 ) THEN
C
            LCST = L(J)
            JJ = L(J) + 1
            NODE = NPQ(J) - JJ
            E = CORLVL(J)
            INIT = MATCH(J).LT.1
            DO K = 1,MESHR
               WK(K,3) = V(K)
            END DO
C
            EMAX = 1D10
            EMIN = -1D10
            DO ITR = 1,ISTOP
C
               CALL CORADA(E,JJ,WK(1,1),RIN,MATCH(J),G1,G2,NN,WK(1,3),
     &                     DR,XR,MESHR,ISR)
C
               IF ( NN.GT.NODE ) THEN
                  EMAX = MIN(E+EB,EMAX)
                  E = MAX(EMAX*1.25D0,(EMAX+EMIN)*5D-1) - EB
                  IF ( INIT ) MATCH(J) = 0
                  CYCLE
               END IF
               IF ( NN.LT.NODE ) THEN
                  EMIN = MAX(E+EB,EMIN)
                  E = MIN(EMIN*0.75D0,(EMAX+EMIN)*5D-1) - EB
                  IF ( INIT ) MATCH(J) = 0
                  CYCLE
               END IF
               DLT = -WK(MATCH(J),1)*(G1-G2)
               IF ( ABS(DLT).LT.TOL ) GOTO 20
               E = E + DLT*0.1D0
            END DO
            WRITE (*,'(a,i3)') 
     &                     '   ***err in cstate...no convergence for j='
     &                     ,J
            WRITE (*,'(1x,6f12.5)') (WK(K,1),K=1,MESHR,10)
            STOP
C
 20         CONTINUE
            DO K = 1,MESHR
               WK(K,4) = V(K) - WK(K,3)
               WK(K,1) = WK(K,1)/XR(K)**2
               WK(K,2) = WK(K,1)
            END DO
C
            CORLVL(J) = E
            RORG(J) = 0D0
            ECOR_LT(1+LCST) = ECOR_LT(1+LCST) + E*CONFIG(J)
            WRITE (*,*) 'ECORE ',J,E
C
            IF ( RIN.LT.0.98D0 ) WRITE (6,*)
     &                                   'WARNING from <CSTATE>:   RIN='
     &                                  ,RIN,'  for J=',J
C
            IF ( E.LT.EBTM ) THEN
               CONFIG(J) = ABS(CONFIG(J))
C
               ANC = ANC + CONFIG(J)*RIN
               TINT = TINT + CONFIG(J)*(1D0-RIN)
               DO K = 1,MESHR - 1
C                 RO(K) = RO(K) + CONFIG(J)*WK(K,1)*XR(K)**2
                  WK(K,1) = WK(K,1)/RIN
                  RO(K) = RO(K) + CONFIG(J)*WK(K,1)
               END DO
               CALL EXTORG(RORG(J),WK,XR)
               IF ( ABS(RORG(J)).LT.1D-20 ) RORG(J) = 0D0
               RORG(J) = CONFIG(J)*RORG(J)
C
            ELSE
               CONFIG(J) = -ABS(CONFIG(J))
            END IF
            IF ( ABS(E-EBTM).LT.SMALL )
     &            STOP 'in <CSTATE>: core level near ebtm found'
         END IF
      END DO
      RO(MESHR) = TINT
      END
C*==corada.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE CORADA(E,J,R,RIN,KMATCH,G1,G2,NODE,V,DR,XR,MESHR,ISR)
C----------------------------------------------------------------------
C     Solves schroedinger equation for core states.
C     coded by H.Akai, 1983, Juelich
C----------------------------------------------------------------------
      IMPLICIT NONE
C*--CORADA304
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 C1,C2,C3,C4,C5,C6
      PARAMETER (C1=274.0720442D0,C2=1D0/C1**2,C3=24D0/9D0,C4=19D0/9D0,
     &           C5=5D0/9D0,C6=1D0/9D0)
C
C Dummy arguments
C
      REAL*8 E,G1,G2,RIN
      INTEGER ISR,J,KMATCH,MESHR,NODE
      REAL*8 DR(MESHR),R(MESHR),V(MESHR),XR(MESHR)
C
C Local variables
C
      REAL*8 BL,C,DELT,DLT,DP(3),DP1,DP2,DP3,DQ(3),DQ1,DQ2,DQ3,EE,EK2,
     &       ER,EREL,ES,GM,H(10),P,P1,P2,P3,PM,PMATCH,PS,Q,Q1,Q2,Q3,QM,
     &       QMATCH,RED,RED1,RED2,RM,RM1,RM2,RM3,S,S1,S2,S3,SMALL,XZERO,
     &       ZERO
      INTEGER I,ITR,ITRMX,K,KB,KBACK,KK,L
      LOGICAL SN1,SN2,SRA
C
C*** End of declarations rewritten by SPAG
C
      DATA C/C1/,ITRMX/5/,SMALL/1D-2/
      SRA = ISR.EQ.1
      L = J - 1
      BL = DBLE(L*J)
C     ES = V(MESHR) + E
      ES = E
      DO K = 1,MESHR - 1
         KB = K
         IF ( V(K)*XR(K)**2+BL.LT.0D0 ) EXIT
      END DO
C
C     ---- expansion of type p=x**(l+1)+... ---
      DO K = 1,3
         EK2 = (V(K)-ES)*XR(K)
         RM = XR(K)
         IF ( SRA ) RM = RM - EK2*C2
         EK2 = BL/RM + EK2
         DELT = DR(K)/XR(K)
         P = XR(K)**J
         Q = DBLE(L)*P/RM
         DP(K) = (P+RM*Q)*DELT/C3
         DQ(K) = (EK2*P-Q)*DELT/C3
         R(K) = P**2
         IF ( SRA ) R(K) = R(K)*(1D0+BL/(RM*C)**2) + (Q/C)**2
      END DO
      S1 = 5D-1*R(1)*DR(1) + 2D0*R(2)*DR(2) + R(3)*DR(3)
      NODE = 0
      SN1 = P.GT.0D0
      IF ( SRA ) THEN
C
C
C     ---- adams-multon for sra case ----
C              ( forward )
         DO K = 4,MESHR - 4
            EK2 = (V(K)-ES)*XR(K)
            RM = XR(K) - EK2*C2
            EK2 = BL/RM + EK2
            DLT = C3*XR(K)/DR(K)
            GM = (1D0-DLT)/EK2
            PM = P + C4*DP(3) - C5*DP(2) + C6*DP(1)
            QM = Q + C4*DQ(3) - C5*DQ(2) + C6*DQ(1)
            Q = (GM*QM-PM)/(RM+GM*(DLT+1D0))
            P = (PM+RM*Q)/(DLT-1D0)
            DP(1) = DP(2)
            DP(2) = DP(3)
            DQ(1) = DQ(2)
            DQ(2) = DQ(3)
            DP(3) = P + RM*Q
            DQ(3) = EK2*P - Q
            P = P*DLT
            Q = Q*DLT
            SN2 = P.GT.0D0
            IF ( SN1 .NEQV. SN2 ) NODE = NODE + 1
            SN1 = SN2
            IF ( KMATCH.LT.1 .AND. K.GT.KB .AND. EK2.GT.0D0 ) GOTO 100
            IF ( K.EQ.KMATCH ) GOTO 100
            R(K) = P**2*(1D0+BL/(RM*C)**2) + (Q/C)**2
            PS = DR(K)*R(K)
            S1 = S1 + PS
            IF ( MOD(K,2).EQ.0 ) S1 = S1 + PS
         END DO
      ELSE
C
C     ---- adams-multon for nrl case ----
C              ( forward )
         DO K = 4,MESHR - 4
            EK2 = (V(K)-ES)*XR(K)
            RM = XR(K)
            EK2 = BL/RM + EK2
            DLT = C3*XR(K)/DR(K)
            GM = (1D0-DLT)/EK2
            PM = P + C4*DP(3) - C5*DP(2) + C6*DP(1)
            QM = Q + C4*DQ(3) - C5*DQ(2) + C6*DQ(1)
            Q = (GM*QM-PM)/(RM+GM*(DLT+1D0))
            P = (PM+RM*Q)/(DLT-1D0)
            DP(1) = DP(2)
            DP(2) = DP(3)
            DQ(1) = DQ(2)
            DQ(2) = DQ(3)
            DP(3) = P + RM*Q
            DQ(3) = EK2*P - Q
            P = P*DLT
            Q = Q*DLT
            SN2 = P.GT.0D0
            IF ( SN1 .NEQV. SN2 ) NODE = NODE + 1
            SN1 = SN2
            IF ( KMATCH.LT.1 .AND. K.GT.KB .AND. EK2.GT.0D0 ) GOTO 100
            IF ( K.EQ.KMATCH ) GOTO 100
            R(K) = P**2
            PS = DR(K)*R(K)
            S1 = S1 + PS
            IF ( MOD(K,2).EQ.0 ) S1 = S1 + PS
         END DO
      END IF
C
      S1 = S1 - 2D0*PS
      K = MESHR - 4
 100  CONTINUE
      KMATCH = K
      PMATCH = P
      QMATCH = Q
C
      XZERO = SQRT(DABS(2D3/ES))
      DO K = KMATCH + 3,MESHR - 1,4
         KBACK = K + 1 - MOD(K,2)
         IF ( XR(K).GT.XZERO ) GOTO 200
      END DO
      KBACK = MESHR - 1
C
C     --- following code gives the natural boundary condition, which is
C         suitable for a single muffin-tin potential in a free space.
C         This, however, causes some problem for positive energies.
C         In order to overcome the difficulty, I introduce somewhat
C         artificcial boundary condition for positive, or negative
C         but nearly zero, energies.
 200  CONTINUE
      EE = ES
      EREL = EE
      IF ( SRA ) EREL = EE + (EE/C)**2
      EREL = 5.D-1*(EREL-SQRT(EREL**2+SMALL**2))
      ER = XR(KBACK)*SQRT(-EREL)
      H(1) = EXP(-ER)/ER
      H(2) = -H(1)
      DO I = 1,J
         H(I+2) = H(I) - DBLE(2*I-1)*H(I+1)/ER
      END DO
      P1 = XR(KBACK)*H(J+1)
      Q1 = ER*(H(J)-DBLE(J)*H(J+1)/ER)*EE/EREL
C     ----- next statement gives normalization in the entire space---
      S3 = 0D0
      IF ( KBACK.GE.MESHR-1 ) S3 = 5D-1*XR(MESHR-1)
     &                             **3*(H(J)*H(J+2)-H(J+1)**2)
C     --- normalization special
C     s3=0d0
C     p1=0d0
C     q1=0d0
      EK2 = (V(KBACK)-ES)*XR(KBACK)
      RM1 = XR(KBACK)
      IF ( SRA ) RM1 = RM1 - EK2*C2
      EK2 = BL/RM1 + EK2
C
C     ---- self-starting formula ----
C              ( backward )
      DELT = -DR(KBACK)/XR(KBACK)
      DP1 = (P1+RM1*Q1)*DELT
      DQ1 = (EK2*P1-Q1)*DELT
      P2 = P1 + DP1
      Q2 = Q1 + DQ1
      EK2 = (V(KBACK-1)-ES)*XR(KBACK-1)
      RM2 = XR(KBACK-1)
      IF ( SRA ) RM2 = RM2 - EK2*C2
      EK2 = BL/RM2 + EK2
      DELT = -DR(KBACK-1)/XR(KBACK-1)
      DP2 = (P2+RM2*Q2)*DELT
      DQ2 = (EK2*P2-Q2)*DELT
      P2 = P1 + 5D-1*(DP1+DP2)
      Q2 = Q1 + 5D-1*(DQ1+DQ2)
      P3 = P1 + 2D0*DP2
      Q3 = Q1 + 2D0*DQ2
      DO ITR = 1,ITRMX
         EK2 = (V(KBACK-1)-ES)*XR(KBACK-1)
         RM2 = XR(KBACK-1)
         IF ( SRA ) RM2 = RM2 - EK2*C2
         EK2 = BL/RM2 + EK2
         DELT = -DR(KBACK-1)/XR(KBACK-1)
         DP2 = (P2+RM2*Q2)*DELT
         DQ2 = (EK2*P2-Q2)*DELT
         EK2 = (V(KBACK-2)-ES)*XR(KBACK-2)
         RM3 = XR(KBACK-2)
         IF ( SRA ) RM3 = RM3 - EK2*C2
         EK2 = BL/RM3 + EK2
         DELT = -DR(KBACK-2)/XR(KBACK-2)
         DP3 = (P3+RM3*Q3)*DELT
         DQ3 = (EK2*P3-Q3)*DELT
         P2 = P1 + (5D0*DP1+8D0*DP2-DP3)/12D0
         Q2 = Q1 + (5D0*DQ1+8D0*DQ2-DQ3)/12D0
         P3 = P1 + (DP1+4D0*DP2+DP3)/3D0
         Q3 = Q1 + (DQ1+4D0*DQ2+DQ3)/3D0
      END DO
C
      IF ( SRA ) THEN
         R(KBACK) = P1**2*(1D0+BL/(RM1*C)**2) + (Q1/C)**2
         R(KBACK-1) = P2**2*(1D0+BL/(RM2*C)**2) + (Q2/C)**2
         R(KBACK-2) = P3**2*(1D0+BL/(RM3*C)**2) + (Q3/C)**2
      ELSE
         R(KBACK) = P1**2
         R(KBACK-1) = P2**2
         R(KBACK-2) = P3**2
      END IF
      DP(1) = DP1/C3
      DP(2) = DP2/C3
      DP(3) = DP3/C3
      DQ(1) = DQ1/C3
      DQ(2) = DQ2/C3
      DQ(3) = DQ3/C3
      P = P3
      Q = Q3
      S2 = 5D-1*R(KBACK)*DR(KBACK) + 2D0*R(KBACK-1)*DR(KBACK-1)
     &     + R(KBACK-2)*DR(KBACK-2)
C
C     ---- adams-multon ----
C          ( backward )
      IF ( SRA ) THEN
C
C     ---- adams-multon for sra ----
C             ( backward )
         DO KK = 4,KBACK - KMATCH + 1
            K = KBACK + 1 - KK
            EK2 = (V(K)-ES)*XR(K)
            RM = XR(K) - EK2*C2
            EK2 = BL/RM + EK2
            DLT = -C3*XR(K)/DR(K)
            GM = (1D0-DLT)/EK2
            PM = P + C4*DP(3) - C5*DP(2) + C6*DP(1)
            QM = Q + C4*DQ(3) - C5*DQ(2) + C6*DQ(1)
            Q = (GM*QM-PM)/(RM+GM*(DLT+1D0))
            P = (PM+RM*Q)/(DLT-1D0)
            DP(1) = DP(2)
            DP(2) = DP(3)
            DQ(1) = DQ(2)
            DQ(2) = DQ(3)
            DP(3) = P + RM*Q
            DQ(3) = EK2*P - Q
            P = P*DLT
            Q = Q*DLT
            R(K) = P**2*(1D0+BL/(RM*C)**2) + (Q/C)**2
            PS = DR(K)*R(K)
            S2 = S2 + PS
            IF ( MOD(K,2).EQ.0 ) S2 = S2 + PS
         END DO
      ELSE
C
C     ---- adams-multon for nrl ----
C             ( backward )
         DO KK = 4,KBACK - KMATCH + 1
            K = KBACK + 1 - KK
            EK2 = (V(K)-ES)*XR(K)
            RM = XR(K)
            EK2 = BL/RM + EK2
            DLT = -C3*XR(K)/DR(K)
            GM = (1D0-DLT)/EK2
            PM = P + C4*DP(3) - C5*DP(2) + C6*DP(1)
            QM = Q + C4*DQ(3) - C5*DQ(2) + C6*DQ(1)
            Q = (GM*QM-PM)/(RM+GM*(DLT+1D0))
            P = (PM+RM*Q)/(DLT-1D0)
            DP(1) = DP(2)
            DP(2) = DP(3)
            DQ(1) = DQ(2)
            DQ(2) = DQ(3)
            DP(3) = P + RM*Q
            DQ(3) = EK2*P - Q
            P = P*DLT
            Q = Q*DLT
            R(K) = P**2
            PS = DR(K)*R(K)
            S2 = S2 + PS
            IF ( MOD(K,2).EQ.0 ) S2 = S2 + PS
         END DO
      END IF
C
      G1 = Q/P
      G2 = QMATCH/PMATCH
      RED = (PMATCH/P)**2
      IF ( S2*RED/S1.LT.1D-18 ) RED = 0D0
      S3 = S3*RED
      S = 2D0*(S1+S2*RED)/3D0 + S3
      RIN = 1D0 - S3/S
      RED1 = 1D0/S
      RED2 = RED1*RED
      ZERO = 1D-36/(ABS(RED2)+1D-36)
      DO K = 1,KMATCH - 1
         R(K) = R(K)*RED1
      END DO
      DO K = KMATCH,KBACK
         IF ( ABS(R(K)).LT.ZERO ) R(K) = 0D0
         R(K) = R(K)*RED2
      END DO
      DO K = KBACK + 1,MESHR
         R(K) = 0D0
      END DO
      END
C*==extorg.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE EXTORG(ORG,A,X)
C-----------------------------------------------------------------------
C     Extrapolate charge density toward the origin.
C     coded by H.Akai, 1983, Juelich
C-----------------------------------------------------------------------
      IMPLICIT NONE
C*--EXTORG633
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ORG
      REAL*8 A(3),X(3)
C
C Local variables
C
      REAL*8 P,Q
C
C*** End of declarations rewritten by SPAG
C
      P = ((A(1)-A(2))/(X(1)-X(2))-(A(3)-A(2))/(X(3)-X(2)))/(X(1)-X(3))
      Q = (A(1)-A(2))/(X(1)-X(2)) - P*(X(1)-X(2))
      ORG = P*X(2)**2 - Q*X(2) + A(2)
      END
