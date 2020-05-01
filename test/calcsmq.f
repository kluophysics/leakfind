C*==calcsmq.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE CALCSMQ(SDIA,SMDIA,SOFF,SMOFF,QDIA,QMDIA,QOFF,QMOFF)
C   ********************************************************************
C   *                                                                  *
C   *   Q-COEFFICIENTS    -  USED TO CALCULATE THE EFG                 *
C   *                                                                  *
C   *   Q(K,K',MUE) =                                                  *
C   *      (-1)**(MUE-1/2)*SQRT*[(2L+1)(2L'+1)]/4PI * C(L,L',2;0,0)    *
C   *      * [ C(K,MUE,2) * C(K',MUE,2) * C(L,L',2;MUE+1/2,-MUE-1/2)   *
C   *        - C(K,MUE,2) * C(K',MUE,2) * C(L,L',2;MUE+1/2,-MUE-1/2)]  *
C   *                                                                  *
C   *   S-COEFFICIENTS    -  USED TO CALCULATED THE SPIN-MAGN          *
C   *                        HYPERFINE TERM                            *
C   *                                                                  *
C   *   WITH L = L'  AND  K=K' / -K'-1                                 *
C   *                                                                  *
C   *   ..DIA/..OFF ARE THE ELEMENTS FOR  K=K'/K=-K'-1                 *
C   *   I  NUMBERS THE  Q's   and   S's   COLUMN-WISE                  *
C   *   STARTING WITH COLUMN  1                                        *
C   *                                                                  *
C   *   CHECK-OPTION:    PRINT Q-COEFFICIENTS                          *
C   *                                                                  *
C   *   THE PREFACTORS FOR  Q  AND  S  ARE INCLUDED                    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:CGC,NLMAX,NKMMAX
      USE MOD_CONSTANTS,ONLY:A0_CGS,MB_CGS,PI
      IMPLICIT NONE
C*--CALCSMQ30
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 VFSHF,QCGS
      PARAMETER (VFSHF=2.0D0*MB_CGS/(A0_CGS*A0_CGS*A0_CGS),
     &           QCGS=(2*4*PI/5.0D0)/A0_CGS**3)
C
C Dummy arguments
C
      REAL*8 QDIA(NKMMAX),QMDIA(NKMMAX),QMOFF(NKMMAX),QOFF(NKMMAX),
     &       SDIA(NKMMAX),SMDIA(NKMMAX),SMOFF(NKMMAX),SOFF(NKMMAX)
C
C Local variables
C
      REAL*8 C3DN,C3UP,CC,CM3DN,CM3UP,J(2),L,LB(2),MJ,Q,QM,QQ1,QQ2,S,SM,
     &       VF,VFM
      REAL*8 CGC_RACAH,GAUNT_CYLM
      LOGICAL CHECK
      INTEGER I,IKM1,IKM2,IL,IMKM1,IMKM2,K1,K2,KAP(2),MUEM05,N,N2,NSOL
      INTEGER IKAPMUE
C
C*** End of declarations rewritten by SPAG
C
C PREFACTOR FOR THE SPIN MAGN HF - TERM
C PREFACTOR FOR THE QUADRUPOLAR  PARAMETER  Q   (IN 1/CM**3)
C
      DO I = 1,NKMMAX
         QDIA(I) = 0.0D0
         QMDIA(I) = 0.0D0
         QOFF(I) = 0.0D0
         QMOFF(I) = 0.0D0
         SDIA(I) = 0.0D0
         SMDIA(I) = 0.0D0
         SOFF(I) = 0.0D0
         SMOFF(I) = 0.0D0
      END DO
C
      CHECK = .TRUE.
C
      IF ( CHECK ) WRITE (6,99004) NLMAX
C
      DO IL = 1,NLMAX
         L = IL - 1
C
         DO MUEM05 = -IL,(IL-1)
            MJ = 1.0D0*MUEM05 + 0.5D0
C
            CC = (-1.0D0)**MUEM05/(4*PI)
C
            VF = CC*(2*L+1)*CGC_RACAH(L,L,2.0D0,0.0D0,0.0D0,0.0D0)
            C3UP = CGC_RACAH(L,L,2.0D0,(MJ-0.5D0),(-MJ+0.5D0),0.0D0)
            C3DN = CGC_RACAH(L,L,2.0D0,(MJ+0.5D0),(-MJ-0.5D0),0.0D0)
C
            KAP(1) = NINT(-L-1)
            J(1) = L + 0.5D0
            LB(1) = L + 1.0D0
            KAP(2) = NINT(+L)
            J(2) = L - 0.5D0
            LB(2) = L - 1.0D0
C
            IF ( ABS(MJ).GT.L ) THEN
               NSOL = 1
            ELSE
               NSOL = 2
            END IF
C
            DO K2 = 1,NSOL
               IKM2 = IKAPMUE(KAP(K2),MUEM05)
               IMKM2 = IKAPMUE(-KAP(K2),MUEM05)
               I = IKM2
               IF ( I.GT.NKMMAX ) THEN
                  WRITE (6,99001) I,NKMMAX
                  STOP
               END IF
C
               DO K1 = 1,NSOL
                  IKM1 = IKAPMUE(KAP(K1),MUEM05)
                  IMKM1 = IKAPMUE(-KAP(K1),MUEM05)
C
                  CM3UP = CGC_RACAH(LB(K1),LB(K2),2.0D0,(MJ-0.5D0),
     &                    (-MJ+0.5D0),0.0D0)
                  CM3DN = CGC_RACAH(LB(K1),LB(K2),2.0D0,(MJ+0.5D0),
     &                    (-MJ-0.5D0),0.0D0)
C
                  VFM = CC*DSQRT((2*LB(K1)+1)*(2*LB(K2)+1))
     &                  *CGC_RACAH(LB(K1),LB(K2),2.0D0,0.0D0,0.0D0,
     &                  0.0D0)
C
                  Q = VF*(CGC(IKM1,2)*CGC(IKM2,2)*C3UP-CGC(IKM1,1)
     &                *CGC(IKM2,1)*C3DN)
C
                  QM = VFM*(CGC(IMKM1,2)*CGC(IMKM2,2)*CM3UP-CGC(IMKM1,1)
     &                 *CGC(IMKM2,1)*CM3DN)
C
                  QQ1 = CGC(IKM1,2)*CGC(IKM2,2)
     &                  *GAUNT_CYLM(NINT(L),MUEM05,NINT(L),MUEM05,2,0)
     &                  + CGC(IKM1,1)*CGC(IKM2,1)
     &                  *GAUNT_CYLM(NINT(L),MUEM05+1,NINT(L),MUEM05+1,2,
     &                  0)
                  QQ2 = CGC(IMKM1,2)*CGC(IMKM2,2)
     &                  *GAUNT_CYLM(NINT(LB(K1)),MUEM05,NINT(LB(K2)),
     &                  MUEM05,2,0) + CGC(IMKM1,1)*CGC(IMKM2,1)
     &                  *GAUNT_CYLM(NINT(LB(K1)),MUEM05+1,NINT(LB(K2)),
     &                  MUEM05+1,2,0)
                  IF ( ABS(QQ1-QQ2).GT.1D-8 ) THEN
                     WRITE (6,*) ' QQ1 <> QQ2 ',QQ1,QQ2
                     WRITE (6,*) ' L,K1,K2,MUEM05 ',L,K1,K2,MUEM05
                  END IF
C
C     no off-diagonal elements (corrected) mb 12.1.93
C
                  S = -DSQRT(4.0D0*PI/5.0D0)
     &                *(CGC(IKM1,1)*CGC(IKM2,1)*GAUNT_CYLM(NINT(L),
     &                NINT(MJ+0.5D0),2,0,NINT(L),NINT(MJ+0.5D0))
     &                -CGC(IKM1,2)*CGC(IKM2,2)
     &                *GAUNT_CYLM(NINT(L),NINT(MJ-0.5D0),2,0,NINT(L),
     &                NINT(MJ-0.5D0)))
                  S = S + 3.0D0*DSQRT(2.0D0*PI/15.0D0)
     &                *(-CGC(IKM1,1)*CGC(IKM2,2)
     &                *GAUNT_CYLM(NINT(L),NINT(MJ+0.5D0),2,1,NINT(L),
     &                NINT(MJ-0.5D0))+CGC(IKM1,2)*CGC(IKM2,1)
     &                *GAUNT_CYLM(NINT(L),NINT(MJ-0.5D0),2,-1,NINT(L),
     &                NINT(MJ+0.5D0)))
                  N = NINT(LB(K1))
                  N2 = NINT(LB(K2))
                  SM = -DSQRT(4.0D0*PI/5.0D0)
     &                 *(CGC(IMKM2,1)*CGC(IMKM1,1)*GAUNT_CYLM(N2,
     &                 NINT(MJ+0.5D0),2,0,N,NINT(MJ+0.5D0))-CGC(IMKM2,2)
     &                 *CGC(IMKM1,2)*GAUNT_CYLM(N2,NINT(MJ-0.5D0),2,0,N,
     &                 NINT(MJ-0.5D0)))
                  SM = SM + 3.0D0*DSQRT(2.0D0*PI/15.0D0)
     &                 *(-CGC(IMKM2,1)*CGC(IMKM1,2)
     &                 *GAUNT_CYLM(N2,NINT(MJ+0.5D0),2,1,N,
     &                 NINT(MJ-0.5D0))+CGC(IMKM2,2)*CGC(IMKM1,1)
     &                 *GAUNT_CYLM(N2,NINT(MJ-0.5D0),2,-1,N,
     &                 NINT(MJ+0.5D0)))
                  IF ( K1.EQ.K2 ) THEN
                     QDIA(I) = Q
                     QMDIA(I) = QM
                     SDIA(I) = S
                     SMDIA(I) = SM
                     IF ( CHECK ) WRITE (6,99003) NINT(L),NINT(2*MJ),
     &                    KAP(K1),KAP(K2),I,QDIA(I),QMDIA(I),SDIA(I),
     &                    SMDIA(I),NINT(2*J(K1)),NINT(2*J(K2)),
     &                    NINT(LB(K1)),NINT(LB(K2))
                  ELSE
                     QOFF(I) = Q
                     QMOFF(I) = QM
                     SOFF(I) = S
                     SMOFF(I) = SM
                     IF ( CHECK ) WRITE (6,99003) NINT(L),NINT(2*MJ),
     &                    KAP(K1),KAP(K2),I,QOFF(I),QMOFF(I),SOFF(I),
     &                    SMOFF(I),NINT(2*J(K1)),NINT(2*J(K2)),
     &                    NINT(LB(K1)),NINT(LB(K2))
                  END IF
C
               END DO
            END DO
C
         END DO
      END DO
      WRITE (6,99002) QCGS,VFSHF
C
C                                         ADD THE VARIOUS PRE-FACTORS
      DO I = 1,NKMMAX
         QDIA(I) = QCGS*QDIA(I)
         QMDIA(I) = QCGS*QMDIA(I)
         SDIA(I) = VFSHF*SDIA(I)
         SMDIA(I) = VFSHF*SMDIA(I)
C
         QOFF(I) = QCGS*QOFF(I)
         QMOFF(I) = QCGS*QMOFF(I)
         SOFF(I) = VFSHF*SOFF(I)
         SMOFF(I) = VFSHF*SMOFF(I)
      END DO
C
99001 FORMAT (//,10X,'STOP IN <CALCSMQ>  ',/,10X,'  INDEX I =',I3,
     &        ' > ARRAYSIZE =  ',I3)
99002 FORMAT (/,5X,'QCGS  = ',E12.5,' [1/cm**3]',/,5X,'VFSHF = ',E12.5,
     &        ' [erg/G*cm**3]',/)
99003 FORMAT (I6,I3,'/2',I4,2I5,4F11.5,2(I3,'/2'),2I4)
99004 FORMAT (/,5X,
     &        'QUAD- and SPIN-MAGN.-hyperfine-coefficients for NL = ',
     &        I3,//,5X,'L  MUE KAP1 KAP2 IKM   Q(K1,K2) Q(-K1,-K2)',3X,
     &        'S(K1,K2)  S(-K1,-K2)  J1   J2  LB1 LB2')
C
      END
