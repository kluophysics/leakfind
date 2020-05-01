C*==posanilambda.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE POSANILAMBDA(RHOCHR_EL,RHO2NS_EL)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the positron life time parameter   LAMBDA             *
C   *                                                                  *
C   *  ASA:  RHOCHR = 4*pi* RHO                                        *
C   *  FP:   RHO2NS = 4*pi* RHO[l,m] * r**2                            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:NRMAX,FLMSF,R,LMISF,KLMSF,NSF,JRCUT,JRWS,NPAN,
     &    FULLPOT,R2DRDI
      USE MOD_ANGMOM,ONLY:NL
      USE MOD_TYPES,ONLY:NLMFPMAX,NTMAX,NLMFPT,KLMFP,IMT,RHO2NS,RHOCHR,
     &    ITBOT,ITTOP,NAT,CONC
      USE MOD_CONSTANTS,ONLY:PI,CONST_4PI,C_CGS,A0_CGS,RE_CGS,SQRT_4PI
      IMPLICIT NONE
C*--POSANILAMBDA19
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NGAU,NMAX_GAUGRID
      PARAMETER (NGAU=20,NMAX_GAUGRID=NGAU*NGAU)
C
C Dummy arguments
C
      REAL*8 RHO2NS_EL(NRMAX,NLMFPMAX,NTMAX,3),RHOCHR_EL(NRMAX,NTMAX)
C
C Local variables
C
      REAL*8 ANGINT,AUX,CTET,GFACTOR(:),NORM,PHI,PREFAC,QPOS,QPOS_T(:),
     &       RATE,RATET(:),RHAT(3),RHO_EL(:),RHO_POS(:),RINT(:),SFN(:),
     &       STET,WGAU(:),W_GAUGRID(:),X,XGAU(:),Y_GAUGRID(:,:)
      INTEGER IA_ERR,IM,IPHI,IR,IRCRIT,IRMTIN,IRSF,IRTOP,ISF,IT,ITET,
     &        I_GAUGRID,JSF,LM,LM_SF,NLMSF,NLSF,N_GAUGRID
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE RHO_POS,RHO_EL,GFACTOR,RINT,QPOS_T
      ALLOCATABLE W_GAUGRID,Y_GAUGRID,RATET,SFN
C
      ALLOCATABLE XGAU,WGAU
C
      ALLOCATE (RATET(NTMAX),RINT(NRMAX),QPOS_T(NTMAX))
C
      WRITE (6,99002)
C
C=======================================================================
C                                  ASA
C=======================================================================
C
      IF ( .NOT.FULLPOT ) THEN
C
         ALLOCATE (RHO_POS(NRMAX),RHO_EL(NRMAX),GFACTOR(NRMAX))
C
C-------------------------------- normalize positron charge density to 1
         QPOS = 0D0
         DO IT = ITBOT,ITTOP
C
            IM = IMT(IT)
            IRTOP = JRWS(IM)
            RINT(1:IRTOP) = RHOCHR(1:IRTOP,IT)*R2DRDI(1:IRTOP,IM)
C
            CALL RRADINT(IM,RINT,QPOS_T(IT))
C
            QPOS = QPOS + NAT(IT)*CONC(IT)*QPOS_T(IT)
C
         END DO
C
         WRITE (6,99005) QPOS
         IF ( QPOS.GT.1D-8 ) THEN
            X = 1D0/QPOS
            RHOCHR(1:NRMAX,ITBOT:ITTOP) = X*RHOCHR(1:NRMAX,ITBOT:ITTOP)
         ELSE
            STOP '<POSANILAMBDA>:   positron charge is 0 '
         END IF
C
         RATE = 0D0
         DO IT = ITBOT,ITTOP
C
            IM = IMT(IT)
            IRTOP = JRWS(IM)
C
            DO IR = 1,IRTOP
               RHO_POS(IR) = RHOCHR(IR,IT)/CONST_4PI
               RHO_EL(IR) = RHOCHR_EL(IR,IT)/CONST_4PI
            END DO
C
            CALL POSANI_GFACTOR(RHO_EL,GFACTOR,IRTOP)
C
            DO IR = 1,IRTOP
               RINT(IR) = RHO_POS(IR)*RHO_EL(IR)*GFACTOR(IR)
     &                    *R2DRDI(IR,IM)
            END DO
C
            CALL RRADINT(IM,RINT,RATET(IT))
C
            RATET(IT) = RATET(IT)*CONST_4PI
C
            RATE = RATE + NAT(IT)*CONC(IT)*RATET(IT)
C
         END DO
C
      ELSE
C=======================================================================
C                        FULL POTENTIAL
C=======================================================================
C
C ----------------------------------------------- construct angular mesh
C
         NLSF = 4*(NL-1) + 1
         NLMSF = NLSF**2
C
         N_GAUGRID = NMAX_GAUGRID
C
         ALLOCATE (W_GAUGRID(NMAX_GAUGRID))
         ALLOCATE (Y_GAUGRID(NLMSF,NMAX_GAUGRID),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP 'ALLOC: POSANILAMBDA -> Y_GAUGRID'
C
         ALLOCATE (XGAU(NGAU),WGAU(NGAU))
C
         CALL GAULEG(-1D0,+1D0,XGAU,WGAU,NGAU)
C
         I_GAUGRID = 0
         DO ITET = 1,NGAU
            CTET = XGAU(ITET)
            STET = SQRT(1D0-CTET*CTET)
            DO IPHI = 1,NGAU
               PHI = 2*PI*IPHI/NGAU
               RHAT(1) = STET*COS(PHI)
               RHAT(2) = STET*SIN(PHI)
               RHAT(3) = CTET
C
               I_GAUGRID = I_GAUGRID + 1
C
               CALL CALC_RHPLM(RHAT(1),RHAT(2),RHAT(3),
     &                         Y_GAUGRID(1,I_GAUGRID),NLSF-1,NLMSF)
C
               W_GAUGRID(I_GAUGRID) = WGAU(ITET)*2*PI/NGAU
C
            END DO
         END DO
C
         ALLOCATE (RHO_POS(NMAX_GAUGRID),RHO_EL(NMAX_GAUGRID))
         ALLOCATE (GFACTOR(NMAX_GAUGRID),SFN(NMAX_GAUGRID))
C
C-------------------------------- normalize positron charge density to 1
         QPOS = 0D0
C
         DO IT = ITBOT,ITTOP
            IM = IMT(IT)
C
            IRCRIT = JRCUT(NPAN(IM),IM)
            IRMTIN = JRCUT(1,IM)
            RINT(1:IRCRIT) = 0D0
C
            DO IR = 1,IRMTIN
               NORM = R(IR,IM)**2
               RINT(IR) = RHO2NS(IR,1,IT,1)*R2DRDI(IR,IM)/NORM*SQRT_4PI
            END DO
C
            DO LM = 1,NLMFPT(IT)
               IF ( KLMFP(LM,IT).NE.0 ) THEN
C
                  DO ISF = 1,NSF(IM)
                     LM_SF = LMISF(ISF,IM)
                     IF ( LM_SF.EQ.LM ) THEN
                        JSF = ISF
                        GOTO 5
                     END IF
                  END DO
                  CYCLE
 5                CONTINUE
                  ISF = JSF
C
                  DO IR = IRMTIN + 1,IRCRIT
                     NORM = R(IR,IM)**2
                     IRSF = IR - IRMTIN
                     RINT(IR) = RINT(IR) + RHO2NS(IR,LM,IT,1)
     &                          *FLMSF(IRSF,ISF,IM)*R2DRDI(IR,IM)/NORM
                  END DO
               END IF
            END DO
C
C-------------------------------------------- perform radial integration
C
            CALL RRADINT(IM,RINT,AUX)
C
            QPOS = QPOS + NAT(IT)*CONC(IT)*AUX
C
         END DO
C
         WRITE (6,99005) QPOS
         IF ( QPOS.GT.1D-8 ) THEN
            X = 1D0/QPOS
            DO IT = ITBOT,ITTOP
C
               IRCRIT = JRCUT(NPAN(IM),IM)
               DO LM = 1,NLMFPT(IT)
                  IF ( KLMFP(LM,IT).NE.0 ) RHO2NS(1:IRCRIT,LM,IT,1)
     &                 = X*RHO2NS(1:IRCRIT,LM,IT,1)
               END DO
            END DO
         ELSE
            STOP '<POSANILAMBDA>:   positron charge is < 1D-8 '
         END IF
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
         RATE = 0D0
         DO IT = ITBOT,ITTOP
C
            IM = IMT(IT)
            IRCRIT = JRCUT(NPAN(IM),IM)
            IRMTIN = JRCUT(1,IM)
C
C------------------------------------------------- loop over radial mesh
            DO IR = 1,IRCRIT
C
C---------- generate the densities and shape function on an angular mesh
C
               RHO_POS(1:N_GAUGRID) = 0.0D0
               RHO_EL(1:N_GAUGRID) = 0.0D0
C
               DO LM = 1,NLMFPT(IT)
                  IF ( KLMFP(LM,IT).NE.0 ) THEN
C
                     DO I_GAUGRID = 1,N_GAUGRID
                        RHO_POS(I_GAUGRID) = RHO_POS(I_GAUGRID)
     &                     + RHO2NS(IR,LM,IT,1)*Y_GAUGRID(LM,I_GAUGRID)
C
                        RHO_EL(I_GAUGRID) = RHO_EL(I_GAUGRID)
     &                     + RHO2NS_EL(IR,LM,IT,1)
     &                     *Y_GAUGRID(LM,I_GAUGRID)
                     END DO
C
                  END IF
               END DO
C
               NORM = R(IR,IM)**2
               DO I_GAUGRID = 1,N_GAUGRID
                  RHO_POS(I_GAUGRID) = RHO_POS(I_GAUGRID)/NORM
                  RHO_EL(I_GAUGRID) = RHO_EL(I_GAUGRID)/NORM
               END DO
C
               CALL POSANI_GFACTOR(RHO_EL,GFACTOR,N_GAUGRID)
C
               IF ( IR.LE.IRMTIN ) THEN
C
                  SFN(1:N_GAUGRID) = 1D0
C
               ELSE
C
                  SFN(1:N_GAUGRID) = 0D0
C
                  IRSF = IR - IRMTIN
C
                  DO ISF = 1,NSF(IM)
                     LM = LMISF(ISF,IM)
                     IF ( KLMSF(LM,IM).EQ.1 ) THEN
C
                        DO I_GAUGRID = 1,N_GAUGRID
                           SFN(I_GAUGRID) = SFN(I_GAUGRID)
     &                        + FLMSF(IRSF,ISF,IM)
     &                        *Y_GAUGRID(LM,I_GAUGRID)
                        END DO
C
                     END IF
                  END DO
C
               END IF
C
C------------------ perform angular integration via Gaussian integration
C
               ANGINT = 0D0
               DO I_GAUGRID = 1,N_GAUGRID
                  ANGINT = ANGINT + RHO_EL(I_GAUGRID)*RHO_POS(I_GAUGRID)
     &                     *GFACTOR(I_GAUGRID)*SFN(I_GAUGRID)
     &                     *W_GAUGRID(I_GAUGRID)
               END DO
C
               RINT(IR) = ANGINT*R2DRDI(IR,IM)
C
            END DO
C------------------------------------------- loop over radial mesh - END
C
C-------------------------------------------- perform radial integration
C
            CALL RRADINT(IM,RINT,RATET(IT))
C
            RATE = RATE + NAT(IT)*CONC(IT)*RATET(IT)
C
         END DO
C
      END IF
C=======================================================================
C
      PREFAC = PI*C_CGS*RE_CGS**2/(A0_CGS**3)
C
      RATET(ITBOT:ITTOP) = PREFAC*RATET(ITBOT:ITTOP)
      RATE = PREFAC*RATE
C
      IF ( RATE.GE.1D-8 ) THEN
         WRITE (6,99001)
         DO IT = ITBOT,ITTOP
            WRITE (6,99003) IT,RATET(IT),1D0/RATET(IT)
         END DO
         WRITE (6,99004) RATE,1D0/RATE,1D+12/RATE
      ELSE
         WRITE (6,*) 'ERROR: Positron annihilation rate ',RATE
      END IF
C
99001 FORMAT (8X,'IT    rate  [1/s]     life time [s]    [ps]'/)
99002 FORMAT (/,1X,79('*'),/,33X,'<POSANILAMBDA>',/,79('*'),//,10X,
     &        'calculating the positron life time',/)
99003 FORMAT (I10,E15.5,1X,E15.5)
99004 FORMAT (10X,2(E15.5,1X),4X,F8.3)
99005 FORMAT (/,10X,'total positronic charge from charge density ',
     &        F12.6,/,10X,'charge density will be normalized to 1 ',/)
      END
