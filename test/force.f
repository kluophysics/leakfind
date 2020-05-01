C*==force.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FORCE(CMNTMTT,VNEW,BNEW,F_LMT,FHF_LMT)
C   ********************************************************************
C   *                                                                  *
C   *   calculate the atomic forces                                    *
C   *                                                                  *
C   *   STEP 1: force on nucleus according to Hellmann - Feynman       *
C   *           theorem from non spherical charge density at the       *
C   *           nucleur site                                           *
C   *                                                                  *
C   *   STEP 2: calculate core correction to get force on atom         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_ANGMOM,ONLY:NL,NSPIN
      USE MOD_RMESH,ONLY:NRMAX,R,JRMT,DRDI_W_RADINT
      USE MOD_TYPES,ONLY:Z,IMT,NTMAX,ITBOT,ITTOP,NLMFPMAX,RHO2NS,
     &    RHOCHRC,RHOSPNC
      USE MOD_CONSTANTS,ONLY:CONST_4PI,CONST_4PIOV3,SQRT_4PIOV3
      USE MOD_FILES,ONLY:IFILBUILDBOT,WRBUILDBOT
      IMPLICIT NONE
C*--FORCE23
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='FORCE')
C
C Dummy arguments
C
      REAL*8 BNEW(NRMAX,NLMFPMAX,NTMAX),CMNTMTT(NLMFPMAX,NTMAX),
     &       FHF_LMT(2:4,NTMAX),F_LMT(2:4,NTMAX),
     &       VNEW(NRMAX,NLMFPMAX,NTMAX)
C
C Local variables
C
      REAL*8 DVS(NRMAX),FC_LMT(2:4,NTMAX),FHF_LMT1(2:4),FHF_LMT2(2:4),
     &       F_IR,RHOCSR2(NRMAX),RTOP,SPNWGT,V1(NRMAX),VINT1,VS(NRMAX)
      INTEGER IM,IR,IRMTIN,IS,IT,LM,M
C
C*** End of declarations rewritten by SPAG
C
      IF ( NL.LT.2 ) THEN
         WRITE (6,FMT=99001)
         STOP
      END IF
C
C-----------------------------------------------------------------------
C      STEP 1: force on nucleus according to Hellmann - Feynman
C              theorem from non spherical charge density at the
C              nucleur site
C-----------------------------------------------------------------------
C
      DO IT = ITBOT,ITTOP
C
         IM = IMT(IT)
         IRMTIN = JRMT(IM)
         RTOP = R(IRMTIN,IM)
C
         DO M = -1,1
            LM = 2 + M + 1
C
C---> integrate with integration scheme (Simpson)
C
            VINT1 = 0D0
            DO IR = 1,IRMTIN
               F_IR = RHO2NS(IR,LM,IT,1)/(R(IR,IM)*R(IR,IM))
               VINT1 = VINT1 + F_IR*DRDI_W_RADINT(IR,IM)
            END DO
C
            FHF_LMT1(LM) = 2.0D0*VINT1
C
C---> use coulomb potential to determine extra atomic contribution
C
            FHF_LMT2(LM) = VNEW(IRMTIN,LM,IT)/(CONST_4PIOV3*RTOP)
     &                     - 2.0D0*CMNTMTT(LM,IT)/(RTOP**3)
C
C---> total Hellman-Feynman force
C
            FHF_LMT(LM,IT) = (FHF_LMT1(LM)+FHF_LMT2(LM))*Z(IT)
C
         END DO
      END DO
C
C-----------------------------------------------------------------------
C                  STEP 2: calculate core correction
C-----------------------------------------------------------------------
C
      DO IT = ITBOT,ITTOP
C
         IM = IMT(IT)
         IRMTIN = JRMT(IM)
         RTOP = R(IRMTIN,IM)
C
         DO M = -1,1
            LM = 2 + M + 1
C
C---> initialize v1
C
            V1(1:IRMTIN) = 0.0D0
C
            DO IS = 1,NSPIN
C
               IF ( IREL.GE.2 ) THEN
                  SPNWGT = NINT((IS-1.5D0)*2D0)
               ELSE
                  SPNWGT = 0D0
               END IF
C
               VS(:) = 0D0
               IF ( IREL.LE.1 ) THEN
                  DO IR = 1,IRMTIN
                     VS(IR) = VNEW(IR,LM,IT)
                     RHOCSR2(IR) = RHOCHRC(IR,IT)/CONST_4PI
                  END DO
               ELSE
                  DO IR = 1,IRMTIN
                     VS(IR) = VNEW(IR,LM,IT) + SPNWGT*BNEW(IR,LM,IT)
                     RHOCSR2(IR) = (RHOCHRC(IR,IT)+SPNWGT*RHOSPNC(IR,IT)
     &                             )/CONST_4PI
                  END DO
               END IF
C
C---> determine the derivative of the potential
C
               CALL RDIFFER(IM,VS,DVS)
C
               DO IR = 1,IRMTIN
                  V1(IR) = V1(IR) + RHOCSR2(IR)
     &                     *(2.0D0*VS(IR)/R(IR,IM)+DVS(IR))
               END DO
C
            END DO
C
C---> integrate with integration scheme (Simpson)
C
            VINT1 = 0D0
            DO IR = 1,IRMTIN
               VINT1 = VINT1 + V1(IR)*DRDI_W_RADINT(IR,IM)
            END DO
C
            FHF_LMT(LM,IT) = SQRT_4PIOV3*FHF_LMT(LM,IT)
            FC_LMT(LM,IT) = -SQRT_4PIOV3*VINT1
            F_LMT(LM,IT) = FHF_LMT(LM,IT) + FC_LMT(LM,IT)
C
         END DO
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
         IF ( WRBUILDBOT ) WRITE (IFILBUILDBOT,99002)
     &                            ROUTINE(1:LEN_TRIM(ROUTINE)),IT,
     &                            FHF_LMT(2:4,IT),FC_LMT(2:4,IT)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
      END DO
C
99001 FORMAT (10x,'error stop in subroutine <FORCE> :',
     &        ' the charge density has to contain non spherical',
     &        ' contributions up to l=1 at least !')
99002 FORMAT ('# BUILDBOT: ',A,': forces (LM=2..4)  FHF and FC ',
     &        'for IT =',I5,/,(1PE22.14))
      END
