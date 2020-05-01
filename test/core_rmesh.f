C*==core_rmesh.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CORE_RMESH
C   ********************************************************************
C   *                                                                  *
C   *  set up expanded radial mesh for core states                     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:NRMAX,JRWS,JRCRI,DRDI,R,NRCMAX,EXPDX_ASA,
     &    FULLPOT,NPAN,JRCUT,RMESHTYPE,R_COR,DRDI_COR,DRDIOVR_COR,NM
      IMPLICIT NONE
C*--CORE_RMESH12
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CORE_RMESH')
C
C Local variables
C
      REAL*8 DR,RAT
      INTEGER IM,IR,IR1,IR2,IRTOP,NRC
C
C*** End of declarations rewritten by SPAG
C
      DO IM = 1,NM
C
C=======================================================================
         IF ( RMESHTYPE.EQ.'EXPONENTIAL ' ) THEN
C
C--------------------------------------------------------------- FULLPOT
            IF ( FULLPOT ) THEN
C
               IRTOP = JRCRI(IM)
C
               NRC = 2*IRTOP
               IF ( NRC.GT.NRCMAX )
     &               CALL STOP_MESSAGE(ROUTINE,'NRC > NRCMAX')
C
               DO IR = 1,IRTOP
                  R_COR(IR,IM) = R(IR,IM)
                  DRDI_COR(IR,IM) = DRDI(IR,IM)
                  DRDIOVR_COR(IR,IM) = DRDI_COR(IR,IM)/R_COR(IR,IM)
               END DO
C
C------------------------- continue last panel with equidistant steps DR
               IR2 = IRTOP
               IR1 = JRCUT(NPAN(IM)-1,IM) + 1
               DR = (R(IR2,IM)-R(IR1,IM))/DBLE(IR2-IR1)
               IF ( ABS(1D0-(R(IR2,IM)-R(IR2-1,IM))/DR).GT.1D-6 )
     &              CALL STOP_MESSAGE(ROUTINE,'DR not consistent')
C
               DO IR = (IRTOP+1),NRC
                  R_COR(IR,IM) = R_COR(IR-1,IM) + DR
                  DRDI_COR(IR,IM) = DRDI_COR(IRTOP,IM)
                  DRDIOVR_COR(IR,IM) = DRDI_COR(IR,IM)/R_COR(IR,IM)
               END DO
C
            ELSE
C------------------------------------------------------------------- ASA
C
               IRTOP = JRWS(IM)
C
               NRC = 2*IRTOP
               IF ( NRC.GT.NRCMAX )
     &               CALL STOP_MESSAGE(ROUTINE,'NRC > NRCMAX')
C
               DO IR = 1,NRMAX
                  R_COR(IR,IM) = R(IR,IM)
                  DRDI_COR(IR,IM) = DRDI(IR,IM)
                  DRDIOVR_COR(IR,IM) = DRDI_COR(IR,IM)/R_COR(IR,IM)
               END DO
C
               RAT = EXPDX_ASA(IM)
C
               DO IR = (NRMAX+1),NRC
                  R_COR(IR,IM) = RAT*R_COR(IR-1,IM)
                  DRDI_COR(IR,IM) = (RAT-1.0D0)*R_COR(IR-1,IM)
                  DRDIOVR_COR(IR,IM) = DRDI_COR(IR,IM)/R_COR(IR,IM)
               END DO
            END IF
C
C=======================================================================
         ELSE IF ( RMESHTYPE.EQ.'JUELICH     ' ) THEN
C
            CALL STOP_MESSAGE(ROUTINE,'TODO: RMESHTYPE.EQ.JUELICH')
C
C=======================================================================
         ELSE
            CALL STOP_MESSAGE(ROUTINE,'check RMESHTYPE ')
         END IF
C=======================================================================
C
      END DO
C
      END
