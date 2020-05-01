C*==cradint.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CRADINT(IM,F,FINT)
C   ********************************************************************
C   *                                                                  *
C   * integrates a COMPLEX FUNCTION F(i) tabulated on the mesh  IM     *
C   * OUTWARDS  using the 3 point Simpson method                       *
C   *                                                                  *
C   *                      r_top                                       *
C   *               FINT = S      F(r') dr'                            *
C   *                      0                                           *
C   *                                                                  *
C   * - works for ASA (r_top = r_WS) and FULLPOT (r_top = r_crit)      *
C   * - FULLPOT:  at each kink the integration is restarted            *
C   * - F(i) has to include the weight  r^2 * dr/di                    *
C   * - the number of mesh points in each panel can be even or odd     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:NRMAX,FULLPOT,JRCRI,JRWS,W_RADINT
      IMPLICIT NONE
C*--CRADINT21
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 FINT
      INTEGER IM
      COMPLEX*16 F(NRMAX)
C
C Local variables
C
      INTEGER IR,IRTOP
C
C*** End of declarations rewritten by SPAG
C
C-----------------------------------------------------------------------
      IF ( FULLPOT ) THEN
         IRTOP = JRCRI(IM)
      ELSE
         IRTOP = JRWS(IM)
      END IF
C-----------------------------------------------------------------------
C
      FINT = 0.0D0
C
      DO IR = 1,IRTOP
C
         FINT = FINT + F(IR)*W_RADINT(IR,IM)
C
      END DO
C
      END
C*==cradint_r.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CRADINT_R(IM,F,FINT)
C   ********************************************************************
C   *                                                                  *
C   * integrates a COMPLEX FUNCTION F(i) tabulated on the mesh  IM     *
C   * OUTWARDS  using the 3 point Simpson method                       *
C   *                                                                  *
C   *                         r                                        *
C   *               FINT(r) = S   F(r') dr'                            *
C   *                         0                                        *
C   *                                                                  *
C   * - works for ASA (r_top = r_WS) and FULLPOT (r_top = r_crit)      *
C   * - FULLPOT:  at each kink the integration is restarted            *
C   * - F(i) has to include the weight  r^2 * dr/di                    *
C   * - the number of mesh points in each panel can be even or odd     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:NPAN,NRMAX,JRCUT
      IMPLICIT NONE
C*--CRADINT_R85
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CRADINT_R')
      REAL*8 R12,R3
      PARAMETER (R12=1D0/12D0,R3=1D0/3D0)
C
C Dummy arguments
C
      INTEGER IM
      COMPLEX*16 F(NRMAX),FINT(NRMAX)
C
C Local variables
C
      COMPLEX*16 FINTREF,FIP1_2,FIP1_4,FI_2,FI_5
      INTEGER IPAN,IR,IR_END,IR_START
C
C*** End of declarations rewritten by SPAG
C
      FINT(1) = 0.0D0
C
C-----------------------------------------------------------------------
C     loop over panels
C-----------------------------------------------------------------------
C
      DO IPAN = 1,NPAN(IM)
C
         IR_END = JRCUT(IPAN,IM)
         IR_START = JRCUT(IPAN-1,IM) + 1
         IF ( IR_START+2.GT.IR_END ) CALL STOP_MESSAGE(ROUTINE,
     &        'at least 3 mesh points in panel required')
C
         IF ( IPAN.EQ.1 ) THEN
            FINT(IR_START) = 0.0D0
         ELSE
            FINT(IR_START) = FINT(IR_START-1)
         END IF
         FINTREF = FINT(IR_START)
C
         IF ( MOD((IR_END-IR_START),2).NE.0 ) THEN
            IR = IR_START
            FIP1_2 = F(IR+1) + F(IR+1)
            FIP1_4 = FIP1_2 + FIP1_2
            FI_2 = F(IR) + F(IR)
            FI_5 = FI_2 + FI_2 + F(IR)
            FINT(IR+1) = FINTREF + (FI_5+FIP1_4+FIP1_4-F(IR+2))*R12
            FINTREF = FINT(IR+1)
            IR_START = IR_START + 1
         END IF
C
         DO IR = IR_START,IR_END - 2,2
            FIP1_2 = F(IR+1) + F(IR+1)
            FIP1_4 = FIP1_2 + FIP1_2
            FI_2 = F(IR) + F(IR)
            FI_5 = FI_2 + FI_2 + F(IR)
            FINT(IR+1) = FINTREF + (FI_5+FIP1_4+FIP1_4-F(IR+2))*R12
            FINT(IR+2) = FINTREF + (F(IR)+FIP1_4+F(IR+2))*R3
            FINTREF = FINT(IR+2)
         END DO
C
      END DO
C
      END
