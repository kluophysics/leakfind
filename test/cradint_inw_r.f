C*==cradint_inw_r.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CRADINT_INW_R(IM,F,FINT)
C   ********************************************************************
C   *                                                                  *
C   * integrates a COMPLEX FUNCTION F(i) tabulated on the mesh  IM     *
C   * INWARDS   using the 3 point Simpson method                       *
C   *                                                                  *
C   *                         r_top                                    *
C   *               FINT(r) = S      F(r') dr'                         *
C   *                         r                                        *
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
C*--CRADINT_INW_R21
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CRADINT_INW_R')
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
      COMPLEX*16 FINTREF,FIP1_2,FIP1_4,FIP2_2,FIP2_5,FI_2,FI_5
      INTEGER IPAN,IR,IR_END,IR_START
      LOGICAL N_MESH_POINTS_IS_EVEN
C
C*** End of declarations rewritten by SPAG
C
C-----------------------------------------------------------------------
C     loop over panels
C-----------------------------------------------------------------------
C
      DO IPAN = NPAN(IM),1, - 1
C
         IR_START = JRCUT(IPAN,IM)
         IR_END = JRCUT(IPAN-1,IM) + 1
         IF ( IR_END+2.GT.IR_START ) CALL STOP_MESSAGE(ROUTINE,
     &        'at least 3 mesh points in panel required')
C
         IF ( IPAN.EQ.NPAN(IM) ) THEN
            FINT(IR_START) = 0.0D0
         ELSE
            FINT(IR_START) = FINT(IR_START+1)
         END IF
         FINTREF = FINT(IR_START)
C
         IF ( MOD((IR_END-IR_START),2).NE.0 ) THEN
            N_MESH_POINTS_IS_EVEN = .TRUE.
            IR_END = IR_END + 1
         ELSE
            N_MESH_POINTS_IS_EVEN = .FALSE.
         END IF
C
         DO IR = IR_START - 2,IR_END, - 2
            FIP1_2 = F(IR+1) + F(IR+1)
            FIP1_4 = FIP1_2 + FIP1_2
            FIP2_2 = F(IR+2) + F(IR+2)
            FIP2_5 = FIP2_2 + FIP2_2 + F(IR+2)
            FINT(IR+1) = FINTREF + (-F(IR)+FIP1_4+FIP1_4+FIP2_5)*R12
            FINT(IR) = FINTREF + (F(IR)+FIP1_4+F(IR+2))*R3
            FINTREF = FINT(IR)
         END DO
C
         IF ( N_MESH_POINTS_IS_EVEN ) THEN
            IR = IR_END - 1
            FI_2 = F(IR) + F(IR)
            FI_5 = FI_2 + FI_2 + F(IR)
            FIP1_2 = F(IR+1) + F(IR+1)
            FIP1_4 = FIP1_2 + FIP1_2
            FINT(IR) = FINTREF + (FI_5+FIP1_4+FIP1_4-F(IR+2))*R12
            FINTREF = FINT(IR)
         END IF
C
      END DO
C
      END
