C*==cdiffer.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE CDIFFER(IM,F,FP)
C   *******************************************************************
C   *                                                                 *
C   *   differentiate COMPLEX function F(r)                           *
C   *                                                                 *
C   *        FP = dF/dr = dI/dr * dF/dI                               *
C   *                                                                 *
C   *******************************************************************
C
      USE MOD_RMESH,ONLY:NPAN,NRMAX,FULLPOT,JRCUT,JRWS,DRDI
      IMPLICIT NONE
C*--CDIFFER13
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CDIFFER')
C
C Dummy arguments
C
      INTEGER IM
      COMPLEX*16 F(NRMAX),FP(NRMAX)
C
C Local variables
C
      INTEGER IPAN,IR,N1,N2
C
C*** End of declarations rewritten by SPAG
C
C---------------------------------------- ensure proper settings for ASA
      IF ( .NOT.FULLPOT ) THEN
         NPAN(IM) = 1
         JRCUT(0,IM) = 0
         JRCUT(1,IM) = JRWS(IM)
      ELSE
         IF ( NPAN(IM).EQ.1 )
     &         CALL STOP_MESSAGE(ROUTINE,'NPAN(IM) = 1  for FULLPOT')
      END IF
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C     loop over panels
C-----------------------------------------------------------------------
C
      DO IPAN = 1,NPAN(IM)
         N1 = JRCUT(IPAN-1,IM) + 1
         N2 = JRCUT(IPAN,IM)
C
C--------------------- forward difference at the beginning of the table
C
         FP(N1) = ((2.D0*F(N1+3)+18.D0*F(N1+1))-(9.D0*F(N1+2)+11.D0*F(N1
     &            )))/6.0D0
         FP(N1+1) = ((2.D0*F(N1+4)+18.D0*F(N1+2))-(9.D0*F(N1+3)+11.D0*F(
     &              N1+1)))/6.0D0
C
C---------------------- central difference at the interior of the table
C
         DO IR = N1 + 2,N2 - 2
            FP(IR) = ((F(IR-2)+8.D0*F(IR+1))-(8.D0*F(IR-1)+F(IR+2)))
     &               /12.0D0
         END DO
C
C-------------------------- backward difference at the end of the table
C
         FP(N2) = ((11.D0*F(N2)+9.D0*F(N2-2))-(18.D0*F(N2-1)+2.D0*F(N2-3
     &            )))/6.0D0
         FP(N2-1) = ((11.D0*F(N2-1)+9.D0*F(N2-3))-(18.D0*F(N2-2)+2.D0*F(
     &              N2-4)))/6.0D0
C
      END DO
C
C-----------------------------------------------------------------------
C                 FP = dF/dr = dI/dr * dF/dI
C-----------------------------------------------------------------------
C
      DO IR = 1,JRCUT(NPAN(IM),IM)
         FP(IR) = FP(IR)/DRDI(IR,IM)
      END DO
C
      END
