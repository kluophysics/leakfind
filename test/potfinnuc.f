C*==potfinnuc.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE POTFINNUC(VT)
C   ********************************************************************
C   *                                                                  *
C   *   change to finite nucleus potential                             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:R,JRWS,NRMAX,FINITE_NUCLEUS
      USE MOD_TYPES,ONLY:TXT_T,Z,NTMAX,IMT,NT
      IMPLICIT NONE
C*--POTFINNUC12
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 VT(NRMAX,NTMAX)
C
C Local variables
C
      INTEGER IM,IT,J,JRNUC
      REAL*8 POT_NUC,RNUCTAB
      REAL*8 RJ,RNUC,ZZ
C
C*** End of declarations rewritten by SPAG
C
      FINITE_NUCLEUS = .TRUE.
C
      WRITE (6,99001)
C
      DO IT = 1,NT
C
         IM = IMT(IT)
C
         ZZ = -VT(1,IT)*R(1,IM)/2.0D0
C
         IF ( ABS(Z(IT)-ZZ).GT.1 ) THEN
C
            WRITE (6,99002) IT
C
         ELSE
C
            RNUC = RNUCTAB(Z(IT))
            JRNUC = 1
            ZZ = DBLE(Z(IT))
C
            DO J = 1,JRWS(IM)
               RJ = R(J,IM)
               IF ( RNUC.LT.RJ ) EXIT
               VT(J,IT) = VT(J,IT) + 2D0*ZZ/RJ + POT_NUC(RJ,ZZ)
               JRNUC = J
            END DO
C
            WRITE (6,99003) IT,TXT_T(IT),Z(IT),JRNUC,RNUC,R(1,IM)
C
         END IF
C
      END DO
C
99001 FORMAT (/,1X,79('*'),/,16X,
     &        'set potential for finite nucleus in <POTFINNUC>',/,1X,
     &        79('*'),//,15X,'type           Z    JRNUC       RNUC',10X,
     &        'R(1)')
99002 FORMAT (5X,'####### potential for type ',I3,
     &        ' is already according to finite nucleus')
99003 FORMAT (10X,I5,3X,A,I5,I8,2F14.7)
      END
