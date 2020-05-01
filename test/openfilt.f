C*==openfilt.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE OPENFILT(DATSET,LDATSET,TXT_T,LTXT_T,IT,EXT,LEXT,
     &                    FILNAM,LFN,KOPEN,IFIL,TEXT,LTEXT,NTMAX)
C   ********************************************************************
C   *                                                                  *
C   *   create a filename for atom type   IT                           *
C   *                                                                  *
C   *   {DATSET_}TXT_T(IT)EXT    {...}  optional                       *
C   *                                                                  *
C   *   KOPEN =   X: create filename                                   *
C   *           > 0:     + OPEN file IFIL                              *
C   *             1:     + OPEN file IFIL + print info TEXT            *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--OPENFILT17
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*80 DATSET,FILNAM
      CHARACTER*(*) EXT,TEXT
      INTEGER IFIL,IT,KOPEN,LDATSET,LEXT,LFN,LTEXT,NTMAX
      INTEGER LTXT_T(NTMAX)
      CHARACTER*(*) TXT_T(NTMAX)
C
C Local variables
C
      CHARACTER*28 FMTSTR
C
C*** End of declarations rewritten by SPAG
C
      IF ( LDATSET.NE.0 ) THEN
         FILNAM = DATSET(1:LDATSET)//'_'//TXT_T(IT)(1:LTXT_T(IT))
         LFN = LDATSET + 1 + LTXT_T(IT)
      ELSE
         FILNAM = TXT_T(IT)(1:LTXT_T(IT))
         LFN = LTXT_T(IT)
      END IF
C
      FILNAM = FILNAM(1:LFN)//EXT(1:LEXT)
      LFN = LFN + LEXT
C
      IF ( KOPEN.GT.0 ) THEN
C
         IF ( IFIL.GT.99 ) THEN
            FMTSTR = '(/,10X,A,'': ('',I3,'') '',A) '
         ELSE
            FMTSTR = '(/,10X,A,'':  ('',I2,'') '',A)'
         END IF
C
         OPEN (UNIT=IFIL,FILE=FILNAM)
C
         IF ( KOPEN.GT.1 ) WRITE (6,FMT=FMTSTR) TEXT(1:LTEXT),IFIL,
     &                            FILNAM(1:LFN)
C
      END IF
C
      END
