C*==input_find_section.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE INPUT_FIND_SECTION(KW,STOP_KEY)
C   ********************************************************************
C   *                                                                  *
C   *  find keyword KW(1:LKW) starting at column 1 in file  IFILINP    *
C   *  combine this and all follwing lines to a string CMD00(1:LCMD)   *
C   *  until a non-blank character in the 1 column or EOF is found     *
C   *  CMDUC contains an  upper case - copy of CMD00                   *
C   *                                                                  *
C   *  the KW is expected to be upper case                             *
C   *                                                                  *
C   *  only  LLMAX  characters per input line are accepted             *
C   *                                                                  *
C   *  lines starting with # are comment - lines                       *
C   *  input after # is interpreted as  comment                        *
C   *  the character '=' is replaced by a blank                        *
C   *                                                                  *
C   *  a blank is added to  CMDUC/CMD00 at the end     >>>  <FINDKW>   *
C   *                                                                  *
C   *  STOP_KEY = 0     CONTINUE if SECTION is not found               *
C   *             1     STOP     if SECTION is not found               *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_FILES,ONLY:LCMDMAX,LCMD,CMD00,CMDUC,IFILINP,FOUND_SECTION
      IMPLICIT NONE
C*--INPUT_FIND_SECTION28
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='INPUT_FIND_SECTION')
      INTEGER LLMAX
      PARAMETER (LLMAX=9999) !REC increased
C
C Dummy arguments
C
      CHARACTER*(*) KW
      INTEGER STOP_KEY
C
C Local variables
C
      INTEGER I,LKW,LL
      CHARACTER*(LLMAX) LINE00,LINEUC
C
C*** End of declarations rewritten by SPAG
C
      FOUND_SECTION = .FALSE.
C
      REWIND IFILINP
      LCMD = 0
      LKW = LEN(KW)
      CMD00(1:LEN(CMD00)) = ' '
C
C-------------------------------------------------------- serach KEYWORD
 100  CONTINUE
      READ (IFILINP,'(A)',END=300) LINE00
      CALL STRING_REPLACE_TAB(LINE00)
C
      LINEUC = LINE00
C
      CALL STRING_CONVERT_TO_UC(LINEUC(1:MIN(LKW,LLMAX)))
      IF ( LINEUC(1:LKW).NE.KW(1:LKW) ) GOTO 100
C
      FOUND_SECTION = .TRUE.
C
C------------------------------------------------------------- set CMD00
      LL = LEN_TRIM(LINEUC)
      LCMD = LL - LKW + 1
      CMD00 = LINE00(MIN(LKW+1,LL):LL)//' '
C
C------------------------------------------------ read continuation line
C
 200  CONTINUE
      READ (IFILINP,'(A)',END=300) LINE00
C
      CALL STRING_REPLACE_TAB(LINE00)
C
      IF ( LINE00(1:1).EQ.'#' ) GOTO 200
C
      IF ( LINE00(1:1).EQ.' ' ) THEN
         LL = LEN_TRIM(LINE00)
         IF ( LL.NE.0 ) THEN
C
            CALL STRING_TRIM_LEFT(LINE00)
            LL = LEN_TRIM(LINE00)
            I = INDEX(LINE00(1:LL),'#')
C
            IF ( I.NE.1 ) THEN
               IF ( I.GT.1 ) LL = I - 1
C
               IF ( LCMD+LL.GT.LCMDMAX ) THEN
                  WRITE (6,'(A,A)') 'CMD :',CMD00(1:LCMD)
                  WRITE (6,'(A,A)') '+  LINE :',LINE00(1:LL)
                  WRITE (6,'(A,I4)') 'TOO LONG !  LCMDMAX=',LCMDMAX
               END IF
C
               LCMD = LEN_TRIM(CMD00)
               CMD00 = CMD00(1:LCMD)//' '//LINE00(1:LL)
               LCMD = LCMD + 1 + LL
C
            END IF
C
            GOTO 200
C
         END IF
C
      END IF
C-----------------------------------------------------------------------
C
 300  CONTINUE
      IF ( FOUND_SECTION ) THEN
C
         CMDUC = CMD00
         CALL STRING_CONVERT_TO_UC(CMDUC)
C
C-------------------------------------------------------- modify SECTION
         LCMD = LCMD + 1
         CMD00(LCMD:LCMD) = ' '
         CMDUC(LCMD:LCMD) = ' '
         CALL EXCHANGE_CHARACTERS('=',' ',CMD00)
         CALL EXCHANGE_CHARACTERS('=',' ',CMDUC)
         CALL EXCHANGE_CHARACTERS('(','{',CMDUC)
         CALL EXCHANGE_CHARACTERS(')','}',CMDUC)
         CALL EXCHANGE_CHARACTERS('[','{',CMDUC)
         CALL EXCHANGE_CHARACTERS(']','}',CMDUC)
         CALL EXCHANGE_CHARACTERS(';',',',CMDUC)
C
C----------------------------------------------- STOP if STOP_KEY is set
C
      ELSE IF ( STOP_KEY.EQ.1 ) THEN
C
         WRITE (6,99001) KW(1:LKW)
         CALL STOP_MESSAGE(ROUTINE,'requested SECTION not found')
C
      END IF
C
99001 FORMAT (//,10X,'expected SECTION in input file: ',A,//)
C
      END
C*==section_find_keyword.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SECTION_FIND_KEYWORD(KW,FOUND)
C   ********************************************************************
C   *                                                                  *
C   *  find the keyword   KW(1:LKW)  in  LINE(1:LL)                    *
C   *  a blank is assumed to be before and after the KW in LINE        *
C   *  LINE is assumed to be upper case                                *
C   *  IPOS = 0    KW is not found                                     *
C   *  IPOS > 0    position of first character after KW//' '           *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_FILES,ONLY:LCMD,CMDUC,FOUND_KEYWORD,IPOS_KEYWORD
      IMPLICIT NONE
C*--SECTION_FIND_KEYWORD173
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SECTION_FIND_KEYWORD')
C
C Dummy arguments
C
      LOGICAL FOUND
      CHARACTER*(*) KW
C
C Local variables
C
      INTEGER I,IPOS,LKW,LKWP
      CHARACTER*82 KWP
C
C*** End of declarations rewritten by SPAG
C
      LKW = LEN_TRIM(KW)
      LCMD = LEN(CMDUC)
C
      IF ( LKW.GT.80 ) CALL STOP_MESSAGE(ROUTINE,'LKW > 80')
C
      KWP = ' '//KW(1:LKW)//' '
      LKWP = LKW + 2
C
      CALL STRING_CONVERT_TO_UC(KWP)
C
      I = INDEX(CMDUC(1:LCMD),KWP(1:LKWP))
C
      IF ( I.EQ.0 ) THEN
         IPOS = 0
         FOUND = .FALSE.
      ELSE
         IPOS = MIN(I+LKWP,LCMD)
         FOUND = .TRUE.
      END IF
C
      IPOS_KEYWORD = IPOS
C
      FOUND_KEYWORD = FOUND
C
      END
C*==section_set_integer.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SECTION_SET_INTEGER(STRX,X,DEFAULT,STOP_KEY)
C   ********************************************************************
C   *                                                                  *
C   *   find the string  STRX in  present section  CMDUC  and          *
C   *   update the   INTEGER   X    if STRX  was found                 *
C   *                                                                  *
C   *   set FOUND_INTEGER = .TRUE.                       ---> RETURN   *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *   STRX not found:                                                *
C   *                                                                  *
C   *   DEFAULT <> 9999     use supplied default value                 *
C   *                       FOUND_INTEGER = .FALSE.      ---> RETURN   *
C   *                                                                  *
C   *   DEFAULT  = 9999                                                *
C   *   STOP_KEY = 1                                     ---> STOP     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:CMDUC,STRINP,FOUND_INTEGER,IPOS_KEYWORD
      IMPLICIT NONE
C*--SECTION_SET_INTEGER256
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SECTION_SET_INTEGER')
C
C Dummy arguments
C
      INTEGER DEFAULT,STOP_KEY,X
      CHARACTER*(*) STRX
C
C Local variables
C
      LOGICAL FOUND,UPDATE
      INTEGER IPOS,LSTRX
C
C*** End of declarations rewritten by SPAG
C
      LSTRX = LEN_TRIM(STRX)
C
      CALL SECTION_FIND_KEYWORD(STRX,UPDATE)
C
      IPOS = IPOS_KEYWORD
C
      IF ( UPDATE ) THEN
         CALL STRING_GET_INTEGER(X,STRINP,CMDUC,IPOS,FOUND)
         IF ( .NOT.FOUND ) STOP 'no integer supplied for update !! '
      END IF
C
      FOUND_INTEGER = UPDATE
C
C=======================================================================
      IF ( FOUND_INTEGER ) RETURN
C=======================================================================
C
C
C-----------------------------------------------------------------------
C                     use DEFAULT if supplied
C-----------------------------------------------------------------------
C
      IF ( DEFAULT.NE.9999 ) THEN
C
         X = DEFAULT
C
         RETURN
C
      END IF
C
C-----------------------------------------------------------------------
C          NO SETTING POSSIBLE --- STOP if STOP_KEY is set
C-----------------------------------------------------------------------
C
      IF ( STOP_KEY.EQ.1 ) THEN
C
         WRITE (6,99001) STRX(1:LSTRX)
         CALL STOP_MESSAGE(ROUTINE,'requested INTEGER not found')
C
      END IF
C
99001 FORMAT (//,10X,'expected INTEGER in input file: ',A,//)
C
      END
C*==section_set_real.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SECTION_SET_REAL(STRX,X,DEFAULT,STOP_KEY)
C   ********************************************************************
C   *                                                                  *
C   *   find the string  STRX in  present section  CMDUC  and          *
C   *   update the   REAL      X    if STRX  was found                 *
C   *                                                                  *
C   *   set FOUND_REAL    = .TRUE.                       ---> RETURN   *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *   STRX not found:                                                *
C   *                                                                  *
C   *   DEFAULT <> 9999.0   use supplied default value                 *
C   *                       FOUND_REAL    = .FALSE.      ---> RETURN   *
C   *                                                                  *
C   *   DEFAULT  = 9999.0                                              *
C   *   STOP_KEY = 1                                     ---> STOP     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:CMDUC,STRINP,FOUND_REAL,IPOS_KEYWORD
      IMPLICIT NONE
C*--SECTION_SET_REAL358
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SECTION_SET_REAL')
C
C Dummy arguments
C
      REAL*8 DEFAULT,X
      INTEGER STOP_KEY
      CHARACTER*(*) STRX
C
C Local variables
C
      LOGICAL FOUND,UPDATE
      INTEGER IPOS,LSTRX
C
C*** End of declarations rewritten by SPAG
C
      LSTRX = LEN_TRIM(STRX)
C
      CALL SECTION_FIND_KEYWORD(STRX,UPDATE)
C
      IPOS = IPOS_KEYWORD
C
      IF ( UPDATE ) THEN
         CALL STRING_GET_REAL(X,STRINP,CMDUC,IPOS,FOUND)
         IF ( .NOT.FOUND ) STOP 'no  real supplied for update !! '
      END IF
C
      FOUND_REAL = UPDATE
C
C=======================================================================
      IF ( FOUND_REAL ) RETURN
C=======================================================================
C
C-----------------------------------------------------------------------
C                     use DEFAULT if supplied
C-----------------------------------------------------------------------
C
      IF ( ABS(DEFAULT-9999D0).GT.1D-5 ) THEN
C
         X = DEFAULT
C
         RETURN
C
      END IF
C
C-----------------------------------------------------------------------
C          NO SETTING POSSIBLE --- STOP if STOP_KEY is set
C-----------------------------------------------------------------------
C
      IF ( STOP_KEY.EQ.1 ) THEN
C
         WRITE (6,99001) STRX(1:LSTRX)
         CALL STOP_MESSAGE(ROUTINE,'requested REAL not found')
C
      END IF
C
99001 FORMAT (//,10X,'expected REAL in input file: ',A,//)
C
      END
C*==section_set_string.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SECTION_SET_STRING(STRX,X,DEFAULT,STOP_KEY)
C   ********************************************************************
C   *                                                                  *
C   *   find the string  STRX in  present section  CMDUC  and          *
C   *   update the   STRING    X    if STRX  was found                 *
C   *                                                                  *
C   *   set FOUND_STRING  = .TRUE.                       ---> RETURN   *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *   STRX not found:                                                *
C   *                                                                  *
C   *   DEFAULT <> '9999'   use supplied default value                 *
C   *                       FOUND_STRING  = .FALSE.      ---> RETURN   *
C   *                                                                  *
C   *   DEFAULT  = '9999'                                              *
C   *   STOP_KEY = 1                                     ---> STOP     *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_FILES,ONLY:CMD00,FOUND_STRING,IPOS_KEYWORD
      IMPLICIT NONE
C*--SECTION_SET_STRING461
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SECTION_SET_STRING')
C
C Dummy arguments
C
      CHARACTER*(*) DEFAULT,STRX,X
      INTEGER STOP_KEY
C
C Local variables
C
      INTEGER IPOS,LSTRX,LX
      LOGICAL UPDATE
C
C*** End of declarations rewritten by SPAG
C
      LSTRX = LEN_TRIM(STRX)
C
      CALL SECTION_FIND_KEYWORD(STRX,UPDATE)
C
      IPOS = IPOS_KEYWORD
C
      IF ( UPDATE ) THEN
         CALL STRING_GET_STRING(X,LX,CMD00,IPOS)
         IF ( LX.EQ.0 ) STOP 'no  string supplied for update !! '
      END IF
C
      FOUND_STRING = UPDATE
C
C=======================================================================
      IF ( FOUND_STRING ) RETURN
C=======================================================================
C
C-----------------------------------------------------------------------
C                     use DEFAULT if supplied
C-----------------------------------------------------------------------
C
      IF ( DEFAULT.NE.'9999' ) THEN
C
         X = DEFAULT
C
         RETURN
C
      END IF
C
C-----------------------------------------------------------------------
C          NO SETTING POSSIBLE --- STOP if STOP_KEY is set
C-----------------------------------------------------------------------
C
      IF ( STOP_KEY.EQ.1 ) THEN
C
         WRITE (6,99001) STRX(1:LSTRX)
         CALL STOP_MESSAGE(ROUTINE,'requested STRING not found')
C
      END IF
C
99001 FORMAT (//,10X,'expected STRING in input file: ',A,//)
C
      END
C*==section_set_integer_array.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SECTION_SET_INTEGER_ARRAY(STRX,X,NX,NXMAX,K_REQUIRED,
     &   DEFAULT,STOP_KEY)
C   ********************************************************************
C   *                                                                  *
C   *   find the string  STRX in  present section  CMDUC  and          *
C   *   update the   INTEGER ARRAY  X    if STRX  was found            *
C   *                                                                  *
C   *   K_REQUIRED = 0     number NX of values found is returned       *
C   *                      if NX = 1: copy value to whole array        *
C   *                                                                  *
C   *   K_REQUIRED = 1     NX values (input) have to be found          *
C   *                      if less than NX values are found the        *
C   *                      program stops                               *
C   *                                                                  *
C   *   set FOUND_INTEGER_ARRAY = .TRUE.                 ---> RETURN   *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *   STRX not found:                                                *
C   *                                                                  *
C   *   DEFAULT <> 9999     use supplied default value for whole array *
C   *                       FOUND_INTEGER_ARRAY = .FALSE.---> RETURN   *
C   *                                                                  *
C   *   DEFAULT  = 9999                                                *
C   *   STOP_KEY = 1                                     ---> STOP     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:CMDUC,FOUND_INTEGER_ARRAY,IPOS_KEYWORD,STRINP
      IMPLICIT NONE
C*--SECTION_SET_INTEGER_ARRAY570
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='INPUT_SET_INTEGER_ARRAY')
C
C Dummy arguments
C
      INTEGER DEFAULT,K_REQUIRED,NX,NXMAX,STOP_KEY
      CHARACTER*(*) STRX
      INTEGER X(*)
C
C Local variables
C
      CHARACTER*1 C,CDEL
      LOGICAL FOUND,UPDATE
      INTEGER I,IPOS,ITOP,J,LCMD,NFOUND,NX_REQUIRED
C
C*** End of declarations rewritten by SPAG
C
      IF ( K_REQUIRED.EQ.1 ) NX_REQUIRED = NX
C
      LCMD = LEN_TRIM(CMDUC) + 1
C
      CALL SECTION_FIND_KEYWORD(STRX,UPDATE)
C
      IPOS = IPOS_KEYWORD
C
      NFOUND = 0
C
      IF ( UPDATE ) THEN
C
         CDEL = ' '
 50      CONTINUE
         C = CMDUC(IPOS:IPOS)
         IF ( C.NE.' ' ) THEN
            IF ( C.EQ.'{' ) CDEL = '}'
            IF ( C.EQ.'[' ) CDEL = ']'
            IF ( C.EQ.'(' ) CDEL = ')'
         ELSE IF ( IPOS.LT.LCMD ) THEN
            IPOS = IPOS + 1
            GOTO 50
         ELSE
            WRITE (6,*) ROUTINE//':  IPOS > LCMD '
            WRITE (6,*) ' CMDUC   ',CMDUC(1:LCMD)
            WRITE (6,*) '         ',(' ',I=1,IPOS-1),'^'
            WRITE (6,*) ' STR     ',STRX
         END IF
C
         IF ( CDEL.EQ.' ' ) THEN
            ITOP = 1
         ELSE
            ITOP = NXMAX
         END IF
C
         DO I = 1,ITOP
            CALL STRING_GET_INTEGER(X(I),STRINP,CMDUC,IPOS,FOUND)
            IF ( CMDUC(IPOS:IPOS).EQ.CDEL .OR. (ITOP.EQ.1) ) THEN
               NFOUND = I
               DO J = NFOUND + 1,NXMAX
                  X(J) = X(NFOUND)
               END DO
               EXIT
            END IF
         END DO
C
      END IF
C
      FOUND_INTEGER_ARRAY = UPDATE
C
C=======================================================================
      IF ( FOUND_INTEGER_ARRAY ) THEN
C
         IF ( (K_REQUIRED.EQ.1) .AND. (NFOUND.NE.NX_REQUIRED) ) THEN
            WRITE (6,99002) NFOUND,NX_REQUIRED
            CALL STOP_MESSAGE(ROUTINE,
     &                        'requested number of values not found')
         END IF
C
         IF ( K_REQUIRED.EQ.0 ) NX = NFOUND
C
         IF ( NX.EQ.1 ) X(2:NXMAX) = X(1)
C
         RETURN
      END IF
C=======================================================================
C
C
C-----------------------------------------------------------------------
C                     use DEFAULT if supplied
C-----------------------------------------------------------------------
C
      IF ( DEFAULT.NE.9999 ) THEN
C
         X(1:NXMAX) = DEFAULT
C
         RETURN
C
      END IF
C
C-----------------------------------------------------------------------
C          NO SETTING POSSIBLE --- STOP if STOP_KEY is set
C-----------------------------------------------------------------------
C
      IF ( STOP_KEY.EQ.1 ) THEN
C
         WRITE (6,99001) STRX(1:LEN_TRIM(STRX))
         CALL STOP_MESSAGE(ROUTINE,'requested INTEGER ARRAY not found')
C
      END IF
C
99001 FORMAT (//,10X,'expected INTEGER in input file: ',A,//)
99002 FORMAT (//,10X,'number of values found     NX          ',I5,/,10X,
     &        'number of values required  NX_REQUIRED ',I5,/)
      END
C*==section_set_real_array.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SECTION_SET_REAL_ARRAY(STRX,X,NX,NXMAX,K_REQUIRED,
     &                                  DEFAULT,STOP_KEY)
C   ********************************************************************
C   *                                                                  *
C   *   find the string  STRX in  present section  CMDUC  and          *
C   *   update the   REAL ARRAY     X    if STRX  was found            *
C   *                                                                  *
C   *   K_REQUIRED = 0     number NX of values found is returned       *
C   *                      if NX = 1: copy value to whole array        *
C   *                                                                  *
C   *   K_REQUIRED = 1     NX values (input) have to be found          *
C   *                      if less than NX values are found the        *
C   *                      program stops                               *
C   *                                                                  *
C   *   set FOUND_REAL_ARRAY = .TRUE.                    ---> RETURN   *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *   STRX not found:                                                *
C   *                                                                  *
C   *   DEFAULT <> 9999     use supplied default value for whole array *
C   *                       FOUND_REAL_ARRAY = .FALSE.   ---> RETURN   *
C   *                                                                  *
C   *   DEFAULT  = 9999                                                *
C   *   STOP_KEY = 1                                     ---> STOP     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:CMDUC,FOUND_REAL_ARRAY,IPOS_KEYWORD,STRINP
      IMPLICIT NONE
C*--SECTION_SET_REAL_ARRAY733
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='INPUT_SET_REAL_ARRAY')
C
C Dummy arguments
C
      REAL*8 DEFAULT
      INTEGER K_REQUIRED,NX,NXMAX,STOP_KEY
      CHARACTER*(*) STRX
      REAL*8 X(*)
C
C Local variables
C
      CHARACTER*1 C,CDEL
      LOGICAL FOUND,UPDATE
      INTEGER I,IPOS,ITOP,J,LCMD,NFOUND,NX_REQUIRED
C
C*** End of declarations rewritten by SPAG
C
      IF ( K_REQUIRED.EQ.1 ) NX_REQUIRED = NX
C
      LCMD = LEN_TRIM(CMDUC) + 1
C
      CALL SECTION_FIND_KEYWORD(STRX,UPDATE)
C
      IPOS = IPOS_KEYWORD
C
      NFOUND = 0
C
      IF ( UPDATE ) THEN
C
         CDEL = ' '
 50      CONTINUE
         C = CMDUC(IPOS:IPOS)
         IF ( C.NE.' ' ) THEN
            IF ( C.EQ.'{' ) CDEL = '}'
            IF ( C.EQ.'[' ) CDEL = ']'
            IF ( C.EQ.'(' ) CDEL = ')'
         ELSE IF ( IPOS.LT.LCMD ) THEN
            IPOS = IPOS + 1
            GOTO 50
         ELSE
            WRITE (6,*) ROUTINE//':  IPOS > LCMD '
            WRITE (6,*) ' CMDUC   ',CMDUC(1:LCMD)
            WRITE (6,*) '         ',(' ',I=1,IPOS-1),'^'
            WRITE (6,*) ' STR     ',STRX
         END IF
C
         IF ( CDEL.EQ.' ' ) THEN
            ITOP = 1
         ELSE
            ITOP = NXMAX
         END IF
C
         DO I = 1,ITOP
            CALL STRING_GET_REAL(X(I),STRINP,CMDUC,IPOS,FOUND)
            IF ( CMDUC(IPOS:IPOS).EQ.CDEL .OR. (ITOP.EQ.1) ) THEN
               NFOUND = I
               DO J = NFOUND + 1,NXMAX
                  X(J) = X(NFOUND)
               END DO
               EXIT
            END IF
         END DO
C
      END IF
C
      FOUND_REAL_ARRAY = UPDATE
C
C=======================================================================
      IF ( FOUND_REAL_ARRAY ) THEN
C
         IF ( (K_REQUIRED.EQ.1) .AND. (NFOUND.NE.NX_REQUIRED) ) THEN
            WRITE (6,99002) NFOUND,NX_REQUIRED
            CALL STOP_MESSAGE(ROUTINE,
     &                        'requested number of values not found')
         END IF
C
         IF ( K_REQUIRED.EQ.0 ) NX = NFOUND
C
         IF ( NX.EQ.1 ) X(2:NXMAX) = X(1)
C
         RETURN
      END IF
C=======================================================================
C
C
C-----------------------------------------------------------------------
C                     use DEFAULT if supplied
C-----------------------------------------------------------------------
C
      IF ( ABS(DEFAULT-9999D0).GT.1D-5 ) THEN
C
         X(1:NXMAX) = DEFAULT
C
         RETURN
C
      END IF
C
C-----------------------------------------------------------------------
C          NO SETTING POSSIBLE --- STOP if STOP_KEY is set
C-----------------------------------------------------------------------
C
      IF ( STOP_KEY.EQ.1 ) THEN
C
         WRITE (6,99001) STRX(1:LEN_TRIM(STRX))
         CALL STOP_MESSAGE(ROUTINE,'requested REAL ARRAY not found')
C
      END IF
C
99001 FORMAT (//,10X,'expected REAL in input file: ',A,//)
99002 FORMAT (//,10X,'number of values found     NX          ',I5,/,10X,
     &        'number of values required  NX_REQUIRED ',I5,/)
      END
C*==section_set_string_array.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SECTION_SET_STRING_ARRAY(STRX,X,NX,NXMAX,K_REQUIRED,
     &   DEFAULT,STOP_KEY)
C   ********************************************************************
C   *                                                                  *
C   *   find the string  STRX in  present section  CMDUC  and          *
C   *   update the   STRING ARRAY     X    if STRX  was found          *
C   *                                                                  *
C   *   K_REQUIRED = 0     number NX of values found is returned       *
C   *                      if NX = 1: copy value to whole array        *
C   *                                                                  *
C   *   K_REQUIRED = 1     NX values (input) have to be found          *
C   *                      if less than NX values are found the        *
C   *                      program stops                               *
C   *                                                                  *
C   *   set FOUND_STRING_ARRAY = .TRUE.                  ---> RETURN   *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *   STRX not found:                                                *
C   *                                                                  *
C   *   DEFAULT <> 9999     use supplied default value for whole array *
C   *                       FOUND_STRING_ARRAY = .FALSE. ---> RETURN   *
C   *                                                                  *
C   *   DEFAULT  = 9999                                                *
C   *   STOP_KEY = 1                                     ---> STOP     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:CMDUC,FOUND_STRING_ARRAY,IPOS_KEYWORD
      IMPLICIT NONE
C*--SECTION_SET_STRING_ARRAY897
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='INPUT_SET_STRING_ARRAY')
C
C Dummy arguments
C
      CHARACTER*(*) DEFAULT,STRX
      INTEGER K_REQUIRED,NX,NXMAX,STOP_KEY
      CHARACTER*(*) X(*)
C
C Local variables
C
      CHARACTER*1 C,CDEL
      INTEGER I,IPOS,ITOP,J,LCMD,LX,NFOUND,NX_REQUIRED
      LOGICAL UPDATE
C
C*** End of declarations rewritten by SPAG
C
      IF ( K_REQUIRED.EQ.1 ) NX_REQUIRED = NX
C
      LCMD = LEN_TRIM(CMDUC) + 1
C
      CALL SECTION_FIND_KEYWORD(STRX,UPDATE)
C
      IPOS = IPOS_KEYWORD
C
      NFOUND = 0
C
      IF ( UPDATE ) THEN
C
         CDEL = ' '
 50      CONTINUE
         C = CMDUC(IPOS:IPOS)
         IF ( C.NE.' ' ) THEN
            IF ( C.EQ.'{' ) CDEL = '}'
            IF ( C.EQ.'[' ) CDEL = ']'
            IF ( C.EQ.'(' ) CDEL = ')'
         ELSE IF ( IPOS.LT.LCMD ) THEN
            IPOS = IPOS + 1
            GOTO 50
         ELSE
            WRITE (6,*) ROUTINE//':  IPOS > LCMD '
            WRITE (6,*) ' CMDUC   ',CMDUC(1:LCMD)
            WRITE (6,*) '         ',(' ',I=1,IPOS-1),'^'
            WRITE (6,*) ' STR     ',STRX
         END IF
C
         IF ( CDEL.EQ.' ' ) THEN
            ITOP = 1
         ELSE
            ITOP = NXMAX
         END IF
C
         DO I = 1,ITOP
            CALL STRING_GET_STRING(X(I),LX,CMDUC,IPOS)
            IF ( CMDUC(IPOS:IPOS).EQ.CDEL .OR. (ITOP.EQ.1) ) THEN
               NFOUND = I
               DO J = NFOUND + 1,NXMAX
                  X(J) = X(NFOUND)
               END DO
               EXIT
            END IF
         END DO
C
      END IF
C
      FOUND_STRING_ARRAY = UPDATE
C
C=======================================================================
      IF ( FOUND_STRING_ARRAY ) THEN
C
         IF ( (K_REQUIRED.EQ.1) .AND. (NFOUND.NE.NX_REQUIRED) ) THEN
            WRITE (6,99002) NFOUND,NX_REQUIRED
            CALL STOP_MESSAGE(ROUTINE,
     &                        'requested number of values not found')
         END IF
C
         IF ( K_REQUIRED.EQ.0 ) NX = NFOUND
C
         IF ( NX.EQ.1 ) X(2:NXMAX) = X(1)
C
         RETURN
      END IF
C=======================================================================
C
C
C-----------------------------------------------------------------------
C                     use DEFAULT if supplied
C-----------------------------------------------------------------------
C
      IF ( DEFAULT.NE.'9999' ) THEN
C
         X(1:NXMAX) = DEFAULT
C
         RETURN
C
      END IF
C
C-----------------------------------------------------------------------
C          NO SETTING POSSIBLE --- STOP if STOP_KEY is set
C-----------------------------------------------------------------------
C
      IF ( STOP_KEY.EQ.1 ) THEN
C
         WRITE (6,99001) STRX(1:LEN_TRIM(STRX))
         CALL STOP_MESSAGE(ROUTINE,'requested STRING ARRAY not found')
C
      END IF
C
99001 FORMAT (//,10X,'expected STRING in input file: ',A,//)
99002 FORMAT (//,10X,'number of values found     NX          ',I5,/,10X,
     &        'number of values required  NX_REQUIRED ',I5,/)
      END
C*==string_get_integer.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE STRING_GET_INTEGER(I,STR,STRING,IPOS,FOUND)
C   ********************************************************************
C   *                                                                  *
C   *   get the next integer      I      in  STRING(IPOS:LSTRING)      *
C   *   delimiters:  blank , ; { } [ ] ( )  ignored                    *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--STRING_GET_INTEGER1039
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL FOUND
      INTEGER I,IPOS
      CHARACTER*(*) STR,STRING
C
C Local variables
C
      INTEGER LSTR
C
C*** End of declarations rewritten by SPAG
C
      CALL STRING_GET_STRING(STR,LSTR,STRING,IPOS)
C
      IF ( LSTR.EQ.0 ) THEN
         FOUND = .FALSE.
      ELSE
         READ (STR,FMT=*,ERR=100) I
         FOUND = .TRUE.
      END IF
      RETURN
 100  CONTINUE
      WRITE (6,*) 
     &     ' TROUBLE with conversion to INTEGER in <STRING_GET_INTEGER>'
      WRITE (6,*) ' STR = ',STR(1:LSTR)
      STOP
      END
C*==string_get_real.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE STRING_GET_REAL(R,STR,STRING,IPOS,FOUND)
C   ********************************************************************
C   *                                                                  *
C   *   get the next real number       R      in  STRING(IPOS:LSTRING) *
C   *   delimiters:  blank , ; { } [ ] ( )  ignored                    *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--STRING_GET_REAL1091
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL FOUND
      INTEGER IPOS
      REAL*8 R
      CHARACTER*(*) STR,STRING
C
C Local variables
C
      INTEGER LSTR
C
C*** End of declarations rewritten by SPAG
C
      CALL STRING_GET_STRING(STR,LSTR,STRING,IPOS)
C
      IF ( LSTR.EQ.0 ) THEN
         FOUND = .FALSE.
      ELSE
         READ (STR,FMT=*,ERR=100) R
         FOUND = .TRUE.
      END IF
      RETURN
 100  CONTINUE
      WRITE (6,*) 
     &           ' TROUBLE with conversion to REAL in <STRING_GET_REAL>'
      WRITE (6,*) ' STR = ',STR(1:LSTR)
      STOP
      END
C*==string_get_string.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE STRING_GET_STRING(STR,LSTR,STRING,IPOS)
C   ********************************************************************
C   *                                                                  *
C   *   get the next string  STR  in STRING(IPOS:LSTRING)              *
C   *   delimiters:  blank , ; { } [ ] ( )  ignored                    *
C   *   push IPOS to position of next non-blank character in STRING    *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--STRING_GET_STRING1145
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IPOS,LSTR
      CHARACTER*(*) STR,STRING
C
C Local variables
C
      INTEGER I,I0,J,LEND,LSTRING,LSTRMAX
C
C*** End of declarations rewritten by SPAG
C
      LSTRING = LEN(STRING)
      LSTRMAX = LEN(STR)
C
      DO I = 1,LSTRMAX
         STR(I:I) = ' '
      END DO
C
      I0 = IPOS
      LSTR = 0
      DO I = I0,LSTRING
         IF ( STRING(I:I).NE.' ' .AND. STRING(I:I).NE.'=' .AND. 
     &        STRING(I:I).NE.',' .AND. STRING(I:I).NE.';' .AND. 
     &        STRING(I:I).NE.'{' .AND. STRING(I:I).NE.'[' .AND. 
     &        STRING(I:I).NE.'(' ) THEN
            DO J = I,LSTRING
               IF ( STRING(J:J).EQ.' ' .OR. STRING(J:J).EQ.',' .OR. 
     &              STRING(J:J).EQ.';' .OR. STRING(J:J).EQ.'}' .OR. 
     &              STRING(J:J).EQ.']' .OR. STRING(J:J).EQ.')' ) EXIT
               LSTR = LSTR + 1
               IF ( LSTR.GT.LSTRMAX ) THEN
                  WRITE (6,'(A,/,A,A,/,A,I3)')
     &                    ' ERROR in <STRING_GET_STRING> ',
     &                   ' STRING:             ',STRING(1:LSTRING),
     &                   ' STR too long LSTRMAX:',LSTRMAX
                  STOP
               END IF
               STR(LSTR:LSTR) = STRING(J:J)
               IPOS = J + 1
            END DO
            GOTO 100
         END IF
      END DO
      RETURN
C
 100  CONTINUE
      LEND = LEN_TRIM(STRING)
      IF ( STRING(IPOS:IPOS).EQ.' ' .AND. IPOS.LT.LEND ) THEN
         IPOS = IPOS + 1
         GOTO 100
      END IF
C
      END
C*==string_convert_to_uc.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE STRING_CONVERT_TO_UC(STRING)
C   ********************************************************************
C   *                                                                  *
C   *  convert characters to upper case in  STRING(I1:I2)              *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--STRING_CONVERT_TO_UC1222
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*(*) STRING
C
C Local variables
C
      INTEGER I,IA,INC,IZ,J
C
C*** End of declarations rewritten by SPAG
C
      IA = ICHAR('a')
      IZ = ICHAR('z')
      INC = ICHAR('A') - IA
C
      DO I = 1,LEN_TRIM(STRING)
         J = ICHAR(STRING(I:I))
         IF ( J.GE.IA .AND. J.LE.IZ ) STRING(I:I) = CHAR(J+INC)
      END DO
      END
C*==string_convert_to_lc.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE STRING_CONVERT_TO_LC(STRING)
C   ********************************************************************
C   *                                                                  *
C   *  convert characters to lower case in  STRING                     *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--STRING_CONVERT_TO_LC1265
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*(*) STRING
C
C Local variables
C
      INTEGER I,IA,INC,IZ,J
C
C*** End of declarations rewritten by SPAG
C
      IA = ICHAR('A')
      IZ = ICHAR('Z')
      INC = ICHAR('a') - IA
C
      DO I = 1,LEN_TRIM(STRING)
         J = ICHAR(STRING(I:I))
         IF ( J.GE.IA .AND. J.LE.IZ ) STRING(I:I) = CHAR(J+INC)
      END DO
      END
C*==string_trim_left.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE STRING_TRIM_LEFT(STRING)
C   ********************************************************************
C   *                                                                  *
C   *  remove leading blanks in   STRING                               *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--STRING_TRIM_LEFT1307
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*(*) STRING
C
C Local variables
C
      INTEGER I,J,K,LL,LX
C
C*** End of declarations rewritten by SPAG
C
      LX = LEN(STRING)
      LL = LEN_TRIM(STRING)
C
      DO I = 1,LL
         IF ( STRING(I:I).NE.' ' ) THEN
            IF ( I.EQ.1 ) RETURN
            K = 0
            DO J = I,LL
               K = K + 1
               STRING(K:K) = STRING(J:J)
            END DO
            LL = LL - (I-1) + 1
            STRING(LL:LX) = ' '
            RETURN
         END IF
      END DO
      END
C*==string_add_n.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE STRING_ADD_N(STRING,N)
C   ********************************************************************
C   *                                                                  *
C   *  add the integer N to     STRING                                 *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--STRING_ADD_N1358
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='STRING_ADD_N')
C
C Dummy arguments
C
      INTEGER N
      CHARACTER*(*) STRING
C
C Local variables
C
      INTEGER I,I0,J,LL,LLMAX,M,NDIG
      REAL*8 X
C
C*** End of declarations rewritten by SPAG
C
      LLMAX = LEN(STRING)
      LL = LEN_TRIM(STRING)
C
      IF ( N.LT.0 ) THEN
         M = -N
         LL = LL + 1
         STRING(LL:LL) = '-'
      ELSE
         M = N
      END IF
C
      DO I = 1,10000
         IF ( M.LT.10 ) THEN
            NDIG = I
            EXIT
         END IF
         M = M/10
      END DO
      X = DBLE(ABS(N))
      I0 = ICHAR('0')
      DO I = NDIG,1, - 1
         LL = LL + 1
         IF ( LL.GT.LLMAX )
     &         CALL STOP_MESSAGE(ROUTINE,'string gets too long:  '//
     &        STRING)
         J = INT(X/10.0D0**(I-1)) - INT(X/10.0D0**I)*10
         STRING(LL:LL) = CHAR(I0+J)
      END DO
      END
C*==input_get_int2str.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE INPUT_GET_INT2STR(MYSTR,STRL,POS,INTG)
C   ********************************************************************
C   *                                                                  *
C   *  INPUT:     intg   ... integer to write into mystr               *
C   *             strl   ... length of mystr                           *
C   *             pos    ... rightmost position of number              *
C   *                                                                  *
C   *  IN/OUT  :  mystr                                                *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--INPUT_GET_INT2STR1435
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER INTG,POS,STRL
      CHARACTER*(*) MYSTR
C
C Local variables
C
      CHARACTER*4 FMTSTR
      INTEGER NDIGITS
C
C*** End of declarations rewritten by SPAG
C
      FMTSTR = '(I1)'
C
      IF ( INTG.LT.0 ) THEN
         WRITE (6,'("int2str:: intg lesser than 0 -- STOP  ")')
         STOP
      END IF
C
      IF ( INTG.EQ.0 ) THEN
         NDIGITS = 1
      ELSE
         NDIGITS = INT(LOG10(DBLE(INTG))) + 1
      END IF
C
      WRITE (FMTSTR(3:3),'(I1)') NDIGITS
C
      IF ( POS.GT.STRL ) THEN
         WRITE (6,'("int2str:: pos larger the strl -- STOP  ")')
         STOP
      END IF
C
      IF ( NDIGITS.GT.POS ) THEN
         WRITE (6,'("int2str:: intg too large  -- STOP")')
         STOP
      END IF
C
      WRITE (MYSTR(POS-NDIGITS+1:POS),FMTSTR) INTG
C
      END
C*==readkwint.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE READKWINT(IFILINP,KW,I,DEF,IWRI,IFLAG)
C   ********************************************************************
C   *                                                                  *
C   *  find the keyword   KW (STR*10)  in  file   IFILINP              *
C   *  at the beginning of a line and read following  INTEGER I        *
C   *  DEF        default value                                        *
C   *  IWRI  = 1  write result to chanel IWRI                          *
C   *          0  no output                                            *
C   *  IFLAG = 1  stop if KW not found                                 *
C   *          0  continue even if KW not found and set to default     *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--READKWINT1505
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DEF,I,IFILINP,IFLAG,IWRI
      CHARACTER*10 KW
C
C Local variables
C
      CHARACTER*80 LINE
      CHARACTER*10 STR10
C
C*** End of declarations rewritten by SPAG
C
      REWIND IFILINP
 100  CONTINUE
      READ (IFILINP,'(A)',ERR=200,END=200) LINE
      CALL STRING_REPLACE_TAB(LINE)
      CALL STRING_TRIM_LEFT(LINE)
      IF ( INDEX(LINE,KW).NE.1 ) GOTO 100
C
      READ (LINE,*,ERR=300,END=300) STR10,I
      IF ( IWRI.NE.0 ) WRITE (IWRI,99001) KW,IFILINP,STR10,I
      RETURN
C
 200  CONTINUE
      IF ( IFLAG.NE.0 ) THEN
         WRITE (6,99002) KW,IFILINP
         STOP
      ELSE
         I = DEF
         IF ( IWRI.NE.0 ) WRITE (IWRI,99001) KW,IFILINP,STR10,I
         REWIND IFILINP
         RETURN
      END IF
C
 300  CONTINUE
      WRITE (6,99003) IFILINP,KW,LINE
      STOP
99001 FORMAT (2X,'keyword ',A,' file',I4,' INTEGER  input: ',A,I5)
99002 FORMAT (//,1X,79('*'),/,2X,'keyword  ',A,
     &        ' not found in file IFILINP',I4,/,2X,
     &        'the STOP flag was set indicating input obligatory',/,
     &        79('*'),/)
99003 FORMAT (//,1X,79('*'),/,2X,'error reading file',I4,
     &        ' for keyword ',A,/,A,/,1X,79('*'),/)
      END
C*==readkwreal.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE READKWREAL(IFILINP,KW,A,DEF,IWRI,IFLAG)
C   ********************************************************************
C   *                                                                  *
C   *  find the keyword   KW (STR*10)  in  file   IFILINP              *
C   *  at the beginning of a line and read following  REAL    A        *
C   *  DEF        default value                                        *
C   *  IWRI  = 1  write result to chanel IWRI                          *
C   *          0  no output                                            *
C   *  IFLAG = 1  stop if KW not found                                 *
C   *          0  continue even if KW not found and set to default     *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--READKWREAL1580
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 A,DEF
      INTEGER IFILINP,IFLAG,IWRI
      CHARACTER*10 KW
C
C Local variables
C
      CHARACTER*80 LINE
      CHARACTER*10 STR10
C
C*** End of declarations rewritten by SPAG
C
      REWIND IFILINP
 100  CONTINUE
      READ (IFILINP,'(A)',ERR=200,END=200) LINE
      CALL STRING_REPLACE_TAB(LINE)
      CALL STRING_TRIM_LEFT(LINE)
      IF ( INDEX(LINE,KW).NE.1 ) GOTO 100
C
      READ (LINE,*,ERR=300,END=300) STR10,A
      IF ( IWRI.NE.0 ) WRITE (IWRI,99001) KW,IFILINP,STR10,A
      RETURN
C
 200  CONTINUE
      IF ( IFLAG.NE.0 ) THEN
         WRITE (6,99002) KW,IFILINP
         STOP
      ELSE
         A = DEF
         IF ( IWRI.NE.0 ) WRITE (IWRI,99001) KW,IFILINP,STR10,A
         REWIND IFILINP
         RETURN
      END IF
C
 300  CONTINUE
      WRITE (6,99003) IFILINP,KW
      STOP
99001 FORMAT (2X,'keyword ',A,' file',I4,' REAL     input: ',A,E15.8)
99002 FORMAT (//,1X,79('*'),/,2X,'keyword  ',A,
     &        ' not found in file IFILINP',I4,/,2X,
     &        'the STOP flag was set indicating input obligatory',/,
     &        79('*'),/)
99003 FORMAT (//,1X,79('*'),/,2X,'error reading file',I4,
     &        ' for keyword ',A,/,1X,79('*'),/)
      END
C*==readkwlog.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE READKWLOG(IFILINP,KW,A,DEF,IWRI,IFLAG)
C   ********************************************************************
C   *                                                                  *
C   *  find the keyword   KW (STR*10)  in  file   IFILINP                 *
C   *  at the beginning of a line and read following  LOGICAL A        *
C   *  DEF        default value                                        *
C   *  IWRI  = 1  write result to chanel IWRI                          *
C   *          0  no output                                            *
C   *  IFLAG = 1  stop if KW not found                                 *
C   *          0  continue even if KW not found and set to default     *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--READKWLOG1656
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL A,DEF
      INTEGER IFILINP,IFLAG,IWRI
      CHARACTER*10 KW
C
C Local variables
C
      CHARACTER*80 LINE
      CHARACTER*10 STR10
C
C*** End of declarations rewritten by SPAG
C
      REWIND IFILINP
 100  CONTINUE
      READ (IFILINP,'(A)',ERR=200,END=200) LINE
      CALL STRING_REPLACE_TAB(LINE)
      CALL STRING_TRIM_LEFT(LINE)
      IF ( INDEX(LINE,KW).NE.1 ) GOTO 100
C
      READ (LINE,*,ERR=300,END=300) STR10,A
      IF ( IWRI.NE.0 ) WRITE (IWRI,99001) KW,IFILINP,STR10,A
      RETURN
C
 200  CONTINUE
      IF ( IFLAG.NE.0 ) THEN
         WRITE (6,99002) KW,IFILINP
         STOP
      ELSE
         A = DEF
         IF ( IWRI.NE.0 ) WRITE (IWRI,99001) KW,IFILINP,STR10,A
         REWIND IFILINP
         RETURN
      END IF
C
 300  CONTINUE
      WRITE (6,99003) IFILINP,KW
      STOP
99001 FORMAT (2X,'keyword ',A,' file',I4,' LOGICAL  input: ',A,L5)
99002 FORMAT (//,1X,79('*'),/,2X,'keyword  ',A,
     &        ' not found in file IFILINP',I4,/,2X,
     &        'the STOP flag was set indicating input obligatory',/,
     &        79('*'),/)
99003 FORMAT (//,1X,79('*'),/,2X,'error reading file',I4,
     &        ' for keyword ',A,/,1X,79('*'),/)
      END
C*==readkwstr.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE READKWSTR(IFILINP,KW,STR,DEF,L,IWRI,IFLAG)
C   ********************************************************************
C   *                                                                  *
C   *  find the keyword   KW (STR*10)  in  file   IFILINP                 *
C   *  at the beginning of a line and read following  STRING STR*(L)   *
C   *  DEF        default value                                        *
C   *  IWRI  = 1  write result to chanel IWRI                          *
C   *          0  no output                                            *
C   *  IFLAG = 1  stop if KW not found                                 *
C   *          0  continue even if KW not found and set to default     *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--READKWSTR1732
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*(*) DEF,STR
      INTEGER IFILINP,IFLAG,IWRI,L
      CHARACTER*10 KW
C
C Local variables
C
      CHARACTER*80 LINE
      CHARACTER*10 STR10
      CHARACTER*70 STR70
C
C*** End of declarations rewritten by SPAG
C
      REWIND IFILINP
 100  CONTINUE
      READ (IFILINP,'(A)',ERR=200,END=200) LINE
      CALL STRING_REPLACE_TAB(LINE)
      CALL STRING_TRIM_LEFT(LINE)
      IF ( INDEX(LINE,KW).NE.1 ) GOTO 100
C
      READ (LINE,*,ERR=300,END=300) STR10,STR70
      STR = STR70(1:MIN(L,70))
      IF ( IWRI.NE.0 ) WRITE (IWRI,99001) KW,IFILINP,STR10,STR
      RETURN
C
 200  CONTINUE
      IF ( IFLAG.NE.0 ) THEN
         WRITE (6,99002) KW,IFILINP
         STOP
      ELSE
         STR = DEF
         IF ( IWRI.NE.0 ) WRITE (IWRI,99001) KW,IFILINP,STR10,STR
         REWIND IFILINP
         RETURN
      END IF
C
 300  CONTINUE
      WRITE (6,99003) IFILINP,KW
      STOP
99001 FORMAT (2X,'keyword ',A,' file',I4,' STRING   input: ',2A)
99002 FORMAT (//,1X,79('*'),/,2X,'keyword  ',A,
     &        ' not found in file IFILINP',I4,/,2X,
     &        'the STOP flag was set indicating input obligatory',/,
     &        79('*'),/)
99003 FORMAT (//,1X,79('*'),/,2X,'error reading file',I4,
     &        ' for keyword ',A,/,1X,79('*'),/)
      END
C*==readkwrarr.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE READKWRARR(IFILINP,KW,A,L,DEF,IWRI,IFLAG)
C   ********************************************************************
C   *                                                                  *
C   *  find the keyword   KW (STR*10)  in  file   IFILINP                 *
C   *  at the beginning of a line and read following  REAL ARRAY A(L)  *
C   *  DEF        default value                                        *
C   *  IWRI  = 1  write result to chanel IWRI                          *
C   *          0  no output                                            *
C   *  IFLAG = 1  stop if KW not found                                 *
C   *          0  continue even if KW not found and set to default     *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--READKWRARR1810
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 DEF
      INTEGER IFILINP,IFLAG,IWRI,L
      CHARACTER*10 KW
      REAL*8 A(L)
C
C Local variables
C
      CHARACTER*80 LINE
      CHARACTER*10 STR10
C
C*** End of declarations rewritten by SPAG
C
      REWIND IFILINP
 100  CONTINUE
      READ (IFILINP,'(A)',ERR=200,END=200) LINE
      CALL STRING_REPLACE_TAB(LINE)
      CALL STRING_TRIM_LEFT(LINE)
      IF ( INDEX(LINE,KW).NE.1 ) GOTO 100
C
      READ (LINE,*,ERR=300,END=300) STR10,A
      IF ( IWRI.NE.0 ) WRITE (IWRI,99001) KW,IFILINP,STR10,A
      RETURN
C
 200  CONTINUE
      IF ( IFLAG.NE.0 ) THEN
         WRITE (6,99002) KW,IFILINP
         STOP
      ELSE
         A = DEF
         IF ( IWRI.NE.0 ) WRITE (IWRI,99001) KW,IFILINP,STR10,A
         REWIND IFILINP
         RETURN
      END IF
C
 300  CONTINUE
      WRITE (6,99003) IFILINP,KW
      STOP
99001 FORMAT (2X,'keyword ',A,' file',I4,' REAL     input: ',A,(20E15.8)
     &        )
99002 FORMAT (//,1X,79('*'),/,2X,'keyword  ',A,
     &        ' not found in file IFILINP',I4,/,2X,
     &        'the STOP flag was set indicating input obligatory',/,
     &        79('*'),/)
99003 FORMAT (//,1X,79('*'),/,2X,'error reading file',I4,
     &        ' for keyword ',A,/,1X,79('*'),/)
      END
C*==is_letter.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      FUNCTION IS_LETTER(C)
C   ********************************************************************
C   *                                                                  *
C   *  checks whether character  C  is letter                          *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--IS_LETTER1882
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*1 C
      LOGICAL IS_LETTER
C
C Local variables
C
      LOGICAL IS_LC_LETTER,IS_UC_LETTER
C
C*** End of declarations rewritten by SPAG
C
      IS_LETTER = IS_UC_LETTER(C) .OR. IS_LC_LETTER(C)
C
      END
C*==is_uc_letter.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      FUNCTION IS_UC_LETTER(C)
C   ********************************************************************
C   *                                                                  *
C   *  checks whether character  C  is Upper Case letter               *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:ICHAR_UCA,ICHAR_UCZ
      IMPLICIT NONE
C*--IS_UC_LETTER1921
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*1 C
      LOGICAL IS_UC_LETTER
C
C Local variables
C
      INTEGER IC
C
C*** End of declarations rewritten by SPAG
C
      IC = ICHAR(C)
C
      IS_UC_LETTER = (ICHAR_UCA.LE.IC) .AND. (IC.LE.ICHAR_UCZ)
C
      END
C*==is_lc_letter.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      FUNCTION IS_LC_LETTER(C)
C   ********************************************************************
C   *                                                                  *
C   *  checks whether character  C  is Lower Case letter               *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:ICHAR_LCA,ICHAR_LCZ
      IMPLICIT NONE
C*--IS_LC_LETTER1962
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*1 C
      LOGICAL IS_LC_LETTER
C
C Local variables
C
      INTEGER IC
C
C*** End of declarations rewritten by SPAG
C
      IC = ICHAR(C)
C
      IS_LC_LETTER = (ICHAR_LCA.LE.IC) .AND. (IC.LE.ICHAR_LCZ)
C
      END
C*==is_number.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      FUNCTION IS_NUMBER(C)
C   ********************************************************************
C   *                                                                  *
C   *  checks whether character  C  is a number                        *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:ICHAR_0,ICHAR_9
      IMPLICIT NONE
C*--IS_NUMBER2003
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*1 C
      LOGICAL IS_NUMBER
C
C Local variables
C
      INTEGER IC
C
C*** End of declarations rewritten by SPAG
C
      IC = ICHAR(C)
C
      IS_NUMBER = (ICHAR_0.LE.IC) .AND. (IC.LE.ICHAR_9)
C
      END
C*==is_name_character.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      FUNCTION IS_NAME_CHARACTER(C)
C   ********************************************************************
C   *                                                                  *
C   *  checks whether character  C  is a character allowed for         *
C   *  a name e.g. of a subroutine:  A-Z, a-z, 0-9 and '_'             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:ICHAR_UNDERSCORE
      IMPLICIT NONE
C*--IS_NAME_CHARACTER2045
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*1 C
      LOGICAL IS_NAME_CHARACTER
C
C Local variables
C
      INTEGER IC
      LOGICAL IS_LETTER,IS_NUMBER
C
C*** End of declarations rewritten by SPAG
C
      IC = ICHAR(C)
C
      IS_NAME_CHARACTER = IS_LETTER(C) .OR. IS_NUMBER(C) .OR. 
     &                    (IC.EQ.ICHAR_UNDERSCORE)
C
      END
C*==symbolkm.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      FUNCTION SYMBOLKM(KAP,MM05)
C   ********************************************************************
C   *                                                                  *
C   *  return term-symbol for quantum numbers     (KAPPA,MJ-1/2)       *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--SYMBOLKM2087
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER KAP,MM05
      CHARACTER*20 SYMBOLKM
C
C Local variables
C
      CHARACTER FUNTXTL
      INTEGER J2,L,LL,M2
      CHARACTER*20 STR
C
C*** End of declarations rewritten by SPAG
C
      L = ABS(KAP) + (SIGN(1,KAP)-1)/2
      J2 = 2*ABS(KAP) - 1
      M2 = 2*MM05 + 1
C
      STR = FUNTXTL(L)//'!s'
      CALL STRING_ADD_N(STR,J2)
      LL = LEN_TRIM(STR)
      IF ( M2.GT.0 ) THEN
         STR = STR(1:LL)//'/2;+'
      ELSE
         STR = STR(1:LL)//'/2;-'
      END IF
      CALL STRING_ADD_N(STR,ABS(M2))
      LL = LEN_TRIM(STR)
      SYMBOLKM = STR(1:LL)//'/2!N'
      END
C*==funtxtl.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      FUNCTION FUNTXTL(L)
C   ********************************************************************
C   *                                                                  *
C   *  return term-symbol for quantum number      L                    *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--FUNTXTL2139
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER L
      CHARACTER*1 FUNTXTL
C
C Local variables
C
      CHARACTER*1 TL(0:2)
C
C*** End of declarations rewritten by SPAG
C
      DATA TL/'s','p','d'/
C
      IF ( L.LE.2 ) THEN
         FUNTXTL = TL(L)
      ELSE
         FUNTXTL = CHAR(ICHAR('f')+L-3)
      END IF
      END
C*==funtxtj.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      FUNCTION FUNTXTJ(J)
C   ********************************************************************
C   *                                                                  *
C   *  return term-symbol for quantum number      J                    *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--FUNTXTJ2182
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 J
      CHARACTER*4 FUNTXTJ
C
C Local variables
C
      INTEGER J2
      CHARACTER*2 STR
C
C*** End of declarations rewritten by SPAG
C
      J2 = NINT(2*J)
C
      WRITE (STR,FMT='(I2)') J2
C
      IF ( J2.GE.10 ) THEN
         FUNTXTJ = STR//'/2'
      ELSE
         FUNTXTJ = STR(2:2)//'/2'
      END IF
C
      END
C*==funtxtmj.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      FUNCTION FUNTXTMJ(MJ)
C   ********************************************************************
C   *                                                                  *
C   *  return term-symbol for quantum number      MJ                   *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--FUNTXTMJ2229
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 MJ
      CHARACTER*5 FUNTXTMJ
C
C Local variables
C
      INTEGER ABSMJ2
      CHARACTER*2 STR
      CHARACTER*1 VZ
C
C*** End of declarations rewritten by SPAG
C
      ABSMJ2 = NINT(2*ABS(MJ))
C
      WRITE (STR,FMT='(I2)') ABSMJ2
C
      IF ( MJ.GT.0D0 ) THEN
         VZ = '+'
      ELSE
         VZ = '-'
      END IF
C
      IF ( ABSMJ2.GE.10 ) THEN
         FUNTXTMJ = VZ//STR(1:2)//'/2'
      ELSE
         FUNTXTMJ = VZ//STR(2:2)//'/2'
      END IF
C     write(6,*) '## FUNTXTMJ ## RESULT ',MJ,  FUNTXTMJ
C
      END
C*==txtmj2num.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE TXTMJ2NUM(TXTMJ,MJ,OK)
C   ********************************************************************
C   *                                                                  *
C   *  convert  MJ  from text to numerical form                        *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--TXTMJ2NUM2284
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 MJ
      LOGICAL OK
      CHARACTER*5 TXTMJ
C
C Local variables
C
      CHARACTER*5 FUNTXTMJ
      INTEGER I,MJM05
      REAL*8 MJTRY
      CHARACTER*5 STR,STRTRY
C
C*** End of declarations rewritten by SPAG
C
      OK = .FALSE.
C
C--------------------------------------------- add sign (+) if necessary
      IF ( TXTMJ(1:1).NE.'+' .AND. TXTMJ(1:1).NE.'-' ) THEN
         STR = '+'//TXTMJ(1:4)
      ELSE
         STR = TXTMJ(1:5)
      END IF
C
C------------------------------------------------- check for ending '/2'
      I = INDEX(STR,'/2')
      IF ( I.EQ.0 ) RETURN
      IF ( I.EQ.3 ) STR = STR(1:4)//' '
C
C---------------------------------------------- try mj = -61/2 ... +61/2
      DO MJM05 = -31,30
         MJTRY = DBLE(MJM05) + 0.5D0
         STRTRY = FUNTXTMJ(MJTRY)
         IF ( STR.EQ.STRTRY ) THEN
            MJ = MJTRY
            OK = .TRUE.
            RETURN
         END IF
      END DO
C
      WRITE (6,*) '############## <TXTMJ2NUM> failed to covert  STR=',
     &            STR
      STOP
      END
C*==getqnmj.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE GETQNMJ(SUB,CMDUC,LCMD,L,MJ,UDT)
C   ********************************************************************
C   *                                                                  *
C   *  get quantum number   mj                                         *
C   *                                                                  *
C   *  if the string MJ=xy/2   with X= blank,+,-  and  y=1,3,5,...     *
C   *  is not found in the input CMDUC, the default MJ=+1/2 is assumed *
C   *  this case is indicated by  UDT = .F.                            *
C   *                                                                  *
C   *  if the input is inconsistent, the program is stopped            *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--GETQNMJ2358
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*(*) CMDUC,SUB
      INTEGER L,LCMD
      REAL*8 MJ
      LOGICAL UDT
C
C Local variables
C
      LOGICAL OK
      CHARACTER*5 TMJ
C
C*** End of declarations rewritten by SPAG
C
C------------------------------------------------------------ scan input
      CALL SECTION_SET_STRING('MJ',TMJ,'9999',0)
C
C------------------------------------------- convert text to real number
      IF ( UDT ) THEN
         CALL TXTMJ2NUM(TMJ,MJ,OK)
         IF ( .NOT.OK ) THEN
            WRITE (6,99001) SUB,CMDUC(1:LCMD)
            STOP 'in <GETQNMJ>'
         END IF
C
      ELSE
C----------------------------------------------------------- set default
         MJ = 0.5D0
C
      END IF
C
C----------------------------------------------------------- final check
      IF ( ABS(MJ).GT.(L+0.5000001D0) ) THEN
         WRITE (6,99002) ' MJ = ',MJ,' too large'
         MJ = SIGN(1D0,MJ)*(L+0.5D0)
         WRITE (6,99002) ' MJ reduced to ',MJ
      END IF
C
      RETURN
99001 FORMAT (/,1X,79('#'),/,10x,'ERROR reading mj quantum number',/,
     &        10X,'invoked by ',A,8X,'example for input: mj=+3/2 ',/,
     &        10X,'input lines read:',/,A)
99002 FORMAT (5x,30('#'),'  WARNING ',A,F5.1,A)
      END
C*==string_replace_tab.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE STRING_REPLACE_TAB(STRING)
C   ********************************************************************
C   *                                                                  *
C   *  replace horizontal TABs  by blanks                              *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--STRING_REPLACE_TAB2425
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*(*) STRING
C
C Local variables
C
      INTEGER I,L
C
C*** End of declarations rewritten by SPAG
C
      L = LEN_TRIM(STRING)
C
      DO I = 1,L
         IF ( IACHAR(STRING(I:I)).EQ.9 ) STRING(I:I) = ' '
      END DO
      END
C*==integer_to_string.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      FUNCTION INTEGER_TO_STRING(I)
C   ********************************************************************
C   *                                                                  *
C   *  write integer I to a string of length 10                        *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--INTEGER_TO_STRING2464
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER I
      CHARACTER*10 INTEGER_TO_STRING
C
C*** End of declarations rewritten by SPAG
C
      WRITE (INTEGER_TO_STRING,'(I10)') I
C
      END
C*==string_replace_blanks.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE STRING_REPLACE_BLANKS(STRING)
C   ********************************************************************
C   *                                                                  *
C   *  replace blanks in   STRING  by "_"                               *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--STRING_REPLACE_BLANKS2494
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*(*) STRING
C
C Local variables
C
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,LEN_TRIM(STRING)
         IF ( STRING(I:I).EQ.' ' ) STRING(I:I) = '_'
      END DO
      END
C*==exchange_characters.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE EXCHANGE_CHARACTERS(C1,C2,STRING)
C   ********************************************************************
C   *                                                                  *
C   *  exchange all characters  C1  by C2 in  STRING(1:LSTR)           *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--EXCHANGE_CHARACTERS2532
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER C1,C2
      CHARACTER*(*) STRING
C
C Local variables
C
      INTEGER I,LSTR
C
C*** End of declarations rewritten by SPAG
C
      LSTR = LEN(STRING)
C
 100  CONTINUE
      I = INDEX(STRING(1:LSTR),C1)
      IF ( I.NE.0 ) THEN
         STRING(I:I) = C2
         GOTO 100
      END IF
C
      END
C*==section_get_core_level_info.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SECTION_GET_CORE_LEVEL_INFO(CL,N,L)
C   ********************************************************************
C   *                                                                  *
C   *  search the string   CX = XX  and set N and L according to XX    *
C   *                                                                  *
C   *  DEFAULT:  CL=2P  implying N=2 L=1                               *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_FILES,ONLY:FOUND_QUANTUM_NUMBER
      IMPLICIT NONE
C*--SECTION_GET_CORE_LEVEL_INFO2581
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SECTION_GET_CORE_LEVEL_INFO')
C
C Dummy arguments
C
      CHARACTER*2 CL
      INTEGER L,N
C
C*** End of declarations rewritten by SPAG
C
      FOUND_QUANTUM_NUMBER = .FALSE.
C
      CALL SECTION_SET_STRING('CL',CL,'2P',0)
C
      N = ICHAR(CL(1:1)) - ICHAR('1') + 1
C
      L = 9999
C
      SELECT CASE (CL(2:2))
C
      CASE ('s','S')
C
         L = 0
C
      CASE ('p','P')
C
         L = 1
C
      CASE ('d','D')
C
         L = 2
C
      CASE ('f','F')
C
         L = 3
C
      CASE DEFAULT
C
         CALL STOP_MESSAGE(ROUTINE,'no proper input value for L')
C
      END SELECT
C
      IF ( N.LT.1 .OR. N.GT.6 )
     &      CALL STOP_MESSAGE(ROUTINE,'n quantum number not found')
      IF ( L.EQ.9999 )
     &      CALL STOP_MESSAGE(ROUTINE,'l quantum number not found')
C
      END
C*==section_get_quantum_number_mj.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SECTION_GET_QUANTUM_NUMBER_MJ(MJ,L)
C   ********************************************************************
C   *                                                                  *
C   *  get quantum number   mj                                         *
C   *                                                                  *
C   *  if the string MJ=xy/2   with X= blank,+,-  and  y=1,3,5,...     *
C   *  is not found in the input CMDUC, the default MJ=+1/2 is assumed *
C   *  this case is indicated by  UDT = .F.                            *
C   *                                                                  *
C   *  if the input is inconsistent, the program is stopped            *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_FILES,ONLY:FOUND_QUANTUM_NUMBER,FOUND_STRING,CMDUC
      IMPLICIT NONE
C*--SECTION_GET_QUANTUM_NUMBER_MJ2663
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SECTION_GET_QUANTUM_NUMBER_MJ')
C
C Dummy arguments
C
      INTEGER L
      REAL*8 MJ
C
C Local variables
C
      LOGICAL OK
      CHARACTER*5 TMJ
C
C*** End of declarations rewritten by SPAG
C
      FOUND_QUANTUM_NUMBER = .FALSE.
C
C------------------------------------------------------------ scan input
      CALL SECTION_SET_STRING('MJ',TMJ,'9999',0)
C
C------------------------------------------- convert text to real number
      IF ( FOUND_STRING ) THEN
C
         CALL TXTMJ2NUM(TMJ,MJ,OK)
C
         IF ( .NOT.OK ) THEN
            WRITE (6,99001) CMDUC(1:LEN_TRIM(CMDUC))
            CALL STOP_MESSAGE(ROUTINE,'quantum number mj not found')
         END IF
C
         FOUND_QUANTUM_NUMBER = .TRUE.
C
      ELSE
C----------------------------------------------------------- set default
         MJ = 0.5D0
C
      END IF
C
C----------------------------------------------------------- final check
      IF ( ABS(MJ).GT.(L+0.5000001D0) ) THEN
         WRITE (6,99002) ' MJ = ',MJ,' too large'
         MJ = SIGN(1D0,MJ)*(L+0.5D0)
         WRITE (6,99002) ' MJ reduced to ',MJ
      END IF
C
      RETURN
C
99001 FORMAT (/,1X,79('#'),/,10x,'ERROR reading mj quantum number',/,
     &        10X,'example for input: mj=+3/2 ',/,10X,
     &        'input lines read:',/,A)
99002 FORMAT (5x,30('#'),'  WARNING ',A,F5.1,A)
      END
