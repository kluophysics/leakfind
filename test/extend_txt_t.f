C*==extend_txt_t.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE EXTEND_TXT_T
C **********************************************************************
C *                                                                    *
C *  add an extension to  TXT_T  if there are more types with same Z   *
C *                                                                    *
C **********************************************************************
      USE MOD_TYPES,ONLY:NT,TXT_T,LTXT_T,Z
      IMPLICIT NONE
C*--EXTEND_TXT_T10
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      INTEGER I0,I1,I2,IC2,ICAL,ICAU,ICZL,ICZU,IEXT,IT,JT,KDONE(NT)
C
C*** End of declarations rewritten by SPAG
C
      ICAL = ICHAR('a')
      ICAU = ICHAR('A')
      ICZL = ICHAR('z')
      ICZU = ICHAR('Z')
C
C----------------------------------- reset type names to chemical symbol
      DO IT = 1,NT
         TXT_T(IT) = TXT_T(IT)(1:2)//'  '
         IC2 = ICHAR(TXT_T(IT)(2:2))
C
         IF ( ((IC2.LT.ICAL) .OR. (IC2.GT.ICZL)) .AND. 
     &        ((IC2.LT.ICAU) .OR. (IC2.GT.ICZU)) ) THEN
            TXT_T(IT)(2:2) = ' '
            LTXT_T(IT) = 1
         ELSE
            LTXT_T(IT) = 2
         END IF
         KDONE(IT) = 0
      END DO
C
      DO IT = 1,NT
         I0 = LTXT_T(IT)
         I1 = LTXT_T(IT) + 1
         I2 = I1 + 1
         IEXT = 1
         IF ( KDONE(IT).NE.1 ) THEN
            DO JT = 1,NT
               IF ( (IT.NE.JT) .AND. 
     &              (TXT_T(IT)(1:I0).EQ.TXT_T(JT)(1:I0)) .AND. 
     &              (Z(IT).EQ.Z(JT)) ) THEN
                  TXT_T(IT)(I1:I2) = '_1'
                  IEXT = IEXT + 1
                  TXT_T(JT)(1:I1) = TXT_T(IT)(1:I1)
                  CALL STRING_ADD_N(TXT_T(JT),IEXT)
                  LTXT_T(JT) = LEN_TRIM(TXT_T(JT))
                  KDONE(JT) = 1
               END IF
            END DO
         END IF
C
         IF ( IEXT.NE.1 ) LTXT_T(IT) = LTXT_T(IT) + 2
         KDONE(IT) = 1
C
      END DO
C
      END
