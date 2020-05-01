C*==wrhead.f    processed by SPAG 6.70Rc at 15:55 on 19 Dec 2016
      SUBROUTINE WRHEAD(IFIL,FILNAM,KW,NE)
C   ********************************************************************
C   *                                                                  *
C   *   open and write  head for data-file    IFIL,FILNAM              *
C   *                                                                  *
C   ********************************************************************
      USE MOD_KSPACE,ONLY:NKTAB
      USE MOD_ENERGY,ONLY:EFERMI
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_FILES,ONLY:INFO,LINFO,SYSTEM,LSYSTEM,TITLE
      USE MOD_SITES,ONLY:IQAT,IQBOT,IQTOP
      USE MOD_TYPES,ONLY:TXT_T,CONC,NAT,ITBOT,ITTOP
      USE MOD_ANGMOM,ONLY:NLQ
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*(*) FILNAM,KW
      INTEGER IFIL,NE
C
C Local variables
C
      REAL*8 EFERMIOUT
      INTEGER IA,IQ,IT
C
C*** End of declarations rewritten by SPAG
C
      OPEN (IFIL,FILE=FILNAM(1:LEN_TRIM(FILNAM)))
C
      WRITE (IFIL,99001) 'KEYWORD   ',KW(1:LEN_TRIM(KW))
      WRITE (IFIL,99001) 'TITLE     ',TITLE(1:LEN_TRIM(TITLE))
      WRITE (IFIL,99001) 'SYSTEM    ',SYSTEM(1:LSYSTEM)
      WRITE (IFIL,99003) 'NQ_eff    ',(IQTOP-IQBOT+1)
      WRITE (IFIL,99003) 'NT_eff    ',(ITTOP-ITBOT+1)
      WRITE (IFIL,99003) '          '
      WRITE (IFIL,99003) 'NE        ',NE
      WRITE (IFIL,99003) 'IREL      ',IREL
      IF ( EFERMI.LT.-9D0 ) THEN
         EFERMIOUT = 0D0
      ELSE
         EFERMIOUT = EFERMI
      END IF
      WRITE (IFIL,99004) 'EFERMI    ',EFERMIOUT
      WRITE (IFIL,99002) 'INFO      ',INFO(1:LINFO),NKTAB
      WRITE (IFIL,99002) '          '
C
      WRITE (IFIL,'(''   IQ  NLQ '')')
      DO IQ = IQBOT,IQTOP
         WRITE (IFIL,FMT='(2I5)') IQ,NLQ(IQ)
      END DO
C
      WRITE (IFIL,'(''   IT       TXT_T      CONC  NAT IQAT'')')
C
      DO IT = ITBOT,ITTOP
         WRITE (IFIL,99005) IT,TXT_T(IT),CONC(IT),NAT(IT),
     &                      ((IQAT(IA,IT)-IQBOT+1),IA=1,NAT(IT))
      END DO
99001 FORMAT (A10,A)
99002 FORMAT (A10,A,I5)
99003 FORMAT (A10,I10)
99004 FORMAT (A10,F10.5)
99005 FORMAT (I5,1X,A,6X,F10.5,I5,10I5,:,/,(41X,10I5))
      END
C*==rdhead.f    processed by SPAG 6.70Rc at 15:55 on 19 Dec 2016
      SUBROUTINE RDHEAD(IFIL)
C **********************************************************************
C *                                                                    *
C *       read the standard head of the input data file                *
C *                                                                    *
C **********************************************************************
C
      USE MOD_FILES,ONLY:SYSTEM,LSYSTEM
      USE MOD_TYPES,ONLY:NT
      USE MOD_SITES,ONLY:NQ
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IFIL
C
C Local variables
C
      INTEGER IQ,IT
C
C*** End of declarations rewritten by SPAG
C
      READ (IFIL,*)
      READ (IFIL,*)
      READ (IFIL,*)
      READ (IFIL,*)
      READ (IFIL,*)
      READ (IFIL,*)
      READ (IFIL,*)
      READ (IFIL,*)
      READ (IFIL,*)
      READ (IFIL,*)
C
C ----------------------------------------------------------------------
      READ (IFIL,*)
C
      DO IQ = 1,NQ
         READ (IFIL,*)
      END DO
C ----------------------------------------------------------------------
      READ (IFIL,*)
      DO IT = 1,NT
         READ (IFIL,*)
      END DO
C
      WRITE (6,*) ' SYSTEM:     ',SYSTEM(1:LSYSTEM)
      WRITE (6,*) '<RDHEAD>  passed '
C
      END
