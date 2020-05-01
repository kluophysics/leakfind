C*==lmatstruct.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LMATSTRUCT(STR,A,N,M,MLIN,MCOL,NQ,K_OUTPUT)
C   ********************************************************************
C   *                                                                  *
C   *   writes structure of   LOGICAL NxN   matrix   A                 *
C   *                                                                  *
C   *   M           is the actual array - size used for   A            *
C   *   MLIN/COL    MODE for line and column indexing                  *
C   *               0: plain, 1: (l,ml), 2: (l,ml,ms), 3: (kap,mue)    *
C   *   TOL         tolerance for difference                           *
C   *   NQ          if NQ > 1  all IQ-JQ-blocks are treated one by one *
C   *   K_OUTPUT    output contral                                     *
C   *               -6, +6:  output to standard output                 *
C   *               else:    the matrix will be printed to a file      *
C   *                        with name temp_STR.tmp (blanks in STR are *
C   *                        replaced by underscores)                  *
C   *                        if the file exists the new matrices will  *
C   *                        be appended                               *
C   *                                                                  *
C   *   any changes should be done in all  *MATSTRUCT  routines !!!    *
C   *                                                                  *
C   *   LMFT >= 6 + 2*(l_max+1)**2 + 2*4*(l_max+1) + 7                 *
C   *   setting works up to l_max = 8                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IOTMP4
      IMPLICIT NONE
C*--LMATSTRUCT29
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='LMATSTRUCT')
C
C Dummy arguments
C
      INTEGER K_OUTPUT,M,MCOL,MLIN,N,NQ
      CHARACTER*(*) STR
      LOGICAL A(M,M)
C
C Local variables
C
      LOGICAL B(:,:),DIAGONAL,F_AVAILABLE,SYMMETRIC
      CHARACTER*3 BTXT(0:3)
      CHARACTER*1 CTAB(:)
      CHARACTER*250 FILENAME,FMT1,FMT2,FMT3,FMT4
      INTEGER I,I1,IA_ERR,IC0,IFIL,IL,ILSEP(20),IQ,ISL,IW(:),J,J0,JP,JQ,
     &        K,L3,LF,MM,N1,N2,N3,NK,NM,NM1,NM2,NM3,NNON0,NQEFF,NSL
C
C*** End of declarations rewritten by SPAG
C
      DATA BTXT/'   ','LM ','LMS','KM '/
C
      ALLOCATABLE B,CTAB,IW
C
      ALLOCATE (B(N,N),CTAB(0:N*N),IW(M),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'B')
C
      DIAGONAL = .TRUE.
      SYMMETRIC = .TRUE.
C
      FILENAME = ' '
      FILENAME(1:LEN_TRIM(STR)) = STR(1:LEN_TRIM(STR))
      CALL STRING_REPLACE_BLANKS(FILENAME)
      FILENAME = 'temp_'//FILENAME(1:LEN_TRIM(FILENAME))//'.tmp'
C
      IF ( ABS(K_OUTPUT).EQ.6 ) THEN
         IFIL = 6
      ELSE
         WRITE (6,*) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &               ' :: writing matrix to file: ',STR(1:LEN_TRIM(STR))
         IFIL = IOTMP4
         I = LEN_TRIM(FILENAME)
         INQUIRE (FILE=FILENAME(1:I),EXIST=F_AVAILABLE)
         IF ( F_AVAILABLE ) THEN
            OPEN (UNIT=IFIL,FILE=FILENAME(1:I),POSITION='APPEND')
         ELSE
            OPEN (UNIT=IFIL,FILE=FILENAME(1:I),STATUS='NEW')
         END IF
      END IF
C
C---------------------------------------------------------------- header
      IC0 = ICHAR('0')
      N3 = N/100
      N2 = N/10 - N3*10
      N1 = N - N2*10 - N3*100
C
      IF ( N.LE.18 ) THEN
         FMT1 = '(8X,I3,''|'','
         FMT2 = '( 9X,''--|'','
         FMT3 = '( 9X,'' #|'','
         FMT4 = '( 9X,''  |'','
      ELSE
         FMT1 = '(   I4,''|'','
         FMT2 = '( 2X,''--|'','
         FMT3 = '( 2X,'' #|'','
         FMT4 = '( 2X,''  |'','
      END IF
C
      LF = 11
      L3 = 11
      IF ( MCOL.EQ.0 ) THEN
         FMT1 = FMT1(1:LF)//CHAR(IC0+N3)//CHAR(IC0+N2)//CHAR(IC0+N1)
     &          //'( 2A1),''|'',I3)'
         FMT2 = FMT2(1:LF)//CHAR(IC0+N3)//CHAR(IC0+N2)//CHAR(IC0+N1)
     &          //'(''--''),''|'',I3)'
         FMT3 = FMT3(1:LF)//'60(2X,I2))'
         FMT4 = FMT4(1:LF)//'60(I2,2X))'
         LF = 21
      ELSE
         IF ( MCOL.EQ.1 ) THEN
            NK = NINT(SQRT(DBLE(N)))
         ELSE IF ( MCOL.EQ.2 ) THEN
            NK = NINT(SQRT(DBLE(N/2)))
         ELSE IF ( MCOL.EQ.3 ) THEN
            NK = 2*NINT(SQRT(DBLE(N/2))) - 1
         END IF
         DO K = 1,NK
            IF ( MCOL.LE.2 ) THEN
               NM = 2*K - 1
            ELSE
               NM = 2*((K+1)/2)
            END IF
            NM2 = NM/10
            NM1 = NM - NM2*10
            NM3 = NM/2
            FMT1 = FMT1(1:LF)//CHAR(IC0+NM2)//CHAR(IC0+NM1)
     &             //'( 2A1),''|'','
            FMT2 = FMT2(1:LF)//CHAR(IC0+NM2)//CHAR(IC0+NM1)
     &             //'(''--''),''|'','
C
            IF ( MCOL.LE.2 ) THEN
               DO MM = 1,NM
                  IF ( MOD(MM,2).EQ.MOD(K,2) ) THEN
                     FMT3 = FMT3(1:L3)//'2X,'
                     FMT4 = FMT4(1:L3)//'I2,'
                  ELSE
                     FMT3 = FMT3(1:L3)//'I2,'
                     FMT4 = FMT4(1:L3)//'2X,'
                  END IF
                  L3 = L3 + 3
               END DO
               FMT3 = FMT3(1:L3)//'''|'','
               FMT4 = FMT4(1:L3)//'''|'','
               L3 = L3 + 4
            ELSE
               FMT3 = FMT3(1:LF)//CHAR(IC0+NM3)//'(2X,I2),''|'','
               FMT4 = FMT4(1:LF)//CHAR(IC0+NM3)//'(I2,2X),''|'','
               L3 = L3 + 13
            END IF
            LF = LF + 13
         END DO
         IF ( MCOL.EQ.2 ) THEN
            FMT1 = FMT1(1:LF)//FMT1(12:LF)
            FMT2 = FMT2(1:LF)//FMT2(12:LF)
C
            FMT3 = FMT3(1:L3)//FMT4(12:L3)
            FMT4 = FMT4(1:L3)//FMT3(12:L3)
            LF = 2*LF - 11
            L3 = 2*L3 - 11
         END IF
         FMT1 = FMT1(1:LF)//'I3)'
         FMT2 = FMT2(1:LF)//'I3)'
         FMT3 = FMT3(1:L3)//'I3)'
         FMT4 = FMT4(1:L3)//'I3)'
      END IF
      IF ( MLIN.EQ.0 ) THEN
         NSL = 1
         ILSEP(1) = N
      ELSE IF ( MLIN.EQ.1 ) THEN
         NSL = NINT(SQRT(DBLE(N)))
         DO IL = 1,NSL
            ILSEP(IL) = IL**2
         END DO
      ELSE IF ( MLIN.EQ.2 ) THEN
         NSL = NINT(SQRT(DBLE(N/2)))
         DO IL = 1,NSL
            ILSEP(IL) = IL**2
         END DO
         DO IL = 1,NSL
            ILSEP(NSL+IL) = ILSEP(NSL) + IL**2
         END DO
         NSL = 2*NSL
      ELSE IF ( MLIN.EQ.3 ) THEN
         NSL = 2*NINT(SQRT(DBLE(N/2))) - 1
         ILSEP(1) = 2
         DO K = 2,NSL
            ILSEP(K) = ILSEP(K-1) + 2*((K+1)/2)
         END DO
      END IF
C
      IF ( L3.GT.250 ) CALL STOP_MESSAGE(ROUTINE,'L3 > 250')
C
      WRITE (IFIL,99001) STR(1:(LEN_TRIM(STR)))
C
      NQEFF = MAX(1,NQ)
      IF ( N*NQEFF.GT.M ) CALL STOP_MESSAGE(ROUTINE,'N*NQ > M')
C
C==================================================================== IQ
      DO IQ = 1,NQEFF
C==================================================================== JQ
         DO JQ = 1,NQEFF
C
C----------------------------------------------------- copy matrix block
C
            J0 = N*(JQ-1)
            DO J = 1,N
               I1 = N*(IQ-1) + 1
               JP = J0 + J
               B(1:N,J) = A(I1:(I1+N-1),JP)
            END DO
C
            IF ( NQ.GT.1 ) WRITE (IFIL,99002) BTXT(MLIN),BTXT(MCOL),IQ,
     &                            JQ
C
            WRITE (IFIL,FMT3) (I,I=2,N,2)
            WRITE (IFIL,FMT4) (I,I=1,N,2)
            WRITE (IFIL,FMT=FMT2)
C------------------------------------------------------------ header end
            NNON0 = 0
            CTAB(0) = ' '
            CTAB(1) = 'T'
C
            DO I = 1,N
               DO J = 1,N
                  IF ( B(I,J) ) THEN
                     IF ( I.NE.J ) DIAGONAL = .FALSE.
                     IF ( .NOT.B(J,I) ) SYMMETRIC = .FALSE.
                     NNON0 = NNON0 + 1
                     IW(J) = 1
                  ELSE
                     IW(J) = 0
                  END IF
               END DO
C------------------------------------------------------------ write line
               WRITE (IFIL,FMT=FMT1) I,(' ',CTAB(IW(J)),J=1,N),I
C
               DO ISL = 1,NSL
                  IF ( I.EQ.ILSEP(ISL) ) WRITE (IFIL,FMT=FMT2)
               END DO
            END DO
C
C------------------------------------------------------------------ foot
C
            WRITE (IFIL,FMT4) (I,I=1,N,2)
            WRITE (IFIL,FMT3) (I,I=2,N,2)
C
            IF ( K_OUTPUT.GT.0 ) THEN
               WRITE (IFIL,99002) NNON0
C
               IF ( .NOT.(DIAGONAL) ) THEN
C
                  IF ( SYMMETRIC ) WRITE (IFIL,99003) 'SYMMETRIC '
C
               ELSE
C
                  WRITE (IFIL,99003) 'DIAGONAL  '
C
               END IF
            ELSE
               WRITE (IFIL,*) ' '
            END IF
C
         END DO
C==================================================================== JQ
      END DO
C==================================================================== IQ
C
      IF ( IFIL.EQ.IOTMP4 ) CLOSE (IOTMP4)
      DEALLOCATE (B,CTAB,IW,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC')
C
99001 FORMAT (/,8X,A,/)
99002 FORMAT (/,8X,I5,' elements set',/)
99003 FORMAT (8X,'the matrix is  ',A)
      END
