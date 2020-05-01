C*==dmft_fft_ttof2.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C============================================================================
C
      SUBROUTINE DMFT_FFT_TTOF2(FT,DTIME,NN,NOM,NS)
      IMPLICIT NONE
C*--DMFT_FFT_TTOF26
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      DOUBLE PRECISION DTIME
      INTEGER NN,NOM,NS
      COMPLEX*16 FT(NN,NOM,NS)
C
C Local variables
C
      COMPLEX*16 F(:)
      INTEGER IL,IS,K
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE F
C
C*** End of declarations rewritten by SPAG
C
C
C
      ALLOCATE (F(NOM))
      DO IS = 1,NS
         DO IL = 1,NN
            DO K = 1,NOM
               F(K) = FT(IL,K,IS)
            END DO
            CALL DMFT_DFOUR1(F,NOM,-1)
C------ This is standart FFT for DEC-alpha
C     CALL ZFFT('C','C','F',f,f,nom,1)
C
            DO K = 1,NOM
               FT(IL,K,IS) = 0.5D0*DTIME*F(K)
            END DO
C
         END DO
      END DO
C
      END
C*==dmft_fft_ftot2.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C============================================================================
C
      SUBROUTINE DMFT_FFT_FTOT2(FT,DOMEGA,NN,NOM,NS)
      IMPLICIT NONE
C*--DMFT_FFT_FTOT264
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 PI
      PARAMETER (PI=3.141592653589793238462643D0)
C
C Dummy arguments
C
      DOUBLE PRECISION DOMEGA
      INTEGER NN,NOM,NS
      COMPLEX*16 FT(NN,NOM,NS)
C
C Local variables
C
      COMPLEX*16 F(:)
      INTEGER IL,IS,K
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE F
C
C*** End of declarations rewritten by SPAG
C
C
C
      ALLOCATE (F(NOM))
      DO IS = 1,NS
         DO IL = 1,NN
C
            DO K = 1,NOM
               F(K) = FT(IL,K,IS)
            END DO
            CALL DMFT_DFOUR1(F,NOM,1)
C------ This is standart FFT for DEC-alpha
C     CALL ZFFT('C','C','B',f,f,nom,1)
C
            DO K = 1,NOM
               FT(IL,K,IS) = DOMEGA/PI*F(K)
C------ "Back" FFT for DECalpha contains "normalization"!
C        ft(il,k,is)=nom*domega/pi*f(k)
            END DO
C
         END DO
      END DO
C
      END
C*==dmft_fft_ttof4.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
C=============================================================================
C
      SUBROUTINE DMFT_FFT_TTOF4(FT,DTIME,NNNN,NOM,NS)
      IMPLICIT NONE
C*--DMFT_FFT_TTOF4134
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      DOUBLE PRECISION DTIME
      INTEGER NNNN,NOM,NS
      COMPLEX*16 FT(NNNN,NOM,NS)
C
C Local variables
C
      COMPLEX*16 F(:)
      INTEGER IL,IS,K
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE F
C
C*** End of declarations rewritten by SPAG
C
C
C
      ALLOCATE (F(NOM))
      DO IS = 1,NS
         DO IL = 1,NNNN
C
            DO K = 1,NOM
               F(K) = FT(IL,K,IS)
            END DO
            CALL DMFT_DFOUR1(F,NOM,-1)
C------ This is standart FFT for DEC-alpha
C     CALL ZFFT('C','C','F',f,f,nom,1)
C
            DO K = 1,NOM
               FT(IL,K,IS) = 0.5D0*DTIME*F(K)
            END DO
C
         END DO
      END DO
C
      END
C*==dmft_fft_ftot4.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C============================================================================
C
      SUBROUTINE DMFT_FFT_FTOT4(FT,DOMEGA,NNNN,NOM,NS)
      IMPLICIT NONE
C*--DMFT_FFT_FTOT4193
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 PI
      PARAMETER (PI=3.141592653589793238462643D0)
C
C Dummy arguments
C
      DOUBLE PRECISION DOMEGA
      INTEGER NNNN,NOM,NS
      COMPLEX*16 FT(NNNN,NOM,NS)
C
C Local variables
C
      COMPLEX*16 F(:)
      INTEGER IL,IS,K
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE F
C
C*** End of declarations rewritten by SPAG
C
C
C
      ALLOCATE (F(NOM))
      DO IS = 1,NS
         DO IL = 1,NNNN
C
            DO K = 1,NOM
               F(K) = FT(IL,K,IS)
            END DO
            CALL DMFT_DFOUR1(F,NOM,1)
C------ This is standart FFT for DEC-alpha
C     CALL ZFFT('C','C','B',f,f,nom,1)
C
            DO K = 1,NOM
               FT(IL,K,IS) = DOMEGA/PI*F(K)
C------ "Back" FFT for DECalpha contains "normalization"!
C        ft(il,k,is)=nom*domega/pi*f(k)
            END DO
C
         END DO
      END DO
C
      END
C*==dmft_fft_ttof4m.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
C=============================================================================
C
      SUBROUTINE DMFT_FFT_TTOF4M(FT,DTIME,NNNN,NOM)
      IMPLICIT NONE
C*--DMFT_FFT_TTOF4M263
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      DOUBLE PRECISION DTIME
      INTEGER NNNN,NOM
      COMPLEX*16 FT(NNNN,NOM)
C
C Local variables
C
      COMPLEX*16 F(:)
      INTEGER IL,K
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE F
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATE (F(NOM))
      DO IL = 1,NNNN
C
         DO K = 1,NOM
            F(K) = FT(IL,K)
         END DO
         CALL DMFT_DFOUR1(F,NOM,-1)
C------ This is standart FFT for DEC-alpha
C     CALL ZFFT('C','C','F',f,f,nom,1)
C
         DO K = 1,NOM
            FT(IL,K) = 0.5D0*DTIME*F(K)
         END DO
C
      END DO
C
      END
C*==dmft_fft_ftot4m.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C============================================================================
C
      SUBROUTINE DMFT_FFT_FTOT4M(FT,DOMEGA,NNNN,NOM)
      IMPLICIT NONE
C*--DMFT_FFT_FTOT4M319
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 PI
      PARAMETER (PI=3.141592653589793238462643D0)
C
C Dummy arguments
C
      DOUBLE PRECISION DOMEGA
      INTEGER NNNN,NOM
      COMPLEX*16 FT(NNNN,NOM)
C
C Local variables
C
      COMPLEX*16 F(:)
      INTEGER IL,K
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE F
C
C*** End of declarations rewritten by SPAG
C
C
C
C
      ALLOCATE (F(NOM))
      DO IL = 1,NNNN
C
         DO K = 1,NOM
            F(K) = FT(IL,K)
         END DO
         CALL DMFT_DFOUR1(F,NOM,1)
C------ This is standart FFT for DEC-alpha
C     CALL ZFFT('C','C','B',f,f,nom,1)
C
         DO K = 1,NOM
            FT(IL,K) = DOMEGA/PI*F(K)
C------ "Back" FFT for DECalpha contains "normalization"!
C        ft(il,k)=nom*domega/pi*f(k)
         END DO
C
      END DO
C
      END
C*==dmft_fft_ftot4s.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
C=============================================================================
C
      SUBROUTINE DMFT_FFT_FTOT4S(FT,DOMEGA,NNNN,NOM,NS)
      IMPLICIT NONE
C*--DMFT_FFT_FTOT4S388
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 PI
      PARAMETER (PI=3.141592653589793238462643D0)
C
C Dummy arguments
C
      DOUBLE PRECISION DOMEGA
      INTEGER NNNN,NOM,NS
      COMPLEX*16 FT(NNNN,NOM,NS+1)
C
C Local variables
C
      COMPLEX*16 F(:)
      INTEGER IL,IS,K
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE F
C
C*** End of declarations rewritten by SPAG
C
C
C
C
      ALLOCATE (F(NOM))
      DO IS = 1,NS + 1
         DO IL = 1,NNNN
C
            DO K = 1,NOM
               F(K) = FT(IL,K,IS)
            END DO
            CALL DMFT_DFOUR1(F,NOM,1)
C------ This is standart FFT for DEC-alpha
C     CALL ZFFT('C','C','B',f,f,nom,1)
C
            DO K = 1,NOM
               FT(IL,K,IS) = DOMEGA/PI*F(K)
C------ "Back" FFT for DECalpha contains "normalization"!
C        ft(il,k,is)=nom*domega/pi*f(k)
            END DO
C
         END DO
      END DO
C
      END
C*==dmft_dfour1.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
C=============================================================================
      SUBROUTINE DMFT_DFOUR1(CDATAM,NN,ISIG)
      IMPLICIT NONE
C*--DMFT_DFOUR1458
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ISIG,NN
      COMPLEX*16 CDATAM(NN)
C
C Local variables
C
      DOUBLE PRECISION DATAM(:),TEMPI,TEMPR,THETA,WI,WPI,WPR,WR,WTEMP
      INTEGER I,ISTEP,J,M,MMAX,N
C
C*** End of declarations rewritten by SPAG
C
C
C
C Local variables
C
      ALLOCATABLE DATAM
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATE (DATAM(2*NN))
      N = 2*NN
      DO I = 1,N - 1,2
         DATAM(I) = DREAL(CDATAM(I))
         DATAM(I+1) = DIMAG(CDATAM(I))
      END DO
C
      J = 1
      DO I = 1,N,2
         IF ( J.GT.I ) THEN
            TEMPR = DATAM(J)
            TEMPI = DATAM(J+1)
            DATAM(J) = DATAM(I)
            DATAM(J+1) = DATAM(I+1)
            DATAM(I) = TEMPR
            DATAM(I+1) = TEMPI
         END IF
         M = N/2
 50      CONTINUE
         IF ( (M.GE.2) .AND. (J.GT.M) ) THEN
            J = J - M
            M = M/2
            GOTO 50
         END IF
         J = J + M
      END DO
      MMAX = 2
 100  CONTINUE
      IF ( N.GT.MMAX ) THEN
         ISTEP = 2*MMAX
         THETA = 6.28318530717959D0/(ISIG*MMAX)
         WPR = -2.D0*SIN(0.5D0*THETA)**2
         WPI = SIN(THETA)
         WR = 1.D0
         WI = 0.D0
         DO M = 1,MMAX,2
            DO I = M,N,ISTEP
               J = I + MMAX
               TEMPR = WR*DATAM(J) - WI*DATAM(J+1)
               TEMPI = WR*DATAM(J+1) + WI*DATAM(J)
               DATAM(J) = DATAM(I) - TEMPR
               DATAM(J+1) = DATAM(I+1) - TEMPI
               DATAM(I) = DATAM(I) + TEMPR
               DATAM(I+1) = DATAM(I+1) + TEMPI
            END DO
            WTEMP = WR
            WR = WR*WPR - WI*WPI + WR
            WI = WI*WPR + WTEMP*WPI + WI
         END DO
         MMAX = ISTEP
         GOTO 100
      END IF
      DO I = 1,N - 1,2
         CDATAM(I) = DCMPLX(DATAM(I),DATAM(I+1))
      END DO
      END
C======================================================================
C
