C*==dmft_sptflex.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_SPTFLEX(UUU,UJJ,TEMP,GOM,SOM,MINOM,OMEGA,DTIME,
     &                        DOMEGA,NSPIN,NLM,LOP,NLMSQ,NLMQU,NOM,
     &                        DMFTDBLC,IPRINT,IFLEX)
C
C***********************************************************
C********** SP-T-FLEX  for  N-band spin-polarised  *********
C********** SP-Tmatrix: A.Lichtenstein (Nijmegen)  *********
C********** M.Katsnelson (IFM)                     *********
C********** J. Minar + L.Chioncel KKR implementation *******
C********** 14.09.2005: JM rewritten               *********
C***********************************************************
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--DMFT_SPTFLEX15
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*4 DMFTDBLC
      REAL*8 DOMEGA,DTIME,TEMP,UJJ,UUU
      INTEGER IFLEX,IPRINT,LOP,NLM,NLMQU,NLMSQ,NOM,NSPIN
      COMPLEX*16 GOM(NLM,NLM,NOM,NSPIN),SOM(NLM,NLM,NOM,NSPIN)
      INTEGER MINOM(NOM)
      REAL*8 OMEGA(NOM)
C
C Local variables
C
      INTEGER IOM,IS,LM,LMP
      REAL*8 TIMEE,TIMES
      COMPLEX*16 UC(:,:),UD(:,:),UM(:,:),US(:,:),UT(:,:)
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE UC,UD,UM,US,UT
C============================================================
C     Bosonic GF set to zero
C============================================================
      DO IS = 1,NSPIN
         DO IOM = 1,NOM
            IF ( MODULO(IOM,2).NE.0 ) THEN
               DO LMP = 1,NLM
                  DO LM = 1,NLM
                     GOM(LM,LMP,IOM,IS) = C0
                  END DO
               END DO
            END IF
         END DO
      END DO
      SOM(1:NLM,1:NLM,1:NOM,1:NSPIN) = C0
C============================================================
C     Setup 4-index U matrix
C============================================================
      ALLOCATE (UC(NLMSQ,NLMSQ),UD(NLMSQ,NLMSQ),UM(NLMSQ,NLMSQ))
      ALLOCATE (US(NLMSQ,NLMSQ),UT(NLMSQ,NLMSQ))
      UC(1:NLMSQ,1:NLMSQ) = C0
      UD(1:NLMSQ,1:NLMSQ) = C0
      UM(1:NLMSQ,1:NLMSQ) = C0
      US(1:NLMSQ,1:NLMSQ) = C0
      UT(1:NLMSQ,1:NLMSQ) = C0
      IF ( IPRINT.GT.0 ) WRITE (6,*) 
     &                            'Calculation of the 4-index U matrix:'
C
      CALL DMFT_VERTEX(UUU,UJJ,UC,UD,UM,US,UT,NLM,LOP,IPRINT)
C
      IF ( IPRINT.GT.0 ) THEN
         WRITE (6,*) 'Ueff (Ry)=',UUU,'Jeff (Ry)=',UJJ
         WRITE (6,*) 'Uij ='
         WRITE (6,'(5f9.4)') (DREAL(UC(LM,LM)),LM=1,NLMSQ)
      END IF
C============================================================
C     Impurity solver
C============================================================
      CALL CPU_TIME(TIMES)
C
      CALL DMFT_SIG_SPTFLEX(GOM,SOM,TEMP,MINOM,OMEGA,DTIME,DOMEGA,UC,
     &                      DMFTDBLC,UUU,UJJ,NLM,NLMSQ,NLMQU,NOM,NSPIN,
     &                      IPRINT,IFLEX)
C
      CALL CPU_TIME(TIMEE)
      IF ( IPRINT.GT.0 ) WRITE (6,99001) TIMEE - TIMES
      RETURN
99001 FORMAT (/,5X,'execution time for DMFT solver   ',F14.3,' secs',/)
      END
