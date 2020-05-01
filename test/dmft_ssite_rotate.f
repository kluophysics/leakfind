C*==dmft_ssite_rotate.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_SSITE_ROTATE(TSST)
C   ********************************************************************
C   *                                                                  *
C   *  rotate type dependent ssite-matrices from the                   *
C   *  GLOBAL to the LOCAL frame                                       *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:IREL,KMROT
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,WKM1,NKMQ
      USE MOD_SITES,ONLY:IQAT,DROTQ,MAGROT_Q
      USE MOD_TYPES,ONLY:ITBOT,ITTOP,NTMAX
      IMPLICIT NONE
C*--DMFT_SSITE_ROTATE15
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='DMFT_SSITE_ROTATE')
C
C Dummy arguments
C
      COMPLEX*16 TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      INTEGER IQ,IT,M,N
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
      IF ( KMROT.EQ.0 ) RETURN
C
      M = NKMMAX
C
      LOOP_IT:DO IT = ITBOT,ITTOP
C
         IQ = IQAT(1,IT)
C
         IF ( .NOT.MAGROT_Q(IQ) ) CYCLE LOOP_IT
C
         IF ( IREL.NE.2 ) THEN
            N = NKMQ(IQ)
         ELSE
            N = NKM
         END IF
C
C------------------------------------------------------------------- TSST
         CALL ROTATE(TSST(1,1,IT),'G->L',WKM1,N,DROTQ(1,1,IQ),M)
         TSST(1:N,1:N,IT) = WKM1(1:N,1:N)
C
      END DO LOOP_IT
C
      END
