C*==runssite.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFILSS,GETIRRSOL,
     &                    ERYD,P,IPRINT,TSST,MSST,SSST,MEZZ,MEZJ,
     &                    ORBPOLSS)
C   ********************************************************************
C   *                                                                  *
C   *   driver to run the various single site routines                 *
C   *   according to the calculation mode                              *
C   *                                                                  *
C   *   CHECK:  set  NPAN = 1, JRWS=JRMT=JRCUR(1),  VNST=BNST=0        *
C   *           run   <FPSSITE>  AND  <SSITE>                          *
C   *           or  <FPNRSSITE>  AND  <NRSSITE>  according to IREL     *
C   *           this restricts the potential regime to the mt-regime   *
C   *           as a consquence all results should be the same         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_RMESH,ONLY:JRNSMIN,NMMAX,NRMAX,NPAN,JRCUT,JRCRI,FULLPOT,
     &    JRWS
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX,NKM
      USE MOD_TYPES,ONLY:NTMAX,NLMFPMAX,BNST,VNST
      IMPLICIT NONE
C*--RUNSSITE24
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      LOGICAL CHECK
      PARAMETER (CHECK=.FALSE.)
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='RUNSSITE')
C
C Dummy arguments
C
      LOGICAL CALCINT,GETIRRSOL
      COMPLEX*16 ERYD,P
      INTEGER IFILSS,IPRINT,IWRIRRWF,IWRREGWF
      CHARACTER*10 ORBPOLSS
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),SSST(NKMMAX,NKMMAX,NTMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      INTEGER IM
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
      IF ( IREL.LE.2 ) THEN
C
C=======================================================================
C             NON- or SCALAR-relativistic calculation
C=======================================================================
C
         IF ( FULLPOT ) THEN
C
            IF ( CHECK ) THEN
               IPRINT = 5
               DO IM = 1,NMMAX
                  NPAN(IM) = 1
                  JRCRI(IM) = JRCUT(1,IM)
                  JRWS(IM) = JRCUT(1,IM)
               END DO
               CALL RINIT((NRMAX-JRNSMIN+1)*NLMFPMAX*NTMAX,VNST)
               CALL RINIT((NRMAX-JRNSMIN+1)*NLMFPMAX*NTMAX,BNST)
               OPEN (6,FILE='zzz_FPNRSSITE.out')
            END IF
C
            CALL FPNRSSITE(IWRREGWF,IWRIRRWF,IFILSS,GETIRRSOL,ERYD,P,
     &                     IPRINT,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOLSS)
C
            IF ( CHECK ) THEN
               CLOSE (6)
               OPEN (6,FILE='zzz_NRSSITE.out')
C
               CALL NRSSITE(IWRREGWF,IWRIRRWF,IFILSS,ERYD,P,IPRINT,TSST,
     &                      MSST,SSST,MEZZ,MEZJ)
               CLOSE (6)
               WRITE (6,*) 'output redirected to zzz_...out - compare:'
               WRITE (6,*) 'tkdiff zzz_FPNRSSITE.out zzz_NRSSITE.out'
               CALL STOP_REGULAR(ROUTINE,'check completed')
            END IF
C
         ELSE
C
            CALL NRSSITE(IWRREGWF,IWRIRRWF,IFILSS,ERYD,P,IPRINT,TSST,
     &                   MSST,SSST,MEZZ,MEZJ)
         END IF
C
C=======================================================================
C                 FULLY - relativistic calculation
C=======================================================================
C
      ELSE IF ( FULLPOT ) THEN
C
         IF ( CHECK ) THEN
            IPRINT = 5
            DO IM = 1,NMMAX
               NPAN(IM) = 1
               JRCRI(IM) = JRCUT(1,IM)
               JRWS(IM) = JRCUT(1,IM)
            END DO
            CALL RINIT((NRMAX-JRNSMIN+1)*NLMFPMAX*NTMAX,VNST)
            CALL RINIT((NRMAX-JRNSMIN+1)*NLMFPMAX*NTMAX,BNST)
            OPEN (6,FILE='zzz_FPNRSSITE.out')
         END IF
C
         CALL FPSSITE(IWRREGWF,IWRIRRWF,IFILSS,GETIRRSOL,ERYD,P,IPRINT,
     &                TSST,MSST,SSST,MEZZ,MEZJ,ORBPOLSS)
C
C
         IF ( CHECK ) THEN
            CLOSE (6)
            OPEN (6,FILE='zzz_NRSSITE.out')
C
            CALL SSITE(IWRREGWF,IWRIRRWF,IFILSS,CALCINT,GETIRRSOL,ERYD,
     &                 P,IPRINT,NKM,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOLSS)
C
            CLOSE (6)
            WRITE (6,*) 'output redirected to zzz_...out - compare:'
            WRITE (6,*) 'tkdiff zzz_FPNRSSITE.out zzz_NRSSITE.out'
            CALL STOP_REGULAR(ROUTINE,'check completed')
         END IF
C
      ELSE
C
         CALL SSITE(IWRREGWF,IWRIRRWF,IFILSS,CALCINT,GETIRRSOL,ERYD,P,
     &              IPRINT,NKM,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOLSS)
C
      END IF
C
      END
