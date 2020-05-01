C*==tau_drive.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TAU_DRIVE(IECURR,IPRINTBAND,ERYD,P,TSSQ,MSSQ,TSST,MSST,
     &                     TAUQ,ICPAFLAG,CPACHNG,ITCPA,ICPACONV,PHASK)
C   ********************************************************************
C   *                                                                  *
C   *    driver routine to call the appropriate routine to solve       *
C   *    the multiple scattering problem according to   KKRMODE        *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:NQ,NQMAX,ITOQ,NOQ
      USE MOD_ANGMOM,ONLY:NKMMAX,NKMQ
      USE MOD_TYPES,ONLY:NTMAX,CONC
      USE MOD_CALCMODE,ONLY:IREL,KMROT,KKRMODE,THERMAL_VIBRA_FLUCT
      USE MOD_KSPACE,ONLY:IBZINT,NPTMIN,NPTMAX
      IMPLICIT NONE
C*--TAU_DRIVE17
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='TAU_DRIVE')
C
C Dummy arguments
C
      REAL*8 CPACHNG
      COMPLEX*16 ERYD,P
      INTEGER ICPACONV,ICPAFLAG,IECURR,IPRINTBAND,ITCPA
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           PHASK(IECURR),TAUQ(NKMMAX,NKMMAX,NQMAX),
     &           TSSQ(NKMMAX,NKMMAX,NQMAX),TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      INTEGER IQ,IT,KTAUIJ,N
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
C***********************************************************************
      IF ( THERMAL_VIBRA_FLUCT ) THEN
         IF ( KKRMODE(1:12).NE.'STANDARD-KKR' )
     &        CALL STOP_MESSAGE(ROUTINE,
     &        'THERMAL_VIBRA_FLUCT but KKRMODE <> STANDARD-KKR')
         IF ( IREL.EQ.2 ) CALL STOP_MESSAGE(ROUTINE,
     &        'THERMAL_VIBRA_FLUCT but IREL.EQ. 2')
C         IF ( IBZINT.NE.2 ) CALL STOP_MESSAGE(ROUTINE,
C     &        'THERMAL_VIBRA_FLUCT but IBZINT.NE.2 ')
      END IF
C***********************************************************************
C
C
C=======================================================================
C                   REAL SPACE  type calculation
C=======================================================================
C
C***********************************************************************
      IF ( KKRMODE(1:10).EQ.'REAL-SPACE' ) THEN
C***********************************************************************
C
         IF ( IBZINT.NE.0 ) CALL STOP_MESSAGE(ROUTINE,
     &        'KKRMODE = REAL-SPACE  and  IBZINT <> 0')
C
         KTAUIJ = 0
C
         CALL TAU_REAL_SPACE(KTAUIJ,IPRINTBAND,ERYD,TSSQ,MSSQ,TAUQ,TSST,
     &                       MSST)
C
C
C=======================================================================
C                STANDARD-KKR  type calculation
C=======================================================================
C
C***********************************************************************
      ELSE IF ( KKRMODE(1:12).EQ.'STANDARD-KKR' ) THEN
C***********************************************************************
C
C=======================================================================
C                   BZ - integration using  WEYL - scheme
C=======================================================================
C
         IF ( IBZINT.EQ.1 ) THEN
C
C-----------------------------------------------------------------------
C               copy t(IT) = TAU(IQ) to get single site DOS
C-----------------------------------------------------------------------
C
            IF ( NPTMIN.EQ.0 .OR. NPTMAX.EQ.0 ) THEN
C
               DO IQ = 1,NQ
                  IT = ITOQ(1,IQ)
                  N = NKMQ(IQ)
                  TAUQ(1:N,1:N,IQ) = TSST(1:N,1:N,IT)
               END DO
C
            ELSE
C
               CALL TAU_STD_KWEYL(ERYD,P,IPRINTBAND,ICPAFLAG,ITCPA,
     &                            ICPACONV,CPACHNG,TSST,MSST,TSSQ,MSSQ,
     &                            TAUQ)
C
            END IF
C
C=======================================================================
C               BZ - integration using  SPECIAL POINTS
C=======================================================================
C
         ELSE IF ( IBZINT.EQ.2 ) THEN
C
            IF ( IREL.NE.2 ) THEN
C
               CALL TAU_STD_KPOINTS(ICPAFLAG,CPACHNG,ERYD,P,IPRINTBAND,
     &                              ITCPA,ICPACONV,CONC,NOQ,ITOQ,PHASK,
     &                              IECURR,NTMAX,TSST,MSST,TSSQ,MSSQ,
     &                              TAUQ)
C
C-----------------------------------------------------------------------
C          BZ - integration using  SPECIAL POINTS  for  SPIN SPRIRALs
C-----------------------------------------------------------------------
C
            ELSE IF ( KMROT.GE.3 ) THEN
C
               CALL TAU_STD_KSPIRAL(ICPAFLAG,CPACHNG,ERYD,P,IPRINTBAND,
     &                              ITCPA,ICPACONV,TSST,MSST,TSSQ,MSSQ,
     &                              TAUQ)
C
C-----------------------------------------------------------------------
C               BZ - integration using  SPECIAL POINTS
C               for NON- or SCALAR relativistic case
C-----------------------------------------------------------------------
C
            ELSE
C
               CALL TAU_STD_KSPSREL(ICPAFLAG,CPACHNG,ERYD,P,IPRINTBAND,
     &                              ITCPA,ICPACONV,PHASK,IECURR,TSST,
     &                              MSST,TSSQ,MSSQ,TAUQ)
            END IF
C
C=======================================================================
C               BZ - integration using  TETRAHEDRON METHOD
C=======================================================================
C
         ELSE IF ( IBZINT.EQ.3 ) THEN
C
            CALL TAU_STD_KTETRA(ICPAFLAG,CPACHNG,ERYD,P,IPRINTBAND,
     &                          ITCPA,ICPACONV,PHASK,IECURR,TSST,MSST,
     &                          TSSQ,MSSQ,TAUQ)
C
C=======================================================================
C               BZ - integration using  ZOOM-IN METHOD
C=======================================================================
C
         ELSE IF ( IBZINT.EQ.4 ) THEN
C
            WRITE (6,*) '###############################'
            WRITE (6,*) '###############################'
            WRITE (6,*) '##### UNDER CONSTRUCTION  #####'
            WRITE (6,*) '###############################'
            WRITE (6,*) '###############################'
            STOP 'TAU_STD_KZOOM deactivated'
C            CALL TAU_STD_KZOOM(ICPAFLAG,CPACHNG,ERYD,P,IPRINTBAND,TAUQ,
C     &                         ITCPA,ICPACONV,MSSQ,MSST)
C
         ELSE
            WRITE (6,99001) IBZINT
            CALL STOP_MESSAGE(ROUTINE,
     &                  'KKRMODE = STANDARD-KKR and IBZINT out of range'
     &                  )
         END IF
C
C
C=======================================================================
C                    TB-KKR  type calculation
C=======================================================================
C
C***********************************************************************
      ELSE IF ( KKRMODE(1:6).EQ.'TB-KKR' ) THEN
C***********************************************************************
C
         CALL TAU_TB(ERYD,P,ICPAFLAG,CPACHNG,ITCPA,ICPACONV,TSST,MSST,
     &               TSSQ,MSSQ,TAUQ)
C
C
C=======================================================================
C                 EMBEDDED-CLUSTER  type calculation
C=======================================================================
C
C***********************************************************************
      ELSE IF ( KKRMODE(1:16).EQ.'EMBEDDED-CLUSTER' ) THEN
C***********************************************************************
C
         CALL TAUIJ_EMBEDDED(IPRINTBAND,ERYD,ICPAFLAG,CPACHNG,ITCPA,
     &                       ICPACONV,IECURR,TSST,MSST,TSSQ,MSSQ,TAUQ)
C
C***********************************************************************
      ELSE
C
         CALL STOP_MESSAGE(ROUTINE,'KKRMODE not set properly')
C
      END IF
C***********************************************************************
C
99001 FORMAT (//,1X,79('#'),/,10X,'IBZINT=',i4,/)
      END
