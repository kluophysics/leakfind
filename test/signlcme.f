C*==signlcme.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIGNLCME(ERYDA,ERYDB,IFILA,IFILB,IQCPA,NQNLCPA,NDIMCLU,
     &                    IND0QCLU,NCFG,PCFG,W1HAT,MSSTA,MSSTB,MUEHATA,
     &                    MUEHATB,OMEGAHATA,OMEGAHATB,MAQQAB,MAQQBA,
     &                    MBQQAB,MBQQBA,MCQQAB,MCQQBA,MDQQAB,MDQQBA,
     &                    MTAB,MTBA)
C   ********************************************************************
C   *                                                                  *
C   *   - calculate the SIGMA matrix elements                          *
C   *   - transform the matrix elements for different E-arguments      *
C   *                                                                  *
C   *   the l-expansion is set via   NKM   for all atomic types to NL  *
C   *                                                                  *
C   *   ALL matrix elements refer to the local frame of reference      *
C   *                                                                  *
C   *   NOTE:  this routine supplies the matrix elements MAQQAB ....   *
C   *          for the NLCPA - case assuming                           *
C   *          - disorder only on site  IQCPA                          *
C   *          - IQCLUREP: site for which matrix elements are taken    *
C   *                                                                  *
C   *    THIS ROUTINE IS RESTRICTED TO NQMAX = 1 AT THE MOMENT !!!!!   *
C   *                                                                  *
C   ********************************************************************
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,NKMQ,WKM1,MZAZB,MZBZA,MIRR_2,
     &    MIRR_3,MIRR_4,NPOLMAX
      USE MOD_SITES,ONLY:ITOQ,NQMAX,IQAT,IQBOT_CHI,IQTOP_CHI
      USE MOD_TYPES,ONLY:CTL,NTMAX,NT
      USE MOD_FILES,ONLY:IPRINT,FOUND_SECTION
      USE MOD_CONSTANTS,ONLY:C1,C0
      IMPLICIT NONE
C*--SIGNLCME31
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIGNLCME')
      LOGICAL ME_CC_BRA_RWF
      PARAMETER (ME_CC_BRA_RWF=.FALSE.)
      INTEGER IWME
      PARAMETER (IWME=0)
C
C Dummy arguments
C
      COMPLEX*16 ERYDA,ERYDB
      INTEGER IFILA,IFILB,IQCPA,NCFG,NDIMCLU,NQNLCPA
      INTEGER IND0QCLU(NQNLCPA)
      COMPLEX*16 MAQQAB(NKMMAX,NKMMAX,3,NQNLCPA,NQNLCPA),
     &           MAQQBA(NKMMAX,NKMMAX,3,NQNLCPA,NQNLCPA),
     &           MBQQAB(NKMMAX,NKMMAX,3,NQNLCPA,NQNLCPA),
     &           MBQQBA(NKMMAX,NKMMAX,3,NQNLCPA,NQNLCPA),
     &           MCQQAB(NKMMAX,NKMMAX,3,NQNLCPA,NQNLCPA),
     &           MCQQBA(NKMMAX,NKMMAX,3,NQNLCPA,NQNLCPA),
     &           MDQQAB(NKMMAX,NKMMAX,3,NQNLCPA,NQNLCPA),
     &           MDQQBA(NKMMAX,NKMMAX,3,NQNLCPA,NQNLCPA),
     &           MSSTA(NKMMAX,NKMMAX,NTMAX),MSSTB(NKMMAX,NKMMAX,NTMAX),
     &           MTAB(NKMMAX,NKMMAX,3,NTMAX),MTBA(NKMMAX,NKMMAX,3,NTMAX)
     &           ,MUEHATA(NDIMCLU,NDIMCLU),MUEHATB(NDIMCLU,NDIMCLU),
     &           OMEGAHATA(NDIMCLU,NDIMCLU),OMEGAHATB(NDIMCLU,NDIMCLU),
     &           W1HAT(NDIMCLU,NDIMCLU)
      REAL*8 PCFG(NCFG)
C
C Local variables
C
      REAL*8 C,T
      COMPLEX*16 CWGT,DMATA(:,:,:,:),DMATB(:,:,:,:),DTILA(:,:,:,:),
     &           DTILB(:,:,:,:),TAUIMPA(:,:),TAUIMPB(:,:),W2HAT(:,:)
      INTEGER I1,IA_ERR,ICFG,IOCCREP,IQ,IQCLU,IQCLUREP,IT,ITREP,J,J0,
     &        JQCLU,M,MUE,N,NN
      LOGICAL INITIALIZE
      CHARACTER*3 MEFORM
      INTEGER NLCPACONF
      SAVE MEFORM
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
C
      ALLOCATABLE DMATA,DMATB,DTILA,DTILB
      ALLOCATABLE TAUIMPA,TAUIMPB,W2HAT
C
      ALLOCATE (W2HAT(NDIMCLU,NDIMCLU))
      ALLOCATE (TAUIMPA(NDIMCLU,NDIMCLU))
      ALLOCATE (TAUIMPB(NDIMCLU,NDIMCLU))
      ALLOCATE (DMATA(NKMMAX,NKMMAX,NQNLCPA,NQNLCPA))
      ALLOCATE (DMATB(NKMMAX,NKMMAX,NQNLCPA,NQNLCPA))
      ALLOCATE (DTILA(NKMMAX,NKMMAX,NQNLCPA,NQNLCPA))
      ALLOCATE (DTILB(NKMMAX,NKMMAX,NQNLCPA,NQNLCPA),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'alloc DTILB')
      ALLOCATE (MZAZB(NKMMAX,NKMMAX,NPOLMAX))
      ALLOCATE (MZBZA(NKMMAX,NKMMAX,NPOLMAX))
      ALLOCATE (MIRR_2(NKMMAX,NKMMAX,NPOLMAX,NPOLMAX))
      ALLOCATE (MIRR_3(NKMMAX,NKMMAX,NPOLMAX,NPOLMAX))
      ALLOCATE (MIRR_4(NKMMAX,NKMMAX,NPOLMAX,NPOLMAX))
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MZAZB')
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
      IF ( INITIALIZE ) THEN
C
         WRITE (6,99001) ROUTINE(1:LEN_TRIM(ROUTINE))
C
         M = NKMMAX
C
         CALL INPUT_FIND_SECTION('TASK',0)
C
C        MEFORM = 'NAB'
         MEFORM = 'ADA'
         IF ( FOUND_SECTION ) CALL SECTION_SET_STRING('ME',MEFORM,
     &        '9999',0)
C
         WRITE (6,'(A10,A)') '   MEFORM ',MEFORM
C
         IF ( (MEFORM.NE.'NAB') .AND. (MEFORM.NE.'ADA') )
     &        CALL STOP_MESSAGE(ROUTINE,'MEFORM not found')
C
C ------------------------------------ calculate angular matrix elements
C
         CALL AME_INIT(MEFORM,6)
C
         INITIALIZE = .FALSE.
C
      END IF
C=======================================================================
C
C   ********************************************************************
C   *                                                                  *
C   * read wave function and calculate matrix elements of momentum     *
C   * operator used to determine specific resistance                   *
C   *                                                                  *
C   ********************************************************************
C
      DO IT = 1,NT
         M = NKMMAX
         N = NKM
         IQ = IQAT(1,IT)
C
         WRITE (6,99002) IT
C
         C = CTL(IT,1)
C
C-------------------------------------------------- energies E_a and E_b
C
         CALL ME_DRIVE(1,NKM,IFILB,ERYDB,1,NKM,IFILA,ERYDA,
     &                 ME_CC_BRA_RWF,IT,MTAB(1,1,1,IT),MTBA(1,1,1,IT),
     &                 MIRR_2,MIRR_3,MIRR_4,C,2,MEFORM,IWME)
C
C-----------------------------------------------------------------------
C           transformation from spherical to cartesian coordinates
C              from now on        ipol = 1,2,3  ==  (x),(y),(z)
C-----------------------------------------------------------------------
C
         CALL CMAT_CONVERT_POLAR(MTAB(1,1,1,IT),'S>C')
C
         CALL CMAT_CONVERT_POLAR(MTBA(1,1,1,IT),'S>C')
C
C-----------------------------------------------------------------------
         IF ( IPRINT.EQ.5 ) THEN
            T = 1D-8
            WRITE (6,99003) MEFORM
            DO MUE = 1,3
               WRITE (6,99006) IT,MUE
               CALL CMATSTRUCT('MTAB  ',MTAB(1,1,MUE,IT),N,M,3,3,0,T,6)
               CALL CMATSTRUCT('MTBA  ',MTBA(1,1,MUE,IT),N,M,3,3,0,T,6)
            END DO
C
         END IF
C-----------------------------------------------------------------------
C
      END DO
C                                                                     IT
C   ********************************************************************
C
C ======================================================================
C     calculate configurational average of matrix element MAQQAB etc.
C ======================================================================
C
      NN = NKMMAX*NKMMAX*3*NQNLCPA*NQNLCPA
C
      CALL CINIT(NN,MAQQAB)
      CALL CINIT(NN,MBQQAB)
      CALL CINIT(NN,MCQQAB)
      CALL CINIT(NN,MDQQAB)
C
      CALL CINIT(NN,MAQQBA)
      CALL CINIT(NN,MBQQBA)
      CALL CINIT(NN,MCQQBA)
      CALL CINIT(NN,MDQQBA)
C
      DO IQ = IQBOT_CHI,IQTOP_CHI
         IF ( IQ.EQ.IQCPA ) THEN
            DO ICFG = 1,NCFG
C
C-------------------------------- calculate projection matrices D and D~
C    from   TAU_imp = D TAU^ = TAU^ D~  and  Omega^ = mue^ - TAU^**-1
C
               N = NKMQ(IQCPA)
               M = NDIMCLU
C
               CALL NLCPAPROJ(ICFG,IQCPA,MSSTA,OMEGAHATA,TAUIMPA,W1HAT,
     &                        NQNLCPA,N,IND0QCLU,ITOQ,NDIMCLU,NQMAX,
     &                        NTMAX,NKMMAX)
C
               W1HAT(1:M,1:M) = MUEHATA(1:M,1:M) - OMEGAHATA(1:M,1:M)
               CALL CMATMUL(M,M,TAUIMPA,W1HAT,W2HAT)
               DO IQCLU = 1,NQNLCPA
                  I1 = IND0QCLU(IQCLU) + 1
                  DO JQCLU = 1,NQNLCPA
                     J0 = IND0QCLU(JQCLU)
                     DO J = 1,N
                        CALL ZCOPY(N,W2HAT(I1,J0+J),1,
     &                             DMATA(1,J,IQCLU,JQCLU),1)
                     END DO
                  END DO
               END DO
C
               CALL CMATMUL(M,M,W1HAT,TAUIMPA,W2HAT)
               DO IQCLU = 1,NQNLCPA
                  I1 = IND0QCLU(IQCLU) + 1
                  DO JQCLU = 1,NQNLCPA
                     J0 = IND0QCLU(JQCLU)
                     DO J = 1,N
                        CALL ZCOPY(N,W2HAT(I1,J0+J),1,
     &                             DTILA(1,J,IQCLU,JQCLU),1)
                     END DO
                  END DO
               END DO
C
               CALL NLCPAPROJ(ICFG,IQCPA,MSSTB,OMEGAHATB,TAUIMPB,W1HAT,
     &                        NQNLCPA,N,IND0QCLU,ITOQ,NDIMCLU,NQMAX,
     &                        NTMAX,NKMMAX)
C
               W1HAT(1:M,1:M) = MUEHATB(1:M,1:M) - OMEGAHATB(1:M,1:M)
C
               CALL CMATMUL(M,M,TAUIMPB,W1HAT,W2HAT)
               DO IQCLU = 1,NQNLCPA
                  I1 = IND0QCLU(IQCLU) + 1
                  DO JQCLU = 1,NQNLCPA
                     J0 = IND0QCLU(JQCLU)
                     DO J = 1,N
                        CALL ZCOPY(N,W2HAT(I1,J0+J),1,
     &                             DMATB(1,J,IQCLU,JQCLU),1)
                     END DO
                  END DO
               END DO
C
               CALL CMATMUL(M,M,W1HAT,TAUIMPB,W2HAT)
               DO IQCLU = 1,NQNLCPA
                  I1 = IND0QCLU(IQCLU) + 1
                  DO JQCLU = 1,NQNLCPA
                     J0 = IND0QCLU(JQCLU)
                     DO J = 1,N
                        CALL ZCOPY(N,W2HAT(I1,J0+J),1,
     &                             DTILB(1,J,IQCLU,JQCLU),1)
                     END DO
                  END DO
               END DO
C
C-----------------------------------------------------------------------
               DO IQCLUREP = 1,NQNLCPA
C
                  CWGT = PCFG(ICFG)/DBLE(NQNLCPA)
C
                  IOCCREP = NLCPACONF(ICFG,NQNLCPA,IQCLUREP,1,2)
                  ITREP = ITOQ(IOCCREP,IQCPA)
C
                  IF ( IPRINT.EQ.5 ) WRITE (6,99004) IQ
C
                  M = NKMMAX
                  DO IQCLU = 1,NQNLCPA
                     DO JQCLU = 1,NQNLCPA
                        DO MUE = 1,3
C
C----------------------------------------------------- MBAR type A and B
C
                           CALL ZGEMM('N','N',N,N,N,C1,
     &                                DTILA(1,1,IQCLU,IQCLUREP),M,
     &                                MTAB(1,1,MUE,ITREP),M,C0,WKM1,M)
C
                           CALL ZGEMM('N','N',N,N,N,CWGT,WKM1,M,
     &                                DMATB(1,1,IQCLUREP,JQCLU),M,C1,
     &                                MAQQAB(1,1,MUE,IQCLU,JQCLU),M)
C
                           CALL ZGEMM('N','C',N,N,N,CWGT,WKM1,M,
     &                                DTILB(1,1,JQCLU,IQCLUREP),M,C1,
     &                                MBQQAB(1,1,MUE,IQCLU,JQCLU),M)
C
C
                           CALL ZGEMM('N','N',N,N,N,C1,
     &                                DTILB(1,1,IQCLU,IQCLUREP),M,
     &                                MTBA(1,1,MUE,ITREP),M,C0,WKM1,M)
C
                           CALL ZGEMM('N','N',N,N,N,CWGT,WKM1,M,
     &                                DMATA(1,1,IQCLUREP,JQCLU),M,C1,
     &                                MAQQBA(1,1,MUE,IQCLU,JQCLU),M)
C
                           CALL ZGEMM('N','C',N,N,N,CWGT,WKM1,M,
     &                                DTILA(1,1,JQCLU,IQCLUREP),M,C1,
     &                                MBQQBA(1,1,MUE,IQCLU,JQCLU),M)
C
C----------------------------------------------------- MBAR type C and D
C
                           CALL ZGEMM('C','N',N,N,N,C1,
     &                                DMATA(1,1,IQCLUREP,IQCLU),M,
     &                                MTAB(1,1,MUE,ITREP),M,C0,WKM1,M)
C
                           CALL ZGEMM('N','N',N,N,N,CWGT,WKM1,M,
     &                                DMATB(1,1,IQCLUREP,JQCLU),M,C1,
     &                                MCQQAB(1,1,MUE,IQCLU,JQCLU),M)
C
                           CALL ZGEMM('N','C',N,N,N,CWGT,WKM1,M,
     &                                DTILB(1,1,JQCLU,IQCLUREP),M,C1,
     &                                MDQQAB(1,1,MUE,IQCLU,JQCLU),M)
C
C
                           CALL ZGEMM('C','N',N,N,N,C1,
     &                                DMATB(1,1,IQCLUREP,IQCLU),M,
     &                                MTBA(1,1,MUE,ITREP),M,C0,WKM1,M)
C
                           CALL ZGEMM('N','N',N,N,N,CWGT,WKM1,M,
     &                                DMATA(1,1,IQCLUREP,JQCLU),M,C1,
     &                                MCQQBA(1,1,MUE,IQCLU,JQCLU),M)
C
                           CALL ZGEMM('N','C',N,N,N,CWGT,WKM1,M,
     &                                DTILA(1,1,JQCLU,IQCLUREP),M,C1,
     &                                MDQQBA(1,1,MUE,IQCLU,JQCLU),M)
C
C-------------------------------------------------------------------- IO
C
                        END DO
                     END DO
                  END DO
C
               END DO
C
            END DO
C-----------------------------------------------------------------------
            IF ( IPRINT.EQ.5 ) THEN
               DO IQCLU = 1,NQNLCPA
                  DO JQCLU = 1,NQNLCPA
                     DO MUE = 1,3
                        WRITE (6,99005) MUE,IQCLU,JQCLU
                        T = 1D-8
                        CALL CMATSTRUCT('MAQQAB ',
     &                                  MAQQAB(1,1,MUE,IQCLU,JQCLU),N,M,
     &                                  3,3,0,T,6)
                        CALL CMATSTRUCT('MBQQAB ',
     &                                  MBQQAB(1,1,MUE,IQCLU,JQCLU),N,M,
     &                                  3,3,0,T,6)
                        CALL CMATSTRUCT('MCQQAB ',
     &                                  MCQQAB(1,1,MUE,IQCLU,JQCLU),N,M,
     &                                  3,3,0,T,6)
                        CALL CMATSTRUCT('MDQQAB ',
     &                                  MDQQAB(1,1,MUE,IQCLU,JQCLU),N,M,
     &                                  3,3,0,T,6)
                     END DO
                  END DO
               END DO
            END IF
C-----------------------------------------------------------------------
         END IF
      END DO
C
99001 FORMAT (//,1X,79('*'),/,35X,'<',A,'>',/,1X,79('*'),/)
99002 FORMAT (//,' J - matrix elements for component IT=',I3)
99003 FORMAT (/,10X,'electric dipole - matrix elements evaluated',
     &        ' using the  ',A,' * A - form ',/)
99004 FORMAT (/,'JBAR-Matrices IQ=',I2,/,44('='))
99005 FORMAT (/,'direction (xyz)=(123)',i2,5x,' IQCLU=',i2,' JQCLU=',i2)
99006 FORMAT (/,'     IT =',i3,' MUE =',i3)
      END
