C*==linresp_magnet.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_MAGNET
C   ********************************************************************
C   *                                                                  *
C   *   calculate STATIC susceptibility                                *
C   *                                                                  *
C   *  indexing of operators and observables                           *
C   *                                                                  *
C   *  number of particles   1        1 IDOS                           *
C   *  spin moment           b s_z    2 ISPN                           *
C   *  orbital moment        b l_z    3 IORB   ____ NPERT              *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_FILES,ONLY:IPRINT,IFILBUILDBOT,WRBUILDBOT
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,IND0Q,NKKR
      USE MOD_SITES,ONLY:NQMAX,NQ,IQBOT_CHI,IQTOP_CHI
      USE MOD_TYPES,ONLY:NTMAX,NLMFPMAX,NT,CONC,NAT,TXT_T
      USE MOD_MPI,ONLY:MPI,MPI_ID
      USE MOD_LINRESP,ONLY:HZ_PERT_LMTP,T0Z,T1Z,TIJZ,TZ,K_PERT_LMTP,
     &    CHI_TO,NZ12MAX,ITTA,ITTB,ITTC,ITTD,ITTQ1,ITTQ2,NTTJ,JTT1,JTT2,
     &    JTTX,WTTJ,ITTMAX,JTTMAX,NTKTKMAX,NTKTKLIN,CHIZ,CHI_DK,TZ_STD,
     &    TZ_DK
      IMPLICIT NONE
C*--LINRESP_MAGNET26
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='LINRESP_MAGNET')
      LOGICAL USE_TAU_DK
      PARAMETER (USE_TAU_DK=.TRUE.)
C
C Local variables
C
      INTEGER I,IA_ERR,IC,IQ,IT,J,JORB,JORB1,JORB2,JQ,JSPN,JSPN1,JSPN2,
     &        JT,L1,L2,L3,L4,M
      REAL*8 MAGMOM_SUM,MAGMOM_TOT,RHO2NSX(:,:,:,:),TIME1,TIME2
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE RHO2NSX
C
C ======================================================================
C
      CALL LINRESP_INIT('MAGNET    ')
C
      JSPN1 = 2
      JORB1 = 5
      JSPN2 = 4
      JORB2 = 7
C
C ======================================================================
C
C----------------------------------------------------- variables for CHI
C
      ALLOCATE (RHO2NSX(NRMAX,NLMFPMAX,NTMAX,3))
C
      WRITE (6,99001)
C
C=======================================================================
C
      CALL CPU_TIME(TIME1)
C
      IF ( MPI ) CALL DRV_MPI_BARRIER
C
      IF ( USE_TAU_DK ) THEN
C
C ======================================================================
C         set up index table for BZ integration for NO SYMMETRY case
C ======================================================================
C
         IF ( ALLOCATED(WTTJ) ) THEN
            DEALLOCATE (WTTJ,JTT1,JTT2,JTTX)
            DEALLOCATE (NTTJ,ITTA,ITTB,ITTC,ITTD,ITTQ1,ITTQ2)
         END IF
C
C------------------------------- consider ALL TAU(k)*TAU(k) combinations
         NTKTKMAX = NKM**4
C
         ITTMAX = NTKTKMAX*NQ*NQ
C
         ALLOCATE (NTTJ(ITTMAX))
         ALLOCATE (ITTA(ITTMAX),ITTB(ITTMAX),ITTC(ITTMAX))
         ALLOCATE (ITTD(ITTMAX),ITTQ1(ITTMAX),ITTQ2(ITTMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate ITTD')
C
         JTTMAX = ITTMAX
C
         ALLOCATE (WTTJ(JTTMAX))
         ALLOCATE (JTT1(JTTMAX),JTT2(JTTMAX),JTTX(JTTMAX))
C
         NTKTKLIN = ITTMAX
C
         I = 0
         DO IQ = IQBOT_CHI,IQTOP_CHI
            DO JQ = IQBOT_CHI,IQTOP_CHI
               DO L1 = 1,NKM
                  DO L4 = 1,NKM
                     DO L2 = 1,NKM
                        DO L3 = 1,NKM
                           I = I + 1
                           ITTA(I) = L1
                           ITTB(I) = L2
                           ITTC(I) = L3
                           ITTD(I) = L4
                           ITTQ1(I) = IQ
                           ITTQ2(I) = JQ
                           NTTJ(I) = 1
C
                           J = I
                           JTT1(J) = (IND0Q(JQ)+L2-1)*NKKR + IND0Q(IQ)
     &                               + L1
                           JTT2(J) = (IND0Q(IQ)+L4-1)*NKKR + IND0Q(JQ)
     &                               + L3
                           JTTX(J) = (IND0Q(JQ)+L3-1)*NKKR + IND0Q(IQ)
     &                               + L4
                           WTTJ(J) = 1D0
C
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
C
C ======================================================================
C                         re-allocate CHIZ
C ======================================================================
C
         DEALLOCATE (CHIZ)
C
         NZ12MAX = 2
         M = NKMMAX
C
         ALLOCATE (CHIZ(M*M*NQMAX,M*M*NQMAX,NZ12MAX),STAT=IA_ERR)
C
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate CHIZ')
C
C
C ======================================================================
C                   initialize k-mesh for FULL BZ
C ======================================================================
C
         CALL INIT_MOD_TAUIJ_KMESH
C
      END IF
C
C ======================================================================
C                         allocate CHI_DK
C ======================================================================
C
      M = NKMMAX
C
      ALLOCATE (CHI_DK(M*M*NQMAX,M*M*NQMAX,3),STAT=IA_ERR)
C
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate CHI_DK')
C
C=======================================================================
C
      HZ_PERT_LMTP(:,:,:,2) = HZ_PERT_LMTP(:,:,:,1)
      HZ_PERT_LMTP(:,:,:,3) = HZ_PERT_LMTP(:,:,:,1)
      K_PERT_LMTP(:,:,2) = K_PERT_LMTP(:,:,1)
      K_PERT_LMTP(:,:,3) = K_PERT_LMTP(:,:,1)
C
      CHI_TO(:,:) = 0D0
      CHI_TO(:,:) = 0D0
      CHI_TO(:,:) = 0D0
C
      CALL LINRESP_MAGNET_ELOOP(RHO2NSX,USE_TAU_DK)
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      IF ( MPI_ID.EQ.0 ) THEN
C
C ======================================================================
C   calculate atom type-resolved  SPIN  and  ORBITAL  magnetic moment
C ======================================================================
C-----------------------------------------------------------------------
         IF ( IPRINT.GE.1 ) THEN
            JSPN = JSPN1
            JORB = JORB1
C
            WRITE (6,99003) (IT,DIMAG(T1Z(IT,1,JSPN)),DIMAG(T0Z(IT,1,
     &                      JSPN)),IT=1,NT)
C
            DO IT = 1,NT
               WRITE (6,'('' IT='',I2)') IT
               WRITE (6,99004) (JT,'Z',DIMAG(TIJZ(IT,JT,1,JSPN)),JT=1,
     &                         NT)
            END DO
C
            WRITE (6,99005) (IT,DIMAG(T1Z(IT,1,JORB)),DIMAG(T0Z(IT,1,
     &                      JORB)),IT=1,NT)
C
            DO IT = 1,NT
               WRITE (6,'('' IT='',I2)') IT
               WRITE (6,99006) (JT,DIMAG(TIJZ(IT,JT,1,JORB)),JT=1,NT)
            END DO
         END IF
C
C ======================================================================
C                           print out
C ======================================================================
C
         TZ(:,:,:) = -TZ(:,:,:)
C
         MAGMOM_TOT = 0D0
C
         DO IT = 1,NT
C
C------------------------------------------------------------------- SPN
            CALL LINRESP_MAGNET_POL(TZ,IT,1)
C ----------------------------------------------------------------------
C SWITCHING ON SHOULD CORRELATE WITH: subroutine LINRESP_MAGNET_STANDARD
C together with  WKM3(IKM1,IKM2) = HAZ_ZZ_T !!!!!!!!!!!!!!!!!!!!
            CALL LINRESP_MAGNET_POL(TZ_STD,IT,1)
C ----------------------------------------------------------------------
C
C------------------------------------------------------------------- ORB
            CALL LINRESP_MAGNET_POL(TZ,IT,4)
C ----------------------------------------------------------------------
C SWITCHING ON SHOULD CORRELATE WITH: subroutine LINRESP_MAGNET_STANDARD
C together with  WKM3(IKM1,IKM2) = HAZ_ZZ_T !!!!!!!!!!!!!!!!!!!!
            CALL LINRESP_MAGNET_POL(TZ_STD,IT,4)
C ----------------------------------------------------------------------
C
            CALL LINRESP_MAGNET_POLXX(TZ,IT)
            CALL LINRESP_MAGNET_POLXX(TZ_STD,IT)
C
            DO IC = 1,3
C
C----------------------------------------------------- LAN  alf  A and B
C
               CALL LINRESP_MAGNET_POL(TZ_DK(1,1,1,IC),IT,16)
C
C----------------------------------------------------- LAN  nab  A and B
C
               CALL LINRESP_MAGNET_POL(TZ_DK(1,1,1,IC),IT,19)
C
            END DO
C
            WRITE (6,99008) IT,TXT_T(IT)
C
            JSPN = JSPN1
            JORB = JORB1
            MAGMOM_SUM = DIMAG(TZ(IT,1,JSPN)) + DIMAG(TZ(IT,1,JORB))
C
            MAGMOM_TOT = MAGMOM_TOT + NAT(IT)*CONC(IT)*MAGMOM_SUM
C
C         WRITE (6,99010) 'MAGMOM S+O :',
C     &                   (DIMAG(TZ(IT,1,JSPN)+TZ(IT,1,JSPN+3)),
C     &                   JSPN=JSPN1,JSPN2)
C
C-----------------------------------------------------------------------
C
            WRITE (6,99010) 'STD    S   :',
     &                      (DIMAG(TZ_STD(IT,1,JSPN)),JSPN=JSPN1,JSPN2)
            WRITE (6,99009) 'G*G    S   :',
     &                      (DIMAG(TZ(IT,1,JSPN)),JSPN=JSPN1,JSPN2)
C
C-----------------------------------------------------------------------
C
            WRITE (6,99010) 'STD    O   :',
     &                      (DIMAG(TZ_STD(IT,1,JORB)),JORB=JORB1,JORB2)
            WRITE (6,99009) 'G*G    O   :',
     &                      (DIMAG(TZ(IT,1,JORB)),JORB=JORB1,JORB2)
C
C-----------------------------------------------------------------------
C
            WRITE (6,99010) 'STD r x a A:',DIMAG(TZ_STD(IT,1,15)),
     &                      DIMAG(TZ_STD(IT,1,10)),
     &                      DIMAG(TZ_STD(IT,1,11))
C
            WRITE (6,99009) 'G*G r x a A:',DIMAG(TZ(IT,1,15)),
     &                      DIMAG(TZ(IT,1,10)),DIMAG(TZ(IT,1,11))
C
C-----------------------------------------------------------------------
C
            WRITE (6,99010) 'STD r x a B:',DIMAG(TZ_STD(IT,1,13)),
     &                      DIMAG(TZ_STD(IT,1,14)),
     &                      DIMAG(TZ_STD(IT,1,09))
C
            WRITE (6,99009) 'G*G r x a B:',DIMAG(TZ(IT,1,13)),
     &                      DIMAG(TZ(IT,1,14)),DIMAG(TZ(IT,1,09))
C
C-----------------------------------------------------------------------
C
            WRITE (6,99009) 'LAN  alf  A:',DIMAG(TZ_DK(IT,1,19,2)),
     &                      DIMAG(TZ_DK(IT,1,17,3)),
     &                      DIMAG(TZ_DK(IT,1,18,1))
C
            WRITE (6,99009) 'LAN  alf  B:',DIMAG(TZ_DK(IT,1,18,3)),
     &                      DIMAG(TZ_DK(IT,1,19,1)),
     &                      DIMAG(TZ_DK(IT,1,17,2))
C
C-----------------------------------------------------------------------
C
            WRITE (6,99009) 'LAN  nab  A:',DIMAG(TZ_DK(IT,1,22,2)),
     &                      DIMAG(TZ_DK(IT,1,20,3)),
     &                      DIMAG(TZ_DK(IT,1,21,1))
C
            WRITE (6,99009) 'LAN  nab  B:',DIMAG(TZ_DK(IT,1,21,3)),
     &                      DIMAG(TZ_DK(IT,1,22,1)),
     &                      DIMAG(TZ_DK(IT,1,20,2))
C
            WRITE (6,'(/,1X,79(''-''))')
C
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
            IF ( WRBUILDBOT ) THEN
C
               WRITE (IFILBUILDBOT,99011) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                'STD    S   :',IT,
     &                (DIMAG(TZ_STD(IT,1,JSPN)),JSPN=JSPN1,JSPN2)
C
               WRITE (IFILBUILDBOT,99011) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                'G*G    S   :',IT,
     &                (DIMAG(TZ(IT,1,JSPN)),JSPN=JSPN1,JSPN2)
C
               WRITE (IFILBUILDBOT,99011) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                'STD    O   :',IT,
     &                (DIMAG(TZ_STD(IT,1,JORB)),JORB=JORB1,JORB2)
C
               WRITE (IFILBUILDBOT,99011) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                'G*G    O   :',IT,
     &                (DIMAG(TZ(IT,1,JORB)),JORB=JORB1,JORB2)
C
               WRITE (IFILBUILDBOT,99011) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                'STD r x a A:',IT,DIMAG(TZ_STD(IT,1,15)),
     &                DIMAG(TZ_STD(IT,1,10)),DIMAG(TZ_STD(IT,1,11))
C
               WRITE (IFILBUILDBOT,99011) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                'G*G r x a A:',IT,DIMAG(TZ(IT,1,15)),
     &                DIMAG(TZ(IT,1,10)),DIMAG(TZ(IT,1,11))
C
               WRITE (IFILBUILDBOT,99011) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                'STD r x a B:',IT,DIMAG(TZ_STD(IT,1,13)),
     &                DIMAG(TZ_STD(IT,1,14)),DIMAG(TZ_STD(IT,1,09))
C
C
               WRITE (IFILBUILDBOT,99011) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                'G*G r x a B:',IT,DIMAG(TZ(IT,1,13)),
     &                DIMAG(TZ(IT,1,14)),DIMAG(TZ(IT,1,09))
C
               WRITE (IFILBUILDBOT,99011) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                'LAN  alf  A:',IT,DIMAG(TZ_DK(IT,1,19,2)),
     &                DIMAG(TZ_DK(IT,1,17,3)),DIMAG(TZ_DK(IT,1,18,1))
C
               WRITE (IFILBUILDBOT,99011) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                'LAN  alf  B:',IT,DIMAG(TZ_DK(IT,1,18,3)),
     &                DIMAG(TZ_DK(IT,1,19,1)),DIMAG(TZ_DK(IT,1,17,2))
C
               WRITE (IFILBUILDBOT,99011) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                'LAN  nab  A:',IT,DIMAG(TZ_DK(IT,1,22,2)),
     &                DIMAG(TZ_DK(IT,1,20,3)),DIMAG(TZ_DK(IT,1,21,1))
C
               WRITE (IFILBUILDBOT,99011) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                'LAN  nab  B:',IT,DIMAG(TZ_DK(IT,1,21,3)),
     &                DIMAG(TZ_DK(IT,1,22,1)),DIMAG(TZ_DK(IT,1,20,2))
C
            END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
C-----------------------------------------------------------------------
         END DO
C
         IF ( NT.EQ.1 ) THEN
            WRITE (6,'(/,1X,79(''=''))')
         ELSE
            WRITE (6,99007) MAGMOM_TOT
         END IF
C
C ======================================================================
C
         CALL CPU_TIME(TIME2)
C
         WRITE (6,99002) TIME2 - TIME1
C
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      CALL STOP_REGULAR(ROUTINE,' ')
C
99001 FORMAT (///,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*        *     *   ***    ****   *     *  *****  *****       *'
     &  ,/,10X,
     &  '*        **   **  *   *  *    *  **    *  *        *         *'
     &  ,/,10X,
     &  '*        * * * *  *   *  *       * *   *  *        *         *'
     &  ,/,10X,
     &  '*        *  *  *  *****  *   **  *  *  *  ****     *         *'
     &  ,/,10X,
     &  '*        *     *  *   *  *    *  *   * *  *        *         *'
     &  ,/,10X,
     &  '*        *     *  *   *  *    *  *    **  *        *         *'
     &  ,/,10X,
     &  '*        *     *  *   *   ****   *     *  *****    *         *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99002 FORMAT (/,' CPU - TIME  USED      ',F10.3,/)
99003 FORMAT (/,8(:,' IT=',I2,' T1SSZ   =',1E13.5,' T0SSZ   =',1E13.5,/)
     &        )
99004 FORMAT (8(:,' JT=',I2,' TIJSS',A,'  =',1E13.5,/))
99005 FORMAT (/,16(:,' IT=',I2,' T1LLZS  =',1E13.5,' T0LLZS  =',1E13.5,
     &        /))
99006 FORMAT (32(:,' JT=',I2,' TIJLLZS =',1E13.5,/))
99007 FORMAT (/,1X,79('='),//,10X,'total magnetic moment',F15.6,' m_B',
     &        //,1X,79('='),/)
99008 FORMAT (/,1X,79('='),//,10X,'results for atom type  IT=',I2,2X,A)
C$$$99009 FORMAT (/,10X,A12,3F15.6,5X)
C$$$99010 FORMAT (//,10X,A12,3F15.6,5X)
C------------------------------------------- TMP: Temporary used format
99009 FORMAT (/,10X,A12,3F20.12,5X)
99010 FORMAT (//,10X,A12,3F20.12,5X)
99011 FORMAT ('# BUILDBOT: ',A,2X,A12,'  for IT =',I5,/,(1PE22.14))
C----------------------------------------------------------------------
      END
C*==linresp_magnet_pol.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
C
      SUBROUTINE LINRESP_MAGNET_POL(T,IT,IPERT0)
C   ********************************************************************
C   *                                                                  *
C   *     change polarisation from SPHERICAL to CARTESIAN              *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:SQRT_2,CI
      USE MOD_LINRESP,ONLY:NPERT,NOBSE
      USE MOD_TYPES,ONLY:NTMAX
      IMPLICIT NONE
C*--LINRESP_MAGNET_POL444
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IPERT0,IT
      COMPLEX*16 T(NTMAX,NOBSE,NPERT)
C
C Local variables
C
      COMPLEX*16 AUX_0,AUX_M,AUX_P,AUX_X,AUX_Y
C
C*** End of declarations rewritten by SPAG
C
C      CHARACTER*40 ROUTINE
C      PARAMETER (ROUTINE='LINRESP_MAGNET_POL')
C
      AUX_M = T(IT,1,IPERT0+1)
      AUX_0 = T(IT,1,IPERT0+2)
      AUX_P = T(IT,1,IPERT0+3)
C
C-----------------------------------------------------------------------
C  To be checked:
C  alpha and nabla operators have different conventions for polariztion
C-----------------------------------------------------------------------
      AUX_X = (+AUX_P+AUX_M)/SQRT_2
      AUX_Y = CI*(-AUX_P+AUX_M)/SQRT_2
C$$$      AUX_X = (-AUX_P+AUX_M)/SQRT_2
C$$$      AUX_Y = CI*(AUX_P+AUX_M)/SQRT_2
C
      T(IT,1,IPERT0+1) = AUX_X
      T(IT,1,IPERT0+2) = AUX_Y
      T(IT,1,IPERT0+3) = AUX_0
C
      END
C*==linresp_magnet_polxx.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
C
      SUBROUTINE LINRESP_MAGNET_POLXX(T,IT)
C   ********************************************************************
C   *                                                                  *
C   *     change polarisation from SPHERICAL to CARTESIAN              *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:SQRT_2,CI,PI
      USE MOD_LINRESP,ONLY:NPERT,NOBSE
      USE MOD_TYPES,ONLY:NTMAX
      IMPLICIT NONE
C*--LINRESP_MAGNET_POLXX504
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IT
      COMPLEX*16 T(NTMAX,NOBSE,NPERT)
C
C Local variables
C
      REAL*8 NXY,NZ
      COMPLEX*16 TX_0,TX_M,TX_P,TX_X,TX_Y,TX_Z,TY_0,TY_M,TY_P,TY_X,TY_Y,
     &           TY_Z,TZ_0,TZ_M,TZ_P,TZ_X,TZ_Y,TZ_Z
C
C*** End of declarations rewritten by SPAG
C
C      CHARACTER*40 ROUTINE
C      PARAMETER (ROUTINE='LINRESP_MAGNET_POL')
C
      NXY = SQRT(4D0*PI/3D0)
      NZ = SQRT(4D0*PI/3D0)
C
      TX_M = T(IT,1,10)*NXY
      TY_M = T(IT,1,08)*NXY
      TZ_M = T(IT,1,09)*NZ
C
      TX_0 = T(IT,1,13)*NXY
      TY_0 = T(IT,1,11)*NXY
      TZ_0 = T(IT,1,12)*NZ
C
      TX_P = T(IT,1,16)*NXY
      TY_P = T(IT,1,14)*NXY
      TZ_P = T(IT,1,15)*NZ
C
      TX_X = (-TX_P+TX_M)/SQRT_2
      TX_Y = CI*(TX_P+TX_M)/SQRT_2
      TX_Z = TX_0
C
      TY_X = (-TY_P+TY_M)/SQRT_2
      TY_Y = CI*(TY_P+TY_M)/SQRT_2
      TY_Z = TY_0
C
      TZ_X = (-TZ_P+TZ_M)/SQRT_2
      TZ_Y = CI*(TZ_P+TZ_M)/SQRT_2
      TZ_Z = TZ_0
C
      T(IT,1,08) = TX_X
      T(IT,1,09) = TY_X
      T(IT,1,10) = TZ_X
C
      T(IT,1,11) = TX_Y
      T(IT,1,12) = TY_Y
      T(IT,1,13) = TZ_Y
C
      T(IT,1,14) = TX_Z
      T(IT,1,15) = TY_Z
      T(IT,1,16) = TZ_Z
C
      END
