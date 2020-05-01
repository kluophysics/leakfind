C*==linresp_magnet_eloop.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_MAGNET_ELOOP(RHO2NSX,USE_TAU_DK)
C   ********************************************************************
C   *                                                                  *
C   *  run the E-loop for linear response calculations                 *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LINRESP,ONLY:NZ12MAX,RHO2NS_GG,T0Z,T1Z,TZ,TIJZ,TZ_STD,
     &    NPERT,NOBSE,TZ_DK,T1Z_DK,T0Z_DK
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_CPA,ONLY:CPATOL,NCPA
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,WKM1,MEZJ,MEZZ,TSSQ,MSSQ,MSST,SSST,
     &    TAUQ,TAUT,TSST
      USE MOD_SITES,ONLY:NQMAX,IQAT
      USE MOD_TYPES,ONLY:NT,NTMAX,NLMFPMAX
      USE MOD_CALCMODE,ONLY:ORBPOL,TASK
      USE MOD_ENERGY,ONLY:ETAB,NEMAX,NETAB,WETAB,PHASK
      USE MOD_FILES,ONLY:IFILCBWF
      USE MOD_MPI,ONLY:MPI,MPI_ID,MPI_KLOOP,MPI_ELOOP
      IMPLICIT NONE
C*--LINRESP_MAGNET_ELOOP23
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='LINRESP_MAGNET_ELOOP')
C
C Dummy arguments
C
      LOGICAL USE_TAU_DK
      REAL*8 RHO2NSX(NRMAX,NLMFPMAX,NTMAX,3)
C
C Local variables
C
      LOGICAL CALCINT,GETIRRSOL
      REAL*8 CPACHNG,CPACHTAB(NEMAX),TIME1,TIME2
      COMPLEX*16 CWKE(:),DMATT(:,:,:),DTILT(:,:,:),ERYD,P,
     &           TAUQZ(NKMMAX,NKMMAX,NQMAX,NZ12MAX),WE
      INTEGER I,IA_ERR,ICPACONV,ICPAFLAG,IE,IECPAFAIL(NEMAX),IOBSE,
     &        IPERT,IPROCE(:),IQ,IT,ITCPA,IWRI,IWRIRRWF,IWRREGWF,LWKE,
     &        LWKEMAX,NCPAFAIL
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DMATT,DTILT,IPROCE,CWKE
C
      ALLOCATE (DMATT(NKMMAX,NKMMAX,NTMAX),IPROCE(NEMAX))
      ALLOCATE (DTILT(NKMMAX,NKMMAX,NTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate DTILT')
C
C ======================================================================
C
      IWRI = 6
C
      IPRINT = 1
      CALCINT = .TRUE.
      GETIRRSOL = .TRUE.
      IWRREGWF = 1
      IWRIRRWF = 1
C
C=======================================================================
C
      NCPAFAIL = 0
C
      CALL CPU_TIME(TIME1)
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      IF ( .NOT.(MPI) ) THEN
         MPI_KLOOP = .FALSE.
         MPI_ELOOP = .FALSE.
      ELSE IF ( MPI ) THEN
         MPI_ELOOP = .TRUE.
         MPI_KLOOP = .FALSE.
      END IF
C
      CALL MPI_DISTRIBUTE(IPROCE,NETAB(1),MPI_ELOOP,'E')
C
      IF ( MPI ) CALL DRV_MPI_BARRIER
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      RHO2NSX(:,:,:,:) = 0.0D0
      RHO2NS_GG(:,:,:,:,:) = 0.0D0
C
      T0Z(:,:,:) = 0.0D0
      T1Z(:,:,:) = 0.0D0
      TZ(:,:,:) = 0.0D0
      TIJZ(:,:,:,:) = 0.0D0
      TZ_DK(:,:,:,:) = 0.0D0
      T1Z_DK(:,:,:,:) = 0.0D0
      T0Z_DK(:,:,:,:) = 0.0D0
C
      TZ_STD(:,:,:) = 0.0D0
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop    START
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
      DO IE = 1,NETAB(1)
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         IF ( MPI_ID.EQ.IPROCE(IE) ) THEN
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
            ERYD = ETAB(IE,1)
            WE = WETAB(IE,1)
C
            ICPAFLAG = 0
            CPACHNG = 0.0D0
C
            CALL CINIT(NKMMAX*NKMMAX*NTMAX,TAUT)
C
C ===================================== solve SS - differential equation
C
            CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFILCBWF,GETIRRSOL,
     &                    ERYD,P,IPRINT,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
            CALL MSSINIT(TSST,MSST,TSSQ,MSSQ)
C
C ======================================================================
C
            IF ( NCPA.GT.0 ) CALL TAU_DRIVE(IE,IPRINT,ERYD,P,TSSQ,MSSQ,
     &           TSST,MSST,TAUQ,ICPAFLAG,CPACHNG,ITCPA,ICPACONV,PHASK)
C
            IF ( USE_TAU_DK ) THEN
C
               CALL LINRESP_KLOOP_DK(ERYD,P,TAUQ,TAUQZ,MSSQ)
C
            ELSE
C
               CALL SIGKLOOPS(ERYD,P,ERYD,P,TAUQ,TAUQ,TAUQZ,MSSQ,MSSQ)
C
            END IF
C
C ======================================================================
C
            CALL PROJTAU(ICPAFLAG,CPACHNG,ERYD,MSST,MSSQ,TAUQZ,TAUT)
C
            DO IT = 1,NT
C
               IQ = IQAT(1,IT)
C
               CALL GETDMAT(TAUQ(1,1,IQ),DMATT(1,1,IT),DTILT(1,1,IT),
     &                      WKM1,NKM,MSSQ(1,1,IQ),MSST(1,1,IT),NKMMAX)
            END DO
C
C ======================================================================
C
            CALL LINRESP_MAGNET_PERTURB(IFILCBWF,IFILCBWF,ERYD)
C
C ======================================================================
C
            CALL LINRESP_MAGNET_RESPONSE(IFILCBWF,MSST,TAUT,DMATT,DTILT,
     &         IFILCBWF,MSST,TAUT,DMATT,DTILT,MEZZ,MEZJ,WE,ERYD)
C
C
            DO I = 1,3
               IF ( TASK.EQ.'MAGNET    ' ) THEN
                  IOBSE = 1
                  IF ( I.EQ.1 ) THEN
                     IPERT = 1
                  ELSE IF ( I.EQ.2 ) THEN
                     IPERT = 4
                  ELSE
                     IPERT = 7
                  END IF
               ELSE
                  IOBSE = I
                  IPERT = I
               END IF
               WRITE (*,*) 'T0Z ',IE,I,T0Z(1,IOBSE,IPERT)
               WRITE (*,*) 'TZ  ',IE,I,TZ(1,IOBSE,IPERT)
            END DO
C
C ======================================================================
C
            IF ( ICPAFLAG.NE.0 ) THEN
               NCPAFAIL = NCPAFAIL + 1
               CPACHTAB(NCPAFAIL) = CPACHNG
               IECPAFAIL(NCPAFAIL) = IE
            END IF
C
C ======================================================================
            IF ( IPRINT.GE.3 ) CALL DUMPTAU(IE,ERYD,IWRI,MSST,MSSQ,TAUT,
     &           TAUQ)
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      END DO
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop    END
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C             collect results for energy points of a E-path
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      IF ( MPI ) THEN
C
         CALL DRV_MPI_BARRIER
C
         LWKEMAX = NTMAX*NTMAX*NOBSE*NPERT
         ALLOCATE (CWKE(LWKEMAX))
C
         LWKE = NTMAX*NOBSE*NPERT
C
         CALL DRV_MPI_REDUCE_C(T0Z(1,1,1),CWKE(1),LWKE)
         CALL DRV_MPI_REDUCE_C(T1Z(1,1,1),CWKE(1),LWKE)
         CALL DRV_MPI_REDUCE_C(TZ(1,1,1),CWKE(1),LWKE)
         CALL DRV_MPI_REDUCE_C(TZ_STD(1,1,1),CWKE(1),LWKE)
         CALL DRV_MPI_REDUCE_C(TIJZ(1,1,1,1),CWKE(1),LWKEMAX)
C
         CALL DRV_MPI_BARRIER
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      CALL CPU_TIME(TIME2)
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      IF ( MPI_ID.EQ.0 ) THEN
C
         WRITE (IWRI,99001) TIME2 - TIME1
C
         IF ( NCPAFAIL.NE.0 ) THEN
            WRITE (IWRI,99002) CPATOL,NCPAFAIL,
     &                         (IECPAFAIL(IE),DREAL(ETAB(IECPAFAIL(IE),
     &                         1)),CPACHTAB(IE),IE=1,NCPAFAIL)
            WRITE (IWRI,'(1X,79(''*''),/)')
         ELSE IF ( NCPA.NE.0 ) THEN
            WRITE (IWRI,99003)
         END IF
C
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
99001 FORMAT (/,' CPU - TIME  USED      ',F10.3,/)
99002 FORMAT (/,1X,79('*'),/,' tolerance for CPA-cycle:',F15.7,/,
     &        ' CPA not converged for',I3,' energies:',/,
     &        3(' E:',I3,F7.4,' C:',F8.5,:,2X))
99003 FORMAT (/,1X,79('*'),/,25X,'no problems with','  CPA-cycle !',/,
     &        1X,79('*'),/)
      END
C*==linresp_magnet_standard.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_MAGNET_STANDARD(IT,WE,TMAT)
C   ********************************************************************
C   *                                                                  *
C   *        calculate expectation values via the standard scheme      *
C   *                                                                  *
C   *                    <A> = -1/PI Im Trace A G                      *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LINRESP,ONLY:NPERT,TZ_STD,HAZ_ZZ_T,HCZ_ZJ_T
      USE MOD_ANGMOM,ONLY:NKMMAX,NKMQ,WKM1,WKM2,WKM3
      USE MOD_SITES,ONLY:IQAT
      USE MOD_CONSTANTS,ONLY:PI
      USE MOD_TYPES,ONLY:IMT
      USE MOD_RMESH,ONLY:JRCRI
      IMPLICIT NONE
C*--LINRESP_MAGNET_STANDARD284
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IT
      COMPLEX*16 WE
      COMPLEX*16 TMAT(NKMMAX,NKMMAX)
C
C Local variables
C
      COMPLEX*16 CMATTRC
      INTEGER IKM1,IKM2,IM,IME,IPERT,IQ,IRCRIT,M,N
C
C*** End of declarations rewritten by SPAG
C
      IQ = IQAT(1,IT)
      IM = IMT(IT)
      IRCRIT = JRCRI(IM)
C
C--------------------- site diagonal multiple scattering matrix G or TAU
C
      M = NKMMAX
      N = NKMQ(IQ)
C
C-----------------------------------------------------------------------
C       calculate standard expectation values for site IQ -- type IT
C-----------------------------------------------------------------------
C
      DO IPERT = 1,NPERT
         IF ( IPERT.EQ.1 ) THEN
            IME = 1
         ELSE IF ( IPERT.LE.4 ) THEN
            IME = 2
         ELSE IF ( IPERT.LE.7 ) THEN
            IME = 3
         ELSE
            IME = IPERT
         END IF
C
         DO IKM1 = 1,N
            DO IKM2 = 1,N
               WKM3(IKM1,IKM2) = HAZ_ZZ_T(IRCRIT,IKM1,IKM2,IT,IPERT)
               WKM2(IKM1,IKM2) = HCZ_ZJ_T(1,IKM1,IKM2,IT,IPERT)
            END DO
         END DO
C
C?????????????????????????????????????????????????????????? TO BE CHECKED
C should be coherent with changes in: subr. LINRESP_MAGNET
C         IF ( IME.LE.3 ) THEN
C         WKM3(1:N,1:N) = MEZZ(1:N,1:N,IT,IME)
C         WKM2(1:N,1:N) = MEZJ(1:N,1:N,IT,IME)
C         END IF
C
         WKM1(1:N,1:N) = MATMUL(WKM3(1:N,1:N),TMAT(1:N,1:N))
C
         WKM2(1:N,1:N) = WKM1(1:N,1:N) - WKM2(1:N,1:N)
C
         TZ_STD(1,1,IPERT) = TZ_STD(1,1,IPERT) - (WE/PI)
     &                       *CMATTRC(N,M,WKM2)
C
      END DO
C
      END
C*==linresp_magnet_resp_sum.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_MAGNET_RESP_SUM(IT,MEZZ,MEZJ,TMATA,TMATB,
     &   THT_DIA_P,THT_OFF_TP,D0,T0,D1,T1,D,T,DIJ,TIJ,WE,IPERT1,IPERT2)
C   ********************************************************************
C   *                                                                  *
C   *   evaluate the 4 contributions to  Delta G  for type IT          *
C   *                                                                  *
C   *   perform partial summations and integrate D* -> T*              *
C   *                                                                  *
C   *   CPA: the type projection for all quantitiers has been done     *
C   *                                                                  *
C   *   depending on GF_CONV_RH alle quantities are defined            *
C   *   in the RH- or ZJ-convention: e.g. TMAT = TAUT or G00[IT]       *
C   *                                                                  *
C   *   THT_DIA_P: site diagonal perturbation                         *
C   *               TMAT(IT) * H(IT) * TMAT(IT)                        *
C   *                                                                  *
C   *   THT_OFF_TP: site off-diagonal perturbation                     *
C   *        sum{JQ}  TMAT(IQ,IT;JQ,JT) * H(JT) * TMAT(JQ,JT;IQ,IT)    *
C   *                                                                  *
C   ********************************************************************
      USE MOD_LINRESP,ONLY:NPERT,NOBSE,MZBZA_O,MIRR2_OP,MIRR3_OP,
     &    MIRR4_OP,DOBS_TO
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,NKMQ,WKM1,WKM2,NMEMAX
      USE MOD_SITES,ONLY:IQAT
      USE MOD_CONSTANTS,ONLY:C0,C1,PI
      USE MOD_TYPES,ONLY:NTMAX,ITBOT,ITTOP
      IMPLICIT NONE
C*--LINRESP_MAGNET_RESP_SUM388
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IPERT1,IPERT2,IT
      COMPLEX*16 WE
      COMPLEX*16 D(NTMAX,NOBSE,NPERT),D0(NTMAX,NOBSE,NPERT),
     &           D1(NTMAX,NOBSE,NPERT),DIJ(NTMAX,NTMAX,NOBSE,NPERT),
     &           MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),T(NTMAX,NOBSE,NPERT),
     &           T0(NTMAX,NOBSE,NPERT),T1(NTMAX,NOBSE,NPERT),
     &           THT_DIA_P(NKMMAX,NKMMAX,NPERT),
     &           THT_OFF_TP(NKMMAX,NKMMAX,NTMAX,NPERT),
     &           TIJ(NTMAX,NTMAX,NOBSE,NPERT),TMATA(NKMMAX,NKMMAX),
     &           TMATB(NKMMAX,NKMMAX)
C
C Local variables
C
      COMPLEX*16 CMATTRC
      COMPLEX*16 D1DIA,D1OFF,D2,D3,D4
      INTEGER I,IME,IOBSE,IPERT,IQ,J,JT,M,N
C
C*** End of declarations rewritten by SPAG
C
      IQ = IQAT(1,IT)
      M = NKMMAX
      N = NKMQ(IQ)
C
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      DO IOBSE = 1,NOBSE
C
C-----------------------------------------------------------------------
C       calculate standard densities DOS for type IT energy A
C-----------------------------------------------------------------------
C
C------------ in doubt: the equivalency of the indices has to be checked
         IME = IOBSE
C
         WKM1(1:N,1:N) = MATMUL(MEZZ(1:N,1:N,IT,IME),TMATA(1:N,1:N))
C
         WKM2(1:N,1:N) = -(WKM1(1:N,1:N)-MEZJ(1:N,1:N,IT,IME))/PI
C
         DOBS_TO(IT,IOBSE) = CMATTRC(N,M,WKM2)
C
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
         DO IPERT = IPERT1,IPERT2
C
C-----------------------------------------------------------------------
C   calculate  G * G = Delta_G_off(1) + Delta_G_dia(1)
C                    + Delta_G(2) + Delta_G(3) + Delta_G(4)
C-----------------------------------------------------------------------
C
C------------------------- Delta G_off(1) site off-diagonal contribution
C
C
            D1OFF = 0D0
            DO JT = ITBOT,ITTOP
C
               WKM1(1:N,1:N) = MATMUL(MZBZA_O(1:N,1:N,IOBSE),THT_OFF_TP(
     &                         1:N,1:N,JT,IPERT))
C
               DIJ(IT,JT,IOBSE,IPERT) = CMATTRC(NKM,NKMMAX,WKM1)
               D1OFF = D1OFF + DIJ(IT,JT,IOBSE,IPERT)
            END DO
C
C
C----------------------------- Delta G_dia(1) site diagonal contribution
C
            WKM1(1:N,1:N) = MATMUL(MZBZA_O(1:N,1:N,IOBSE),THT_DIA_P(1:N,
     &                      1:N,IPERT))
C
            D1DIA = CMATTRC(NKM,NKMMAX,WKM1)
C
C------------------------------------------------------------ Delta G(2)
C
            CALL ZGEMM('N','N',N,N,N,C1,MIRR2_OP(1,1,IOBSE,IPERT),M,
     &                 TMATB(1,1),M,C0,WKM1,M)
C
            D2 = -CMATTRC(NKM,NKMMAX,WKM1)
C
C------------------------------------------------------------ Delta G(3)
C
            WKM2(1:N,1:N) = TRANSPOSE(TMATA(1:N,1:N))
C
            CALL ZGEMM('N','N',N,N,N,C1,WKM2,M,MIRR3_OP(1,1,IOBSE,IPERT)
     &                 ,M,C0,WKM1,M)
C
            D3 = -CMATTRC(NKM,NKMMAX,WKM1)
C
C------------------------------------------------------------ Delta G(4)
C
            D4 = 0D0
            DO I = 1,NKM
               DO J = 1,NKM
                  D4 = D4 + MIRR4_OP(I,J,IOBSE,IPERT)
               END DO
            END DO
C
C-----------------------------------------------------------------------
C               perform summations and integrations
C-----------------------------------------------------------------------
C
            DO JT = ITBOT,ITTOP
               TIJ(IT,JT,IOBSE,IPERT) = TIJ(IT,JT,IOBSE,IPERT)
     &                                  + DIJ(IT,JT,IOBSE,IPERT)*WE/PI
               IF ( IT.EQ.JT ) TIJ(IT,JT,IOBSE,IPERT)
     &              = TIJ(IT,JT,IOBSE,IPERT) + D1DIA*WE/PI
            END DO
C
            D1(IT,IOBSE,IPERT) = D1OFF + D1DIA
C
            T1(IT,IOBSE,IPERT) = T1(IT,IOBSE,IPERT) + D1(IT,IOBSE,IPERT)
     &                           *WE/PI
C
            D0(IT,IOBSE,IPERT) = D2 + D3 + D4
C
            T0(IT,IOBSE,IPERT) = T0(IT,IOBSE,IPERT) + D0(IT,IOBSE,IPERT)
     &                           *WE/PI
C
            D(IT,IOBSE,IPERT) = D0(IT,IOBSE,IPERT) + D1(IT,IOBSE,IPERT)
            T(IT,IOBSE,IPERT) = T0(IT,IOBSE,IPERT) + T1(IT,IOBSE,IPERT)
C
         END DO
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
C
      END DO
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
C
      END
