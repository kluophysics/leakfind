C*==linresp_eloop.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_ELOOP(MEZJ,MEZZ,TSSQ,MSSQ,MSST,SSST,TAUQ,TAUT,
     &                         TSST,RHO2NSX,PHASK)
C   ********************************************************************
C   *                                                                  *
C   *  run the E-loop for linear response calculations                 *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LINRESP,ONLY:NZ12MAX,RHO2NS_GG,LINRESP_CHECK_SUM_RULES,
     &    T0Z,T0X,T1Z,T1X,TZ,TX,TIJZ,TIJX,QVEC_PERT_EQ_0VEC
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_CPA,ONLY:CPATOL,NCPA
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX,NKM,WKM1
      USE MOD_SITES,ONLY:NQMAX,IQAT
      USE MOD_TYPES,ONLY:NT,NTMAX,NLMFPMAX
      USE MOD_CALCMODE,ONLY:ORBPOL,IREL
      USE MOD_ENERGY,ONLY:ETAB,NEMAX,NETAB,WETAB
      USE MOD_FILES,ONLY:IFILCBWF
      USE MOD_MPI,ONLY:MPI,MPI_ID,MPI_KLOOP,MPI_ELOOP
      IMPLICIT NONE
C*--LINRESP_ELOOP23
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='LINRESP_ELOOP')
C
C Dummy arguments
C
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           PHASK(NEMAX),SSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX),TAUT(NKMMAX,NKMMAX,NTMAX),
     &           TSSQ(NKMMAX,NKMMAX,NQMAX),TSST(NKMMAX,NKMMAX,NTMAX)
      REAL*8 RHO2NSX(NRMAX,NLMFPMAX,NTMAX,3)
C
C Local variables
C
      LOGICAL CALCINT,GETIRRSOL
      REAL*8 CPACHNG,CPACHTAB(NEMAX),TIME1,TIME2
      COMPLEX*16 DMATT(:,:,:),DTILT(:,:,:),ERYD,P,
     &           TAUQZ(NKMMAX,NKMMAX,NQMAX,NZ12MAX),WE
      INTEGER I,IA_ERR,ICPACONV,ICPAFLAG,IE,IECPAFAIL(NEMAX),IQ,IT,
     &        ITCPA,IWRI,IWRIRRWF,IWRREGWF,NCPAFAIL
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DMATT,DTILT
C
      ALLOCATE (DMATT(NKMMAX,NKMMAX,NTMAX))
      ALLOCATE (DTILT(NKMMAX,NKMMAX,NTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate DTILT')
C
C ======================================================================
C
      WRITE (*,*) 'LINRESP_CHECK_SUM_RULES ',LINRESP_CHECK_SUM_RULES
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
         MPI_ELOOP = .FALSE.
         MPI_KLOOP = .TRUE.
      END IF
C
      IF ( MPI ) CALL DRV_MPI_BARRIER
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      RHO2NSX(:,:,:,:) = 0.0D0
      RHO2NS_GG(:,:,:,:,:) = 0.0D0
C
      T0Z(:,:,:) = 0.0D0
      T0X(:,:,:) = 0.0D0
      T1Z(:,:,:) = 0.0D0
      T1X(:,:,:) = 0.0D0
      TZ(:,:,:) = 0.0D0
      TX(:,:,:) = 0.0D0
      TIJZ(:,:,:,:) = 0.0D0
      TIJX(:,:,:,:) = 0.0D0
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop    START
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
      DO IE = 1,NETAB(1)
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
     &                 ERYD,P,IPRINT,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
         CALL MSSINIT(TSST,MSST,TSSQ,MSSQ)
C
C ======================================================================
C
         IF ( NCPA.GT.0 ) CALL TAU_DRIVE(IE,IPRINT,ERYD,P,TSSQ,MSSQ,
     &        TSST,MSST,TAUQ,ICPAFLAG,CPACHNG,ITCPA,ICPACONV,PHASK)
C
         IF ( IREL.EQ.2 ) THEN
C
            CALL SIGKLOOPS_SPSREL(ERYD,P,ERYD,P,TAUQ,TAUQ,TAUQZ,MSSQ,
     &                            MSSQ)
C
C
         ELSE IF ( QVEC_PERT_EQ_0VEC ) THEN
C
            CALL SIGKLOOPS(ERYD,P,ERYD,P,TAUQ,TAUQ,TAUQZ,MSSQ,MSSQ)
C
         ELSE
C
            CALL LINRESP_KLOOPS_Q(ERYD,P,ERYD,P,TAUQ,TAUQ,TAUQZ,MSSQ,
     &                            MSSQ)
C
C
         END IF
C
C ======================================================================
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         IF ( MPI_ID.EQ.0 ) THEN
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
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
            CALL LINRESP_PERTURB(IFILCBWF,IFILCBWF)
C
C ======================================================================
C
            IF ( LINRESP_CHECK_SUM_RULES ) THEN
               WE = 1D0
               RHO2NSX(:,:,:,:) = 0D0
               RHO2NS_GG(:,:,:,:,:) = 0D0
            END IF
C
C ======================================================================
C
            CALL LINRESP_RESPONSE(IE,IFILCBWF,MSST,TAUT,DMATT,DTILT,
     &                            IFILCBWF,MSST,TAUT,DMATT,DTILT,MEZZ,
     &                            MEZJ,WE,RHO2NSX)
C
C
            DO I = 1,3
               WRITE (*,*) 'T0Z ',IE,I,T0Z(1,I,I)
               WRITE (*,*) 'T0X ',IE,I,T0X(1,I,I)
               WRITE (*,*) 'TZ  ',IE,I,TZ(1,I,I)
               WRITE (*,*) 'TX  ',IE,I,TX(1,I,I)
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
      CALL CPU_TIME(TIME2)
C
      WRITE (IWRI,99001) TIME2 - TIME1
C
      IF ( NCPAFAIL.NE.0 ) THEN
         WRITE (IWRI,99002) CPATOL,NCPAFAIL,
     &                      (IECPAFAIL(IE),DREAL(ETAB(IECPAFAIL(IE),1)),
     &                      CPACHTAB(IE),IE=1,NCPAFAIL)
         WRITE (IWRI,'(1X,79(''*''),/)')
      ELSE IF ( NCPA.NE.0 ) THEN
         WRITE (IWRI,99003)
      END IF
C
99001 FORMAT (/,' CPU - TIME  USED      ',F10.3,/)
99002 FORMAT (/,1X,79('*'),/,' tolerance for CPA-cycle:',F15.7,/,
     &        ' CPA not converged for',I3,' energies:',/,
     &        3(' E:',I3,F7.4,' C:',F8.5,:,2X))
99003 FORMAT (/,1X,79('*'),/,25X,'no problems with','  CPA-cycle !',/,
     &        1X,79('*'),/)
      END
C*==linresp_sumrule.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_SUMRULE(IT,IE,TMAT,MEZZ,MEZJ,RHO2NSX,D0Z,D1Z)
C   ********************************************************************
C   *                                                                  *
C   *   check the sum rule  d G / d E = - G G                          *
C   *                                                                  *
C   ********************************************************************
      USE MOD_LINRESP,ONLY:NPERT,NOBSE,RHO2NS_GG,
     &    LINRESP_CHECK_SUM_RULES,TZ,T0X,TIJX,CHI_TO
      USE MOD_RMESH,ONLY:DRDI,NRMAX,NRSFTOT,NSF,LMISF,FLMSF,JRMT,JRCRI,R
      USE MOD_ENERGY,ONLY:ETAB
      USE MOD_ANGMOM,ONLY:NKMMAX,NMEMAX,NKMQ,WKM1,WKM2
      USE MOD_SITES,ONLY:IQAT
      USE MOD_CONSTANTS,ONLY:C0,PI,SQRT_4PI
      USE MOD_TYPES,ONLY:NT,NTMAX,NLMFPMAX,IMT,KLMFP
      IMPLICIT NONE
C*--LINRESP_SUMRULE260
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IE,IT
      COMPLEX*16 D0Z(NTMAX,NOBSE,NPERT),D1Z(NTMAX,NOBSE,NPERT),
     &           MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),TMAT(NKMMAX,NKMMAX)
      REAL*8 RHO2NSX(NRMAX,NLMFPMAX,NTMAX,3)
C
C Local variables
C
      COMPLEX*16 CMATTRC
      REAL*8 DDOS_GG,ICHR,RHO2NS0_GG(:,:,:,:,:),RHO2NSX0(:,:,:,:),
     &       RHO2NSX_G(:,:,:,:),RHO2NS_GG_INT(:,:,:,:,:),RHOINC,RINT1(:)
     &       ,RSQ,WRINT
      COMPLEX*16 DELE,DELTA_G(:),DELTA_G0(:),DOS0_T(:),DOSINC,DOS_T(:),
     &           DOS_T_GG(:),SUM_X
      INTEGER ILMRHO,IM,IOBSE,IPERT,IQ,IR,IRCRIT,IRMTIN,IRSF,ISF,JT,LM,
     &        M,N,NLMRHO
      LOGICAL INITIALIZE
      SAVE DELTA_G,DELTA_G0,DOS0_T,DOS_T,DOS_T_GG,RHO2NS0_GG,RHO2NSX0,
     &     RHO2NSX_G,RHO2NS_GG_INT
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE RINT1
      ALLOCATABLE RHO2NSX0,RHO2NS0_GG,RHO2NS_GG_INT,RHO2NSX_G
      ALLOCATABLE DOS_T,DOS0_T,DOS_T_GG,DELTA_G,DELTA_G0
      ALLOCATE (RINT1(NRMAX))
C
      DATA INITIALIZE/.TRUE./
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
C
      IF ( INITIALIZE ) THEN
C
         ALLOCATE (DOS_T(NTMAX),DOS0_T(NTMAX),DOS_T_GG(NTMAX))
         ALLOCATE (DELTA_G0(NTMAX),DELTA_G(NTMAX))
         ALLOCATE (RHO2NSX0(NRMAX,NLMFPMAX,NTMAX,3))
         ALLOCATE (RHO2NSX_G(NRMAX,NLMFPMAX,NTMAX,3))
         ALLOCATE (RHO2NS0_GG(NRMAX,NLMFPMAX,NTMAX,NOBSE,NPERT))
         ALLOCATE (RHO2NS_GG_INT(NRMAX,NLMFPMAX,NTMAX,NOBSE,NPERT))
C
         DELTA_G0(:) = C0
         DOS_T_GG(:) = C0
         RHO2NSX0(:,:,:,:) = 0D0
         RHO2NS0_GG(:,:,:,:,:) = 0D0
         RHO2NS_GG_INT(:,:,:,:,:) = 0D0
C
         INITIALIZE = .FALSE.
C
      END IF
C=======================================================================
C
      IPERT = 1
      IOBSE = 1
C
      IQ = IQAT(1,IT)
      IM = IMT(IT)
      IRMTIN = JRMT(IM)
      IRCRIT = JRCRI(IM)
C
C--------------------- site diagonal multiple scattering matrix G or TAU
C
      M = NKMMAX
      N = NKMQ(IQ)
C
C
      IF ( LINRESP_CHECK_SUM_RULES ) THEN
C-----------------------------------------------------------------------
C             calculate standard DOS for site IQ -- type IT
C-----------------------------------------------------------------------
C
         WKM1(1:N,1:N) = MATMUL(MEZZ(1:N,1:N,IT,1),TMAT(1:N,1:N))
C
         WKM2(1:N,1:N) = -(WKM1(1:N,1:N)-MEZJ(1:N,1:N,IT,1))/PI
C
         DOS_T(IT) = CMATTRC(N,M,WKM2)
C
C-----------------------------------------------------------------------
C   calculate  G * G = Delta_G_off(1) + Delta_G_dia(1)
C                    + Delta_G(2) + Delta_G(3) + Delta_G(4)
C-----------------------------------------------------------------------
C
C------------------------- Delta G_off(1) site off-diagonal contribution
C----------------------------- Delta G_dia(1) site diagonal contribution
C
         DELTA_G(IT) = D1Z(1,1,IT) + D0Z(1,1,IT)
C
         IF ( IE.EQ.1 ) THEN
C
            DOS_T_GG(IT) = 0D0
            DOS0_T(IT) = DOS_T(IT)
            RHO2NS0_GG(:,:,:,:,IT) = 0D0
            RHO2NSX0(:,:,IT,:) = RHO2NSX(:,:,IT,:)
C
         ELSE
C
            DELE = ETAB(IE,1) - ETAB(IE-1,1)
C
            DOSINC = 0.5D0*(DELTA_G(IT)+DELTA_G0(IT))*DELE/PI
            DOS_T_GG(IT) = DOS_T_GG(IT) + DOSINC
C
            WRINT = DREAL(0.5D0*DELE/PI)
            ILMRHO = 0
            DO LM = 1,NLMFPMAX
               IF ( KLMFP(LM,IT).NE.0 ) THEN
                  ILMRHO = ILMRHO + 1
                  NLMRHO = ILMRHO
                  DO IOBSE = 1,NOBSE
                     IPERT = IOBSE
                     DO IR = 1,IRCRIT
                        RHOINC = WRINT*(RHO2NS0_GG(IR,LM,IT,IOBSE,IPERT)
     &                           +RHO2NS_GG(IR,LM,IT,IOBSE,IPERT))
                        RHO2NS_GG_INT(IR,ILMRHO,IT,IOBSE,IPERT)
     &                     = RHO2NS_GG_INT(IR,ILMRHO,IT,IOBSE,IPERT)
     &                     + RHOINC
C
                        RHO2NSX_G(IR,ILMRHO,IT,IOBSE)
     &                     = RHO2NSX(IR,LM,IT,IOBSE)
     &                     - RHO2NSX0(IR,LM,IT,IOBSE)
                     END DO
                  END DO
               END IF
            END DO
C
            WRITE (100+IT,99002) DREAL(ETAB(IE,1)),
     &                           DIMAG(DOS_T(IT)-DOS0_T(IT)),
     &                           DIMAG(DOS_T_GG(IT))
C
            WRITE (6,99001) DREAL(ETAB(IE,1)),IT,DOS_T(IT),DOS_T(IT)
     &                      - DOS0_T(IT),DOS_T_GG(IT),DOS_T_GG(IT)
     &                      - (DOS_T(IT)-DOS0_T(IT))
C
            DO IR = 1,IRCRIT
               WRITE (200+IT,99003) R(IR,IM),
     &                              (RHO2NSX_G(IR,ILMRHO,IT,1),ILMRHO=1,
     &                              NLMRHO),
     &                              (RHO2NS_GG_INT(IR,ILMRHO,1,1,IT),
     &                              ILMRHO=1,NLMRHO)
            END DO
C
         END IF
C
         DELTA_G0(IT) = DELTA_G(IT)
         RHO2NS0_GG(:,:,:,:,IT) = RHO2NS_GG(:,:,:,:,IT)
C
      END IF
C
C=======================================================================
C             check energy resolved induced densities
C=======================================================================
C
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      DO IOBSE = 1,NOBSE
         IPERT = IOBSE
C
         DO LM = 1,NLMFPMAX
            DO IR = 1,IRCRIT
               RSQ = R(IR,IM)*R(IR,IM)
               RHO2NSX(IR,LM,IT,1) = RHO2NS_GG(IR,LM,IT,IOBSE,IPERT)
     &                               *RSQ*DRDI(IR,IM)
            END DO
         END DO
C
C----- convolute RHO2NS_GG with shape functions to get spherical density
C
         DO IR = 1,IRMTIN
            RINT1(IR) = RHO2NSX(IR,1,IT,1)*SQRT_4PI
         END DO
C
         RINT1((IRMTIN+1):IRCRIT) = 0D0
         DO ISF = 1,NSF(IM)
            LM = LMISF(ISF,IM)
            IF ( LM.LE.NLMFPMAX ) THEN
C
               DO IRSF = 1,NRSFTOT(IM)
                  IR = IRSF + IRMTIN
C
                  RINT1(IR) = RINT1(IR) + FLMSF(IRSF,ISF,IM)
     &                        *RHO2NSX(IR,LM,IT,1)
C
               END DO
C
            END IF
         END DO
C
C ------------------------------------------- check integrated densities
C
         CALL RRADINT(IM,RINT1,ICHR)
C
         IF ( LINRESP_CHECK_SUM_RULES ) THEN
            DDOS_GG = DIMAG(DELTA_G(IT)/PI)
         ELSE
            DDOS_GG = DIMAG(TZ(IT,IOBSE,IPERT))
         END IF
C
         SUM_X = T0X(IT,IOBSE,IPERT)
         DO JT = 1,NT
            SUM_X = SUM_X + TIJX(IT,JT,IOBSE,IPERT)
         END DO
C
         DDOS_GG = DDOS_GG + DIMAG(SUM_X)*CHI_TO(IT,IOBSE)
C
         DDOS_GG = DDOS_GG
C
         WRITE (6,99004) IT,IOBSE,DDOS_GG,ICHR,DDOS_GG/ICHR
C
      END DO
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
C
99001 FORMAT (/,' E = ',F7.3,' IT = ',I2,/,'  G      (STD)',2F12.6,/,
     &        '  G - G0 (STD)',2F12.6,/,'  Int dE G * G',2F12.6,4X,
     &        'Diff ',2F12.6)
99002 FORMAT (20F12.5)
99003 FORMAT (20E14.5)
99004 FORMAT (/,' IT ',I3,2X,'OBS:',I3,2X,'DNS',F20.10,/,18X,'INT',
     &        F20.10,5X,'RAT',F20.14)
      END
C*==linresp_resp_sum.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_RESP_SUM(IT,MEZZ,MEZJ,TMATA,TMATB,THT_DIA_P,
     &                            THT_OFF_TP,D0,T0,D1,T1,D,T,DIJ,TIJ,WE)
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
C*--LINRESP_RESP_SUM523
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IT
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
         DO IPERT = 1,NPERT
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
