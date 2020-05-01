C*==chrdns_norm.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHRDNS_NORM(RENORMALIZE,SCLNOS,TOTDOS,TOTNOS,TOTNOSX,
     &                       IEPANEL,IEPATH,IECURR,OBS_LT,OBS_T,OBS_LTX,
     &                       OBS_TX,RHOCHRX,RHOSPNX,RHOORBX,RHO2NSX,
     &                       CMNTTX)
C   ********************************************************************
C   *                                                                  *
C   * - renormalize GF-related quantities                              *
C   * - add the PATH-related results (*_LTX etc) to the global arrays  *
C   *                                                                  *
C   *   SCF_CHECK_SPLITSS:  print single site and back scattering      *
C   *                       contribution in case of  SPLITSS           *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ENERGY,ONLY:SPLITSS,ETAB,NETAB
      USE MOD_ANGMOM,ONLY:NOBSMAX,IDOS,ISMT,IOMT,IHFF,IBND,NLMAX
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_RMESH,ONLY:NRMAX,JRWS,FULLPOT,JRCRI
      USE MOD_TYPES,ONLY:ITBOT,ITTOP,CMNTT,CONC,NAT,IMT,NTMAX,RHOCHR,
     &    RHOSPN,RHOORB,RHO2NS,NLMFPMAX,KLMFP,TXT_T,NLT,NT,OBS_X
      USE MOD_SITES,ONLY:NLQMAD,NLMQMAD
      USE MOD_SCF,ONLY:SCF_CHECK_SPLITSS
      IMPLICIT NONE
C*--CHRDNS_NORM25
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IECURR,IEPANEL,IEPATH
      LOGICAL RENORMALIZE
      REAL*8 SCLNOS,TOTDOS,TOTNOS,TOTNOSX
      REAL*8 CMNTTX(NLMFPMAX,NTMAX),OBS_LT(0:3,NOBSMAX,NLMAX,NTMAX),
     &       OBS_LTX(0:3,NOBSMAX,NLMAX,NTMAX),OBS_T(0:3,NOBSMAX,NTMAX),
     &       OBS_TX(0:3,NOBSMAX,NTMAX),RHO2NSX(NRMAX,NLMFPMAX,NTMAX,3),
     &       RHOCHRX(NRMAX,NTMAX),RHOORBX(NRMAX,NTMAX),
     &       RHOSPNX(NRMAX,NTMAX)
C
C Local variables
C
      INTEGER I1,IL,IM,IR,IRTOP,IT,IW,LM,NDENS
C
C*** End of declarations rewritten by SPAG
C
C=======================================================================
C                     normalize GF-related quantities
C=======================================================================
C
      IF ( RENORMALIZE ) THEN
C
         WRITE (6,99001) SCLNOS
C
         TOTNOSX = TOTNOSX*SCLNOS
C
         DO IT = ITBOT,ITTOP
C
            OBS_TX(:,:,IT) = OBS_TX(:,:,IT)*SCLNOS
            OBS_LTX(:,:,:,IT) = OBS_LTX(:,:,:,IT)*SCLNOS
C
            IF ( NLQMAD.NE.1 ) CMNTTX(1:NLMQMAD,IT)
     &           = CMNTTX(1:NLMQMAD,IT)*SCLNOS
C
         END DO
C
C--------------------------------------------------------------- FULLPOT
         IF ( FULLPOT ) THEN
C
            DO IT = ITBOT,ITTOP
               IM = IMT(IT)
               IRTOP = JRCRI(IM)
C
               DO LM = 1,NLMFPMAX
                  IF ( KLMFP(LM,IT).NE.0 ) THEN
                     DO I1 = 1,3
                        DO IR = 1,IRTOP
                           RHO2NSX(IR,LM,IT,I1) = RHO2NSX(IR,LM,IT,I1)
     &                        *SCLNOS
                        END DO
                     END DO
                  END IF
               END DO
C
            END DO
C
C------------------------------------------------------------------- ASA
         ELSE
C
            DO IT = ITBOT,ITTOP
               IM = IMT(IT)
               IRTOP = JRWS(IM)
C
               DO IR = 1,IRTOP
                  RHOCHRX(IR,IT) = RHOCHRX(IR,IT)*SCLNOS
                  RHOSPNX(IR,IT) = RHOSPNX(IR,IT)*SCLNOS
                  RHOORBX(IR,IT) = RHOORBX(IR,IT)*SCLNOS
               END DO
C
            END DO
C
         END IF
C
      END IF
C
C???????????????????????????????????????????????????????????????????????
C       SCF_CHECK_SPLITSS:  print single site and back scattering
C                           contribution in case of  SPLITSS
C???????????????????????????????????????????????????????????????????????
      IF ( SPLITSS .AND. SCF_CHECK_SPLITSS ) THEN
C
         IW = 6
C
         WRITE (IW,99008) IEPANEL,IEPATH,IECURR
C
         OBS_X(:,:) = 0D0
C
         DO IT = ITBOT,ITTOP
C
            OBS_X(:,:) = OBS_X(:,:) + OBS_TX(:,:,IT)*CONC(IT)*NAT(IT)
C
            WRITE (IW,99002) NETAB(1),ETAB(NETAB(1),1),IT,TXT_T(IT)
C
            WRITE (IW,99003) (OBS_LTX(0,IDOS,IL,IT),OBS_LTX(0,ISMT,IL,IT
     &                       ),OBS_LTX(0,IOMT,IL,IT),
     &                       OBS_LTX(0,IHFF,IL,IT)*1D-3,IL=1,
     &                       MIN(3,NLT(IT)))
            IF ( NLT(IT).GT.3 ) WRITE (IW,99004)
     &                                 (OBS_LTX(0,IDOS,IL,IT),OBS_LTX(0,
     &                                 ISMT,IL,IT),OBS_LTX(0,IOMT,IL,IT)
     &                                 ,OBS_LTX(0,IHFF,IL,IT)*1D-3,IL=4,
     &                                 NLT(IT))
C
            WRITE (IW,99005) OBS_TX(0,IDOS,IT),OBS_TX(0,ISMT,IT),
     &                       OBS_TX(0,IOMT,IT),(OBS_TX(0,IHFF,IT)*1D-3)
            IF ( NT.GT.1 ) WRITE (IW,99006) OBS_TX(0,IBND,IT)
C
            IF ( IT.LT.ITTOP ) THEN
               WRITE (IW,'(1X,79(''-''))')
            ELSE
               WRITE (IW,99007) TOTDOS,OBS_X(0,IDOS),OBS_X(0,ISMT),
     &                          OBS_X(0,IOMT),OBS_X(0,IBND)
C
            END IF
C
         END DO
C
      END IF
C???????????????????????????????????????????????????????????????????????
C
C=======================================================================
C              add contribution along path and reset arrays
C=======================================================================
C
      TOTNOS = TOTNOS + TOTNOSX
C
      DO IT = ITBOT,ITTOP
C
         OBS_T(:,:,IT) = OBS_T(:,:,IT) + OBS_TX(:,:,IT)
         OBS_LT(:,:,:,IT) = OBS_LT(:,:,:,IT) + OBS_LTX(:,:,:,IT)
C
         IF ( NLQMAD.NE.1 ) CMNTT(1:NLMQMAD,IT) = CMNTT(1:NLMQMAD,IT)
     &        + CMNTTX(1:NLMQMAD,IT)
C
      END DO
C
C--------------------------------------------------------------- FULLPOT
      IF ( FULLPOT ) THEN
C
         DO IT = ITBOT,ITTOP
            IM = IMT(IT)
            IRTOP = JRCRI(IM)
C
            IF ( IREL.LE.1 ) THEN
               DO LM = 1,NLMFPMAX
                  IF ( KLMFP(LM,IT).NE.0 ) RHO2NSX(1:IRTOP,LM,IT,1)
     &                 = 2*RHO2NSX(1:IRTOP,LM,IT,1)
               END DO
               NDENS = 1
            ELSE IF ( IREL.EQ.2 ) THEN
               NDENS = 2
            ELSE
               NDENS = 3
            END IF
C
            DO LM = 1,NLMFPMAX
               IF ( KLMFP(LM,IT).NE.0 ) RHO2NS(1:IRTOP,LM,IT,1:NDENS)
     &              = RHO2NS(1:IRTOP,LM,IT,1:NDENS)
     &              + RHO2NSX(1:IRTOP,LM,IT,1:NDENS)
            END DO
C
         END DO
C
         RHO2NSX(:,:,:,:) = 0.0D0
C
C------------------------------------------------------------------- ASA
      ELSE
C
         DO IT = ITBOT,ITTOP
            IM = IMT(IT)
            IRTOP = JRWS(IM)
C
            IF ( IREL.LE.1 ) THEN
               DO IR = 1,IRTOP
                  RHOCHR(IR,IT) = RHOCHR(IR,IT) + 2*RHOCHRX(IR,IT)
                  RHOSPN(IR,IT) = 0.0D0
                  RHOORB(IR,IT) = 0.0D0
               END DO
            ELSE
               DO IR = 1,IRTOP
                  RHOCHR(IR,IT) = RHOCHR(IR,IT) + RHOCHRX(IR,IT)
                  RHOSPN(IR,IT) = RHOSPN(IR,IT) + RHOSPNX(IR,IT)
                  RHOORB(IR,IT) = RHOORB(IR,IT) + RHOORBX(IR,IT)
               END DO
            END IF
C
         END DO
C
         RHOCHRX(:,:) = 0.0D0
         RHOSPNX(:,:) = 0.0D0
         RHOORBX(:,:) = 0.0D0
C
      END IF
C
C-----------------------------------------------------------------------
C
      TOTNOSX = 0D0
C
C--------------------------------- moments, hyperfine field, band energy
C
      OBS_TX(:,:,:) = 0D0
      OBS_LTX(:,:,:,:) = 0D0
      CMNTTX(:,:) = 0D0
C
Cc         IF ( LLOYD ) THEN
Ccc            SCLEBAND = EBANDLD/EBAND
Ccc            SCLEBAND = 1
Cc            SCLEBAND = SCLNOS
Cc            EBAND = 0D0
Cc            DO IT = ITBOT,ITTOP
Cc               EBND_T(IT) = EBND_T(IT)*SCLEBAND
Cc               EBAND = EBAND + EBND_T(IT)*CONC(IT)*NAT(IT)
Cc            END DO
Ccc            IF ( ABS(EBANDLD/EBAND-1D0).GT.1D-8 )
Ccc     &            STOP 'in <CHRDNS_NORM>:  SCLEBAND'
Cc         END IF
C
C-----------------------------------------------------------------------
C
99001 FORMAT (' Lloyd scaling     ',F12.8)
99002 FORMAT (/,I4,' E=',2F7.4,10X,'IT=',I2,2X,A,:,/,15X,
     &        'DOS  [1/Ry]  |  m_spin  [m_B]  |  m_orb   [m_B]  |',
     &        '   B_tot   [kG]',/,' INT(DE) crystal  ',F8.3,10X,F8.3,
     &        10X,F8.3,10X,F8.1,/,' TOTAL   crystal  ',F8.3,10X,F8.3,
     &        10X,F8.3,10X,F8.1)
99003 FORMAT ('         DOS      NOS              m_spin',
     &        '             m_orb    B_val   ',/,'  s ',9X,F9.4,10X,
     &        F9.4,10X,F9.5,F8.2,:,/,'  p ',9X,F9.4,10X,F9.4,10X,F9.5,
     &        F8.2,:,/,'  d ',9X,F9.4,10X,F9.4,10X,F9.5,F8.2)
99004 FORMAT ('  f ',9X,F9.4,10X,F9.4,10X,F9.5,F8.2,:,/,'  g ',9X,F9.4,
     &        10X,F9.4,10X,F9.5,F8.2)
99005 FORMAT (' sum',9X,F9.4,10X,F9.4,10X,F9.5,F8.2)
99006 FORMAT (' E_band',F19.8,' [Ry]')
99007 FORMAT (' ',79('-'),/,' TOT',2F9.4,10X,F9.4,10X,F9.5,/,' E_band',
     &        F19.8,' [Ry]',/,' ',79('='))
99008 FORMAT (/,' <CHRDNS_NORM>:  ',/,10X,
     &        'type resolved single site results for',/,10X,'IEPANEL =',
     &        I3,3X,'IEPATH =',I3,3X,'(1 = SS, 2 = BS)',3X,'IECURR =',
     &        I3)
      END
