C*==rmeshasa.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE RMESHASA
C   ********************************************************************
C   *                                                                  *
C   *   generate the radial meshes for       ASA                       *
C   *                                                                  *
C   *   DX or BRMSH will be reset                                      *
C   *   as the formatted input may be inaccurate                       *
C   *                                                                  *
C   *   RMESHTYPE = EXPONENTIAL                                        *
C   *                                                                  *
C   *   fix mesh using  R(1), RWS and JRWS - adjust DX if necessary    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:NM,NMMAX,DX,R,DRDI,R2DRDI,RMESHTYPE,RWS,RMT,
     &    JRWS,JRMT,JRCUT,NPAN,NRMAX,ARMSH,BRMSH,DRDIOVR,ARMSH_ASA,
     &    BRMSH_ASA,DX_ASA,EXPDX_ASA
      IMPLICIT NONE
C*--RMESHASA20
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      REAL*8 A1,B1,DXCALC,EA,EXPDX,RWSTST
      INTEGER IM,IR
C
C*** End of declarations rewritten by SPAG
C
      IF ( .NOT.ALLOCATED(DX_ASA) ) THEN
         ALLOCATE (ARMSH_ASA(NMMAX),BRMSH_ASA(NMMAX))
         ALLOCATE (DX_ASA(NMMAX),EXPDX_ASA(NMMAX))
      END IF
C
C=======================================================================
      IF ( RMESHTYPE.EQ.'EXPONENTIAL ' ) THEN
C
         DO IM = 1,NM
C
            DXCALC = LOG(RWS(IM)/R(1,IM))/DBLE(JRWS(IM)-1)
            IF ( ABS(DX(IM)-DXCALC).GT.1D-7 ) WRITE (6,99002) IM,DX(IM),
     &           R(1,IM),DXCALC
C
            EXPDX = DEXP(DX(IM))
C
            DO IR = 2,NRMAX
               R(IR,IM) = R(IR-1,IM)*EXPDX
            END DO
            DO IR = 1,NRMAX
               DRDI(IR,IM) = R(IR,IM)*DX(IM)
               R2DRDI(IR,IM) = R(IR,IM)*R(IR,IM)*R(IR,IM)*DX(IM)
               DRDIOVR(IR,IM) = DRDI(IR,IM)/R(IR,IM)
            END DO
C
            DX_ASA(IM) = DX(IM)
            EXPDX_ASA(IM) = EXPDX
         END DO
C
C-----------------------------------------------------------------------
C                          print mesh info
C-----------------------------------------------------------------------
         WRITE (6,99003) '       DX'
         DO IM = 1,NM
            IR = JRMT(IM)
            WRITE (6,99004) IM,IR,R(IR,IM),RMT(IM),JRWS(IM),RWS(IM),
     &                      DX(IM)
         END DO
C
C=======================================================================
      ELSE IF ( RMESHTYPE.EQ.'JUELICH     ' ) THEN
         DO IM = 1,NM
C
            RWSTST = BRMSH(IM)*(EXP(ARMSH(IM)*JRWS(IM))-1.0D0)
            IF ( ABS(RWSTST-RWS(IM)).GT.1D-6 ) THEN
               WRITE (6,99001)
               WRITE (6,*) 'input data inconsistent'
               WRITE (6,*) ' IM               ',IM
               WRITE (6,*) ' RWS     JRWS     ',RWS(IM),JRWS(IM)
               WRITE (6,*) ' RWS(JRWS,A,B)    ',RWSTST
            END IF
C
C drop first mesh point r=0 of juelich format
C
            BRMSH(IM) = RWS(IM)/(EXP(ARMSH(IM)*JRWS(IM))-1.0D0)
C
            A1 = ARMSH(IM)
            B1 = BRMSH(IM)
            ARMSH_ASA(IM) = ARMSH(IM)
            BRMSH_ASA(IM) = BRMSH(IM)
C
            DO IR = 1,JRWS(IM)
               EA = EXP(A1*DBLE(IR))
               R(IR,IM) = B1*(EA-1.0D0)
               DRDI(IR,IM) = A1*B1*EA
               R2DRDI(IR,IM) = R(IR,IM)*R(IR,IM)*DRDI(IR,IM)
               DRDIOVR(IR,IM) = DRDI(IR,IM)/R(IR,IM)
            END DO
C
         END DO
C
C-----------------------------------------------------------------------
C                          print mesh info
C-----------------------------------------------------------------------
         WRITE (6,99003) '     ARMSH      BRMSH'
         DO IM = 1,NM
            IR = JRMT(IM)
            WRITE (6,99005) IM,IR,R(IR,IM),RMT(IM),JRWS(IM),RWS(IM),
     &                      ARMSH(IM),BRMSH(IM)
         END DO
C
C=======================================================================
      ELSE
         STOP 'in <RMESHASA>: check RMESHTYPE in POTFILE '
      END IF
C
C=======================================================================
C       ensure settings according to ASA for integration routines
C=======================================================================
C
      DO IM = 1,NM
         NPAN(IM) = 1
         JRCUT(0,IM) = 0
         JRCUT(1,IM) = JRWS(IM)
      END DO
C
C=======================================================================
C     weights for radial integration on the basis of the simpson rule
C=======================================================================
C
      CALL CORE_RMESH
C
      CALL GET_W_RADINT
C
      RETURN
C
99001 FORMAT (/,' ##### TROUBLE in <RMESHASA> ',51('#'))
99002 FORMAT (/,60('*'),/,10X,'WARNING from <RMESHASA>',/,10X,
     &        ' MESH     ',I15,/,10X,' DX read  ',2F15.10,/,10X,
     &        ' DX calc  ',F15.10,/,60('*'),/)
99003 FORMAT (/,10X,'radial mesh parameters:',//,10X,
     &        'IM    JRMT   R(JRMT)    RMT      JRWS    RWS',A)
99004 FORMAT (8X,I4,I7,2F10.5,I7,F10.5,2F14.9)
99005 FORMAT (8X,I4,I7,2F10.5,I7,F10.5,2F11.8)
      END
C*==rmeshfp.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE RMESHFP(NRPAN,XRSF,DRSF,SCLM)
C   ********************************************************************
C   *                                                                  *
C   *      generate the radial meshes for       FULLPOT                *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LATTICE,ONLY:ALAT
      USE MOD_RMESH,ONLY:NM,DX,R,DRDI,R2DRDI,RMESHTYPE,RWS,RMT,JRWS,
     &    JRMT,NRMAX,ARMSH,BRMSH,WINTLM,DRDIOVR,JRCRI,JRCUT,JRNS1,FLMSF,
     &    LMISF,NRSFTOT,NSF,NRNS,NPAN,NMMAX,NRSFMAX,JRNSMIN,NPANMAX,
     &    R2DRDI_W_RADINT
      USE MOD_ANGMOM,ONLY:L_LM,M_LM
      IMPLICIT NONE
C*--RMESHFP169
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 DRSF(NRSFMAX,NMMAX),SCLM(NMMAX),XRSF(NRSFMAX,NMMAX)
      INTEGER NRPAN(NPANMAX,NMMAX)
C
C Local variables
C
      REAL*8 A1,B1,EA,EXPDX,RMTTST
      INTEGER IM,IPAN,IR,IROFF,IRSF,ISF,J,JNMT
C
C*** End of declarations rewritten by SPAG
C
C=======================================================================
      IF ( RMESHTYPE.EQ.'EXPONENTIAL ' ) THEN
C
         DO IM = 1,NM
C
            EXPDX = DEXP(DX(IM))
C
            DO IR = 2,JRMT(IM)
               R(IR,IM) = R(IR-1,IM)*EXPDX
            END DO
            DO IR = 1,JRMT(IM)
               DRDI(IR,IM) = R(IR,IM)*DX(IM)
               R2DRDI(IR,IM) = R(IR,IM)*R(IR,IM)*DRDI(IR,IM)
               DRDIOVR(IR,IM) = DRDI(IR,IM)/R(IR,IM)
            END DO
         END DO
C
C=======================================================================
      ELSE IF ( RMESHTYPE.EQ.'JUELICH     ' ) THEN
C
C drop first mesh point r=0 of original juelich format
C
         DO IM = 1,NM
C
            RMTTST = BRMSH(IM)*(EXP(ARMSH(IM)*JRMT(IM))-1.0D0)
            IF ( ABS(RMTTST-RMT(IM)).GT.1D-6 ) THEN
               WRITE (6,99001)
               WRITE (6,*) 'input data inconsistent'
               WRITE (6,*) ' IM               ',IM
               WRITE (6,*) ' RMT     JRMT     ',RMT(IM),JRMT(IM)
               WRITE (6,*) ' RMT(JRMT,A,B)    ',RMTTST
            END IF
C
            BRMSH(IM) = RMT(IM)/(EXP(ARMSH(IM)*JRMT(IM))-1.0D0)
C
            A1 = ARMSH(IM)
            B1 = BRMSH(IM)
C
            DO J = 1,JRMT(IM)
               EA = EXP(A1*DBLE(J))
               R(J,IM) = B1*(EA-1.0D0)
               DRDI(J,IM) = A1*B1*EA
               R2DRDI(J,IM) = R(J,IM)*R(J,IM)*DRDI(J,IM)
               DRDIOVR(IR,IM) = DRDI(IR,IM)/R(IR,IM)
            END DO
         END DO
C
C=======================================================================
      ELSE
         STOP 'in <RMESHFP>: check RMESHTYPE in POTFILE '
      END IF
C
C-----------------------------------------------------------------------
C     fill cell-type depending mesh points in the non-muffin-tin-region
C-----------------------------------------------------------------------
C
      DO IM = 1,NM
         DO JNMT = 1,NRSFTOT(IM)
            J = JNMT + JRMT(IM)
            R(J,IM) = SCLM(IM)*ALAT*XRSF(JNMT,IM)
            DRDI(J,IM) = SCLM(IM)*ALAT*DRSF(JNMT,IM)
            R2DRDI(J,IM) = R(J,IM)*R(J,IM)*DRDI(J,IM)
         END DO
C
         JRCUT(0,IM) = 0
         JRCUT(1,IM) = JRMT(IM)
         DO IPAN = 2,NPAN(IM)
            JRCUT(IPAN,IM) = JRCUT(IPAN-1,IM) + NRPAN(IPAN,IM)
         END DO
         JRCRI(IM) = JRCUT(NPAN(IM),IM)
         IF ( JRCRI(IM).GT.NRMAX ) THEN
            WRITE (6,99001)
            WRITE (6,99006) 'MESH    ',IM
            WRITE (6,99006) 'JRCRI   ',JRCRI(IM),'  > NRMAX=',NRMAX
         END IF
C
         JRWS(IM) = JRCRI(IM)
         RWS(IM) = R(JRWS(IM),IM)
         NRNS(IM) = JRCUT(NPAN(IM),IM) - JRNS1(IM)
      END DO
C
      WRITE (6,99002)
      DO IM = 1,NM
         WRITE (6,99005) ' '
         WRITE (6,99005) 'mesh',IM
         IF ( RMESHTYPE.EQ.'EXPONENTIAL ' ) THEN
            WRITE (6,99008) DX(IM)
         ELSE
            WRITE (6,99009) ARMSH(IM),BRMSH(IM)
         END IF
         WRITE (6,99003) 'JRNSMIN ',JRNSMIN,R(JRNSMIN,IM)
         WRITE (6,99003) 'JRNS1   ',JRNS1(IM),R(JRNS1(IM),IM)
         WRITE (6,99004) 'JRMT    ',JRMT(IM),R(JRMT(IM),IM),RMT(IM)
         DO IPAN = 0,NPAN(IM)
            J = JRCUT(IPAN,IM)
            IF ( IPAN.EQ.0 ) THEN
               WRITE (6,99005) 'IPAN',IPAN,J
            ELSE IF ( IPAN.LT.NPAN(IM) ) THEN
               WRITE (6,99005) 'IPAN',IPAN,J,R(J,IM),R(J+1,IM)
            ELSE
               WRITE (6,99005) 'IPAN',IPAN,J,R(J,IM)
            END IF
         END DO
         WRITE (6,99004) 'JRWS    ',JRWS(IM),R(JRWS(IM),IM),RWS(IM)
         WRITE (6,99004) 'JRCRI   ',JRCRI(IM),R(JRCRI(IM),IM)
         WRITE (6,99004) 'NRNS    ',NRNS(IM)
         WRITE (6,99005) ' '
         WRITE (6,99005) 'shape functions:'
         WRITE (6,99006) 'NSF ',NSF(IM),'NRSFTOT ',NRSFTOT(IM)
         WRITE (6,99007) (LMISF(ISF,IM),L_LM(LMISF(ISF,IM)),
     &                   M_LM(LMISF(ISF,IM)),ISF=1,NSF(IM))
      END DO
C
C
      CALL CORE_RMESH
C
      CALL GET_W_RADINT
C
C=======================================================================
C                   set up integration weights  WINTLM
C=======================================================================
C
      DO IM = 1,NM
C
         IROFF = JRCUT(1,IM)
C
         DO ISF = 1,NSF(IM)
C
            DO IRSF = 1,NRSFTOT(IM)
               IR = IRSF + IROFF
               WINTLM(IRSF,ISF,IM) = FLMSF(IRSF,ISF,IM)
     &                               *R2DRDI_W_RADINT(IR,IM)
            END DO
C
         END DO
C
      END DO
C
      RETURN
C
99001 FORMAT (/,' ##### TROUBLE in <RMESHFP>  ',51('#'))
99002 FORMAT (/,10X,'radial mesh parameters:')
99003 FORMAT (10X,12X,4(A,I5,:,5X,'R =',F12.8,5X))
99004 FORMAT (10X,12X,A,I5,5X,:,'R =',F12.8,5X,F12.8)
99005 FORMAT (10X,A,I4,:,4X,'JRCUT   ',I5,:,5X,'R =',F12.8,5X,F12.8)
99006 FORMAT (10X,A,I4,4X,A,I5)
99007 FORMAT (10X,'LM  ',5(I4,' (',I2,',',I3,')'),:,/,
     &        (14X,5(I4,' (',I2,',',I3,')')))
99008 FORMAT (39X,'DX =',F12.8)
99009 FORMAT (40X,'A =',F12.8,/,40X,'B =',F12.8)
      END
C*==calc_vol_ac.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE CALC_VOL_AC
C   ********************************************************************
C   *                                                                  *
C   *  set the volume for the atomic cells                             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:NRMAX,NMMAX,NPAN,JRCUT,JRMT,R2DRDI,RMT,JRWS,
     &    RWS,NM,R,SPHERCELL,VOL_IS_ASA_NUM,VOL_IS_ASA_ANA,
     &    VOL_IS_FP_NUM,VOL_AC,VOL_AC_ASA,VOL_AC_FP,FULLPOT,FLMSF,
     &    R2DRDI_W_RADINT
      USE MOD_LATTICE,ONLY:VOLUC,ALAT
      USE MOD_CONSTANTS,ONLY:PI,SQRT_4PI
      USE MOD_SITES,ONLY:IMQ,NQHOST
      IMPLICIT NONE
C*--CALC_VOL_AC363
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CALC_VOL_AC')
C
C Local variables
C
      REAL*8 DDOT,YLAG
      INTEGER IM,IQ,IR1,IRCRIT,NR
      REAL*8 SUM1,U1(NRMAX),U1MT
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (VOL_AC(NMMAX))
      IF ( .NOT.ALLOCATED(VOL_AC_ASA) ) THEN
         ALLOCATE (VOL_AC_ASA(NMMAX),VOL_AC_FP(NMMAX))
         ALLOCATE (VOL_IS_FP_NUM(NMMAX))
         ALLOCATE (VOL_IS_ASA_NUM(NMMAX),VOL_IS_ASA_ANA(NMMAX))
      END IF
C
      VOL_AC(:) = 0D0
      VOL_AC_FP(:) = 0D0
      VOL_AC_ASA(:) = 0D0
      VOL_IS_FP_NUM(:) = 0D0
      VOL_IS_ASA_NUM(:) = 0D0
      VOL_IS_ASA_ANA(:) = 0D0
C
      WRITE (6,99006) ROUTINE(1:LEN_TRIM(ROUTINE))
C
      IF ( .NOT.FULLPOT ) THEN
C
C-----------------------------------------------------------------------
C         calculate cell-related interstitial volume for  ASA
C-----------------------------------------------------------------------
C
         DO IM = 1,NM
C
            CALL RRADINT_R(IM,R2DRDI(1,IM),U1)
C
            U1MT = YLAG(RMT(IM),R(1,IM),U1,0,3,JRWS(IM))
C
            VOL_IS_ASA_NUM(IM) = 4D0*PI*(U1(JRWS(IM))-U1MT)
            VOL_IS_ASA_ANA(IM) = 4D0*PI*(RWS(IM)**3-RMT(IM)**3)/3D0
C
            VOL_AC_ASA(IM) = 4D0*PI*RWS(IM)**3/3D0
C
            IF ( ABS(1D0-VOL_IS_ASA_NUM(IM)/VOL_IS_ASA_ANA(IM))
     &           .GT.1D-5 ) WRITE (6,*) ROUTINE(1:LEN_TRIM(ROUTINE)),IM,
     &                                  VOL_IS_ASA_NUM(IM),
     &                                  VOL_IS_ASA_ANA(IM)
C
            WRITE (6,99003) IM,VOL_IS_ASA_NUM(IM),VOL_AC_ASA(IM)
            VOL_AC(IM) = VOL_AC_ASA(IM)
C
         END DO
C
      ELSE
C
C-----------------------------------------------------------------------
C         calculate cell-related interstitial volume for  FULLPOT
C-----------------------------------------------------------------------
C
         WRITE (6,99002)
C
         DO IM = 1,NM
C
            IRCRIT = JRCUT(NPAN(IM),IM)
C
            IF ( SPHERCELL ) THEN
C
               VOL_AC_FP(IM) = 4D0*PI*R(IRCRIT,IM)**3/3D0
C
            ELSE
C
               NR = IRCRIT - JRMT(IM)
               IR1 = JRMT(IM) + 1
C
               VOL_IS_FP_NUM(IM) = SQRT_4PI*DDOT(NR,FLMSF(1,1,IM),1,
     &                             R2DRDI_W_RADINT(IR1,IM),1)
C
               VOL_AC_FP(IM) = 4D0*PI*RMT(IM)**3/3D0 + VOL_IS_FP_NUM(IM)
C
            END IF
C
            WRITE (6,99003) IM,VOL_IS_FP_NUM(IM),VOL_AC_FP(IM)
C
            VOL_AC(IM) = VOL_AC_FP(IM)
C
         END DO
C
      END IF
C
C ======================================================================
C                     check consistency of results
C ======================================================================
C
      SUM1 = 0D0
      DO IQ = 1,NQHOST
         IM = IMQ(IQ)
         SUM1 = SUM1 + VOL_AC(IM)
      END DO
C
      WRITE (6,99004) 'SUM  ',SUM1
      WRITE (6,99005) 'VOLUC',VOLUC*ALAT**3
      IF ( ABS(1D0-VOLUC*ALAT**3/SUM1).GT.1D-5 ) THEN
         WRITE (6,99001) ROUTINE(1:LEN_TRIM(ROUTINE))
         STOP
      END IF
C
C ======================================================================
99001 FORMAT (/,10X,'<',A,'>: volume of unit cell not properly',/,10X,
     &        'reproduced by integration over shape functions',/)
99002 FORMAT (/,10X,'atomic interstitial and cell volumes '//,10X,'IM',
     &        8X,'IS (FP)         AC (FP)')
99003 FORMAT (I12,3X,2F16.8)
99004 FORMAT (10X,A,16X,F16.8)
99005 FORMAT (10X,A,16X,F16.8,/)
99006 FORMAT (/,1X,79('*'),/,34X,'<',A,'>',/,14X,
     &        'setting the volume VOL_AC for the atomic cells',/,1X,
     &        79('*'),/)
C
      END
C*==get_w_radint.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE GET_W_RADINT
C   ********************************************************************
C   *                                                                  *
C   *  call  GET_INTEGRATION_WEIGHTS  for each mesh IM to get          *
C   *  W_RADINT(IR,IM)  for the coded integration scheme  (SIMPSON)    *
C   *                                                                  *
C   *  set:     DRDI_W_RADINT(IR,IM) =   DRDI(IR,IM)*W_RADINT(IR,IM)   *
C   *         R2DRDI_W_RADINT(IR,IM) = R2DRDI(IR,IM)*W_RADINT(IR,IM)   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:NPAN,NRMAX,NMMAX,JRCUT,NM,W_RADINT,
     &    DRDI_W_RADINT,R2DRDI_W_RADINT,R2DRDI,DRDI
      IMPLICIT NONE
C*--GET_W_RADINT517
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      INTEGER IM,IR
C
C*** End of declarations rewritten by SPAG
C
      IF ( ALLOCATED(W_RADINT) ) THEN
         DEALLOCATE (W_RADINT,DRDI_W_RADINT)
         DEALLOCATE (R2DRDI_W_RADINT)
      END IF
C
      ALLOCATE (W_RADINT(NRMAX,NMMAX))
      ALLOCATE (DRDI_W_RADINT(NRMAX,NMMAX),R2DRDI_W_RADINT(NRMAX,NMMAX))
C
      W_RADINT(:,:) = 0D0
      DRDI_W_RADINT(:,:) = 0D0
      R2DRDI_W_RADINT(:,:) = 0D0
C
C-----------------------------------------------------------------------
C     loop over radial meshes
C-----------------------------------------------------------------------
C
      LOOP_IM:DO IM = 1,NM
C
         CALL GET_INTEGRATION_WEIGHTS(NPAN(IM),JRCUT(0,IM),NRMAX,
     &                                W_RADINT(1,IM))
C
         DO IR = 1,JRCUT(NPAN(IM),IM)
            DRDI_W_RADINT(IR,IM) = DRDI(IR,IM)*W_RADINT(IR,IM)
            R2DRDI_W_RADINT(IR,IM) = R2DRDI(IR,IM)*W_RADINT(IR,IM)
         END DO
C
      END DO LOOP_IM
C
      END
C*==get_integration_weights.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE GET_INTEGRATION_WEIGHTS(NPAN,JRCUT,NR,W)
C   ********************************************************************
C   *                                                                  *
C   * this subroutine supplies weights w(i) for an radial integration  *
C   *                                                                  *
C   *         r_max                                                    *
C   *     I = S      f(r') dr'   =  Sum_{i=1}^{i_max} w(i) * F(i)      *
C   *         0                                                        *
C   *                                                                  *
C   * with F(i) = f(r) * dr/di                                         *
C   *                                                                  *
C   * the variables NPAN and JRCUT specify the various panels          *
C   * of the integration regime.                                       *
C   *                                                                  *
C   * ASA: NPAN = 1, JRCUT(0)=0, JRCUT(1)=JRWS                         *
C   *                                                                  *
C   * HERE: the standard Simpson scheme has been implemented           *
C   * implying within a panel the weights  1 4 2 4 2 .... 2 4 1        *
C   * in case of an EVEN number of radial mesh points in a panel       *
C   * the first interval will be treated separatly via Simpson         *
C   *                                                                  *
C   * MAKE SURE THAT THE SAME INTEGRATION SCHEME IS USED EVERYWHERE    *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--GET_INTEGRATION_WEIGHTS591
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='GET_INTEGRATION_WEIGHTS')
      REAL*8 R12,R3
      PARAMETER (R12=1D0/12D0,R3=1D0/3D0)
C
C Dummy arguments
C
      INTEGER NPAN,NR
      INTEGER JRCUT(0:NPAN)
      REAL*8 W(NR)
C
C Local variables
C
      INTEGER IPAN,IR,IR_END,IR_START,IR_START0
      REAL*8 WSUM
C
C*** End of declarations rewritten by SPAG
C
      W(1:NR) = 0D0
C
C-----------------------------------------------------------------------
C     loop over panels
C-----------------------------------------------------------------------
C
      DO IPAN = 1,NPAN
C
         IR_END = JRCUT(IPAN)
         IR_START = JRCUT(IPAN-1) + 1
         IR_START0 = IR_START
         IF ( IR_START+2.GT.IR_END ) CALL STOP_MESSAGE(ROUTINE,
     &        'at least 3 mesh points in panel required')
C
         IF ( MOD((IR_END-IR_START),2).NE.0 ) THEN
            IR = IR_START
            W(IR+0) = W(IR+0) + 5*R12
            W(IR+1) = W(IR+1) + 8*R12
            W(IR+2) = W(IR+2) - 1*R12
            IR_START = IR_START + 1
         END IF
C
         DO IR = IR_START,IR_END - 2,2
            W(IR+0) = W(IR+0) + 1*R3
            W(IR+1) = W(IR+1) + 4*R3
            W(IR+2) = W(IR+2) + 1*R3
         END DO
C
         WSUM = SUM(W(IR_START0:IR_END))
         IF ( ABS(1D0-DFLOAT(IR_END-IR_START0)/WSUM).GT.1D-14 ) THEN
            WRITE (6,99001) IPAN,WSUM,IR_END - IR_START0
            CALL STOP_MESSAGE(ROUTINE,'WSUM != IR_END-IR_START0')
         END IF
C
      END DO
C
99001 FORMAT (/,5('#'),' IPAN      ',I4,/,5('#'),' SUM(w)    ',F20.14,/,
     &        5('#'),' Delta I   ',I5,/)
      END
