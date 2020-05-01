C*==fpasastart.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPASASTART
C   ********************************************************************
C   *                                                                  *
C   *  start FP type calculations using an ASA potential               *
C   *                                                                  *
C   * 08/07/14  VMTZ: sign convention as in ASA                        *
C   ********************************************************************
C
      USE MOD_ENERGY,ONLY:EFERMI
      USE MOD_CALCMODE,ONLY:ORBPOL,KMROT,IREL
      USE MOD_RMESH,ONLY:RMESHTYPE,NRMAX,JRNSMIN,NRSFMAX,NLMSFMAX,
     &    NSFMAX,NPANMAX,NMMAX,LMISF,NSF,KLMSF,ISFLM,FLMSF,NPAN,NRSFTOT,
     &    JRCRI,JRCUT,JRMT,JRNS1,R2DRDI,RMT,JRWS,RWS,BRMSH,ARMSH,DX,NM,
     &    R,SPHERCELL,VOL_IS_ASA_NUM,VOL_IS_ASA_ANA,VOL_IS_FP_NUM,
     &    VOL_AC_ASA,VOL_AC_FP,FULLPOT,R2DRDI_W_RADINT
      USE MOD_TYPES,ONLY:AOPT,ITBOT,ITTOP,NLMFPMAX,NTMAX,NAT,CONC,
     &    LTXT_T,TXT_T,BNST,VNST,BT,VT,LMIFP,NFPT,NLMFPT,IMT,Z,KLMFP
      USE MOD_ANGMOM,ONLY:L_LM,M_LM
      USE MOD_LATTICE,ONLY:VOLUC,ALAT,SYSTEM_DIMENSION,SYSTEM_TYPE
      USE MOD_CONSTANTS,ONLY:PI,SQRT_4PI
      USE MOD_SITES,ONLY:IMQ,NQHOST,NQ_L,NQ_I,NQ_R,NOQ,ITOQ
      IMPLICIT NONE
C*--FPASASTART24
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='FPASASTART')
      REAL*8 RNS1
      PARAMETER (RNS1=0.10D0)
C
C Local variables
C
      REAL*8 AUX,BWSASA,DELR,DRSF(:,:),F00SF,RASA(NRMAX,NMMAX),RAT,
     &       RINTVAV,RINTVOL,RWSASA,SCLM(:),SUM1,SUM2,SUM3,SUM4,
     &       U1(NRMAX),U1MT,U1_IR,U2(NRMAX),VAV,VMTZ_ASA,VOL,VWSASA,X,
     &       XINC,XMT,XRSF(:,:)
      INTEGER IFP,IM,IO,IPAN,IQ,IR,IRCRIT,IRMTIN,IRSF,IS,ISF,IT,
     &        ITBOT_EXTENDED,ITTOP_EXTENDED,J,JRWSASA(NMMAX),LM,NPANEL,
     &        NRPAN(:,:)
      REAL*8 YLAG
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DRSF,SCLM,XRSF,NRPAN
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (VOL_AC_ASA(NMMAX),VOL_AC_FP(NMMAX))
      ALLOCATE (VOL_IS_FP_NUM(NMMAX))
      ALLOCATE (VOL_IS_ASA_NUM(NMMAX),VOL_IS_ASA_ANA(NMMAX))
      ALLOCATE (DRSF(NRSFMAX,NMMAX),SCLM(NMMAX))
      ALLOCATE (XRSF(NRSFMAX,NMMAX),NRPAN(NPANMAX,NMMAX))
C
      WRITE (6,99006) ROUTINE(1:LEN_TRIM(ROUTINE))
C
C-----------------------------------------------------------------------
C     for  2D LIR- (bulk) and  LIV (surface) calculations
C     perform the SAME transformation of the potential for
C     the atoms in the interaction zone as for those of the left bulk
C-----------------------------------------------------------------------
C
      IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' .AND. 
     &     (SYSTEM_TYPE(1:3).EQ.'LIR' .OR. SYSTEM_TYPE(1:3).EQ.'LIV') )
     &     THEN
C
         ITBOT_EXTENDED = ITBOT
         ITTOP_EXTENDED = ITTOP
         DO IQ = NQ_L + 1,NQ_L + NQ_I
            DO IO = 1,NOQ(IQ)
               IT = ITOQ(IO,IQ)
               ITBOT_EXTENDED = MIN(ITBOT_EXTENDED,IT)
               ITTOP_EXTENDED = MAX(ITTOP_EXTENDED,IT)
            END DO
         END DO
C
      ELSE
C
         ITBOT_EXTENDED = ITBOT
         ITTOP_EXTENDED = ITTOP
C
      END IF
C
C-----------------------------------------------------------------------
C                         save  ASA  mesh info
C-----------------------------------------------------------------------
C
      DO IM = 1,NM
C
         RASA(:,IM) = R(:,IM)
         JRWSASA(IM) = JRWS(IM)
C
      END DO
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
         IF ( ABS(1D0-VOL_IS_ASA_NUM(IM)/VOL_IS_ASA_ANA(IM)).GT.1D-5 )
     &        WRITE (6,*) '<FPASASTART>: ',IM,VOL_IS_ASA_NUM(IM),
     &                    VOL_IS_ASA_ANA(IM)
      END DO
C
C***********************************************************************
C***********************************************************************
C            from now on ALL standard radial mesh parameters
C                 refer to the FULL POTENTIAL mesh
C***********************************************************************
C***********************************************************************
C
      FULLPOT = .TRUE.
C
C-----------------------------------------------------------------------
C                         spherical cells
C-----------------------------------------------------------------------
      IF ( SPHERCELL ) THEN
C
         F00SF = SQRT_4PI
C
         DO IM = 1,NM
C
            NPAN(IM) = 2
            NRPAN(2,IM) = 60
            NRSFTOT(IM) = NRPAN(2,IM)
            NSF(IM) = 1
C
            DO LM = 1,NLMSFMAX
               KLMSF(LM,IM) = 0
            END DO
C
            LMISF(1,IM) = 1
            KLMSF(1,IM) = 1
            ISFLM(1,IM) = 1
            DO J = 1,NRSFTOT(IM)
               FLMSF(J,1,IM) = F00SF
            END DO
C
            SCLM(IM) = 1D0
            DELR = (RWS(IM)-RMT(IM))/DBLE(NRSFTOT(IM)-1)
            DO J = 1,NRSFTOT(IM)
               XRSF(J,IM) = (RMT(IM)+DELR*(J-1))/ALAT
               DRSF(J,IM) = DELR/ALAT
            END DO
         END DO
C
      ELSE
C
C-----------------------------------------------------------------------
C                POLYHEDRAL cells  -  read SHAPE FUNCTIONS
C-----------------------------------------------------------------------
C
         CALL SFNREAD(NM,NRPAN,XRSF,DRSF,SCLM)
C
C-----------------------------------------------------------------------
C       rotate the shape functions
C         to the local frames of reference if necessary
C-----------------------------------------------------------------------
C
         IF ( KMROT.NE.0 .AND. IREL.EQ.3 ) CALL FP_SFN_LOCAL
C
      END IF
C
C-----------------------------------------------------------------------
C              generate new FULL POTENTIAL radial mesh
C-----------------------------------------------------------------------
C
      DO IM = 1,NM
C
         RMT(IM) = SCLM(IM)*ALAT*XRSF(1,IM)
C
      END DO
C
      IF ( RMESHTYPE.EQ.'EXPONENTIAL ' ) THEN
C
         DO IM = 1,NM
C
            JRMT(IM) = INT(LOG(RMT(IM)/R(1,IM))/DX(IM)) + 2
            JRMT(IM) = MIN(JRMT(IM),NRMAX-NRSFTOT(IM)-5)
            DX(IM) = LOG(RMT(IM)/R(1,IM))/DBLE(JRMT(IM)-1)
C
            JRNS1(IM) = INT(LOG(RNS1/R(1,IM))/DX(IM)) + 1
C
C------------------------------------------ recover mesh parameter JRCUT
            JRCUT(0,IM) = 0
            JRCUT(1,IM) = JRMT(IM)
            DO IPAN = 2,NPAN(IM)
               JRCUT(IPAN,IM) = JRCUT(IPAN-1,IM) + NRPAN(IPAN,IM)
            END DO
C
         END DO
C
      ELSE IF ( RMESHTYPE.EQ.'JUELICH     ' ) THEN
C
         DO IM = 1,NM
C
            JRMT(IM) = MIN(350,NRMAX-NRSFTOT(IM)-5)
C
C --------------------------- fix ARMSH and BRMSH using Newton algorithm
            RAT = RMT(IM)/R(1,IM)
            X = 1.05D0
 20         CONTINUE
            XMT = X**JRMT(IM)
            XINC = -(XMT-1D0-RAT*(X-1D0))/(JRMT(IM)*XMT/X-RAT)
            IF ( ABS(XINC/X).GT.1D-8 ) THEN
               X = X + XINC
               GOTO 20
            END IF
C
            BRMSH(IM) = RMT(IM)/(XMT-1D0)
            ARMSH(IM) = LOG(XMT)/DBLE(JRMT(IM))
C
            JRNS1(IM) = INT(LOG(RNS1/BRMSH(IM)+1D0)/ARMSH(IM)) + 1
C
         END DO
C
      END IF
C
      DO IM = 1,NM
         JRNS1(IM) = MAX(JRNSMIN,JRNS1(IM))
         JRNS1(IM) = MIN(JRMT(IM),JRNS1(IM))
      END DO
C
      CALL RMESHFP(NRPAN,XRSF,DRSF,SCLM)
C
C-----------------------------------------------------------------------
C         calculate cell-related interstitial volume for  FULLPOT
C-----------------------------------------------------------------------
C
      WRITE (6,99002)
C
      DO IM = 1,NM
C
         NPANEL = NPAN(IM)
         IRMTIN = JRMT(IM)
         IRCRIT = JRCUT(NPANEL,IM)
C
         AUX = 0D0
         DO IR = IRMTIN + 1,IRCRIT
            IRSF = IR - IRMTIN
            AUX = AUX + FLMSF(IRSF,1,IM)*R2DRDI_W_RADINT(IR,IM)
         END DO
         VOL_IS_FP_NUM(IM) = AUX*SQRT_4PI
C
         VOL_AC_FP(IM) = 4D0*PI*RMT(IM)**3/3D0 + VOL_IS_FP_NUM(IM)
C
         WRITE (6,99003) IM,VOL_IS_ASA_NUM(IM),VOL_AC_ASA(IM),
     &                   VOL_IS_FP_NUM(IM),VOL_AC_FP(IM)
      END DO
C
      SUM1 = 0D0
      SUM2 = 0D0
      SUM3 = 0D0
      SUM4 = 0D0
      DO IQ = 1,NQHOST
         IM = IMQ(IQ)
C
         SUM1 = SUM1 + VOL_IS_ASA_NUM(IM)
         SUM2 = SUM2 + VOL_AC_ASA(IM)
C
         SUM3 = SUM3 + VOL_IS_FP_NUM(IM)
         SUM4 = SUM4 + VOL_AC_FP(IM)
      END DO
      WRITE (6,99004) 'SUM  ',SUM1,SUM2,SUM3,SUM4
      WRITE (6,99005) 'VOLUC',VOLUC*ALAT**3,VOLUC*ALAT**3
      IF ( ABS(1D0-VOLUC*ALAT**3/SUM4).GT.1D-5 ) THEN
         WRITE (6,99001) ROUTINE(1:LEN_TRIM(ROUTINE))
         STOP
      END IF
C
C=======================================================================
C                      initialize potentials
C=======================================================================
C
      DO IT = ITBOT_EXTENDED,ITTOP_EXTENDED
C
         IM = IMT(IT)
C
         VNST(:,:,IT) = 0.0D0
         BNST(:,:,IT) = 0.0D0
C
C-----------------------------------------------------------------------
C        interpolate potential functions V and B from ASA to FP mesh
C-----------------------------------------------------------------------
C
         DO IR = 1,JRWSASA(IM)
            U1(IR) = VT(IR,IT) + 2D0*Z(IT)/RASA(IR,IM)
            U2(IR) = BT(IR,IT)
         END DO
C
         IR = JRWSASA(IM)
         RWSASA = RASA(IR,IM)
         VWSASA = VT(IR,IT)
         BWSASA = BT(IR,IT)
         DO IR = 1,JRCRI(IM)
            IF ( R(IR,IM).LE.RWSASA ) THEN
               VT(IR,IT) = YLAG(R(IR,IM),RASA(1,IM),U1,0,3,JRWSASA(IM))
     &                     - 2D0*Z(IT)/R(IR,IM)
               BT(IR,IT) = YLAG(R(IR,IM),RASA(1,IM),U2,0,3,JRWSASA(IM))
            ELSE
               VT(IR,IT) = VWSASA
               BT(IR,IT) = BWSASA
            END IF
         END DO
C
         IF ( ORBPOL(1:6).EQ.'BROOKS' ) THEN
            DO IS = 1,2
               CALL DCOPY(JRWSASA(IM),AOPT(1,IS,IT),1,U1,1)
C
               DO IR = 1,JRCRI(IM)
                  AOPT(IR,IS,IT) = YLAG(R(IR,IM),RASA(1,IM),U1,0,3,
     &                             JRWSASA(IM))
               END DO
            END DO
         END IF
      END DO
C
C-----------------------------------------------------------------------
C                   shift muffin-tin zero
C-----------------------------------------------------------------------
C
      VOL = 0D0
      VAV = 0D0
C
      DO IT = ITBOT,ITTOP
C
         IM = IMT(IT)
C
         NPANEL = NPAN(IM)
         IRMTIN = JRMT(IM)
         IRCRIT = JRCUT(NPANEL,IM)
C
         RINTVOL = 0D0
         RINTVAV = 0D0
         DO IR = IRMTIN + 1,IRCRIT
            IRSF = IR - IRMTIN
C
            U1_IR = FLMSF(IRSF,1,IM)*R2DRDI_W_RADINT(IR,IM)
            RINTVOL = RINTVOL + U1_IR
            RINTVAV = RINTVAV + U1_IR*VT(IR,IT)
         END DO
C
         VOL = VOL + RINTVOL*CONC(IT)*NAT(IT)
         VAV = VAV + RINTVAV*CONC(IT)*NAT(IT)
      END DO
C
      VMTZ_ASA = VAV/VOL
C
      WRITE (6,99011) VMTZ_ASA
C
      DO IT = ITBOT_EXTENDED,ITTOP_EXTENDED
         IM = IMT(IT)
         DO IR = 1,JRCUT(NPAN(IM),IM)
            VT(IR,IT) = VT(IR,IT) - VMTZ_ASA
         END DO
      END DO
C
      IF ( SYSTEM_TYPE(1:16).NE.'EMBEDDED-CLUSTER' ) EFERMI = EFERMI - 
     &     VMTZ_ASA
C
C-----------------------------------------------------------------------
C        convolute spherical  ASA  potential with shape function
C-----------------------------------------------------------------------
C
      WRITE (6,99010)
      DO IT = ITBOT_EXTENDED,ITTOP_EXTENDED
C
         IM = IMT(IT)
         WRITE (6,99007) ' '
         WRITE (6,99008) 'atom type',IT,TXT_T(IT)(1:LTXT_T(IT))
         WRITE (6,99008) 'NFP ',NFPT(IT),'mesh IM ',IM
         WRITE (6,99009) (LMIFP(IFP,IT),L_LM(LMIFP(IFP,IT)),
     &                   M_LM(LMIFP(IFP,IT)),IFP=1,NFPT(IT))
C
         IRMTIN = JRCUT(1,IM)
C
         DO LM = 2,NLMFPT(IT)
            IF ( KLMSF(LM,IM).EQ.1 ) THEN
               ISF = ISFLM(LM,IM)
C
               IF ( ISF.GT.NSFMAX .OR. ISF.LT.1 )
     &              CALL STOP_MESSAGE(ROUTINE,'ISF out of range')
C
               DO IR = IRMTIN + 1,JRCRI(IM)
                  IRSF = IR - IRMTIN
                  VNST(IR,LM,IT) = VT(IR,IT)*FLMSF(IRSF,ISF,IM)
                  BNST(IR,LM,IT) = BT(IR,IT)*FLMSF(IRSF,ISF,IM)
               END DO
C
            END IF
         END DO
C
         VNST(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX) = 0D0
         BNST(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX) = 0D0
C
         DO IR = IRMTIN + 1,JRCRI(IM)
            IRSF = IR - IRMTIN
            VT(IR,IT) = VT(IR,IT)*FLMSF(IRSF,1,IM)/SQRT_4PI
            BT(IR,IT) = BT(IR,IT)*FLMSF(IRSF,1,IM)/SQRT_4PI
         END DO
C
      END DO
C
C ======================================================================
C
C-----------------------------------------------------------------------
C  for  SYSTEM_TYPE = VIV i.e. slab  set VT = 0 for left and right host
C-----------------------------------------------------------------------
C
      IF ( SYSTEM_TYPE(1:3).EQ.'VIV' ) THEN
C
         DO IQ = 1,NQ_L + NQ_I + NQ_R
            IF ( IQ.GT.NQ_L .AND. IQ.LE.(NQ_L+NQ_I) ) CYCLE
            DO IO = 1,NOQ(IQ)
               IT = ITOQ(IO,IQ)
               NFPT(IT) = 1
               NLMFPT(IT) = NLMFPMAX
               KLMFP(:,IT) = 0
               KLMFP(1,IT) = 1
               VT(:,IT) = 0.0D0
               BT(:,IT) = 0.0D0
               VNST(:,:,IT) = 0.0D0
               BNST(:,:,IT) = 0.0D0
            END DO
         END DO
C
      END IF
C
      DEALLOCATE (DRSF,SCLM,XRSF,NRPAN)
C
C ======================================================================
99001 FORMAT (/,10X,'<',A,'>: volume of unit cell not properly',/,10X,
     &        'reproduced by integration over shape functions',/)
99002 FORMAT (/,10X,'atomic interstitial and cell volumes ',
     &        '- only for comparison'//,10X,'IM',8X,
     &        'IS (ASA)        AC (ASA)        IS (FP)         AC (FP)')
99003 FORMAT (I12,3X,4F16.8)
99004 FORMAT (10X,A,4F16.8)
99005 FORMAT (10X,A,2(16X,F16.8),/)
99006 FORMAT (/,1X,79('*'),/,34X,'<',A,'>',/,14X,
     &        'starting  FULL POTENTIAL CALCULATION  using ASA data',/,
     &        1X,79('*'),/)
99007 FORMAT (10X,A,I4,:,4X,'JRCUT   ',I5,5X,'R =',F12.8,5X,F12.8)
99008 FORMAT (10X,A,I4,4X,A,I5)
99009 FORMAT (10X,'LM  ',5(I4,' (',I2,',',I3,')'),:,/,
     &        (14X,5(I4,' (',I2,',',I3,')')))
99010 FORMAT (/,10X,'full potential parameters:')
99011 FORMAT (/,1X,79('*'),/,10X,'shift of muffin-tin zero  VMTZ ',
     &        F10.6,'  for  ASA -> FP',/,1X,79('*'),/)
C
      END
C
