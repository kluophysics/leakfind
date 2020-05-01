C*==gilsum.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE GILSUM(ALF0Q,ALF1QQ_NV,ALF1QQ_VC,ALF0QO,ALF1QOQ_NV,
     &                  ALF1QOQ_VC,I_TEMP_LAT,N_TEMP_LAT,TEMP_LAT)
C   ********************************************************************
C   *                                                                  *
C   *   sum up results                                                 *
C   *   calculate final Gilbert damping parameter ALPHA                *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CPA,ONLY:USENLCPA
      USE MOD_FILES,ONLY:IPRINT,DATSET,LDATSET,DATFIL,LDATFIL,IFILDAT,
     &    IOTMP,LSYSTEM,SYSTEM,IFILBUILDBOT,WRBUILDBOT
      USE MOD_CALCMODE,ONLY:KMROT
      USE MOD_ANGMOM,ONLY:ISMT,IOMT
      USE MOD_SITES,ONLY:NQMAX,NOMAX,MROTQ,NOQ,ITOQ,IQBOT_CHI,IQTOP_CHI
      USE MOD_TYPES,ONLY:CONC,OBS_T,TXT_T,Z
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C*--GILSUM20
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='GILSUM')
C
C Dummy arguments
C
      INTEGER I_TEMP_LAT,N_TEMP_LAT
      REAL*8 TEMP_LAT
      COMPLEX*16 ALF0Q(3,3,NQMAX),ALF0QO(3,3,NQMAX,NOMAX),
     &           ALF1QOQ_NV(3,3,NQMAX,NOMAX,NQMAX),
     &           ALF1QOQ_VC(3,3,NQMAX,NOMAX,NQMAX),
     &           ALF1QQ_NV(3,3,NQMAX,NQMAX),ALF1QQ_VC(3,3,NQMAX,NQMAX)
C
C Local variables
C
      COMPLEX*16 ALF0(3,3),ALF1_NV(3,3),ALF1_VC(3,3),ALFQO_NV(:,:,:,:),
     &           ALFQO_VC(:,:,:,:),ALF_NV(3,3),ALF_VC(3,3),CHK0(3,3),
     &           CHK_NV(3,3),CHK_VC(3,3)
      REAL*8 ALPHA0_AVG(3,3),ALPHAL_NV(3,3),ALPHAL_VC(3,3),ALPHA_NV(3,3)
     &       ,ALPHA_NV_AVG(3,3),ALPHA_VC(3,3),ALPHA_VC_AVG(3,3),
     &       ANVTMP(:,:,:),AVCTMP(:,:,:),GEFF,GOVM,GOVM_QO(:,:),
     &       MTMP(3,3),MUEORB,MUESPN,X,XMAX,XMIN,XTMP(:),YMAX,YMIN,
     &       YTMP(:)
      CHARACTER*80 FILNAM,SYS,YTXT1,YTXT2
      INTEGER IC,IO,IQ,IT,JQ,LFILNAM,LSYS,LYTXT1,LYTXT2,MUE,NCURVES,
     &        NQ_CHI,NUE
      CHARACTER*20 LEG(2)
      SAVE ANVTMP,AVCTMP,XTMP,YTMP
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE XTMP,YTMP,AVCTMP,ANVTMP,ALFQO_NV,ALFQO_VC
      ALLOCATABLE GOVM_QO
C
      ALLOCATE (ALFQO_NV(3,3,NQMAX,NOMAX),ALFQO_VC(3,3,NQMAX,NOMAX))
      ALLOCATE (GOVM_QO(NQMAX,NOMAX))
C
      CALL TRACK_INFO(ROUTINE)
C
C ------------------------------------------  multiply prefactor  (1/PI)
C
      ALF0Q(:,:,:) = -ALF0Q(:,:,:)/PI
      ALF1QQ_NV(:,:,:,:) = -ALF1QQ_NV(:,:,:,:)/PI
      ALF1QQ_VC(:,:,:,:) = -ALF1QQ_VC(:,:,:,:)/PI
C
      ALF0(:,:) = 0D0
      ALF1_NV(:,:) = 0D0
      ALF1_VC(:,:) = 0D0
C
      DO IQ = IQBOT_CHI,IQTOP_CHI
         ALF0(:,:) = ALF0(:,:) + ALF0Q(:,:,IQ)
         DO JQ = IQBOT_CHI,IQTOP_CHI
            ALF1_NV(:,:) = ALF1_NV(:,:) + ALF1QQ_NV(:,:,IQ,JQ)
            ALF1_VC(:,:) = ALF1_VC(:,:) + ALF1QQ_VC(:,:,IQ,JQ)
         END DO
      END DO
C
C ------------------------------------------  multiply prefactor  (1/PI)
C
      ALF0QO(:,:,:,:) = -ALF0QO(:,:,:,:)/PI
      ALF1QOQ_NV(:,:,:,:,:) = -ALF1QOQ_NV(:,:,:,:,:)/PI
      ALF1QOQ_VC(:,:,:,:,:) = -ALF1QOQ_VC(:,:,:,:,:)/PI
C
      ALFQO_NV(:,:,:,:) = ALF0QO(:,:,:,:)
      ALFQO_VC(:,:,:,:) = ALF0QO(:,:,:,:)
C
      CHK0(:,:) = 0D0
      CHK_NV(:,:) = ALF0(:,:)
      CHK_VC(:,:) = ALF0(:,:)
C
      DO IQ = IQBOT_CHI,IQTOP_CHI
         DO IO = 1,NOQ(IQ)
            X = CONC(ITOQ(IO,IQ))
            CHK0(:,:) = CHK0(:,:) + ALF0QO(:,:,IQ,IO)*X
C
            DO JQ = IQBOT_CHI,IQTOP_CHI
               CHK_NV(:,:) = CHK_NV(:,:) + ALF1QOQ_NV(:,:,IQ,IO,JQ)*X
               CHK_VC(:,:) = CHK_VC(:,:) + ALF1QOQ_VC(:,:,IQ,IO,JQ)*X
C
               ALFQO_NV(:,:,IQ,IO) = ALFQO_NV(:,:,IQ,IO)
     &                               + ALF1QOQ_NV(:,:,IQ,IO,JQ)
               ALFQO_VC(:,:,IQ,IO) = ALFQO_VC(:,:,IQ,IO)
     &                               + ALF1QOQ_VC(:,:,IQ,IO,JQ)
            END DO
         END DO
      END DO
C
C***********************************************************************
C
C     calculate total Gilbert damping tensor
C
      DO MUE = 1,3
         DO NUE = 1,3
C
            ALF_NV(MUE,NUE) = ALF1_NV(MUE,NUE) + ALF0(MUE,NUE)
            ALF_VC(MUE,NUE) = ALF1_VC(MUE,NUE) + ALF0(MUE,NUE)
C
C           check consistency of data sets
C
            IF ( ABS(ALF0(MUE,NUE)-CHK0(MUE,NUE)).GT.1D-10 ) THEN
               WRITE (6,99014) 'ALF0 : ',MUE,NUE,ALF0(MUE,NUE)
               WRITE (6,99014) 'CHK0 : ',MUE,NUE,CHK0(MUE,NUE)
            END IF
C
            IF ( ABS(ALF_NV(MUE,NUE)-CHK_NV(MUE,NUE)).GT.1D-10 ) THEN
               WRITE (6,99014) 'ALF_NV : ',MUE,NUE,ALF_NV(MUE,NUE)
               WRITE (6,99014) 'CHK_NV : ',MUE,NUE,CHK_NV(MUE,NUE)
            END IF
C
            IF ( ABS(ALF_VC(MUE,NUE)-CHK_VC(MUE,NUE)).GT.1D-10 ) THEN
               WRITE (6,99014) 'ALF_VC : ',MUE,NUE,ALF_VC(MUE,NUE)
               WRITE (6,99014) 'CHK_VC : ',MUE,NUE,CHK_VC(MUE,NUE)
            END IF
C
C           check if tensor elements are all real
C
            IF ( ABS(DIMAG(ALF_NV(MUE,NUE))).GT.1.D-6 ) WRITE (6,99005)
     &            'NVC',MUE,NUE,ALF_NV(MUE,NUE)
            IF ( ABS(DIMAG(ALF_VC(MUE,NUE))).GT.1.D-6 ) WRITE (6,99005)
     &            ' VC',MUE,NUE,ALF_NV(MUE,NUE)
C
            ALPHA_NV(MUE,NUE) = DREAL(ALF_NV(MUE,NUE))
            ALPHA_VC(MUE,NUE) = DREAL(ALF_VC(MUE,NUE))
C
         END DO
      END DO
C
C=======================================================================
C                    write out total alpha-tensor
C=======================================================================
C
      IF ( IPRINT.GT.5 ) THEN
         WRITE (6,99002)
         WRITE (6,99003)
         CALL RMATSTRUCT('ALPHA total-matrix',ALPHA_NV,3,3,0,0,0,1D-8,6)
C
         WRITE (6,99004)
         CALL RMATSTRUCT('ALPHA total-matrix',ALPHA_VC,3,3,0,0,0,1D-8,6)
      END IF
C
      WRITE (6,99001) ALPHA_NV(1,1),DREAL(ALF0(1,1)),DREAL(ALF1_NV(1,1))
     &                ,ALPHA_VC(1,1),DREAL(ALF0(1,1)),
     &                DREAL(ALF1_VC(1,1))
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
      IF ( WRBUILDBOT ) WRITE (IFILBUILDBOT,99015)
     &                         ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                         ALPHA_NV(1,1),DREAL(ALF0(1,1)),
     &                         DREAL(ALF1_NV(1,1)),ALPHA_VC(1,1),
     &                         DREAL(ALF0(1,1)),DREAL(ALF1_VC(1,1))
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
C-------------------------------------- site and component resolved data
C---------------------------- G/MEFF not printed for empty spheres (Z=0)
C
      WRITE (6,99008)
C
      NQ_CHI = IQTOP_CHI - IQBOT_CHI + 1
      MUESPN = 0.0D0
      MUEORB = 0.0D0
C
      DO IQ = IQBOT_CHI,IQTOP_CHI
         DO IO = 1,NOQ(IQ)
            IT = ITOQ(IO,IQ)
            GEFF = 2D0*(1D0+OBS_T(0,IOMT,IT)/OBS_T(0,ISMT,IT))
            GOVM_QO(IQ,IO) = GEFF/(OBS_T(0,ISMT,IT)+OBS_T(0,IOMT,IT))
C
            IF ( Z(IT).GT.0 ) THEN
               WRITE (6,99009) IQ,IO,IT,TXT_T(IT)(1:4),CONC(IT),
     &                         OBS_T(0,ISMT,IT),OBS_T(0,IOMT,IT),GEFF,
     &                         GOVM_QO(IQ,IO)
            ELSE
               WRITE (6,99009) IQ,IO,IT,TXT_T(IT)(1:4),CONC(IT),
     &                         OBS_T(0,ISMT,IT),OBS_T(0,IOMT,IT),GEFF
            END IF
C
            MUESPN = MUESPN + OBS_T(0,ISMT,IT)*CONC(IT)
            MUEORB = MUEORB + OBS_T(0,IOMT,IT)*CONC(IT)
         END DO
      END DO
C
      GEFF = 2D0*(1D0+MUEORB/MUESPN)
      GOVM = GEFF/(MUESPN+MUEORB)
C
      WRITE (6,99010) MUESPN/DBLE(NQ_CHI),MUEORB/DBLE(NQ_CHI),GEFF,
     &                GOVM*DBLE(NQ_CHI),MUESPN,MUEORB,GEFF,GOVM
C
      WRITE (6,99011)
C
      ALPHA_VC_AVG(:,:) = 0D0
      ALPHA_NV_AVG(:,:) = 0D0
      ALPHA0_AVG(:,:) = 0D0
C
      DO IQ = IQBOT_CHI,IQTOP_CHI
         DO IO = 1,NOQ(IQ)
            IT = ITOQ(IO,IQ)
C
            X = GOVM_QO(IQ,IO)
C
            WRITE (6,99012) IQ,TXT_T(IT)(1:4),
     &                      DREAL(ALFQO_VC(1,1,IQ,IO)*X),
     &                      DREAL(ALFQO_VC(2,2,IQ,IO)*X),
     &                      DREAL(ALFQO_NV(1,1,IQ,IO)*X),
     &                      DREAL(ALFQO_NV(2,2,IQ,IO)*X),
     &                      DREAL(ALF0QO(1,1,IQ,IO)*X),
     &                      DREAL(ALF0QO(2,2,IQ,IO)*X)
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
            IF ( WRBUILDBOT .AND. .NOT.USENLCPA )
     &           WRITE (IFILBUILDBOT,99016) ROUTINE(1:LEN_TRIM(ROUTINE))
     &           ,IQ,IO,IT,DREAL(ALFQO_VC(1,1,IQ,IO)*X),
     &           DREAL(ALFQO_VC(2,2,IQ,IO)*X),
     &           DREAL(ALFQO_NV(1,1,IQ,IO)*X),
     &           DREAL(ALFQO_NV(2,2,IQ,IO)*X),DREAL(ALF0QO(1,1,IQ,IO)*X)
     &           ,DREAL(ALF0QO(2,2,IQ,IO)*X)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
            X = GOVM_QO(IQ,IO)*CONC(IT)/DBLE(NQ_CHI)
C
            ALPHA_VC_AVG(:,:) = ALPHA_VC_AVG(:,:)
     &                          + DREAL(X*ALFQO_VC(:,:,IQ,IO))
            ALPHA_NV_AVG(:,:) = ALPHA_NV_AVG(:,:)
     &                          + DREAL(X*ALFQO_NV(:,:,IQ,IO))
            ALPHA0_AVG(:,:) = ALPHA0_AVG(:,:)
     &                        + DREAL(X*ALF0QO(:,:,IQ,IO))
C
         END DO
      END DO
C
      WRITE (6,99013) ALPHA_VC_AVG(1,1),ALPHA_VC_AVG(2,2),
     &                ALPHA_NV_AVG(1,1),ALPHA_NV_AVG(2,2),
     &                ALPHA0_AVG(1,1),ALPHA0_AVG(2,2),ALPHA_VC(1,1)
     &                *GOVM,ALPHA_VC(2,2)*GOVM,ALPHA_NV(1,1)*GOVM,
     &                ALPHA_NV(2,2)*GOVM,DREAL(ALF0(1,1))*GOVM,
     &                DREAL(ALF0(1,1))*GOVM
C
C------------------- if there is a global rotation of the magnetisation:
C----------------------------------- give the tensors in the local frame
C------------------------------------ use rotation matrix MROTQ for IQ=1
C
      IF ( KMROT.EQ.2 ) THEN
C
         CALL STOP_MESSAGE(ROUTINE,'not working for rotated moment')
C
         CALL DGEMM('N','T',3,3,3,1D0,ALPHA_NV,3,MROTQ,3,0D0,MTMP,3)
         CALL DGEMM('N','N',3,3,3,1D0,MROTQ,3,MTMP,3,0D0,ALPHAL_NV,3)
C
         CALL DGEMM('N','T',3,3,3,1D0,ALPHA_VC,3,MROTQ,3,0D0,MTMP,3)
         CALL DGEMM('N','N',3,3,3,1D0,MROTQ,3,MTMP,3,0D0,ALPHAL_VC,3)
C
         WRITE (6,99006) ALPHAL_NV,ALPHAL_VC
      END IF
C
C=======================================================================
C                   write results to data file  DATFIL
C=======================================================================
C
      IF ( I_TEMP_LAT.EQ.1 ) THEN
         IF ( N_TEMP_LAT.GT.1 ) THEN
            DATFIL = DATSET(1:LDATSET)//'alfa_T.dat'
         ELSE
            DATFIL = DATSET(1:LDATSET)//'alfa_x.dat'
         END IF
         LDATFIL = LDATSET + 10
         OPEN (UNIT=IFILDAT,FILE=DATFIL(1:LDATFIL))
      END IF
C
      IF ( N_TEMP_LAT.EQ.1 ) THEN
         WRITE (IFILDAT,99007) 'CONC(1)',CONC(1),ALPHA_NV,ALPHA_VC
      ELSE
         WRITE (IFILDAT,99007) 'TMP',TEMP_LAT,ALPHA_NV,ALPHA_VC
      END IF
C
C=======================================================================
      IF ( N_TEMP_LAT.EQ.1 ) RETURN
C=======================================================================
C
C=======================================================================
C                create xmgrace file for ALFA(T)
C=======================================================================
      IF ( I_TEMP_LAT.EQ.1 ) ALLOCATE (XTMP(N_TEMP_LAT),YTMP(N_TEMP_LAT)
     &                                 ,ANVTMP(3,3,N_TEMP_LAT),
     &                                 AVCTMP(3,3,N_TEMP_LAT))
      XTMP(I_TEMP_LAT) = TEMP_LAT
      ANVTMP(1:3,1:3,I_TEMP_LAT) = ALPHA_NV(1:3,1:3)
      AVCTMP(1:3,1:3,I_TEMP_LAT) = ALPHA_VC(1:3,1:3)
C
C ----------------------------------------------------------------------
      IF ( I_TEMP_LAT.EQ.N_TEMP_LAT ) THEN
C
         XMIN = XTMP(1)
         XMAX = XTMP(N_TEMP_LAT)
         YMIN = 0D0
         YMAX = 0D0
         DO IC = 1,N_TEMP_LAT
            DO MUE = 1,2
               ANVTMP(MUE,MUE,IC) = ABS(ANVTMP(MUE,MUE,IC))
               AVCTMP(MUE,MUE,IC) = ABS(AVCTMP(MUE,MUE,IC))
               YMAX = MAX(YMAX,ANVTMP(MUE,MUE,IC))
               YMAX = MAX(YMAX,AVCTMP(MUE,MUE,IC))
            END DO
         END DO
C
         YTXT1 = '!xa!0(T) (d.u.)'
         LYTXT1 = 15
         YTXT2 = ' '
         LYTXT2 = 1
         SYS = SYSTEM
         LSYS = LSYSTEM
         CALL XMGRSUBSCRIPTS(SYS,LSYS,80)
C
         CALL XMGRHEAD(DATSET,LDATSET,'alfa_T',6,' ',0,FILNAM,80,
     &                 LFILNAM,IOTMP,1,XMIN,1,XMAX,1,YMIN,0,YMAX,1,YMIN,
     &                 0,YMAX,1,'temperature (T)',15,YTXT1,LYTXT1,YTXT2,
     &                 LYTXT2,'SPR-KKR calculations for '//SYS(1:LSYS),
     &                 25+LSYS,'Gilbert damping parameter !xa!0(T)',34,
     &                 .FALSE.)
C
         NCURVES = 2
C
         CALL XMGRCURVES(IOTMP,1,NCURVES,NCURVES,2,1,0)
C
         LEG(1) = 'NV'
         LEG(2) = 'VC'
C
         CALL XMGRLEGEND(IOTMP,1,NCURVES,NCURVES,LEG,LEG)
C
         YTMP(1:N_TEMP_LAT) = ANVTMP(1,1,1:N_TEMP_LAT)
         CALL XMGRTABLE(0,0,XTMP,YTMP,1D0,N_TEMP_LAT,IOTMP)
C
         YTMP(1:N_TEMP_LAT) = AVCTMP(1,1,1:N_TEMP_LAT)
         CALL XMGRTABLE(0,1,XTMP,YTMP,1D0,N_TEMP_LAT,IOTMP)
C
         WRITE (6,*) '  '
         WRITE (6,*) '   ALFA written to file  ',FILNAM(1:LFILNAM)
         CLOSE (IOTMP)
C
      END IF
C ----------------------------------------------------------------------
C
C=======================================================================
99001 FORMAT (/,10X,'Gilbert damping parameter ALPHA',//,23X,'total ',
     &        12X,'term 0',9X,'term 1 ',/,10X,'(-VC)',3X,F15.10,3X,
     &        2F15.10,/,10X,'(+VC)',3X,F15.10,3X,2F15.10,//,10X,
     &        'to be scaled by factor g_eff / m_eff',/)
99002 FORMAT (/,'ALPHA tensor total:',/,18('='),/)
99003 FORMAT ('without vertex-corrections:')
99004 FORMAT ('including vertex-corrections:')
99005 FORMAT ('ALPHAS_',a3,' not real!!! mue/nue/ALPHA=',2I4,2E13.5)
99006 FORMAT (/,5X,'Gilbert damping tensors for the local frame',/,5X,
     &        'ALPHA  (NV)        ',3E13.5,/,5X,'[.........]        ',
     &        3E13.5,/,23X,3E13.5,/,5X,'ALPHA  (VC)        ',3E13.5,/,
     &        5X,'[.........]        ',3E13.5,/,23X,3E13.5,/)
99007 FORMAT ('#  tabulating ',A,', ALPHAL_NV, ALPHAL_VC',/,F8.3,
     &        6(2X,3E13.5))
99008 FORMAT (/,5X,'magnetic moment information',
     &        '   with  GEFF = 2(SMT+OMT)/SMT',//,5X,'IQ   IO',
     &        '   IT  TXT     CONC    SMT     OMT     GEFF   G/MEFF')
99009 FORMAT (I7,2I5,2X,A,5F8.3)
99010 FORMAT (5X,'average (site)',12X,4F8.3,/,5X,'average (sum) ',12X,
     &        4F8.3)
99011 FORMAT (//,5X,'site and component resolved ALPHA',
     &        ' (including factor g_eff / m_eff)',//,23X,'alfa (+VC)',
     &        12X,'alfa (-VC)',12X,'alfa (0)',/,5X,'IQ  TXT',6X,
     &        '     xx        yy ',2(9X,'xx       yy'))
99012 FORMAT (I7,2X,A,6X,2F10.6,2X,2F9.5,2X,2F9.5)
99013 FORMAT (5X,'average (site)',2F10.6,2X,2F9.5,2X,2F9.5,/,5X,
     &        'average (sum) ',2F10.6,2X,2F9.5,2X,2F9.5,/,80('#'))
99014 FORMAT ('>>>>>> INCONSISTENCY ',A,2I3,2F14.10)
99015 FORMAT ('# BUILDBOT: ',A,
     &        ':  ALFA  TOTAL  term 0  term 1  (-VC)  (+VC)',/,
     &        (1PE22.14))
99016 FORMAT ('# BUILDBOT: ',A,':  ALFA (xx yy)^VC (xx yy)^NV (xx yy)^0'
     &        ,'  for IQ, IO, IT =',3I5,/,(1PE22.14))
      END
