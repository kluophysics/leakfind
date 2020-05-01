C*==signlcsum.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIGNLCSUM(SIG1_NVQ,SIG1_VCQ,SIG0Q,SIGMAAU_NV,
     &                     SIGMAAU_VC,RHOAU_NV,RHOAU_VC,SIGMA_NV,
     &                     SIGMA_VC,RHO_NV,RHO_VC,IPRINT)
C   ********************************************************************
C   *                                                                  *
C   *   sum up results                                                 *
C   *   calculate final resistivity and conductivity                   *
C   *   scalar  and  tensor  form                                      *
C   *                                                                  *
C   *   NOTE: JB gives HBAR etc. in SI units                           *
C   *                                                                  *
C   * 06/10/99  HE  based on JB's routine                              *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:KMROT
      USE MOD_SITES,ONLY:IMQ,NQMAX,MROTQ,IQBOT_CHI,IQTOP_CHI
      USE MOD_RMESH,ONLY:RWS
      USE MOD_CONSTANTS,ONLY:A0_SI,E0_SI,HBAR_SI,M0_SI,PI
      USE MOD_FILES,ONLY:IFILBUILDBOT,WRBUILDBOT
      IMPLICIT NONE
C*--SIGNLCSUM22
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIGNLCSUM')
C
C Dummy arguments
C
      INTEGER IPRINT
      REAL*8 RHOAU_NV,RHOAU_VC,SIGMAAU_NV,SIGMAAU_VC
      REAL*8 RHO_NV(3,3),RHO_VC(3,3),SIGMA_NV(3,3),SIGMA_VC(3,3)
      COMPLEX*16 SIG0Q(3,3,NQMAX),SIG1_NVQ(3,3,NQMAX,NQMAX),
     &           SIG1_VCQ(3,3,NQMAX,NQMAX)
C
C Local variables
C
      REAL*8 CONSI,DIFF,MTMP(3,3),PREFAC,RESIST_NV,RESIST_VC,
     &       RHOL_NV(3,3),RHOL_VC(3,3),SIG0AU,SIG1AU_NV,SIG1AU_VC,
     &       SIGMAL_NV(3,3),SIGMAL_VC(3,3),VUC
      INTEGER IQ,IQ1,IQ2,MUE,NUE
      COMPLEX*16 SIG0(3,3),SIG1_NV(3,3),SIG1_VC(3,3),SIGTOT_NV(3,3),
     &           SIGTOT_VC(3,3)
C
C*** End of declarations rewritten by SPAG
C
C ------------------------- conversion of resistivity result to SI-units
C -------------------------------  get volume of unit cell VUC in  [m^3]
C
      VUC = 0.0D0
      DO IQ = IQBOT_CHI,IQTOP_CHI
         VUC = VUC + RWS(IMQ(IQ))**3
      END DO
      VUC = (4.D0*PI/3.D0)*VUC*A0_SI**3
      PREFAC = (E0_SI*HBAR_SI*A0_SI/M0_SI)**2
      CONSI = PREFAC*4.D0*M0_SI**2/(PI*HBAR_SI**3*VUC)
C
      DO MUE = 1,3
         DO NUE = 1,3
            SIG0(MUE,NUE) = 0D0
            SIG1_NV(MUE,NUE) = 0D0
            SIG1_VC(MUE,NUE) = 0D0
         END DO
      END DO
C
      DO IQ1 = IQBOT_CHI,IQTOP_CHI
         DO MUE = 1,3
            DO NUE = 1,3
               SIG0(MUE,NUE) = SIG0(MUE,NUE) + SIG0Q(MUE,NUE,IQ1)
               DO IQ2 = IQBOT_CHI,IQTOP_CHI
                  SIG1_NV(MUE,NUE) = SIG1_NV(MUE,NUE)
     &                               + SIG1_NVQ(MUE,NUE,IQ1,IQ2)
                  SIG1_VC(MUE,NUE) = SIG1_VC(MUE,NUE)
     &                               + SIG1_VCQ(MUE,NUE,IQ1,IQ2)
               END DO
            END DO
         END DO
      END DO
C
C***********************************************************************
C***********************************************************************
C***********************************************************************
C --------------------------------------- SUPPRESS OFF-DIAGONAL ELEMENTS
C
      DO MUE = 1,3
         DO NUE = 1,3
            IF ( MUE.NE.NUE ) THEN
               SIG0(MUE,NUE) = 0D0
               SIG1_NV(MUE,NUE) = 0D0
               SIG1_VC(MUE,NUE) = 0D0
            END IF
         END DO
      END DO
C***********************************************************************
C***********************************************************************
C***********************************************************************
C
C     calculate total conductivity tensor:
C     ------------------------------------
      DO MUE = 1,3
         DO NUE = 1,3
C
            SIGTOT_NV(MUE,NUE) = SIG1_NV(MUE,NUE) + SIG0(MUE,NUE)
            SIGTOT_VC(MUE,NUE) = SIG1_VC(MUE,NUE) + SIG0(MUE,NUE)
C
C           check if tensor elements are all real:
C           --------------------------------------
            IF ( ABS(DIMAG(SIGTOT_NV(MUE,NUE))).GT.1.D-6 )
     &           WRITE (6,99007) 'NVC',MUE,NUE,SIGTOT_NV(MUE,NUE)
            IF ( ABS(DIMAG(SIGTOT_VC(MUE,NUE))).GT.1.D-6 )
     &           WRITE (6,99007) ' VC',MUE,NUE,SIGTOT_NV(MUE,NUE)
C
            SIGMA_NV(MUE,NUE) = DREAL(SIGTOT_NV(MUE,NUE))
            SIGMA_VC(MUE,NUE) = DREAL(SIGTOT_VC(MUE,NUE))
C
         END DO
      END DO
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
      IF ( WRBUILDBOT ) WRITE (IFILBUILDBOT,99011)
     &                         ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                         ((SIGMA_NV(MUE,NUE),NUE=1,3),MUE=1,3),
     &                         ((SIGMA_VC(MUE,NUE),NUE=1,3),MUE=1,3)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
C     write out total sigma and rho-tensors
C     -------------------------------------
      WRITE (6,99001)
C
      WRITE (6,99002)
      CALL RMATSTRUCT('Sigma total-matrix',SIGMA_NV,3,3,0,0,0,1D-8,6)
C
      CALL RMATINV(3,3,SIGMA_NV,RHO_NV)
C
      CALL RMATSTRUCT('Resistivity-matrix',RHO_NV,3,3,0,0,0,1D-8,6)
C
      WRITE (6,99003)
      CALL RMATSTRUCT('Sigma total-matrix',SIGMA_VC,3,3,0,0,0,1D-8,6)
C
      CALL RMATINV(3,3,SIGMA_VC,RHO_VC)
C
      CALL RMATSTRUCT('Resistivity-matrix',RHO_VC,3,3,0,0,0,1D-8,6)
C
C------------------- if there is a global rotation of the magnetisation:
C----------------------------------- give the tensors in the local frame
C------------------------------------ use rotation matrix MROTQ for IQ=1
C
      IF ( KMROT.EQ.2 ) THEN
         CALL DGEMM('N','T',3,3,3,1D0,SIGMA_NV,3,MROTQ,3,0D0,MTMP,3)
         CALL DGEMM('N','N',3,3,3,1D0,MROTQ,3,MTMP,3,0D0,SIGMAL_NV,3)
C
         CALL DGEMM('N','T',3,3,3,1D0,SIGMA_VC,3,MROTQ,3,0D0,MTMP,3)
         CALL DGEMM('N','N',3,3,3,1D0,MROTQ,3,MTMP,3,0D0,SIGMAL_VC,3)
C
         CALL DGEMM('N','T',3,3,3,1D0,RHO_NV,3,MROTQ,3,0D0,MTMP,3)
         CALL DGEMM('N','N',3,3,3,1D0,MROTQ,3,MTMP,3,0D0,RHOL_NV,3)
C
         CALL DGEMM('N','T',3,3,3,1D0,RHO_VC,3,MROTQ,3,0D0,MTMP,3)
         CALL DGEMM('N','N',3,3,3,1D0,MROTQ,3,MTMP,3,0D0,RHOL_VC,3)
C
         WRITE (6,99010) SIGMAL_NV,SIGMAL_VC,RHOL_NV,RHOL_VC
      END IF
C
C     remark sigtot_vc/nvc is overwritten in sigmaout
C
C     sum up diagonal elements:
C     -------------------------
      SIG0AU = 0.D0
      SIG1AU_NV = 0.D0
      SIG1AU_VC = 0.D0
      RHOAU_NV = 0.D0
      RHOAU_VC = 0.D0
      DO MUE = 1,3
         SIG0AU = SIG0AU + DREAL(SIG0(MUE,MUE))/3D0
         SIG1AU_NV = SIG1AU_NV + DREAL(SIG1_NV(MUE,MUE))/3D0
         SIG1AU_VC = SIG1AU_VC + DREAL(SIG1_VC(MUE,MUE))/3D0
         RHOAU_NV = RHOAU_NV + RHO_NV(MUE,MUE)/3D0
         RHOAU_VC = RHOAU_VC + RHO_VC(MUE,MUE)/3D0
      END DO
      SIGMAAU_NV = SIG1AU_NV + SIG0AU
      SIGMAAU_VC = SIG1AU_VC + SIG0AU
C
C        calculate resistivities in SI-units
C        -----------------------------------
      RESIST_NV = RHOAU_NV*1.D8/CONSI
      RESIST_VC = RHOAU_VC*1.D8/CONSI
      IF ( RESIST_NV.LT.0 ) WRITE (6,99009)
      IF ( RESIST_VC.LT.0 ) WRITE (6,99009)
C
      DIFF = 100.D0*(SIGMAAU_VC-SIGMAAU_NV)/SIGMAAU_NV
C
      IF ( IPRINT.NE.0 ) THEN
         WRITE (6,99004)
         WRITE (6,99008) SIGMAAU_NV,SIGMAAU_NV*CONSI,
     &                   1.D8/(SIGMAAU_NV*CONSI),RESIST_NV,' (no VC)  '
         WRITE (6,99005)
         WRITE (6,99008) SIGMAAU_VC,SIGMAAU_VC*CONSI,
     &                   1.D8/(SIGMAAU_VC*CONSI),RESIST_VC,' (with VC)'
         WRITE (6,99006) DIFF
         IF ( DIFF.LT.0D0 ) WRITE (6,*) 'THIS IS STRANGE,PLEASE CHECK!!'
      END IF
99001 FORMAT (/,'Sigmatensor total:',/,18('='),/)
99002 FORMAT ('without vertex-corrections:')
99003 FORMAT ('including vertex-corrections:')
99004 FORMAT (/,'Results without Vertex-corrections:',/,35('='))
99005 FORMAT (/,'Results including Vertex-corrections:',/,37('='))
99006 FORMAT ('Inclusion of vertex-correction increases',
     &        ' conductivity by ',f8.4,'%',/)
99007 FORMAT ('Sigma_',a3,' not real!!! mue/nue/sigma=',2I4,2E13.5)
99008 FORMAT (/,'Conductivity in a.u.:          ',f18.5,/,
     &        'Conductivity in [1/(Ohm*m)]    ',f18.5,/,
     &        'Inverse Conductivity [muOhm.cm]',f18.5,/,
     &        'Resistivity in [muOhm.cm]      ',f18.5,a10)
99009 FORMAT (/,45('!'),/,' Negative resistivity...bad luck... ',
     &        'try again',/,45('!'))
99010 FORMAT (/,5X,'galvano-magnetic tensors for the local frame',/,5X,
     &        'SIGMA  (NV)        ',3E13.5,/,5X,'[1/(Ohm*m)]        ',
     &        3E13.5,/,23X,3E13.5,/,5X,'SIGMA  (VC)        ',3E13.5,/,
     &        5X,'[1/(Ohm*m)]        ',3E13.5,/,23X,3E13.5,/,5X,
     &        'RHO    (NV)        ',3E13.5,/,5X,'[muOhm.cm]         ',
     &        3E13.5,/,23X,3E13.5,/,5X,'RHO    (VC)        ',3E13.5,/,
     &        5X,'[muOhm.cm]         ',3E13.5,/,23X,3E13.5,/)
99011 FORMAT ('# BUILDBOT: ',A,':  SIG tensor (a.u.) NV -- VC',/,
     &        (1PE22.14))
      END
