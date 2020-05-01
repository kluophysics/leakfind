C*==photonvec.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE PHOTONVEC(AA,THQ,FIQ,ALQ,DELQ,OMHAR,ICIRC,IDREH)
C     /****************************************************************/
C     # purpose      : determine vectorpotential                       *
C                      of the electromagnetic field                    *
C                                                                      *
C     # parameter:                                                     *
C       input:                                                         *
C          angles:                                                     *
C          thq=thetaq, fiq= phiq, alq=alpha, delq=delta                *
C          (alpha=polarisation angle, delta=polarisation phase)        *
C          icirc  = 0 elliptically polarized                           *
C             = 1 linear polarized                                     *
C             = 2 circularly polarized                                 *
C             = 3 polarization changed by metal-optics (fresnel)       *
C          note: idreh=0 results always in linearly polarized light    *
C       output:                                                        *
C         aa                                                           *
C                                                                      *
C     # attention:                                                     *
C       it is assumed that the z-axis is the positive surface normal ! *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:PI,CIMAG
      USE MOD_SPEC_QXQY,ONLY:QX,QY
      USE MOD_SPEC_OUTPUT,ONLY:NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALQ,DELQ,FIQ,OMHAR,THQ
      INTEGER ICIRC,IDREH
      COMPLEX*16 AA(3)
C
C Local variables
C
      REAL*8 AL,AXZ,DE,PH,QB,QZ,RD,S0,S1,S2,S3,TH
      COMPLEX*16 COSBXZ,EPSI,NIK,SINBXZ,SP,SS,TP,TS
C
C*** End of declarations rewritten by SPAG
C
      RD = PI/180.D0
C
      TH = THQ*RD
      PH = FIQ*RD
      AL = ALQ*RD
C
C
C     set phase fixed for pure linear or circular polarization
C     just to check, because it should be already corrected in input
      IF ( ICIRC.EQ.1 ) DELQ = 0.D0
      IF ( ICIRC.EQ.2 ) DELQ = 90.D0
C     these ifs should be removed later
C
      DE = DBLE(IDREH)*DELQ*RD
      QB = OMHAR/137.036D0
C
C     vector-potential and propagation vector in vacuum
      IF ( (ABS(THQ).LT.0.1D0 .OR. ABS(THQ-180.0D0).LT.0.1D0) .AND. 
     &     ICIRC.EQ.2 ) THEN
         AA(1) = SQRT(0.5D0)
         AA(2) = DBLE(IDREH)*SQRT(0.5D0)*CIMAG
         AA(3) = 0.D0
C
         QX = 0.D0
         QY = 0.D0
         QZ = QB
         IF ( ABS(THQ-180.0D0).LT.0.1D0 ) THEN
            QZ = -QZ
            AA(2) = -AA(2)
         END IF
      ELSE
         AA(1) = COS(TH)*COS(PH)*COS(AL) - SIN(PH)*SIN(AL)
     &           *CDEXP(CIMAG*DE)
         AA(2) = COS(TH)*SIN(PH)*COS(AL) + COS(PH)*SIN(AL)
     &           *CDEXP(CIMAG*DE)
         AA(3) = -SIN(TH)*COS(AL)
C
         QX = SIN(TH)*COS(PH)*QB
         QY = SIN(TH)*SIN(PH)*QB
         QZ = COS(TH)*QB
      END IF
C
C     calculate with fresnel optics
      IF ( ICIRC.GE.3 ) THEN
C
C         the values of the complex index of refraction nik
C         may be read in from a file
C
C         example for non-unit index of refraction
C         nik = dcmplx(-2.66d0,5.13d0) ! 4.3 eV
C         nik = dcmplx(-2.51d0,5.01d0) ! 4.4 eV
C         nik = dcmplx(-2.36d0,4.90d0) ! 4.5 eV
C         nik = dcmplx(-2.22d0,4.82d0) ! 4.6 eV
C         nik = dcmplx(-2.13d0,4.74d0) ! 4.7 eV
C         nik = dcmplx(-2.04d0,4.65d0) ! 4.8 eV
C         nik = dcmplx(-1.95d0,4.56d0) ! 4.9 eV
C         nik = dcmplx(-1.87d0,4.49d0) ! 5.0 eV
C
         NIK = DCMPLX(0.8D0,0.3D0)
         EPSI = NIK*NIK
C
         AXZ = ATAN2(ABS(AA(1)),ABS(AA(3)))
         SINBXZ = SIN(AXZ)/NIK
         COSBXZ = CDSQRT(1.D0-SINBXZ*SINBXZ)
C
         TS = 2.D0*COS(AXZ)/(COS(AXZ)+CDSQRT(EPSI-SIN(AXZ)**2))
         TP = 2.D0*NIK*COS(AXZ)/(EPSI*COS(AXZ)+CDSQRT(EPSI-SIN(AXZ)**2))
C
         AA(1) = -COSBXZ*TP*AA(1)
         AA(2) = TS*AA(2)
         AA(3) = SINBXZ*TP*AA(3)
      END IF
C
      WRITE (NOUT1,99001) QX,QY,QZ
      WRITE (NOUT1,99002)
      WRITE (NOUT1,99003) AA(1),AA(2),AA(3)
C
      SP = CDSQRT(AA(1)**2+AA(3)**2)
      SS = AA(2)
      S0 = ABS(SP)**2 + ABS(SS)**2
      S1 = ABS(SP)**2 - ABS(SS)**2
      S2 = 2.D0*DBLE(SP*DCONJG(SS))
      S3 = -2.D0*DIMAG(SP*DCONJG(SS))
      WRITE (NOUT1,99004) S0,100.D0*S1/S0,100.D0*S2/S0,100.D0*S3/S0
      RETURN
C
99001 FORMAT (1x,'vacuum-wavevector of the photonfield:',3(1x,f12.8))
99002 FORMAT (1x,'vector potential  of the photonfield:')
99003 FORMAT (1x,'aa(i=1,3)= ',3(1x,'(',2F12.8,')'))
99004 FORMAT (1x,'stookes-vector    of the photonfield',1x,e12.5,
     &        3(1x,f8.3,'%'))
      END
C*==detmat.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE DETMAT(AA,BB,VLM,AMATX_A,AMATY_A,AMATV_A,AMAT_A,
     &                  AMATX_NA,AMATY_NA,AMATV_NA,AMAT_NA,AMAT_AT,
     &                  AMAT_NAT,NATL,LAYS)
C     /****************************************************************/
C     # purpose      : prepares the calculation of the                 *
C                      angular matrix elements                         *
C                                                                      *
C     # calls the subroutines:                                         *
C       kmatpfp      cind                                              *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,ML,MQD,XMAXE,NFULLPOT,EPS12,CZERO,
     &    NTPHOMAX
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      USE MOD_SPEC_RINDC,ONLY:NRL,NRM
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER LAYS
      COMPLEX*16 AA(3),AMAT_A(XMAXE,4,NFULLPOT,2),
     &           AMAT_AT(XMAXE,4,NFULLPOT,2),AMAT_NA(XMAXE,4,NFULLPOT),
     &           AMAT_NAT(XMAXE,4,NFULLPOT),
     &           VLM(LAYSM,NATLM,MQD,NRMAX,NTPHOMAX)
      INTEGER AMATV_A(XMAXE,4,NFULLPOT,2),AMATV_NA(XMAXE,4,NFULLPOT),
     &        AMATX_A(MQD+1,4,NFULLPOT,2),AMATX_NA(MQD+1,4,NFULLPOT),
     &        AMATY_A(MQD+1,4,NFULLPOT,2),AMATY_NA(MQD+1,4,NFULLPOT),
     &        NATL(LAYSM)
      REAL*8 BB(3)
C
C Local variables
C
      INTEGER A,AMX(:,:,:),AMY(:,:),ANGULAR_TYPE,I,IK,IKS,IM,IMS,ISI,J,
     &        K,KS,L1,LK,LKS,LMP,M1,NSI
      COMPLEX*16 AMV(:,:),AMVT(:,:),AXB(3),DPP(:,:),YLM(3),YLM1(3)
      CHARACTER*3 ANAME(8)
      REAL*8 BFIELD,MU,MUS
      INTEGER KAPMUE2K
      EXTERNAL KAPMUE2K
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE AMV,AMX,AMY,DPP,AMVT
      ALLOCATE (AMV(XMAXE,8),AMX(MQD,2,8),AMY(XMAXE,8),DPP(MQD,MQD))
      ALLOCATE (AMVT(XMAXE,8))
C
C*** End of declarations rewritten by SPAG
C
C     full potential angular matrix elements:
C     =======================================
C     a dependend matrix elements are d and f:
C     ========================================
C     not a-dependend matrixelements are a and e:
C     ===========================================
C
      DATA ANAME/'d++','d--','f++','f--','a+-','a-+','e++','e--'/
C
      BFIELD = SQRT(BB(1)**2+BB(2)**2+BB(3)**2)
C
      CALL KMATPFP(AA,BB,BFIELD,VLM,YLM,YLM1,AXB,AMATX_A,AMATY_A,
     &             AMATV_A,AMAT_A,AMATX_NA,AMATY_NA,AMATV_NA,AMAT_NA,
     &             NATL,LAYS)
C
      CALL KMATPFPT(AA,BB,BFIELD,VLM,YLM,YLM1,AXB,AMATX_A,AMATY_A,
     &              AMATV_A,AMAT_AT,AMATX_NA,AMATY_NA,AMATV_NA,AMAT_NAT,
     &              NATL,LAYS)
C
C     some output follows
      IF ( IP.GT.1 ) THEN
         WRITE (NOUT1,'(a21)') 'vectorpotential used:'
         WRITE (NOUT1,'(''aa(1)'',2f12.8)') AA(1)
         WRITE (NOUT1,'(''aa(2)'',2f12.8)') AA(2)
         WRITE (NOUT1,'(''aa(3)'',2f12.8)') AA(3)
C
         WRITE (NOUT1,'(a33)') 'ylms used for d++, d--, f++, f--:'
         WRITE (NOUT1,'(''ylm1(1)'',2f12.8)') YLM1(1)
         WRITE (NOUT1,'(''ylm1(2)'',2f12.8)') YLM1(2)
         WRITE (NOUT1,'(''ylm1(3)'',2f12.8)') YLM1(3)
C
         WRITE (NOUT1,'(a23)') 'ylms used for a+-, a-+:'
         WRITE (NOUT1,'(''ylm(1)'',2f12.8)') YLM(1)
         WRITE (NOUT1,'(''ylm(2)'',2f12.8)') YLM(2)
         WRITE (NOUT1,'(''ylm(3)'',2f12.8)') YLM(3)
C
         IF ( BFIELD.GT.EPS12 ) THEN
            WRITE (NOUT1,'(a15)') 'magnetic field:'
            WRITE (NOUT1,'(''b(1)'',2f12.8)') BB(1)
            WRITE (NOUT1,'(''b(2)'',2f12.8)') BB(2)
            WRITE (NOUT1,'(''b(3)'',2f12.8)') BB(3)
            WRITE (NOUT1,'(''b   '',2f12.8)') BFIELD
C
            WRITE (NOUT1,'(a6)') 'a x b:'
            WRITE (NOUT1,'(''axb(1)'',2f12.8)') AXB(1)
            WRITE (NOUT1,'(''axb(2)'',2f12.8)') AXB(2)
            WRITE (NOUT1,'(''axb(3)'',2f12.8)') AXB(3)
         END IF
      END IF
C
      IF ( IP.GE.3 ) THEN
         WRITE (NOUT1,99001)
         WRITE (NOUT1,99002)
         L1 = 0
         M1 = 0
         DO LMP = 1,NFULLPOT
            IF ( NFULLPOT.NE.1 ) THEN
               L1 = NRL(LMP)
               M1 = NRM(LMP)
            END IF
            DO A = 1,2
               CALL CIND(LMP,A,AMATX_A,AMATY_A,AMATV_A,AMAT_A,AMATX_NA,
     &                   AMATY_NA,AMATV_NA,AMX,AMY,AMV,AMAT_AT,AMAT_NAT,
     &                   AMVT)
               DO ANGULAR_TYPE = 1,4
                  WRITE (NOUT1,99003) LMP,A,ANGULAR_TYPE
                  DO I = 1,MQD
                     DO J = 1,MQD
                        DPP(I,J) = CZERO
                     END DO
                  END DO
                  DO K = 1,MQD
                     DO IK = AMX(K,1,ANGULAR_TYPE),AMX(K,2,ANGULAR_TYPE)
                        KS = AMY(IK,ANGULAR_TYPE)
                        DPP(K,KS) = AMV(IK,ANGULAR_TYPE)
                     END DO
                  END DO
C
                  DO LK = 1,ML
                     DO NSI = 1,2
                        ISI = (-1)**(NSI+1)
                        K = ISI*(LK-1) + (ISI-1)/2
                        IF ( K.NE.0 ) THEN
                           DO LKS = 1,3
                              IF ( LKS.EQ.1 ) KS = K - 1
                              IF ( LKS.EQ.2 ) KS = -K
                              IF ( LKS.EQ.3 ) KS = K + 1
                              IF ( KS.NE.0 .AND. KS.GE.-ML .AND. 
     &                             KS.LE.ML-1 ) THEN
                                 DO IM = -(2*ABS(K)-1),2*ABS(K) - 1,2
                                    MU = IM/2.D0
                                    IK = KAPMUE2K(K,MU,ML)
                                    DO IMS = -(2*ABS(KS)-1),2*ABS(KS)
     &                                 - 1,2
                                       MUS = IMS/2.D0
                                       IKS = KAPMUE2K(KS,MUS,ML)
                                       IF ( CDABS(DPP(IK,IKS))
     &                                    .GT.1D-16 )
     &                                    WRITE (NOUT1,99005)
     &                                    ANAME(ANGULAR_TYPE),K,IM,KS,
     &                                    IMS,IK,IKS,L1,M1,(IMS+IM)/2,A,
     &                                    DPP(IK,IKS)
                                    END DO
                              !mus
                                 END DO
                              END IF
                              !mu
                           END DO
                        END IF
                              !lks
                     END DO   !nsi
                  END DO      !lk
               END DO         !angular_type
            END DO            !a
C
C         matrix elements not depending on a
            DO ANGULAR_TYPE = 5,8
               WRITE (NOUT1,99004) LMP,ANGULAR_TYPE
               DO I = 1,MQD
                  DO J = 1,MQD
                     DPP(I,J) = CZERO
                  END DO
               END DO
               DO K = 1,MQD
                  DO IK = AMX(K,1,ANGULAR_TYPE),AMX(K,2,ANGULAR_TYPE)
                     KS = AMY(IK,ANGULAR_TYPE)
                     DPP(K,KS) = AMV(IK,ANGULAR_TYPE)
                  END DO
               END DO
C
               DO LK = 1,ML
                  DO NSI = 1,2
                     ISI = (-1)**(NSI+1)
                     K = ISI*(LK-1) + (ISI-1)/2
                     IF ( K.NE.0 ) THEN
                        DO LKS = 1,3
                           IF ( LKS.EQ.1 ) KS = K - 1
                           IF ( LKS.EQ.2 ) KS = -K
                           IF ( LKS.EQ.3 ) KS = K + 1
                           IF ( KS.NE.0 .AND. KS.GE.-ML .AND. 
     &                          KS.LE.ML-1 ) THEN
                              DO IM = -(2*ABS(K)-1),2*ABS(K) - 1,2
                                 MU = IM/2.D0
                                 IK = KAPMUE2K(K,MU,ML)
                                 DO IMS = -(2*ABS(KS)-1),2*ABS(KS) - 1,2
                                    MUS = IMS/2.D0
                                    IKS = KAPMUE2K(KS,MUS,ML)
                                    IF ( CDABS(DPP(IK,IKS)).GT.1D-16 )
     &                                 WRITE (NOUT1,99006)
     &                                 ANAME(ANGULAR_TYPE),K,IM,KS,IMS,
     &                                 IK,IKS,L1,M1,(IMS+IM)/2,
     &                                 DPP(IK,IKS)
                                 END DO
                              END DO
                           END IF
                        END DO
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END IF
C
      RETURN
C
99001 FORMAT ('# non zero angular matrix elements:')
99002 FORMAT (1x,'m+-((k,mu) -> (ks,mus),...)')
99003 FORMAT ('# lmp=',i2,'  a=',i2,'  atype=',i2)
99004 FORMAT ('# lmp=',i2,'  atype=',i2)
C
99005 FORMAT (1x,a3,'( (',2I3,') -> (',2I3,')',1x,i2,'->',i2,2x,'l=',i2,
     &        1x,'m=',i2,1x,'mq=',i2,1x,'a=',i1,')',2F14.9)
99006 FORMAT (1x,a3,'( (',2I3,') -> (',2I3,')',1x,i2,'->',i2,2x,'l=',i2,
     &        1x,'m=',i2,1x,'mq=',i2,')',4x,2F14.9)
C
      END
C*==kmatpfp.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE KMATPFP(AA,BB,BFIELD,VLM,YLM,YLM1,AXB,AMATX_A,AMATY_A,
     &                   AMATV_A,AMAT_A,AMATX,AMATY,AMATV,AMAT,NATL,
     &                   LAYS)
C     /****************************************************************/
C     # purpose      : calculate the angular matrix elements           *
C                      d(k,mue,ks,mues)   = dmatp                      *
C                      d(-k,mue,-ks,mues) = dmatm                      *
C                      a(k,mue,-ks,mues)  = amatpm                     *
C                      a(-k,mue,ks,mues)  = amatmp                     *
C                      f(k,mue,ks,mues)   = fmatp,                     *
C                      f(-k,mue,-ks,mues) = fmatm                      *
C       note         : 1  dmatp(i,j)  = -dmatm(i,j)                    *
C                      2  due to storage reduction the matrizes are    *
C                         stored as vectors together with index field  *
C                                                                      *
C     # calls the subroutines and functions:                           *
C       clmmgo      clpmgo          blm                                *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,ML,MQD,RSTEP,XMAXE,NFULLPOT,PI,
     &    EPS12,CZERO,CIMAG,NTPHOMAX
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_SITES,ONLY:NOQ
      USE MOD_TRANSPHO_LAYPOT,ONLY:IQ_SPR_LAYAT
      USE MOD_SPEC_RINDC,ONLY:REL
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 BFIELD
      INTEGER LAYS
      COMPLEX*16 AA(3),AMAT(XMAXE,4,NFULLPOT),AMAT_A(XMAXE,4,NFULLPOT,2)
     &           ,AXB(3),VLM(LAYSM,NATLM,MQD,NRMAX,NTPHOMAX),YLM(3),
     &           YLM1(3)
      INTEGER AMATV(XMAXE,4,NFULLPOT),AMATV_A(XMAXE,4,NFULLPOT,2),
     &        AMATX(MQD+1,4,NFULLPOT),AMATX_A(MQD+1,4,NFULLPOT,2),
     &        AMATY(MQD+1,4,NFULLPOT),AMATY_A(MQD+1,4,NFULLPOT,2),
     &        NATL(LAYSM)
      REAL*8 BB(3)
C
C Local variables
C
      INTEGER A,CT(4),CTX(4),CTX_A(4,2),CT_A(4,2),HM11,HM12,HM2,HM31,
     &        HM32,I,IMU,IMUS,IN1,IN2,INDEX1,INDEX2,IQ,J,K,KC,KS,L,L1,
     &        L2,LAST_IN1(4),LAST_IN1_A(4,2),LM,LMP,LS,LSM,M1,M2,MS,
     &        VLM_COM(:)
      COMPLEX*16 AB0,ABM,ABP,B0,BM,BP,CGGA,D0,DM1,DP1,MAT(4),MAT_A(4,2)
      REAL*8 BLM,CGDO,CGUP
      REAL*8 CLM(:),CLP(:),CMM1,CMM2,CMM3,CMM4,CMP1,CMP2,CMP3,CMP4,CPM1,
     &       CPM2,CPM3,CPM4,CPP1,CPP2,CPP3,CPP4,DUM1,DUM2,GA1,GA2,GA3,
     &       GA4,IPMAT,MUE,MUES,PF,RFP,RPF,RPMAT,RTF,X(:,:,:,:)
      INTEGER LM2LMP
      EXTERNAL BLM,CGDO,CGUP,LM2LMP
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE X,CLM,CLP,VLM_COM
      ALLOCATE (X(0:ML-1,-ML-1:ML-1,-1:1,2),CLM(4*ML*ML))
      ALLOCATE (CLP(4*ML*ML),VLM_COM(MQD))
C
      RFP = SQRT(4.D0*PI)/137.036
      PF = 4.D0*PI/(3.D0*137.036)
      RPF = SQRT(3.D0/(4.D0*PI))
      RTF = SQRT(3.D0)/2.D0
C
C     calculate spherical harmonics ylm(l=1; m=0,+-1)(aa)
C     note: y11  = -(x + iy) / sqrt(2)
C           y10  =   z
C           y1-1 =  (x - iy) / sqrt(2)
C
C     some factors are left out (~sqrt(3/4pi)):
C     =========================================
C     all (x+-iy) components should be divided by root(2), see above
      YLM(3) = -(AA(1)+CIMAG*AA(2))/2.D0     !/ sqrt(2.d0)
      YLM(2) = AA(3)/2.D0
      YLM(1) = (AA(1)-CIMAG*AA(2))/2.D0      !/ sqrt(2.d0)
C
Ccghf whats the difference to dm1..d0 beside the factors ?
C     seems these are really the ylm despite a factor root(1/pi)
      YLM1(3) = -(AA(1)+CIMAG*AA(2))*RTF/SQRT(2.D0)
      YLM1(2) = AA(3)*RTF
      YLM1(1) = (AA(1)-CIMAG*AA(2))*RTF/SQRT(2.D0)
C
C     definition of some frequently used constants:
C     =============================================
      DM1 = YLM(3)
      D0 = YLM(2)
      DP1 = YLM(1)
C
C     for non-magnetic case
      AXB(1) = CZERO
      AXB(2) = CZERO
      AXB(3) = CZERO
      BP = CZERO
      BM = CZERO
      B0 = CZERO
      ABM = CZERO
      ABP = CZERO
      AB0 = CZERO
C     for magnetic case
      IF ( BFIELD.GT.EPS12 ) THEN
         BP = (BB(1)+CIMAG*BB(2))     !/ sqrt(2.d0)
         BM = (BB(1)-CIMAG*BB(2))     !/ sqrt(2.d0)
         B0 = BB(3)
C
C     vector-product of aa and bb:
C     ============================
         AXB(1) = AA(2)*BB(3) - AA(3)*BB(2)
         AXB(2) = AA(3)*BB(1) - AA(1)*BB(3)
         AXB(3) = AA(1)*BB(2) - AA(2)*BB(1)
         ABP = (AXB(1)+CIMAG*AXB(2))    !/ sqrt(2.d0)
         ABM = (AXB(1)-CIMAG*AXB(2))    !/ sqrt(2.d0)
         AB0 = AXB(3)
      END IF
C
C     init arrays:
C     ============
      DO I = 1,MQD + 1
         DO J = 1,4
            DO K = 1,NFULLPOT
               AMATX(I,J,K) = 0
               AMATY(I,J,K) = 0
               DO L = 1,2
                  AMATX_A(I,J,K,L) = 0
                  AMATY_A(I,J,K,L) = 0
               END DO
            END DO
         END DO
      END DO
C
      DO I = 1,XMAXE
         DO J = 1,4
            DO K = 1,NFULLPOT
               AMATV(I,J,K) = 0
               AMAT(I,J,K) = CZERO
               DO L = 1,2
                  AMATV_A(I,J,K,L) = 0
                  AMAT_A(I,J,K,L) = CZERO
               END DO
            END DO
         END DO
      END DO
C
C     calculate clebsch-gordon coefficients c(-+kappa,mue,+-1/2):
C     ===========================================================
      CALL CLMMGO(CLM,ML)
      CALL CLPMGO(CLP,ML)
C
C     calculate chi(a,l,m,ms) coefficients:
C     ==================================
      DO I = 0,ML - 1
         DO J = -ML - 1,ML - 1
            DO K = -1,1
               DO L = 1,2
                  X(I,J,K,L) = 0.D0
               END DO
            END DO
         END DO
      END DO
C
      IF ( NFULLPOT.EQ.1 ) THEN
         DO M2 = -1,1
C             chi(0,0,ms, a=1)
            X(0,0,M2,1) = 0.D0
C             chi(0,0,ms, a=2)
            X(0,0,M2,2) = RPF*SQRT(1.D0/3.D0)
         END DO
      ELSE
         DO L1 = 0,ML - 1
            DO M1 = -L1,L1
               DO M2 = -1,1
C             chi(l,m,ms, a=1)
                  DUM1 = DBLE((L1+M1)*(L1-M1)*(2*L1-M2*(M1+M2*(L1+1))))
                  DUM2 = DBLE((2*L1+1)*(2*L1-1)*2*(L1+M1*M2))
C                  IF ( DUM1.EQ.0.D0 ) THEN
                  IF ( ABS(DUM1).LE.1.0D-16 ) THEN
                     X(L1,M1,M2,1) = 0.D0
                  ELSE
                     X(L1,M1,M2,1) = ((-1.D0)**(M2))*RPF*SQRT(DUM1/DUM2)
                  END IF
C             chi(l,m,ms, a=2)
                  DUM1 = DBLE((L1+1+M1)*(L1+1-M1)
     &                   *(2*(L1+1)+M2*(M1-M2*L1)))
                  DUM2 = DBLE((2*L1+1)*(2*L1+3)*2*((L1+1)-M1*M2))
C                  IF ( DUM1.EQ.0.D0 ) THEN
                  IF ( ABS(DUM1).LE.1.0D-16 ) THEN
                     X(L1,M1,M2,2) = 0.D0
                  ELSE
                     X(L1,M1,M2,2) = RPF*SQRT(DUM1/DUM2)
                  END IF
               END DO
            END DO
         END DO
      END IF
C
C     find non zero components of potential and store in vlm_com:
C     ===========================================================
      DO I = 1,MQD
         VLM_COM(I) = 0
      END DO
C
      DO I = 1,LAYS
         DO J = 1,NATL(I)
            IQ = IQ_SPR_LAYAT(J,I)
            DO KC = 1,NOQ(IQ)
               DO K = 1,NFULLPOT
C                  IF ( CDABS(VLM(I,J,REL(K,1),RSTEP-10,KC)).NE.0.0 )
C     &                 VLM_COM(K) = 1
                  IF ( ABS(CDABS(VLM(I,J,REL(K,1),RSTEP-10,KC)))
     &                 .GT.1.0D-16 ) VLM_COM(K) = 1
               END DO
            END DO
         END DO
      END DO
C
C     calculate fp-angular matrix elements:
C     =====================================
C     if (ip.gt.4) write(nout1,9990)
C
      INDEX1 = 0
C
C     for every l1,m1:
C     ================
      DO L1 = 0,ML - 1
         DO M1 = -L1,L1
            LMP = LM2LMP(L1,M1)
C
            IF ( VLM_COM(LMP).NE.0 ) THEN
               DO I = 1,4
                  CT(I) = 0
                  CTX(I) = 0
                  DO J = 1,2
                     CT_A(I,J) = 0
                     CTX_A(I,J) = 0
                  END DO
               END DO
C         for every kappa,mue:
C         ====================
               INDEX1 = 0
               DO K = -ML,ML - 1
                  IF ( K.NE.0 ) THEN
                     DO IMU = -(2*ABS(K)-1),(2*ABS(K)-1),2
                        MUE = DBLE(IMU)*0.5D0
                        INDEX1 = INDEX1 + 2
                        INDEX2 = 0
C             for every kappa',mue':
C             ======================
                        DO KS = -ML,ML - 1
                           IF ( KS.NE.0 ) THEN
                              DO IMUS = -(2*ABS(KS)-1),(2*ABS(KS)-1),2
                                 MUES = DBLE(IMUS)*0.5D0
                                 INDEX2 = INDEX2 + 2
                                 IN1 = INDEX1/2
                                 IN2 = INDEX2/2
C                 some coefficients:
C                 ==================
C
                                 CPP1 = CLP(INDEX1)*CLP(INDEX2-1)
                                 CPP2 = CLP(INDEX1-1)*CLP(INDEX2)
                                 CPP3 = CLP(INDEX1)*CLP(INDEX2)
                                 CPP4 = CLP(INDEX1-1)*CLP(INDEX2-1)
                                 CMM1 = CLM(INDEX1)*CLM(INDEX2-1)
                                 CMM2 = CLM(INDEX1-1)*CLM(INDEX2)
                                 CMM3 = CLM(INDEX1)*CLM(INDEX2)
                                 CMM4 = CLM(INDEX1-1)*CLM(INDEX2-1)
                                 CPM1 = CLP(INDEX1)*CLM(INDEX2-1)
                                 CPM2 = CLP(INDEX1-1)*CLM(INDEX2)
                                 CPM3 = CLP(INDEX1)*CLM(INDEX2)
                                 CPM4 = CLP(INDEX1-1)*CLM(INDEX2-1)
                                 CMP1 = CLM(INDEX1)*CLP(INDEX2-1)
                                 CMP2 = CLM(INDEX1-1)*CLP(INDEX2)
                                 CMP3 = CLM(INDEX1)*CLP(INDEX2)
                                 CMP4 = CLM(INDEX1-1)*CLP(INDEX2-1)
C
C                 calculate corresponding l (l):
C                 ==============================
                                 IF ( K.LT.0 ) THEN
                                    L = -K - 1
                                    LM = -K
                                 ELSE
                                    L = K
                                    LM = K - 1
                                 END IF
                                 IF ( KS.LT.0 ) THEN
                                    LS = -KS - 1
                                    LSM = -KS
                                 ELSE
                                    LS = KS
                                    LSM = KS - 1
                                 END IF
C
                                 DO I = 1,4
                                    DO J = 1,2
                                       MAT_A(I,J) = CZERO
                                    END DO
                                 END DO
                                 DO A = 1,2
                                    L2 = L1 + (-1)**A
                                    IF ( L2.GT.0 ) THEN
                                       DO MS = -1,1
C                             calculate d++:
C                             ==============
                                         HM2 = M1 + MS
C
                                         HM11 = INT(MUE+0.5D0)
                                         HM12 = INT(MUE-0.5D0)
C
                                         HM31 = INT(MUES-0.5D0)
                                         HM32 = INT(MUES+0.5D0)
C
C                             gaunt coefficients
                                         GA1 = BLM(L,HM11,L2,HM2,LS,
     &                                      HM31)
                                         GA2 = BLM(L,HM12,L2,HM2,LS,
     &                                      HM32)
C
                                         CPP1 = CGDO(K,MUE)
     &                                      *CGUP(KS,MUES)
                                         CPP2 = CGUP(K,MUE)
     &                                      *CGDO(KS,MUES)
C
                                         CGGA = PF*YLM1(MS+2)
     &                                      *X(L1,M1,MS,A)
     &                                      *(CPP1*GA1-CPP2*GA2)
C
                                         MAT_A(1,A) = MAT_A(1,A) + CGGA
C
C                             calculate d-- (d-- = -d++):
C                             ===========================
                                         MAT_A(2,A) = -MAT_A(1,A)
C
C                             for the magnetic case:
C                             ======================
C                             calculate f++:
C                             ==============
                                         GA3 = BLM(L,HM11,L2,HM2,LS,
     &                                      HM32)
                                         GA4 = BLM(L,HM12,L2,HM2,LS,
     &                                      HM31)
                                         MAT_A(3,A) = MAT_A(3,A)
     &                                      - PF*YLM1(MS+2)
     &                                      *X(L1,M1,MS,A)
     &                                      *(CPP4*BP*GA4-CPP3*BM*GA3-
     &                                      CPP2*B0*GA2-CPP1*B0*GA1)
C                             calculate f--:
C                             ==============
                                         GA1 = BLM(LM,HM11,L2,HM2,LSM,
     &                                      HM31)
                                         GA2 = BLM(LM,HM12,L2,HM2,LSM,
     &                                      HM32)
                                         GA3 = BLM(LM,HM11,L2,HM2,LSM,
     &                                      HM32)
                                         GA4 = BLM(LM,HM12,L2,HM2,LSM,
     &                                      HM31)
                                         MAT_A(4,A) = MAT_A(4,A)
     &                                      - PF*YLM1(MS+2)
     &                                      *X(L1,M1,MS,A)
     &                                      *(CMM4*BP*GA4-CMM3*BM*GA3-
     &                                      CMM2*B0*GA2-CMM1*B0*GA1)
                                       END DO
                                    END IF
                                 END DO
C
C
C                 calculate a+-:
C                 ==============
                                 MAT(1)
     &                              = -RFP*(CPM4*DM1*BLM(L,HM12,L1,M1,
     &                              LSM,HM31)
     &                              -CPM2*D0*BLM(L,HM12,L1,M1,LSM,HM32)
     &                              -CPM3*DP1*BLM(L,HM11,L1,M1,LSM,HM32)
     &                              -CPM1*D0*BLM(L,HM11,L1,M1,LSM,HM31))
C                 calculate a-+:
C                 ==============
                                 MAT(2)
     &                              = -RFP*(CMP4*DM1*BLM(LM,HM12,L1,M1,
     &                              LS,HM31)
     &                              -CMP2*D0*BLM(LM,HM12,L1,M1,LS,HM32)
     &                              -CMP3*DP1*BLM(LM,HM11,L1,M1,LS,HM32)
     &                              -CMP1*D0*BLM(LM,HM11,L1,M1,LS,HM31))
C
C                 for the magnetic case:
C                 ======================
C                 calculate e++:
C                 ==============
                                 MAT(3)
     &                              = -RFP*(CPP4*ABP*BLM(L,HM12,L1,M1,
     &                              LS,HM31)
     &                              -CPP2*AB0*BLM(L,HM12,L1,M1,LS,HM32)
     &                              -CPP3*ABM*BLM(L,HM11,L1,M1,LS,HM32)
     &                              -CPP1*AB0*BLM(L,HM11,L1,M1,LS,HM31))
C                 calculate e--:
C                 ==============
                                 MAT(4)
     &                              = -RFP*(CMM4*ABP*BLM(LM,HM12,L1,M1,
     &                              LSM,HM31)
     &                              -CMM2*AB0*BLM(LM,HM12,L1,M1,LSM,
     &                              HM32)
     &                              -CMM3*ABM*BLM(LM,HM11,L1,M1,LSM,
     &                              HM32)
     &                              -CMM1*AB0*BLM(LM,HM11,L1,M1,LSM,
     &                              HM31))
C
C                 store fp-angular matrixelements in vector-fields:
C                 =================================================
C                 store d++,d--,f++,f-- for a=1,2:
C                 ================================
                                 DO I = 1,4
                                    DO A = 1,2
                                       IF ( CDABS(MAT_A(I,A)).GT.EPS12 )
     &                                    THEN
                                         CT_A(I,A) = CT_A(I,A) + 1
                                         IF ( LAST_IN1_A(I,A).NE.IN1 )
     &                                      THEN
                                         LAST_IN1_A(I,A) = IN1
                                         CTX_A(I,A) = CTX_A(I,A) + 1
                                         AMATX_A(CTX_A(I,A),I,LMP,A)
     &                                      = IN1
                                         AMATY_A(CTX_A(I,A),I,LMP,A)
     &                                      = CT_A(I,A)
                                         END IF
                                         AMATV_A(CT_A(I,A),I,LMP,A)
     &                                      = IN2
C                             amat_a(ct_a(i,a),i,lmp,a)  = mat_a(i,a)
                                         RPMAT = DBLE(MAT_A(I,A))
                                         IPMAT = DIMAG(MAT_A(I,A))
                                         IF ( ABS(RPMAT).LT.EPS12 )
     &                                      RPMAT = 0.D0
                                         IF ( ABS(IPMAT).LT.EPS12 )
     &                                      IPMAT = 0.D0
C
                                         AMAT_A(CT_A(I,A),I,LMP,A)
     &                                      = DCMPLX(RPMAT,IPMAT)
C        if (i.eq.1) write (9,1111) 'ang',k,ks,mue,mues,a,
C    1               ct_a(i,a),mat_a(1,a)
C1111    format(1x,2i4,2x,2f5.2,2x,2i3,2x,2e14.7)
                                       END IF
                                    END DO
                                 END DO
C
C                 store a+-,a-+,e++,e--:
C                 ======================
                                 DO I = 1,4
                                    IF ( CDABS(MAT(I)).GT.EPS12 ) THEN
                                       CT(I) = CT(I) + 1
                                       IF ( LAST_IN1(I).NE.IN1 ) THEN
                                         LAST_IN1(I) = IN1
                                         CTX(I) = CTX(I) + 1
                                         AMATX(CTX(I),I,LMP) = IN1
                                         AMATY(CTX(I),I,LMP) = CT(I)
                                       END IF
                                       AMATV(CT(I),I,LMP) = IN2
                                       RPMAT = DBLE(MAT(I))
                                       IPMAT = DIMAG(MAT(I))
                                       IF ( ABS(RPMAT).LT.EPS12 )
     &                                    RPMAT = 0.D0
                                       IF ( ABS(IPMAT).LT.EPS12 )
     &                                    IPMAT = 0.D0
                                       AMAT(CT(I),I,LMP)
     &                                    = DCMPLX(RPMAT,IPMAT)
C        if (i.eq.1) write (9,*) 'ang',k,ks,mue,mues,
C    1               dcmplx(rpmat,ipmat)
C1111    format(1x,2i4,2x,2f5.2,2x,2e14.7)
                                    END IF
                                 END DO
C
C             enddo mue':
C             ===========
                              END DO
                           END IF
C             enddo k':
C             =========
                        END DO
C
                        DO I = 1,4
                           DO A = 1,2
C                     amatv_a(mqd+1,i,lmp,a) = ct_a(i,a)+1
                              AMATY_A(MQD+1,I,LMP,A) = CT_A(I,A) + 1
                           END DO
C                 amatv(mqd+1,i,lmp) = ct(i)+1
                           AMATY(MQD+1,I,LMP) = CT(I) + 1
                        END DO
C
C         enddo mue:
C         ==========
                     END DO
                     DO I = 1,4
                        DO A = 1,2
                           AMATX_A(MQD+1,I,LMP,A) = CTX_A(I,A)
                           AMATY_A(CTX_A(I,A)+1,I,LMP,A) = CT_A(I,A) + 1
                        END DO
                        AMATX(MQD+1,I,LMP) = CTX(I)
                        AMATY(CTX(I)+1,I,LMP) = CT(I) + 1
                     END DO
                  END IF
C         enddo k:
C         ========
               END DO
            END IF
C     enddo m1:
C     =========
         END DO
C     enddo l1:
C     =========
      END DO
C
C
      END
C*==kmatpfpt.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE KMATPFPT(AA,BB,BFIELD,VLM,YLM,YLM1,AXB,AMATX_A,AMATY_A,
     &                    AMATV_A,AMAT_A,AMATX,AMATY,AMATV,AMAT,NATL,
     &                    LAYS)
C     /****************************************************************/
C     # purpose      : calculate the angular matrix elements           *
C                      d(k,mue,ks,mues)   = dmatp                      *
C                      d(-k,mue,-ks,mues) = dmatm                      *
C                      a(k,mue,-ks,mues)  = amatpm                     *
C                      a(-k,mue,ks,mues)  = amatmp                     *
C                      f(k,mue,ks,mues)   = fmatp,                     *
C                      f(-k,mue,-ks,mues) = fmatm                      *
C       note         : 1  dmatp(i,j)  = -dmatm(i,j)                    *
C                      2  due to storage reduction the matrizes are    *
C                         stored as vectors together with index field  *
C                                                                      *
C     # calls the subroutines and functions:                           *
C       clmmgo      clpmgo          blm                                *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,ML,MQD,RSTEP,XMAXE,NFULLPOT,PI,
     &    EPS12,CZERO,CIMAG,NTPHOMAX
      USE MOD_SITES,ONLY:NOQ
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_TRANSPHO_LAYPOT,ONLY:IQ_SPR_LAYAT
      USE MOD_SPEC_RINDC,ONLY:REL
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 BFIELD
      INTEGER LAYS
      COMPLEX*16 AA(3),AMAT(XMAXE,4,NFULLPOT),AMAT_A(XMAXE,4,NFULLPOT,2)
     &           ,AXB(3),VLM(LAYSM,NATLM,MQD,NRMAX,NTPHOMAX),YLM(3),
     &           YLM1(3)
      INTEGER AMATV(XMAXE,4,NFULLPOT),AMATV_A(XMAXE,4,NFULLPOT,2),
     &        AMATX(MQD+1,4,NFULLPOT),AMATX_A(MQD+1,4,NFULLPOT,2),
     &        AMATY(MQD+1,4,NFULLPOT),AMATY_A(MQD+1,4,NFULLPOT,2),
     &        NATL(LAYSM)
      REAL*8 BB(3)
C
C Local variables
C
      INTEGER A,CT(4),CTX(4),CTX_A(4,2),CT_A(4,2),HM11,HM12,HM2,HM31,
     &        HM32,I,IMU,IMUS,IN1,IN2,INDEX1,INDEX2,IQ,J,K,KC,KS,L,L1,
     &        L2,LAST_IN1(4),LAST_IN1_A(4,2),LM,LMP,LS,LSM,M1,M2,MS,
     &        VLM_COM(:)
      COMPLEX*16 AB0,ABM,ABP,B0,BM,BP,CGGA,D0,DM1,DP1,MAT(4),MAT_A(4,2)
      REAL*8 BLM,CGDO,CGUP
      REAL*8 CLM(:),CLP(:),CMM1,CMM2,CMM3,CMM4,CMP1,CMP2,CMP3,CMP4,CPM1,
     &       CPM2,CPM3,CPM4,CPP1,CPP2,CPP3,CPP4,DUM1,DUM2,GA1,GA2,GA3,
     &       GA4,IPMAT,MUE,MUES,PF,RFP,RPF,RPMAT,RTF,X(:,:,:,:)
      INTEGER LM2LMP
      EXTERNAL BLM,CGDO,CGUP,LM2LMP
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE X,CLM,CLP,VLM_COM
      ALLOCATE (X(0:ML-1,-ML-1:ML-1,-1:1,2),CLM(4*ML*ML))
      ALLOCATE (CLP(4*ML*ML),VLM_COM(MQD))
C
C*** End of declarations rewritten by SPAG
C
      RFP = SQRT(4.D0*PI)/137.036
      PF = 4.D0*PI/(3.D0*137.036)
      RPF = SQRT(3.D0/(4.D0*PI))
      RTF = SQRT(3.D0)/2.D0
C
C     calculate spherical harmonics ylm(l=1; m=0,+-1)(aa)
C     note: y11  = -(x + iy) / sqrt(2)
C           y10  =   z
C           y1-1 =  (x - iy) / sqrt(2)
C
C     some factors are left out (~sqrt(3/4pi)):
C     =========================================
C     all (x+-iy) components should be divided by root(2), see above
      YLM(3) = -(AA(1)+CIMAG*AA(2))/2.D0     !/ sqrt(2.d0)
      YLM(2) = AA(3)/2.D0
      YLM(1) = (AA(1)-CIMAG*AA(2))/2.D0      !/ sqrt(2.d0)
C
Ccghf whats the difference to dm1..d0 beside the factors ?
C     seems these are really the ylm despite a factor root(1/pi)
      YLM1(3) = -(AA(1)+CIMAG*AA(2))*RTF/SQRT(2.D0)
      YLM1(2) = AA(3)*RTF
      YLM1(1) = (AA(1)-CIMAG*AA(2))*RTF/SQRT(2.D0)
C
C     definition of some frequently used constants:
C     =============================================
      DM1 = YLM(3)
      D0 = YLM(2)
      DP1 = YLM(1)
C
C     for non-magnetic case
      AXB(1) = CZERO
      AXB(2) = CZERO
      AXB(3) = CZERO
      BP = CZERO
      BM = CZERO
      B0 = CZERO
      ABM = CZERO
      ABP = CZERO
      AB0 = CZERO
C     for magnetic case
      IF ( BFIELD.GT.EPS12 ) THEN
         BP = (BB(1)+CIMAG*BB(2))     !/ sqrt(2.d0)
         BM = (BB(1)-CIMAG*BB(2))     !/ sqrt(2.d0)
         B0 = BB(3)
C
C     vector-product of aa and bb:
C     ============================
         AXB(1) = AA(2)*BB(3) - AA(3)*BB(2)
         AXB(2) = AA(3)*BB(1) - AA(1)*BB(3)
         AXB(3) = AA(1)*BB(2) - AA(2)*BB(1)
         ABP = (AXB(1)+CIMAG*AXB(2))    !/ sqrt(2.d0)
         ABM = (AXB(1)-CIMAG*AXB(2))    !/ sqrt(2.d0)
         AB0 = AXB(3)
      END IF
C
C     init arrays:
C     ============
      DO I = 1,MQD + 1
         DO J = 1,4
            DO K = 1,NFULLPOT
               AMATX(I,J,K) = 0
               AMATY(I,J,K) = 0
               DO L = 1,2
                  AMATX_A(I,J,K,L) = 0
                  AMATY_A(I,J,K,L) = 0
               END DO
            END DO
         END DO
      END DO
C
      DO I = 1,XMAXE
         DO J = 1,4
            DO K = 1,NFULLPOT
               AMATV(I,J,K) = 0
               AMAT(I,J,K) = CZERO
               DO L = 1,2
                  AMATV_A(I,J,K,L) = 0
                  AMAT_A(I,J,K,L) = CZERO
               END DO
            END DO
         END DO
      END DO
C
C     calculate clebsch-gordon coefficients c(-+kappa,mue,+-1/2):
C     ===========================================================
      CALL CLMMGO(CLM,ML)
      CALL CLPMGO(CLP,ML)
C
C     calculate chi(a,l,m,ms) coefficients:
C     ==================================
      DO I = 0,ML - 1
         DO J = -ML - 1,ML - 1
            DO K = -1,1
               DO L = 1,2
                  X(I,J,K,L) = 0.D0
               END DO
            END DO
         END DO
      END DO
C
      IF ( NFULLPOT.EQ.1 ) THEN
         DO M2 = -1,1
C             chi(0,0,ms, a=1)
            X(0,0,M2,1) = 0.D0
C             chi(0,0,ms, a=2)
            X(0,0,M2,2) = RPF*SQRT(1.D0/3.D0)
         END DO
      ELSE
         DO L1 = 0,ML - 1
            DO M1 = -L1,L1
               DO M2 = -1,1
C             chi(l,m,ms, a=1)
                  DUM1 = DBLE((L1+M1)*(L1-M1)*(2*L1-M2*(M1+M2*(L1+1))))
                  DUM2 = DBLE((2*L1+1)*(2*L1-1)*2*(L1+M1*M2))
C                  IF ( DUM1.EQ.0.D0 ) THEN
                  IF ( ABS(DUM1).LE.1.0D-16 ) THEN
                     X(L1,M1,M2,1) = 0.D0
                  ELSE
                     X(L1,M1,M2,1) = ((-1.D0)**(M2))*RPF*SQRT(DUM1/DUM2)
                  END IF
C             chi(l,m,ms, a=2)
                  DUM1 = DBLE((L1+1+M1)*(L1+1-M1)
     &                   *(2*(L1+1)+M2*(M1-M2*L1)))
                  DUM2 = DBLE((2*L1+1)*(2*L1+3)*2*((L1+1)-M1*M2))
C                  IF ( DUM1.EQ.0.D0 ) THEN
                  IF ( ABS(DUM1).LE.1.0D-16 ) THEN
                     X(L1,M1,M2,2) = 0.D0
                  ELSE
                     X(L1,M1,M2,2) = RPF*SQRT(DUM1/DUM2)
                  END IF
               END DO
            END DO
         END DO
      END IF
C
C     find non zero components of potential and store in vlm_com:
C     ===========================================================
      DO I = 1,MQD
         VLM_COM(I) = 0
      END DO
C
      DO I = 1,LAYS
         DO J = 1,NATL(I)
            IQ = IQ_SPR_LAYAT(J,I)
            DO KC = 1,NOQ(IQ)
               DO K = 1,NFULLPOT
C                  IF ( CDABS(VLM(I,J,REL(K,1),RSTEP-10,KC)).NE.0.0 )
C     &                 VLM_COM(K) = 1
                  IF ( ABS(CDABS(VLM(I,J,REL(K,1),RSTEP-10,KC)))
     &                 .GT.1.0D-16 ) VLM_COM(K) = 1
               END DO
            END DO
         END DO
      END DO
C
C     calculate fp-angular matrix elements:
C     =====================================
C     if (ip.gt.4) write(nout1,9990)
C
      INDEX1 = 0
C
C     for every l1,m1:
C     ================
      DO L1 = 0,ML - 1
         DO M1 = -L1,L1
            LMP = LM2LMP(L1,M1)
C
            IF ( VLM_COM(LMP).NE.0 ) THEN
               DO I = 1,4
                  CT(I) = 0
                  CTX(I) = 0
                  DO J = 1,2
                     CT_A(I,J) = 0
                     CTX_A(I,J) = 0
                  END DO
               END DO
C         for every kappas,mues:
C         ====================
               INDEX1 = 0
               DO KS = -ML,ML - 1
                  IF ( KS.NE.0 ) THEN
                     DO IMUS = -(2*ABS(KS)-1),(2*ABS(KS)-1),2
                        MUES = DBLE(IMUS)*0.5D0
                        INDEX1 = INDEX1 + 2
                        INDEX2 = 0
C             for every kappa',mue':
C             ======================
                        DO K = -ML,ML - 1
                           IF ( K.NE.0 ) THEN
                              DO IMU = -(2*ABS(K)-1),(2*ABS(K)-1),2
                                 MUE = DBLE(IMU)*0.5D0
                                 INDEX2 = INDEX2 + 2
                                 IN1 = INDEX1/2
                                 IN2 = INDEX2/2
C                 some coefficients:
C                 ==================
C
                                 CPP1 = CLP(INDEX1)*CLP(INDEX2-1)
                                 CPP2 = CLP(INDEX1-1)*CLP(INDEX2)
                                 CPP3 = CLP(INDEX1)*CLP(INDEX2)
                                 CPP4 = CLP(INDEX1-1)*CLP(INDEX2-1)
                                 CMM1 = CLM(INDEX1)*CLM(INDEX2-1)
                                 CMM2 = CLM(INDEX1-1)*CLM(INDEX2)
                                 CMM3 = CLM(INDEX1)*CLM(INDEX2)
                                 CMM4 = CLM(INDEX1-1)*CLM(INDEX2-1)
                                 CPM1 = CLP(INDEX1)*CLM(INDEX2-1)
                                 CPM2 = CLP(INDEX1-1)*CLM(INDEX2)
                                 CPM3 = CLP(INDEX1)*CLM(INDEX2)
                                 CPM4 = CLP(INDEX1-1)*CLM(INDEX2-1)
                                 CMP1 = CLM(INDEX1)*CLP(INDEX2-1)
                                 CMP2 = CLM(INDEX1-1)*CLP(INDEX2)
                                 CMP3 = CLM(INDEX1)*CLP(INDEX2)
                                 CMP4 = CLM(INDEX1-1)*CLP(INDEX2-1)
C
C                 calculate corresponding l (l):
C                 ==============================
                                 IF ( K.LT.0 ) THEN
                                    L = -K - 1
                                    LM = -K
                                 ELSE
                                    L = K
                                    LM = K - 1
                                 END IF
                                 IF ( KS.LT.0 ) THEN
                                    LS = -KS - 1
                                    LSM = -KS
                                 ELSE
                                    LS = KS
                                    LSM = KS - 1
                                 END IF
C
                                 DO I = 1,4
                                    DO J = 1,2
                                       MAT_A(I,J) = CZERO
                                    END DO
                                 END DO
                                 DO A = 1,2
                                    L2 = L1 + (-1)**A
                                    IF ( L2.GT.0 ) THEN
                                       DO MS = -1,1
C                             calculate d++:
C                             ==============
                                         HM2 = M1 + MS
C
                                         HM11 = INT(MUE+0.5D0)
                                         HM12 = INT(MUE-0.5D0)
C
                                         HM31 = INT(MUES-0.5D0)
                                         HM32 = INT(MUES+0.5D0)
C
C                             gaunt coefficients
                                         GA1 = BLM(L,HM11,L2,HM2,LS,
     &                                      HM31)
                                         GA2 = BLM(L,HM12,L2,HM2,LS,
     &                                      HM32)
C
                                         CPP1 = CGDO(K,MUE)
     &                                      *CGUP(KS,MUES)
                                         CPP2 = CGUP(K,MUE)
     &                                      *CGDO(KS,MUES)
C
                                         CGGA = PF*YLM1(MS+2)
     &                                      *X(L1,M1,MS,A)
     &                                      *(CPP1*GA1-CPP2*GA2)
C
                                         MAT_A(1,A) = MAT_A(1,A) + CGGA
C
C                             calculate d-- (d-- = -d++):
C                             ===========================
                                         MAT_A(2,A) = -MAT_A(1,A)
C
C                             for the magnetic case:
C                             ======================
C                             calculate f++:
C                             ==============
                                         GA3 = BLM(L,HM11,L2,HM2,LS,
     &                                      HM32)
                                         GA4 = BLM(L,HM12,L2,HM2,LS,
     &                                      HM31)
                                         MAT_A(3,A) = MAT_A(3,A)
     &                                      - PF*YLM1(MS+2)
     &                                      *X(L1,M1,MS,A)
     &                                      *(CPP4*BP*GA4-CPP3*BM*GA3-
     &                                      CPP2*B0*GA2-CPP1*B0*GA1)
C                             calculate f--:
C                             ==============
                                         GA1 = BLM(LM,HM11,L2,HM2,LSM,
     &                                      HM31)
                                         GA2 = BLM(LM,HM12,L2,HM2,LSM,
     &                                      HM32)
                                         GA3 = BLM(LM,HM11,L2,HM2,LSM,
     &                                      HM32)
                                         GA4 = BLM(LM,HM12,L2,HM2,LSM,
     &                                      HM31)
                                         MAT_A(4,A) = MAT_A(4,A)
     &                                      - PF*YLM1(MS+2)
     &                                      *X(L1,M1,MS,A)
     &                                      *(CMM4*BP*GA4-CMM3*BM*GA3-
     &                                      CMM2*B0*GA2-CMM1*B0*GA1)
                                       END DO
                                    END IF
                                 END DO
C
C
C                 calculate a+-:
C                 ==============
                                 MAT(1)
     &                              = -RFP*(CPM4*DM1*BLM(L,HM12,L1,M1,
     &                              LSM,HM31)
     &                              -CPM2*D0*BLM(L,HM12,L1,M1,LSM,HM32)
     &                              -CPM3*DP1*BLM(L,HM11,L1,M1,LSM,HM32)
     &                              -CPM1*D0*BLM(L,HM11,L1,M1,LSM,HM31))
C                 calculate a-+:
C                 ==============
                                 MAT(2)
     &                              = -RFP*(CMP4*DM1*BLM(LM,HM12,L1,M1,
     &                              LS,HM31)
     &                              -CMP2*D0*BLM(LM,HM12,L1,M1,LS,HM32)
     &                              -CMP3*DP1*BLM(LM,HM11,L1,M1,LS,HM32)
     &                              -CMP1*D0*BLM(LM,HM11,L1,M1,LS,HM31))
C
C                 for the magnetic case:
C                 ======================
C                 calculate e++:
C                 ==============
                                 MAT(3)
     &                              = -RFP*(CPP4*ABP*BLM(L,HM12,L1,M1,
     &                              LS,HM31)
     &                              -CPP2*AB0*BLM(L,HM12,L1,M1,LS,HM32)
     &                              -CPP3*ABM*BLM(L,HM11,L1,M1,LS,HM32)
     &                              -CPP1*AB0*BLM(L,HM11,L1,M1,LS,HM31))
C                 calculate e--:
C                 ==============
                                 MAT(4)
     &                              = -RFP*(CMM4*ABP*BLM(LM,HM12,L1,M1,
     &                              LSM,HM31)
     &                              -CMM2*AB0*BLM(LM,HM12,L1,M1,LSM,
     &                              HM32)
     &                              -CMM3*ABM*BLM(LM,HM11,L1,M1,LSM,
     &                              HM32)
     &                              -CMM1*AB0*BLM(LM,HM11,L1,M1,LSM,
     &                              HM31))
C
C                 store fp-angular matrixelements in vector-fields:
C                 =================================================
C                 store d++,d--,f++,f-- for a=1,2:
C                 ================================
                                 DO I = 1,4
                                    DO A = 1,2
                                       IF ( CDABS(MAT_A(I,A)).GT.EPS12 )
     &                                    THEN
                                         CT_A(I,A) = CT_A(I,A) + 1
                                         IF ( LAST_IN1_A(I,A).NE.IN1 )
     &                                      THEN
                                         LAST_IN1_A(I,A) = IN1
                                         CTX_A(I,A) = CTX_A(I,A) + 1
                                         AMATX_A(CTX_A(I,A),I,LMP,A)
     &                                      = IN1
                                         AMATY_A(CTX_A(I,A),I,LMP,A)
     &                                      = CT_A(I,A)
                                         END IF
                                         AMATV_A(CT_A(I,A),I,LMP,A)
     &                                      = IN2
C                             amat_a(ct_a(i,a),i,lmp,a)  = mat_a(i,a)
                                         RPMAT = DBLE(MAT_A(I,A))
                                         IPMAT = DIMAG(MAT_A(I,A))
                                         IF ( ABS(RPMAT).LT.EPS12 )
     &                                      RPMAT = 0.D0
                                         IF ( ABS(IPMAT).LT.EPS12 )
     &                                      IPMAT = 0.D0
C
                                         AMAT_A(CT_A(I,A),I,LMP,A)
     &                                      = -DCMPLX(RPMAT,IPMAT)
C      if (i.eq.1) write (8,1111) 'angT',k,ks,mue,mues,a,
C    1             ct_a(i,a),-mat_a(1,a)
C1111  format(1x,2i4,2x,2f5.2,2x,2i3,2x,2e14.7)
                                       END IF
                                    END DO
                                 END DO
C
C                 store a+-,a-+,e++,e--:
C                 ======================
                                 DO I = 1,4
                                    IF ( CDABS(MAT(I)).GT.EPS12 ) THEN
                                       CT(I) = CT(I) + 1
                                       IF ( LAST_IN1(I).NE.IN1 ) THEN
                                         LAST_IN1(I) = IN1
                                         CTX(I) = CTX(I) + 1
                                         AMATX(CTX(I),I,LMP) = IN1
                                         AMATY(CTX(I),I,LMP) = CT(I)
                                       END IF
                                       AMATV(CT(I),I,LMP) = IN2
                                       RPMAT = DBLE(MAT(I))
                                       IPMAT = DIMAG(MAT(I))
                                       IF ( ABS(RPMAT).LT.EPS12 )
     &                                    RPMAT = 0.D0
                                       IF ( ABS(IPMAT).LT.EPS12 )
     &                                    IPMAT = 0.D0
                                       AMAT(CT(I),I,LMP)
     &                                    = -DCMPLX(RPMAT,IPMAT)
                                    END IF
                                 END DO
C
C             enddo mue':
C             ===========
                              END DO
                           END IF
C             enddo k':
C             =========
                        END DO
C
                        DO I = 1,4
                           DO A = 1,2
C                     amatv_a(mqd+1,i,lmp,a) = ct_a(i,a)+1
                              AMATY_A(MQD+1,I,LMP,A) = CT_A(I,A) + 1
                           END DO
C                 amatv(mqd+1,i,lmp) = ct(i)+1
                           AMATY(MQD+1,I,LMP) = CT(I) + 1
                        END DO
C
C         enddo mue:
C         ==========
                     END DO
                     DO I = 1,4
                        DO A = 1,2
                           AMATX_A(MQD+1,I,LMP,A) = CTX_A(I,A)
                           AMATY_A(CTX_A(I,A)+1,I,LMP,A) = CT_A(I,A) + 1
                        END DO
                        AMATX(MQD+1,I,LMP) = CTX(I)
                        AMATY(CTX(I)+1,I,LMP) = CT(I) + 1
                     END DO
                  END IF
C         enddo k:
C         ========
               END DO
            END IF
C     enddo m1:
C     =========
         END DO
C     enddo l1:
C     =========
      END DO
C
C
      END
C*==clpmgo.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C
      SUBROUTINE CLPMGO(CLP,NML)
C     /****************************************************************/
C     # purpose       : clebsch gordon coefficients for large component*
C                       c(kappa,mue,+-1/2)                             *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:MLSP
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NML
      REAL*8 CLP(MLSP)
C
C Local variables
C
      INTEGER II1,INDEX1,K,KAP1,S
      REAL*8 MUE
C
C*** End of declarations rewritten by SPAG
C
      INDEX1 = 0
      DO K = -NML,NML - 1
         IF ( K.NE.0 ) THEN
C           DO MUE = -(ABS(K)-0.5D0),ABS(K) - 0.5D0
            KAP1 = 2*ABS(K)
            DO II1 = 1,KAP1
               MUE = -ABS(K) + 0.5D0 + (II1*1.D0-1.0D0)
               DO S = 1, - 1, - 2
                  INDEX1 = INDEX1 + 1
                  IF ( S.EQ.-1 ) THEN
                     CLP(INDEX1) = SQRT((K+MUE+0.5D0)/DBLE(2*K+1))
                  ELSE
                     CLP(INDEX1) = (-K/ABS(K))
     &                             *SQRT((K-MUE+0.5D0)/DBLE(2*K+1))
                  END IF
               END DO
            END DO
         END IF
      END DO
C
      END
C*==clmmgo.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE CLMMGO(CLM,NML)
C     /****************************************************************/
C     # purpose:        clebsch gordon coefficients for small component*
C                       c(-kappa,mue,+-1/2)                            *
C     /****************************************************************/
CC
      USE MOD_SPEC,ONLY:MLSP
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NML
      REAL*8 CLM(MLSP)
C
C Local variables
C
      INTEGER II1,INDEX1,K,KAP1,S
      REAL*8 MUE
C
C*** End of declarations rewritten by SPAG
C
      INDEX1 = 0
      DO K = -NML,NML - 1
         IF ( K.NE.0 ) THEN
C
C           DO MUE = -(ABS(K)-0.5D0),ABS(K) - 0.5D0
C
            KAP1 = 2*ABS(K)
            DO II1 = 1,KAP1
               MUE = -ABS(K) + 0.5D0 + (II1*1.D0-1.0D0)
C
               DO S = 1, - 1, - 2
                  INDEX1 = INDEX1 + 1
                  IF ( S.EQ.-1 ) THEN
                     CLM(INDEX1) = SQRT((K-MUE-0.5D0)/DBLE(2*K-1))
                  ELSE
                     CLM(INDEX1) = (K/ABS(K))
     &                             *SQRT((K+MUE-0.5D0)/DBLE(2*K-1))
                  END IF
               END DO
            END DO
         END IF
      END DO
C
      END
C*==lm2lmp.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      INTEGER FUNCTION LM2LMP(L,M)
C     /****************************************************************/
C     # purpose       : return index into vlm (potentials) for given   *
C                       l and m quantumnumber                          *
C     /****************************************************************/
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER L,M
C
C*** End of declarations rewritten by SPAG
C
      LM2LMP = L**2 + 1 + L + M
      END
C*==rdipolem.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE RDIPOLEM(LAY,ATOM,VLM,RDIP1M,RDIP2M,CLIGHT,OMHAR,ALPHA,
     &                    IO)
C
C     # purpose      : determine radial parts of dipole operator       *
C                      for cpa-typ nq                                  *
C                      in acceleration representation for munich mesh  *
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MQD,RSTEP,NFULLPOT,EPS12,CZERO,
     &    NTPHOMAX
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_SPEC_MESH,ONLY:RADMESHM,RADMESH3M
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALPHA,CLIGHT,OMHAR
      INTEGER ATOM,IO,LAY
      COMPLEX*16 RDIP1M(NFULLPOT,2,NRMAX,2,NTPHOMAX),
     &           RDIP2M(NFULLPOT,NRMAX,2,NTPHOMAX),
     &           VLM(LAYSM,NATLM,MQD,NRMAX,NTPHOMAX)
C
C Local variables
C
      COMPLEX*16 B(9),CB,CV,DB,DB1H,DB2H,DB4H,DV,DV1H,DV2H,DV4H,DVLMB,
     &           DVLMV,V(9)
      INTEGER I,J,K,L,LMP,M,N,RL
      REAL*8 RAD3M(:),RADM(:)
C
C*** End of declarations rewritten by SPAG
C
C
C
C     # note:                                                          *
C         rdip1m(lm,a=1,r,1,nq) = fac*(d/dr + (l+1)/r) v(r)            *
C         rdip1m(lm,a=2,r,1,nq) = fac*(d/dr -   l/r  ) v(r)            *
C         rdip1m(lm,a=1,r,2,nq) = fac*(d/dr + (l+1)/r) b(r)            *
C         rdip1m(lm,a=2,r,2,nq) = fac*(d/dr -   l/r  ) b(r)            *
C         rdip2m(lm,r,1,nq) = fac * omega/c * v(r)                     *
C         rdip2m(lm,r,2,nq) = fac * omega/c * b(r)                     *
C                                                                      *
C         fac = -2*i*c/((2*e+omega)*omega)                             *
C         (this factor will be multiplied later in the program)        *
C                                                                      *
C         the derivative of v(r) and b(r) is determined                *
C         by a richardson extrapolation up to the order of o**6        *
C         since the derivative of v(r) and b(r) is very                *
C         steep, the derivative of r*v(r) and r*b(r) is                *
C         calculated; the derivatives of v(r) and b(r)                 *
C         are then given by:                                           *
C         v'(r) = 1/r*( (r*v(r))' - v(r))                              *
C         b'(r) = 1/r*( (r*b(r))' - b(r))                              *
C                                                                      *
C     # calls the subroutines:   detvb       ind                       *
C
C
      ALLOCATABLE radm,rad3m
      ALLOCATE (RADM(NRMAX),RAD3M(NRMAX))
C
      DO I = 1,NFULLPOT
         DO J = 1,2
            DO K = 1,NRMAX
               DO N = 1,IO
                  RDIP2M(I,K,J,N) = CZERO
                  DO L = 1,2
                     RDIP1M(I,J,K,L,N) = CZERO
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      DO I = 1,NRMAX
         RADM(I) = RADMESHM(I,ATOM,LAY)
         RAD3M(I) = RADMESH3M(I,ATOM,LAY)
      END DO
C     determine sum v and b:
C
      CALL DETVBM(LAY,ATOM,VLM,RDIP2M,IO)
C
      DO LMP = 1,NFULLPOT
         IF ( CDABS(RDIP2M(LMP,RSTEP-20,1,IO)).GT.EPS12 ) THEN
            DO RL = 1,RSTEP
C              calculate the 8 values which are necessary for
C              determining the derivative:
C
               DO I = 1,9
                  IF ( (RL-5+I).LE.0 ) THEN
                     IF ( LMP.EQ.1 .AND. DBLE(RDIP2M(LMP,1,1,IO))
     &                    .LT.(-100.D0) ) THEN
                        V(I) = RDIP2M(LMP,1,1,IO)*RADM(1)
                        B(I) = RDIP2M(LMP,1,2,IO)*RADM(1)
                     ELSE
                        V(I) = RDIP2M(LMP,1,1,IO)*RADM(RSTEP)
     &                         *DEXP(DBLE(RL-5+I-RSTEP)*ALPHA)
                        B(I) = RDIP2M(LMP,1,2,IO)*RADM(RSTEP)
     &                         *DEXP(DBLE(RL-5+I-RSTEP)*ALPHA)
                     END IF
                  ELSE IF ( (RL-5+I).GT.RSTEP ) THEN
                     V(I) = RDIP2M(LMP,RSTEP,1,IO)*RADM(RSTEP)
     &                      *DEXP(DBLE(RL-5+I-RSTEP)*ALPHA)
                     B(I) = RDIP2M(LMP,RSTEP,2,IO)*RADM(RSTEP)
     &                      *DEXP(DBLE(RL-5+I-RSTEP)*ALPHA)
                  ELSE
                     V(I) = RDIP2M(LMP,RL-5+I,1,IO)*RADM(RL-5+I)
                     B(I) = RDIP2M(LMP,RL-5+I,2,IO)*RADM(RL-5+I)
                  END IF
               END DO
C
C              calculate derivative of vlmv:
C
               DV1H = (V(6)-V(4))/2.D0
               DV2H = (V(7)-V(3))/4.D0
               DV4H = (V(9)-V(1))/8.D0
               DV = (64.D0*DV1H-20.D0*DV2H+DV4H)
               DV = DV/(45.D0*ALPHA*RADM(RL))
               DVLMV = (DV-RDIP2M(LMP,RL,1,IO))/RADM(RL)
C
C              calculate derivative of vlmb:
C
               DB1H = (B(6)-B(4))/2.D0
               DB2H = (B(7)-B(3))/4.D0
               DB4H = (B(9)-B(1))/8.D0
               DB = (64.D0*DB1H-20.D0*DB2H+DB4H)
               DB = DB/(45.D0*ALPHA*RADM(RL))
               DVLMB = (DB-RDIP2M(LMP,RL,2,IO))/RADM(RL)
C
C              calculate dipole operator:
C
               CALL IND(LMP,L,M)
               RDIP1M(LMP,1,RL,1,IO) = DVLMV + (L+1)/RADM(RL)
     &                                 *RDIP2M(LMP,RL,1,IO)
               RDIP1M(LMP,2,RL,1,IO) = DVLMV - L/RADM(RL)
     &                                 *RDIP2M(LMP,RL,1,IO)
               RDIP1M(LMP,1,RL,2,IO) = DVLMB + (L+1)/RADM(RL)
     &                                 *RDIP2M(LMP,RL,2,IO)
               RDIP1M(LMP,2,RL,2,IO) = DVLMB - L/RADM(RL)
     &                                 *RDIP2M(LMP,RL,2,IO)
            END DO
         END IF
      END DO
C
C     multiply coefficients to v and b:
C
      CV = OMHAR/CLIGHT
      CB = OMHAR/CLIGHT
C
      DO I = 1,NFULLPOT
         DO J = 1,RSTEP
            RDIP2M(I,J,1,IO) = CV*RAD3M(J)*RDIP2M(I,J,1,IO)
            RDIP2M(I,J,2,IO) = CB*RAD3M(J)*RDIP2M(I,J,2,IO)
C
            RDIP1M(I,1,J,1,IO) = RAD3M(J)*RDIP1M(I,1,J,1,IO)
            RDIP1M(I,1,J,2,IO) = RAD3M(J)*RDIP1M(I,1,J,2,IO)
            RDIP1M(I,2,J,1,IO) = RAD3M(J)*RDIP1M(I,2,J,1,IO)
            RDIP1M(I,2,J,2,IO) = RAD3M(J)*RDIP1M(I,2,J,2,IO)
         END DO
      END DO
C
      END
C*==detvbm.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE DETVBM(LAY,ATOM,VLM,RDIP2,IO)
C
C     # purpose      : determine v and b as sum und difference of
C                      vlmup and vlmdown
C
      USE MOD_SPEC
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_TYPES,ONLY:NTMAX
      USE MOD_SPEC_RINDC,ONLY:REL
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATOM,IO,LAY
      COMPLEX*16 RDIP2(NFULLPOT,NRMAX,2,NTMAX),
     &           VLM(LAYSM,NATLM,MQD,NRMAX,NTPHOMAX)
C
C Local variables
C
      INTEGER LMP,LMP1,LMP2,RL
      REAL*8 Y00
C
C*** End of declarations rewritten by SPAG
C
C
      Y00 = 1.D0/DSQRT(4.D0*PI)
C
C
C     the potential vlm has a factor 1/y00
C     multiplied to it, which has to be removed:
C
      DO LMP = 1,NFULLPOT
         LMP1 = REL(LMP,1)
         LMP2 = REL(LMP,2)
         DO RL = 1,RSTEP
            RDIP2(LMP,RL,1,IO) = Y00*(VLM(LAY,ATOM,LMP1,RL,IO)+VLM(LAY,
     &                           ATOM,LMP2,RL,IO))/2.D0
            RDIP2(LMP,RL,2,IO) = Y00*(VLM(LAY,ATOM,LMP1,RL,IO)-VLM(LAY,
     &                           ATOM,LMP2,RL,IO))/2.D0
         END DO
      END DO
C
      END
C*==ind.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE IND(JJ,LLC,MM)
C
      USE MOD_SPEC,ONLY:LL
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER JJ,LLC,MM
C
C*** End of declarations rewritten by SPAG
C
C      LLC = AINT(SQRT(DBLE(JJ-1)))
      LLC = INT(SQRT(DBLE(JJ-1)))
      MM = JJ + LLC - ((LL+1)**2)
C
      END
C*==cgup.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      REAL*8 FUNCTION CGUP(K,MU)
C     /****************************************************************/
C     # purpose      : Clebsch Gordan for spinor up component          *
C     /****************************************************************/
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER K
      REAL*8 MU
C
C Local variables
C
      REAL*8 SIG
C
C*** End of declarations rewritten by SPAG
C
      CGUP = 0.D0
      IF ( K.EQ.0 ) RETURN
      IF ( INT(2.*ABS(MU)).GT.(2*ABS(K)-1) ) RETURN
C
      SIG = -DBLE(K/ABS(K))
      CGUP = SIG*SQRT((DBLE(K)-MU+0.5D0)/DBLE(2*K+1))
C
      END
C*==cgdo.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      REAL*8 FUNCTION CGDO(K,MU)
C     /****************************************************************/
C     # purpose      : Clebsch Gordan for spinor down component        *
C     /****************************************************************/
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER K
      REAL*8 MU
C
C*** End of declarations rewritten by SPAG
C
      CGDO = 0.D0
      IF ( K.EQ.0 ) RETURN
      IF ( INT(2.*ABS(MU)).GT.(2*ABS(K)-1) ) RETURN
C
      CGDO = SQRT((DBLE(K)+MU+0.5D0)/DBLE(2*K+1))
C
      END
C*==cgsp.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      REAL*8 FUNCTION CGSP(K,MU,S)
C     /****************************************************************/
C     # purpose      : Clebsch Gordan for spinors in k,mu notation     *
C                      coupling of l,ml,s,ms -> j,mj                   *
C                      here: half-integer mu                           *
C     # parameter:                                                     *
C                      s = +1 up   - component  (ms=+1/2)              *
C                      s = -1 down - componenet (ms=-1/2)              *
C     /****************************************************************/
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER K,S
      REAL*8 MU
C
C Local variables
C
      REAL*8 SIG
C
C*** End of declarations rewritten by SPAG
C
      CGSP = 0.D0
      IF ( K.EQ.0 ) RETURN
      IF ( INT(2.*ABS(MU)).GT.(2*ABS(K)-1) ) RETURN
C
      SIG = -DBLE(K/ABS(K))
      IF ( S.EQ.1 ) CGSP = SIG*SQRT((DBLE(K)-MU+0.5D0)/DBLE(2*K+1))
      IF ( S.EQ.-1 ) CGSP = SQRT((DBLE(K)+MU+0.5D0)/DBLE(2*K+1))
C
      END
