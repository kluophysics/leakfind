C*==fscat.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FSCAT(IPRINT,TSST,MSST)
C   ********************************************************************
C   *                                                                  *
C   * program to calculate the relativistic general scattering         *
C   * amplitude matrix  f_scat                                         *
C   *                                                                  *
C   ********************************************************************
C   *                                                                  *
C   * 11/10/97  HE                                                     *
C   ********************************************************************
C
      USE MOD_ENERGY,ONLY:ETAB,NEMAX,EFERMI,NETAB
      USE MOD_CALCMODE,ONLY:ORBPOL
      USE MOD_SITES,ONLY:NQ
      USE MOD_ANGMOM,ONLY:NL,NMEMAX,NKMMAX,NLMAX,CGC,KAPTAB,NMUETAB,
     &    LTAB,NKMQ,NLQ,NKM,NK
      USE MOD_FILES,ONLY:DATSET,LSYSTEM,SYSTEM,LDATSET
      USE MOD_TYPES,ONLY:NTMAX,NLT,NT,LTXT_T,TXT_T
      USE MOD_CONSTANTS,ONLY:C0,CI,PI,RY_EV
      IMPLICIT NONE
C*--FSCAT22
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NLEGMAX,NGEOMAX
      PARAMETER (NLEGMAX=50,NGEOMAX=10)
C
C Dummy arguments
C
      INTEGER IPRINT
      COMPLEX*16 MSST(NKMMAX,NKMMAX,NTMAX),TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      REAL*8 AA,BBABSTAB(:,:,:,:,:),BBPHATAB(:,:,:,:,:),COSTET,DIR(3),
     &       DIRFIN(3),DIRINI(3,NGEOMAX),EPLOT(:),PHI,PLEG(:),XJ,XMJ,
     &       XMS,YMAXABST(:),YMAXPHAT(:),YMINABST(:),YMINPHAT(:)
      COMPLEX*16 BB,CYFIN(:,:,:),CYINI(:,:,:),ERYD,FT,
     &           MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),P,
     &           SSST(NKMMAX,NKMMAX,NTMAX),YLMFIN(:,:),YLMINI(:,:)
      CHARACTER*80 FILNAM
      INTEGER I,IA_ERR,IC,IE,IFIL,IGEO,IK,IKM,IKM1,IKM2,IMS,IMS1,IMS2,
     &        IQ,IT,L,LCHPMAX,LFILNAM,LM,MJM05,ML,NCURVES,NE,NGEO,
     &        NLMCHPMAX
      CHARACTER*20 LEG(NLEGMAX)
      CHARACTER*10 STR10
C
C*** End of declarations rewritten by SPAG
C
C
C
C Local variables
C
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE PLEG,CYFIN,CYINI,YLMFIN,YLMINI,EPLOT
      ALLOCATABLE BBABSTAB,BBPHATAB,YMAXABST,YMAXPHAT,YMINABST,YMINPHAT
C
      ALLOCATE (YMAXABST(NT),YMAXPHAT(NT),YMINABST(NT),YMINPHAT(NT))
      ALLOCATE (BBABSTAB(NEMAX,2,2,NGEOMAX,NT))
      ALLOCATE (BBPHATAB(NEMAX,2,2,NGEOMAX,NT),EPLOT(NEMAX))
      ALLOCATE (CYINI(NKMMAX,2,NGEOMAX))
      ALLOCATE (CYFIN(NKMMAX,2,NGEOMAX),PLEG(NLMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:fscat -> YLMINI'
C
C-----------------------------------------------------------------------
C               initialize complex spherical harmonics
C-----------------------------------------------------------------------
      LCHPMAX = NL - 1
      NLMCHPMAX = (LCHPMAX+1)**2
C
      ALLOCATE (YLMFIN(NLMCHPMAX,NGEOMAX))
      ALLOCATE (YLMINI(NLMCHPMAX,NGEOMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:symrot -> Q_CHP'
C
C-----------------------------------------------------------------------
C
      IF ( NGEOMAX.LE.0 ) STOP '<FSCAT> -> NGEOMAX'
      WRITE (6,99002)
C
      DO IT = 1,NT
         NLT(IT) = NL
      END DO
      DO IQ = 1,NQ
         NLQ(IQ) = NL
         NKMQ(IQ) = NKM
      END DO
C
      NGEO = 3
C   GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
      DO IGEO = 1,NGEO
C
         IF ( IGEO.EQ.1 ) THEN
            DIRINI(1,IGEO) = 0D0
            DIRINI(2,IGEO) = 0D0
            DIRINI(3,IGEO) = 1D0
         ELSE IF ( IGEO.EQ.2 ) THEN
            DIRINI(1,IGEO) = 1D0
            DIRINI(2,IGEO) = 0D0
            DIRINI(3,IGEO) = 0D0
         ELSE IF ( IGEO.EQ.3 ) THEN
            DIRINI(1,IGEO) = 1D0/SQRT(3D0)
            DIRINI(2,IGEO) = 1D0/SQRT(3D0)
            DIRINI(3,IGEO) = 1D0/SQRT(3D0)
         END IF
         DO I = 1,3
            DIRFIN(I) = -DIRINI(I,IGEO)
         END DO
C
         DIR(1:3) = DIRINI(1:3,IGEO)
         CALL RVECNORM(3,DIR)
C
         CALL CALC_CHPLM(DIR(1),DIR(2),DIR(3),YLMINI(1,IGEO),LCHPMAX,
     &                   NLMCHPMAX)
C
         DIR(1:3) = DIRFIN(1:3)
         CALL RVECNORM(3,DIR)
C
         CALL CALC_CHPLM(DIR(1),DIR(2),DIR(3),YLMFIN(1,IGEO),LCHPMAX,
     &                   NLMCHPMAX)
C
         IKM = 0
         DO IK = 1,NK
            L = LTAB(IK)
            XJ = DBLE(LTAB(IK)) - 0.5D0*SIGN(1,KAPTAB(IK))
            DO MJM05 = NINT(-XJ-0.5D0),NINT(XJ-0.5D0)
               XMJ = DBLE(MJM05) + 0.5D0
               IKM = IKM + 1
C
               DO IMS = 1,2
                  XMS = -1.5D0 + DBLE(IMS)
                  ML = NINT(XMJ-XMS)
                  LM = L*(L+1) + L + ML + 1
C
                  CYINI(IKM,IMS,IGEO) = CGC(IKM,IMS)
     &                                  *DCONJG(YLMINI(LM,IGEO))
C
                  CYFIN(IKM,IMS,IGEO) = CGC(IKM,IMS)*YLMFIN(LM,IGEO)
C
               END DO
            END DO
         END DO
C
      END DO
C   GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
C
C
      COSTET = -1D0
C
      CALL LEGENDRE(PLEG,NLMAX-1,COSTET)
C
      DO IT = 1,NT
         YMINABST(IT) = 0D0
         YMAXABST(IT) = 0D0
         YMINPHAT(IT) = 0D0
         YMAXPHAT(IT) = 0D0
      END DO
C
      NE = NETAB(1)
C
      DO IE = 1,NE
C
         ERYD = ETAB(IE,1)
         EPLOT(IE) = DREAL(ERYD-EFERMI)*RY_EV
C
         CALL SSITE(0,0,6,.FALSE.,.FALSE.,ERYD,P,IPRINT,NKM,TSST,MSST,
     &              SSST,MEZZ,MEZJ,ORBPOL)
C
C   GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
         DO IGEO = 1,NGEO
C
            DO IT = 1,NT
               IFIL = 30 + IT*IGEO
C
               FT = C0
C
               I = 0
C
               DO IK = 1,NK
C
                  L = LTAB(IK)
                  I = I + NMUETAB(IK)
C
                  FT = FT - TSST(I,I,IT)*DBLE(2*L+1)*PLEG(L+1)/2D0
C
                  IFIL = 20 + IT
C
               END DO
C
               DO IMS1 = 1,2
                  DO IMS2 = 1,2
C
                     BB = C0
C
                     DO IKM1 = 1,NKM
                        DO IKM2 = 1,NKM
C
                           BB = BB + CYINI(IKM1,IMS1,IGEO)
     &                          *CYFIN(IKM2,IMS2,IGEO)
     &                          *TSST(IKM1,IKM2,IT)
C
                        END DO
                     END DO
C
                     BB = -4*PI*BB
                     AA = CDABS(BB)
                     PHI = DREAL(LOG(BB/AA)/CI)
                     IF ( ABS(AA).GT.1D-8 ) THEN
                        BBABSTAB(IE,IMS1,IMS2,IGEO,IT) = AA
                        BBPHATAB(IE,IMS1,IMS2,IGEO,IT) = PHI
                     ELSE
                        BBABSTAB(IE,IMS1,IMS2,IGEO,IT) = 0D0
                        BBPHATAB(IE,IMS1,IMS2,IGEO,IT) = 0D0
                     END IF
                     YMINABST(IT) = MIN(YMINABST(IT),AA)
                     YMAXABST(IT) = MAX(YMAXABST(IT),AA)
                     YMINPHAT(IT) = MIN(YMINPHAT(IT),PHI)
                     YMAXPHAT(IT) = MAX(YMAXPHAT(IT),PHI)
C
                  END DO
               END DO
C
            END DO
         END DO
C   GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
      END DO
C
      IFIL = 7
C
      NCURVES = 4
      LEG(1) = '!DN!DN'
      LEG(2) = '!DN!UP'
      LEG(3) = '!UP!DN'
      LEG(4) = '!UP!UP'
C
      DO IT = 1,NT
         DO IGEO = 1,NGEO
C
            WRITE (STR10,'(A3,I1)') 'GEO',IGEO
C
            CALL XMGRHEAD(DATSET,LDATSET,STR10,4,TXT_T(IT),LTXT_T(IT),
     &                    FILNAM,80,LFILNAM,IFIL,2,EPLOT(1),1,EPLOT(NE),
     &                    1,0D0,0,YMAXABST(IT),1,YMINPHAT(IT),1,
     &                    YMAXPHAT(IT),1,'energy (eV)',11,
     &                    '|B!sm!ss!N!sm!m{1}!ss!M{1}!s''!N(E)|',35,
     &                    '!xF!0!sm!ss!N!sm!m{1}!ss!M{1}!s''!N(E)',37,
     &                    'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM)
     &                    ,25+LSYSTEM,
     &                    'scattering amplitude matrix B of '//TXT_T(IT)
     &                    (1:LTXT_T(IT)),(33+LTXT_T(IT)),.FALSE.)
C
            CALL XMGRCURVES(IFIL,2,NCURVES,NCURVES,2,1,0)
C
            CALL XMGRLEGEND(IFIL,2,NCURVES,NCURVES,LEG,LEG)
C
            IC = -1
            DO IMS1 = 1,2
               DO IMS2 = 1,2
                  IC = IC + 1
                  CALL XMGRTABLE(0,IC,EPLOT,
     &                           BBABSTAB(1,IMS1,IMS2,IGEO,IT),1.0D0,NE,
     &                           IFIL)
                  CALL XMGRTABLE(1,IC,EPLOT,
     &                           BBPHATAB(1,IMS1,IMS2,IGEO,IT),1.0D0,NE,
     &                           IFIL)
               END DO
            END DO
C
            WRITE (6,*) ' '
            WRITE (6,99001) IGEO,(DIRINI(I,IGEO),I=1,3)
            WRITE (6,*) 
     &               '   scattering amplitude matrix B written to file '
     &               ,FILNAM(1:LFILNAM)
            WRITE (6,*) ' '
            CLOSE (IFIL)
C
         END DO
      END DO
C
C ======================================================================
C
      DEALLOCATE (PLEG,CYFIN,CYINI,YLMFIN,YLMINI,EPLOT)
      DEALLOCATE (BBABSTAB,BBPHATAB,YMAXABST,YMAXPHAT,YMINABST,YMINPHAT)
C
      STOP
C
99001 FORMAT ('    for geometry  IGEO =',I2,'  ->p(ini) = (',2(F5.2,',')
     &        ,F5.2,')')
99002 FORMAT (' '//,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*       ******          ****    ****     **    *******       *'
     &  ,/,10X,
     &  '*       *              *    *  *    *   *  *      *          *'
     &  ,/,10X,
     &  '*       *              *       *       *    *     *          *'
     &  ,/,10X,
     &  '*       *****    ***    ****   *       ******     *          *'
     &  ,/,10X,
     &  '*       *                   *  *       *    *     *          *'
     &  ,/,10X,
     &  '*       *              *    *  *    *  *    *     *          *'
     &  ,/,10X,
     &  '*       *               ****    ****   *    *     *          *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
C
      END
C*==legendre.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE LEGENDRE(P,LMAX,X)
C   ********************************************************************
C   *                                                                  *
C   * program to calculate the ordinary Legendre polynoms P_l(x)       *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--LEGENDRE333
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER LMAX
      REAL*8 X
      REAL*8 P(*)
C
C Local variables
C
      INTEGER IL,L
C
C*** End of declarations rewritten by SPAG
C
C
      P(1) = 1D0
      P(2) = X
C
      DO L = 2,LMAX
         IL = L + 1
C
         P(IL) = (DBLE(2*L-1)*X*P(IL-1)-DBLE(L-1)*P(IL-2))/DBLE(L)
C
      END DO
      END
