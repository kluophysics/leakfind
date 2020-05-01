C*==rwfbs_norm.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE RWFBS_NORM(IT,IRTOP,GET_ALL_SPHER_SOL,C,ERYD,P,IS,CV0,
     &                      PRL,QRL,PHL,QHL,TMT0L,SMT0L,GAMMA,IGAMMA,
     &                      RGAM,RIGAM,GF_CONV_RH,PJL)
C   ********************************************************************
C   *                                                                  *
C   *      GET_ALL_SPHER_SOL controls which solutions to the           *
C   *      single site problem for a complex spherical potential       *
C   *      functions V(r) == CV0  will be returned                     *
C   *      all solutions up to NLT(IT) will e calculated for type IT   *
C   *                                                                  *
C   *      GET_ALL_SPHER_SOL = FALSE                                   *
C   *                                                                  *
C   *     unnormalized   regular   PH                                  *
C   *                                                                  *
C   *      only up to IRTOP  --  the solution will be continued        *
C   *      by solver for non-spherical  potential                      *
C   *                                                                  *
C   *                                                                  *
C   *      GET_ALL_SPHER_SOL = TRUE                                    *
C   *                                                                  *
C   *       normalized   regular   PR -> J - ip SUM H * t              *
C   *       normalized irregular   PH -> H                             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:R,NRMAX,DRDI,JRCUT,NPAN
      USE MOD_ANGMOM,ONLY:NL,NLMAX
      USE MOD_TYPES,ONLY:NTMAX,IMT,Z,NLT
      USE MOD_CALCMODE,ONLY:IREL
      IMPLICIT NONE
C*--RWFBS_NORM32
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      COMPLEX*16 CI
      PARAMETER (CI=(0.0D0,1.0D0))
      LOGICAL WRONSKI
      PARAMETER (WRONSKI=.FALSE.)
C
C Dummy arguments
C
      REAL*8 C
      COMPLEX*16 ERYD,P
      LOGICAL GET_ALL_SPHER_SOL,GF_CONV_RH
      INTEGER IRTOP,IS,IT
      COMPLEX*16 CV0(NRMAX),PHL(NRMAX,NL),PJL(NRMAX,NL),PRL(NRMAX,NL),
     &           QHL(NRMAX,NL),QRL(NRMAX,NL),SMT0L(NLMAX),TMT0L(NLMAX)
      REAL*8 GAMMA(0:NLMAX,1:NTMAX),IGAMMA(0:NLMAX,1:NTMAX),
     &       RGAM(1:NRMAX,0:NLMAX,1:NTMAX),
     &       RIGAM(1:NRMAX,0:NLMAX,1:NTMAX)
C
C Local variables
C
      COMPLEX*16 ALPHAL,CHL,CHLP1,CJL,CJLP1,CRSQ,DPRIRTOP,H_NORM,NORM,
     &           PLFAC,QJL(:,:),R_NORM,TLP,TSST0,W,WRON,X,Y,ZTOP
      COMPLEX*16 CJLZ,CNLZ
      REAL*8 CSQR,DROVR
      INTEGER IL,IM,IR,L,NSTEP,NVIEW
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE QJL
      ALLOCATE (QJL(NRMAX,NL))
      CSQR = C*C
C
      IM = IMT(IT)
C
      PLFAC = P
C
CLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      DO IL = 1,NLT(IT)
C
         L = IL - 1
C
         PLFAC = PLFAC*DBLE(2*L+1)/P
C
C --------------------------------------- calculate radial wave function
C        unnormalized   regular   (PH,QH)    (temporary)
C          normalized irregular   (PJ,QJ) -> J
C
         CALL RWFBS(GET_ALL_SPHER_SOL,C,ERYD,P,L,CV0,Z(IT),R(1,IM),
     &              DRDI(1,IM),IRTOP,JRCUT(0,IM),NPAN(IM),PHL(1,IL),
     &              QHL(1,IL),DPRIRTOP,PJL(1,IL),QJL(1,IL),NRMAX,
     &              GAMMA(L,IT),IGAMMA(L,IT),IREL,GF_CONV_RH)
C
C=======================================================================
         IF ( GET_ALL_SPHER_SOL ) THEN
C
            ZTOP = P*R(IRTOP,IM)
            DROVR = DRDI(IRTOP,IM)/R(IRTOP,IM)
C
C --------------------------------------------------- calculate t-matrix
C
            CJL = CJLZ(L,ZTOP)
            CHL = CJL + CI*CNLZ(L,ZTOP)
C
            CJLP1 = CJLZ(L+1,ZTOP)
            CHLP1 = CJLP1 + CI*CNLZ(L+1,ZTOP)
C
            X = (DBLE(L)/R(IRTOP,IM))*CHL - P*CHLP1
            Y = (DBLE(L)/R(IRTOP,IM))*CJL - P*CJLP1
C
            W = (DPRIRTOP/(PHL(IRTOP,IL)*DROVR)+GAMMA(L,IT)-1.0D0)
     &          /R(IRTOP,IM)
C
            TLP = -CI*((CJL*W-Y)/(CHL*W-X))
C
            TSST0 = TLP/P
C
            TMT0L(IL) = TSST0
C
            ALPHAL = (CJL-CI*CHL*TLP)*R(IRTOP,IM)
     &               /(PHL(IRTOP,IL)*R(IRTOP,IM)**GAMMA(L,IT))
C
            SMT0L(IL) = ALPHAL*PLFAC/TSST0
C
C --------------------------------------------- normalize wave functions
C          normalized   regular   (PR,QR) -> J - ip SUM H * t
C          normalized irregular   (PH,QH) -> H
C
C NOTE:  the wave functions PHL and QHL contain the prefactor  ip
C
            R_NORM = (CJL-CI*P*CHL*TSST0)
     &               /(PHL(IRTOP,IL)*RGAM(IRTOP,L,IT)/R(IRTOP,IM))
            H_NORM = 1D0/(-CI*P*TSST0)
C
            DO IR = 1,IRTOP
C
               NORM = R_NORM*RGAM(IR,L,IT)
C
               PRL(IR,IL) = PHL(IR,IL)*NORM
               QRL(IR,IL) = QHL(IR,IL)*NORM
C
               PJL(IR,IL) = PJL(IR,IL)*RIGAM(IR,L,IT)
               QJL(IR,IL) = QJL(IR,IL)*RIGAM(IR,L,IT)
C
               PHL(IR,IL) = (PRL(IR,IL)-PJL(IR,IL))*H_NORM
               QHL(IR,IL) = (QRL(IR,IL)-QJL(IR,IL))*H_NORM
            END DO
C
C---------------------------------------------- check wronskian relation
C
            IF ( WRONSKI ) THEN
               CRSQ = (1.0D0+ERYD/CSQR)*P*CI
C
               WRITE (6,99001) IT,L,IS,ERYD
               NSTEP = 20
               NVIEW = 10
               IR = 0
 10            CONTINUE
               IF ( IR.LT.NVIEW .OR. IR.GE.(IRTOP-NVIEW) ) THEN
                  IR = IR + 1
               ELSE IF ( IR.LT.(IRTOP-NVIEW-NSTEP) ) THEN
                  IR = IR + NSTEP
               ELSE
                  IR = IRTOP - NVIEW
               END IF
               IF ( IR.LE.IRTOP ) THEN
                  WRON = QRL(IR,IL)*PHL(IR,IL) - PRL(IR,IL)*QHL(IR,IL)
                  WRITE (6,'(3I4,10F18.14)') IT,L,IR,WRON*CRSQ
                  GOTO 10
               END IF
            END IF
C
C ----------------------------------------------------------------------
C
         END IF
C=======================================================================
      END DO
CLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
      IF ( WRONSKI ) STOP 'in <NRSSITE>:  check of Wronskian done'
99001 FORMAT (/,' wronski relation for IT =',I3,' L =',I3,' IS =',I3,
     &        ' E =',2F10.5)
      END
