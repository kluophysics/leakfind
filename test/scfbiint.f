C*==scfbiint.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SCFBIINT(WGT,ZG,ZF,IT,LAM1,LAM2,IKMCB,LMAX,NSOL,IRTOP,
     &                    RINT,RPW,ABITNEW,TL,TG,SL,SG)
C   ********************************************************************
C   *                                                                  *
C   *  auxilary routine to calculate the vector potential due to the   *
C   *                       BREIT INTERACTION                          *
C   *                                                                  *
C   *  26/01/95  HE                                                    *
C   ********************************************************************
      USE MOD_ANGMOM,ONLY:AMEBI1,NLABIMAX,NLMAX,GBIG,GBIL,AMEBI2,NKMMAX
      USE MOD_RMESH,ONLY:R2DRDI,NRMAX
      USE MOD_TYPES,ONLY:IMT,NTMAX,NCPLWFMAX
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C*--SCFBIINT16
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCFBIINT')
      LOGICAL INCLUDE_RETARDATION_TERM
      PARAMETER (INCLUDE_RETARDATION_TERM=.TRUE.)
      REAL*8 PREMAG,PRERET,WOD,WEV
      PARAMETER (PREMAG=(1.0D0/PI)*(4.0D0*PI),PRERET=(1.0D0/PI)
     &           *(4.0D0*PI)**2/3.0D0,WOD=1.0D0/3.0D0,WEV=1.0D0/2.0D0)
C
C Dummy arguments
C
      INTEGER IRTOP,IT,LAM1,LAM2,LMAX,NSOL
      COMPLEX*16 WGT
      REAL*8 ABITNEW(NRMAX,NLABIMAX,-1:+1,NTMAX),RINT(NRMAX,2),
     &       RPW(NRMAX,2*NLMAX),SG(NRMAX,2),SL(NRMAX,2),TG(NRMAX,2),
     &       TL(NRMAX,2)
      INTEGER IKMCB(2)
      COMPLEX*16 ZF(NRMAX,NCPLWFMAX,NKMMAX),ZG(NRMAX,NCPLWFMAX,NKMMAX)
C
C Local variables
C
      REAL*8 AFG,AFGG,AFGL,AGF,AGFG,AGFL,FAC,FACGBI,GG2MGG1,GL2MGL1,
     &       IMWFG,IMWGF,PREL,WR
      INTEGER IKMA,IKMB,ILA,ILR,IM,IR,KA,KB,LA,LR,M,MP,N
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
      IM = IMT(IT)
C
      DO KA = 1,NSOL
         IKMA = IKMCB(KA)
C
         DO KB = 1,NSOL
            IKMB = IKMCB(KB)
C=======================================================================
            DO LR = 1,(2*LMAX+1),2
               ILR = (LR+1)/2
C
C ---------------------------------------------------- set up integrands
               DO IR = 1,IRTOP
                  IMWGF = DIMAG(WGT*ZG(IR,KA,LAM1)*ZF(IR,KB,LAM2))
                  IMWFG = DIMAG(WGT*ZF(IR,KA,LAM1)*ZG(IR,KB,LAM2))
                  WR = R2DRDI(IR,IM)*RPW(IR,LR)
                  TL(IR,1) = IMWGF*WR
                  TL(IR,2) = IMWFG*WR
                  WR = R2DRDI(IR,IM)/RPW(IR,LR+1)
                  TG(IR,1) = IMWGF*WR
                  TG(IR,2) = IMWFG*WR
               END DO
C -------------------------------------------- evaluate radial integrals
               DO N = 1,2
                  SL(1,N) = 0.0D0
                  SG(1,N) = 0.0D0
                  DO IR = 3,IRTOP,2
                     SL(IR,N) = SL(IR-2,N)
     &                          + WOD*(TL(IR-2,N)+4.0D0*TL(IR-1,N)
     &                          +TL(IR,N))
                     SL(IR-1,N) = SL(IR,N) - WEV*(TL(IR-1,N)+TL(IR,N))
C
                     SG(IR,N) = SG(IR-2,N)
     &                          + WOD*(TG(IR-2,N)+4.0D0*TG(IR-1,N)
     &                          +TG(IR,N))
                     SG(IR-1,N) = SG(IR,N) - WEV*(TG(IR-1,N)+TG(IR,N))
                  END DO
               END DO
C ------------------------------------ add prefactor 1/r**(L+1) and r**L
               DO N = 1,2
                  DO IR = 1,IRTOP
                     SL(IR,N) = SL(IR,N)/RPW(IR,LR+1)
                     SG(IR,N) = (SG(IRTOP,N)-SG(IR,N))*RPW(IR,LR)
                     RINT(IR,N) = SL(IR,N) + SG(IR,N)
                  END DO
               END DO
C ------------------------------------------- determine vector potential
               DO M = -1, + 1
C ........................................................ MAGNETIC TERM
                  ILA = ILR
                  LA = LR
                  PREL = PREMAG/DBLE(2*LA+1)
                  AGF = PREL*AMEBI1(IKMA,IKMB,ILA,-M)
                  AFG = PREL*AMEBI2(IKMA,IKMB,ILA,-M)
C     WRITE (6,'(a,4i3,2e15.6)') ' BI ' ,ikma,ikmb,ila,m, AGF,   AFG
                  DO IR = 1,IRTOP
                     ABITNEW(IR,ILA,M,IT) = ABITNEW(IR,ILA,M,IT)
     &                  + RINT(IR,1)*AGF - RINT(IR,2)*AFG
                  END DO
C
C ..................................................... RETARDATION TERM
                  IF ( INCLUDE_RETARDATION_TERM ) THEN
                     DO MP = -1, + 1
                        ILA = ILR
                        FAC = (-1)**(M+MP)*PRERET
                        AGF = FAC*AMEBI1(IKMA,IKMB,ILA,MP)
                        AFG = FAC*AMEBI2(IKMA,IKMB,ILA,MP)
                        GL2MGL1 = GBIL(2,ILA,0,M,MP)
     &                            - GBIL(1,ILA,0,M,MP)
                        GG2MGG1 = GBIG(2,ILA,0,M,MP)
     &                            - GBIG(1,ILA,0,M,MP)
                        AGFL = AGF*GL2MGL1
                        AFGL = AFG*GL2MGL1
                        AGFG = AGF*GG2MGG1
                        AFGG = AFG*GG2MGG1
                        DO IR = 1,IRTOP
                           ABITNEW(IR,ILA,M,IT) = ABITNEW(IR,ILA,M,IT)
     &                        + SL(IR,1)*AGFL - SL(IR,2)*AFGL + SG(IR,1)
     &                        *AGFG - SG(IR,2)*AFGG
                        END DO
C
                        ILA = ILR + 1
                        IF ( ILA.LE.(LMAX+1) ) THEN
                           FACGBI = FAC*GBIL(2,ILR,1,M,MP)
                           AGFL = FACGBI*AMEBI1(IKMA,IKMB,ILR,MP)
                           AFGL = FACGBI*AMEBI2(IKMA,IKMB,ILR,MP)
                           DO IR = 1,IRTOP
                              ABITNEW(IR,ILA,M,IT)
     &                           = ABITNEW(IR,ILA,M,IT) + SL(IR,1)
     &                           *AGFL - SL(IR,2)*AFGL
                           END DO
                        END IF
C
                        ILA = ILR
                        IF ( (ILA-1).GE.1 ) THEN
                           FACGBI = -FAC*GBIL(1,ILA,1,M,MP)
                           AGFL = FACGBI*AMEBI1(IKMA,IKMB,ILA-1,MP)
                           AFGL = FACGBI*AMEBI2(IKMA,IKMB,ILA-1,MP)
                           DO IR = 1,IRTOP
                              ABITNEW(IR,ILA,M,IT)
     &                           = ABITNEW(IR,ILA,M,IT) + SL(IR,1)
     &                           *AGFL - SL(IR,2)*AFGL
                           END DO
                        END IF
C
                        ILA = ILR
                        IF ( (ILA+1).LE.(LMAX+1) ) THEN
                           FACGBI = FAC*GBIG(2,ILA,1,M,MP)
                           AGFG = FACGBI*AMEBI1(IKMA,IKMB,ILA+1,MP)
                           AFGG = FACGBI*AMEBI2(IKMA,IKMB,ILA+1,MP)
                           DO IR = 1,IRTOP
                              ABITNEW(IR,ILA,M,IT)
     &                           = ABITNEW(IR,ILA,M,IT) + SG(IR,1)
     &                           *AGFG - SG(IR,2)*AFGG
                           END DO
                        END IF
C
                        ILA = ILR - 1
                        IF ( ILA.GE.1 ) THEN
                           FACGBI = -FAC*GBIG(1,ILR,1,M,MP)
                           AGFG = FACGBI*AMEBI1(IKMA,IKMB,ILR,MP)
                           AFGG = FACGBI*AMEBI2(IKMA,IKMB,ILR,MP)
                           DO IR = 1,IRTOP
                              ABITNEW(IR,ILA,M,IT)
     &                           = ABITNEW(IR,ILA,M,IT) + SG(IR,1)
     &                           *AGFG - SG(IR,2)*AFGG
                           END DO
                        END IF
                     END DO
                  END IF
C ..................................................... RETARDATION TERM
C
               END DO
C
            END DO
C=======================================================================
         END DO
      END DO
      END
