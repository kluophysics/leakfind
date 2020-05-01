C*==scfbipot.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SCFBIPOT(EFCORRECT,MODE,SHFTEF,ABITNEW,WE0,TAUT,TSST)
C   ********************************************************************
C   *                                                                  *
C   *     subroutine to calculate the vector potential due to the      *
C   *                   BREIT INTERACTION                              *
C   *                                                                  *
C   *     evaluation of the contribution involving the irregular       *
C   *     solution has been suppressed  >>  cirr                       *
C   *                                                                  *
C   *  26/01/95  HE                                                    *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IFILCBWF
      USE MOD_RMESH,ONLY:R,JRWS,NRMAX
      USE MOD_TYPES,ONLY:IMT,NTMAX,ITBOT,ITTOP,NCPLWFMAX,IKMCPLWF
      USE MOD_ANGMOM,ONLY:NL,IKM1LIN,IKM2LIN,NLMAX,NLABIMAX,NKMMAX,
     &    NCPLWF
      IMPLICIT NONE
C*--SCFBIPOT20
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCFBIPOT')
C
C Dummy arguments
C
      LOGICAL EFCORRECT
      INTEGER MODE
      REAL*8 SHFTEF
      COMPLEX*16 WE0
      REAL*8 ABITNEW(NRMAX,NLABIMAX,-1:+1,NTMAX)
      COMPLEX*16 TAUT(NKMMAX,NKMMAX,NTMAX),TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      INTEGER I1,I2,IA_ERR,IKMCB(2),IL,ILA,IM,IMLAST,IR,IRTOP,IT,K1,K2,
     &        KAP1,KAP2,L,LAM1,LAM2,LIN,LMAX,M,MJM05,NSOL
      INTEGER IKAPMUE
      COMPLEX*16 JF(:,:,:),JG(:,:,:),WE,WREG,ZF(:,:,:),ZG(:,:,:)
      REAL*8 MJ,RINT(:,:),RPW(:,:),SG(:,:),SL(:,:),TG(:,:),TL(:,:)
      CHARACTER*60 OUTPUT
C
C*** End of declarations rewritten by SPAG
C
C
C
C Dummy arguments
C
C
C Local variables
C
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE RINT,SG,SL,TG,TL,ZF,ZG,JF,JG,RPW
C
      ALLOCATE (RINT(NRMAX,2),SG(NRMAX,2),SL(NRMAX,2))
      ALLOCATE (TG(NRMAX,2),TL(NRMAX,2),RPW(NRMAX,2*NLMAX))
      ALLOCATE (JF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JG(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZG(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZG')
C
      IF ( EFCORRECT ) THEN
         WE = DCMPLX(SHFTEF,0.0D0)
      ELSE
         WE = WE0
      END IF
C
      IMLAST = 0
C
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = ITBOT,ITTOP
C
         IM = IMT(IT)
         IRTOP = JRWS(IM)
C
         CALL WAVFUN_READ_REL(IFILCBWF,IT,1,ZG,ZF,JG,JF,IRTOP,NCPLWF,
     &                        IKMCPLWF)
C
         IF ( IM.NE.IMLAST ) THEN
            DO IR = 1,NRMAX
               RPW(IR,1) = R(IR,IM)
               DO IL = 2,2*NLMAX
                  RPW(IR,IL) = RPW(IR,IL-1)*R(IR,IM)
               END DO
            END DO
            IMLAST = IM
         END IF
C
         LMAX = NL - 1
         LIN = 0
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
         DO L = 0,LMAX
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            KAP1 = -L - 1
            KAP2 = L
            IF ( L.EQ.0 ) KAP2 = KAP1
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
            DO MJM05 = -L - 1, + L
               MJ = DBLE(MJM05) + 0.5D0
C           no coupling for:  abs(mue)= j   +  j=l+1/2 == kap=-l-1
               IF ( ABS(MJ).GT.DBLE(L) ) THEN
                  NSOL = 1
               ELSE
                  NSOL = 2
               END IF
C------------------------------------------------------------------------
               IKMCB(1) = IKAPMUE(KAP1,NINT(MJ-0.5D0))
               IKMCB(2) = IKAPMUE(KAP2,NINT(MJ-0.5D0))
C------------------------------------------------------------------------
C
               DO K1 = 1,NSOL
                  LAM1 = IKMCB(K1)
                  DO K2 = 1,NSOL
                     LAM2 = IKMCB(K2)
                     LIN = LIN + 1
                     I1 = IKM1LIN(LIN)
                     I2 = IKM2LIN(LIN)
                     IF ( MODE.EQ.1 ) THEN
                        WREG = WE*(TAUT(I1,I2,IT)-TSST(I1,I2,IT))
                     ELSE
                        WREG = WE*TSST(I1,I2,IT)
                     END IF
C
                     CALL SCFBIINT(WREG,ZG,ZF,IT,LAM1,LAM2,IKMCB,LMAX,
     &                             NSOL,IRTOP,RINT,RPW,ABITNEW,TL,TG,SL,
     &                             SG)
C
                  END DO
               END DO
C
            END DO
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
C
         END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
C-----------------------------------------------------------------------
      IF ( .NOT.EFCORRECT ) RETURN
C-----------------------------------------------------------------------
C
      DO ILA = 1,3
         IF ( MODE.EQ.1 ) THEN
            OUTPUT = 'abs'//CHAR(ICHAR('1')-1+ILA)
         ELSE
            OUTPUT = 'ass'//CHAR(ICHAR('1')-1+ILA)
         END IF
         OPEN (60,FILE=OUTPUT)
         WRITE (6,*) ' writing vector potential to  ',OUTPUT
         DO IR = 1,IRTOP
            WRITE (60,'(a,2i5,f10.7,3e20.12,3x,e15.6)') ' A ',ILA,IR,
     &             R(IR,1),(ABITNEW(IR,ILA,M,1),M=-1,+1),
     &             (ABITNEW(IR,ILA,-1,1)+ABITNEW(IR,ILA,+1,1))
         END DO
         CLOSE (60)
      END DO
C
      DEALLOCATE (RINT,SG,SL,TG,TL,ZF,ZG,RPW,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC')
C
      END
