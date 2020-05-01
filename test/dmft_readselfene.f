C*==dmft_readselfene.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_READSELFENE(ERYDC,EFERMI,KSELF,SEVT,SEBT,SEVNST,
     &                            SEBNST,IMT,JRNS1,JRCRI,JRNSMIN,NTMAX,
     &                            NLMFPMAX,NRMAX,FULLPOT)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
      USE MOD_FILES,ONLY:IOTMP
      USE MOD_CONSTANTS,ONLY:C0,PI,RY_EV
      IMPLICIT NONE
C*--DMFT_READSELFENE12
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NSMAX,NREPMAX
      PARAMETER (NSMAX=2,NREPMAX=2)
C
C Dummy arguments
C
      REAL*8 EFERMI
      COMPLEX*16 ERYDC
      LOGICAL FULLPOT
      INTEGER JRNSMIN,NLMFPMAX,NRMAX,NTMAX
      INTEGER IMT(NTMAX),JRCRI(NTMAX),JRNS1(NTMAX),KSELF(NTMAX)
      COMPLEX*16 SEBNST(JRNSMIN:NRMAX,NLMFPMAX,NTMAX),SEBT(NRMAX,NTMAX),
     &           SEVNST(JRNSMIN:NRMAX,NLMFPMAX,NTMAX),SEVT(NRMAX,NTMAX)
C
C Local variables
C
      REAL*8 A,B00I,B00R,B40I,B40R,B44I,B44R,C,D,E,ENE(1000),ERYD,
     &       ESELF(1000),ESELFI(NREPMAX,NSMAX,NTMAX),
     &       ESELFR(NREPMAX,NSMAX,NTMAX),
     &       ESELFTABI(1000,NREPMAX,NSMAX,NTMAX),
     &       ESELFTABR(1000,NREPMAX,NSMAX,NTMAX),ESELFTMP(1000),V00I,
     &       V00R,V40I,V40R,V44I,V44R,VDN00I,VDN00R,VDN40I,VDN40R,
     &       VDN44I,VDN44R,VUP00I,VUP00R,VUP40I,VUP40R,VUP44I,VUP44R,Y00
      INTEGER IDN,IE,IEG,IIE,IM,IN,IR,IREP,IS,ISE,IT,IT2G,IUP,LINP,
     &        NESELF
      CHARACTER*80 INPFILE,MYIT
      LOGICAL READSELF
      REAL*8 YLAG
C
C*** End of declarations rewritten by SPAG
C
C---------- Gaunt coefficients  <2m|l'm'|2m> for real sperical harmonics
C
C 2-2   2-2   0 0      0.2820947917738782
C 2-2   2-2   4 0      0.0402992559676969
C 2-2   2-2   4 4     -0.2384136135044482
C
C 2-1   2-1   0 0      0.2820947917738782
C 2-1   2-1   4 0     -0.1611970238707876
C
C 2 0   0 0   2 0      0.2820947917738782
C 2 0   2 0   4 0      0.2417955358061813
C
C 2 1   2 1   0 0      0.2820947917738782
C 2 1   2 1   4 0     -0.1611970238707876
C
C 2 2   2 2   0 0      0.2820947917738782
C 2 2   2 2   4 0      0.0402992559676969
C 2 2   2 2   4 4      0.2384136135044482
C
      ERYD = DREAL(ERYDC)
      A = 0.2820947917738782D0
      C = -0.2384136135044482D0
      D = -0.1611970238707876D0
      E = 0.2417955358061813D0
C
      Y00 = SQRT(1D0/(4D0*PI))
C
      READSELF = .FALSE.
C
      DO IT = 1,NTMAX
         IM = IMT(IT)
         DO IR = 1,JRCRI(IM)
            SEVT(IR,IT) = C0
            SEBT(IR,IT) = C0
         END DO
         DO IR = JRNS1(IM),JRCRI(IM)
            DO IN = 1,NLMFPMAX
               SEVNST(IR,IN,IT) = C0
               SEBNST(IR,IN,IT) = C0
            END DO
         END DO
      END DO
C------------------------------------------------ loop over atomic types
      DO IT = 1,NTMAX
         IF ( IT.LT.10 ) THEN
            WRITE (MYIT,'(I1)') IT
            LINP = 9
         ELSE
            WRITE (MYIT,'(I2)') IT
            LINP = 10
         END IF
         INPFILE(1:LINP) = 'selfene_'//MYIT
         INQUIRE (FILE=INPFILE(1:LINP),EXIST=READSELF)
         IF ( READSELF ) THEN
            OPEN (UNIT=IOTMP,FILE=INPFILE(1:LINP),STATUS='old')
            KSELF(IT) = 1
            DO IE = 1,20000
               READ (IOTMP,*,END=20,ERR=100) ESELF(IE),
     &               ((ESELFTABR(IE,IREP,IS,IT),ESELFTABI(IE,IREP,IS,IT)
     &               ,IREP=1,NREPMAX),IS=1,NSMAX)
               ESELF(IE) = ESELF(IE)/RY_EV
               DO IREP = 1,2
                  DO IS = 1,2
                     ESELFTABR(IE,IREP,IS,IT) = ESELFTABR(IE,IREP,IS,IT)
     &                  /RY_EV
                     ESELFTABI(IE,IREP,IS,IT) = ESELFTABI(IE,IREP,IS,IT)
     &                  /RY_EV
                  END DO
               END DO
               NESELF = IE
            END DO
 20         CONTINUE
            CLOSE (IOTMP)
C            DO IE=1,NESELF
C                ESELFTABR(IE,1,1,1)=-ESELF(IE)*0.05
C                ESELFTABR(IE,2,1,1)=-ESELF(IE)*0.05
Cc                ESELFTABR(IE,1,1,1)=0.0d0
Cc                ESELFTABR(IE,2,1,1)=0.0d0
C                ESELFTABR(IE,1,2,1)=0.0d0
C                ESELFTABR(IE,2,2,1)=0.0d0
C            END DO
         END IF
      END DO
      DO IT = 1,NTMAX
         DO IS = 1,2
            DO IREP = 1,2
               ESELFR(IREP,IS,IT) = 0.0D0
               ESELFI(IREP,IS,IT) = 0.0D0
            END DO
         END DO
      END DO
C
C--------------------------------------- INTERPOLATION OF THE SELFENERGY
C
      DO IIE = 1,NESELF
         ENE(IIE) = ESELF(IIE) + EFERMI
      END DO
C
      IF ( ERYD.LT.ENE(1) .OR. ERYD.GT.ENE(NESELF) ) THEN
         WRITE (6,*) 'NO SELFENERGY FOR ENERGY',DREAL(ERYDC)
         WRITE (6,*) ERYD,ENE(1),ENE(NESELF)
         DO IR = 1,JRCRI(IM)
            SEVT(IR,IT) = C0
            SEBT(IR,IT) = C0
         END DO
         IF ( FULLPOT ) THEN
            DO IR = JRNS1(IM),JRCRI(IM)
               SEVNST(IR,21,IT) = C0
               SEVNST(IR,25,IT) = C0
               SEBNST(IR,21,IT) = C0
               SEBNST(IR,25,IT) = C0
            END DO
         END IF
         RETURN
      END IF
      DO IT = 1,NTMAX
         IM = IMT(IT)
         IF ( KSELF(IT).EQ.1 ) THEN
            DO IREP = 1,2
               DO IS = 1,2
                  DO ISE = 1,NESELF
                     ESELFTMP(ISE) = ESELFTABR(ISE,IREP,IS,IT)
                  END DO
                  ESELFR(IREP,IS,IT) = YLAG(ERYD,ENE,ESELFTMP,0,3,1000)
                  DO ISE = 1,NESELF
                     ESELFTMP(ISE) = ESELFTABI(ISE,IREP,IS,IT)
                  END DO
                  ESELFI(IREP,IS,IT) = YLAG(ERYD,ENE,ESELFTMP,0,3,1000)
               END DO
            END DO
C
C
            IT2G = 1
            IEG = 2
            IDN = 1
            IUP = 2
            VDN40R = (ESELFR(IT2G,IDN,IT)-ESELFR(IEG,IDN,IT))/(D-E)
            VDN44R = VDN40R*(D-E)/(2*C)
            VDN00R = (ESELFR(IT2G,IDN,IT)-D*VDN40R)/A
            VUP40R = (ESELFR(IT2G,IUP,IT)-ESELFR(IEG,IUP,IT))/(D-E)
            VUP44R = VUP40R*(D-E)/(2*C)
            VUP00R = (ESELFR(IT2G,IUP,IT)-D*VUP40R)/A
C
            V00R = (VDN00R+VUP00R)*0.5D0*Y00
            V40R = (VDN40R+VUP40R)*0.5D0
            V44R = (VDN44R+VUP44R)*0.5D0
C
            B00R = (VDN00R-VUP00R)*0.5D0*Y00
            B40R = (VDN40R-VUP40R)*0.5D0
            B44R = (VDN44R-VUP44R)*0.5D0
C
            VDN40I = (ESELFI(IT2G,IDN,IT)-ESELFI(IEG,IDN,IT))/(D-E)
            VDN44I = VDN40I*(D-E)/(2*C)
            VDN00I = (ESELFI(IT2G,IDN,IT)-D*VDN40I)/A
            VUP40I = (ESELFI(IT2G,IUP,IT)-ESELFI(IEG,IUP,IT))/(D-E)
            VUP44I = VUP40I*(D-E)/(2*C)
            VUP00I = (ESELFI(IT2G,IUP,IT)-D*VUP40I)/A
C
            V00I = (VDN00I+VUP00I)*0.5D0*Y00
            V40I = (VDN40I+VUP40I)*0.5D0
            V44I = (VDN44I+VUP44I)*0.5D0
C
            B00I = (VDN00I-VUP00I)*0.5D0*Y00
            B40I = (VDN40I-VUP40I)*0.5D0
            B44I = (VDN44I-VUP44I)*0.5D0
C
            DO IR = 1,JRCRI(IM)
               SEVT(IR,IT) = DCMPLX(V00R,V00I)
               SEBT(IR,IT) = DCMPLX(B00R,B00I)
            END DO
            IF ( FULLPOT ) THEN
               DO IR = JRNS1(IM),JRCRI(IM)
                  SEVNST(IR,21,IT) = DCMPLX(V40R,V40I)
                  SEVNST(IR,25,IT) = DCMPLX(V44R,V44I)
                  SEBNST(IR,21,IT) = DCMPLX(B40R,B40I)
                  SEBNST(IR,25,IT) = DCMPLX(B44R,B44I)
               END DO
            END IF
         END IF
      END DO
      RETURN
 100  CONTINUE
      STOP 'ERROR READING SELFENE FILE'
      END
