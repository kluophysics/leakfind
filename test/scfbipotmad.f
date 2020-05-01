C*==scfbipotmad.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SCFBIPOTMAD(ABITNEW)
C
C   ********************************************************************
C   *                                                                  *
C   *     subroutine to calculate the vector potential due to the      *
C   *                   BREIT INTERACTION                              *
C   *                                                                  *
C   *     evaluation of the lattice contribution                       *
C   *                                                                  *
C   *  26/01/95  HE                                                    *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0,C_AU,SQRT_2
      USE MOD_RMESH,ONLY:R,JRWS,NRMAX
      USE MOD_ANGMOM,ONLY:ISMT,IOMT,NLABIMAX,NLMAX,NL
      USE MOD_TYPES,ONLY:IMT,NTMAX,CONC,ITBOT,ITTOP,OBS_T
      USE MOD_SITES,ONLY:NQMAX,NQ,NOQ,ITOQ,IQAT,ABIMAD,IQBOT,IQTOP
      IMPLICIT NONE
C*--SCFBIPOTMAD20
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCFBIPOTMAD')
C
C Dummy arguments
C
      REAL*8 ABITNEW(NRMAX,NLABIMAX,-1:+1,NTMAX)
C
C Local variables
C
      REAL*8 EOVC_AU,MNT_Q(:),RPW(:,:)
      INTEGER I,IA_ERR,IL,ILA,IM,IMLAST,IO,IQ,IR,IT,JQ,JTOP,LA,LMA,LMAX,
     &        M,MA
      CHARACTER*60 OUTPUT
      COMPLEX*16 QSUM
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE MNT_Q,RPW
C
      EOVC_AU = SQRT_2/C_AU
C
C-----------------------------------------------------------------------
C                 deal with lattice contribution
C-----------------------------------------------------------------------
C
      ALLOCATE (MNT_Q(NQMAX))
      ALLOCATE (RPW(NRMAX,2*NLMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: RPW')
C
      MNT_Q(1:NQMAX) = 0.0D0
C
C----------------------------------------------------- moment on site IQ
      DO IQ = IQBOT,IQTOP
         DO IO = 1,NOQ(IQ)
            IT = ITOQ(IO,IQ)
            MNT_Q(IQ) = MNT_Q(IQ) + CONC(IT)
     &                  *(OBS_T(0,ISMT,IT)+OBS_T(0,IOMT,IT))
         END DO
      END DO
C
      IMLAST = 0
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = ITBOT,ITTOP
         IM = IMT(IT)
         JTOP = JRWS(IM)
         IQ = IQAT(1,IT)
         LMAX = NL - 1
C
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
C ------------------------------------------- determine vector potential
         DO LA = 1,(2*LMAX+1),2
            ILA = (LA+1)/2
            DO M = -1, + 1
               MA = M
               LMA = LA*LA + LA + MA + 1
C
               QSUM = C0
               DO JQ = 1,NQ
                  QSUM = QSUM + ABIMAD(IQ,JQ,M,LMA)*MNT_Q(JQ)
               END DO
               QSUM = EOVC_AU*QSUM
C
               DO I = 1,JTOP
                  ABITNEW(I,ILA,M,IT) = ABITNEW(I,ILA,M,IT)
     &                                  + REAL(QSUM*RPW(I,LA))
               END DO
            END DO
         END DO
C
      END DO
C
      DO ILA = 1,3
         OUTPUT = 'ass_mad'//CHAR(ICHAR('1')-1+ILA)
         OPEN (60,FILE=OUTPUT)
         WRITE (6,*) ' writing vector potential to  ',OUTPUT
         DO I = 1,JTOP
            WRITE (60,'(a,2i5,f10.7,3e15.6,3x,e15.6)') ' A ',ILA,I,
     &             R(I,1),(ABITNEW(I,ILA,M,1),M=-1,+1),
     &             (ABITNEW(I,ILA,-1,1)+ABITNEW(I,ILA,+1,1))
         END DO
         CLOSE (60)
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
C
      END
