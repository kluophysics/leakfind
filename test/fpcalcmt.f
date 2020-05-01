C*==fpcalcmt.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPCALCMT(IMT,Z,RMT,RMTNEW,ALAT,DRDI,A,B,IRWS,R,KSHAPE)
      IMPLICIT NONE
C*--FPCALCMT4
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 A,ALAT,B,RMT,RMTNEW,Z
      INTEGER IMT,IRWS,KSHAPE
      REAL*8 DRDI(*),R(*)
C
C Local variables
C
      REAL*8 DRD1,DRDWS,RIMT,RIMTM1,RNUC,RWS
      INTEGER IDELTA,IH,IMTL,IRWSM2
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C     this subroutine calculates imt and rmt(cal-rmt)
C                     and prints some informations about the used meshes
C        imtl = maximum number of meshpoints generating a radius
C               less or equal than rmt
C        imt  = number of meshpoint generating a new mt-radius closer th
C               mt-radius than every ather meshpoint
C***********************************************************************
C
      IF ( KSHAPE.EQ.0 ) THEN
C th
C        RIMT = LOG(RMT/B+1.D0)/A + 1.D0
         RIMT = LOG(RMT/B+1.D0)/A - 1.0D0
         IMTL = INT(RIMT)
         IRWSM2 = IRWS - 2
         IDELTA = INT((RIMT-IMTL)*2)
         IF ( IDELTA.EQ.0 ) IMT = IMTL
         IF ( IDELTA.GT.0 ) IMT = IMTL + 1
C th
C        RIMTM1 = REAL(IMT-1)
         RIMTM1 = DBLE(IMT+1)
         RMTNEW = B*EXP(A*RIMTM1) - B
C
         IF ( IMT.GT.IRWSM2 ) THEN
            WRITE (6,FMT=99001)
            STOP 'in calcmt'
C
         END IF
C
C        IF (MOD(IMT,2).EQ.0) THEN
C          WRITE (6,FMT=*) ' error stop in calrmt - imt = ',IMT,
C     +      ' has to be odd to get proper core charge ! '
C          STOP 'in calcmt'
C
C        END IF
C
      END IF
C
      IH = IRWS/2
      DRD1 = DRDI(1)
      DRDWS = DRDI(IRWS)
C----- nucleus radius rnuc in bohr's radii
      RNUC = 2.2677022D-5*(2.D0*Z)**(1.0D0/3.0D0)
C-----
      RWS = R(IRWS)
      WRITE (6,FMT=99002) Z,A,B,RNUC,R(2),IH,R(IH),DRD1,DRDWS
      WRITE (6,FMT=99003) IRWS,IMT,RWS,RMT,RMTNEW,ALAT
C
99001 FORMAT (1x,'potentials need more meshpoints',/,50('*'))
99002 FORMAT (' rmesh  z=',f5.2,'  a=',f6.4,'  b=',f8.6,'  rnuc=',f10.8,
     &        '  r(2)=',f10.8,'  r(',i3,')=',f6.4,'   drdi(1)=',f10.8,
     &        '   drdi(irws)=',f8.6)
99003 FORMAT (1x,' irws=',i6,' imt=',i6,' rws=',f12.8,' rmt=',f12.8,
     &        ' rmtnew=',f12.8,' alat=',f12.8)
      END
