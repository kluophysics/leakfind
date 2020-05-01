C*==sfnmtmesh.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SFNMTMESH(IM,NPAN,MESHN,NM,XRN,DRN,NFUN,FLMSF,
     &                     RMTREDFILL,RMTRED0,NRSFMAX,NPANMAX,NLMSFMAX)
C   ********************************************************************
C   *                                                                  *
C   *   <SFNMTMESH>       adds one extra pannel inside the             *
C   *   muffin-tin sphere to allow lattice relaxations.                *
C   *                                                                  *
C   *   NRADD : number of points added inside the MT radius            *
C   *    = 0 : no extra pannel added                                   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:SQRT_4PI
      IMPLICIT NONE
C*--SFNMTMESH16
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IM,MESHN,NFUN,NLMSFMAX,NPAN,NPANMAX,NRSFMAX
      REAL*8 RMTRED0,RMTREDFILL
      REAL*8 DRN(NRSFMAX),FLMSF(NRSFMAX,NLMSFMAX),XRN(NRSFMAX)
      INTEGER NM(NPANMAX)
C
C Local variables
C
      REAL*8 DIST,STEP
      INTEGER IFUN,IPAN,IR,N1,N2,NRADD
C
C*** End of declarations rewritten by SPAG
C
      IF ( RMTREDFILL-RMTRED0.LT.-1D-8 ) THEN
         WRITE (6,99004) IM,RMTREDFILL,RMTRED0
         STOP 'in <SFNCREATE>'
      END IF
C
C=======================================================================
C                     add a panel if necessary
C=======================================================================
      IF ( ABS(RMTREDFILL-RMTRED0).GT.1D-7 ) THEN
C
         STEP = 0.95D0*(XRN(MESHN)-XRN(1))/DBLE(MESHN-1)
         DIST = XRN(1) - RMTRED0
         NRADD = INT(DIST/STEP) + 1
C
         IF ( DIST.LT.1.0D-5 ) THEN
            WRITE (6,*) 'DIST',DIST
            WRITE (6,*) 'Error from MTMESH '
            WRITE (6,*) 'MT-radius > minimum shape radius '
            WRITE (6,*) 'Your MT-Radius .....',RMTRED0
            WRITE (6,*) 'Shape Radius .......',XRN(1)
            STOP
         END IF
         IF ( NPAN+1.GT.NPANMAX ) THEN
            WRITE (6,FMT=*) ' npan , npanmax ',NPAN + 1,NPANMAX
            STOP
         END IF
         IF ( MESHN+NRADD.GT.NRSFMAX ) THEN
            WRITE (6,99001) MESHN,NRADD,NRSFMAX,STEP,DIST,XRN(1),RMTRED0
            STOP
         END IF
C
C-------------------------------------------- shift data 1 panel upwards
         DO IPAN = NPAN,1, - 1
            NM(IPAN+1) = NM(IPAN)
         END DO
         DO IR = MESHN,1, - 1
            XRN(IR+NRADD) = XRN(IR)
            DRN(IR+NRADD) = DRN(IR)
         END DO
C
         DO IFUN = 1,NFUN
            DO IR = MESHN,1, - 1
               FLMSF(IR+NRADD,IFUN) = FLMSF(IR,IFUN)
            END DO
         END DO
C
         NPAN = NPAN + 1
         MESHN = MESHN + NRADD
         NM(1) = NRADD
C
         STEP = DIST/(NRADD-1)
         DO IR = 1,NRADD
            XRN(IR) = RMTRED0 + STEP*(IR-1)
            DRN(IR) = STEP
         END DO
C
         DO IR = 1,NRADD
            FLMSF(IR,1) = SQRT_4PI
            DO IFUN = 2,NLMSFMAX
               FLMSF(IR,IFUN) = 0.0D0
            END DO
         END DO
C
         WRITE (6,99002)
         N2 = 0
         DO IPAN = 1,NPAN
            N1 = N2 + 1
            N2 = N2 + NM(IPAN)
            WRITE (6,99003) IPAN,XRN(N1),XRN(N2),NM(IPAN)
         END DO
         WRITE (6,*)
C
      END IF
C
C=======================================================================
C
      RETURN
99001 FORMAT (/,'<SFNMTMESH>  MESHN+NRADD > NRSFMAX ',/,10X,
     &        'MESHN     = ',I10,/,10X,'NRADD     = ',I10,/,10X,
     &        'NRSFMAX   = ',I10,/,10X,'STEP      = ',F20.10,/,10X,
     &        'DIST      = ',F20.10,/,10X,'XRN(1)    = ',F20.10,/,10X,
     &        'RMTRED0   = ',F20.10,/)
99002 FORMAT (/,24X,'updated radial mesh',//,15X,'IPAN',6X,'from',10X,
     &        'to',10X,'POINTS',/)
99003 FORMAT (13X,I5,2F14.7,I10)
99004 FORMAT (/,' ##### TROUBLE in <SFNMTMESH> ',80('#'),/,10X,
     &        'for mesh  IM = ',I3,/,10X,'RMT(FILL) =',F10.6,
     &        ' < RMT(INP) =',F10.6,/,80('#'),/)
      END
C*==sfnwrite.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SFNWRITE(IWSF,NPAN,MESHN,NM,XRN,DRN,NFUN,FLMSF,LMISF_S,
     &                    NRSFMAX,NPANMAX,NLMSFMAX)
C   ********************************************************************
C   *                                                                  *
C   *     write the shape function to file IWSF for mesh  IM           *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--SFNWRITE144
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IWSF,MESHN,NFUN,NLMSFMAX,NPAN,NPANMAX,NRSFMAX
      REAL*8 DRN(NRSFMAX),FLMSF(NRSFMAX,NLMSFMAX),XRN(NRSFMAX)
      INTEGER LMISF_S(NLMSFMAX),NM(NPANMAX)
C
C Local variables
C
      INTEGER IFUN,IPAN,IR
C
C*** End of declarations rewritten by SPAG
C
      WRITE (IWSF,FMT=99001) NPAN,MESHN
      WRITE (IWSF,FMT=99001) (NM(IPAN),IPAN=1,NPAN)
      WRITE (IWSF,FMT=99002) (XRN(IR),DRN(IR),IR=1,MESHN)
      WRITE (IWSF,FMT=99001) NFUN
      DO IFUN = 1,NFUN
         WRITE (IWSF,FMT=99001) LMISF_S(IFUN)
         WRITE (IWSF,FMT=99002) (FLMSF(IR,IFUN),IR=1,MESHN)
      END DO
      RETURN
C
99001 FORMAT (16I5)
99002 FORMAT (4D20.12)
      END
