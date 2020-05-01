C*==spec_convertpot.f    processed by SPAG 6.70Rc at 17:01 on  7 Mar 2017
      SUBROUTINE SPEC_CONVERTPOT()
C
      USE MOD_TRANSPHO_LAYPOT,ONLY:IQ_SPR_LAYAT
      USE MOD_RMESH,ONLY:NRMAX,JRWS,R,DX,FULLPOT,JRCRI
      USE MOD_SITES,ONLY:NQ,NQ_R,NQ_L,ITOQ,NOQ
      USE MOD_TYPES,ONLY:NTMAX,IMT,CONC,Z,BT,VT,LTXT_T,TXT_T
      USE MOD_SPEC
      USE MOD_MPI,ONLY:MPI
      USE MOD_FILES,ONLY:IFILSPECSTR,IPRINT
      USE MOD_LATTICE,ONLY:SYSTEM_DIMENSION
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      REAL*8 ALPHAM(:,:),AX,AY,BX,BY,POS(:,:,:),RBS(:,:),SPA,
     &       VTMP(:,:,:,:,:),ZDIST,ZP(:,:,:)
      INTEGER CIV,IA,IAN,IAT(:,:),ILATT,IM,IMAG,IPOTUD,IQ,IR,IRUMP(:),
     &        ISTACK,IT,ITS,J,JTOP,K,LAYS,NA,NATL(:),NROT,POTKOM
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE IAT,RBS,POS,NATL,IRUMP,ALPHAM
      ALLOCATABLE VTMP, ZP
      ALLOCATE (IAT(NATLM,LAYSM),RBS(LAYSM,NATLM))
      ALLOCATE (POS(3,NATLM,LAYSM),NATL(LAYSM),IRUMP(LAYSM))
      ALLOCATE (ALPHAM(LAYSM,NATLM))
C
      ALLOCATE (VTMP(NRMAX,2,NATLM,LAYSM,NTMAX),ZP(NATLM,LAYSM,NTMAX))
C
C     *****************************************************************
C     We have only one layer (lays = 1) with only one (natl(...) = 1)
C     muffin-tin atom (potkom = 1) that should be at first site
C     (pos(...) = (z,0,0) with z perpendicular to surface)
C     *****************************************************************
      REWIND (IFILSPECSTR)
C
C
C     input data describing layer types in the crystal
C
      READ (IFILSPECSTR,'(4i5)') ILATT,CIV,NROT,ISTACK
C     lets make spag happy
      IF ( .FALSE. ) WRITE (6,*) ILATT,CIV,NROT,ISTACK
C
C
C     lays : no. of different layer types in photoemission
C
      READ (IFILSPECSTR,'(e14.6)') SPA
      IF ( .FALSE. ) WRITE (6,*) SPA
C
C     ax,ay,bx,by = real space lattice vectors c
C
      READ (IFILSPECSTR,'(2e14.6)') AX,AY
      READ (IFILSPECSTR,'(2e14.6)') BX,BY
C     lets make spag happy
      IF ( .FALSE. ) WRITE (6,*) AX,AY,BX,BY
C
C     input data describing layer types in the crystal
C     lays : no. of different layer types
C
      READ (IFILSPECSTR,'(i5)') LAYS
C
C     number of potential components
      POTKOM = 1
C
      IF ( IPRINT.GE.1 ) WRITE (6,99001) LAYS
      IF ( LAYS.GT.LAYSM ) THEN
         WRITE (6,99002) LAYS,LAYSM
         STOP
      END IF
C
      ALLOCATE (IQ_SPR_LAYAT(NATLM,LAYSM))
      DO IT = 1,LAYS
         DO IAN = 1,NATLM
            IQ_SPR_LAYAT(IAN,IT) = 1
         END DO
      END DO
C
      DO IT = 1,LAYS
C
C        for each layer type
C        natl : no. of atoms
C        ia   : types of atoms in locations 1,2..natl
C        pos  : disposition of atoms in locations 1,2..natl
C               along (z,a,b) in units of (spa,!a!,!b!)
C               with respect to layer unit cell origin
C               immediately projected along (z,x,y) in real*8
C               space units
C
         READ (IFILSPECSTR,'(i5)') NATL(IT)
C
         NA = NATL(IT)
         IF ( IPRINT.GE.1 ) WRITE (6,99003) IT,NA
         IF ( NA.GT.NATLM ) THEN
            WRITE (*,99004) IT,NA,NATLM
            STOP
         END IF
C
         DO IA = 1,NA
C
            READ (IFILSPECSTR,*) IAT(IA,IT),IQ_SPR_LAYAT(IA,IT)
            IF ( SYSTEM_DIMENSION(1:2).NE.'3D' ) THEN
               IF ( IQ_SPR_LAYAT(IA,IT).LE.NQ_L .OR. IQ_SPR_LAYAT(IA,IT)
     &              .GE.(NQ-NQ_R+1) ) THEN
                  WRITE (6,*) 
     &              '<SPEC_CONVERTPOT>: In the case of 2D calculations,'
                  WRITE (6,*) 'please avoid to used sites bellow NQ_L'
               END IF
            END IF
C
C            iat(ia,it) = 1
C
            READ (IFILSPECSTR,'(3e14.6)') (POS(J,IA,IT),J=1,3)
C
            IF ( IPRINT.GE.1 ) WRITE (6,99005) IA,IAT(IA,IT),
     &                                (POS(J,IA,IT),J=1,3)
         END DO
C
C        rumpled layer has different components of vector normal to layer
C
         IRUMP(IT) = 0
         ZDIST = POS(1,1,IT)
         DO IA = 2,NA
            IF ( ABS(ZDIST-POS(1,IA,IT)).GT.0.D0 ) IRUMP(IT) = 1
         END DO
         IF ( IPRINT.GE.1 ) THEN
            IF ( IRUMP(IT).EQ.0 ) WRITE (6,99006) IT
            IF ( IRUMP(IT).EQ.1 ) WRITE (6,99007) IT
         END IF
C
      END DO
      WRITE (6,99008)
      WRITE (6,99009)
      WRITE (6,99010)
      DO IT = 1,LAYS
         NA = NATL(IT)
         DO IAN = 1,NA
            IQ = IQ_SPR_LAYAT(IAN,IT)
            WRITE (6,99011) IT,IAN,IQ,
     &                      (ITOQ(K,IQ),TXT_T(ITOQ(K,IQ))(1:LTXT_T
     &                      (ITOQ(K,IQ))),CONC(ITOQ(K,IQ)),K=1,NOQ(IQ))
         END DO
      END DO
C Calculate number of types needed for Photoemission
      NTPHOMAX = 0
      DO IT = 1,LAYS
         NA = NATL(IT)
         DO IAN = 1,NA
            IQ = IQ_SPR_LAYAT(IAN,IT)
            NTPHOMAX = MAX(NTPHOMAX,NOQ(IQ))
         END DO
      END DO
C
      WRITE (6,99008)
C
C
      DO IT = 1,LAYS
         NA = NATL(IT)
         DO IAN = 1,NA
            IQ = IQ_SPR_LAYAT(IAN,IT)
            DO K = 1,NOQ(IQ)
               ITS = ITOQ(K,IQ)
               ZP(IAN,IT,K) = DFLOAT(Z(ITS))
            END DO
         END DO
      END DO
C
C     CONVERT SPRKKR POTENTIAL TO PHOTOEMISSION FORM
C
C
      IPOTUD = 2
      IMAG = 1
C
      DO IT = 1,LAYS
         NA = NATL(IT)
         DO IAN = 1,NA
            IQ = IQ_SPR_LAYAT(IAN,IT)
            DO K = 1,NOQ(IQ)
               ITS = ITOQ(K,IQ)
               IM = IMT(ITS)
               IF ( FULLPOT ) THEN
                  JTOP = JRCRI(IM)
               ELSE
                  JTOP = JRWS(IM)
               END IF
C
               DO IR = 1,JTOP
                  VTMP(IR,1,IAN,IT,K) = (R(IR,IM)*(VT(IR,ITS)+BT(IR,ITS)
     &                                  )/2.0D0)
                  VTMP(IR,2,IAN,IT,K) = (R(IR,IM)*(VT(IR,ITS)-BT(IR,ITS)
     &                                  )/2.0D0)
               END DO
C
               RBS(IT,IAN) = R(JTOP,IM)
               ALPHAM(IT,IAN) = DX(IM)
C
               DO IR = 1,JTOP
                  IF ( VTMP(IR,1,IAN,IT,K)-VTMP(IR,2,IAN,IT,K).LT.0.D0 )
     &                 IMAG = MAX(2,IMAG)
                  IF ( VTMP(IR,1,IAN,IT,K)-VTMP(IR,2,IAN,IT,K).GT.0.D0 )
     &                 IMAG = MAX(3,IMAG)
               END DO
            END DO
         END DO
      END DO
      CALL POTWRITE_RSLAB(NATL,LAYS,ZP,IPOTUD,POTKOM,ALPHAM,VTMP,RBS,
     &                    JTOP)
      IF ( MPI ) CALL DRV_MPI_BARRIER
C
      RETURN
C
99001 FORMAT (///,' NO.OF DIFFERENT LAYER TYPES FOR PHOTOEMISSION',/,
     &        10X,'LAYS',/,10X,i3,/)
99002 FORMAT (' ERROR IN CONVERT: LAYS=',i3,' >  LAYSM=',i3,/)
99003 FORMAT (' LAYER TYPE ',i3,' CONTAINS ',i3,' ATOMS',/)
99004 FORMAT (' LAYER TYPE ',i3,' NATL=',i3,'   NATLM=',i3,/)
99005 FORMAT (' ATOM AT LOCATION ',i3,' IS OF TYPE ',i3,//,
     &        ' DISPLACEMENT FROM UNIT CELL ORIGIN=',3E12.4,/)
99006 FORMAT (' LAYER NUMBER',i3,' IS COPLANAR',/)
99007 FORMAT (' LAYER NUMBER',i3,' IS RUMPLED',/)
99008 FORMAT (//,1X,79('*'))
99009 FORMAT (/,10X,'*************************************************',
     &        /,10X,'**********   TABLE OF SPR-KKR SITES   ***********',
     &        /,10X,'********** used for SPEC calculations ***********',
     &        /,10X,'*************************************************')
99010 FORMAT (5X,'LAYER ',6X,'ATOM IN LAYER ',5X,'SPR-KKR IQ ')
99011 FORMAT (6X,I3,10X,I3,12X,I5,/,(28X,'IT=',I5,3X,A8,3x,'conc=',F8.3)
     &        )
C
      END
C*==potwrite_rslab.f    processed by SPAG 6.70Rc at 17:01 on  7 Mar 2017
      SUBROUTINE POTWRITE_RSLAB(NATL,LAYS,Z,IPOTUD,POTKOM,ALPHA,VTMP,
     &                          RBS,JTOP)
C     /****************************************************************/
C     # purpose      : reads the values of the potential from sprkkr
C                      and writes the in_potm.inp file for rslab
C                      format
C     /****************************************************************/
C
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_SPEC,ONLY:LAYSM,NATLM
      USE MOD_TYPES,ONLY:NTMAX
      USE MOD_SITES,ONLY:NOQ
      USE MOD_TRANSPHO_LAYPOT,ONLY:IQ_SPR_LAYAT
      USE MOD_FILES,ONLY:IPRINT,IFILSPECPOT,FOUND_SECTION
      USE MOD_MPI,ONLY:MPI,MPI_ID
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IPOTUD,JTOP,LAYS,POTKOM
      REAL*8 ALPHA(LAYSM,NATLM),RBS(LAYSM,NATLM),
     &       VTMP(NRMAX,2,NATLM,LAYSM,NTMAX),Z(NATLM,LAYSM,NTMAX)
      INTEGER NATL(LAYSM)
C
C Local variables
C
      REAL*8 EB(3)
      LOGICAL FOUND
      INTEGER I,IAN,IQ,ISPINVAL,IT,K,KC,L,NA,NIN
      CHARACTER*20 INFILE
C
C*** End of declarations rewritten by SPAG
C
      IF ( MPI_ID.EQ.0 ) THEN
         INFILE(1:10) = 'in_pot.inp'
         INQUIRE (FILE=INFILE(1:10),EXIST=FOUND)
         IF ( FOUND ) THEN
            OPEN (IFILSPECPOT,FILE=INFILE(1:10),STATUS='REPLACE')
         ELSE
            OPEN (IFILSPECPOT,FILE=INFILE(1:10),STATUS='NEW')
         END IF
         REWIND (IFILSPECPOT)
C
         DO I = 1,3
            EB(I) = 0.0D0
         END DO
         EB(3) = -1.0D0
C     eb writes the x,y,z component of unity b-field-vector
C
         CALL INPUT_FIND_SECTION('SPEC',0)
C
         IF ( FOUND_SECTION ) CALL SECTION_SET_REAL_ARRAY('BFIELD',EB,
     &        NIN,3,0,9999D0,0)
C
         WRITE (IFILSPECPOT,99001) EB(1),EB(2),EB(3)
         IF ( IPRINT.GE.1 ) WRITE (6,*) 'Direction of B-Field is',EB(1),
     &                                  EB(2),EB(3)
C
C     writes the potential values for a (non)spherical potential
C
         DO IT = 1,LAYS
            WRITE (IFILSPECPOT,99002) IT
C
            NA = NATL(IT)
            DO IAN = 1,NA
               IQ = IQ_SPR_LAYAT(IAN,IT)
C
               DO KC = 1,NOQ(IQ)
                  WRITE (IFILSPECPOT,99004) ALPHA(IT,IAN)
                  WRITE (IFILSPECPOT,99003) IAN,KC
                  WRITE (IFILSPECPOT,99004) RBS(IT,IAN)
C
                  WRITE (IFILSPECPOT,99002) 2*POTKOM
C
                  DO K = 1,POTKOM
                     IF ( IPOTUD.GT.2 ) STOP
C
C                write potential:
C
                     DO ISPINVAL = 1,IPOTUD
                        WRITE (IFILSPECPOT,99005) 1,ISPINVAL
                        WRITE (IFILSPECPOT,99006)
     &                         (VTMP(L,ISPINVAL,IAN,IT,KC)+Z(IAN,IT,KC),
     &                         L=1,JTOP)
                     END DO
C
                  END DO
               END DO
            END DO
         END DO
C
         CALL FLUSH(IFILSPECPOT)
         IF ( MPI ) CLOSE (IFILSPECPOT)
C
      END IF
      IF ( MPI ) CALL DRV_MPI_BARRIER
      IF ( MPI ) THEN
         INFILE(1:10) = 'in_pot.inp'
         OPEN (IFILSPECPOT,FILE=INFILE(1:10),STATUS='OLD')
         REWIND (IFILSPECPOT)
      END IF
C
      RETURN
C
99001 FORMAT (3E14.7)
99002 FORMAT (i4)
99003 FORMAT (2I4)
99004 FORMAT (e14.7)
99005 FORMAT (1x,2I3)
99006 FORMAT (5E14.7)
C
      END
