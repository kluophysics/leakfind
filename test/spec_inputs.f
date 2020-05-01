C*==openfiles.f    processed by SPAG 6.70Rc at 08:30 on 19 Apr 2017
      SUBROUTINE OPENFILES(INPVER,INFILE,LU,LFIL)
C   ********************************************************************
C   *                                                                  *
C   *  open files for input / output for spectroscopy calculations     *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER INPVER
      CHARACTER*80 INFILE(5)
      INTEGER LFIL(5),LU(5)
C
C Local variables
C
      LOGICAL EXISTS
      INTEGER NOUT1
C
C*** End of declarations rewritten by SPAG
C
      NOUT1 = 6
      IF ( INPVER.EQ.0 ) THEN
C     BAR AND INPUT FILE
         INQUIRE (FILE=INFILE(3)(1:LFIL(3)),EXIST=EXISTS)
         IF ( .NOT.EXISTS ) THEN
            WRITE (6,99001) INFILE(3)(1:LFIL(3))
            WRITE (NOUT1,99001) INFILE(3)(1:LFIL(3))
            STOP
         ELSE
            OPEN (UNIT=LU(3),FILE=INFILE(3)(1:LFIL(3)))
         END IF
         INQUIRE (FILE=INFILE(2)(1:LFIL(2)),EXIST=EXISTS)
         IF ( .NOT.EXISTS ) THEN
            WRITE (6,99001) INFILE(2)(1:LFIL(2))
            WRITE (NOUT1,99001) INFILE(2)(1:LFIL(2))
            STOP
         ELSE
            OPEN (UNIT=LU(2),FILE=INFILE(2)(1:LFIL(2)))
         END IF
      END IF
C
      RETURN
C
99001 FORMAT (2x,'** error file:',2x,a15,1x,'does not exist')
      END
C*==inpstru.f    processed by SPAG 6.70Rc at 08:30 on 19 Apr 2017
      SUBROUTINE INPSTRU(NATL,POS,IRUMP,ISEQ,LAYS,LAYB,IAT,BRUSEP,EMACH,
     &                   AD1,AD2,ILATT,CIV,NROT,ISTACK,BARABU,BARABD)
C   ********************************************************************
C   *                                                                  *
C   *  read in structural file: old format                             *
C   *                                                                  *
C   ********************************************************************
C     /****************************************************************/
C     purpose:            input of structure data                      *
C                                                                      *
C     subroutines called: rotsep                                       *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,PI
      USE MOD_SPEC_GEOM,ONLY:SPA,SEP,AR1,AR2,RAR,TV
      USE MOD_FILES,ONLY:IFILSPECSTR,IFILSPECOU3
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 BARABD,BARABU,EMACH
      INTEGER CIV,ILATT,ISTACK,LAYB,LAYS,NROT
      REAL*8 AD1(2),AD2(2),BRUSEP(3),POS(3,NATLM,LAYSM)
      INTEGER IAT(NATLM,LAYSM),IRUMP(LAYSM),ISEQ(0:LAYSM),NATL(LAYSM)
C
C Local variables
C
      REAL*8 AX,AY,BX,BY,DIST(:,:),DISTX,DISTY,DISTZ,RTV,X,Y,Z
      INTEGER I,IA,IATOM,IJDIST(:,:,:),IP,IS,IT,J,JATOM,N,NA,NDIST,
     &        NFIND,NOUT1
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DIST,IJDIST
      ALLOCATE (DIST(3,NATLM*NATLM*LAYSM),IJDIST(NATLM,NATLM,LAYSM))
C
C      common /output/ip, nout1
      NOUT1 = IFILSPECOU3
      IP = 1
C
C     mark that a new input file is used for batch runs !
      WRITE (NOUT1,*) '** start new input **'
      WRITE (NOUT1,*) '====================='
C
      REWIND IFILSPECSTR
C
C     input of crystalline structure parameters
C     ilatt = 1: bcc
C             2: fcc
C             3: hcp
C             4..etc: unknown at present
C     civ     symmetry of the surface =1: c1v,
C                                     =2: c2v, etc
C             be aware that the hcp(0001) surface atoms belong to c3v and not c6v
C             because of the 3-fold screw axis
C             using c6v for fcc(111) or hcp(0001) will
C             give an average over all possible stacking orders
C     nrot:   not yet defined,
C             to be used for a rotation of the initial orientation of the structure
C     istack: stacking order in fcc, hcp or dhcp lattices (rotates the shift vector by 60?)
C             =1 ab, abc, or abac for hcp, fcc, or dhcp, respectively
C             =2 ac, acb, or acab for hcp, fcc, or dhcp, respectively
C
      READ (IFILSPECSTR,'(4i5)') ILATT,CIV,NROT,ISTACK
      IF ( CIV.LT.1 ) CIV = 1
      IF ( ISTACK.LT.1 .OR. ISTACK.GT.2 ) ISTACK = 1
      WRITE (NOUT1,'(4i5)') ILATT,CIV,NROT,ISTACK
C
C     input for crystal structure parameters
C     spa = lattice constant in bohr
      READ (IFILSPECSTR,'(e14.6)') SPA
C
C     ax,ay,bx,by = real space lattice vectors c
      READ (IFILSPECSTR,'(2e14.6)') AX,AY
      READ (IFILSPECSTR,'(2e14.6)') BX,BY
C
      AX = AX*SPA
      AY = AY*SPA
      BX = BX*SPA
      BY = BY*SPA
      TV = ABS(AX*BY-AY*BX)
      AR1(1) = AX
      AR1(2) = AY
      AR2(1) = BX
      AR2(2) = BY
      AD1(1) = AX
      AD1(2) = AY
      AD2(1) = BX
      AD2(2) = BY
      RTV = 2.D0*PI/TV
      RAR(1) = BY*RTV
      RAR(2) = -BX*RTV
      RAR(3) = -AY*RTV
      RAR(4) = AX*RTV
      IF ( IP.GT.0 ) WRITE (NOUT1,99003) SPA,AR1(1),AR1(2),AR2(1),AR2(2)
     &                      ,(RAR(I),I=1,4)
C
C     input data describing layer types in the crystal
C     lays : no. of different layer types
C
      READ (IFILSPECSTR,'(i5)') LAYS
      IF ( IP.GT.0 ) WRITE (NOUT1,99004) LAYS
      IF ( LAYS.GT.LAYSM ) THEN
         WRITE (NOUT1,99005) LAYS,LAYSM
         STOP
      END IF
      DO IT = 1,LAYS
C
C      for each layer type
C      natl : no. of atoms
C      ia   : types of atoms in locations 1,2..natl
C      pos  : disposition of atoms in locations 1,2..natl
C             along (z,a,b) in units of (spa,!a!,!b!)
C             with respect to layer unit cell origin
C             immediately projected along (z,x,y) in real*8 space units
C
         READ (IFILSPECSTR,'(i5)') NATL(IT)
         NA = NATL(IT)
         IF ( IP.GT.0 ) WRITE (NOUT1,99006) IT,NA
         IF ( NA.GT.NATLM ) THEN
            WRITE (6,99007) I,NA,NATLM
            STOP
         END IF
         DO IA = 1,NA
            READ (IFILSPECSTR,'(i5)') IAT(IA,IT)
            READ (IFILSPECSTR,'(3e14.6)') (POS(J,IA,IT),J=1,3)
            IF ( IP.GT.0 ) WRITE (NOUT1,99008) IA,IAT(IA,IT),
     &                            (POS(J,IA,IT),J=1,3)
            POS(1,IA,IT) = POS(1,IA,IT)*SPA
            X = POS(2,IA,IT)*AX + POS(3,IA,IT)*BX
            Y = POS(2,IA,IT)*AY + POS(3,IA,IT)*BY
            POS(2,IA,IT) = X
            POS(3,IA,IT) = Y
         END DO
C
C       rumpled layer has different components of vector normal to layer
C
         IRUMP(IT) = 0
         Z = POS(1,1,IT)
         DO IA = 2,NA
            IF ( ABS(Z-POS(1,IA,IT)).GT.EMACH ) IRUMP(IT) = 1
         END DO
         IF ( IRUMP(IT).EQ.0 .AND. IP.GT.0 ) WRITE (NOUT1,99009) IT
         IF ( IRUMP(IT).EQ.1 .AND. IP.GT.0 ) WRITE (NOUT1,99010) IT
      END DO
C
C     determine the smallest set of interatom intralayer
C     vector distances for calculation of lattice sums
C
      NDIST = 1
      DIST(1,NDIST) = 0.D0
      DIST(2,NDIST) = 0.D0
      DIST(3,NDIST) = 0.D0
C
C     separation of atom from itself is zero for all atoms in all layers
C
      DO IT = 1,LAYS
         DO IATOM = 1,NATL(IT)
            IJDIST(IATOM,IATOM,IT) = NDIST
         END DO
      END DO
      DO IT = 1,LAYS
         DO I = 1,NATL(IT) - 1
            DO J = I + 1,NATL(IT)
               DISTX = POS(2,J,IT) - POS(2,I,IT)
               DISTY = POS(3,J,IT) - POS(3,I,IT)
               DISTZ = POS(1,J,IT) - POS(1,I,IT)
               NFIND = 0
               DO N = 2,NDIST,2
C                  IF ( DISTX.EQ.DIST(2,N) .AND. DISTY.EQ.DIST(3,N) .AND.
C     &                 DISTZ.EQ.DIST(1,N) ) THEN
                  IF ( ABS(DISTX-DIST(2,N)).LT.1.0D-16 .AND. 
     &                 ABS(DISTY-DIST(3,N)).LT.1.0D-16 .AND. 
     &                 ABS(DISTZ-DIST(1,N)).LT.1.0D-16 ) THEN
                     IJDIST(I,J,IT) = N
                     IJDIST(J,I,IT) = N + 1
                     NFIND = 1
                  END IF
C                  IF ( DISTX.EQ.(-DIST(2,N)) .AND. DISTY.EQ.(-DIST(3,N))
C     &                 .AND. DISTZ.EQ.(-DIST(1,N)) ) THEN
                  IF ( ABS(DISTX+DIST(2,N)).LT.1.0D-16 .AND. 
     &                 ABS(DISTY+DIST(3,N)).LT.1.0D-16 .AND. 
     &                 ABS(DISTZ+DIST(1,N)).LT.1.0D-16 ) THEN
                     IJDIST(I,J,IT) = N + 1
                     IJDIST(J,I,IT) = N
                     NFIND = 1
                  END IF
               END DO
               IF ( NFIND.EQ.0 ) THEN
                  NDIST = NDIST + 1
                  IJDIST(I,J,IT) = NDIST
                  DIST(1,NDIST) = DISTZ
                  DIST(2,NDIST) = DISTX
                  DIST(3,NDIST) = DISTY
                  NDIST = NDIST + 1
                  IJDIST(J,I,IT) = NDIST
                  DIST(1,NDIST) = -DISTZ
                  DIST(2,NDIST) = -DISTX
                  DIST(3,NDIST) = -DISTY
               END IF
            END DO
         END DO
      END DO
      IF ( IP.GT.0 ) WRITE (NOUT1,99015) NDIST
      IF ( IP.GT.0 ) WRITE (NOUT1,99016) ((DIST(J,N),J=1,3),N=1,NDIST)
      DO IT = 1,LAYS
         IF ( IP.GT.0 ) WRITE (NOUT1,99017) IT
         DO IATOM = 1,NATL(IT)
            DO JATOM = 1,NATL(IT)
               IF ( IP.GT.0 ) WRITE (NOUT1,99018) IATOM,JATOM,
     &                               IJDIST(IATOM,JATOM,IT)
            END DO
         END DO
      END DO
C
C        input data describing sequence of layers  starting from
C        the embedding plane and moving into the bulk
C        lays   : no. of layers in overlayer + bulk repeat unit sequence
C        layb   : first layer of bulk repeat unit in sequence
C
      READ (IFILSPECSTR,'(2i5)') LAYS,LAYB
      IF ( IP.GT.0 ) WRITE (NOUT1,99011) LAYS,LAYB
      IF ( LAYS.GT.LAYSM ) THEN
         WRITE (6,99012) LAYS,LAYSM
         STOP
      END IF
C
C        sep(i,0) : displacement of origin of unit cell in first
C                   layer from origin of
C                   barrier plane, along (z,a,b) in units of
C                   (spa,!a!,!b!). immediately projected along
C                   (z,a,b) in real space units
C        iseq(0) = 0 defines the surface barrier type
C
      READ (IFILSPECSTR,'(3e14.6)') (SEP(I,0),I=1,3)
C
C JM: Barier parameters are taken from input or in_bar.inp
      SEP(3,0) = 0.0D0
      SEP(1,0) = BARABU
      SEP(2,0) = BARABD
C
      ISEQ(0) = 0
      DO IS = 1,LAYS
C
C        for each layer of atoms in the sequence
C        iseq : layer type
C        sep  : displacement of origin of unit cell in next layer into
C               the bulk from unit cell origin in current layer,
C               along (z,a,b) in units of (spa,!a!,!b!). immediately
C               projected along (z,x,y) in real space units.
C
         READ (IFILSPECSTR,99001) ISEQ(IS)
         READ (IFILSPECSTR,99002) (SEP(J,IS),J=1,3)
      END DO
C
      SEP(1,0) = SEP(1,0)*SPA
      SEP(2,0) = SEP(2,0)*SPA
C
      DO IS = 1,LAYS
         SEP(1,IS) = SEP(1,IS)*SPA
         X = SEP(2,IS)*AX + SEP(3,IS)*BX
         Y = SEP(2,IS)*AY + SEP(3,IS)*BY
         SEP(2,IS) = X
         SEP(3,IS) = Y
         IF ( IP.GT.0 ) WRITE (NOUT1,99013) IS,ISEQ(IS),
     &                         (SEP(J,IS),J=1,3)
      END DO
C
C     brusep : displacement of origin of unit cell in next bulk repeat
C              unit
C
      BRUSEP(1) = 0.D0
      BRUSEP(2) = 0.D0
      BRUSEP(3) = 0.D0
      DO IS = LAYB,LAYS
         BRUSEP(1) = BRUSEP(1) + SEP(1,IS)
         BRUSEP(2) = BRUSEP(2) + SEP(2,IS)
         BRUSEP(3) = BRUSEP(3) + SEP(3,IS)
      END DO
      WRITE (NOUT1,99014) (BRUSEP(J),J=1,3)
C
      IF ( ISTACK.EQ.2 ) THEN
         CALL ROTSEP(ILATT,CIV,0,ISTACK,LAYS,LAYB,ISEQ,SEP,BRUSEP)
C         reset istack
         ISTACK = 1
      END IF
C
C     close input file
C     close(IFILSPECSTR)
      RETURN
C
99001 FORMAT (20I5)
99002 FORMAT (5E14.6)
99003 FORMAT (1x,'lattice constant      :',1F14.7,/1x,
     &        'real basis 1          :',2F10.5,/1x,
     &        'real basis 2          :',2F10.5,/1x,
     &        'reziprocel basis 1    :',2F10.5,/1x,
     &        'reziprocel basis 2    :',2F10.5)
99004 FORMAT (' no.of different layer types lays',i3,/)
99005 FORMAT (' error in inpstru: lays=',i3,' >  laysm=',i3)
99006 FORMAT (' layer type ',i3,' contains ',i3,' atoms')
99007 FORMAT (' layer type ',i3,' natl=',i3,'   natlm=',i3)
99008 FORMAT (' atom at location ',i3,' is of type ',i3,/,
     &        ' displacement from unit cell origin=',3E12.4)
99009 FORMAT (' layer number',i3,' is coplanar',/)
99010 FORMAT (' layer number',i3,' is rumpled',/)
99011 FORMAT (' no. of layers in overlayer and ',
     &        'bulk repeat unit sequence lays=',i3,/,
     &        ' first layer of bulk repeat unit in sequence layb=',i3)
99012 FORMAT (' lays=',i3,'   laysm=',i3)
99013 FORMAT (' layer no. ',i3,' is of type ',i3,/,
     &        ' displacement to next layer in sequence=',3E12.4)
99014 FORMAT (' bulkrepeat unit= ',3E12.4)
99015 FORMAT (/,' number of different distance vectors connecting',/,
     &        ' atoms in layers: ndist',i3)
99016 FORMAT (' distance vectors(z,x,y):',3E12.4,' au')
99017 FORMAT (' distance vector indices for layer it=',i3)
99018 FORMAT (' atom i=',i3,' to atom j=',i3,' ndist=',i3)
      END
C*==rotsep.f    processed by SPAG 6.70Rc at 08:30 on 19 Apr 2017
      SUBROUTINE ROTSEP(ILATT,CIV,NROT,ISTACK,LAYS,LAYB,ISEQ,SEP,BRUSEP)
C     /****************************************************************/
C     purpose:            rotates the structure                        *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,PI
      USE MOD_SPEC_OUTPUT,ONLY:NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER CIV,ILATT,ISTACK,LAYB,LAYS,NROT
      REAL*8 BRUSEP(3),SEP(3,0:LAYSM)
      INTEGER ISEQ(0:LAYSM)
C
C Local variables
C
      INTEGER IS,J
      REAL*8 X,XPI,Y
C
C*** End of declarations rewritten by SPAG
C
      XPI = 0.D0
      IF ( CIV.NE.0 ) XPI = 2.D0*PI*DBLE(NROT-1)/DBLE(CIV)
      IF ( ISTACK.EQ.2 ) XPI = PI/3.D0
C
C     rotation
      IF ( CIV.GT.0 .OR. ISTACK.EQ.2 ) THEN
         DO IS = 0,LAYS
            X = SEP(2,IS)
            Y = SEP(3,IS)
            SEP(2,IS) = X*COS(XPI) - Y*SIN(XPI)
            SEP(3,IS) = X*SIN(XPI) + Y*COS(XPI)
         END DO
      END IF
      WRITE (NOUT1,*) '======================================='
      WRITE (NOUT1,*) ' lattice             ilatt =',ILATT
      IF ( ISTACK.EQ.2 ) THEN
         WRITE (NOUT1,*) ' change stacking order from ab to ac:'
      ELSE
         WRITE (NOUT1,*) ' layer rotation   irot =',NROT
      END IF
      WRITE (NOUT1,99001) XPI*180.D0/PI
      DO IS = 0,LAYS
         WRITE (NOUT1,99002) IS,ISEQ(IS),(SEP(J,IS),J=1,3)
      END DO
C
C     brusep : displacement of origin of unit cell in next bulk repeat
C              unit
C
      BRUSEP(1) = 0.D0
      BRUSEP(2) = 0.D0
      BRUSEP(3) = 0.D0
      DO IS = LAYB,LAYS
         BRUSEP(1) = BRUSEP(1) + SEP(1,IS)
         BRUSEP(2) = BRUSEP(2) + SEP(2,IS)
         BRUSEP(3) = BRUSEP(3) + SEP(3,IS)
      END DO
      WRITE (NOUT1,99003) (BRUSEP(J),J=1,3)
C
      RETURN
C
99001 FORMAT (2x,'rotation by:',f7.1,' deg')
99002 FORMAT (2x,'layer no. ',i3,' is of type ',i3,/,
     &        ' displacement to next layer in sequence=',3E12.4)
99003 FORMAT (' bulkrepeat unit= ',3E12.4)
      END
C*==potfullm2.f    processed by SPAG 6.70Rc at 08:30 on 19 Apr 2017
      SUBROUTINE POTFULLM2(LAYS,NATL,ZM,IFM,VLM)
C
C     # purpose   : gives the values of the munich potential
C                   full: vlm(laysm,natlm,mlq,rstep,ntype)
C                   para: vrm(rstep,laysm,natlm,ntype)
C                   nrmax = number of mesh points for munich
C
      USE MOD_RMESH,ONLY:NRMAX,JRWS,R,FULLPOT,JRCRI
      USE MOD_SITES,ONLY:ITOQ,NOQ
      USE MOD_TYPES,ONLY:NTMAX,IMT
      USE MOD_TRANSPHO_LAYPOT,ONLY:IQ_SPR_LAYAT
      USE MOD_SPEC
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      USE MOD_SPEC_POTLM,ONLY:MESH,EB,BB,ALPHA
      USE MOD_SPEC_RINDC,ONLY:REL
      USE MOD_FILES,ONLY:IFILSPECPOT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IFM,LAYS
      INTEGER NATL(LAYSM)
      COMPLEX*16 VLM(LAYSM,NATLM,MQD,NRMAX,NTPHOMAX)
      REAL*8 ZM(NATLM,LAYSM,NTMAX)
C
C Local variables
C
      INTEGER ATOM,I,IAN,IM,IQ,IT,ITS,J,JTOP,K,KA,KC,L,LAYPOT,LMC,M,N,
     &        NA,POTKOM,SPIN
      REAL*8 BFIELD,POT(:),RBS(:,:),VIN(:,:),Y00
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE RBS,VIN,POT
      ALLOCATE (RBS(LAYSM,NATLM),VIN(NRMAX,2),POT(NRMAX))
C
      Y00 = 1.D0/DSQRT(4.0D0*PI)
C
      REWIND IFILSPECPOT
C
C     init potentials:
C
      DO I = 1,LAYSM
         DO J = 1,NATLM
            DO K = 1,MQD
               DO M = 1,NRMAX
                  DO N = 1,NTPHOMAX
                     VLM(I,J,K,M,N) = CZERO
                  END DO
               END DO
            END DO
         END DO
      END DO
C
C
C     *** eb reads the x,y,z component of unity b-field-vector ***
C
      READ (IFILSPECPOT,99002) EB(1),EB(2),EB(3)
      IF ( IP.GT.0 ) WRITE (*,99002) EB(1),EB(2),EB(3)
      DO I = 1,3
         BB(I) = EB(I)
      END DO
C
      BFIELD = SQRT(EB(1)**2+EB(2)**2+EB(3)**2)
      IF ( BFIELD.GT.0.0001D0 ) IFM = 1
      IF ( BFIELD.GT.1.0D0 ) THEN
         WRITE (*,99001) BFIELD
         WRITE (NOUT1,99001) BFIELD
      END IF
C
C           This is a small change for fp-spleed
C           Until now, only spherical part of in_pot
C           is read in
      DEALLOCATE (MESH)
      ALLOCATE (MESH(LAYS,NATLM,NRMAX))
      DO IT = 1,LAYS
C
         READ (IFILSPECPOT,99003) LAYPOT
         IF ( IP.GT.0 ) WRITE (6,99003) LAYPOT
C
         DO IAN = 1,NATL(IT)
            IQ = IQ_SPR_LAYAT(IAN,IT)
            ITS = ITOQ(1,IQ)
            IM = IMT(ITS)
C
            IF ( FULLPOT ) THEN
               JTOP = JRCRI(IM)
            ELSE
               JTOP = JRWS(IM)
            END IF
C
C
C
            DO I = 1,JTOP
               MESH(IT,IAN,I) = R(I,IM)
            END DO
C
            DO KA = 1,NOQ(IQ)
               READ (IFILSPECPOT,99005) ALPHA(IT,IAN)
               READ (IFILSPECPOT,99004) ATOM,KC
               READ (IFILSPECPOT,99005) RBS(IT,IAN)
               IF ( IP.GT.0 ) THEN
                  WRITE (6,99005) ALPHA(IT,IAN)
                  WRITE (6,99004) ATOM,KC
                  WRITE (6,99005) RBS(IT,IAN)
               END IF
C
C
C              number of potential components:
C
               READ (IFILSPECPOT,99003) POTKOM
               IF ( POTKOM.GT.2 .AND. NFULLPOT.EQ.1 ) THEN
                  WRITE (*,*) 'you are about to read in nonspherical',
     &                        ' components to the potential.'
                  WRITE (*,*) 'in this case nfullpot should',
     &                        'be set to ',ML**2
                  STOP
               END IF
C
               DO I = 1,NRMAX
                  VIN(I,1) = 0.D0
                  VIN(I,2) = 0.D0
               END DO
C
               DO K = 1,POTKOM
C
C                 nonrelativistic index of potential component, spin (1=down):
C
                  READ (IFILSPECPOT,99006) LMC,SPIN
C
                  IF ( SPIN.GT.2 ) STOP
C
                  READ (IFILSPECPOT,99007) (POT(L),L=1,JTOP)
C
                  IF ( LMC.EQ.1 ) THEN
C                    spherical components:
C
                     DO L = 1,JTOP
                        VIN(L,SPIN) = POT(L)
                        VLM(IT,ATOM,REL(LMC,SPIN),L,KA)
     &                     = (DCMPLX(VIN(L,SPIN),0.D0)-ZM(IAN,IT,KA))
     &                     /(MESH(IT,IAN,L)*Y00)
                     END DO
                  ELSE
C
C                    nonspherical components:
C
                     DO L = 1,JTOP
                        VLM(IT,ATOM,REL(LMC,SPIN),L,KA)
     &                     = DCMPLX(POT(L),0.D0)/Y00
                     END DO
                  END IF
               END DO
C
            END DO
         END DO
      END DO
C
      IF ( IP.GT.0 ) THEN
C
C        print magnetization
C
         WRITE (NOUT1,99008) (EB(I),I=1,3)
C
C        print fullpot vlm
C
         DO IT = 1,LAYS
            NA = NATL(IT)
            DO IAN = 1,NA
               IQ = IQ_SPR_LAYAT(IAN,IT)
               DO KA = 1,NOQ(IQ)
                  WRITE (NOUT1,99009) IT,IAN,KA
                  DO K = 1,MQD
                     IF ( CDABS(VLM(IT,IAN,K,1,KA)).GT.EPS12 ) THEN
                        WRITE (NOUT1,99010) K,MESH(IT,IAN,1),
     &                         VLM(IT,IAN,K,1,KA)*MESH(IT,IAN,1)
     &                         *Y00 + ZM(IAN,IT,KA)
                        WRITE (NOUT1,99010) K,MESH(IT,IAN,JTOP),
     &                         VLM(IT,IAN,K,JTOP,KA)*MESH(IT,IAN,JTOP)
     &                         *Y00 + ZM(IAN,IT,KA)
                     END IF
                  END DO
               END DO
            END DO
         END DO
C
      END IF
C
      RETURN
C
99001 FORMAT ('* attention, super-ferromagnetic case b>1 =',f14.7)
99002 FORMAT (3E14.7)
99003 FORMAT (i4)
99004 FORMAT (2I4)
99005 FORMAT (e14.7)
99006 FORMAT (1x,2I3)
99007 FORMAT (5E14.7)
99008 FORMAT ('magnetization eb=',3(1x,e14.7))
99009 FORMAT ('non zero full potential vlm for layer',i3,'atom: ',i3,2x,
     &        'type',2x,i3)
99010 FORMAT (1x,i3,3(1x,e14.7))
      END
C*==input.f    processed by SPAG 6.70Rc at 08:30 on 19 Apr 2017
      SUBROUTINE INPUT(EMESH,ESTEP,EFEV,WFEV,VIH,VIL,OMEV,EREAL,VPREV,
     &                 THQ,FIQ,ALQ,NOUT1,IP,NREL,NPOL,NSPLEE,NSPIN,SPOL,
     &                 LAYER1,GANZ,LANZ1,LANZ2,TMIN,TMAX,PMIN,PMAX,NT,
     &                 NP,TYP,ISTR,TEMP,DTEMP,MASS,POL0,POL0L,IDORA,
     &                 IDREH,ICIRC,IFSP,Q1,Q2,Q3,Q4,IBLOCH,CORE,ESTART,
     &                 EEND,ASYM,XMATOFF,PKSCAN,PKEMIN,PKEMAX,ASIG,
     &                 FROMSIGMA,TOSIGMA,EMIN,EMAX,VBH,VBL,VBSTP,NATL,
     &                 LAYS,DELQ,MCD,USEEULER,NPARA,IXAS,LOWHM,LOIHM,
     &                 LOEF,IUSEDOS,IDCALC,IGCONV,FWHM,FTEMP)
C
C   ********************************************************************
C   *                                                                  *
C   *  read in input file: old format                                  *
C   *                                                                  *
C   ********************************************************************
C     /****************************************************************/
C     purpose       : read program control parameters
C     /****************************************************************/
C
      USE MOD_TYPES,ONLY:
      USE MOD_TRANSPHO_LAYPOT,ONLY:IQ_SPR_LAYAT
      USE MOD_SITES,ONLY:NOQ
      USE MOD_ENERGY,ONLY:EFERMI
      USE MOD_CONSTANTS,ONLY:RY_EV
      USE MOD_SPEC,ONLY:LAYSM,MPW
      USE MOD_FILES,ONLY:IFILSPECINP
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALQ,ASYM,DELQ,DTEMP,EEND,EFEV,EMAX,EMESH,EMIN,EREAL,ESTART,
     &       ESTEP,FIQ,FROMSIGMA,FTEMP,FWHM,LOEF,LOIHM,LOWHM,MASS,OMEV,
     &       PKEMAX,PKEMIN,PMAX,PMIN,TEMP,THQ,TMAX,TMIN,TOSIGMA,VBH,VBL,
     &       VBSTP,VIH,VIL,VPREV,WFEV
      INTEGER ASIG,CORE,GANZ,IBLOCH,ICIRC,IDCALC,IDORA,IDREH,IFSP,
     &        IGCONV,IP,IUSEDOS,IXAS,LANZ1,LANZ2,LAYER1,LAYS,MCD,NOUT1,
     &        NP,NPARA,NPOL,NREL,NSPIN,NSPLEE,NT,PKSCAN,SPOL,TYP,
     &        USEEULER,XMATOFF
      COMPLEX*16 Q1,Q2,Q3,Q4
      INTEGER ISTR(2),NATL(LAYSM)
      REAL*8 POL0(3),POL0L(3)
C
C Local variables
C
      INTEGER IAN,IO,IQ,IT,ITMP,LINESREAD
      REAL*8 RTMP
C
C*** End of declarations rewritten by SPAG
C
C     ALLOCATABLE z
C     ALLOCATE(z(natlm,laysm,NTMAX))
      REWIND (IFILSPECINP)
C
C     emesh = number of points in energy mesh
C     estep = next point in energy mesh (stepwidth)
C     efev  = fermi energy in ev
C             startenergy for ups and aes calculations
C     vih   = imaginary part of the potential in ev (final state)
C     vil   = imaginary part of the potential in ev (initial state)
C     omev  = photon energy in ev
C     ereal = controls absorption in the mtp (ereal = 1)
C     vprev = inner potential of the bulk crystal in ev
C     thq,fiq =   polar and azimuthal angles of photon beam
C                 normal incidence thq=180 !!
C     alq   = alignment of polarization vector or pol.ellipsis
C     delq  = phase shift between real and imaginary part of e-vector,
C             delq=90 for circular polarized light !!
C
      LINESREAD = 1
      READ (IFILSPECINP,99001,ERR=100) EMESH,ESTEP,EFEV,VIH,VIL
      LINESREAD = LINESREAD + 1
      READ (IFILSPECINP,99001,ERR=100) OMEV,EREAL,VPREV,THQ,FIQ
      EREAL = EFERMI*RY_EV
      VPREV = VPREV + EREAL
      LINESREAD = LINESREAD + 1
      READ (IFILSPECINP,99002,ERR=100) ALQ,DELQ
      LINESREAD = LINESREAD + 1
C
      WFEV = VPREV - EFEV
C     z(natlm,laysm,ntype)  = atomic numbers of the atoms
C     rn = muffin-tin-radius for the corresponding atom
C
      DO IT = 1,LAYS
         DO IAN = 1,NATL(IT)
            IQ = IQ_SPR_LAYAT(IAN,IT)
            DO IO = 1,NOQ(IQ)
               READ (IFILSPECINP,99002,ERR=100) RTMP
            END DO
         END DO
      END DO
C
      DO IT = 1,LAYS
         DO IAN = 1,NATL(IT)
            READ (IFILSPECINP,99002,ERR=100) RTMP
         END DO
      END DO
C     This line only to make spag happy
      IF ( .FALSE. ) WRITE (6,*) RTMP
C
C
C     nout1   = standard outputkanal (usually 4)
C     ip      = controls amount of output
C             = -1 minimal
C             = 0  few
C             = 1  more essential
C             = 2..4  even more, most, all
C             = 5  more than enaugh
C     nrel    = 0 nonrelativistic calculation,
C             = 1 relativistic calculation
C     idora   controls calculation in upsrun
C             in upsrun (ibloch=4 etc)
C             = 0 paramagnetc calc. (use solvedirac),
C             = 1 ferromagnetic calc. (use doublerad)
C             in bands (ibloch = 1)
C             = 0 use rslab,
C             = 1 use fullf
C     nsplee  = 0 spleed calculation,
C             = 1 photoemission calculation
C     nspin   = 2: has to be 2 (but is principally not longer in use)
C     spol    = 1 unpolarized calculation,
C             = 2 spin-polarized calculation
C     layer1  = number of bulk layer in a photoemission calculation
C     ganz    = number of reciprocal lattice vectors
C     lanz1   = | number of layerdoublings for final state
C     lanz2   = | number of layerdoublings for initial state
C     ibloch  = 0  band structure calc.
C                  calculate core-level if core = 1
C             = 1  spleed calculation,
C             = 2  xps (corelevel) calc.
C             = 3  auger calc.
C             = 4  ups (valence-band) calc. (arups)
C             = 5  xes calc.
C             = 6  xas calc.
C             = 7  angular integrated ups
C             = 8  secondary emission spectrum
C             = 9  2ppe => two photon photoemission (work in progress)
C            >= 10 special purpose (see rslab)
C
      READ (IFILSPECINP,99003,ERR=100) ITMP,IP,NREL,IDORA,NSPLEE,NSPIN,
     &                                 SPOL,LAYER1,GANZ,LANZ1,LANZ2,
     &                                 IBLOCH
      LINESREAD = LINESREAD + 1
C     This line only to make spag happy
      IF ( .FALSE. ) WRITE (6,*) ITMP
C
C     tmin,tmax     define the range of polar angles
C     pmin,pmax     define the range of azimut angles
C     temp        = temperature parameter
C     dtemp       = debeye temperature
C     mass        = atomic mass in units of the electron mass
C                   matom/me= 1.822843e+3
C
      READ (IFILSPECINP,99001,ERR=100) TMIN,TMAX,PMIN,PMAX,TEMP
      LINESREAD = LINESREAD + 1
      READ (IFILSPECINP,99001,ERR=100) DTEMP,MASS
      LINESREAD = LINESREAD + 1
C
C     pol0,pol0l = initial pol. ,initial pol. in the laboratory system
C
      READ (IFILSPECINP,99001,ERR=100) POL0(1),POL0(2),POL0(3),POL0L(1),
     &                                 POL0L(2)
      LINESREAD = LINESREAD + 1
      READ (IFILSPECINP,99001,ERR=100) POL0L(3)
      LINESREAD = LINESREAD + 1
C
C     nt,np  numbers of angular values for a rotation diagram
C            nt: polar,  np: azimuth
C     typ    crystal coordinats in splout, xpsrun, or upsrun
C           =0: i(e) diagram,
C           =1: rotation diagram,         -> phi scan
C           =2: scattering-angle diagram  -> theta scan
C           =3: orthonormal projection    | 3,4 only for angular resolved
C           =4: stereographic projection  | pe (ups, xps) note: nt=np-> nx,ny
C     istr   beam number (h,k)
C
      READ (IFILSPECINP,99004,ERR=100) NT,NP,TYP,ISTR(1),ISTR(2)
      LINESREAD = LINESREAD + 1
C
C     dimension for angular scan is actually too small (mpw=1)
C     but other fields have to be cleaned up before larger dimensions will work
      IF ( NT.GT.MPW .OR. NP.GT.MPW ) THEN
         WRITE (*,*) '*** dimension for angular scan too small ***'
         WRITE (*,*) 'nt, np =',NT,NP
         WRITE (*,*) '***     change dimensions in parms.h     ***'
         STOP
      END IF
C
C     # parameter controlling the photon polarisation and dichroism
C     npol    controls the polarization and dichroism
C             = 0 unpolarized and p-s dichroism for the calculation
C             = 1 p-pol or rcp or elliptical (depends on icirc, etc)
C             = 2 s-pol or lcp or elliptical (depends on icirc, etc)
C             = 3 dichroism (cdad, ldad)
C     icirc   controls polarisation / ellipticity of the photons
C             = 0 elliptically pol. light:  alq, delq arbitrary
C             = 1 linear pol. light:        alq arbitrary, delq=0
C             = 2 circular pol. light:      alq=45, delq = 90
C                 (should not depend on alq, but the matrixelements become symmetric)
C             = 3 (icirc=2 with fresnel optics )
C             = 4 fresnel optics  (for all kind of polarizations like icirc=0)
C     idreh   controls helicity of the photons
C             = 1, -1: sigma+, sigma- => right or left circular polarization
C               (note: lcp, rcp are exchanged in some books)
C             = 0 linearly polarized (equals icirc=1)
C     ifsp    = 0 fixed photon azimuth angle
C             = 1 variable (equal to electron) azimuth angle
C
      READ (IFILSPECINP,99005,ERR=100) NPOL,ICIRC,IDREH,IFSP
      LINESREAD = LINESREAD + 1
C
C     q1,q2,q3,q4 = amplitudes of the photoelectron used in spin
C                   polarized calculations
C
      READ (IFILSPECINP,99001,ERR=100) Q1,Q2,Q3,Q4
      LINESREAD = LINESREAD + 1
C
C     # parameters for all core level calculations (xps, xas, xes, aes)
C     core    = 1 calculate core wavefunctions (xps, aes, etc.)
C                 works also with ibloch = 1
C                 calculate corelevel instead of bands
C             = 2 read core wavefunctions from file
C                 has to be 1 for mcd calculations
C     estart  = start energy for xps, aes, etc in ev
C     eend    = stop energy for xps, aes, etc in ev
C               estart, estop define the energy intervall where
C               core level are searched
C
      READ (IFILSPECINP,99006,ERR=100) CORE,ESTART,EEND
      LINESREAD = LINESREAD + 1
C
C     # parameters for spleed, ups, xps, xas, xes
C     mcd         = 0 no mcd
C                 = 1 calculate regular mcd for ups, xps, xas, or xes
C                     (invert b: b -> -b)
C                 = 2 in plane 90deg mcd (example: bx -> by)
C     xmatoff     = 0 use xmatrix
C                 = 1 turn off xmatrix
C     useeuler    = 0 calculate wave functions using full bfield
C                 = 1 rotate bfield and calculate for bz
C                 (note: useeuler will be changed by selecteuler,
C                        if inadequately used)
C     npara       = 0 switch off selection rules
C                 = 1 switch on selection rules
C
      READ (IFILSPECINP,99007,ERR=100) MCD,XMATOFF,USEEULER,NPARA
      LINESREAD = LINESREAD + 1
C
C     asym    = asymmetry parameter for doniach-sunjic lineshape
      READ (IFILSPECINP,99012,ERR=100) ASYM
      LINESREAD = LINESREAD + 1
C
C     # parameters for energy scans
C     pkscan  use peakscan for calculation of energy (result depends on ibloch)
C                 (pkscan is not in use for bands or spleed, ibloch=0,1)
C             = 0 in xps, xes, xas, and direct:
C                    calculate with fixed energymesh in [pkemin,pkemax];
C                 in band calculate with fixed energymesh from 0
C                 in ups calculate with fixed energymesh from efev
C                 in aes calculate with fixed energymesh from efev
C             = 1 in xps calculate mesh of energypoints in the peakregion
C                    finer than in the valley region;
C                 in xes, xas use energyintervall determined by
C                    core-level binding energies and dos
C                 in band, ups calculate with fixed energymesh in [pkmin,pkmax]
C                 in directrun: use energymesh from dos
C                    mainly used if called from secondaryrun;
C             = 2 in xps calculate intensities only for the exact eigenvalues within
C                    the energy-interval [pkemin,pkemax];
C                 in band calculate with fixed stepwidth in [pkemin,pkemax];
C                 in xes, xas use omhar read from inputfile
C     pkemin      |
C     pkemax      | [pkmin,pkmax] energy interval for calculation of energy-vector
C
      READ (IFILSPECINP,99009,ERR=100) PKSCAN,PKEMIN,PKEMAX
      LINESREAD = LINESREAD + 1
C
C     # parameters for all core level calculations (xps, xas, xes, aes)
C     asig        adjust self energy from fromsigma (emin) to tosigma (emax),linearly
C                 = 0 constant, use sigma from vil
C                 = 1 adjust sigma linearly
C                 = 2 adjust sigma = f(j,mj)
C                 asig=2 is new, see subroutine detg
C                 asig=1 should be corrected if out of range, depending on pkscan
C     emin        |
C     emax        | energy intervalls for self energy interpolation
C     fromsigma   |
C     tosigma     | first and last value  for self energy interpolation
      READ (IFILSPECINP,99010,ERR=100) ASIG,EMIN,EMAX,FROMSIGMA,TOSIGMA
      LINESREAD = LINESREAD + 1
C
C     # parameters for aes
C     vbstp   = stepwidth for valence band in aes
C     vbh     = width of valence band for aes
C     vbl     = limit for valence band in aes
      READ (IFILSPECINP,99011,ERR=100) VBSTP,VBH,VBL
      LINESREAD = LINESREAD + 1
C
C     # parameters for xas, xes, aes
C     ixas    = method for xas
C                 =0 incoherent, =1 coherent
C     lowhm   = half width at half maximum of greens function
C     loihm   = width about ef
C     loef    = width and height of greens-f. at ef
C               the width of the accompanied lorentzian varies
C               between loef(e=ef=0) and lowhm*pi/2 (|e|>>ef)
C               loihmc determines how fast the transition
C               between these values occurs.
C               lowhm = 0.25d0:  seems to work fine with xas
C               loihm = 1.0d0 :  should be approximately 1.+-0.1 for xas
C               otherwise structures of dos are shifted too much !!!
C               loef = 0.1d0  :  seems to work fine with xas
      READ (IFILSPECINP,99013,ERR=100) IXAS,LOWHM,LOIHM,LOEF
      LINESREAD = LINESREAD + 1
C
C     # parameters for secondaryrun
C     iusedos    = 0  use direct emission spectrum
C                = 1  use dos only (whithout matrix elements)
C     idcalc     = 0  do not run direct emission
C                     read direct emission spectrum from file (r41, l51)
C                     or use dos only (iusedos=1)
C                = 1  calculate direct emission spectrum
      READ (IFILSPECINP,99008,ERR=100) IUSEDOS,IDCALC
      LINESREAD = LINESREAD + 1
C
C     # parameters for dirac fermi and gauss convolution
C     igconv  parameter for broadening of spectra in experiment
C             =0: off,
C             =1: convolute with gaussian (all spectra)
C             =2: multiply by occupied fermi (ups)
C             =3: multiply by unoccupied fermi (ipes)
C             =4: multiply by fermi (2) and convolute with gaussian (ups)
C             =5: multiply by fermi (3) and convolute with gaussian (ipes)
C     fwhm    = full width at half maximum for gauss convolution
C               is also used in cgauss (called by det_dos)
C               <1e-3 switch off gauss broadening in cgauss (det_dos)
C     ftemp   temperature in fermi-dirac distribution
C             ftemp is also used in det_dos
      READ (IFILSPECINP,99013,ERR=100) IGCONV,FWHM,FTEMP
      LINESREAD = LINESREAD + 1
C
C     # parameters for xas, xes, aes
C     dosunit = logical unit numbers for files containing the dos
C               for each atom and layer
      DO IT = 1,LAYS
C        read (IFILSPECINP, 990,err=300) (dosunit(ian,it),ian=1,natl(it))
         LINESREAD = LINESREAD + 1
      END DO
C
C     check and correct for wrong inputs
C     ==================================
C
C     correct input for use of dos if idcalc=1
      IF ( IDCALC.EQ.1 ) IUSEDOS = 0
C
C     correct input for non-scattering methods
      IF ( IBLOCH.EQ.5 .OR. IBLOCH.EQ.6 .OR. IBLOCH.EQ.7 .OR. 
     &     IBLOCH.EQ.8 ) THEN
         LAYER1 = 1
         GANZ = 1
         LANZ1 = 1
         LANZ2 = 1
      END IF
      IF ( IBLOCH.EQ.8 .AND. MCD.NE.0 .AND. (NPOL.EQ.0 .OR. NPOL.EQ.3) )
     &     THEN
         WRITE (*,*) 
     &     ' mcd only for fixed polarisation possible in secondaryrun !'
         WRITE (*,*) ' program stopped, check method in input   !'
         WRITE (NOUT1,*) 
     &     ' mcd only for fixed polarisation possible in secondaryrun !'
         WRITE (NOUT1,*) ' program stopped, check method in input   !'
         STOP
      END IF
C
      IF ( NSPIN.NE.2 ) THEN
         WRITE (*,*) ' nspin is reset to 2 !'
         WRITE (NOUT1,*) ' nspin is reset to 2 !'
         NSPIN = 2
      END IF
C
C     correct the photon polarisation such that npol, icirc, etc. fit together
      IF ( NPOL.EQ.0 ) ICIRC = 1
      SELECT CASE (ICIRC)
      CASE (1,4)
C         note that s and p are not defined a priori !
C         for thq=0,180 they cannot be distinguished, but s=ey, p=ex
C         for thq/=0 the following gives  s, p for phiq=0, 180,
C                                 and     p, s for phiq=90, 270,
C                                 the notation s,p will depend also on alq for other phq
C         however, the construction of unpolarised light is always correct !
C
         DELQ = 0.D0
         IDREH = 0
         IF ( NPOL.EQ.0 .OR. NPOL.EQ.1 ) THEN
            ALQ = 0.D0
         ELSE IF ( NPOL.EQ.2 ) THEN
            ALQ = 90.D0
         END IF
      CASE (2,3)
         ALQ = 45.D0
         DELQ = 90.D0
         IF ( NPOL.EQ.1 .OR. NPOL.EQ.3 ) THEN
            IDREH = +1
         ELSE
            IDREH = -1
         END IF
      CASE DEFAULT
         IF ( IDREH.NE.0 ) IDREH = IDREH/IABS(IDREH)
      END SELECT
C
C     write input data
      WRITE (NOUT1,*) 'input data:'
      WRITE (NOUT1,*) '==========='
      WRITE (NOUT1,99001) EMESH,ESTEP,EFEV,VIH,VIL,OMEV,EREAL,VPREV,THQ,
     &                    FIQ
      WRITE (NOUT1,99002) ALQ,DELQ
      DO IT = 1,LAYS
         DO IAN = 1,NATL(IT)
            IQ = IQ_SPR_LAYAT(IAN,IT)
            DO IO = 1,NOQ(IQ)
               WRITE (NOUT1,99001,ERR=100) 0.0D0
            END DO
            LINESREAD = LINESREAD + 1
         END DO
      END DO
      DO IT = 1,LAYS
         DO IAN = 1,NATL(IT)
            WRITE (NOUT1,99001) 0.0D0
         END DO
      END DO
      WRITE (NOUT1,99003) NOUT1,IP,NREL,IDORA,NSPLEE,NSPIN,SPOL,LAYER1,
     &                    GANZ,LANZ1,LANZ2,IBLOCH
      WRITE (NOUT1,99001) TMIN,TMAX,PMIN,PMAX,TEMP,DTEMP,MASS
      WRITE (NOUT1,99001) POL0(1),POL0(2),POL0(3),POL0L(1),POL0L(2),
     &                    POL0L(3)
      WRITE (NOUT1,99004) NT,NP,TYP,ISTR(1),ISTR(2)
      WRITE (NOUT1,99005) NPOL,ICIRC,IDREH
      WRITE (NOUT1,99001) Q1,Q2,Q3,Q4
      WRITE (NOUT1,99006) CORE,ESTART,EEND
      WRITE (NOUT1,99007) MCD,XMATOFF,USEEULER,NPARA
      WRITE (NOUT1,99012) ASYM
      WRITE (NOUT1,99009) PKSCAN,PKEMIN,PKEMAX
      WRITE (NOUT1,99010) ASIG,EMIN,EMAX,FROMSIGMA,TOSIGMA
      WRITE (NOUT1,99011) VBSTP,VBH,VBL
      WRITE (NOUT1,99013) IXAS,LOWHM,LOIHM,LOEF
      WRITE (NOUT1,99008) IUSEDOS,IDCALC
      WRITE (NOUT1,99013) IGCONV,FWHM,FTEMP
C
      IF ( XMATOFF.EQ.1 ) WRITE (NOUT1,*) 'xmatrix not in use'
      WRITE (NOUT1,*) 'input data end'
      WRITE (NOUT1,*) '=============='
C
      RETURN
C
 100  CONTINUE
      WRITE (*,*) 'error in line ',LINESREAD,' reading input file!!!'
      STOP
C
99001 FORMAT (5E14.7)
99002 FORMAT (3E14.7,i3)
99003 FORMAT (12I3)
99004 FORMAT (5I3)
99005 FORMAT (4I3)
99006 FORMAT (1x,i1,2E14.7)
99007 FORMAT (4I3)
99008 FORMAT (2I3)
99009 FORMAT (i3,2E14.7)
99010 FORMAT (i3,4E14.7)
99011 FORMAT (3E14.7)
99012 FORMAT (1E14.7)
99013 FORMAT (i2,3E14.7)
      END
C*==spec_input.f    processed by SPAG 6.70Rc at 08:30 on 19 Apr 2017
      SUBROUTINE SPEC_INPUT(EMESH,ESTEP,EFEV,VIH,VIL,OMEV,EREAL,VPREV,
     &                      THQ,FIQ,ALQ,NOUT1,IP,NREL,NPOL,NSPLEE,NSPIN,
     &                      SPOL,LAYER1,GANZ,LANZ1,LANZ2,TMIN,TMAX,PMIN,
     &                      PMAX,NT,NP,TYP,ISTR,TEMP,DTEMP,MASS,POL0,
     &                      POL0L,IDORA,IDREH,ICIRC,IFSP,Q1,Q2,Q3,Q4,
     &                      IBLOCH,CORE,ESTART,EEND,ASYM,XMATOFF,PKSCAN,
     &                      PKEMIN,PKEMAX,ASIG,FROMSIGMA,TOSIGMA,
     &                      EMININP,EMAXINP,VBH,VBL,VBSTP,NATL,LAYS,
     &                      DELQ,MCD,USEEULER,NPARA,IXAS,LOWHM,LOIHM,
     &                      LOEF,IUSEDOS,IDCALC,IGCONV,FWHM,FTEMP,
     &                      FIXTHQ,ITXRAY,LCXRAY,NCXRAY,IE_ICST,PSPIN,
     &                      POL0_VEC_TILT,POL0_INITIAL)
C
C   ********************************************************************
C   *                                                                  *
C   *  read in SPR-KKR input file                                      *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:NTMAX
      USE MOD_TRANSPHO_LAYPOT,ONLY:IQ_SPR_LAYAT
      USE MOD_SITES,ONLY:NOQ
      USE MOD_SPEC,ONLY:LAYSM,MPW
      USE MOD_FILES,ONLY:IPRINT,FOUND_SECTION,FOUND_REAL,FOUND_STRING,
     &    FOUND_REAL_ARRAY,FOUND_INTEGER_ARRAY,N_FOUND
      USE MOD_ENERGY,ONLY:NETAB,EMIN,EMAX,EFERMI,EWORK
      USE MOD_CONSTANTS,ONLY:RY_EV,C0,C1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALQ,ASYM,DELQ,DTEMP,EEND,EFEV,EMAXINP,EMESH,EMININP,EREAL,
     &       ESTART,ESTEP,FIQ,FIXTHQ,FROMSIGMA,FTEMP,FWHM,LOEF,LOIHM,
     &       LOWHM,MASS,OMEV,PKEMAX,PKEMIN,PMAX,PMIN,TEMP,THQ,TMAX,TMIN,
     &       TOSIGMA,VBH,VBL,VBSTP,VIH,VIL,VPREV
      INTEGER ASIG,CORE,GANZ,IBLOCH,ICIRC,IDCALC,IDORA,IDREH,IE_ICST,
     &        IFSP,IGCONV,IP,ITXRAY,IUSEDOS,IXAS,LANZ1,LANZ2,LAYER1,
     &        LAYS,MCD,NOUT1,NP,NPARA,NPOL,NREL,NSPIN,NSPLEE,NT,PKSCAN,
     &        SPOL,TYP,USEEULER,XMATOFF
      LOGICAL POL0_VEC_TILT
      COMPLEX*16 Q1,Q2,Q3,Q4
      INTEGER ISTR(2),LCXRAY(NTMAX),NATL(LAYSM),NCXRAY(NTMAX)
      REAL*8 POL0(3),POL0L(3),POL0_INITIAL(3),PSPIN(3)
C
C Local variables
C
      CHARACTER*2 CL,POL_E,POL_PH
      LOGICAL COORD_SYS_POL0_POL0L,UDL,UDP,UDT
      INTEGER IAN,IINP(2),IO,IQ,IT
      REAL*8 RINP(2),RINP3(3)
C
C*** End of declarations rewritten by SPAG
C
C---------------------------------------------------------------------
C Energy dependent input parameters
C---------------------------------------------------------------------
C
C     emesh = number of points in energy mesh
C     estep = next point in energy mesh (stepwidth)
C     efev  = fermi energy in ev
C             startenergy for ups and aes calculations
C     vih   = imaginary part of the potential in ev (final state)
C     vil   = imaginary part of the potential in ev (initial state)
C     vprev = inner potential of the bulk crystal in ev
C     ereal = controls absorption in the mtp (ereal = 1)
C---------------------------------------------------------------------
C     Default values (in eV):
C---------------------------------------------------------------------
      IF ( NSPLEE.NE.0 ) THEN
         EMESH = 200
         ESTEP = (EFERMI-0.5D0)*RY_EV/REAL(EMESH-1.0D0)
         EFEV = EFERMI*RY_EV
         VIH = 1.0D0
         VIL = 0.1D0
      ELSE IF ( NSPLEE.EQ.0 ) THEN
         EMIN = 11.0/RY_EV
         EMAX = 60/RY_EV
         EMESH = 200
         ESTEP = (EMAX-EMIN)*RY_EV/REAL(EMESH-1.0D0)
         EFEV = 1.0D0
         VIH = 0.2D0
         VIL = 0.2D0
         EREAL = 1.0D0
      END IF
C---------------------------------------------------------------------
C
      CALL INPUT_FIND_SECTION('ENERGY',1)
C
      EMESH = NETAB(1)
      IF ( EMESH.EQ.1 ) THEN
         ESTEP = 0.0D0
      ELSE IF ( NSPLEE.NE.0 ) THEN
C
         ESTEP = (EMIN-EMAX)*RY_EV/REAL(EMESH-1.D0)
C
      ELSE IF ( NSPLEE.EQ.0 ) THEN
C
         CALL SECTION_SET_REAL('EMINEV',EMIN,9999D0,0)
         CALL SECTION_SET_REAL('EMAXEV',EMAX,9999D0,0)
C
         ESTEP = (EMAX-EMIN)/REAL(EMESH-1.0D0)
C
      END IF
C
      IF ( IBLOCH.EQ.1 ) THEN
         EFEV = EMIN
      ELSE
         EFEV = EMAX*RY_EV
      END IF
C
      EREAL = EFERMI*RY_EV
C
      CALL SECTION_SET_REAL('IMV_FIN',VIH,9999D0,0)
      IF ( .NOT.FOUND_REAL ) THEN
         CALL SECTION_SET_REAL('VIH',VIH,9999D0,0)
         IF ( .NOT.FOUND_REAL ) THEN
            CALL SECTION_SET_REAL('IMV_FIN_EV',VIH,9999D0,0)
            IF ( FOUND_REAL ) VIH = VIH/RY_EV
         ELSE
            VIH = VIH/RY_EV
         END IF
      END IF
C
      CALL SECTION_SET_REAL('IMV_INI',VIH,9999D0,0)
      IF ( .NOT.FOUND_REAL ) THEN
         CALL SECTION_SET_REAL('VIL',VIL,9999D0,0)
         IF ( .NOT.FOUND_REAL ) THEN
            CALL SECTION_SET_REAL('IMV_INI_EV',VIL,9999D0,0)
            IF ( FOUND_REAL ) VIL = VIL/RY_EV
         ELSE
            VIL = VIL/RY_EV
         END IF
      END IF
      VIL = VIL*RY_EV
      VIH = VIH*RY_EV
C
C        read in energy for NOT.SPLEED
C
      IF ( IBLOCH.NE.1 ) THEN
C
         CALL SECTION_SET_REAL('EWORK',VPREV,9999D0,0)
         UDT = FOUND_REAL
         IF ( .NOT.UDT ) THEN
            CALL SECTION_SET_REAL('VPREV',VPREV,9999D0,0)
            UDP = FOUND_REAL
            IF ( .NOT.UDP ) THEN
               CALL SECTION_SET_REAL('EWORK_EV',VPREV,9999D0,0)
               UDL = FOUND_REAL
               IF ( UDL ) VPREV = VPREV/RY_EV
            ELSE
               VPREV = VPREV/RY_EV
            END IF
         END IF
C-
         IF ( UDT .OR. UDP .OR. UDL ) THEN
            WRITE (*,*) 'reading workfunction from input'
            EWORK = VPREV
         ELSE
            WRITE (*,*) 'workfunction from scf calculation'
         END IF
C
         VPREV = VPREV*RY_EV
C
         VPREV = VPREV + EREAL
C
      END IF
C
C        read in energy for SPLEED
C
      IF ( IBLOCH.EQ.1 ) THEN
C
         CALL SECTION_SET_REAL('EWORK',EWORK,9999D0,0)
         UDT = FOUND_REAL
         IF ( .NOT.UDT ) THEN
            CALL SECTION_SET_REAL('EWORK_EV',EWORK,9999D0,0)
            UDL = FOUND_REAL
            IF ( UDL ) EWORK = EWORK/RY_EV
         END IF
C
         IF ( UDL ) THEN
            WRITE (*,*) 'reading workfunction from input'
         ELSE
            WRITE (*,*) 'workfunction from scf calculation'
         END IF
C
         VPREV = EWORK + EFERMI
C
         VPREV = VPREV*RY_EV
C
      END IF
C
C
      EWORK = EWORK*RY_EV
      IF ( EWORK.LT.0.0D0 ) STOP 
     &                      '<SPEC_INPUT>: Please specify work function'
C
C---------------------------------------------------------------------
C Geometry related input parameters (ELECTRONS)
C---------------------------------------------------------------------
C     tmin,tmax     define the range of polar angles
C     pmin,pmax     define the range of azimut angles
C     nt,np  numbers of angular values for a rotation diagram
C            nt: polar,  np: azimuth
C     typ    crystal coordinats in splout, xpsrun, or upsrun
C           =0: i(e) diagram,
C           =1: rotation diagram,         -> phi scan
C           =2: scattering-angle diagram  -> theta scan
C           =3: orthonormal projection    | 3,4 only for angular resolved
C           =4: stereographic projection  | pe (ups, xps) note: nt=np-> nx,ny
C     istr   beam number (h,k)
C
C     pol0,pol0l = initial pol. ,initial pol. in the laboratory system
C
C     q1,q2,q3,q4 = amplitudes of the photoelectron used in spin
C                   polarized calculations
C
C
C     spol    = 1 unpolarized calculation,
C             = 2 spin-polarized calculation
C             = 4 spin density matrix calculations
C---------------------------------------------------------------------
C     Default values:
C---------------------------------------------------------------------
      IF ( NSPLEE.NE.0 ) THEN
         TMIN = 0.0D0
         TMAX = 0.0D0
         PMIN = 0.0D0
         PMAX = 0.0D0
      ELSE IF ( NSPLEE.EQ.0 ) THEN
         TMIN = 45.0D0
         TMAX = 45.0D0
         PMIN = 270.0D0
         PMAX = 270.0D0
      END IF
C
      NT = 1
      NP = 1
      TYP = 1
C
      ISTR(1) = 0
      ISTR(2) = 0
C
      IF ( NSPLEE.NE.0 ) THEN
         POL0(1:3) = 0.0D0
         POL0L(1:3) = 0.0D0
      ELSE IF ( NSPLEE.EQ.0 ) THEN
         POL0(1:3) = 0.0D0
         POL0_INITIAL(1:3) = 0.0D0
         POL0L(1:3) = 0.0D0
      END IF
C
      SPOL = 2
      IF ( IBLOCH.EQ.1 ) SPOL = 1
C
      POL_E = 'PZ'
      COORD_SYS_POL0_POL0L = .TRUE.
      Q1 = C1
      Q2 = C0
      Q3 = C0
      Q4 = C1
C
      CALL INPUT_FIND_SECTION('SPEC_EL',0)
C
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_SET_REAL_ARRAY('THETA',RINP,N_FOUND,2,0,9999D0,0)
         IF ( FOUND_REAL_ARRAY ) THEN
            TMIN = RINP(1)
            TMAX = RINP(2)
         END IF
         CALL SECTION_SET_REAL_ARRAY('PHI',RINP,N_FOUND,2,0,9999D0,0)
         IF ( FOUND_REAL_ARRAY ) THEN
            PMIN = RINP(1)
            PMAX = RINP(2)
         END IF
C
         CALL SECTION_SET_INTEGER('SPOL',SPOL,9999,0)
         PSPIN(1:3) = 0.0D0
         PSPIN(3) = 1.0D0
         CALL SECTION_SET_REAL_ARRAY('PSPIN',RINP3,N_FOUND,3,0,9999D0,0)
         IF ( FOUND_REAL_ARRAY ) THEN
            PSPIN(1) = RINP3(1)
            PSPIN(2) = RINP3(2)
            PSPIN(3) = RINP3(3)
         END IF
         CALL SECTION_SET_INTEGER('NT',NT,9999,0)
         CALL SECTION_SET_INTEGER('NP',NP,9999,0)
         IF ( ABS(TMIN-TMAX).GT.1D-5 .AND. NT.EQ.1 ) THEN
            STOP '<SPEC_INPUT>: FOR ANGULAR SCAN SPECIFY NT'
         ELSE IF ( ABS(TMIN-TMAX).LE.1D-5 .AND. NT.GT.1 ) THEN
            STOP 
     &          '<SPEC_INPUT>: FOR ANGULAR SCAN SPECIFY RANGE OF ANGLES'
         END IF
         IF ( ABS(PMIN-PMAX).GT.1D-5 .AND. NP.EQ.1 ) THEN
            STOP '<SPEC_INPUT>: FOR ANGULAR SCAN SPECIFY NP'
         ELSE IF ( ABS(PMIN-PMAX).LE.1D-5 .AND. NP.GT.1 ) THEN
C
            CALL SECTION_SET_REAL('BETA1',RINP(1),9999D0,0)
            IF ( .NOT.FOUND_REAL ) THEN
               STOP 
     &          '<SPEC_INPUT>: FOR ANGULAR SCAN SPECIFY RANGE OF ANGLES'
            ELSE
               WRITE (6,*) 'Sample Tilt considered'
            END IF
C
            CALL SECTION_FIND_KEYWORD('POL_TILT',UDT)
C
            IF ( UDT ) THEN
               WRITE (*,*) ''
               WRITE (*,*) 'TILT OF ELECTRON-VECTOR POL0 CONSIDERED'
               WRITE (*,*) ''
               POL0_VEC_TILT = .TRUE.
            END IF
C
         END IF
C
         CALL SECTION_SET_INTEGER('TYP',TYP,9999,0)
C
         CALL SECTION_SET_INTEGER_ARRAY('ISTR',ISTR,N_FOUND,2,0,9999,0)
C
         CALL SECTION_SET_REAL_ARRAY('POL0',POL0,N_FOUND,3,0,9999D0,0)
C
         IF ( FOUND_REAL_ARRAY .AND. 
     &        (ABS(POL0(1))+ABS(POL0(2))+ABS(POL0(3))).GT.0.0D0 ) THEN
            WRITE (*,*) ''
            WRITE (*,*) 
     &             '        >> POL0 coordinate system chosen <<        '
            WRITE (*,*) 
     &            '>> NO transformation to crystal coordinate system <<'
            WRITE (*,*) ''
            POL0_INITIAL = POL0
                              !POL0_INITIAL for TILT calculation
         END IF
         CALL SECTION_SET_REAL_ARRAY('POL0L',POL0L,N_FOUND,3,0,9999D0,0)
         IF ( FOUND_REAL_ARRAY .AND. 
     &        (ABS(POL0L(1))+ABS(POL0L(2))+ABS(POL0L(3))).GT.0.0D0 )
     &        THEN
            WRITE (*,*) ''
            WRITE (*,*) '        >> POL0L coordinate system chosen <<'
            WRITE (*,*) 
     &         '>> transformation to crystal coordinate system POL0G <<'
            WRITE (*,*) ''
         END IF
         IF ( (ABS(POL0(1))+ABS(POL0(2))+ABS(POL0(3))).GT.1.0D-16 .AND. 
     &        (ABS(POL0L(1))+ABS(POL0L(2))+ABS(POL0L(3))).GT.1.0D-16 )
     &        THEN
            WRITE (*,*) '>> ERR: POL0 and POL0L chosen <<'
            STOP
         ELSE IF ( ABS(POL0(1))+ABS(POL0(2))+ABS(POL0(3))
     &             .LT.1.0D-16 .AND. 
     &             (ABS(POL0L(1))+ABS(POL0L(2))+ABS(POL0L(3)))
     &             .LT.1.0D-16 ) THEN
            WRITE (*,*) ''
            WRITE (*,*) '>>  POL0 and POL0L not chosen <<'
            WRITE (*,*) '>>    POL_E must be defined   <<'
            WRITE (*,*) ''
            COORD_SYS_POL0_POL0L = .FALSE.
         END IF
C
         CALL SECTION_SET_STRING('POL_E',POL_E,'9999',0)
C
         IF ( .NOT.FOUND_STRING ) THEN
            CALL SECTION_SET_REAL_ARRAY('Q1',RINP,N_FOUND,2,0,9999D0,0)
            IF ( FOUND_REAL_ARRAY ) Q1 = DCMPLX(RINP(1),RINP(2))
            CALL SECTION_SET_REAL_ARRAY('Q2',RINP,N_FOUND,2,0,9999D0,0)
            IF ( FOUND_REAL_ARRAY ) Q2 = DCMPLX(RINP(1),RINP(2))
            CALL SECTION_SET_REAL_ARRAY('Q3',RINP,N_FOUND,2,0,9999D0,0)
            IF ( FOUND_REAL_ARRAY ) Q3 = DCMPLX(RINP(1),RINP(2))
            CALL SECTION_SET_REAL_ARRAY('Q4',RINP,N_FOUND,2,0,9999D0,0)
            IF ( FOUND_REAL_ARRAY ) Q4 = DCMPLX(RINP(1),RINP(2))
         ELSE IF ( POL_E(1:2).EQ.'PZ' ) THEN
            Q1 = C1
            Q2 = C0
            Q3 = C0
            Q4 = C1
         ELSE IF ( POL_E(1:2).EQ.'PX' ) THEN
            Q1 = DCMPLX(0.707D0,0.0D0)
            Q2 = DCMPLX(0.707D0,0.0D0)
            Q3 = DCMPLX(0.707D0,0.0D0)
            Q4 = DCMPLX(-0.707D0,0.0D0)
         ELSE IF ( POL_E(1:2).EQ.'PY' ) THEN
            Q1 = DCMPLX(0.707D0,0.0D0)
            Q2 = DCMPLX(0.0D0,0.707D0)
            Q3 = DCMPLX(0.707D0,0.0D0)
            Q4 = DCMPLX(0.0D0,-0.707D0)
         ELSE
            STOP '<SPEC_INPUT>: POL_E not known, check input file'
         END IF
C
         IF ( .NOT.UDT .AND. COORD_SYS_POL0_POL0L ) THEN
            WRITE (*,*) ''
            WRITE (*,*) 
     &             '>>           POL_E not considered                <<'
            WRITE (*,*) 
     &             '>> POL0 or POL0L must be defined in SPLEED input <<'
            WRITE (*,*) ''
         ELSE IF ( .NOT.UDT .AND. COORD_SYS_POL0_POL0L ) THEN
            WRITE (*,*) ''
            WRITE (*,*) '>> ERR: POL_E, POL0, POL0L not considered <<'
            WRITE (*,*) '>>      use at least one definition!      <<'
            WRITE (*,*) ''
            STOP
         ELSE IF ( UDT .AND. COORD_SYS_POL0_POL0L ) THEN
            WRITE (*,*) ''
            WRITE (*,*) '>> ERR: POL_E and POL0 (POL0L) considered <<'
            WRITE (*,*) '>>       use POL_E or POL0 (POL0L)        <<'
            WRITE (*,*) ''
            STOP
         END IF
C
      END IF
C---------------------------------------------------------------------
C Geometry related input parameters (PHOTONS)
C---------------------------------------------------------------------
C
C     omev  = photon energy in ev
C     thq,fiq =   polar and azimuthal angles of photon beam
C                 normal incidence thq=180 !!
C     alq   = alignment of polarization vector or pol.ellipsis
C     delq  = phase shift between real and imaginary part of e-vector,
C             delq=90 for circular polarized light !!
C
C     # parameter controlling the photon polarisation and dichroism
C     npol    controls the polarization and dichroism
C             = 0 unpolarized and p-s dichroism for the calculation
C             = 1 p-pol or rcp or elliptical (depends on icirc, etc)
C             = 2 s-pol or lcp or elliptical (depends on icirc, etc)
C             = 3 dichroism (cdad, ldad)
C     icirc   controls polarisation / ellipticity of the photons
C             = 0 elliptically pol. light:  alq, delq arbitrary
C             = 1 linear pol. light:        alq arbitrary, delq=0
C             = 2 circular pol. light:      alq=45, delq = 90
C      (should not depend on alq, but the matrixelements become symmetric)
C             = 3 (icirc=2 with fresnel optics )
C             = 4 fresnel optics
C                 (for all kind of polarizations like icirc=0)
C     idreh   controls helicity of the photons
C             = 1, -1: sigma+, sigma- => right or left
C                              circular polarization
C               (note: lcp, rcp are exchanged in some books)
C             = 0 linearly polarized (equals icirc=1)
C     ifsp    = 0 fixed photon azimuth angle
C             = 1 variable (equal to electron) azimuth angle
C---------------------------------------------------------------------
C     Default values:
C---------------------------------------------------------------------
      IF ( NSPLEE.NE.0 ) THEN
         OMEV = 21.0D0
         THQ = 45.0D0
         FIQ = 90.0D0
         ALQ = 45.0D0
         DELQ = 0.0D0
      ELSE IF ( NSPLEE.EQ.0 ) THEN
         OMEV = 0.0
         THQ = 65
               !dummy for LEED (only for PES)
         FIQ = 90
               !dummy for LEED (only for PES)
         ALQ = 45
               !dummy for LEED (only for PES)
         DELQ = 90
                !dummy for LEED (only for PES)
      END IF
C
      NPOL = 1
      ICIRC = 1
      IDREH = 0
      IFSP = 0
C
      CALL INPUT_FIND_SECTION('SPEC_PH',0)
C
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_SET_REAL('THETA',THQ,9999D0,0)
         CALL SECTION_SET_REAL('PHI',FIQ,9999D0,0)
         CALL SECTION_SET_REAL('EPHOT',OMEV,9999D0,0)
C
         CALL SECTION_SET_REAL('ALQ',ALQ,9999D0,0)
         CALL SECTION_SET_REAL('DELQ',DELQ,9999D0,0)
C
         CALL SECTION_SET_STRING('POL_P',POL_PH,'9999',0)
C
         IF ( .NOT.FOUND_STRING ) THEN
            CALL SECTION_SET_INTEGER('NPOL',NPOL,9999,0)
            CALL SECTION_SET_INTEGER('ICIRC',ICIRC,9999,0)
            CALL SECTION_SET_INTEGER('IDREH',IDREH,9999,0)
         ELSE
            IF ( POL_PH(1:1).EQ.'P' ) THEN
               NPOL = 1
               ICIRC = 1
               IDREH = 1
            ELSE IF ( POL_PH(1:1).EQ.'S' ) THEN
               NPOL = 2
               ICIRC = 1
               IDREH = 1
            ELSE IF ( POL_PH(1:2).EQ.'C+' ) THEN
               NPOL = 1
               ICIRC = 2
               IDREH = 1
            ELSE IF ( POL_PH(1:2).EQ.'C-' ) THEN
               NPOL = 2
               ICIRC = 2
               IDREH = 1
            ELSE
               STOP '<SPEC_INPUT>: POL_P not known, check input file'
            END IF
            CALL SECTION_SET_INTEGER('IFSP',IFSP,9999,0)
            IF ( IFSP.EQ.1 ) THEN
               CALL SECTION_SET_REAL('THETA_FIX',FIXTHQ,9999D0,0)
               IF ( .NOT.FOUND_REAL ) FIXTHQ = THQ
            END IF
         END IF
C
      END IF
C---------------------------------------------------------------------
C Structure related input parameters
C---------------------------------------------------------------------
C
C     layer1  = number of bulk layer in a photoemission calculation
C     ganz    = number of reciprocal lattice vectors
C     lanz1   = | number of layerdoublings for final state
C     lanz2   = | number of layerdoublings for initial state
C---------------------------------------------------------------------
      IF ( NSPLEE.NE.0 ) THEN
         LAYER1 = 30
      ELSE IF ( NSPLEE.EQ.0 ) THEN
         LAYER1 = 1
      END IF
      GANZ = 30
      LANZ1 = 10
      LANZ1 = 10
C
      CALL INPUT_FIND_SECTION('SPEC_STR',0)
C
      IF ( FOUND_SECTION ) THEN
         N_FOUND = 2
         CALL SECTION_SET_INTEGER_ARRAY('N_LAYDBL',IINP,N_FOUND,2,0,
     &                                  9999,0)
         IF ( FOUND_INTEGER_ARRAY ) THEN
            LANZ2 = IINP(2)
            LANZ1 = IINP(1)
         END IF
         CALL SECTION_SET_INTEGER('NLAT_G_VEC',GANZ,9999,0)
         CALL SECTION_SET_INTEGER('N_LAYER',LAYER1,9999,0)
      END IF
C---------------------------------------------------------------------
C General input parameters
C---------------------------------------------------------------------
C
C     nout1   = standard outputkanal
C     ip      = controls amount of output
C             = -1 minimal
C             = 0  few
C             = 1  more essential
C             = 2..4  even more, most, all
C             = 5  more than enaugh
C     nrel    = 0 nonrelativistic calculation,
C             = 1 relativistic calculation
C     idora   controls calculation in upsrun
C             in upsrun (ibloch=4 etc)
C             = 0 paramagnetc calc. (use solvedirac),
C             = 1 ferromagnetic calc. (use doublerad)
C             in bands (ibloch = 1)
C             = 0 use rslab,
C             = 1 use fullf
C     nsplee  = 0 spleed calculation,
C             = 1 photoemission calculation
C     nspin   = 2: has to be 2 (but is principally not longer in use)
C     THESE PARAMETERS ARE NOT RELEVANT ANYMORE
      NREL = 1
      IDORA = 1
      NSPIN = 2
      TEMP = 0.0D0
      DTEMP = 0.0D0
      MASS = 0.0D0
C
C---------------------------------------------------------------------
      IP = IPRINT
C
      CALL INPUT_FIND_SECTION('SPEC',0)
C
      IF ( FOUND_SECTION ) CALL SECTION_SET_INTEGER('SPOL',SPOL,9999,0)
C
C     dimension for angular scan is actually too small (mpw=1)
C     but other fields have to be cleaned up before
C     larger dimensions will work
      IF ( NT.GT.MPW .OR. NP.GT.MPW ) THEN
         WRITE (6,*) '*** dimension for angular scan too small ***'
         WRITE (6,*) 'nt, np =',NT,NP
         WRITE (6,*) '***     change dimensions in parms.h     ***'
         STOP
      END IF
C
C---------------------------------------------------------------------
C     From here on all parameters are dummy and not yet
C     implemented.
C
C---------------------------------------------------------------------
C     # parameters for all core level calculations (xps, xas, xes, aes)
C     core    = 1 calculate core wavefunctions (xps, aes, etc.)
C                 works also with ibloch = 1
C                 calculate corelevel instead of bands
C             = 2 read core wavefunctions from file
C                 has to be 1 for mcd calculations
C     estart  = start energy for xps, aes, etc in ev
C     eend    = stop energy for xps, aes, etc in ev
C               estart, estop define the energy intervall where
C               core level are searched
C
      CORE = 1
      ESTART = -3.0D0
      EEND = -500.0D0
      ITXRAY = 1
      IE_ICST = 1
      LCXRAY(1:NTMAX) = 0
      NCXRAY(1:NTMAX) = 1
C
      CALL INPUT_FIND_SECTION('SPEC_CL',0)
C
      IF ( FOUND_SECTION ) THEN
C
         ITXRAY = 1
         CALL SECTION_SET_INTEGER('IT',ITXRAY,9999,1)
C
         CALL SECTION_SET_STRING('CL',CL,'2P',0)
C
         IF ( FOUND_STRING ) THEN
            NCXRAY(ITXRAY) = ICHAR(CL(1:1)) - ICHAR('1') + 1
            IF ( CL(2:2).EQ.'S' ) LCXRAY(ITXRAY) = 0
            IF ( CL(2:2).EQ.'P' ) LCXRAY(ITXRAY) = 1
            IF ( CL(2:2).EQ.'D' ) LCXRAY(ITXRAY) = 2
            IF ( CL(2:2).EQ.'F' ) LCXRAY(ITXRAY) = 3
         END IF
C
         CALL SECTION_SET_INTEGER('IE_ICST',IE_ICST,9999,0)
C
      END IF
C
C
C     # parameters for spleed, ups, xps, xas, xes
C     mcd         = 0 no mcd
C                 = 1 calculate regular mcd for ups, xps, xas, or xes
C                     (invert b: b -> -b)
C                 = 2 in plane 90deg mcd (example: bx -> by)
C     xmatoff     = 0 use xmatrix
C                 = 1 turn off xmatrix
C     useeuler    = 0 calculate wave functions using full bfield
C                 = 1 rotate bfield and calculate for bz
C                 (note: useeuler will be changed by selecteuler,
C                        if inadequately used)
C     npara       = 0 switch off selection rules
C                 = 1 switch on selection rules
C
      MCD = 1
C
      XMATOFF = 0
      USEEULER = 0
      NPARA = 1
C
C     asym    = asymmetry parameter for doniach-sunjic lineshape
      ASYM = 0.0D0
C
C     # parameters for energy scans
C     pkscan  use peakscan for calculation of energy (result depends on
C     ibloch)
C                 (pkscan is not in use for bands or spleed, ibloch=0,1)
C             = 0 in xps, xes, xas, and direct:
C                    calculate with fixed energymesh in [pkemin,pkemax];
C                 in band calculate with fixed energymesh from 0
C                 in ups calculate with fixed energymesh from efev
C                 in aes calculate with fixed energymesh from efev
C             = 1 in xps calculate mesh of energypoints in the peakregion
C                    finer than in the valley region;
C                 in xes, xas use energyintervall determined by
C                    core-level binding energies and dos
C     in band, ups calculate with fixed energymesh in [pkmin,pkmax]
C                 in directrun: use energymesh from dos
C                    mainly used if called from secondaryrun;
C     = 2 in xps calculate intensities only for the exact eigenvalues
C     within
C                    the energy-interval [pkemin,pkemax];
C     in band calculate with fixed stepwidth in [pkemin,pkemax];
C     in xes, xas use omhar read from inputfile
C     pkemin      |
C     pkemax      | [pkmin,pkmax] energy interval for calculation of
C     energy-vector
C
C
      PKSCAN = 0
      PKEMIN = -0.3D0
      PKEMAX = 0.5D0
C
C     # parameters for all core level calculations (xps, xas, xes, aes)
C     asig        adjust self energy from fromsigma (emin) to tosigma
C     (emax),linearly
C     = 0 constant, use sigma from vil
C     = 1 adjust sigma linearly
C     = 2 adjust sigma = f(j,mj)
C     asig=2 is new, see subroutine detg
C     asig=1 should be corrected if out of range, depending on pkscan
C     emin        |
C     emax        | energy intervalls for self energy interpolation
C     fromsigma   |
C     tosigma     | first and last value  for self energy interpolation
      ASIG = 0
      EMININP = -44.D0
      EMAXINP = -36.0D0
      FROMSIGMA = 1.0D0
      TOSIGMA = 0.5D0
C
C     # parameters for aes
C     vbstp   = stepwidth for valence band in aes
C     vbh     = width of valence band for aes
C     vbl     = limit for valence band in aes
C
      VBSTP = 0.05D0
      VBH = 9.2D0
      VBL = 0.1D0
C
C     # parameters for xas, xes, aes
C     ixas    = method for xas
C                 =0 incoherent, =1 coherent
C     lowhm   = half width at half maximum of greens function
C     loihm   = width about ef
C     loef    = width and height of greens-f. at ef
C               the width of the accompanied lorentzian varies
C               between loef(e=ef=0) and lowhm*pi/2 (|e|>>ef)
C               loihmc determines how fast the transition
C               between these values occurs.
C               lowhm = 0.25d0:  seems to work fine with xas
C               loihm = 1.0d0 :  should be approximately 1.+-0.1 for xas
C               otherwise structures of dos are shifted too much !!!
C               loef = 0.1d0  :  seems to work fine with xas
      IXAS = 0
      LOWHM = 0.1D0
      LOIHM = 1.0D0
      LOEF = 0.01D0
C
C     # parameters for secondaryrun
C     iusedos    = 0  use direct emission spectrum
C                = 1  use dos only (whithout matrix elements)
C     idcalc     = 0  do not run direct emission
C                     read direct emission spectrum from file (r41, l51)
C                     or use dos only (iusedos=1)
C                = 1  calculate direct emission spectrum
      IUSEDOS = 0
      IDCALC = 0
C
C     # parameters for dirac fermi and gauss convolution
C     igconv  parameter for broadening of spectra in experiment
C             =0: off,
C             =1: convolute with gaussian (all spectra)
C             =2: multiply by occupied fermi (ups)
C             =3: multiply by unoccupied fermi (ipes)
C             =4: multiply by fermi (2) and convolute with gaussian (ups)
C             =5: multiply by fermi (3) and convolute with gaussian (ipes)
C     fwhm    = full width at half maximum for gauss convolution
C               is also used in cgauss (called by det_dos)
C               <1e-3 switch off gauss broadening in cgauss (det_dos)
C     ftemp   temperature in fermi-dirac distribution
C             ftemp is also used in det_dos
C
      IGCONV = 2
      FWHM = 0.1D0
      FTEMP = 300D0
C
C---------------------------------------------------------------------
C
C     # parameters for xas, xes, aes
C     dosunit = logical unit numbers for files containing the dos
C               for each atom and layer
C        read (IFILSPECINP, 990,err=300) (dosunit(ian,it),ian=1,natl(it))
C
C
C     check and correct for wrong inputs
C     ==================================
C
C     correct input for use of dos if idcalc=1
      IF ( IDCALC.EQ.1 ) IUSEDOS = 0
C
C     correct input for non-scattering methods
      IF ( IBLOCH.EQ.5 .OR. IBLOCH.EQ.6 .OR. IBLOCH.EQ.7 .OR. 
     &     IBLOCH.EQ.8 ) THEN
         LAYER1 = 1
         GANZ = 1
         LANZ1 = 1
         LANZ2 = 1
      END IF
      IF ( IBLOCH.EQ.8 .AND. MCD.NE.0 .AND. (NPOL.EQ.0 .OR. NPOL.EQ.3) )
     &     THEN
         WRITE (*,*) 
     &     ' mcd only for fixed polarisation possible in secondaryrun !'
         WRITE (*,*) ' program stopped, check method in input   !'
         WRITE (NOUT1,*) 
     &     ' mcd only for fixed polarisation possible in secondaryrun !'
         WRITE (NOUT1,*) ' program stopped, check method in input   !'
         STOP
      END IF
C
      IF ( NSPIN.NE.2 ) THEN
         WRITE (*,*) ' nspin is reset to 2 !'
         WRITE (NOUT1,*) ' nspin is reset to 2 !'
         NSPIN = 2
      END IF
C
C     correct the photon polarisation such that npol, icirc, etc. fit together
      IF ( NPOL.EQ.0 ) ICIRC = 1
      SELECT CASE (ICIRC)
      CASE (1,4)
C         note that s and p are not defined a priori !
C         for thq=0,180 they cannot be distinguished, but s=ey, p=ex
C         for thq/=0 the following gives  s, p for phiq=0, 180,
C                                 and     p, s for phiq=90, 270,
C                                 the notation s,p will depend also on alq for other phq
C         however, the construction of unpolarised light is always correct !
C
         DELQ = 0.D0
         IDREH = 0
         IF ( NPOL.EQ.0 .OR. NPOL.EQ.1 ) THEN
            ALQ = 0.D0
         ELSE IF ( NPOL.EQ.2 ) THEN
            ALQ = 90.D0
         END IF
      CASE (2,3)
         ALQ = 45.D0
         DELQ = 90.D0
         IF ( NPOL.EQ.1 .OR. NPOL.EQ.3 ) THEN
            IDREH = +1
         ELSE
            IDREH = -1
         END IF
      CASE DEFAULT
         IF ( IDREH.NE.0 ) IDREH = IDREH/IABS(IDREH)
      END SELECT
C
C     write input data
      WRITE (NOUT1,*) 'input data:'
      WRITE (NOUT1,*) '==========='
      WRITE (NOUT1,99001) REAL(EMESH),ESTEP,EFEV,VIH,VIL,OMEV,EREAL,
     &                    VPREV,THQ,FIQ
      WRITE (NOUT1,99002) ALQ,DELQ
      DO IT = 1,LAYS
         DO IAN = 1,NATL(IT)
            IQ = IQ_SPR_LAYAT(IAN,IT)
            DO IO = 1,NOQ(IQ)
               WRITE (NOUT1,99001) 0.0D0
            END DO
         END DO
      END DO
      DO IT = 1,LAYS
         DO IAN = 1,NATL(IT)
            WRITE (NOUT1,99001) 0.0D0
         END DO
      END DO
      WRITE (NOUT1,99003) NOUT1,IP,NREL,IDORA,NSPLEE,NSPIN,SPOL,LAYER1,
     &                    GANZ,LANZ1,LANZ2,IBLOCH
      WRITE (NOUT1,99001) TMIN,TMAX,PMIN,PMAX,TEMP,DTEMP,MASS
      WRITE (NOUT1,99001) POL0(1),POL0(2),POL0(3),POL0L(1),POL0L(2),
     &                    POL0L(3)
      WRITE (NOUT1,99004) NT,NP,TYP,ISTR(1),ISTR(2)
      WRITE (NOUT1,99005) NPOL,ICIRC,IDREH
      WRITE (NOUT1,99001) Q1,Q2,Q3,Q4
      WRITE (NOUT1,99006) CORE,ESTART,EEND
      WRITE (NOUT1,99007) MCD,XMATOFF,USEEULER,NPARA
      WRITE (NOUT1,99012) ASYM
      WRITE (NOUT1,99009) PKSCAN,PKEMIN,PKEMAX
      WRITE (NOUT1,99010) ASIG,EMININP,EMAXINP,FROMSIGMA,TOSIGMA
      WRITE (NOUT1,99011) VBSTP,VBH,VBL
      WRITE (NOUT1,99013) IXAS,LOWHM,LOIHM,LOEF
      WRITE (NOUT1,99008) IUSEDOS,IDCALC
      WRITE (NOUT1,99013) IGCONV,FWHM,FTEMP
C      DO IT = 1,LAYS
C         WRITE (NOUT1,99018) (DOSUNIT(IAN,IT),IAN=1,NATL(IT))
C      END DO
C
      IF ( XMATOFF.EQ.1 ) WRITE (NOUT1,*) 'xmatrix not in use'
      WRITE (NOUT1,*) 'input data end'
      WRITE (NOUT1,*) '=============='
C----------------------------------------------------------------------
      WRITE (6,99014)
      WRITE (6,99015) 'energy path (eV):'
      WRITE (6,99017) 'Nr. of energies: ',REAL(EMESH)
      IF ( NSPLEE.NE.0 ) THEN
         WRITE (6,99017) 'E_min :          ',EMIN*RY_EV  !Ausgabe in eV
         WRITE (6,99017) 'E_max :          ',EMAX*RY_EV  !Ausgabe in eV
      ELSE
         WRITE (6,99017) 'E_min :          ',EMIN  !Ausgabe in eV
         WRITE (6,99017) 'E_max :          ',EMAX  !Ausgabe in eV
      END IF
      WRITE (6,99017) 'EWORK :          ',EWORK  !Ausgabe in eV
      WRITE (6,99017) 'Inner Pot :      ',VPREV  !Ausgabe in eV
      WRITE (6,99017) 'Im(V) Ini.:      ',VIL    !Ausgabe in eV
      WRITE (6,99017) 'Im(V) Fin.:      ',VIH    !Ausgabe in eV
C
      IF ( NSPLEE.NE.0 ) THEN
         WRITE (6,99015) 'Geometry parameters (Light):'
         WRITE (6,99017) 'Photon energy:   ',OMEV
         WRITE (6,99018) 'Polarisation:    ',POL_PH
         WRITE (6,99017) 'Polar angle:     ',THQ
         WRITE (6,99017) 'Azimuthal angle: ',FIQ
      ELSE IF ( NSPLEE.EQ.0 ) THEN
         WRITE (6,*) ''
         WRITE (6,*) '--------------------------------------------'
         WRITE (6,*) 'Geometry parameters Light obsolet for SPLEED'
         WRITE (6,*) '--------------------------------------------'
         WRITE (6,*) ''
      END IF
C
      IF ( IFSP.EQ.1 ) THEN
         WRITE (6,99015) 
     &           'Azimuthal angle between light and electrons is fixed.'
         WRITE (6,99017) 'Angle Ele.-Ph. : ',FIXTHQ
      END IF
C
      WRITE (6,99015) 'Geometry parameters (Electrons):'
      IF ( NT.EQ.1 .AND. NP.EQ.1 ) THEN
         WRITE (6,99017) 'Polar angle:     ',TMIN
         WRITE (6,99017) 'Azimuthal angle: ',PMIN
      ELSE
         WRITE (6,99016) 'Nr. of polar ang:',NT
         WRITE (6,99017) 'Polar angle range',TMIN,TMAX
         WRITE (6,99016) 'Nr. of azim. ang:',NP
         WRITE (6,99017) 'Azim. angle range',PMIN,PMAX
      END IF
      IF ( SPOL.EQ.4 ) THEN
         WRITE (6,*) 'Spin density matrix used:'
         WRITE (6,99017) 'Vector of spin polarisation: ',PSPIN(1),
     &                   PSPIN(2),PSPIN(3)
      ELSE IF ( COORD_SYS_POL0_POL0L ) THEN
         WRITE (*,*) '         POL_E not used'
      ELSE IF ( .NOT.COORD_SYS_POL0_POL0L ) THEN
         WRITE (6,99018) 'Spin polarisation: ',POL_E
      END IF
C
C
      WRITE (6,99015) 'Geometry parameters (Structure):'
      WRITE (6,99016) 'Nr. of reciprocal lattice vectors:    ',GANZ
      WRITE (6,99016) 'Total Nr. of layers:                  ',LAYER1
      WRITE (6,99016) 'Nr. of layerdoublings: final state:   ',LANZ1
      WRITE (6,99016) 'Nr. of layerdoublings: initial state: ',LANZ2
C
      WRITE (6,99014)
C----------------------------------------------------------------------
C
      RETURN
C
99001 FORMAT (5E14.7)
99002 FORMAT (3E14.7,i3)
99003 FORMAT (12I3)
99004 FORMAT (5I3)
99005 FORMAT (4I3)
99006 FORMAT (1x,i1,2E14.7)
99007 FORMAT (4I3)
99008 FORMAT (2I3)
99009 FORMAT (i3,2E14.7)
99010 FORMAT (i3,4E14.7)
99011 FORMAT (3E14.7)
99012 FORMAT (1E14.7)
99013 FORMAT (i2,3E14.7)
99014 FORMAT (/,1X,79('-'),//)
99015 FORMAT (//,10X,A)
99016 FORMAT (10X,A,5I10)
99017 FORMAT (10X,A,5F10.6)
99018 FORMAT (10X,A,A)
      END
C*==spec_inpstru.f    processed by SPAG 6.70Rc at 08:30 on 19 Apr 2017
      SUBROUTINE SPEC_INPSTRU(STRFIL,IFILSPECSTR,USE_CRY_PRIM_VECS,
     &                        MILLER_INDICES,IQSURF,TOLZ)
C
C   ********************************************************************
C   *                                                                  *
C   *  Create automatically structure file                             *
C   *  for spectroscopy calculations                                   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LATTICE,ONLY:ALAT,ABAS,SYSTEM_TYPE,SYSTEM_DIMENSION,
     &    BRAVAIS,BDBINV,BBAS,VOLUC
      USE MOD_SITES,ONLY:NQMAX,ITOQ,NOQ,NQ_L,NQ_I,QBAS,NQ
      USE MOD_TYPES,ONLY:Z,NTMAX
      USE MOD_FILES,ONLY:IPRINT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IFILSPECSTR,IQSURF
      CHARACTER*80 STRFIL
      REAL*8 TOLZ
      LOGICAL USE_CRY_PRIM_VECS
      INTEGER MILLER_INDICES(3)
C
C Local variables
C
      REAL*8 AB(:,:),ABAS_NEW(3,3),ADAINV(:,:),ADAMAT(:,:),AV(:),B1(3),
     &       B2(3),BASXY(:,:,:),DEL(:,:),DELR(:,:),DELXY(:,:,:),
     &       DELXYTOP(3),DELZ(:),QBASC(:,:),QBASNEW(:,:),R1,R2,R3,RTMP,
     &       RTMP2(0:2),RTMP3(:),RV(:),TOL,WRK2X2(:,:),Z_LAY(:)
      REAL*8 DDOT
      INTEGER I,IERROR,IFIL,ILAY,ILAYTOP,IO,IOCNTR(:),IQ,IQCNTR(:),
     &        IQ_IQNEW(:),IQ_JQ_Z(:,:),IQ_JQ_ZO(:,:),ISRF,IT,ITMP,
     &        ITMP2(:),J,JQ,K,L,L1,L2,LAYS,NATL(:),NIO_LAY(:),NLAYL,
     &        NLAYR,NLAYTOT
      LOGICAL IQDONE(:),LAYISVAC(:),RUMPLED
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE Z_LAY,IQ_JQ_Z,IQDONE,NIO_LAY,LAYISVAC
      ALLOCATABLE NATL,DELZ,DELXY,BASXY,IQCNTR
      ALLOCATABLE QBASC,AB,RV,AV,ADAMAT,ADAINV,WRK2X2
      ALLOCATABLE IOCNTR,IQ_JQ_ZO,DEL,QBASNEW,IQ_IQNEW
      ALLOCATABLE RTMP3,ITMP2,DELR
      ALLOCATE (Z_LAY(NQMAX),IQ_JQ_Z(NQMAX,NQMAX))
      ALLOCATE (DELR(NQMAX,NQMAX))
      ALLOCATE (IQDONE(NQMAX),NIO_LAY(NQMAX),LAYISVAC(NQMAX))
      ALLOCATE (NATL(NQMAX),DELZ(0:NQMAX))
      ALLOCATE (DELXY(0:NQMAX,0:NQMAX,0:2),DEL(0:NQMAX,1:2))
      ALLOCATE (BASXY(0:NQMAX,0:NQMAX,0:2))
      ALLOCATE (IQCNTR(0:NQMAX))
      ALLOCATE (IOCNTR(0:NQMAX))
      ALLOCATE (QBASC(3,NQMAX))
      ALLOCATE (RV(2),AV(2),ADAINV(2,2),WRK2X2(2,2),AB(2,2))
      ALLOCATE (ADAMAT(2,2),IQ_JQ_ZO(0:NQMAX,NQMAX))
      TOL = 1D-8
C     Tolarance for considering two layers to be rumpled
      IF ( ABS(TOLZ-99999.0D0).LT.1.0D-16 ) TOLZ = 0.1D0
C
      IQDONE(:) = .FALSE.
      LAYISVAC(:) = .FALSE.
      RUMPLED = .FALSE.
C
C
C----------------------------------------------------------------------
      IF ( SYSTEM_DIMENSION(1:2).EQ.'3D' ) THEN
         ALLOCATE (QBASNEW(3,NQMAX))
         ALLOCATE (IQ_IQNEW(NQMAX))
         ALLOCATE (ITMP2(NQMAX))
C
         IERROR = 0
         QBASNEW = 0.0D0
         ABAS_NEW = 0.0D0
         ILAYTOP = 0
C
         CALL CREATE_3D_SURFACE(6,BRAVAIS,USE_CRY_PRIM_VECS,
     &                          MILLER_INDICES,ABAS,BBAS,BDBINV,NQ,QBAS,
     &                          NOQ,ITOQ,VOLUC,IERROR,NQMAX,NTMAX,
     &                          ABAS_NEW,QBASNEW,IQ_IQNEW,IPRINT)
C
         DO I = 1,3
            DO J = 1,3
               IF ( ABS(ABAS_NEW(I,J)).LT.TOL ) ABAS_NEW(I,J) = 0.0D0
            END DO
         END DO
C
C
         IF ( IERROR.NE.0 ) STOP 
     &               '<SPEC_INPUTS>: problems to create 2D from 3D str.'
C
         ILAY = 0
         NIO_LAY(:) = 1
         NLAYL = 0
         NLAYR = 0
         DO IQ = 1,NQ
            IF ( .NOT.IQDONE(IQ) ) THEN
               ILAY = ILAY + 1
               Z_LAY(ILAY) = QBASNEW(3,IQ)
               IQ_JQ_Z(ILAY,1) = IQ
               IO = 1
               DO JQ = IQ,NQ
                  IF ( IQ.NE.JQ .AND. .NOT.IQDONE(JQ) ) THEN
                     B1(1) = QBASNEW(1,JQ) + ABAS_NEW(1,3)
                     B1(2) = QBASNEW(2,JQ) + ABAS_NEW(2,3)
                     B1(3) = QBASNEW(3,JQ) + ABAS_NEW(3,3)
                     B2(1) = QBASNEW(1,JQ) - ABAS_NEW(1,3)
                     B2(2) = QBASNEW(2,JQ) - ABAS_NEW(2,3)
                     B2(3) = QBASNEW(3,JQ) - ABAS_NEW(3,3)
                     R1 = QBASNEW(3,IQ) - QBASNEW(3,JQ)
                     R2 = QBASNEW(3,IQ) - B1(3)
                     R3 = QBASNEW(3,IQ) - B2(3)
                     IF ( ABS(R1).LE.TOLZ ) THEN
                        IO = IO + 1
                        IQDONE(JQ) = .TRUE.
                        IQ_JQ_Z(ILAY,IO) = JQ
                        NIO_LAY(ILAY) = IO
                        DELR(ILAY,IO) = QBASNEW(3,JQ) - QBASNEW(3,IQ)
                     ELSE IF ( ABS(R2).LE.TOLZ .OR. ABS(R3).LE.TOLZ )
     &                         THEN
                        IO = IO + 1
                        IQDONE(JQ) = .TRUE.
                        IQ_JQ_Z(ILAY,IO) = JQ
                        NIO_LAY(ILAY) = IO
                        IF ( ABS(R2).LE.TOLZ ) QBASNEW(1:3,JQ) = B1(1:3)
                        IF ( ABS(R3).LE.TOLZ ) QBASNEW(1:3,JQ) = B2(1:3)
                        DELR(ILAY,IO) = QBASNEW(3,JQ) - QBASNEW(3,IQ)
                     END IF
                     IF ( ABS(DELR(ILAY,IO)).LE.1D-4 ) DELR(ILAY,IO)
     &                    = 0.0D0
                     IF ( ABS(DELR(ILAY,IO)).GT.TOL ) RUMPLED = .TRUE.
                  END IF
               END DO
C
            END IF
         END DO
         NLAYTOT = ILAY
         IF ( IPRINT.GE.1 ) WRITE (6,99003) 
     &                           'Found total numer of layers (in KKR):'
     &                           ,NLAYTOT
C
C     Sort layers with increasing Z
C     Layers are shifted in order to have  IQSURF at Z=0.0
C     If some layers are above top most IQSURF layer, shift
C     if by lattice vector
C
         ALLOCATE (RTMP3(NQMAX))
         DO K = 1,NLAYTOT - 1
            DO L = K + 1,NLAYTOT
               IF ( Z_LAY(K).GT.Z_LAY(L) ) THEN
                  RTMP = Z_LAY(K)
                  Z_LAY(K) = Z_LAY(L)
                  Z_LAY(L) = RTMP
C
                  ITMP = NIO_LAY(K)
                  NIO_LAY(K) = NIO_LAY(L)
                  NIO_LAY(L) = ITMP
C
                  ITMP2(1:NQMAX) = IQ_JQ_Z(K,1:NQMAX)
                  RTMP3(1:NQMAX) = DELR(K,1:NQMAX)
                  IQ_JQ_Z(K,1:NQMAX) = IQ_JQ_Z(L,1:NQMAX)
                  IQ_JQ_Z(L,1:NQMAX) = ITMP2(1:NQMAX)
                  DELR(K,1:NQMAX) = DELR(L,1:NQMAX)
                  DELR(L,1:NQMAX) = RTMP3(1:NQMAX)
               END IF
            END DO
         END DO
         RTMP = 9999999D0
         I = 1
         DO ILAY = 1,NLAYTOT
            IF ( ABS(Z_LAY(ILAY)).LT.TOL ) Z_LAY(ILAY) = 0.0D0
            IF ( IPRINT.GE.1 ) WRITE (6,*) ILAY,Z_LAY(ILAY)
            IF ( Z_LAY(ILAY).LE.RTMP ) THEN
               J = ILAY
               RTMP = Z_LAY(ILAY)
            END IF
            DO I = 1,NIO_LAY(ILAY)
               IQ = IQ_JQ_Z(ILAY,I)
               IF ( IQ_IQNEW(IQ).EQ.IQSURF ) ILAYTOP = ILAY
            END DO
         END DO
         R1 = Z_LAY(ILAYTOP)
         DO ILAY = 1,NLAYTOT
            Z_LAY(ILAY) = Z_LAY(ILAY) - R1
         END DO
C
         IF ( IPRINT.GE.1 ) WRITE (6,99003) 'Layer with lowest Z:',J,
     &                             ILAYTOP
         DO ILAY = 1,NLAYTOT
            IF ( Z_LAY(ILAY).LT.0.0D0 ) THEN
               Z_LAY(ILAY) = Z_LAY(ILAY) + ABAS_NEW(3,3)
               DO I = 1,NIO_LAY(ILAY)
                  IQ = IQ_JQ_Z(ILAY,I)
                  QBASNEW(1,IQ) = QBASNEW(1,IQ) + ABAS_NEW(1,3)
                  QBASNEW(2,IQ) = QBASNEW(2,IQ) + ABAS_NEW(2,3)
               END DO
            END IF
         END DO
C
         DO K = 1,NLAYTOT - 1
            DO L = K + 1,NLAYTOT
               IF ( Z_LAY(K).GT.Z_LAY(L) ) THEN
                  RTMP = Z_LAY(K)
                  Z_LAY(K) = Z_LAY(L)
                  Z_LAY(L) = RTMP
C
                  ITMP = NIO_LAY(K)
                  NIO_LAY(K) = NIO_LAY(L)
                  NIO_LAY(L) = ITMP
C
                  ITMP2(1:NQMAX) = IQ_JQ_Z(K,1:NQMAX)
                  RTMP3(1:NQMAX) = DELR(K,1:NQMAX)
                  IQ_JQ_Z(K,1:NQMAX) = IQ_JQ_Z(L,1:NQMAX)
                  IQ_JQ_Z(L,1:NQMAX) = ITMP2(1:NQMAX)
                  DELR(K,1:NQMAX) = DELR(L,1:NQMAX)
                  DELR(L,1:NQMAX) = RTMP3(1:NQMAX)
               END IF
            END DO
         END DO
C
         IF ( IPRINT.GE.1 ) THEN
            WRITE (6,99002) 'Shifted ILAYTOP to Z=0.0:'
            DO ILAY = 1,NLAYTOT
               DO I = 1,NIO_LAY(ILAY)
                  IQ = IQ_JQ_Z(ILAY,I)
                  WRITE (6,99004) '',QBASNEW(1,IQ),QBASNEW(2,IQ),
     &                            Z_LAY(ILAY),DELR(ILAY,I),
     &                            DFLOAT(IQ_IQNEW(IQ))
               END DO
               WRITE (6,99004) ''
            END DO
         END IF
C
         AB(1,1) = ABAS_NEW(1,1)
         AB(1,2) = ABAS_NEW(1,2)
         AB(2,1) = ABAS_NEW(2,1)
         AB(2,2) = ABAS_NEW(2,2)
C
         ADAMAT = 0.0D0
         DO I = 1,2
            DO J = 1,2
               ADAMAT(I,J) = DDOT(2,ABAS_NEW(1,I),1,ABAS_NEW(1,J),1)
            END DO
         END DO
C
         WRK2X2(1:2,1:2) = ADAMAT(1:2,1:2)
C
         CALL RMATINV(2,2,WRK2X2,ADAINV)
C
         IF ( IPRINT.GE.1 ) WRITE (6,99002)
     &                              ' Crystalographic coordinates:'
         DO ILAY = 1,NLAYTOT
C convert to crystalographic coordinates
            DO IO = 1,NIO_LAY(ILAY)
               IQ = IQ_JQ_Z(ILAY,IO)
               AV(1:2) = QBASNEW(1:2,IQ)
               RV(1:2) = 0.0D0
               CALL RVECEXPAND2D(AV,AB,ADAINV,RV)
               QBASNEW(1:2,IQ) = RV(1:2)
               IF ( IPRINT.GE.1 ) WRITE (6,99004) '',QBASNEW(1,IQ),
     &              QBASNEW(2,IQ)
            END DO
            IF ( IPRINT.GE.1 ) WRITE (6,99003)
         END DO
         AV(1) = ABAS_NEW(1,3)
         AV(2) = ABAS_NEW(2,3)
         RV(1:2) = 0.0D0
         CALL RVECEXPAND2D(AV,AB,ADAINV,RV)
         ABAS_NEW(1,3) = RV(1)
         ABAS_NEW(2,3) = RV(2)
         IFIL = IFILSPECSTR
         IF ( IPRINT.GE.1 ) WRITE (6,*) 'ABAS_NEW(3,:) Cryst.',
     &                                  ABAS_NEW(1:3,3)
         LAYS = NLAYTOT
         NATL(1:LAYS) = NIO_LAY(1:LAYS)
C
         J = 1
         IQCNTR(:) = 0
         IOCNTR(:) = 0
         DELXYTOP(1:2) = 0.0D0
         L1 = 1
         L2 = LAYS
         DO ILAY = L1,L2
C
C     Find distance to the origin for each atom in layer
C     and find out which atom is sitting at origin
C     In the case that top in the most layer no atom at origin is found
C     all atoms in all layers are shifted in order to have one atom at
C     top most layer in origin.
C     If such is not found make corresponding shift of coordinate system
C
            DO IO = 1,NIO_LAY(ILAY)
               IQ = IQ_JQ_Z(ILAY,IO)
               DELXY(J,IO,0) = SQRT(QBASNEW(1,IQ)**2+QBASNEW(2,IQ)**2)
               IF ( DELXY(J,IO,0).LE.TOL ) THEN
                  IQCNTR(J) = IQ
                  IOCNTR(J) = IO
               END IF
            END DO
            IF ( IQCNTR(J).EQ.0 .AND. ILAY.EQ.1 ) THEN
               RTMP = 999999D0
               DO IO = 1,NIO_LAY(ILAY)
                  RTMP = MIN(RTMP,DELXY(J,IO,0))
               END DO
               I = 0
               DO IO = 1,NIO_LAY(ILAY)
                  IQ = IQ_JQ_Z(ILAY,IO)
                  IF ( ABS(RTMP-DELXY(J,IO,0)).LE.TOL ) I = IQ
                  IF ( I.EQ.IQ ) EXIT
               END DO
               DELXYTOP(1:2) = QBASNEW(1:2,I)
            END IF
            J = J + 1
         END DO
C
         J = 1
         IF ( IPRINT.GE.1 ) WRITE (6,99002) 
     &               ' Crystalographic coordinates (shifted to origin):'
         DO ILAY = L1,L2
            DO IO = 1,NIO_LAY(ILAY)
               IQ = IQ_JQ_Z(ILAY,IO)
               QBASNEW(1:2,IQ) = QBASNEW(1:2,IQ) - DELXYTOP(1:2)
               DELXY(J,IO,0) = SQRT(QBASNEW(1,IQ)**2+QBASNEW(2,IQ)**2)
               IF ( IPRINT.GE.1 ) WRITE (6,99004) '',QBASNEW(1,IQ),
     &              QBASNEW(2,IQ)
            END DO
            IF ( IPRINT.GE.1 ) WRITE (6,99003)
            J = J + 1
         END DO
C
         AB(1,1) = 1.0D0
         AB(1,2) = 0.0D0
         AB(2,1) = 0.0D0
         AB(2,2) = 1.0D0
C
         IF ( IPRINT.GE.1 ) WRITE (6,99002) 
     &               ' Crystalographic coordinates (folded back      ):'
         DO ILAY = L1,L2
            DO IO = 1,NIO_LAY(ILAY)
               IQ = IQ_JQ_Z(ILAY,IO)
               AV(1:2) = QBASNEW(1:2,IQ)
               CALL KFOLD(AV,RV,AB)
               QBASNEW(1:2,IQ) = RV(1:2)
               IF ( IPRINT.GE.1 ) WRITE (6,99004) '',QBASNEW(1,IQ),
     &              QBASNEW(2,IQ)
            END DO
            IF ( IPRINT.GE.1 ) WRITE (6,99003)
         END DO
C
         J = 1
         IQCNTR(:) = 0
         IOCNTR(:) = 0
C
         DO ILAY = L1,L2
C
            DO IO = 1,NIO_LAY(ILAY)
               IQ = IQ_JQ_Z(ILAY,IO)
               DELXY(J,IO,0) = SQRT(QBASNEW(1,IQ)**2+QBASNEW(2,IQ)**2)
               IF ( DELXY(J,IO,0).LE.TOL ) THEN
                  IQCNTR(J) = IQ
                  IOCNTR(J) = IO
               END IF
            END DO
            IF ( IQCNTR(J).EQ.0 .AND. ILAY.EQ.1 )
     &            STOP 'Top most layer does not have atom at origin'
            IF ( IQCNTR(J).EQ.0 ) THEN
               IF ( IPRINT.GE.1 ) WRITE (6,99002)
     &               'For layer ILAY no atom at orign found:'
               RTMP = 999999D0
               DO IO = 1,NIO_LAY(ILAY)
                  RTMP = MIN(RTMP,DELXY(J,IO,0))
               END DO
               DO IO = 1,NIO_LAY(ILAY)
                  IQ = IQ_JQ_Z(ILAY,IO)
                  IF ( ABS(RTMP-DELXY(J,IO,0)).LE.TOL ) THEN
                     IQCNTR(J) = IQ
                     IOCNTR(J) = IO
                  END IF
                  IF ( IQCNTR(J).EQ.IQ ) EXIT
               END DO
               IF ( IPRINT.GE.1 ) WRITE (6,99003)
     &               'Layer: IQ atom to be shifted to origin',ILAY,
     &              IQCNTR(J)
               DO IO = 1,NIO_LAY(ILAY)
                  IQ = IQ_JQ_Z(ILAY,IO)
                  DELXY(J,IO,1) = QBASNEW(1,IQCNTR(J))
                  DELXY(J,IO,2) = QBASNEW(2,IQCNTR(J))
                  BASXY(J,IO,1) = QBASNEW(1,IQ) - QBASNEW(1,IQCNTR(J))
                  BASXY(J,IO,2) = QBASNEW(2,IQ) - QBASNEW(2,IQCNTR(J))
                  IF ( IPRINT.GE.1 ) WRITE (6,99004) 'New position',
     &                 BASXY(J,IO,1),BASXY(J,IO,2)
               END DO
            ELSE
               DO IO = 1,NIO_LAY(ILAY)
                  IQ = IQ_JQ_Z(ILAY,IO)
                  DELXY(J,IO,1:2) = 0.0D0
                  BASXY(J,IO,1:2) = QBASNEW(1:2,IQ)
               END DO
            END IF
C
            IF ( IQCNTR(J).EQ.0 ) STOP 
     &                         '<SPEC_INPSTRU>: no site found at origin'
C
            J = J + 1
         END DO
C
C
         IF ( IPRINT.GE.1 ) WRITE (6,99002) 
     &               ' Crystalographic coordinates (folded back      ):'
         J = 1
         DO ILAY = L1,L2
            DO IO = 1,NIO_LAY(ILAY)
               IQ = IQ_JQ_Z(ILAY,IO)
               AV(1:2) = BASXY(J,IO,1:2)
               CALL KFOLD(AV,RV,AB)
               BASXY(J,IO,1:2) = RV(1:2)
               IF ( IPRINT.GE.1 ) WRITE (6,99004) '',BASXY(J,IO,1),
     &              BASXY(J,IO,2)
            END DO
            IF ( IPRINT.GE.1 ) WRITE (6,99003)
            J = J + 1
         END DO
C
         J = 1
         DO ILAY = L1,L2
C
C     resort atoms in a layer, atom in origin should be first one
C
            IF ( IOCNTR(J).NE.1 ) THEN
               RTMP2(0:2) = DELXY(J,1,0:2)
               DELXY(J,1,0:2) = DELXY(J,IOCNTR(J),0:2)
               DELXY(J,IOCNTR(J),0:2) = RTMP2(0:2)
C
               RTMP2(0:2) = BASXY(J,1,0:2)
               BASXY(J,1,0:2) = BASXY(J,IOCNTR(J),0:2)
               BASXY(J,IOCNTR(J),0:2) = RTMP2(0:2)
C
               RTMP3(1) = DELR(ILAY,1)
               DELR(ILAY,1) = DELR(ILAY,IOCNTR(J))
               DELR(ILAY,IOCNTR(J)) = RTMP3(1)
C
               I = IQ_JQ_Z(ILAY,1)
               IQ_JQ_Z(ILAY,1) = IQ_JQ_Z(ILAY,IOCNTR(J))
               IQ_JQ_Z(ILAY,IOCNTR(J)) = I
C
            END IF
            DO IO = 1,NIO_LAY(ILAY)
               IQ_JQ_ZO(J,IO) = IQ_JQ_Z(ILAY,IO)
            END DO
            J = J + 1
         END DO
C
C Take care of origin shift going from one to another layer
         J = 1
         DEL(:,:) = 0.0D0
         DO ILAY = L1,L2 - 1
            DEL(J,1) = DELXY(J+1,1,1) - DELXY(J,1,1)
            DEL(J,2) = DELXY(J+1,1,2) - DELXY(J,1,2)
            J = J + 1
         END DO
         J = L2
         DEL(J,1) = ABAS_NEW(1,3) - DELXY(J,1,1)
         DEL(J,2) = ABAS_NEW(2,3) - DELXY(J,1,2)
C
         DO ILAY = 1,LAYS - 1
            DELZ(ILAY) = Z_LAY(ILAY+1) - Z_LAY(ILAY)
         END DO
         DELZ(LAYS) = ABAS_NEW(3,3) - Z_LAY(LAYS)
C
         WRITE (6,99002) 'File for layered structure will be created:'
         WRITE (6,99005) 'STRFIL=',STRFIL
         WRITE (6,99003) 'UNIT  =',IFIL
C
         REWIND IFIL
         WRITE (IFIL,'(4i5)') 1,1,1,1
         WRITE (IFIL,'(e14.6)') ALAT
         WRITE (IFIL,'(2e14.6)') ABAS_NEW(1,1),ABAS_NEW(2,1)
         WRITE (IFIL,'(2e14.6)') ABAS_NEW(1,2),ABAS_NEW(2,2)
         WRITE (IFIL,'(i5)') LAYS
C
         DO ILAY = 1,LAYS
            WRITE (IFIL,'(i5)') NATL(ILAY)
            DO J = 1,NATL(ILAY)
               WRITE (IFIL,'(4i5)') J,IQ_IQNEW(IQ_JQ_Z(ILAY,J))
               WRITE (IFIL,'(3e14.6)') DELR(ILAY,J),BASXY(ILAY,J,1),
     &                                 BASXY(ILAY,J,2)
            END DO
         END DO
         WRITE (IFIL,'(4i5)') LAYS,1
         WRITE (IFIL,'(3e14.6)') 0.0D0,0.0D0,0.0D0
         DO ILAY = 1,LAYS
            WRITE (IFIL,'(i5)') ILAY
            WRITE (IFIL,'(3e14.6)') DELZ(ILAY),DEL(ILAY,1),DEL(ILAY,2)
         END DO
         IF ( RUMPLED ) THEN
            WRITE (6,99001)
            WRITE (6,99002) 'Rumpled layers found and set'
            WRITE (6,99004) 
     &     'Tolerance for automatic detection of rumpled  layers:(a.u.)'
     &     ,TOLZ
            WRITE (6,99001)
         END IF
      ELSE
C  Automatic structure detection for LIV systems
         AB(1,1) = ABAS(1,1)
         AB(1,2) = ABAS(1,2)
         AB(2,1) = ABAS(2,1)
         AB(2,2) = ABAS(2,2)
C
         ADAMAT = 0.0D0
         DO I = 1,2
            DO J = 1,2
               ADAMAT(I,J) = DDOT(2,AB(1,I),1,AB(1,J),1)
            END DO
         END DO
C
         WRK2X2(1:2,1:2) = ADAMAT(1:2,1:2)
C
         CALL RMATINV(2,2,WRK2X2,ADAINV)
C
C     Find out which IQ are sitting on the same layer
C     NLAYTOT - total number of layers with same Z
C     NIO_LAY -- Number of atoms per layer
C     IQ_JQ_Z(ILAY,IO) -- IQ of layer ilay and atom io in layer
C     Z_LAY -- Z of layer ILAY
C     ! Z should increase from L --> I --> R Zone !
C
         ILAY = 0
         NIO_LAY(:) = 1
         NLAYL = 0
         NLAYR = 0
         DO IQ = 1,NQ
            IF ( .NOT.IQDONE(IQ) ) THEN
               ILAY = ILAY + 1
               Z_LAY(ILAY) = QBAS(3,IQ)
               IQ_JQ_Z(ILAY,1) = IQ
               IO = 1
               IF ( IQ.LE.NQ_L ) NLAYL = NLAYL + 1
               IF ( IQ.GT.(NQ_L+NQ_I) ) NLAYR = NLAYR + 1
               DO JQ = IQ,NQ
                  IF ( IQ.NE.JQ ) THEN
                     IF ( ABS(QBAS(3,IQ)-QBAS(3,JQ)).LE.TOL ) THEN
                        IO = IO + 1
                        IQDONE(JQ) = .TRUE.
                        IQ_JQ_Z(ILAY,IO) = JQ
                        NIO_LAY(ILAY) = IO
                     END IF
                  END IF
               END DO
C
            END IF
         END DO
         NLAYTOT = ILAY
         WRITE (6,99003) 'Found total numer of layers (in KKR):',NLAYTOT
         DO ILAY = 1,NLAYTOT
            WRITE (6,99003) 'Layer:',ILAY
            WRITE (6,99003) '     Sites IQ:',
     &                      (IQ_JQ_Z(ILAY,IO),IO=1,NIO_LAY(ILAY))
         END DO
         WRITE (6,99003) 'Layers in Left Bulk:',NLAYL
         WRITE (6,99003) 'Layers in Right Bulk:',NLAYR
C
C     Find out which layer contains only Vacuum
C
         DO ILAY = 1,NLAYTOT
            JQ = NIO_LAY(ILAY)
            J = 0
            DO I = 1,NIO_LAY(ILAY)
               IQ = IQ_JQ_Z(ILAY,I)
               JQ = JQ + (NOQ(IQ)-1)
               DO IO = 1,NOQ(IQ)
                  IT = ITOQ(IO,IQ)
                  IF ( Z(IT).EQ.0 ) J = J + 1
               END DO
            END DO
            IF ( J.EQ.JQ ) LAYISVAC(ILAY) = .TRUE.
         END DO
         WRITE (6,99003) 'Layer is Vacuum                      '
         DO ILAY = 1,NLAYTOT
            WRITE (6,99003) 'Layer:',ILAY
            WRITE (6,99006) '              ',LAYISVAC(ILAY)
         END DO
C
         IF ( SYSTEM_TYPE(1:3).EQ.'LIV' ) THEN
C
            ISRF = 1
            DO ILAY = 1,NLAYTOT
               IF ( .NOT.LAYISVAC(ILAY) ) ISRF = ILAY + 1
            END DO
C
C        First NQ_L layers above Left region are considered to
C        represent bulk potentials.
C
            WRITE (6,99002) 
     &                 'Warning: for spec calculations poptentials from'
            WRITE (6,99003) 'L-zone are not used.'
            WRITE (6,99003) 'Semi-infinite bulk is created from first',
     &                      NLAYL
            WRITE (6,99003) 
     &               'layers above original L-zone by layer double-ing.'
            WRITE (6,99003) 
     &               'Make sure that interaction zone is thick enought.'
            WRITE (6,99003)
C
            LAYS = ISRF - 1 - NLAYL
            L2 = ISRF - 1
            L1 = NLAYL + 1
            WRITE (6,*) 'Considering layers from',L1,'to ',L2
            J = 0
            DO ILAY = L1,L2
               J = J + 1
               NATL(J) = NIO_LAY(ILAY)
            END DO
C
C        Calculate distance between layers
C
            J = 0
            DO ILAY = L1,L2
               DELZ(J) = Z_LAY(ILAY) - Z_LAY(ILAY-1)
               WRITE (6,99004) 'Del_Z',DELZ(J)
               J = J + 1
            END DO
C
C        In every layer one atom has to be at origin
C
C
            WRITE (6,99002) ' Crystalographic coordinates:'
            DO ILAY = L1 - 1,L2
C convert to crystalographic coordinates
               DO IO = 1,NIO_LAY(ILAY)
                  IQ = IQ_JQ_Z(ILAY,IO)
                  AV(1:2) = QBAS(1:2,IQ)
                  RV(1:2) = 0.0D0
                  CALL RVECEXPAND2D(AV,AB,ADAINV,RV)
                  QBASC(1:2,IQ) = RV(1:2)
                  WRITE (6,99004) '',QBASC(1,IQ),QBASC(2,IQ)
               END DO
               WRITE (6,99003)
            END DO
C
            J = 0
            IQCNTR(:) = 0
            IOCNTR(:) = 0
            DELXYTOP(1:2) = 0.0D0
            DO ILAY = L1 - 1,L2
C
C     Find distance to the origin for each atom in layer
C     and find out which atom is sitting at origin
C     In the case that top in the most layer no atom at origin is found
C     all atoms in all layers are shifted in order to have one atom at
C     top most layer in origin.
C     If such is not found make corresponding shift of coordinate system
C
               DO IO = 1,NIO_LAY(ILAY)
                  IQ = IQ_JQ_Z(ILAY,IO)
                  DELXY(J,IO,0) = SQRT(QBASC(1,IQ)**2+QBASC(2,IQ)**2)
                  IF ( DELXY(J,IO,0).LE.TOL ) THEN
                     IQCNTR(J) = IQ
                     IOCNTR(J) = IO
                  END IF
               END DO
               IF ( IQCNTR(J).EQ.0 .AND. ILAY.EQ.L2 ) THEN
                  RTMP = 999999D0
                  DO IO = 1,NIO_LAY(ILAY)
                     RTMP = MIN(RTMP,DELXY(J,IO,0))
                  END DO
                  I = 0
                  DO IO = 1,NIO_LAY(ILAY)
                     IQ = IQ_JQ_Z(ILAY,IO)
                     IF ( ABS(RTMP-DELXY(J,IO,0)).LE.TOL ) I = IQ
                     IF ( I.EQ.IQ ) EXIT
                  END DO
                  DELXYTOP(1:2) = QBASC(1:2,I)
               END IF
               J = J + 1
            END DO
C
            J = 0
            WRITE (6,99002) 
     &               ' Crystalographic coordinates (shifted to origin):'
            DO ILAY = L1 - 1,L2
               DO IO = 1,NIO_LAY(ILAY)
                  IQ = IQ_JQ_Z(ILAY,IO)
                  QBASC(1:2,IQ) = QBASC(1:2,IQ) - DELXYTOP(1:2)
                  DELXY(J,IO,0) = SQRT(QBASC(1,IQ)**2+QBASC(2,IQ)**2)
                  WRITE (6,99004) '',QBASC(1,IQ),QBASC(2,IQ)
               END DO
               WRITE (6,99003)
               J = J + 1
            END DO
C
            AB(1,1) = 1.0D0
            AB(1,2) = 0.0D0
            AB(2,1) = 0.0D0
            AB(2,2) = 1.0D0
C
            WRITE (6,99002) 
     &               ' Crystalographic coordinates (folded back      ):'
            DO ILAY = L1 - 1,L2
               DO IO = 1,NIO_LAY(ILAY)
                  IQ = IQ_JQ_Z(ILAY,IO)
                  AV(1:2) = QBASC(1:2,IQ)
                  CALL KFOLD(AV,RV,AB)
                  QBASC(1:2,IQ) = RV(1:2)
                  WRITE (6,99004) '',QBASC(1,IQ),QBASC(2,IQ)
               END DO
               WRITE (6,99003)
            END DO
C
            J = 0
            IQCNTR(:) = 0
            IOCNTR(:) = 0
C
            DO ILAY = L1 - 1,L2
C
               DO IO = 1,NIO_LAY(ILAY)
                  IQ = IQ_JQ_Z(ILAY,IO)
                  DELXY(J,IO,0) = SQRT(QBASC(1,IQ)**2+QBASC(2,IQ)**2)
                  IF ( DELXY(J,IO,0).LE.TOL ) THEN
                     IQCNTR(J) = IQ
                     IOCNTR(J) = IO
                  END IF
               END DO
               IF ( IQCNTR(J).EQ.0 .AND. ILAY.EQ.L2 )
     &               STOP 'Top most layer does not have atom at origin'
               IF ( IQCNTR(J).EQ.0 ) THEN
                  WRITE (6,99002) 
     &                          'For layer ILAY no atom at orign found:'
                  RTMP = 999999D0
                  DO IO = 1,NIO_LAY(ILAY)
                     RTMP = MIN(RTMP,DELXY(J,IO,0))
                  END DO
                  DO IO = 1,NIO_LAY(ILAY)
                     IQ = IQ_JQ_Z(ILAY,IO)
                     IF ( ABS(RTMP-DELXY(J,IO,0)).LE.TOL ) THEN
                        IQCNTR(J) = IQ
                        IOCNTR(J) = IO
                     END IF
                     IF ( IQCNTR(J).EQ.IQ ) EXIT
                  END DO
                  WRITE (6,99003) 
     &                          'Layer: IQ atom to be shifted to origin'
     &                          ,ILAY,IQCNTR(J)
                  DO IO = 1,NIO_LAY(ILAY)
                     IQ = IQ_JQ_Z(ILAY,IO)
                     DELXY(J,IO,1) = QBASC(1,IQCNTR(J))
                     DELXY(J,IO,2) = QBASC(2,IQCNTR(J))
                     BASXY(J,IO,1) = QBASC(1,IQ) - QBASC(1,IQCNTR(J))
                     BASXY(J,IO,2) = QBASC(2,IQ) - QBASC(2,IQCNTR(J))
                     WRITE (6,99004) 'New position',BASXY(J,IO,1),
     &                               BASXY(J,IO,2)
                  END DO
               ELSE
                  DO IO = 1,NIO_LAY(ILAY)
                     IQ = IQ_JQ_Z(ILAY,IO)
                     DELXY(J,IO,1:2) = 0.0D0
                     BASXY(J,IO,1:2) = QBASC(1:2,IQ)
                  END DO
               END IF
C
               IF ( IQCNTR(J).EQ.0 )
     &               STOP '<SPEC_INPSTRU>: no site found at origin'
C
               J = J + 1
            END DO
C
            WRITE (6,99002) 
     &               ' Crystalographic coordinates (folded back      ):'
            J = 0
            DO ILAY = L1 - 1,L2
               DO IO = 1,NIO_LAY(ILAY)
                  IQ = IQ_JQ_Z(ILAY,IO)
                  AV(1:2) = BASXY(J,IO,1:2)
                  CALL KFOLD(AV,RV,AB)
                  BASXY(J,IO,1:2) = RV(1:2)
                  WRITE (6,99004) '',BASXY(J,IO,1),BASXY(J,IO,2)
               END DO
               WRITE (6,99003)
               J = J + 1
            END DO
            J = 0
            DO ILAY = L1 - 1,L2
C
C     resort atoms in a layer, atom in origin should be first one
C
               IF ( IOCNTR(J).NE.1 ) THEN
                  RTMP2(0:2) = DELXY(J,1,0:2)
                  DELXY(J,1,0:2) = DELXY(J,IOCNTR(J),0:2)
                  DELXY(J,IOCNTR(J),0:2) = RTMP2(0:2)
C
                  RTMP2(0:2) = BASXY(J,1,0:2)
                  BASXY(J,1,0:2) = BASXY(J,IOCNTR(J),0:2)
                  BASXY(J,IOCNTR(J),0:2) = RTMP2(0:2)
                  I = IQ_JQ_Z(ILAY,1)
                  IQ_JQ_Z(ILAY,1) = IQ_JQ_Z(ILAY,IOCNTR(J))
                  IQ_JQ_Z(ILAY,IOCNTR(J)) = I
               END IF
               DO IO = 1,NIO_LAY(ILAY)
                  IQ_JQ_ZO(J,IO) = IQ_JQ_Z(ILAY,IO)
               END DO
               J = J + 1
            END DO
C
C Take care of origin shift going from one to another layer
            J = 0
            DEL(:,:) = 0.0D0
            DO ILAY = L1 - 1,L2 - 1
               DEL(J,1) = DELXY(J,1,1) - DELXY(J+1,1,1)
               DEL(J,2) = DELXY(J,1,2) - DELXY(J+1,1,2)
               J = J + 1
            END DO
C
         ELSE IF ( SYSTEM_TYPE(1:3).EQ.'VIV' ) THEN
            WRITE (6,*) '<SPEC_INPSTRU> VIV system not yet implemented'
         ELSE IF ( SYSTEM_TYPE(1:3).EQ.'LIR' ) THEN
            WRITE (6,*) '<SPEC_INPSTRU> LIR system not yet implemented'
         END IF
C
         IFIL = IFILSPECSTR
         WRITE (6,99002) 'File for layered structure will be created:'
         WRITE (6,99005) 'STRFIL=',STRFIL
         WRITE (6,99003) 'UNIT  =',IFIL
C
         REWIND IFIL
         WRITE (IFIL,'(4i5)') 1,1,1,1
         WRITE (IFIL,'(e14.6)') ALAT
         WRITE (IFIL,'(2e14.6)') ABAS(1,1),ABAS(2,1)
         WRITE (IFIL,'(2e14.6)') ABAS(1,2),ABAS(2,2)
         WRITE (IFIL,'(i5)') LAYS
C
         DO ILAY = (L2-L1+1),1, - 1
            WRITE (IFIL,'(i5)') NATL(ILAY)
            DO J = 1,NATL(ILAY)
               WRITE (IFIL,'(4i5)') J,IQ_JQ_ZO(ILAY,J)
               WRITE (IFIL,'(3e14.6)') 0.0D0,(BASXY(ILAY,J,I),I=1,2)
            END DO
         END DO
         WRITE (IFIL,'(4i5)') (L2-L1+1),(L2-L1+1) - NLAYL + 1
         WRITE (IFIL,'(3e14.6)') 0.0D0,0.0D0,0.0D0
         I = 1
         DO ILAY = (L2-L1),1, - 1
            WRITE (IFIL,'(i5)') I
            WRITE (IFIL,'(3e14.6)') DELZ(ILAY),DEL(ILAY,1),DEL(ILAY,2)
            I = I + 1
         END DO
         WRITE (IFIL,'(i5)') I
         WRITE (IFIL,'(3e14.6)') DELZ(0),DEL(0,1),DEL(0,2)
      END IF
C
C
99001 FORMAT (/,1X,79('-'),//)
99002 FORMAT (//,10X,A)
99003 FORMAT (10X,A,20I10)
99004 FORMAT (10X,A,5F10.6)
99005 FORMAT (10X,A,A)
99006 FORMAT (10X,A,L2)
      END
C*==rvecexpand2d.f    processed by SPAG 6.70Rc at 08:30 on 19 Apr 2017
      SUBROUTINE RVECEXPAND2D(A,B,BBINV,C)
C   ********************************************************************
C   *                                                                  *
C   *   expand A with respect to basis vectors B_i                     *
C   *                                                                  *
C   *                   A = sum(i) c_i * B_i                           *
C   *                                                                  *
C   *   with          c_i =  sum(j)  BBINV(i,j) * (A*B_j)              *
C   *                                                                  *
C   *   A     real vector of dimension 2                               *
C   *   B     real 2x2-matrix containing the basis vectors             *
C   *   BBINV real 2x2-matrix with                                     *
C   *         BBINV = BB^(-1)  and  BB(i,j) = (B_i*B_j)                *
C   *   C     real vector of dimension 2 with expansion coefficients   *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 A(2),B(2,2),BBINV(2,2),C(2)
C
C Local variables
C
      REAL*8 ADOTB(2)
      REAL*8 DDOT
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,2
         ADOTB(I) = DDOT(2,A,1,B(1,I),1)
      END DO
C
      DO I = 1,2
         C(I) = 0D0
         DO J = 1,2
            C(I) = C(I) + BBINV(I,J)*ADOTB(J)
         END DO
      END DO
C
      END
C*==create_3d_surface.f    processed by SPAG 6.70Rc at 08:30 on 19 Apr 2017
C
      SUBROUTINE CREATE_3D_SURFACE(IWR,BRAVAIS,USE_CRY_PRIM_VECS,
     &                             MILLER_INDICES,ABAS,BBAS,BDBINV,NQ,
     &                             QBAS,NOQ,ITOQ,VOLUC,IERROR,NQMAX,
     &                             NTMAX,ABAS_NEW,QBASNEW,IQ_IQNEW,
     &                             IPRINT)
C   ********************************************************************
C   *                                                                  *
C   *  change the 3D basis vectors from the standard setting to a new  *
C   *  setting that has 2 basis vetors parallel to a given surface     *
C   *  specified by the Miller indices (h,k,l)                         *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER BRAVAIS,IERROR,IPRINT,IWR,NQ,NQMAX,NTMAX
      LOGICAL USE_CRY_PRIM_VECS
      REAL*8 VOLUC
      REAL*8 ABAS(3,3),ABAS_NEW(3,3),BBAS(3,3),BDBINV(3,3),QBAS(3,NQ),
     &       QBASNEW(3,NQMAX)
      INTEGER IQ_IQNEW(NQMAX),ITOQ(NTMAX,NQMAX),MILLER_INDICES(3),
     &        NOQ(NQMAX)
C
C Local variables
C
      REAL*8 A1DOTV,A1_NEW,A1_P(3),A1_PHI,A1_PP(3),A1_PPP(3),A1_XY,
     &       ABAS_HKL(3,3),ADAINV_NEW(3,3),ADAMAT_NEW(3,3),AIBJ,AP(3),
     &       AQVEC(3),BBAS_HKL(3,3),BBAS_NEW(3,3),BDBINV_HKL(3,3),
     &       BDBMAT_HKL(3,3),BP(3),CA,CB,CG,MA(3,3),MB(3,3),MBA(3,3),
     &       MG(3,3),MGBA(3,3),N_LNG,N_P(3),N_PHI,N_PP(3),N_PPP(3),
     &       N_TET,N_XY,PVEC(3),QBASTMP(:,:),QVEC_NEW(3),QZ,SA,SB,SG,V,
     &       V23(3),VEC(3),VEC_SEL(:),VMIN,VNORM(:),VNORM_INP(:),VOL,
     &       VOLMIN,VOLUC_HKL,VOLUC_NEW,WRK3X3_HKL(3,3)
      LOGICAL CHECK
      REAL*8 DDOT,DNRM2
      INTEGER I,I1,I2,I3,IMAX,IO,IQ,IQNEW,IQNEW_IQ(:),ITOQNEW(:,:),J,JQ,
     &        K,KQ,MILLER_INDICES_INP(3),NABAS_NEW(3,3),NOQNEW(:),
     &        NVEC(3),NVEC_SEL(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE QBASTMP,IQNEW_IQ,NOQNEW
      ALLOCATABLE ITOQNEW
      ALLOCATABLE VEC_SEL,NVEC_SEL,VNORM,VNORM_INP
C
C
C
C*** End of declarations rewritten by SPAG
C
      IMAX = 20
C
      CHECK = .FALSE.
      ALLOCATE (QBASTMP(3,NQMAX))
      ALLOCATE (NOQNEW(NQMAX))
      ALLOCATE (IQNEW_IQ(NQMAX))
      ALLOCATE (ITOQNEW(NTMAX,NQMAX))
C
      ALLOCATE (VEC_SEL(3),NVEC_SEL(3),VNORM(3),VNORM_INP(3))
C
C-----------------------------------------------------------------------
      IF ( MILLER_INDICES(1).EQ.0 .AND. MILLER_INDICES(2).EQ.0 .AND. 
     &     +MILLER_INDICES(3).EQ.0 ) THEN
         WRITE (6,*) 'all Miller indices are 0 !!!!!!!!!!!!!!!'
         IERROR = 1
         RETURN
      END IF
C
C
      MILLER_INDICES_INP(1:3) = MILLER_INDICES(1:3)
C-----------------------------------------------------------------------
C
      ABAS_HKL(1:3,1:3) = ABAS(1:3,1:3)
      BBAS_HKL(1:3,1:3) = BBAS(1:3,1:3)
C
      IF ( USE_CRY_PRIM_VECS ) THEN
C
         SELECT CASE (BRAVAIS)
C----------------------------------------------------------- P primitive
         CASE (4,8,12)
C
            ABAS_HKL(1:3,1:3) = ABAS(1:3,1:3)
C
C------------------------------------------------------- I body centered
         CASE (6,9,14)
C
            ABAS_HKL(1:3,1) = ABAS(1:3,2) + ABAS(1:3,3)
            ABAS_HKL(1:3,2) = ABAS(1:3,1) + ABAS(1:3,3)
            ABAS_HKL(1:3,3) = ABAS(1:3,1) + ABAS(1:3,2)
C
C------------------------------------------------------- F face centered
         CASE (7,13)
C
            ABAS_HKL(1:3,1) = -ABAS(1:3,1) + ABAS(1:3,2) + ABAS(1:3,3)
            ABAS_HKL(1:3,2) = +ABAS(1:3,1) - ABAS(1:3,2) + ABAS(1:3,3)
            ABAS_HKL(1:3,3) = +ABAS(1:3,1) + ABAS(1:3,2) - ABAS(1:3,3)
C
C----------------------------------------------------- B,C base centered
         CASE (5)
C
            IF ( ABS(ABAS(1,2)).LT.1D-6 .AND. ABS(ABAS(3,2)).LT.1D-6 )
     &           THEN
C                                                            B base (AC)
C
               ABAS_HKL(1:3,1) = +ABAS(1:3,1) + ABAS(1:3,3)
               ABAS_HKL(1:3,2) = +ABAS(1:3,2)
               ABAS_HKL(1:3,3) = -ABAS(1:3,1) + ABAS(1:3,3)
C
            ELSE
C                                                            C base (AB)
C
               ABAS_HKL(1:3,1) = +ABAS(1:3,1) - ABAS(1:3,2)
               ABAS_HKL(1:3,2) = +ABAS(1:3,1) + ABAS(1:3,2)
               ABAS_HKL(1:3,3) = +ABAS(1:3,3)
C
            END IF
C
C------------------------------------------------------------------ else
         CASE (1,2,3,10,11)
C
            ABAS_HKL(1:3,1:3) = ABAS(1:3,1:3)
C
         CASE DEFAULT
C
            STOP 'Bravais lattice not set'
C
         END SELECT
C
         IF ( BRAVAIS.GE.12 .AND. BRAVAIS.LE.14 ) THEN
            ABAS_HKL(1:3,1:3) = 0D0
            DO I = 1,3
               ABAS_HKL(I,I) = 1D0
            END DO
         END IF
C
C---------------------- primitive vectors (BBAS_HKL) of reciprocal space
C
         DO I = 1,3
            I1 = 1 + MOD(I,3)
            I2 = 1 + MOD(I1,3)
            BBAS_HKL(1,I) = ABAS_HKL(2,I1)*ABAS_HKL(3,I2)
     &                      - ABAS_HKL(3,I1)*ABAS_HKL(2,I2)
            BBAS_HKL(2,I) = ABAS_HKL(3,I1)*ABAS_HKL(1,I2)
     &                      - ABAS_HKL(1,I1)*ABAS_HKL(3,I2)
            BBAS_HKL(3,I) = ABAS_HKL(1,I1)*ABAS_HKL(2,I2)
     &                      - ABAS_HKL(2,I1)*ABAS_HKL(1,I2)
         END DO
         VOLUC_HKL = DABS(ABAS_HKL(1,1)*BBAS_HKL(1,1)+ABAS_HKL(2,1)
     &               *BBAS_HKL(2,1)+ABAS_HKL(3,1)*BBAS_HKL(3,1))
C
         BBAS_HKL(1:3,1:3) = BBAS_HKL(1:3,1:3)/VOLUC_HKL
C
         DO I = 1,3
            DO J = 1,3
               BDBMAT_HKL(I,J) = DDOT(3,BBAS_HKL(1,I),1,BBAS_HKL(1,J),1)
            END DO
         END DO
C
         WRK3X3_HKL(1:3,1:3) = BDBMAT_HKL(1:3,1:3)
C
         CALL RMATINV(3,3,WRK3X3_HKL,BDBINV_HKL)
C
         MILLER_INDICES(1:3) = 0
C
         DO K = 1,3
C
            DO I = 1,3
               BP(I) = DDOT(3,BBAS_HKL(1,K),1,BBAS(1,I),1)
            END DO
C
            DO I = 1,3
               AP(I) = 0D0
               DO J = 1,3
                  AP(I) = AP(I) + BDBINV(I,J)*BP(J)
               END DO
               AP(I) = 2*AP(I)
               IF ( (NINT(AP(I))-AP(I)).GT.1D-6 ) WRITE (6,*) '***',K,I,
     &              AP(I)
            END DO
C
            DO I = 1,3
               MILLER_INDICES(K) = MILLER_INDICES(K)
     &                             + MILLER_INDICES_INP(I)*NINT(AP(I))
            END DO
C
         END DO
C
      END IF
C
      CALL RVECLCIB(MILLER_INDICES_INP(1),MILLER_INDICES_INP(2),
     &              MILLER_INDICES_INP(3),BBAS_HKL,VNORM_INP)
C
      CALL RVECLCIB(MILLER_INDICES(1),MILLER_INDICES(2),
     &              MILLER_INDICES(3),BBAS,VNORM)
C
C-----------------------------------------------------------------------
C                            get   A1_new
C-----------------------------------------------------------------------
      IF ( CHECK ) WRITE (6,*) '**************************************'
C
      VMIN = 1D+6
      VEC_SEL(1:3) = -1D+6
C
      DO I1 = -IMAX,IMAX
         NVEC(1) = I1
         DO I2 = -IMAX,IMAX
            NVEC(2) = I2
            DO I3 = -IMAX,IMAX
               NVEC(3) = I3
C
               IF ( (I1*MILLER_INDICES(1)+I2*MILLER_INDICES(2)+I3*
     &              MILLER_INDICES(3)).EQ.0 .AND. 
     &              .NOT.(I1.EQ.0 .AND. I2.EQ.0 .AND. I3.EQ.0) ) THEN
C
                  CALL RVECLCIB(I1,I2,I3,ABAS,VEC)
C
                  V = DNRM2(3,VEC,1)
C
                  IF ( ABS(V-VMIN).LT.1D-8 ) THEN
C
                     IF ( ABS(VEC(1)-VEC_SEL(1)).LT.1D-8 ) THEN
C
                        IF ( ABS(VEC(2)-VEC_SEL(2)).LT.1D-8 ) THEN
C
                           IF ( ABS(VEC(3)-VEC_SEL(3)).LT.1D-8 ) THEN
                              WRITE (6,*) 'for A1: vec = vec_sel !!!!'
                              IERROR = 1
                              RETURN
                           ELSE IF ( VEC(3).GT.VEC_SEL(3) ) THEN
                              NVEC_SEL = NVEC
                              VEC_SEL = VEC
                              VMIN = V
                           END IF
C
                        ELSE IF ( VEC(2).GT.VEC_SEL(2) ) THEN
                           NVEC_SEL = NVEC
                           VEC_SEL = VEC
                           VMIN = V
                        END IF
C
                     ELSE IF ( VEC(1).GT.VEC_SEL(1) ) THEN
                        NVEC_SEL = NVEC
                        VEC_SEL = VEC
                        VMIN = V
                     END IF
C
                  ELSE IF ( V.LT.VMIN ) THEN
                     NVEC_SEL = NVEC
                     VEC_SEL = VEC
                     VMIN = V
                  END IF
C
                  IF ( CHECK ) WRITE (6,99004) I1,I2,I3,VEC,V,VEC_SEL,
     &                                VMIN
C
               END IF
C
            END DO
         END DO
      END DO
C
      ABAS_NEW(1:3,1) = VEC_SEL(1:3)
      NABAS_NEW(1:3,1) = NVEC_SEL(1:3)
      A1_NEW = DNRM2(3,ABAS_NEW(1,1),1)
C
C-----------------------------------------------------------------------
C                            get   A2_new
C-----------------------------------------------------------------------
      IF ( CHECK ) WRITE (6,*) '**************************************'
C
      VMIN = 1D+6
      VEC_SEL(1:3) = -1D+6
C
      DO I1 = -IMAX,IMAX
         NVEC(1) = I1
         DO I2 = -IMAX,IMAX
            NVEC(2) = I2
            DO I3 = -IMAX,IMAX
               NVEC(3) = I3
C
               CALL RVECLCIB(I1,I2,I3,ABAS,VEC)
C
               V = DNRM2(3,VEC,1)
C
               A1DOTV = ABS(DDOT(3,VEC,1,ABAS_NEW(1,1),1))
C
               IF ( (I1*MILLER_INDICES(1)+I2*MILLER_INDICES(2)+I3*
     &              MILLER_INDICES(3)).EQ.0 .AND. ABS(A1DOTV-A1_NEW*V)
     &              .GT.1D-8 ) THEN
C
                  IF ( ABS(V-VMIN).LT.1D-8 ) THEN
C
                     IF ( ABS(VEC(1)-VEC_SEL(1)).LT.1D-8 ) THEN
C
                        IF ( ABS(VEC(2)-VEC_SEL(2)).LT.1D-8 ) THEN
C
                           IF ( ABS(VEC(3)-VEC_SEL(3)).LT.1D-8 ) THEN
                              WRITE (6,*) 'for A3: vec = vec_sel !!!!'
                              IERROR = 1
                              RETURN
                           ELSE IF ( VEC(3).GT.VEC_SEL(3) ) THEN
                              NVEC_SEL = NVEC
                              VEC_SEL = VEC
                              VMIN = V
                           END IF
C
                        ELSE IF ( VEC(2).GT.VEC_SEL(2) ) THEN
                           NVEC_SEL = NVEC
                           VEC_SEL = VEC
                           VMIN = V
                        END IF
C
                     ELSE IF ( VEC(1).GT.VEC_SEL(1) ) THEN
                        NVEC_SEL = NVEC
                        VEC_SEL = VEC
                        VMIN = V
                     END IF
C
                  ELSE IF ( V.LT.VMIN ) THEN
                     NVEC_SEL = NVEC
                     VEC_SEL = VEC
                     VMIN = V
                  END IF
C
                  IF ( CHECK ) WRITE (6,99004) I1,I2,I3,VEC,V,VEC_SEL,
     &                                VMIN
C
               END IF
C
            END DO
         END DO
      END DO
C
      ABAS_NEW(1:3,2) = VEC_SEL(1:3)
      NABAS_NEW(1:3,2) = NVEC_SEL(1:3)
C
C-----------------------------------------------------------------------
C                            get   A3_new
C-----------------------------------------------------------------------
      IF ( CHECK ) WRITE (6,*) '**************************************'
C
      VOLMIN = 1D+6
      VMIN = 1D+6
      VEC_SEL(1:3) = -1D+6
C
      DO I1 = -IMAX,IMAX
         NVEC(1) = I1
         DO I2 = -IMAX,IMAX
            NVEC(2) = I2
            DO I3 = -IMAX,IMAX
               NVEC(3) = I3
C
               IF ( (I1*MILLER_INDICES(1)+I2*MILLER_INDICES(2)+I3*
     &              MILLER_INDICES(3)).NE.0 ) THEN
C
                  CALL RVECLCIB(I1,I2,I3,ABAS,VEC)
C
                  V = DNRM2(3,VEC,1)
C
                  V23(1) = ABAS_NEW(2,2)*VEC(3) - ABAS_NEW(3,2)*VEC(2)
                  V23(2) = ABAS_NEW(3,2)*VEC(1) - ABAS_NEW(1,2)*VEC(3)
                  V23(3) = ABAS_NEW(1,2)*VEC(2) - ABAS_NEW(2,2)*VEC(1)
C
                  VOL = ABS(ABAS_NEW(1,1)*V23(1)+ABAS_NEW(2,1)*V23(2)
     &                  +ABAS_NEW(3,1)*V23(3))
C
                  IF ( ABS(VOL-VOLMIN).LT.1D-8 ) THEN
C
                     IF ( ABS(V-VMIN).LT.1D-8 ) THEN
C
                        IF ( ABS(VEC(3)-VEC_SEL(3)).LT.1D-8 ) THEN
C
                           IF ( ABS(VEC(2)-VEC_SEL(2)).LT.1D-8 ) THEN
C
                              IF ( ABS(VEC(1)-VEC_SEL(1)).LT.1D-8 ) THEN
                                 WRITE (6,*)
     &                                   'for A3: vec = vec_sel !!!!'
                                 IERROR = 1
                                 RETURN
                              ELSE IF ( VEC(1).GT.VEC_SEL(1) ) THEN
                                 NVEC_SEL = NVEC
                                 VEC_SEL = VEC
                                 VMIN = V
                                 VOLMIN = VOL
                              END IF
C
                           ELSE IF ( VEC(2).GT.VEC_SEL(2) ) THEN
                              NVEC_SEL = NVEC
                              VEC_SEL = VEC
                              VMIN = V
                              VOLMIN = VOL
                           END IF
C
                        ELSE IF ( VEC(3).LT.VEC_SEL(3) ) THEN
                           NVEC_SEL = NVEC
                           VEC_SEL = VEC
                           VMIN = V
                           VOLMIN = VOL
                        END IF
C
                     ELSE IF ( V.LT.VMIN ) THEN
                        NVEC_SEL = NVEC
                        VEC_SEL = VEC
                        VMIN = V
                        VOLMIN = VOL
                     END IF
C
                  ELSE IF ( VOL.LT.VOLMIN ) THEN
                     NVEC_SEL = NVEC
                     VEC_SEL = VEC
                     VMIN = V
                     VOLMIN = VOL
                  END IF
C
                  IF ( CHECK ) WRITE (6,99004) I1,I2,I3,VEC,V,VEC_SEL,
     &                                VMIN,VOL,VOLMIN
C
               END IF
C
            END DO
         END DO
      END DO
C
      IF ( DDOT(3,VEC_SEL,1,VNORM,1).LT.0D0 ) THEN
         ABAS_NEW(1:3,3) = -VEC_SEL(1:3)
         NABAS_NEW(1:3,3) = -NVEC_SEL(1:3)
      ELSE
         ABAS_NEW(1:3,3) = VEC_SEL(1:3)
         NABAS_NEW(1:3,3) = NVEC_SEL(1:3)
      END IF
C
C-----------------------------------------------------------------------
C                       adjust sign of A2_new
C-----------------------------------------------------------------------
C
      VEC(1) = ABAS_NEW(2,2)*ABAS_NEW(3,3) - ABAS_NEW(3,2)*ABAS_NEW(2,3)
      VEC(2) = ABAS_NEW(3,2)*ABAS_NEW(1,3) - ABAS_NEW(1,2)*ABAS_NEW(3,3)
      VEC(3) = ABAS_NEW(1,2)*ABAS_NEW(2,3) - ABAS_NEW(2,2)*ABAS_NEW(1,3)
C
      VOLUC_NEW = ABAS_NEW(1,1)*VEC(1) + ABAS_NEW(2,1)*VEC(2)
     &            + ABAS_NEW(3,1)*VEC(3)
C
      IF ( VOLUC_NEW.LT.0D0 ) THEN
         VOLUC_NEW = -VOLUC_NEW
         ABAS_NEW(1:3,2) = -ABAS_NEW(1:3,2)
         NABAS_NEW(1:3,2) = -NABAS_NEW(1:3,2)
      END IF
C
C-----------------------------------------------------------------------
      IF ( CHECK ) WRITE (6,*) '**************************************'
      IF ( IPRINT.GE.1 ) THEN
C
C
         WRITE (6,99005) 'old primitve vectors'
         WRITE (6,99006) 'A1    ',(ABAS(1:3,1)),DNRM2(3,ABAS(1,1),1)
         WRITE (6,99006) 'A2    ',(ABAS(1:3,2)),DNRM2(3,ABAS(1,2),1)
         WRITE (6,99006) 'A3    ',(ABAS(1:3,3)),DNRM2(3,ABAS(1,3),1)
C
         WRITE (6,99005) 'old reciprocal primitve vectors'
         WRITE (6,99006) 'B1    ',(BBAS(1:3,1)),DNRM2(3,BBAS(1,1),1)
         WRITE (6,99006) 'B2    ',(BBAS(1:3,2)),DNRM2(3,BBAS(1,2),1)
         WRITE (6,99006) 'B3    ',(BBAS(1:3,3)),DNRM2(3,BBAS(1,3),1)
C
         WRITE (6,99005) 'new primitve vectors'
         WRITE (6,99007) 'A1    ',(NABAS_NEW(1:3,1)),(ABAS_NEW(1:3,1)),
     &                   DNRM2(3,ABAS_NEW(1,1),1)
         WRITE (6,99007) 'A2    ',(NABAS_NEW(1:3,2)),(ABAS_NEW(1:3,2)),
     &                   DNRM2(3,ABAS_NEW(1,2),1)
         WRITE (6,99007) 'A3    ',(NABAS_NEW(1:3,3)),(ABAS_NEW(1:3,3)),
     &                   DNRM2(3,ABAS_NEW(1,3),1)
C
         WRITE (6,*)
         WRITE (6,99008) MILLER_INDICES_INP,VNORM_INP,MILLER_INDICES,
     &                   VNORM
C
         WRITE (6,99005) 'used old reciprocal primitve vectors'
         WRITE (6,99006) 'B1    ',(BBAS_HKL(1:3,1))
         WRITE (6,99006) 'B2    ',(BBAS_HKL(1:3,2))
         WRITE (6,99006) 'B3    ',(BBAS_HKL(1:3,3))
C
         WRITE (6,99009)
         WRITE (6,99006) 'A1 * N',DDOT(3,ABAS_NEW(1,1),1,VNORM,1)
         WRITE (6,99006) 'A2 * N',DDOT(3,ABAS_NEW(1,2),1,VNORM,1)
         WRITE (6,99006) 'A3 * N',DDOT(3,ABAS_NEW(1,3),1,VNORM,1)
C
         WRITE (6,99010) VOLUC,VOLUC_NEW
      END IF
C
C--------------------------------------------------- check basis vectors
C
      DO I = 1,3
         DO J = 1,3
            AIBJ = DDOT(3,ABAS(1,I),1,BBAS(1,J),1)
            IF ( I.EQ.J ) THEN
               IF ( ABS(AIBJ-1D0).GT.1D-8 ) THEN
                  IERROR = 1
                  WRITE (6,'(2I3,F10.7)') I,J,AIBJ
               END IF
            ELSE IF ( ABS(AIBJ).GT.1D-8 ) THEN
               IERROR = 1
               WRITE (6,'(2i3,f10.7)') I,J,AIBJ
            END IF
         END DO
      END DO
C
C-----------------------------------------------------------------------
C                      find rotation matrix  MGBA
C-----------------------------------------------------------------------
C
      N_LNG = DNRM2(3,VNORM,1)
      N_TET = ACOS(VNORM(3)/N_LNG)
C
      N_XY = SQRT(VNORM(1)**2+VNORM(2)**2)
C
      IF ( ABS(N_XY).LT.1D-8 ) THEN
         N_PHI = 0D0
      ELSE IF ( VNORM(2).GE.0D0 ) THEN
         N_PHI = ACOS(VNORM(1)/N_XY)
      ELSE IF ( VNORM(1).LT.0D0 ) THEN
         N_PHI = PI + ACOS(-VNORM(1)/N_XY)
      ELSE
         N_PHI = 2*PI - ACOS(VNORM(1)/N_XY)
      END IF
C
      CA = COS(N_PHI)
      SA = SIN(N_PHI)
      CB = COS(N_TET)
      SB = SIN(N_TET)
C
      MA(1:3,1:3) = 0D0
      MA(1,1) = +CA
      MA(1,2) = +SA
      MA(2,1) = -SA
      MA(2,2) = +CA
      MA(3,3) = 1D0
C
      MB(1:3,1:3) = 0D0
      MB(1,1) = +CB
      MB(1,3) = -SB
      MB(3,1) = +SB
      MB(3,3) = +CB
      MB(2,2) = 1D0
C
      N_P(1:3) = MATMUL(MA(1:3,1:3),VNORM(1:3))
      A1_P(1:3) = MATMUL(MA(1:3,1:3),ABAS_NEW(1:3,1))
C
      IF ( IPRINT.GE.1 ) THEN
C
         WRITE (6,'(A,3F10.5)') '    '
         WRITE (6,'(A,3F10.5)') '  N     ',VNORM(1:3)
         WRITE (6,'(A,3F10.5)') ' A1     ',ABAS_NEW(1:3,1)
         WRITE (6,'(A,3F10.5)') '    '
         WRITE (6,'(A,3F10.5)') '  N_P   ',N_P
         WRITE (6,'(A,3F10.5)') ' A1_P   ',A1_P
      END IF
C
      N_PP(1:3) = MATMUL(MB(1:3,1:3),N_P(1:3))
      A1_PP(1:3) = MATMUL(MB(1:3,1:3),A1_P(1:3))
C
      IF ( IPRINT.GE.1 ) THEN
C
         WRITE (6,'(A,3F10.5)') '    '
         WRITE (6,'(A,3F10.5)') '  N_PP  ',N_PP
         WRITE (6,'(A,3F10.5)') ' A1_PP  ',A1_PP
      END IF
C
      A1_XY = SQRT(A1_PP(1)**2+A1_PP(2)**2)
C
      IF ( ABS(A1_XY).LT.1D-8 ) THEN
         A1_PHI = 0D0
      ELSE IF ( A1_PP(2).GE.0D0 ) THEN
         A1_PHI = ACOS(A1_PP(1)/A1_XY)
      ELSE IF ( A1_PP(1).LT.0D0 ) THEN
         A1_PHI = PI + ACOS(-A1_PP(1)/A1_XY)
      ELSE
         A1_PHI = 2*PI - ACOS(A1_PP(1)/A1_XY)
      END IF
C
      CG = COS(A1_PHI)
      SG = SIN(A1_PHI)
C
      MG(1:3,1:3) = 0D0
      MG(1,1) = +CG
      MG(1,2) = +SG
      MG(2,1) = -SG
      MG(2,2) = +CG
      MG(3,3) = 1D0
C
      N_PPP(1:3) = MATMUL(MG(1:3,1:3),N_PP(1:3))
      A1_PPP(1:3) = MATMUL(MG(1:3,1:3),A1_PP(1:3))
C
      IF ( IPRINT.GE.1 ) THEN
C
         WRITE (6,'(A,3F10.5)') '    '
         WRITE (6,'(A,3F10.5)') '  N_PPP ',N_PPP
         WRITE (6,'(A,3F10.5)') ' A1_PPP ',A1_PPP
      END IF
C
      MBA(1:3,1:3) = MATMUL(MB(1:3,1:3),MA(1:3,1:3))
      MGBA(1:3,1:3) = MATMUL(MG(1:3,1:3),MBA(1:3,1:3))
C
      N_PPP(1:3) = MATMUL(MGBA(1:3,1:3),VNORM(1:3))
      A1_PPP(1:3) = MATMUL(MGBA(1:3,1:3),ABAS_NEW(1:3,1))
C
      IF ( IPRINT.GE.1 ) THEN
C
         WRITE (6,'(A,3F10.5)') 
     &                       ' ----------------------------------------'
         WRITE (6,'(A,3F10.5)') '  N_PPP ',N_PPP
         WRITE (6,'(A,3F10.5)') ' A1_PPP ',A1_PPP
C
      END IF
C
C-----------------------------------------------------------------------
C                  rotate vectors to new coordinate system
C-----------------------------------------------------------------------
C
      DO I = 1,3
         VEC(1:3) = MATMUL(MGBA(1:3,1:3),ABAS_NEW(1:3,I))
         ABAS_NEW(1:3,I) = VEC(1:3)
      END DO
C
      IF ( IPRINT.GE.1 ) THEN
C
         WRITE (6,99005) 'new primitve vectors ---- rotated'
         WRITE (6,99007) 'A1    ',(NABAS_NEW(1:3,1)),(ABAS_NEW(1:3,1)),
     &                   DNRM2(3,ABAS_NEW(1,1),1)
         WRITE (6,99007) 'A2    ',(NABAS_NEW(1:3,2)),(ABAS_NEW(1:3,2)),
     &                   DNRM2(3,ABAS_NEW(1,2),1)
         WRITE (6,99007) 'A3    ',(NABAS_NEW(1:3,3)),(ABAS_NEW(1:3,3)),
     &                   DNRM2(3,ABAS_NEW(1,3),1)
      END IF
C
      DO I = 1,3
         DO J = 1,3
            ADAMAT_NEW(I,J) = DDOT(3,ABAS_NEW(1,I),1,ABAS_NEW(1,J),1)
         END DO
      END DO
C
      CALL RINVGJ(ADAINV_NEW,ADAMAT_NEW,3,3)
C
C---------------------- primitive vectors (BBAS_NEW) of reciprocal space
C
      DO I = 1,3
         I1 = 1 + MOD(I,3)
         I2 = 1 + MOD(I1,3)
         BBAS_NEW(1,I) = ABAS_NEW(2,I1)*ABAS_NEW(3,I2) - ABAS_NEW(3,I1)
     &                   *ABAS_NEW(2,I2)
         BBAS_NEW(2,I) = ABAS_NEW(3,I1)*ABAS_NEW(1,I2) - ABAS_NEW(1,I1)
     &                   *ABAS_NEW(3,I2)
         BBAS_NEW(3,I) = ABAS_NEW(1,I1)*ABAS_NEW(2,I2) - ABAS_NEW(2,I1)
     &                   *ABAS_NEW(1,I2)
      END DO
      VOLUC = DABS(ABAS_NEW(1,1)*BBAS_NEW(1,1)+ABAS_NEW(2,1)
     &        *BBAS_NEW(2,1)+ABAS_NEW(3,1)*BBAS_NEW(3,1))
C
      BBAS_NEW(1:3,1:3) = BBAS_NEW(1:3,1:3)/VOLUC
C
C
C------------------- basis vectors - rotate and shift into new unit cell
      DO I = 1,3
         IF ( IPRINT.GE.1 ) WRITE (6,'(A,I3,A,3F10.5)') 'I ',I,
     &                             '  ABAS   ',ABAS(1:3,I)
      END DO
C
      DO IQ = 1,NQ
C
         QVEC_NEW(1:3) = MATMUL(MGBA(1:3,1:3),QBAS(1:3,IQ))
         IF ( IPRINT.GE.1 ) THEN
C
            WRITE (6,'(A,I3,A,3F10.5)') 'IQ',IQ
            WRITE (6,'(A,I3,A,3F10.5)') 'IQ',IQ,'  Q_old  ',QBAS(1:3,IQ)
            WRITE (6,'(A,I3,A,3F10.5)') 'IQ',IQ,'  Q_rot  ',QVEC_NEW
         END IF
C
         DO I = 1,3
            AQVEC(I) = DDOT(3,ABAS_NEW(1,I),1,QVEC_NEW,1)
         END DO
         IF ( IPRINT.GE.1 ) WRITE (6,'(A,I3,A,3F10.5)') 'IQ',IQ,
     &                             '  Q*A    ',AQVEC
C
         PVEC(1:3) = MATMUL(ADAINV_NEW(1:3,1:3),AQVEC(1:3))
         IF ( IPRINT.GE.1 ) WRITE (6,'(A,I3,A,3F10.5)') 'IQ',IQ,
     &                             '  P      ',PVEC
C
         DO I = 1,3
            PVEC(I) = PVEC(I) + 1000D0
            PVEC(I) = PVEC(I) - INT(PVEC(I))
            IF ( ABS(PVEC(I)-1D0).LT.1D-8 ) PVEC(I) = 0D0
         END DO
C
         QBASNEW(1:3,IQ) = MATMUL(ABAS_NEW(1:3,1:3),PVEC(1:3))
         IF ( IPRINT.GE.1 ) WRITE (6,'(A,I3,A,3F10.5)') 'IQ',IQ,
     &                             '  Q_new  ',QBASNEW(1:3,IQ)
C
      END DO
C
C-------------------- sort basis vectors QBASNEW with increasing z-value
C
      IQNEW_IQ(1) = 1
      DO IQ = 2,NQ
         QZ = QBASNEW(3,IQ)
         IF ( QZ.LT.QBASNEW(3,IQNEW_IQ(IQ-1)) ) THEN
            DO KQ = 1,(IQ-1)
               IF ( QZ.LT.QBASNEW(3,IQNEW_IQ(KQ)) ) THEN
                  JQ = IQ
                  DO J = KQ,(IQ-1)
                     IQNEW_IQ(JQ) = IQNEW_IQ(JQ-1)
                     JQ = JQ - 1
                  END DO
                  IQNEW_IQ(KQ) = IQ
                  GOTO 100
               END IF
            END DO
         END IF
C
         IQNEW_IQ(IQ) = IQ
 100  END DO
C
      QBASTMP(1:3,1:NQ) = QBASNEW(1:3,1:NQ)
C
      DO IQ = 1,NQ
         IQNEW = IQNEW_IQ(IQ)
         IQ_IQNEW(IQNEW) = IQ
C
         QBASNEW(1:3,IQNEW) = QBASTMP(1:3,IQ)
      END DO
C
      DO IQNEW = 1,NQ
         IQ = IQ_IQNEW(IQNEW)
         IF ( IPRINT.GE.1 ) WRITE (6,*) 'QBASNEW ordered',IQNEW,IQ,
     &                                  QBASNEW(3,IQNEW)
      END DO
C
      DO IQNEW = 1,NQ
C
         IQ = IQ_IQNEW(IQNEW)
C
         NOQNEW(IQNEW) = NOQ(IQ)
         DO IO = 1,NOQ(IQ)
            ITOQNEW(IO,IQNEW) = ITOQ(IO,IQ)
         END DO
C
      END DO
C
      IF ( IPRINT.GE.1 ) THEN
C
C
C=======================================================================
C          write result to      xband_geometry.out
C=======================================================================
C
         WRITE (IWR,99002) 'GEOMETRY:  result for TASK = cr_3D_surface'
         WRITE (IWR,99002) 'IERROR   = ',IERROR
         WRITE (IWR,99002) '   ABAS (new)'
         DO J = 1,3
            WRITE (IWR,99003) J,(ABAS_NEW(I,J),I=1,3)
         END DO
         WRITE (IWR,99002) 'NQ       = ',NQ
         WRITE (IWR,99002) '   IQ     QX(IQ)    QY(IQ)    QZ(IQ)  new'
         DO IQ = 1,NQ
            WRITE (IWR,99003) IQ,(QBASNEW(I,IQ),I=1,3)
         END DO
C
         WRITE (IWR,99002) 'new occupation '
         DO IQ = 1,NQ
            WRITE (IWR,99001) IQ,NOQNEW(IQ)
            WRITE (IWR,99001) (ITOQNEW(IO,IQ),IO=1,NOQNEW(IQ))
         END DO
C
         WRITE (IWR,99002) 'IQ_3DBULK(old)  of  IQ_3DSURF(new) '
         DO IQNEW = 1,NQ
            WRITE (IWR,99001) IQNEW,IQ_IQNEW(IQNEW)
         END DO
      END IF
C
C=======================================================================
C
99001 FORMAT (10I5)
99002 FORMAT (A,3I10,F10.3)
99003 FORMAT (I5,1X,3F18.12)
99004 FORMAT (3I3,3F8.3,2x,f8.3,'  >>> ',3F8.3,2x,3F8.3)
99005 FORMAT (/,10X,A,/)
99006 FORMAT (10X,A,9X,2X,3F12.6,2X,3F12.6)
99007 FORMAT (10X,A,3I3,2X,3F12.6,2X,3F12.6)
99008 FORMAT (10X,'NORM     h  k  l',/,10X,'inp    ',3I3,2X,3F12.6,/,
     &        10X,'old    ',3I3,2X,3F12.6)
99009 FORMAT (/,10X,'dot-products',/)
99010 FORMAT (/,10X,'volume of the unit cell ',/,10X,'old',f26.6,/,10X,
     &        'new',f26.6,/)
C
      END
