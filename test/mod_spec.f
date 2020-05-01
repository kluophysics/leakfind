C*==mod_spec.f    processed by SPAG 6.70Rc at 13:07 on 23 May 2012
      MODULE MOD_SPEC
C   ********************************************************************
C   *                                                                  *
C   *  module to store all parameter needed for SPEC   calculations    *
C   *                                                                  *
C   *  This is former parms.h                                          *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER LAYSM,NATLM,NTCPA,ML,LMSELF,MLM2,MLQ,MLQNAT,MLSP,MLS,MLP,
     &        MQD,MLSQ,NFULLPOT,LL,LG,LH,LB,LX,RSTEP,RSTEPP5,NPM,LAYTM,
     &        MQK,MAXA,NSTATES,MAXIPP,MAXN,MLXPS,MAXEIG,XMAXE,MAXZ,
     &        MAXCOEF,MAXEV,MAXPTS,MAXSPEC,MAXAU,MAXAUVAL,MAXCSTATES,
     &        MAXCORE,NDOS,MAXL_AUGER,MLQ_AUGER,MQD_AUGER,MVB,MEW,MPW,
     &        MHDIM,MHDIM1,MHDIM2,MHSTEP,MHMAX,MHMAX1,MLZP
      PARAMETER (LAYSM=51,NATLM=8,NTCPA=2,ML=4,LMSELF=9,MLM2=(2*ML)-1,
     &           MLQ=ML*ML,MLQNAT=MLQ*NATLM,MLSP=4*MLQ,MLS=2*MLQ,
     &           MLP=ML+1,MQD=2*MLQ,MLSQ=2*MLQ*MLQ,NFULLPOT=1,LL=110,
     &           LG=90,LH=2*LG,LB=4*LG,LX=2*LL,RSTEP=721,
     &           RSTEPP5=RSTEP+5,NPM=3000,LAYTM=1,MQK=MQD*MQD,MAXA=4,
     &           NSTATES=3,MAXIPP=2000,MAXN=3,MLXPS=4,
     &           MAXEIG=2*MLXPS*MLXPS,XMAXE=1301,MAXZ=35000,MAXCOEF=20,
     &           MAXEV=20,MAXPTS=10,MAXSPEC=901,MAXAU=150,
     &           MAXAUVAL=35000,MAXCSTATES=14,
     &           MAXCORE=NATLM*LAYSM*MAXCSTATES**2,NDOS=1100,
     &           MAXL_AUGER=3,MLQ_AUGER=MAXL_AUGER*MAXL_AUGER,
     &           MQD_AUGER=MLQ_AUGER+MLQ_AUGER,MVB=100,MEW=501,MPW=501,
     &           MHDIM=300,MHDIM1=MHDIM+1,MHDIM2=2*MHDIM+1,MHSTEP=5,
     &           MHMAX=MHDIM*MHSTEP,MHMAX1=MHMAX+1,MLZP=36124)
      REAL*8 HARTRE,PI,EPS12,EPS16
      PARAMETER (HARTRE=27.21166D0,
     &           PI=3.1415926535897932384626433832795D0,EPS12=1.0D-12,
     &           EPS16=1.0D-16)
      COMPLEX*16 CZERO,CONE,CHALF,CIMAG
      PARAMETER (CZERO=(0.D0,0.D0),CONE=(1.D0,0.D0),CHALF=(0.5D0,0.D0),
     &           CIMAG=(0.D0,1.D0))
      LOGICAL TRANSP_BAR
C
C Local variables
C
      INTEGER MAXWF,NTPHOMAX,NVFTPHOMAX
      SAVE MAXWF,NTPHOMAX,NVFTPHOMAX,TRANSP_BAR
C
C*** End of declarations rewritten by SPAG
C
C     NOTES:
C
C     for bandstructure calculations:            ml.le.10, emesh.le.2501
C     for photoemission calculations: ll.le.100, ml.le.10, emesh.le.300
C     for leed (spleed) calculations:            ml.le.10
C
C     integer parameters:
C     ===================
C
C     natlm: maximum number of different atomtypes per layer,
C     natlm.ge.1
C     laysm: control parameter for fcc,bcc and hcp-structures,
C     laysm=1 for lattices with single layer
C     laysm=2 for lattice with 2 layer (e.g. hcp-lattice)
C     laysm=3 for lattice with 3 layer (e.g. fcc)
C     laysm=n for lattice with n layer
C     ======================================================
C
C     ml = maximum number of phaseshifts
C     number of orbital momenta ml=(l+1) should be: 5<ml<10
C     ======================================================
C
C     nfullpot=mqd=mls in case of nonspherical potential
C     ==================================================
C
C     ll = maximum number of layers including adsorbate layers
C     lg = maximum number of reziprocal lattice vectors
C     ========================================================
C
C
C     number of states which have to be stored:
C     =========================================
C
C     maximum number of possible
C     couplings in phasefunction ansatz:
C     >= 660, but much more for magnetic case
C     para         maxh   = 8*ml**2 = 4*mqd
C     bz                  = 8*(2*mlq-ml)
C     bx, by, bxy         = 8*(mlq + mlm2**2)
C     bxyz                = 8*(3*mlq-2*ml+1)
C     =======================================
C
C     maximum number of wavefunctions in
C     phasefunction ansatz;
C     for b in z-direction maxwf=4*ml^2
C     for b in x-direction maxwf=4*ml^3 should be enough:
C     4*mqd*ml = 8*ml^3
C     ===================================================
C
C
C     maximum number of nonzero matrixelements in calculation
C     of xps-zmatrix. 10000 should be sufficient for most cases:
C     ==========================================================
C
C     number of different core states per atom:
C     should be 4l+2 for a single l-doublett: e.g.: l=3->14
C     =====================================================
C
C     number of wavefunctions which represent the core states:
C     for b=0                (4l+2) *natlm*laysm
C     for b in z-direction   (8l+2) *natlm*laysm:
C     for b in xy  direction (20l+2)*natlm*laysm
C     for b in xyz direction (24l+2)*natlm*laysm
C     in fullpot much more up to maximum: (4l+2)**2 *natlm*laysm
C                                        = natlm*laysm*maxcstates**2
C     ================================================================
C     maxcore       = 4*24)
C
C     number of energysteps in dos:
C     =============================
C
C     for auger calculations
C     ======================
C
C     maximum number of energypoints over valence band range:
C     =======================================================
C
C     maximum number of different electron angles and energies:
C     =========================================================
C     mpw should be 1 and mew=301 for energy scans
C                         mew >= width of spectrum / resolution
C                    e.g. 251 >= 12eV / 0.05eV (-> res = 50meV)
C                    e.g. 301 >= 30eV / 0.1eV (-> res = 100meV)
C     mew should be 1 and mpw=181...361 for angular scans
C                         mpw >= angular range / resolution
C                    e.g. 181 >= 180deg / 1deg = 90deg / 0.5deg
C
C
C    *           mlzp = (((12*ml + 3)*ml + 7)*ml + 8)
C    *                 *(ml + 1)*ml/60 +1,
C    *           mlzp = (((12*ml + 3)*ml + 7)*ml + 8)
C    *                 *(ml + 1)*ml/60,
C
      END
C*==mod_spec_mesh.f    processed by SPAG 6.70Rc at 13:07 on 23 May 2012
C
      MODULE MOD_SPEC_MESH
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      REAL*8 HM,RADMESH(:,:,:),RADMESH3(:,:,:),RADMESH3M(:,:,:),
     &       RADMESHM(:,:,:)
      SAVE RADMESH,RADMESH3,RADMESH3M,RADMESHM
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE RADMESH,RADMESH3,RADMESH3M,RADMESHM
C
      END
C*==mod_spec_wave.f    processed by SPAG 6.70Rc at 13:07 on 23 May 2012
C
C MODULE FOR WAVEFUNCTIONS
      MODULE MOD_SPEC_WAVE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      COMPLEX*16 WFF(:,:,:,:),WFFM(:,:,:,:,:),WFFMATOM(:,:,:,:,:),
     &           WFFMZ(:,:,:,:,:),WFFX(:,:,:),WFFY(:,:)
      COMPLEX*16 FFI_G1M(:,:,:),FFI_S1M(:,:,:),
     &           FF_G1(:,:,:),FF_G1M(:,:,:),FF_G1MZ(:,:,:),FF_S1(:,:,:),
     &     FF_S1M(:,:,:),FF_S1MZ(:,:,:),FF_S2(:,:,:),FF_S2M(:,:,:)
      INTEGER WFFXM(:,:,:,:),WFFYM(:,:,:)
      SAVE WFF,WFFM,WFFMATOM,WFFMZ,WFFX,WFFXM,WFFY,WFFYM
      SAVE FFI_G1M,FFI_S1M,FF_G1,FF_G1M,FF_G1MZ,FF_S1,FF_S1M,FF_S1MZ,
     &     FF_S2,FF_S2M
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE WFF,WFFX,WFFY,WFFM,WFFMATOM,WFFMZ,WFFXM,WFFYM
      ALLOCATABLE FFI_G1M,FFI_S1M,FF_G1,FF_G1M,FF_G1MZ,FF_S1,FF_S1M,
     &     FF_S1MZ,FF_S2,FF_S2M
C
      END

      MODULE MOD_SPEC_TRANSFO

      INTEGER NB(:),REL(:,:),SB(:)
      SAVE NB,SB,REL

      ALLOCATABLE NB,SB,REL

      END
