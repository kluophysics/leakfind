C*==epath.f    processed by SPAG 6.70Rc at 19:23 on  1 Nov 2014
      SUBROUTINE EPATH(IGRID,EMIN,EMAX,EIMAG,NE,ETAB,WETAB,EILOW,IPRINT,
     &                 NEMAX)
C   ********************************************************************
C   *                                                                  *
C   *  ROUTINE TO CREATE  E - PATH                                     *
C   *                                                                  *
C   *  IGRID= 0    only real energies -- Gauss integration mesh        *
C   *  IGRID= 1    only real energies -- equidistant path              *
C   *  IGRID= 2    rectangular complex path                            *
C   *  IGRID= 3    straight complex path parallel to real axis         *
C   *  IGRID= 4    rectangular complex grid, return to real axis       *
C   *              on log scale - assumes nearest approach to real     *
C   *              axis prior to interpolation is 0.5 mRY.             *
C   *  IGRID= 5    arc in the complex plain                            *
C   *              from point 1 to NE      NE <= 60                    *
C   *  IGRID= 6    standard X-ray mesh from E_Fermi to 3.5 Ry          *
C   *              with variable step size  from  E_Fermi to 1.1 Ry    *
C   *  IGRID= 7    X-ray mesh                                          *
C   *              with Gaussian weights for integration               *
C   *  IGRID= 8    arc in the complex plain                            *
C   *              from point 1 to NE  with logaritmic distribution    *
C   *  IGRID= 9    arc in the complex plain                            *
C   *              from point 1 to NE  with Akai's scheme              *
C   *  IGRID=10    ellipse in the complex plain                        *
C   *              from point 1 to NE, DK                              *
C   *  IGRID=11    rectengular path including temperature              *
C   *                                                                  *
C   *                                                                  *
C   *  EIMAG= distance from real axis                                  *
C   *  EMIN = lowest   real(E)                                         *
C   *  EMAX = highest  real(E)                                         *
C   *                                                                  *
C   *  HE 00/03/17                                                     *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:CI,PI
      USE MOD_ENERGY,ONLY:NEPOL,NEFD1,nefd2,nefd3,necontn,lactive_contn,
     &                    lepath_contn,emin_contn,emax_contn,ime_contn,
     &                    nepol_contn
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='EPATH')
C
C Dummy arguments
C
      REAL*8 EILOW,EIMAG,EMAX,EMIN
      INTEGER IGRID,IPRINT,NE,NEMAX
      COMPLEX*16 ETAB(NEMAX),WETAB(NEMAX)
C
C Local variables
C
      REAL*8 CNST,DE1,DE2,DX,EBAR,EISTEP,EMID,EPS,ERSTEP,ESTEP,PHI(NE),
     &       PHI1,PHI2,R,W(:),X,X0,Y,Y1,Y2,YY(NE),Z(:),Z0
      COMPLEX*16 CSUM
      INTEGER I,IA_ERR,J,NE1,NE10,NE2,NE20,NPARA,NPERP,NPERP1,NPERP2
      integer nevertical
      LOGICAL SETEPS
      real*8 kb
      DATA KB /0.6333659D-5/
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE W,Z
C
      ALLOCATE (W(NEMAX),Z(NEMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: Z')
C
      ESTEP = 0.0D0
      CALL CINIT(NEMAX,WETAB)
C
C---------------------------------------------------------------------
C IGRID = 0                  gauss-legendre quadrature along real axis
C
      IF ( IGRID.EQ.0 ) THEN
C
         IF ( NE.GT.1 ) THEN
C
            CALL GAULEG(EMIN,EMAX,Z,W,NE)
C
            DO I = 1,NE
               ETAB(I) = Z(I)
               WETAB(I) = W(I)
            END DO
C
         ELSE
C
            ETAB(1) = EMIN
            WETAB(1) = 1D0
C
         END IF
C
      ELSE IF ( IGRID.EQ.1 ) THEN
C---------------------------------------------------------------------
C IGRID= 1                      only real energies -- equidistant path
C
         IF ( NE.GT.1 ) ESTEP = (EMAX-EMIN)/DBLE(NE-1)
C
         ETAB(1) = EMIN
         DO I = 2,NE
            ETAB(I) = ETAB(I-1) + ESTEP
         END DO
C
      ELSE IF ( IGRID.EQ.2 ) THEN
C
C---------------------------------------------------------------------
C IGRID= 2                                    rectangular complex path
C
         ESTEP = (EMAX-EMIN+2*EIMAG)/DBLE(NE-1)
         NPERP = INT(EIMAG/ESTEP) + 2
         NPARA = NE - 2*NPERP
         EISTEP = EIMAG/DBLE(NPERP-1)
         ERSTEP = (EMAX-EMIN)/DBLE(NPARA+1)
         ETAB(1) = EMIN
         DO I = 2,NPERP
            ETAB(I) = ETAB(I-1) + EISTEP*CI
         END DO
         DO I = (NPERP+1),(NPERP+NPARA+1)
            ETAB(I) = ETAB(I-1) + ERSTEP
         END DO
         DO I = (NE-NPERP+2),NE
            ETAB(I) = ETAB(I-1) - EISTEP*CI
C
         END DO
C
      ELSE IF ( IGRID.EQ.3 ) THEN
C---------------------------------------------------------------------
C IGRID= 3                 straight complex path parallel to real axis
C
         IF ( NE.GT.1 ) ESTEP = (EMAX-EMIN)/DBLE(NE-1)
         ETAB(1) = EMIN + EIMAG*CI
C
         DO I = 2,NE
            ETAB(I) = ETAB(I-1) + ESTEP
         END DO
C                                    weights for trapez rule
         WETAB(1) = 0.5D0*ESTEP
         WETAB(2:(NE-1)) = ESTEP
         WETAB(NE) = 0.5D0*ESTEP
C
      ELSE IF ( IGRID.EQ.4 ) THEN
C---------------------------------------------------------------------
C IGRID= 4               rectangular complex grid, return to real axis
C
         EILOW = 0.0005D0
         EIMAG = MAX(EIMAG,0.0005D0)
         NPERP1 = INT(NE*0.25)
         NPERP2 = INT(NE*0.5+1)
         NPARA = MAX(1,NE-1-NPERP1-NPERP2)
         NE = 1 + NPERP1 + NPARA + NPERP2
         EISTEP = (EIMAG-EILOW)/DBLE(NPERP1)
         ERSTEP = (EMAX-EMIN)/DBLE(NPARA)
         DX = LOG(EIMAG/EILOW)/DBLE(NPERP2)
         X0 = LOG(EIMAG)
C
         ETAB(1) = EMIN + EILOW*CI
         DO I = 2,(1+NPERP1)
            ETAB(I) = ETAB(I-1) + EISTEP*CI
         END DO
C
         DO I = (1+NPERP1+1),(1+NPERP1+NPARA)
            ETAB(I) = ETAB(I-1) + ERSTEP
         END DO
C
         X = X0
         DO I = (1+NPERP1+NPARA+1),NE
            X = X - DX
            ETAB(I) = DCMPLX(EMAX,EXP(X))
         END DO
C
      ELSE IF ( IGRID.EQ.5 ) THEN
C---------------------------------------------------------------------
C IGRID = 5                complex arc using gauss-legendre quadrature
C
         CALL GAULEG(-1.0D0,1.0D0,Z,W,NE)
C
         CNST = 0.25D0*PI*(EMAX-EMIN)
         EBAR = 0.50D0*(EMAX+EMIN)
C
!         EIMAG = MAX(EIMAG,0.0005D0) !REC shift imaginary part by EIMAG, but let it be zero
         DO I = 1,NE
            Y = 0.5D0*PI*Z(I)
            ETAB(I) = EBAR + 0.5D0*(EMAX-EMIN)
     &                *CI*(CDEXP(DCMPLX(0.0D0,-Y))+EIMAG)
            WETAB(I) = CNST*W(I)*CDEXP(DCMPLX(0.0D0,-Y))
         END DO
      ELSE IF ( IGRID.EQ.6 ) THEN
C
C---------------------------------------------------------------------
C IGRID = 6                         standard mesh for X-ray absorption
C
         IF ( EMAX.LT.5D0 ) THEN
            NE10 = 50
            NE20 = 120
         ELSE
            NE10 = 30
            NE20 = 170
         END IF
C
         NE1 = MAX(2,(NE*NE10)/(NE20+NE10))
         NE2 = NE - NE1
C
         EMID = 1.1D0
C
         DE2 = (EMAX-EMID)/DBLE(NE2)
         DE1 = (EMID-EMIN)/DBLE((NE1-1)**2)
C
         ETAB(1) = DCMPLX(EMIN,EIMAG)
         DO I = 2,NE1
            ETAB(I) = ETAB(1) + DE1*(I-1)**2
         END DO
C
         DO I = NE1 + 1,NE
            ETAB(I) = ETAB(I-1) + DE2
         END DO
C
      ELSE IF ( IGRID.EQ.7 ) THEN
C---------------------------------------------------------------------
C IGRID = 7  mesh for X-ray absorption  and  gauss-legendre quadrature
C
         CALL GAULEG(EMIN,EMAX,Z,W,NE)
C
         DO I = 1,NE
            ETAB(I) = DCMPLX(Z(I),EIMAG)
            WETAB(I) = W(I)
         END DO
C
      ELSE IF ( IGRID.EQ.8 ) THEN
C---------------------------------------------------------------------
C IGRID = 8                 logarithmic E-mesh for contour integration
C
         CALL GAULEG(-1.0D0,1.0D0,Z,W,NE)
C
         PHI1 = PI
         PHI2 = 0.0D0
C        EPS = 0.025D0 ! 0.050D0: Safer value, quick fix by OS
         EPS = 0.050D0
C
         CALL EGRID8_IME(EPS,EMIN,EMAX,NE,SETEPS)
         R = (EMAX-EMIN)/2.0D0
         Z0 = EMIN + R
         Y1 = -LOG((EPS+PHI1)/EPS)
         Y2 = -LOG((EPS+PHI2)/EPS)
C
         DO J = 1,NE
            YY(J) = 0.5D0*(Y2-Y1)*Z(J) + 0.5D0*(Y2+Y1)
            PHI(J) = EPS*(EXP(-YY(J))-1.0D0)
            ETAB(J) = Z0 + R*EXP(CI*PHI(J))
            WETAB(J) = -0.5D0*CI*EPS*EXP(-YY(J))*(Y2-Y1)*(ETAB(J)-Z0)
     &                 *W(J)
         END DO
C
         IF ( SETEPS ) WRITE (*,99004) 'Energy path 8:    EPS =',EPS,
     &                                 '      Im(ETAB(NE)) =',
     &                                 DIMAG(ETAB(NE))
C
      ELSE IF ( IGRID.EQ.9 ) THEN
C---------------------------------------------------------------------
C IGRID = 9                  H. Akai's  E-mesh for contour integration
C
         CALL EPATH_AKAI(EMIN,EMAX,EIMAG,NE,ETAB,WETAB,NEMAX)
C
      ELSE IF ( IGRID.EQ.10 ) THEN
C---------------------------------------------------------------------
C IGRID =10                 DK, ellipse E-mesh for contour integration
C
         CALL EPATH_ELLIPSE(EMIN,EMAX,EIMAG,0.5D0,NE,ETAB,WETAB,NEMAX)
C
C---------------------------------------------------------------------
C
      ELSE IF ( IGRID.EQ.11 ) THEN
C---------------------------------------------------------------------
C IGRID =11                 Fermi-Dirac with Matsubara poles
C---------------------------------------------------------------------
c nepole - number of poles for scf
c nefd1 - number of points in [emin, emin + ci*height]
c nefd2 - number of points in [emin + ci*height, emax - 30*kB*T + ci*height]
c nefd3 - number of points in [emax - 30*kB*T + ci*height, emax + 30*kB*T + ci*height]
c
         npara = nefd1 + nefd2 + nefd3
         CALL EPATH_TEMP(ETAB,WETAB,npara+nepol,EMIN,EMAX,EMAX,EIMAG,
     &        NEPOL,NEFD1,nefd2,nefd3,npara+nepol)
c
c modified by XJQ: add additional energies with zeor weights for 
c                  analytic continuation using rectangular complex grids from emin_contn to emax_contn
c                  when lactive_contn=.true.,
c                  energy points from ne-nepole-necontn+1 to ne will be used for continuation;
c
         if(necontn>0) then
           do i=1,nepol
             etab(ne-i+1)=etab(npara+nepol-i+1)
             wetab(ne-i+1)=wetab(npara+nepol-i+1)
           enddo
           wetab(npara+1:npara+necontn) = (0d0,0d0)
c
           if(ime_contn<1e-8) ime_contn = aimag(etab(npara))
           if(emax_contn .le. emax) 
     &       emax_contn = real(etab(npara)) + ci*ime_contn
c           nevertical = max( int( ( ime_contn / aimag(etab(npara)) ) 
c     &                            * nefd1 + 0.25 ) + 1, 2 )
c           if(necontn-1-2*nevertical .le. 0) stop 'too few necontn'
c
c           eistep = ime_contn / ( 2d0 * nevertical + 1d0 )
c           etab(npara+1) = emin_contn + ci*eistep
c           do i=npara+2,npara+nevertical
c             etab(i) = etab(i-1) + ci*eistep*2d0
c           enddo
c           etab(npara+nevertical+1) = emin_contn + ci*ime_contn
           do i=1,nepol_contn
             etab(npara+i) = emax + ci*(2*i-1)*pi*kb*eimag
           enddo
           if(nepol_contn < necontn) then
             etab(npara+nepol_contn+1) = emin_contn + ci*ime_contn
c
c           erstep = (emax_contn - emin_contn) / (necontn-1-2*nevertical)
             erstep = (emax_contn - emin_contn) / 
     &                (necontn-1-nepol_contn)
c           do i=npara+nevertical+2,npara+necontn-nevertical
             do i=npara+nepol_contn+2,npara+necontn
               etab(i) = etab(i-1) + erstep
             enddo
c           etab(npara+necontn) = emax_contn + ci*eistep
c           do i=npara+necontn-1,npara+necontn-nevertical+1,-1
c             etab(i) = etab(i+1) + ci*eistep*2d0
c           enddo
c           if(.not. lactive_contn) then
c             etab(npara+1:npara+nepol) = etab(ne-nepol+1:ne)
c           endif
           endif
         endif
c end-mod-xjq
c
      END IF
C
      DEALLOCATE (W,Z,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC')
C
      IF ( IPRINT.LT.2 .AND. IGRID.NE.9 ) RETURN
C
      WRITE (6,99001) IGRID
      CSUM = 0D0
      DO I = 1,NE
         CSUM = CSUM + WETAB(I)
         WRITE (6,99002) I,NE,ETAB(I),EMIN,EIMAG,EMAX,WETAB(I)
      END DO
      WRITE (6,99003) CSUM
C
99001 FORMAT (//,2X,'ENERGY INTEGRATION PATH ',I2,' (in Ry)',//,
     &        '    I   NE        ETAB',15X,'EMIN      EIMAG     EMAX'/)
99002 FORMAT (2I5,2F16.8,3X,3F10.5,2F16.8)
99003 FORMAT (/,10X,'sum of weights: ',2F16.8,/)
99004 FORMAT (/,a,f10.5,a,f12.7)
      END
C*==epath_ellipse.f    processed by SPAG 6.70Rc at 19:23 on  1 Nov 2014
      SUBROUTINE EPATH_ELLIPSE(EMIN,EMAX,EIMAG,H,NE,ETAB,WETAB,NEMAX)
C   ********************************************************************
C   *                                                                  *
C   *  ROUTINE TO CREATE  E - PATH, ELLIPSE                            *
C   *                                                                  *
C   *  EMIN = lowest   real(E)                                         *
C   *  EMAX = highest  real(E)                                         *
C   *                                                                  *
C   *  H = B/A                                                         *
C   *        + semi-major axis A and semi-minor axis B                 *
C   *        + if H=0 then it is set: B = EIMAG                        *
C   *        + Note: H=1 gives semicircle                              *
C   *                                                                  *
C   *  DK 14-10-21                                                     *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:PI,CI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='EPATH_ELLIPSE')
C
C Dummy arguments
C
      REAL*8 EIMAG,EMAX,EMIN,H
      INTEGER NE,NEMAX
      COMPLEX*16 ETAB(NEMAX),WETAB(NEMAX)
C
C Local variables
C
      REAL*8 A,B,DRDPHI,PHI,PI2,R,W(:),X,Z(:)
      COMPLEX*16 DGDX,Z0
      INTEGER IA_ERR,IE
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE W,Z
C
      ALLOCATE (W(NEMAX),Z(NEMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: Z')
C
      PI2 = PI/2D0
      Z0 = (EMAX+EMIN)/2D0
C
      A = (EMAX-EMIN)/2D0
      IF ( ABS(H).LT.1D-8 ) THEN
         B = EIMAG
      ELSE
         B = A*H
      END IF
C
C      a= 1d0
C      b= 1d0
C
      CALL GAULEG(-1.0D0,1.0D0,Z,W,NE)
C
      DO IE = 1,NE
         X = Z(IE)
         PHI = PI2*(1-X)
         CALL ELLIPSER(A,B,PHI,R,DRDPHI)
         ETAB(IE) = Z0 + R*CDEXP(CI*PI2*(1-X))
         DGDX = -PI2*CDEXP(CI*PI2*(1-X))*DRDPHI + R*PI2*CDEXP(-CI*PI2*X)
         WETAB(IE) = W(IE)*DGDX
      END DO
C
      X = 1D0 - DREAL(SUM(WETAB))/(2*A)
      IF ( ABS(X).GT.1D-6 ) THEN
         WRITE (6,'(10X,"Rel error in weights",E25.12)') X
         CALL STOP_MESSAGE(ROUTINE,'Error in weights')
      END IF
      END
C*==ellipser.f    processed by SPAG 6.70Rc at 19:23 on  1 Nov 2014
      SUBROUTINE ELLIPSER(A,B,PHI,R,DRDPHI)
C   ********************************************************************
C   *  Ellipse: Polar form relative to center                          *
C   *                                                                  *
C   *  Input:                                                          *
C   *        A : semi-major axis  a                                    *
C   *        B : semi-minor axis  b                                    *
C   *      PHI : angle phi                                             *
C   *                                                                  *
C   *  Output:                                                         *
C   *        R : r(phi)                                                *
C   *   DRDPHI : dr/dphi                                               *
C   *                                                                  *
C   *  DK 14-10-21                                                     *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 A,B,DRDPHI,PHI,R
C
C Local variables
C
      REAL*8 D
C
C*** End of declarations rewritten by SPAG
C
      R = A*B/DSQRT(B**2*(COS(PHI)**2)+A**2*(SIN(PHI)**2))
C
      D = ((B*COS(PHI))**2+(A*SIN(PHI))**2)**(3D0/2D0)
      DRDPHI = A*B*(B+A)*(B-A)*COS(PHI)*SIN(PHI)/D
C
      END
C*==epath_akai.f    processed by SPAG 6.70Rc at 19:23 on  1 Nov 2014
      SUBROUTINE EPATH_AKAI(EMIN,EMAX,EIMAG,NE,ETAB,WETAB,NEMAX)
C   ********************************************************************
C   *                                                                  *
C   *  ROUTINE TO CREATE  E - PATH according to H. Akai's scheme       *
C   *                                                                  *
C   *                                                                  *
C   *  EIMAG= distance from real axis                                  *
C   *  EMIN = lowest   real(E)                                         *
C   *  EMAX = highest  real(E)                                         *
C   *                                                                  *
C   *  HE 00/03/17                                                     *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0,C1,PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='EPATH_AKAI')
      INTEGER NGMAX
      PARAMETER (NGMAX=30)
C
C Dummy arguments
C
      REAL*8 EIMAG,EMAX,EMIN
      INTEGER NE,NEMAX
      COMPLEX*16 ETAB(NEMAX),WETAB(NEMAX)
C
C Local variables
C
      COMPLEX*16 A,B,BNDE,C,DETL(:),E(:),T(:),T0(:),T1(:),T2(:),
     &           WT(:,:,:)
      REAL*8 BETA,EDELT,EF,EW,EWIDTH,EZ,F,H,ONE,R,RC,THETA
      INTEGER I,IE,K,KC,KMX,NG
C
C*** End of declarations rewritten by SPAG
C
C     data h/5d-1/,rc/4d-1/,one/0.99999999d0/
      DATA H/2D-1/,RC/4D-1/,ONE/0.99999999D0/
C
      ALLOCATABLE E,WT,T,T0,T1,T2,DETL
C
C=======================================================================
C     subroutine cemesh(ef,ewidth,edelt,ebtm,e,kmx)
C=======================================================================
C
      EF = EMAX
      EWIDTH = EMAX - EMIN
      EDELT = 1D-3
      KMX = NE
C
      ALLOCATE (E(KMX))
C
C--------------------------------------------------------------------
C     Generate a semielliptic energy contour. Mesh points are located
C     following fermi's distribution function such that they are
C     distributed densely near the real axis.
C     The following are the examples of typical cases used in the past.
C       edelt/3d-4/, r/0.65d0 /, h/2d-1/, rc/4d-1/
C       edelt/1d-3/, r/0.70d0 /, h/2d-1/, rc/4d-1/
C       edelt/1d-3/, r/1.20d0 /, h/2d-1/, rc/4d-1/
C       edelt/3d-3/, r/1.40d0 /, h/2d-1/, rc/4d-1/
C       edelt/1d-4/, r/0.40d0 /, h/5d-1/, rc/4d-1/
C     coded by H.Akai, 1983, Juelich
C     latest version, 30 Nov.1997, Osaka
C--------------------------------------------------------------------
C
      R = EWIDTH/2D0
C
      IF ( EIMAG.GT.0.001D0 ) H = EIMAG/R
C
C     pi=4d0*atan(1d0)
      KC = INT(RC*DBLE(KMX)+5D-1)
      BETA = LOG(PI/(EDELT/H/R)-1D0)/DBLE(KMX-KC)
      F = PI*(EXP(BETA*DBLE(1-KC))+1D0)*ONE
      DO K = 1,KMX
         THETA = F/(EXP(BETA*DBLE(K-KC))+1D0)
         E(K) = EF + R*DCMPLX(COS(THETA)-1D0,H*SIN(THETA))
      END DO
C     EBTM = DBLE(E(1))
C
C=======================================================================
C     subroutine cgnwt(e,wt,ng,kmx,ew,ez)
C=======================================================================
C
      NG = 1
      EW = (EMAX+EMIN)/2D0
      EZ = (EMAX-EMIN)/2D0
C
      ALLOCATE (WT(NG,3,KMX),T(NGMAX),T0(NGMAX),T1(NGMAX),T2(NGMAX))
C
C-----------------------------------------------------------------------
C     Gives weighting function which is used for the contour
C     integration along complex energy.
C     coded by H.Akai, 1983, Juelich
C-----------------------------------------------------------------------
      IF ( NG.GT.NGMAX ) THEN
         WRITE (6,'(//,10X,''NG='',I3,3X,''NGMAX='',I3)') NG,NGMAX
         CALL STOP_MESSAGE(ROUTINE,'NG > NGMAX')
      END IF
C
      DO K = 1,KMX
         T(1) = 1D0
         T(2) = (E(K)-EW)/EZ
         T0(1) = T(2)
         T0(2) = T(2)**2/2D0
         T1(1) = T0(2)
         T1(2) = T(2)**3/3D0
         DO I = 3,NG
            T(I) = 2D0*T(2)*T(I-1) - T(I-2)
            T0(I) = (DBLE(I-4)*T0(I-2)-(2D0-4D0*T0(2))*T(I-1))/DBLE(I)
            T1(I) = (DBLE(I-1)*T0(I-1)-(1D0-2D0*T0(2))*T(I))/DBLE(I+1)
         END DO
         T2(1) = T1(2)
         DO I = 2,NG
            T2(I) = (DBLE(I-1)*T1(I-1)+T0(I)-T(2)*(1D0-2D0*T0(2))*T(I))
     &              /DBLE(I+2)
         END DO
         DO I = 1,NG
            WT(I,1,K) = EZ*T0(I)
            WT(I,2,K) = EZ**2*T1(I) + EW*WT(I,1,K)
            WT(I,3,K) = EZ**3*T2(I) + 2D0*EW*EZ**2*T1(I)
     &                  + EW**2*WT(I,1,K)
         END DO
      END DO
C
C=======================================================================
C     generate the corresponding weights for E-integration
C=======================================================================
C     subroutine banden(detl,e,wt,bnde,kmx,ng)
C=======================================================================
C
      ALLOCATE (DETL(KMX))
C
C-----------------------------------------------------------------------
C     Calculate the band energy by integrating the phase.
C     coded by H.Akai, 1985, Juelich
C-----------------------------------------------------------------------
C
      DO IE = 1,KMX
C
         DETL(:) = C0
         DETL(IE) = C1
C
         BNDE = 0D0
C
         DO K = 1,KMX - 2,2
            A = ((DETL(K+1)-DETL(K))/(E(K+1)-E(K))-(DETL(K+2)-DETL(K+1))
     &          /(E(K+2)-E(K+1)))/(E(K)-E(K+2))
            B = (DETL(K+1)-DETL(K))/(E(K+1)-E(K)) - A*(E(K+1)+E(K))
            C = DETL(K) - (A*E(K)+B)*E(K)
            BNDE = BNDE + A*(WT(1,3,K+2)-WT(1,3,K))
     &             + B*(WT(1,2,K+2)-WT(1,2,K))
     &             + C*(WT(1,1,K+2)-WT(1,1,K))
         END DO
C
         ETAB(IE) = E(IE)
         WETAB(IE) = BNDE
C
      END DO
C
      END
C 20.09.95 *************************************************************
      SUBROUTINE EPATH_TEMP(EZ,DF,NPNT,EBOT,EMU,EFERMI,TK,
     &                  NPOL,NPNT1,NPNT2,NPNT3,IEMXD)
C **********************************************************************
C *                                                                    *
C * This subroutine provides the energy mesh in array EZ and the       *
C * appropriate integration weights in array DF.                       *
C *                                                                    *
C * Poles of the Fermi function C (Matsubara frequencies) and          *
C * a contour in the complex energy are used as described in (????).   *
C *                                                                    *
C * The contour consists of three straight lines with                  *
C * NPNT1, NPNT2, and NPNT3 integration points and is determined by    *
C * the input arguments: EBOT, EMU, TK, and NPOL.                      *
C *                                                                    *
C *            TK   = temperature in K                                 *
C *            EMU  = chemical potential in Ry                         *
C *            EBOT = bottom of contour in Ry                          *
C *            NPOL = number of Matsubara frequencies                  *
C *                                                                    *
C * The three lines are defined by:                                    *
C *                                                                    *
C *  1. the line from EBOT to EBOT+2*NPOL*pi*i*k*TK                    *
C *              with NPNT1 integration points (Gauss-Legendre rule)   *
C *                                                                    *
C *  2. the line from EBOT+2*NPOL*pi*i*k*TK to                         *
C *                   EMU+(2*NPOL*pi*i-30)*k*TK                        *
C *              with NPNT2 integration points (Gauss-Legendre rule)   *
C *                                                                    *
C *  3. the line from EMU+(2*NPOL*pi*i-30)*k*TK to infinity            *
C *              with NPNT3 integration points (Gauss-Fermi-Dirac rule)*
C *                                                                    *
C *  The total number of integration points is given by:               *
C *              NPNT=NPNT1+NPNT2+NPNT3+NPOL                           *
C *                                                                    *
C *  The integration points and weights on three lines are chosen      *
C *  according to Gauss integration rules. Only in third interval      *
C *  the Fermi function matters since exp(x) < 10**(-10) for x < -25.  *
C *                                                                    *
C *  There are two special cases determined by NPOL = 0 and NPOL < 0.  *
C *                                                                    *
C *  a) NPOL = 0 leads to density-of-states calculations               *
C *  with constant integration weights and equally distributed points  *
C *  between EBOT - pi*i*k*TK and EMU - pi*i*k*TK.                     *
C *                                                                    *
C *  The total number of integration points is given by:               *
C *              NPNT=NPNT2                                            *
C *                                                                    *
C *  b) NPOL < 0 is meant for calculations where the Fermi-Dirac       *
C *  function is replaced by a step function with step at EMU. When    *
C *  this option is used no poles of the Fermi-Dirac function are used *
C *  and the contour consists of the three straight lines:             *
C *                                                                    *
C *  1. the line from EBOT to EBOT-2*NPOL*pi*i*k*TK                    *
C *              with NPNT1 integration points (Gauss-Legendre rule)   *
C *                                                                    *
C *  2. the line from EBOT-2*NPOL*pi*i*k*TK to EMU-2*NPOL*pi*i*k*TK    *
C *              with NPNT2 integration points (Gauss-Legendre rule)   *
C *                                                                    *
C *  3. the line from EMU-2*NPOL*pi*i*k*TK to EMU                      *
C *              with NPNT3 integration points (Gauss-Legendre rule)   *
C *                                                                    *
C *  The total number of integration points is given by:               *
C *              NPNT=NPNT1+NPNT2+NPNT3                                *
C *                                                                    *
C **********************************************************************
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION EBOT,EMU,TK,EFERMI
      INTEGER NPNT,NPNT1,NPNT2,NPNT3,NPOL,IEMXD
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX DF(*),EZ(*)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX DE
      DOUBLE PRECISION ER,ETK,KB,PI,RYD
      INTEGER I
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION WI(128),XI(128)
C     ..
C     .. External Subroutines ..
      LOGICAL OPT
      EXTERNAL GAUFD,GAULEG,OPT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DCMPLX
C     ..
C     .. Data statements ..
      DATA PI /3.14159265358979312D0/
      DATA KB /0.6333659D-5/
      DATA RYD/13.6058D0/
C     ..
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (6,'(5X,A,F10.6," (Ry)",8X,A,F10.6," (Ry)")')
     &     'E min = ',EBOT,'Fermi energy = ',EFERMI
      WRITE (6,'(5X,A,F10.6," (Ry)",8X,A,F15.6," (K )",/,5X,62(1H-))')
     &     'E max = ',EMU,'Temperature  = ',TK
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
      ETK = PI*KB*TK
C ======================================================================
      IF (NPOL.EQ.0) THEN
         DE = (EMU-EBOT)
         IF (NPNT2.GT.1) THEN
            DE = DE/(NPNT2-1)
         ELSE
            DE=DCMPLX(1.0D0,0.0D0)
         END IF
         NPNT = 0
         DO 10 I = 1,NPNT2
            NPNT = NPNT + 1
            IF ( NPNT.GT.IEMXD ) THEN 
               WRITE(6,'(/,5X,2A,I4)') 
     &              'Dimension ERROR: Increase IEMXD in inc.p to ',
     &              'at least ',NPNT
               STOP '     < EMESHT1 >'
            END IF
            ER = EBOT + (I-1)*DE
            EZ(NPNT) = DCMPLX(ER,ETK)
            DF(NPNT) = DE
 10      CONTINUE
         WRITE (6,FMT=9000) NPNT,ETK,ETK*RYD
C ------------------------------------------------------------- NPOL > 0
      ELSE IF (NPOL.GT.0) THEN
         CALL GAULEG(-1.0D0,1.0D0,XI,WI,NPNT1)
         DE = NPOL*DCMPLX(0.0D0,ETK)
         NPNT = 0
         DO 20 I = 1,NPNT1
            NPNT = NPNT + 1
            IF ( NPNT.GT.IEMXD ) THEN 
               WRITE(6,'(/,5X,2A,I4)') 
     &              'Dimension ERROR: Increase IEMXD in inc.p to ',
     &              'at least ',NPNT
               STOP '     < EMESHT2 >'
            END IF
            EZ(NPNT) = XI(I)*DE + DE + EBOT
            DF(NPNT) = WI(I)*DE
 20      CONTINUE
         CALL GAULEG(-1.0D0,1.0D0,XI,WI,NPNT2)
         DE = (EMU-30*KB*TK-EBOT)*0.5D0
         DO 30 I = 1,NPNT2
            NPNT = NPNT + 1
            IF ( NPNT.GT.IEMXD ) THEN 
               WRITE(6,'(/,5X,2A,I4)') 
     &              'Dimension ERROR: Increase IEMXD in inc.p to ',
     &              'at least ',NPNT
               STOP '     < EMESHT3 >'
            END IF
            EZ(NPNT) = XI(I)*DE + DE + EBOT + 2*NPOL*DCMPLX(0.0D0,ETK)
            DF(NPNT) = WI(I)*DE
 30      CONTINUE
         CALL GAUFD(XI,WI,NPNT3)
         DE = 30*KB*TK
         DO 40 I = 1,NPNT3
            NPNT = NPNT + 1
            IF ( NPNT.GT.IEMXD ) THEN 
               WRITE(6,'(/,5X,2A,I4)') 
     &              'Dimension ERROR: Increase IEMXD in inc.p to ',
     &              'at least ',NPNT
               STOP '     < EMESHT4 >'
            END IF
            EZ(NPNT) = XI(I)*DE + EMU + 2*NPOL*DCMPLX(0.0D0,ETK)
            DF(NPNT) = WI(I)*DE
 40      CONTINUE
         DO 50 I = NPOL,1,-1
            NPNT = NPNT + 1
            IF ( NPNT.GT.IEMXD ) THEN 
               WRITE(6,'(/,5X,2A,I4)') 
     &              'Dimension ERROR: Increase IEMXD in inc.p to ',
     &              'at least ',NPNT
               STOP '     < EMESHT5 >'
            END IF
            EZ(NPNT) = EMU + (2*I-1)*DCMPLX(0.0D0,ETK)
            DF(NPNT) = -2*DCMPLX(0.0D0,ETK)
 50      CONTINUE
         WRITE(6,9090) NPNT,NPOL,NPNT1,NPNT2,NPNT3
C ------------------------------------------------------------- NPOL < 0
      ELSE
         IF (NPNT1.GT.0) CALL GAULEG(-1.0D0,1.0D0,XI,WI,NPNT1)
         DE = -NPOL*DCMPLX(0.0D0,ETK)
         NPNT = 0
         DO 60 I = 1,NPNT1
            IF ( NPNT.GT.IEMXD ) THEN 
               WRITE(6,'(/,5X,2A,I4)') 
     &              'Dimension ERROR: Increase IEMXD in inc.p to ',
     &              'at least ',NPNT
               STOP '     < EMESHT6 >'
            END IF
            NPNT = NPNT + 1
            EZ(NPNT) = XI(I)*DE + DE + EBOT
            DF(NPNT) = WI(I)*DE
 60      CONTINUE
         CALL GAULEG(-1.0D0,1.0D0,XI,WI,NPNT2)
         DE = (EMU-EBOT)*0.5D0
         DO 70 I = 1,NPNT2
            NPNT = NPNT + 1
            IF ( NPNT.GT.IEMXD ) THEN 
               WRITE(6,'(/,5X,2A,I4)') 
     &              'Dimension ERROR: Increase IEMXD in inc.p to ',
     &              'at least ',NPNT
               STOP '     < EMESHT7 >'
            END IF
            EZ(NPNT) = XI(I)*DE + DE + EBOT - 2*NPOL*DCMPLX(0.0D0,ETK)
            IF (OPT('GF-EF   ')) EZ(NPNT) = EMU + NPOL*DCMPLX(0.0D0,ETK)
            DF(NPNT) = WI(I)*DE
 70      CONTINUE
         IF (NPNT3.GT.0) CALL GAULEG(-1.0D0,1.0D0,XI,WI,NPNT3)
         DE = -NPOL*DCMPLX(0.0D0,ETK)
         DO 80 I = NPNT3,1,-1
            NPNT = NPNT + 1
            IF ( NPNT.GT.IEMXD ) THEN 
               WRITE(6,'(/,5X,2A,I4)') 
     &              'Dimension ERROR: Increase IEMXD in inc.p to ',
     &              'at least ',NPNT
               STOP '     < EMESHT8 >'
            END IF
            EZ(NPNT) = XI(I)*DE + DE + EMU
            DF(NPNT) = -WI(I)*DE
 80      CONTINUE
         WRITE(6,9091) NPNT,-NPOL,NPNT1,NPNT2,NPNT3
      END IF
C ======================================================================
      WRITE(6,*)
      RETURN
 9000 FORMAT (5X,'Density-of-States calculation',/,
     &        5X,'Number of energy points :',I4,4X,'broadening =',
     &        3P,F9.3,' ( mRy )',/,48X,' =',3P,F9.3,' ( meV )')
 9090 FORMAT (5X,'GF integration rectangular contour ( ImE > 0 )',/,
     &     5X,'Number of energy points :',I4,13X,'poles =',I2,/,
     &     23X,'contour: N1 =',I2,', N2 =',I4,', N3 =',I2)
 9091 FORMAT (5X,'GF integration rectangular contour ( ImE < 0 )',/,
     &     5X,'Number of energy points :',I4,13X,'poles =',I2,/,
     &     23X,'contour: N1 =',I2,', N2 =',I4,', N3 =',I2)
      END
c     ************************************************
      SUBROUTINE GAUFD(XI,WI,N)
c     ************************************************
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION WI(*),XI(*)
C     ..
      IF (N.EQ.1) THEN
        XI(1) = -49817229548128141768.D-20
        WI(1) = 10000000000000031192.D-19
        GO TO 10

      END IF

      IF (N.EQ.2) THEN
        XI(1) = -78465071850839016234.D-20
        XI(2) = -20091536266094051757.D-20
        WI(1) = 50923235990870048433.D-20
        WI(2) = 49076764009130263488.D-20
        GO TO 10

      END IF

      IF (N.EQ.3) THEN
        XI(1) = -88288518955458358024.D-20
        XI(2) = -48117621892777473749.D-20
        XI(3) = -88198184413497647625.D-21
        WI(1) = 28858444436509900908.D-20
        WI(2) = 45966895698954759346.D-20
        WI(3) = 25174659864535651667.D-20
        GO TO 10

      END IF

      IF (N.EQ.4) THEN
        XI(1) = -92613063531202843773.D-20
        XI(2) = -64918327008663578157.D-20
        XI(3) = -28982568853420020298.D-20
        XI(4) = -24595209663255169680.D-21
        WI(1) = 18501429405165520392.D-20
        WI(2) = 34614391006511784214.D-20
        WI(3) = 34152482191988153127.D-20
        WI(4) = 12731697396334854188.D-20
        GO TO 10

      END IF

      IF (N.EQ.5) THEN
        XI(1) = -94875333872503463082.D-20
        XI(2) = -74805843506753178608.D-20
        XI(3) = -45504655263391074765.D-20
        XI(4) = -16657582360358973599.D-20
        XI(5) = 27402283545708211900.D-21
        WI(1) = 12939804504572789754.D-20
        WI(2) = 26102400189213290231.D-20
        WI(3) = 30851911091450589451.D-20
        WI(4) = 24746815229701880449.D-20
        WI(5) = 53590689850617620359.D-21
        GO TO 10

      END IF

      IF (N.EQ.6) THEN
        XI(1) = -96204950250095729781.D-20
        XI(2) = -80971428101130972258.D-20
        XI(3) = -57293627456482418171.D-20
        XI(4) = -30755197635518367504.D-20
        XI(5) = -82123839469384988331.D-21
        XI(6) = 83748358371240941581.D-21
        WI(1) = 96268650841705383829.D-21
        WI(2) = 20246201047059595265.D-20
        WI(3) = 26160719441051813381.D-20
        WI(4) = 25781980698475975536.D-20
        WI(5) = 16683001513553609336.D-20
        WI(6) = 15012322156887800205.D-21
        GO TO 10

      END IF

      IF (N.EQ.7) THEN
        XI(1) = -97053934379083423143.D-20
        XI(2) = -85045695849615413757.D-20
        XI(3) = -65665104053460540522.D-20
        XI(4) = -42357896269371657364.D-20
        XI(5) = -19472732441816555564.D-20
        XI(6) = -19669621223691542539.D-21
        XI(7) = 15142830586888806919.D-20
        WI(1) = 74948008822570509041.D-21
        WI(2) = 16170863905729061704.D-20
        WI(3) = 22007120289205973485.D-20
        WI(4) = 23880411919774885893.D-20
        WI(5) = 20952460047488907594.D-20
        WI(6) = 92465405554445737538.D-21
        WI(7) = 24780240009985858690.D-22
        GO TO 10

      END IF

      IF (N.EQ.8) THEN
        XI(1) = -97630544447925725992.D-20
        XI(2) = -87873822716479965943.D-20
        XI(3) = -71736329217593360204.D-20
        XI(4) = -51463306578144813387.D-20
        XI(5) = -29967081434747298359.D-20
        XI(6) = -10763455942936048359.D-20
        XI(7) = 35963113675701677498.D-21
        XI(8) = 23003149140664609750.D-20
        WI(1) = 60394634019629989770.D-21
        WI(2) = 13252509350880929004.D-20
        WI(3) = 18643612522057003210.D-20
        WI(4) = 21413715867867937533.D-20
        WI(5) = 21005092708864293339.D-20
        WI(6) = 16003068683842947897.D-20
        WI(7) = 36159126989806650464.D-21
        WI(8) = 26624765543536915040.D-23
        GO TO 10

      END IF

      IF (N.EQ.9) THEN
        XI(1) = -98041275487012188695.D-20
        XI(2) = -89918326179154863440.D-20
        XI(3) = -76254129548477842110.D-20
        XI(4) = -58579104527384144901.D-20
        XI(5) = -38924212142470946276.D-20
        XI(6) = -19724340764961096691.D-20
        XI(7) = -40039281758884590381.D-21
        XI(8) = 97228170103579374416.D-21
        XI(9) = 31678885353558278864.D-20
        WI(1) = 49992516372028853833.D-21
        WI(2) = 11099301824870447793.D-20
        WI(3) = 15971411690431220541.D-20
        WI(4) = 19037877203046567198.D-20
        WI(5) = 19869087157813151863.D-20
        WI(6) = 17972334325952047726.D-20
        WI(7) = 10203571121909080322.D-20
        WI(8) = 84501828581921130722.D-22
        WI(9) = 21467529556997868476.D-24
        GO TO 10

      END IF

      IF (N.EQ.10) THEN
        XI(1) = -98345122025502045873.D-20
        XI(2) = -91446749996879318119.D-20
        XI(3) = -79700500547314513626.D-20
        XI(4) = -64189534981349313375.D-20
        XI(5) = -46376588343242516012.D-20
        XI(6) = -28030431525349494354.D-20
        XI(7) = -11327091328726333942.D-20
        XI(8) = 17437648086722052805.D-21
        XI(9) = 16877498338102917782.D-20
        XI(10) = 40960465258252015313.D-20
        WI(1) = 42278597323639457484.D-21
        WI(2) = 94666349251635366832.D-21
        WI(3) = 13843777024241956101.D-20
        WI(4) = 16932936699837666261.D-20
        WI(5) = 18398357022114735352.D-20
        WI(6) = 17939886390638648260.D-20
        WI(7) = 14468854182396060463.D-20
        WI(8) = 46026485095922891703.D-21
        WI(9) = 11890402956686871419.D-22
        WI(10) = 14148408460516817666.D-25
        GO TO 10

      END IF

      IF (N.EQ.11) THEN
        XI(1) = -98576901837451635280.D-20
        XI(2) = -92621727156102677473.D-20
        XI(3) = -82389243156123939088.D-20
        XI(4) = -68670708816882492198.D-20
        XI(5) = -52549052940365991088.D-20
        XI(6) = -35349156561982307316.D-20
        XI(7) = -18652071146560858606.D-20
        XI(8) = -45389164233559550280.D-21
        XI(9) = 76984180593432347734.D-21
        XI(10) = 24899533750455431614.D-20
        XI(11) = 50711636785486806957.D-20
        WI(1) = 36383684790132198923.D-21
        WI(2) = 81985364434128201418.D-21
        WI(3) = 12133566247788805356.D-20
        WI(4) = 15122112006362489825.D-20
        WI(5) = 16900090791849557413.D-20
        WI(6) = 17240157268363508589.D-20
        WI(7) = 15745585899461757802.D-20
        WI(8) = 97600157144810676257.D-21
        WI(9) = 12496828256639735424.D-21
        WI(10) = 11876318920871395759.D-23
        WI(11) = 80046822403386311030.D-27
        GO TO 10

      END IF

      IF (N.EQ.12) THEN
        XI(1) = -98758247347129831371.D-20
        XI(2) = -93546465146779806654.D-20
        XI(3) = -84528996754470930223.D-20
        XI(4) = -72299594230844519839.D-20
        XI(5) = -57679398168141327066.D-20
        XI(6) = -41683730779892996801.D-20
        XI(7) = -25514627335790291149.D-20
        XI(8) = -10710838211747769681.D-20
        XI(9) = 12720145729326415607.D-21
        XI(10) = 14540842218988328389.D-20
        XI(11) = 33552500235752414908.D-20
        XI(12) = 60838109964484063119.D-20
        WI(1) = 31765161579790701148.D-21
        WI(2) = 71927618746964313778.D-21
        WI(3) = 10742555378156694842.D-20
        WI(4) = 13578811351554214795.D-20
        WI(5) = 15492042553417744038.D-20
        WI(6) = 16300300254834219520.D-20
        WI(7) = 15784577013790806216.D-20
        WI(8) = 12921482926208917372.D-20
        WI(9) = 46096943233133302568.D-21
        WI(10) = 20030610755774790850.D-22
        WI(11) = 95165705752725893549.D-25
        WI(12) = 40143360822128708729.D-28
        GO TO 10

      END IF

      IF (N.EQ.13) THEN
        XI(1) = -98903182721370020265.D-20
        XI(2) = -94288936524363459773.D-20
        XI(3) = -86261843870640242196.D-20
        XI(4) = -75277808759167753869.D-20
        XI(5) = -61972590294795871779.D-20
        XI(6) = -47139332563986024748.D-20
        XI(7) = -31718188942187627557.D-20
        XI(8) = -16854863011308355787.D-20
        XI(9) = -41195843159851553906.D-21
        XI(10) = 71957380142115164738.D-21
        XI(11) = 22223926926874000328.D-20
        XI(12) = 42682885634093164862.D-20
        XI(13) = 71270930856714354732.D-20
        WI(1) = 28069991026027589482.D-21
        WI(2) = 63803895087070663653.D-21
        WI(3) = 95973484361405430270.D-21
        WI(4) = 12264378189747678145.D-20
        WI(5) = 14213612346123977130.D-20
        WI(6) = 15296686007570952707.D-20
        WI(7) = 15358437552921000921.D-20
        WI(8) = 14007635729175637795.D-20
        WI(9) = 87531230524252970103.D-21
        WI(10) = 12989730151883234012.D-21
        WI(11) = 22351943999969127535.D-23
        WI(12) = 65097139765619073344.D-26
        WI(13) = 18257341724040876662.D-29
        GO TO 10

      END IF

      IF (N.EQ.14) THEN
        XI(1) = -99021130855943209687.D-20
        XI(2) = -94895368426058288869.D-20
        XI(3) = -87686856465753704289.D-20
        XI(4) = -77752669471002194917.D-20
        XI(5) = -65594116901081876554.D-20
        XI(6) = -51841232227159879604.D-20
        XI(7) = -37243750660439082187.D-20
        XI(8) = -22693429290756856295.D-20
        XI(9) = -93940943648510570987.D-21
        XI(10) = 16521198218716065629.D-21
        XI(11) = 13919799114797561344.D-20
        XI(12) = 30521886852802066309.D-20
        XI(13) = 52192337126752562221.D-20
        XI(14) = 81957965081548293179.D-20
        WI(1) = 25060310888021301605.D-21
        WI(2) = 57137272611562033779.D-21
        WI(3) = 86434450014324433897.D-21
        WI(4) = 11141118228632175288.D-20
        WI(5) = 13070790263291078499.D-20
        WI(6) = 14310195071194851995.D-20
        WI(7) = 14737968606274298328.D-20
        WI(8) = 14154903694980505066.D-20
        WI(9) = 11456160782223814050.D-20
        WI(10) = 40466499493397342820.D-21
        WI(11) = 21701008894932486895.D-22
        WI(12) = 19960253076851250807.D-24
        WI(13) = 39376501060604877095.D-27
        WI(14) = 76596142918862399780.D-31
        GO TO 10

      END IF

      IF (N.EQ.15) THEN
        XI(1) = -99118619138431485634.D-20
        XI(2) = -95398089203095832045.D-20
        XI(3) = -88874665207045485764.D-20
        XI(4) = -79832886799647722652.D-20
        XI(5) = -68674462209286747178.D-20
        XI(6) = -55907326778454372362.D-20
        XI(7) = -42138595122137487519.D-20
        XI(8) = -28083407355763995168.D-20
        XI(9) = -14649293944496725019.D-20
        XI(10) = -30865949117072113052.D-21
        XI(11) = 75989566859912966734.D-21
        XI(12) = 21425891814116860148.D-20
        XI(13) = 39280262275215780450.D-20
        XI(14) = 62012182191671475949.D-20
        XI(15) = 92858877219218103945.D-20
        WI(1) = 22570991165870390473.D-21
        WI(2) = 51589746641923392000.D-21
        WI(3) = 78401918844466166239.D-21
        WI(4) = 10176234626640128024.D-20
        WI(5) = 12055819130110177262.D-20
        WI(6) = 13377324647273569326.D-20
        WI(7) = 14041818603247829422.D-20
        WI(8) = 13919569003129657925.D-20
        WI(9) = 12562361445602688222.D-20
        WI(10) = 74852662340708470150.D-21
        WI(11) = 10996744175647251144.D-21
        WI(12) = 25513307315040157893.D-23
        WI(13) = 15270418102934789627.D-25
        WI(14) = 21560859319293022163.D-28
        WI(15) = 30032040385910287756.D-32
        GO TO 10

      END IF

      IF (N.EQ.16) THEN
        XI(1) = -99200289748411473927.D-20
        XI(2) = -95820266378296634182.D-20
        XI(3) = -89876661129475763142.D-20
        XI(4) = -81599671254865832971.D-20
        XI(5) = -71315812647978444249.D-20
        XI(6) = -59440032425488487666.D-20
        XI(7) = -46470396871945791541.D-20
        XI(8) = -32991653294098863600.D-20
        XI(9) = -19716091326862980561.D-20
        XI(10) = -76605243443508959615.D-21
        XI(11) = 26155046503992069925.D-21
        XI(12) = 14307776307824938682.D-20
        XI(13) = 29506185654032182160.D-20
        XI(14) = 48403577800553841578.D-20
        XI(15) = 72091584865612160132.D-20
        XI(16) = 10394188783181811718.D-19
        WI(1) = 20484388078614008045.D-21
        WI(2) = 46916532350372347409.D-21
        WI(3) = 71569877291069983495.D-21
        WI(4) = 93424466379672137196.D-21
        WI(5) = 11156011364306951297.D-20
        WI(6) = 12512553084306063601.D-20
        WI(7) = 13329704953113185969.D-20
        WI(8) = 13510959073859290681.D-20
        WI(9) = 12840858805365359846.D-20
        WI(10) = 10016528657871746742.D-20
        WI(11) = 32102655847303900301.D-21
        WI(12) = 18115418480524121495.D-22
        WI(13) = 24274994772381143993.D-24
        WI(14) = 10371321943363515335.D-26
        WI(15) = 10868941709467004901.D-29
        WI(16) = 11117372791599461059.D-33
        GO TO 10

      END IF

   10 CONTINUE
      RETURN

      END !       SUBROUTINE GAUFD(XI,WI,N)
