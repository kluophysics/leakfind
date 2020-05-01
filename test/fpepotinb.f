C*==fpepotinb.f    processed by SPAG 6.70Rc at 21:35 on 19 Dec 2016
      SUBROUTINE FPEPOTINB(EPOTIN,NSPIN,R2RHO,VT,BT,VNST,BNST)
C   ********************************************************************
C   *                                                                  *
C   * attention:  energy zero      electro static zero                 *
C   *                                                                  *
C   *             since input potential and single particle energies   *
C   *             are using muffin tin zero as zero the energy shift   *
C   *             is cancelled in the kinetic energy contribution !    *
C   *                                                                  *
C   *                                                                  *
C   * calculate the energy of the input potential                      *
C   * the energy for the representive atom i is given by               *
C   *                                                                  *
C   *                           rws                                    *
C   *   epotin(i) = - sqrt(4 pi) {  dr' vm2z(r',i)*r2rho(r',1,i)       *
C   *                            0                                     *
C   *                                                                  *
C   * in case of non spherical input potential one has to add          *
C   *                                                                  *
C   *             rirt                                                 *
C   *        {  -  {  dr' vns(r',lm,i)*r2rho(r',1,lm,1)   }            *
C   *             rmin                                                 *
C   *                                    (summed over lm)              *
C   *                                                                  *
C   * remember:  the non spherical part of the input potential is      *
C   *            different from zero only between r(jrns1) and r(irt)  *
C   *                                                                  *
C   * attention:  vm2z is the spherically averaged input potential,    *
C   *             vns contains the non spherical contribution of the   *
C   *             potential and r2rho(.,1) is the  real charge density *
C   *             times r**2 - vns and r2rho are expanded into         *
C   *             spherical harmonics.                                 *
C   *                                                                  *
C   * remember:   in case of shape corrections  the contribution of    *
C   *             the nuclear potential - 2*z/r has to be explicitly   *
C   *             taken into account between muffin tin sphere and     *
C   *             circum scribed sphere.                               *
C   *             only within the muffin tin sphere this term is       *
C   *             analytically cancelled wtih the contribution of      *
C   *             the coulomb potential                                *
C   *                                                                  *
C   *                                                                  *
C   *             modified for non spherical potential and shape       *
C   *             corrections                                          *
C   *                                                                  *
C   *                           b.drittler   oct. 1989                 *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:IMT,Z,NTMAX,NLMFPMAX,NLFPMAX
      USE MOD_RMESH,ONLY:R,JRNS1,NRMAX,JRCUT,NPAN,JRNSMIN,DRDI_W_RADINT
      USE MOD_TYPES,ONLY:ITBOT,ITTOP,NLFP
      USE MOD_CONSTANTS,ONLY:SQRT_4PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NSPIN
      REAL*8 BNST(JRNSMIN:NRMAX,NLMFPMAX,NTMAX),BT(NRMAX,NTMAX),
     &       EPOTIN(NTMAX),R2RHO(NRMAX,NLMFPMAX,NTMAX,3),
     &       VNST(JRNSMIN:NRMAX,NLMFPMAX,NTMAX),VT(NRMAX,NTMAX)
C
C Local variables
C
      REAL*8 AUX,ENS(:,:),ER(:),R2RHOD,R2RHOU,VD,VM2ZD,VM2ZU,VU
      INTEGER I,IM,IPAN1,IRCRIT,IRMTIN,IT,JRMIN1,L1,LM,M1
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE ENS,ER
C
      ALLOCATE (ENS(0:(NLFPMAX-1),NTMAX),ER(NRMAX))
C
      DO IT = ITBOT,ITTOP
C
         IM = IMT(IT)
         IPAN1 = NPAN(IM)
         IRCRIT = JRCUT(IPAN1,IM)
         IRMTIN = JRCUT(1,IM)
C
C   the variable-names are misleading :  ?????????????
C   Vm2Z(.,ipotu) is the spin-down-component of Vold,
C   Vm2Z(.,ipotd) is the spin-up-component   of Vold
C
C-----------------------------------------------------------------------
C     calculate charge density times INPUT potential
C-----------------------------------------------------------------------
C
         AUX = 0D0
         DO I = 1,IRMTIN
            R2RHOU = (R2RHO(I,1,IT,1)-R2RHO(I,1,IT,NSPIN))/2.0D0
            R2RHOD = (R2RHO(I,1,IT,1)+R2RHO(I,1,IT,NSPIN))/2.0D0
C
            VM2ZU = VT(I,IT) - BT(I,IT) + 2D0*Z(IT)/R(I,IM)
            VM2ZD = VT(I,IT) + BT(I,IT) + 2D0*Z(IT)/R(I,IM)
C
            ER(I) = -(R2RHOU*VM2ZU+R2RHOD*VM2ZD)*SQRT_4PI
            AUX = AUX + ER(I)*DRDI_W_RADINT(I,IM)
         END DO
C
C-----------------------------------------------------------------------
C     remember the form of vm2z between mt sphere and rirc ?????????????
C-----------------------------------------------------------------------
C
         IF ( IPAN1.GT.1 ) THEN
            DO I = IRMTIN + 1,IRCRIT
               R2RHOU = (R2RHO(I,1,IT,1)-R2RHO(I,1,IT,NSPIN))/2.0D0
               R2RHOD = (R2RHO(I,1,IT,1)+R2RHO(I,1,IT,NSPIN))/2.0D0
C
               VU = VT(I,IT) - BT(I,IT)
               VD = VT(I,IT) + BT(I,IT)
C
               ER(I) = -(R2RHOU*VU+R2RHOD*VD)*SQRT_4PI
               AUX = AUX + ER(I)*DRDI_W_RADINT(I,IM)
            END DO
         END IF
C
         EPOTIN(IT) = AUX
         ENS(0,IT) = AUX
C-----------------------------------------------------------------------
C     add non spher. contribution in case of non spher. input potential
C-----------------------------------------------------------------------
         DO L1 = 1,(NLFP-1)
            ENS(L1,IT) = 0.0D0
         END DO
C
         JRMIN1 = JRNS1(IM)
         IF ( JRMIN1.LE.IRMTIN ) THEN
C
            DO L1 = 1,(NLFP-1)
               DO I = 1,NRMAX
                  ER(I) = 0.0D0
               END DO
               DO M1 = -L1,L1
                  LM = L1*(L1+1) + M1 + 1
C
C-----------------------------------------------------------------------
C     calculate charge density times potential
C-----------------------------------------------------------------------
C
                  DO I = JRMIN1,IRCRIT
                     R2RHOU = (R2RHO(I,LM,IT,1)-R2RHO(I,LM,IT,NSPIN))
     &                        /2.0D0
                     R2RHOD = (R2RHO(I,LM,IT,1)+R2RHO(I,LM,IT,NSPIN))
     &                        /2.0D0
C
                     VU = VNST(I,LM,IT) - BNST(I,LM,IT)
                     VD = VNST(I,LM,IT) + BNST(I,LM,IT)
C
                     ER(I) = ER(I) - R2RHOU*VU - R2RHOD*VD
                  END DO
               END DO
C
               AUX = 0D0
               DO I = JRMIN1,IRCRIT
                  AUX = AUX + ER(I)*DRDI_W_RADINT(I,IM)
               END DO
C
               EPOTIN(IT) = EPOTIN(IT) + AUX
               ENS(L1,IT) = AUX
C
            END DO
C
         END IF
C
      END DO
C
      DO IT = ITBOT,ITTOP
         WRITE (6,99001) IT,EPOTIN(IT)
         WRITE (6,99002) (L1,ENS(L1,IT),L1=0,(NLFP-1))
      END DO
C
99001 FORMAT (/,5X,'atom type IT=',I2,'  EPOTIN =',f15.8)
99002 FORMAT ((9X,4(i2,':',f14.6)))
      END
