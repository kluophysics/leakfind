C*==fpecoub.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPECOUB(CMNTMTT,ECOU,NLFP,NSPIN,R2RHO,V,KVMAD,
     &                   NRGNT_VSF_LM,LMRGNT_VSF,RGNT_VSF,NRGNT_VSF)
C   ********************************************************************
C   *                                                                  *
C   *  attention:  energy zero      electro static zero                *
C   *                                                                  *
C   *  calculate the electrostatic potential-energies without the      *
C   *  electron-nuclear interaction in the cell itself.                *
C   *  the energy of the representive atom i is given by               *
C   *                                                                  *
C   *                       rc                                         *
C   *   ECOU(i) =  1/2 (  {  s dr' V(r',LM,i)*R2RHO(r',LM,i,1) }       *
C   *                        0                                         *
C   *                                                                  *
C   *                                    -  Z(i) * VMAD ( ri )     )   *
C   *                                                                  *
C   *                                                                  *
C   *                                      ( {..} = summed over lm )   *
C   *                                                                  *
C   *  V is the coulomb potential of the atom WITHOUT the nuclear      *
C   *          potential of the atom                                   *
C   *  R2RHO(...,1) is the real charge density times r**2              *
C   *                                                                  *
C   *  both expanded into spherical harmonics                          *
C   *                                                                  *
C   *  Z   is the nuclear charge of the atom                           *
C   *                                                                  *
C   *  VMAD ( ri ) is a generalized madelung potential                 *
C   *              = 1/sqrt(4 pi) * V(IRMTIN,1,is)                     *
C   *                - sqrt(4 pi) * 2 * CMNTMTT(1,IT) / rws            *
C   *                                                                  *
C   *                           ( <..> = spherical averaged )          *
C   *                                                                  *
C   *  attention:  this subroutine has to be called bevor the          *
C   *              exchange correlation potential is added to          *
C   *              the potential V.                                    *
C   *              the energy calculated here is splitted into         *
C   *              l-dependent parts to see the l -convergency.        *
C   *                                                                  *
C   *  attention:  in case of shape corrections the contribution of    *
C   *              the coulomb potential the of the nucleus is         *
C   *              analytically cancelled only in the muffin tin       *
C   *              sphere in the interstial region it has to be taken  *
C   *              into account                                        *
C   *                                                                  *
C   *              modified for band structure code                    *
C   *                            b.drittler   jan. 1990                *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:JRCUT,NPAN,R,FLMSF,KLMSF,NRMAX,ISFLM,JRNS1,
     &    DRDI_W_RADINT
      USE MOD_TYPES,ONLY:ITBOT,ITTOP,IMT,NTMAX,NLFPMAX,NLMFPMAX,Z
      USE MOD_CONSTANTS,ONLY:SQRT_4PI
      IMPLICIT NONE
C*--FPECOUB57
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER KVMAD,NLFP,NRGNT_VSF,NSPIN
      REAL*8 CMNTMTT(NLMFPMAX,NTMAX),ECOU(0:(NLFPMAX-1),NTMAX),
     &       R2RHO(NRMAX,NLMFPMAX,NTMAX,3),RGNT_VSF(NRGNT_VSF),
     &       V(NRMAX,NLMFPMAX,NTMAX)
      INTEGER LMRGNT_VSF(NRGNT_VSF,3),NRGNT_VSF_LM(0:NLMFPMAX)
C
C Local variables
C
      REAL*8 DDOT
      REAL*8 ER(:),RHOSP,SN,VM,VMAD
      INTEGER IM,IPAN1,IR,IRCRIT,IRMTIN,IRSF,ISF,ISPIN,IT,J,L,LM,LM2,M
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE ER
C
      ALLOCATE (ER(NRMAX))
      ER(1:NRMAX) = 0.0D0
C
      DO IT = ITBOT,ITTOP
C
         IM = IMT(IT)
         IPAN1 = NPAN(IM)
         IRMTIN = JRCUT(1,IM)
         IRCRIT = JRCUT(IPAN1,IM)
C
         DO L = 0,NLFP - 1
C
            DO IR = 1,IRCRIT
               ER(IR) = 0.0D0
            END DO
C
            DO ISPIN = 1,NSPIN
               IF ( ISPIN.EQ.NSPIN ) THEN
                  SN = 1.0D0
               ELSE
                  SN = -1.0D0
               END IF
C
               DO M = -L,L
                  LM = L*L + L + M + 1
                  DO IR = 1,IRMTIN
                     RHOSP = (R2RHO(IR,LM,IT,1)+SN*R2RHO(IR,LM,IT,NSPIN)
     &                       )/4.0D0
                     IF ( IR.LT.JRNS1(IM) .AND. L.NE.0 ) RHOSP = 0.0D0
                     ER(IR) = ER(IR) + RHOSP*V(IR,LM,IT)
                  END DO
C
C-----------------------------------------------------------------------
C     convolute with shape function
C-----------------------------------------------------------------------
                  DO J = NRGNT_VSF_LM(LM-1) + 1,NRGNT_VSF_LM(LM)
                     LM2 = LMRGNT_VSF(J,2)
                     IF ( KLMSF(LMRGNT_VSF(J,3),IM).GT.0 ) THEN
                        ISF = ISFLM(LMRGNT_VSF(J,3),IM)
                        IF ( LM2.EQ.1 ) THEN
                           DO IR = IRMTIN + 1,IRCRIT
                              IRSF = IR - IRMTIN
                              RHOSP = (R2RHO(IR,LM,IT,1)
     &                                +SN*R2RHO(IR,LM,IT,NSPIN))/2.0D0
C-----------------------------------------------------------------------
C     remember that in the interstial -2z/r has to be taken into account
C-----------------------------------------------------------------------
                              ER(IR) = ER(IR) + RHOSP*RGNT_VSF(J)
     &                                 *FLMSF(IRSF,ISF,IM)
     &                                 *(V(IR,1,IT)/2.0D0-Z(IT)/R(IR,IM)
     &                                 *SQRT_4PI)
                           END DO
C
                        ELSE
C
                           DO IR = IRMTIN + 1,IRCRIT
                              IRSF = IR - IRMTIN
                              RHOSP = (R2RHO(IR,LM,IT,1)
     &                                +SN*R2RHO(IR,LM,IT,NSPIN))/2.0D0
                              ER(IR) = ER(IR) + RHOSP*RGNT_VSF(J)
     &                                 *FLMSF(IRSF,ISF,IM)*V(IR,LM2,IT)
     &                                 /2.0D0
                           END DO
                        END IF
                     END IF
C
                  END DO
C
               END DO
C
            END DO
C
C-----------------------------------------------------------------------
C     now integrate
C-----------------------------------------------------------------------
C
            ECOU(L,IT) = DDOT(IRCRIT,ER(1),1,DRDI_W_RADINT(1,IM),1)
C
         END DO
C
C-----------------------------------------------------------------------
C     calculate the madelung potential
C-----------------------------------------------------------------------
C
         VMAD = V(IRMTIN,1,IT)/SQRT_4PI - SQRT_4PI*2.0D0*CMNTMTT(1,IT)
     &          /R(IRMTIN,IM)
C
C-----------------------------------------------------------------------
C     add to ECOU
C-----------------------------------------------------------------------
C
         ECOU(0,IT) = ECOU(0,IT) - Z(IT)*VMAD/2.0D0
C
C-----------------------------------------------------------------------
C     option to calculate full generalized madelung potential
C                                  rc
C     vm(rn) = vmad +2*sqrt(4*pi)* s  dr*r*rho(lm=1,r)
C                                  0
C-----------------------------------------------------------------------
         IF ( KVMAD.EQ.1 ) THEN
C
            VM = 0D0
            DO IR = 2,IRMTIN
               VM = VM + DRDI_W_RADINT(IR,IM)*R2RHO(IR,1,IT,1)/R(IR,IM)
            END DO
C
            VM = 2.0D0*SQRT_4PI*VM + VMAD
C
            WRITE (6,99001) IT,VM
C
         END IF
C
      END DO
C
99001 FORMAT (13x,'full generalized madelung pot. for atom',1X,I3,1X,
     &        ': ',D16.8)
      END
