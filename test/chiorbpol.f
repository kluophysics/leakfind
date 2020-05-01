C*==chiorbpol.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHIORBPOL(AXCN,GAMMAOP,ORBSQINT,R2DRDI,RPW,Z,ORBPOL,
     &                     SYMT,IT,IRTOP,NLMAX,NRMAX)
C   ********************************************************************
C   *                                                                  *
C   * calculate orbital polarisation induced by an electron            *
C   *                                at the Fermi surface              *
C   *                                                                  *
C   * based on scfoppot.f                                              *
C   * called by CHIGAMMA                                               *
C   *                                                                  *
C   * ORBPOL = BROOKS     M.S.S. Brooks' OP-scheme                     *
C   * OPRPOL = BROOKS-F   DEALS WITH F ORBITALS                        *
C   *                                                                  *
C   * 01/08/99  HF                                                     *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:IMT
      USE MOD_CONSTANTS,ONLY:RY_EV
      IMPLICIT NONE
C*--CHIORBPOL21
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IRTOP,IT,NLMAX,NRMAX,Z
      CHARACTER*10 ORBPOL
      CHARACTER*2 SYMT
      REAL*8 AXCN(NRMAX,2,NLMAX),GAMMAOP(NRMAX,2),ORBSQINT(0:NLMAX,2),
     &       R2DRDI(NRMAX),RPW(NRMAX,NLMAX*2)
C
C Local variables
C
      REAL*8 BRACAH,FCLMB,PRE(2*(NLMAX-1)),RHORSQ,SG(NRMAX),SL(NRMAX),
     &       TG(NRMAX),TL(NRMAX)
      INTEGER I,ICALL,IM,IS,LF,LFMAX,LOP,NLOP
      LOGICAL OP
      CHARACTER*2 RACPAR
C
C*** End of declarations rewritten by SPAG
C
      DATA ICALL/0/
      ICALL = ICALL + 1
C
      CALL RINIT(NRMAX*2*NLMAX,AXCN)
C
      IM = IMT(IT)
C
      IF ( Z.GE.0 ) THEN
         LOP = 2
         RACPAR = 'B '
         PRE(2) = +9.0D0/441.0D0
         PRE(4) = -5.0D0/441.0D0
         LFMAX = 4
         OP = .TRUE.
         IF ( ORBPOL(1:8).EQ.'BROOKS-F' ) THEN
            LOP = 3
            RACPAR = 'E3'
            PRE(2) = +5.0D0/(225.0D0*3.0D0)
            PRE(4) = +6.0D0/(1089.0D0*3.0D0)
            PRE(6) = -91.0D0/(7361.64D0*3.0D0)
            LFMAX = 6
            IF ( Z.LT.57 ) OP = .FALSE.
         END IF
         NLOP = LOP + 1
C
         IF ( OP ) THEN
C -------------------------------------------------- calculate potential
            DO IS = 1,2
               DO LF = 2,LFMAX,2
C -------------------------------------- evaluate first radial integrals
                  DO I = 1,IRTOP
                     RHORSQ = 2.0D0*GAMMAOP(I,IS)*R2DRDI(I)
                     TL(I) = RHORSQ*RPW(I,LF)
                     TG(I) = RHORSQ/RPW(I,LF+1)
                  END DO
C
                  CALL RRADINT_R(IM,TL,SL)
                  CALL RRADINT_R(IM,TG,SG)
C
C ------------------------------------ add prefactor 1/r**(L+1) and r**L
                  DO I = 1,IRTOP
                     SL(I) = SL(I)/RPW(I,LF+1) + (SG(IRTOP)-SG(I))
     &                       *RPW(I,LF)
                  END DO
C ----------------------------------------- check Coulomb integrals F^LF
                  DO I = 1,IRTOP
                     SG(I) = GAMMAOP(I,IS)*SL(I)*R2DRDI(I)
                  END DO
                  CALL RRADINT(IM,SG,FCLMB)
C ----------------------------------------------------- set up potential
                  DO I = 1,IRTOP
                     AXCN(I,IS,NLOP) = AXCN(I,IS,NLOP) - PRE(LF)*SL(I)
                  END DO
C
C ----------------------------------------- check Racah parameter B / E3
                  DO I = 1,IRTOP
                     SG(I) = GAMMAOP(I,IS)*AXCN(I,IS,NLOP)*R2DRDI(I)
                  END DO
                  CALL RRADINT(IM,SG,BRACAH)
                  BRACAH = BRACAH/ORBSQINT(NLOP,IS)
C
                  IF ( ICALL.EQ.1 ) THEN
C                    WRITE (6,*) 'PRE(LF)=',PRE(LF)
                     WRITE (6,99001) LF,FCLMB*1000D0,FCLMB*RY_EV,IT,
     &                               SYMT,IS
                     WRITE (6,99002) RACPAR,BRACAH*1000D0,BRACAH*RY_EV,
     &                               IT,SYMT,IS
                     WRITE (6,99003) ORBSQINT(NLOP,IS),IT,SYMT,IS
                  END IF
C
               END DO
C ======= DO LF=2,LFMAX,2 =======
C
            END DO
C === DO IS=1,2 ===
C
         END IF
C === IF (OP) ===
C
      END IF
C === Z .GE. 0 ===
C
99001 FORMAT (' COULOMB-integral F^',I1,' = ',F8.3,' mRy = ',F8.3,' eV',
     &        ' for  IT=',I2,2X,A,'  MS=',I2)
99002 FORMAT (' RACAH-parameter  ',A2,1X,' = ',F8.3,' mRy = ',F8.3,
     &        ' eV',' for  IT=',I2,2X,A,'  MS=',I2)
99003 FORMAT (' <l_z>              ',1X,' = ',F8.3,18X,' for  IT=',I2,
     &        2X,A,'  MS=',I2)
C
      END
