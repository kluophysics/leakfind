C*==rotmat_rylm.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE ROTMAT_RYLM(NL_LOC,NLMMAX_LOC,ALFDEG,BETDEG,GAMDEG,
     &                       DROT_RYLM)
C   ********************************************************************
C   *                                                                  *
C   *  calculate rotation matrix  DROT_RYLM(ALFDEG,BETDEG,GAMDEG)      *
C   *  for REAL spherical harmonics for any  NL_LOC                    *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_FILES,ONLY:IPRINT
      IMPLICIT NONE
C*--ROTMAT_RYLM14
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALFDEG,BETDEG,GAMDEG
      INTEGER NLMMAX_LOC,NL_LOC
      REAL*8 DROT_RYLM(NLMMAX_LOC,NLMMAX_LOC)
C
C Local variables
C
      COMPLEX*16 DROT(:,:),RC_LM(:,:),W1(:,:),W2(:,:)
      INTEGER I,IA_ERR,J,NLM_LOC
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DROT,W1,W2,RC_LM
C
C ----------------------------------------------------------------------
C        allow to deal with arbitrary angular momentum cut off
C ----------------------------------------------------------------------
C
      NLM_LOC = NL_LOC**2
C
      IF ( NLM_LOC.GT.NLMMAX_LOC ) THEN
         WRITE (6,99001) NL_LOC,NLM_LOC,NLMMAX_LOC
         STOP
      END IF
C
C-----------------------------------------------------------------------
C   transformation matrix  RC_LM = RC:  complex -> real spher. harm.
C-----------------------------------------------------------------------
C
      ALLOCATE (RC_LM(NLM_LOC,NLM_LOC))
C
      CALL CALC_RC(NL_LOC-1,RC_LM,NLM_LOC,1)
C
C ----------------------------------------------------------------------
C         rotation matrix  DROT for COMPLEX spherical harmonics
C ----------------------------------------------------------------------
C
      ALLOCATE (DROT(NLM_LOC,NLM_LOC))
C
      CALL ROTMAT(NL_LOC,1,ALFDEG,BETDEG,GAMDEG,DROT,NLM_LOC)
C
C ----------------------------------------------------------------------
C         rotation matrix  DROT_RYLM for REAL spherical harmonics
C ----------------------------------------------------------------------
      ALLOCATE (W1(NLM_LOC,NLM_LOC),W2(NLM_LOC,NLM_LOC),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc: ROTMAT_YRLM -> LLIM'
C
      DROT_RYLM(:,:) = 0D0
C
      CALL CMATTRANS(NLM_LOC,NLM_LOC,W2,RC_LM,.TRUE.,DROT,'SAS+',W1)
C
      DO J = 1,NLM_LOC
         DO I = 1,NLM_LOC
C
            IF ( ABS(DIMAG(W1(I,J))).GT.1D-8 )
     &            STOP 'in <ROTMAT_YRLM>:  |Im(DROT(I,J)| > 1D-8'
C
            DROT_RYLM(I,J) = DREAL(W1(I,J))
C
            IF ( ABS(DROT_RYLM(I,J)).LT.1D-12 ) DROT_RYLM(I,J) = 0D0
C
         END DO
      END DO
C
C-----------------------------------------------------------------------
      IF ( IPRINT.LT.6 ) RETURN
C-----------------------------------------------------------------------
      WRITE (6,99002) NL_LOC,ALFDEG,BETDEG,GAMDEG
      CALL RMATSTRUCT(' ',DROT_RYLM,NLMMAX_LOC,NLMMAX_LOC,1,1,0,1D-8,6)
C-----------------------------------------------------------------------
99001 FORMAT (/,10X,'ERROR in <ROTMAT_RYLM>:',/,10X,'NL_LOC     =',I4,/,
     &        /,10X,'NLM_LOC    =',I4,/,10X,'NLMMAX_LOC =',I4,/,10X,
     &        'NLM_LOC <> NLMMAX_LOC',/)
99002 FORMAT (/,10X,'rotation matrix w.r.t. REAL spherical harmonics',/,
     &        10X,'for NL_LOC =',I2,'  ALF =',F6.3,'  BET =',F6.3,
     &        '  GAM =',F6.3)
      END
