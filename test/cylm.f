C*==cylm.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CYLM(V,LOMAX,Y)
C
C CALCULATES SINES AND COSINES OF THE POLAR ANGLES OF THE VECTOR V.
C
      IMPLICIT NONE
C*--CYLM7
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CYLM')
C
C Dummy arguments
C
      INTEGER LOMAX
      REAL*8 V(3)
      COMPLEX*16 Y((LOMAX+1)**2)
C
C Local variables
C
      REAL*8 C1L,C2L,CD,CFI,CMFI,CN,CSR,CTH,CYP,FPI,ONE,P(:,:),RF,SFI,
     &       SGNM,SMFI,SNULL,STH,TCTH,XY,XYZ,YI,YR,ZERO
      INTEGER I,IA_ERR,IDWN,L,L1,L2,LM,LM2,M,M1
C
C*** End of declarations rewritten by SPAG
C
      DATA ZERO/0.D0/,SNULL/1.D-10/,ONE/1.D0/,FPI/12.5663706D0/
C
      ALLOCATABLE P
C
      ALLOCATE (P(LOMAX+2,LOMAX+2),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: P')
C
      FPI = 4.D0*ACOS(-1.D0)
      XY = V(1)**2 + V(2)**2
      XYZ = XY + V(3)**2
C
      IF ( XY.GT.SNULL ) THEN
         XY = SQRT(XY)
         XYZ = SQRT(XYZ)
         CTH = V(3)/XYZ
         STH = XY/XYZ
         CFI = V(1)/XY
         SFI = V(2)/XY
      ELSE
         CTH = ONE
         IF ( V(3).LT.ZERO ) CTH = -CTH
         STH = ZERO
         CFI = ONE
         SFI = ZERO
      END IF
      RF = ONE/SQRT(FPI)
      YR = RF
      Y(1) = YR
      I = 1
      P(1,1) = ONE
      P(2,1) = CTH
      C2L = CTH
      TCTH = CTH + CTH
      L1 = 2
      L = 1
 100  CONTINUE
      M = 1
      I = I + L
      IDWN = I + 2
      M1 = 2
      L2 = L
      L = L1
      CMFI = ONE
      SMFI = ZERO
      L1 = L + 1
      LM = L2
      LM2 = L
      CD = ONE
      C2L = C2L + TCTH
      SGNM = ONE
C
C.....RECURSE UPWARD IN L.
C
 200  CONTINUE
      P(L1,M) = (C2L*P(L,M)-LM*P(L2,M))/LM2
C
      C1L = (LM+1)*CTH
      P(L,M1) = ZERO
C
C.....RECURSE UPWARD IN M.
C
      IF ( ABS(STH).GE.SNULL ) P(L,M1) = (C1L*P(L,M)-LM2*P(L1,M))/STH
 300  CONTINUE
      I = I + 1
      IDWN = IDWN - 1
      CSR = SQRT((2*L-ONE)/(FPI*CD))
      CYP = SGNM*CSR*P(L,M)
      YR = CYP*CMFI
      YI = CYP*SMFI
      Y(I) = DCMPLX(YR,YI)
      IF ( IDWN.NE.I ) Y(IDWN) = SGNM*DCMPLX(YR,-YI)
      CN = CMFI
      CMFI = CN*CFI - SFI*SMFI
      SMFI = SFI*CN + SMFI*CFI
      M = M1
      M1 = M + 1
      LM2 = LM2 - 1
      LM = LM + 1
      CD = CD*LM*LM2
      SGNM = -SGNM
      IF ( M.LT.L ) GOTO 200
      IF ( M.EQ.L ) GOTO 300
      IF ( L.LE.LOMAX ) GOTO 100
C
      DEALLOCATE (P,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC')
      END
