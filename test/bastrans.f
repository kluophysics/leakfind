C*==bastrmat.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE BASTRMAT(LMAX,RC,CREL,RREL,NKMMAX)
C   ********************************************************************
C   *                                                                  *
C   *  initialize matrices RC, CREL and RREL for basis transformation  *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   * CREL:    COMPLEX (L,M,S)  to  (KAP,MUE) - representation         *
C   *                 |LAM> = sum[LC] |LC> * CREL(LC,LAM)              *
C   *------------------------------------------------------------------*
C   *   RC:    REAL to  COMPLEX (L,M,S) - representation               *
C   *                 |LC> = sum[LR] |LR> * RC(LR,LC)                  *
C   *------------------------------------------------------------------*
C   * RREL:    REAL (L,M,S)  to  (KAP,MUE) - representation            *
C   *                 |LAM> = sum[LR] |LR> * RREL(LR,LAM)              *
C   * -----------------------------------------------------------------*
C   *    standard indexing:                                            *
C   *    IKM  = L*2*(J+1/2) + J + MUE + 1                              *
C   *    LM   = L*(L+1)     +     M   + 1                              *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:CGC
      USE MOD_CONSTANTS,ONLY:C1,C0
      IMPLICIT NONE
C*--BASTRMAT26
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER LMAX,NKMMAX
      COMPLEX*16 CREL(NKMMAX,NKMMAX),RC(NKMMAX,NKMMAX),
     &           RREL(NKMMAX,NKMMAX)
C
C Local variables
C
      INTEGER IKM,JP05,K,L,LM,LNR,M,MUEM05,MUEP05,NK,NKM,NLM
C
C*** End of declarations rewritten by SPAG
C
      NK = 2*(LMAX+1) + 1
      NLM = (LMAX+1)**2
      NKM = 2*NLM
C
      IF ( NKM.GT.NKMMAX ) STOP '<BASTRMAT>  NKM > NKMMAX'
C
C ----------------------------------------------------------------------
C CREL  transforms from  COMPLEX (L,M,S)  to  (KAP,MUE) - representation
C                 |LAM> = sum[LC] |LC> * CREL(LC,LAM)
C ----------------------------------------------------------------------
      CALL CINIT(NKMMAX*NKMMAX,CREL)
C
      LM = 0
      DO LNR = 0,LMAX
         DO M = -LNR,LNR
            LM = LM + 1
C
            IKM = 0
            DO K = 1,NK
               L = K/2
               IF ( 2*L.EQ.K ) THEN
                  JP05 = L
               ELSE
                  JP05 = L + 1
               END IF
C
               DO MUEM05 = -JP05,(JP05-1)
                  MUEP05 = MUEM05 + 1
                  IKM = IKM + 1
C
                  IF ( L.EQ.LNR ) THEN
                     IF ( MUEP05.EQ.M ) CREL(LM,IKM) = CGC(IKM,1)
                     IF ( MUEM05.EQ.M ) CREL(LM+NLM,IKM) = CGC(IKM,2)
                  END IF
C
               END DO
            END DO
C
         END DO
      END DO
C
C ----------------------------------------------------------------------
C    RC  transforms from  REAL to  COMPLEX (L,M,S) - representation
C                 |LC> = sum[LR] |LR> * RC(LR,LC)
C ----------------------------------------------------------------------
C
      CALL CALC_RC(LMAX,RC,NKMMAX,2)
C
C ----------------------------------------------------------------------
C RREL  transforms from   REAL (L,M,S)  to  (KAP,MUE) - representation
C                 |LAM> = sum[LR] |LR> * RREL(LR,LAM)
C ----------------------------------------------------------------------
C
      CALL ZGEMM('N','N',NKM,NKM,NKM,C1,RC,NKMMAX,CREL,NKMMAX,C0,RREL,
     &           NKMMAX)
C
      END
C*==calc_rc.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE CALC_RC(LMAX,RC,NXMMAX,NSPIN)
C   ********************************************************************
C   *                                                                  *
C   *  RC  transforms from  REAL to  COMPLEX (L,M,S) - representation  *
C   *               |LC> = sum[LR] |LR> * RC(LR,LC)                    *
C   *                                                                  *
C   *  standard indexing:    LM   = L*(L+1) + M + 1                    *
C   *                                                                  *
C   *  NSPIN = 1:  create the standard non-relativistic matrix         *
C   *  NSPIN = 2:  assume a spin-depenent representation and           *
C   *              copy the matrix for ISPIN=1 to ISPIN=2              *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:CI,C1,C0
      IMPLICIT NONE
C*--CALC_RC127
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER LMAX,NSPIN,NXMMAX
      COMPLEX*16 RC(NXMMAX,NXMMAX)
C
C Local variables
C
      INTEGER I,I1,I2,J,L,M,NLM_LOC
      REAL*8 W
C
C*** End of declarations rewritten by SPAG
C
      NLM_LOC = (LMAX+1)**2
C
      IF ( NLM_LOC.GT.NXMMAX ) STOP '<CALC_RC>  NLM_LOC > NXMMAX'
C
C ----------------------------------------------------------------------
C
      RC(:,:) = C0
C
      W = 1.0D0/SQRT(2.0D0)
C
      DO L = 0,LMAX
         DO M = -L,L
            I = L*(L+1) + M + 1
            J = L*(L+1) - M + 1
C
            IF ( M.LT.0 ) THEN
               RC(I,I) = -CI*W
               RC(J,I) = W
            ELSE IF ( M.EQ.0 ) THEN
               RC(I,I) = C1
            ELSE
               RC(I,I) = W*(-1.0D0)**M
               RC(J,I) = CI*W*(-1.0D0)**M
            END IF
         END DO
      END DO
C
C-----------------------------------------------------------------------
C
      IF ( NSPIN.EQ.1 ) RETURN
C
C-----------------------------------------------------------------------
C             create 2nd spin-diagonal block by copying
C-----------------------------------------------------------------------
C
      IF ( 2*NLM_LOC.GT.NXMMAX ) STOP '<CALC_RC>  2*NLM_LOC > NXMMAX'
C
      I1 = NLM_LOC + 1
      I2 = NLM_LOC + NLM_LOC
C
      RC(I1:I2,I1:I2) = RC(1:NLM_LOC,1:NLM_LOC)
C
      END
C*==changerep.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE CHANGEREP(N,M,A,MODE,B)
C   ********************************************************************
C   *                                                                  *
C   *   change the representation of matrix A and store in B           *
C   *   according to MODE:                                             *
C   *                                                                  *
C   *   RLM>REL   non-relat. REAL spher. harm.  >   (kappa,mue)        *
C   *   REL>RLM   (kappa,mue)  > non-relat. REAL spher. harm.          *
C   *   CLM>REL   non-relat. CMPLX. spher. harm.  >   (kappa,mue)      *
C   *   REL>CLM   (kappa,mue)  > non-relat. CMPLX. spher. harm.        *
C   *   RLM>CLM   non-relat. REAL spher. harm.  >  CMPLX. spher. harm. *
C   *   CLM>RLM   non-relat. CMPLX. spher. harm.  >  REAL spher. harm. *
C   *                                                                  *
C   *   the non-relat. representations include the  spin index         *
C   *                                                                  *
C   *   for IPRINT >= 3 the old and new matrices, A and B, are printed *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,WKM_LPK,RC,CREL,RREL
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_CONSTANTS,ONLY:C1,C0
      IMPLICIT NONE
C*--CHANGEREP221
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER M,N
      CHARACTER*7 MODE
      COMPLEX*16 A(M,M),B(M,M)
C
C*** End of declarations rewritten by SPAG
C
C========================================================================
C                         check arguments
C========================================================================
C
      IF ( N.NE.NKM ) STOP '<CHANGEREP>  N <> NKM'
      IF ( M.NE.NKMMAX ) STOP '<CHANGEREP>  M <> NKMMAX'
C
C========================================================================
C                    perform transformation
C========================================================================
C
      IF ( MODE.EQ.'REL>RLM' ) THEN
C
         CALL ZGEMM('N','N',N,N,N,C1,RREL,M,A,M,C0,WKM_LPK,M)
         CALL ZGEMM('N','C',N,N,N,C1,WKM_LPK,M,RREL,M,C0,B,M)
C
      ELSE IF ( MODE.EQ.'RLM>REL' ) THEN
C
         CALL ZGEMM('C','N',N,N,N,C1,RREL,M,A,M,C0,WKM_LPK,M)
         CALL ZGEMM('N','N',N,N,N,C1,WKM_LPK,M,RREL,M,C0,B,M)
C
      ELSE IF ( MODE.EQ.'REL>CLM' ) THEN
C
         CALL ZGEMM('N','N',N,N,N,C1,CREL,M,A,M,C0,WKM_LPK,M)
         CALL ZGEMM('N','C',N,N,N,C1,WKM_LPK,M,CREL,M,C0,B,M)
C
      ELSE IF ( MODE.EQ.'CLM>REL' ) THEN
C
         CALL ZGEMM('C','N',N,N,N,C1,CREL,M,A,M,C0,WKM_LPK,M)
         CALL ZGEMM('N','N',N,N,N,C1,WKM_LPK,M,CREL,M,C0,B,M)
C
      ELSE IF ( MODE.EQ.'CLM>RLM' ) THEN
C
         CALL ZGEMM('N','N',N,N,N,C1,RC,M,A,M,C0,WKM_LPK,M)
         CALL ZGEMM('N','C',N,N,N,C1,WKM_LPK,M,RC,M,C0,B,M)
C
      ELSE IF ( MODE.EQ.'RLM>CLM' ) THEN
C
         CALL ZGEMM('C','N',N,N,N,C1,RC,M,A,M,C0,WKM_LPK,M)
         CALL ZGEMM('N','N',N,N,N,C1,WKM_LPK,M,RC,M,C0,B,M)
C
      ELSE
         WRITE (6,*) ' MODE = ',MODE
         STOP 'in <CHANGEREP>  MODE not allowed'
      END IF
C
C========================================================================
C           print out matrices depending on print level
C========================================================================
C
      IF ( IPRINT.LT.5 ) RETURN
C
C========================================================================
C
      IF ( MODE.EQ.'REL>RLM' ) THEN
         CALL CMATSTRUCT('<CHANGEREP>: REL',A,N,M,3,3,0,1D-8,6)
         CALL CMATSTRUCT('<CHANGEREP>: RLM',B,N,M,2,2,0,1D-8,6)
      ELSE IF ( MODE.EQ.'RLM>REL' ) THEN
         CALL CMATSTRUCT('<CHANGEREP>: RLM',A,N,M,2,2,0,1D-8,6)
         CALL CMATSTRUCT('<CHANGEREP>: REL',B,N,M,3,3,0,1D-8,6)
      ELSE IF ( MODE.EQ.'REL>CLM' ) THEN
         CALL CMATSTRUCT('<CHANGEREP>: REL',A,N,M,3,3,0,1D-8,6)
         CALL CMATSTRUCT('<CHANGEREP>: CLM',B,N,M,2,2,0,1D-8,6)
      ELSE IF ( MODE.EQ.'CLM>REL' ) THEN
         CALL CMATSTRUCT('<CHANGEREP>: CLM',A,N,M,2,2,0,1D-8,6)
         CALL CMATSTRUCT('<CHANGEREP>: REL',B,N,M,3,3,0,1D-8,6)
      ELSE IF ( MODE.EQ.'CLM>RLM' ) THEN
         CALL CMATSTRUCT('<CHANGEREP>: CLM',A,N,M,1,1,0,1D-8,6)
         CALL CMATSTRUCT('<CHANGEREP>: RLM',B,N,M,1,1,0,1D-8,6)
      ELSE IF ( MODE.EQ.'RLM>CLM' ) THEN
         CALL CMATSTRUCT('<CHANGEREP>: RLM',A,N,M,1,1,0,1D-8,6)
         CALL CMATSTRUCT('<CHANGEREP>: CLM',B,N,M,1,1,0,1D-8,6)
      ELSE
         WRITE (6,*) ' MODE = ',MODE
         STOP 'in <CHANGEREP>  MODE not allowed'
      END IF
C
      END
C*==cmattrans.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE CMATTRANS(N,M,W,S,UNITARY,A,KEY,B)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the matrix-matrix transformation                      *
C   *                                                                  *
C   *     B = S  * A * S+     for  KEY = 'SAS+'                        *
C   *     B = S+ * A * S      for  KEY = 'S+AS'                        *
C   *                                                                  *
C   *     for UNITARY = .FALSE. (anti-unitary transformation)          *
C   *                                                                  *
C   *     B = S  * A+ * S+    for  KEY = 'SAS+'                        *
C   *     B = S+ * A+ * S     for  KEY = 'S+AS'                        *
C   *                                                                  *
C   *   A       matrix to be transformed  - unchanged on exit          *
C   *   B       result of transformation                               *
C   *   S       transformation matrix                                  *
C   *   W       work space                                             *
C   *                                                                  *
C   *   A,B,S,W complex  SQUARE  N x N - matrices                      *
C   *   N       dimension  of matrices                                 *
C   *   M       array size of matrices with M >= N                     *
C   *   UNITARY logical  .T. == unitary                                *
C   *                                                                  *
C   *   driver routine using LAPACK routines                           *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0,C1
      IMPLICIT NONE
C*--CMATTRANS349
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*4 KEY
      INTEGER M,N
      LOGICAL UNITARY
      COMPLEX*16 A(M,M),B(M,M),S(M,M),W(M,M)
C
C*** End of declarations rewritten by SPAG
C
C========================================================================
C                        check arguments
C========================================================================
C
      IF ( N.GT.M ) STOP '<CMATTRANS>  N > M'
C
C========================================================================
C                    perform transformation
C========================================================================
C
C
C------------------------------------------------------  B = S  * A * S+
      IF ( KEY.EQ.'SAS+' ) THEN
C
         IF ( UNITARY ) THEN
C
            CALL ZGEMM('N','N',N,N,N,C1,S,M,A,M,C0,W,M)
C
            CALL ZGEMM('N','C',N,N,N,C1,W,M,S,M,C0,B,M)
C
         ELSE
C
            CALL ZGEMM('N','C',N,N,N,C1,S,M,A,M,C0,W,M)
C
            CALL ZGEMM('N','C',N,N,N,C1,W,M,S,M,C0,B,M)
C
         END IF
C
C-------------------------------------------------------  B = S+ * A * S
      ELSE IF ( KEY.EQ.'S+AS' ) THEN
C
         IF ( UNITARY ) THEN
C
            CALL ZGEMM('C','N',N,N,N,C1,S,M,A,M,C0,W,M)
C
            CALL ZGEMM('N','N',N,N,N,C1,W,M,S,M,C0,B,M)
C
         ELSE
C
            CALL ZGEMM('C','C',N,N,N,C1,S,M,A,M,C0,W,M)
C
            CALL ZGEMM('N','N',N,N,N,C1,W,M,S,M,C0,B,M)
C
         END IF
C
      ELSE
C
         STOP '<CMATTRANS>  KEY  not allowed '
C
      END IF
C-----------------------------------------------------------------------
C
      END
