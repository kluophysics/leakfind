C*==rotate.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE ROTATE(T1,MODE,T2,N,ROT,M)
C   ********************************************************************
C   *                                                                  *
C   *   performs the rotation of the matrix  T1  using the rotation-   *
C   *   matrix  ROT, set up by <ROTMAT>                                *
C   *                                                                  *
C   *          T2 = ROT  * T1 * ROT+     IF  MODE = 'L->G'             *
C   *          T2 = ROT+ * T1 * ROT      IF  MODE = 'G->L'             *
C   *                                                                  *
C   *   see:     E.M. ROSE  ELEMENTARY THEORY OF ANGULAR MOMENTUM      *
C   *                                                                  *
C   *   NOTE: T1 is left unchanged                                     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,WKM_LPK
      USE MOD_CONSTANTS,ONLY:C1,C0
      IMPLICIT NONE
C*--ROTATE20
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER M,N
      CHARACTER*4 MODE
      COMPLEX*16 ROT(M,M),T1(M,M),T2(M,M)
C
C Local variables
C
      CHARACTER*1 FL1,FL2
C
C*** End of declarations rewritten by SPAG
C
C========================================================================
C                         check arguments
C========================================================================
C
      IF ( N.NE.NKM ) STOP '<ROTATE>  N <> NKM'
      IF ( M.NE.NKMMAX ) STOP '<ROTATE>  M <> NKMMAX'
C
C========================================================================
C                    perform transformation
C========================================================================
C
      IF ( MODE.EQ.'L->G' ) THEN
         FL1 = 'N'
         FL2 = 'C'
      ELSE IF ( MODE.EQ.'G->L' ) THEN
         FL1 = 'C'
         FL2 = 'N'
      ELSE
         WRITE (6,*) ' MODE = ',MODE
         STOP 'in <ROTATE>  MODE not allowed'
      END IF
C
      CALL ZGEMM(FL1,'N',N,N,N,C1,ROT,M,T1,M,C0,WKM_LPK,M)
C
      CALL ZGEMM('N',FL2,N,N,N,C1,WKM_LPK,M,ROT,M,C0,T2,M)
C
      END
