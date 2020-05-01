C*==ametimrev.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE AMETIMREV(A1,A2,KEY,NKMMAX)
C   ********************************************************************
C   *                                                                  *
C   *  perform time reversal on the angular matrix elements            *
C   *  for the various forms of the interaction operator               *
C   *  KEY = 1: deal with   <  kap,mue |X| +/- kap',mue'>              *
C   *        2: deal with   < -kap,mue |X| +/- kap',mue'>              *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--AMETIMREV12
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER KEY,NKMMAX
      REAL*8 A1(NKMMAX,NKMMAX,3),A2(NKMMAX,NKMMAX,3)
C
C Local variables
C
      REAL*8 F1,F2,J,MJ,SK,W1(:,:),W2(:,:)
      INTEGER IA_ERR,IPOL,L,LAM,LAMP,LAMT,NLAM,NPOL
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE W1,W2
C
      ALLOCATE (W1(NKMMAX,NKMMAX),W2(NKMMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:ametimrev -> W2'
C
      NPOL = 3
      NLAM = NKMMAX
C
      CALL RINIT(NKMMAX*NKMMAX,W1)
      CALL RINIT(NKMMAX*NKMMAX,W2)
C
      DO IPOL = 1,NPOL
         DO LAM = 1,NLAM
            L = INT(SQRT(DBLE(LAM)/2.0D0)-0.0001D0)
            J = L + SIGN(0.5D0,DBLE(LAM-(2*L*(L+1))-0.0001D0))
            MJ = DBLE(LAM) - (L*2*(J+0.5D0)+J+1.0D0)
            SK = 2*L - NINT(2*J)
            LAMT = NINT((J+0.5D0)*2*L+J-MJ+1)
C
            F1 = (-1D0)**NINT(MJ+SK/2.0D0)
            IF ( KEY.EQ.1 ) THEN
               F2 = F1
            ELSE
               F2 = (-1D0)**NINT(MJ-SK/2.0D0)
            END IF
C
            DO LAMP = 1,NLAM
               W1(LAM,LAMP) = F1*A1(LAMT,LAMP,IPOL)
               W2(LAM,LAMP) = F2*A2(LAMT,LAMP,IPOL)
            END DO
C
         END DO
C
         DO LAM = 1,NLAM
            DO LAMP = 1,NLAM
               A1(LAM,LAMP,IPOL) = W1(LAM,LAMP)
               A2(LAM,LAMP,IPOL) = W2(LAM,LAMP)
            END DO
         END DO
C
      END DO
C
      DEALLOCATE (W1,W2,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'dealloc:ametimrev -> W2'
C
      END
