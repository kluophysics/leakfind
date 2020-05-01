C*==dvdrspline.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DVDRSPLINE(V,R,DVDR,N)
C   ********************************************************************
C   *                                                                  *
C   *   CALCULATE DERIVATIVE OF THE POTENTIAL                          *
C   *   DV/DR= ( D(R**2*V)/DR-2.R*V ) / R**2                           *
C   *                                                                  *
C   *   calling J.REDINGER's routine <DERSPL>                          *
C   *                                                                  *
C   * 27/10/94  HE  adopted for SPRKKR-package                         *
C   ********************************************************************
      IMPLICIT NONE
C*--DVDRSPLINE13
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N
      REAL*8 DVDR(N),R(N),V(N)
C
C Local variables
C
      REAL*8 DR2VDR(N),R2V(N)
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
C
C
      DO I = 1,N
         R2V(I) = V(I)*R(I)**2
      END DO
C
      CALL DERSPL(N,R,R2V,DR2VDR)
C
      DO I = 1,N
         DVDR(I) = (DR2VDR(I)-2.D0*R2V(I)/R(I))/R(I)**2
      END DO
      END
