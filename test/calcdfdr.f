C*==calcdfdr.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE CALCDFDR(F,FP,DRDI,N)
C   ********************************************************************
C   *                                                                  *
C   *   fp = d f(r)/dr = d f(i)/di / di/dr   with  drdi = dr / di      *
C   *                                                                  *
C   *   using a 5-point formula                                        *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--CALCDFDR11
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N
      REAL*8 DRDI(N),F(N),FP(N)
C
C Local variables
C
      INTEGER I,NM2
C
C*** End of declarations rewritten by SPAG
C
C =============================================================== N >= 5
      IF ( N.GE.5 ) THEN
C
         NM2 = N - 2
C
C     forward difference at the beginning of the table
C
         FP(1) = ((2D0*F(4)+18D0*F(2))-(9D0*F(3)+11D0*F(1)))/6D0
         FP(2) = ((2D0*F(5)+18D0*F(3))-(9D0*F(4)+11D0*F(2)))/6D0
C
C     central difference at the interior of the table
C
         DO I = 3,NM2
            FP(I) = ((F(I-2)+8D0*F(I+1))-(8D0*F(I-1)+F(I+2)))/12D0
         END DO
C
C     backward difference at the end of the table
C
         FP(N) = ((11D0*F(N)+9D0*F(N-2))-(18D0*F(N-1)+2D0*F(N-3)))/6D0
         FP(N-1) = ((11D0*F(N-1)+9D0*F(N-3))-(18D0*F(N-2)+2D0*F(N-4)))
     &             /6D0
C
         DO I = 1,N
            FP(I) = FP(I)/DRDI(I)
         END DO
C
C ================================================================ N = 4
      ELSE IF ( N.EQ.4 ) THEN
C
         FP(1) = (2D0*F(4)-9D0*F(3)+18D0*F(2)-11D0*F(1))/(6D0*DRDI(1))
         FP(2) = (-F(4)+6D0*F(3)-3D0*F(2)-2D0*F(1))/(6D0*DRDI(2))
         FP(3) = (2D0*F(4)+3D0*F(3)-6D0*F(2)+F(1))/(6D0*DRDI(3))
         FP(4) = (11D0*F(4)-18D0*F(3)+9D0*F(2)-2D0*F(1))/(6D0*DRDI(4))
C
C ================================================================ N = 3
      ELSE IF ( N.EQ.3 ) THEN
C
         FP(1) = (-F(3)+4D0*F(2)-3D0*F(1))/(2D0*DRDI(1))
         FP(2) = (F(3)-F(1))/(2D0*DRDI(2))
         FP(3) = (+3D0*F(3)-4D0*F(2)+F(1))/(2D0*DRDI(3))
C
C ================================================================ N = 2
      ELSE IF ( N.EQ.2 ) THEN
C
         FP(1) = (F(2)-F(1))/DRDI(1)
         FP(2) = FP(1)
C
C ================================================================ N = 1
      ELSE IF ( N.EQ.1 ) THEN
C
         STOP '<CALCDFDR> N=1 !!!!!'
C
      END IF
C
      END
C*==calc_dcfdr.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE CALC_DCFDR(CF,CFP,DRDI,N)
C   ********************************************************************
C   *                                                                  *
C   *   fp = d f(r)/dr = d f(i)/di / di/dr   with  drdi = dr / di      *
C   *                                                                  *
C   *   using a 5-point formula                                        *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--CALC_DCFDR102
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N
      COMPLEX*16 CF(N),CFP(N)
      REAL*8 DRDI(N)
C
C Local variables
C
      INTEGER I,NM2
C
C*** End of declarations rewritten by SPAG
C
C =============================================================== N >= 5
      IF ( N.GE.5 ) THEN
C
         NM2 = N - 2
C
C     forward difference at the beginning of the table
C
         CFP(1) = ((2D0*CF(4)+18D0*CF(2))-(9D0*CF(3)+11D0*CF(1)))/6D0
         CFP(2) = ((2D0*CF(5)+18D0*CF(3))-(9D0*CF(4)+11D0*CF(2)))/6D0
C
C     central difference at the interior of the table
C
         DO I = 3,NM2
            CFP(I) = ((CF(I-2)+8D0*CF(I+1))-(8D0*CF(I-1)+CF(I+2)))/12D0
         END DO
C
C     backward difference at the end of the table
C
         CFP(N) = ((11D0*CF(N)+9D0*CF(N-2))-(18D0*CF(N-1)+2D0*CF(N-3)))
     &            /6D0
         CFP(N-1) = ((11D0*CF(N-1)+9D0*CF(N-3))-(18D0*CF(N-2)+2D0*CF(N-4
     &              )))/6D0
C
         DO I = 1,N
            CFP(I) = CFP(I)/DRDI(I)
         END DO
C
C ================================================================ N = 4
      ELSE IF ( N.EQ.4 ) THEN
C
         CFP(1) = (2D0*CF(4)-9D0*CF(3)+18D0*CF(2)-11D0*CF(1))
     &            /(6D0*DRDI(1))
         CFP(2) = (-CF(4)+6D0*CF(3)-3D0*CF(2)-2D0*CF(1))/(6D0*DRDI(2))
         CFP(3) = (2D0*CF(4)+3D0*CF(3)-6D0*CF(2)+CF(1))/(6D0*DRDI(3))
         CFP(4) = (11D0*CF(4)-18D0*CF(3)+9D0*CF(2)-2D0*CF(1))
     &            /(6D0*DRDI(4))
C
C ================================================================ N = 3
      ELSE IF ( N.EQ.3 ) THEN
C
         CFP(1) = (-CF(3)+4D0*CF(2)-3D0*CF(1))/(2D0*DRDI(1))
         CFP(2) = (CF(3)-CF(1))/(2D0*DRDI(2))
         CFP(3) = (+3D0*CF(3)-4D0*CF(2)+CF(1))/(2D0*DRDI(3))
C
C ================================================================ N = 2
      ELSE IF ( N.EQ.2 ) THEN
C
         CFP(1) = (CF(2)-CF(1))/DRDI(1)
         CFP(2) = CFP(1)
C
C ================================================================ N = 1
      ELSE IF ( N.EQ.1 ) THEN
C
         STOP '<CALC_DCFDR> N=1 !!!!!'
C
      END IF
C
      END
