C*==dmft_padeoz.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_PADEOZ(OMEGA,SIGMAO,EFG,NZM,ZM,SIGMAZ,NOM)
C---- Pade continuations from Matzubara to semicircular
C
C
      USE MOD_CONSTANTS,ONLY:CI
      IMPLICIT NONE
C*--DMFT_PADEOZ8
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 EFG
      INTEGER NOM,NZM
      REAL*8 OMEGA(NOM)
      COMPLEX*16 SIGMAO(NOM),SIGMAZ(NZM),ZM(NZM)
C
C Local variables
C
      logical lnan
      COMPLEX*16 EZ,GE,PZ(:,:),ZOM(:)
      INTEGER IOM,LZ,NOMI
C
C*** End of declarations rewritten by SPAG
C
C
C
C Dummy arguments
C
C
C Local variables
C
C
      ALLOCATABLE PZ,ZOM
      NOMI = NOM / 2
      ALLOCATE (PZ(NOMI,NOMI),ZOM(NOM))
C
C*** End of declarations rewritten by SPAG
C
      ZOM = (0.0D0,0.0D0)
      PZ = (0.0D0,0.0D0)
      DO IOM = 1,NOM/2
         ZOM(IOM) = CI*OMEGA(2*IOM) + EFG
      END DO
      CALL DMFT_PADECOF(SIGMAO,ZOM,PZ,NOM,NOMI,lnan)
      if(lnan) then
         print *,'in padeoz, mats-->epath, a pade coff is NaN'
         stop
      endif
      DO LZ = 1,NZM
         EZ = ZM(LZ)          ! new
         CALL DMFT_GPADE(GE,EZ,ZOM,PZ,NOM,NOMI)
         SIGMAZ(LZ) = GE
      END DO
      END
C*==dmft_padezo.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_PADEZO(EFG,NZM,ZM,GZ,OMEGA,MINOM,GO,NOM)
C
C---- Pade continuations of Green function from semicircular to Matzubara
C
C
      USE MOD_CONSTANTS,ONLY:CI
      IMPLICIT NONE
C*--DMFT_PADEZO71
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 EFG
      INTEGER NOM,NZM
      COMPLEX*16 GO(NOM),GZ(NZM),ZM(NZM)
      INTEGER MINOM(NOM)
      REAL*8 OMEGA(NOM)
C
C Local variables
C
      logical lnan
      COMPLEX*16 EZ,GE,P(:,:)
      INTEGER IOM,MIOM
C
C*** End of declarations rewritten by SPAG
C
C
C
C Dummy arguments
C
C
C Local variables
C
C
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE P
      ALLOCATE (P(NZM,NZM))
C
C
C
      CALL DMFT_PADECOF(GZ,ZM,P,NZM,NZM,lnan)
      if(lnan) then
        print *,'in padezo, epath-->mats, a pade coff is NaN'
        stop
      endif
C
      DO IOM = 1,NOM/4
         MIOM = MINOM(2*IOM)
         EZ = CI*OMEGA(2*IOM) + EFG
                                 ! new
C         ez=-ci*omega(2*iom)+efg  ! new
         CALL DMFT_GPADE(GE,EZ,ZM,P,NZM,NZM)
C          go(miom)=ge
         GO(MIOM) = DCONJG(GE)
C     go(2*iom)=DCONJG(ge)
         GO(2*IOM) = GE
      END DO
C
      END
C*==dmft_padecof.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C     *********************
      SUBROUTINE DMFT_PADECOF(G,Z,P,NOM,N,lnan)
C     *********************
C     Recursion for Pade coefficient (J.Serene)
C
      IMPLICIT NONE
C*--DMFT_PADECOF140
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      logical lnan
      INTEGER N,NOM
      COMPLEX*16 G(NOM),P(N,N),Z(NOM)
      complex*16 tmp1,tmp2
      real*8 pr, pi
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
C
C
C
      DO J = 1,N
         P(1,J) = G(J)
      END DO
C
      DO J = 2,N
         DO I = 2,J
c            P(I,J) = (P(I-1,I-1)-P(I-1,J))/(Z(J)-Z(I-1))
c     &               /(P(I-1,J)+1D-18)
             tmp1 = p(i-1,i-1) / p(i-1,j)
             tmp2 = p(i-1,j) / p(i-1,j)
             p(i,j) = ( tmp1 - tmp2 ) / ( z(j) - z(i-1) )
         END DO
      END DO
      lnan = .false.
      do j=1,n
        pr = real(p(j,j))
        pi = aimag(p(j,j))
        if( (pr .ne. pr) .or. (pi .ne. pi) ) then
          lnan = .true.
          return
        endif
      enddo
      END
C*==dmft_gpade.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
C
      SUBROUTINE DMFT_GPADE(GE,E,Z,P,NZM,N)
C
C     *********************
C
C     Calculation of a Green's function for a given pade-coeff-p(i,j)
C     on the real axis e=e+i0
C
C
      USE MOD_CONSTANTS,ONLY:C0,C1
      IMPLICIT NONE
C*--DMFT_GPADE194
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 E,GE
      INTEGER N,NZM
      COMPLEX*16 P(N,N),Z(NZM)
C
C Local variables
C
      COMPLEX*16 A(:),B(:)
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
C
C
C Dummy arguments
C
C
C Local variables
C
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE A,B
C
      ALLOCATE (A(0:N),B(0:N))
C
C       print*, 'in gpade', 'nzm=',nzm
C       print*, 'in gpade', 'e=',e
C
      A(0) = C0
      A(1) = P(1,1)
      B(0) = C1
      B(1) = C1
      DO I = 1,N - 1
         A(I+1) = A(I) + (E-Z(I))*P(I+1,I+1)*A(I-1)
         B(I+1) = B(I) + (E-Z(I))*P(I+1,I+1)*B(I-1)
      END DO
      GE = A(N)/B(N)
      END
C***************************************************
