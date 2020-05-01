C*==transpho.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE TRANSPHO(P,TSSTM,NKMMAX,NLMAX,TSSTJB)
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NKMMAX,NLMAX
      COMPLEX*16 P
      COMPLEX*16 TSSTJB(NKMMAX,NKMMAX),TSSTM(NKMMAX,NKMMAX)
C
C Local variables
C
      INTEGER I,I11,J,J11,KAP(:),MUD(:),NB(:),NR(:),NRL(:),NRM(:),
     &        REL(:,:),SB(:)
      REAL*8 NRS(:)
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE NB,SB,NR,KAP,REL,MUD,NRL,NRM,NRS
      ALLOCATE (NB(NKMMAX),SB(NKMMAX),NR(NKMMAX),KAP(NKMMAX))
      ALLOCATE (REL(NKMMAX,2),MUD(NKMMAX),NRL(NKMMAX),NRM(NKMMAX))
      ALLOCATE (NRS(NKMMAX))
C
C*** End of declarations rewritten by SPAG
C
      CALL TRANS_HEJB1(NLMAX,NKMMAX,KAP,MUD,NRL,NRM,NRS,NR,REL)
      CALL TRANS_HEJB2(NB,SB,NKMMAX,NLMAX)
C
C     WRITE (12,*) 'SCATTERING MATRIX AT ENERGY',E
      DO I = 1,NKMMAX
         I11 = REL(NB(I),SB(I))
         DO J = 1,NKMMAX
            J11 = REL(NB(J),SB(J))
            TSSTJB(I11,J11) = DCMPLX(0.D0,-1.D0)*TSSTM(I,J)*P
C           IF (CDABS(TSSTJB(I11,J11)).GT.0.D0) THEN
C               WRITE (12,100) I,J,I11,J11,TSSTJB(I11,J11)
C           ENDIF
         END DO
      END DO
C
C 100 FORMAT (4I5,2X,2D15.8)
      END
C*==transphofak.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C
      SUBROUTINE TRANSPHOFAK(P,TSSTM,NKMMAX,TSSTJB)
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NKMMAX
      COMPLEX*16 P
      COMPLEX*16 TSSTJB(NKMMAX,NKMMAX),TSSTM(NKMMAX,NKMMAX)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,NKMMAX
         DO J = 1,NKMMAX
            TSSTJB(I,J) = DCMPLX(0.D0,-1.D0)*TSSTM(I,J)*P
C           IF (CDABS(TSSTJB(I,J)).GT.1.D-10) THEN
C               WRITE (12,100) I,J,TSSTJB(I,J)
C           ENDIF
         END DO
      END DO
C
      END
C*==transphocpa.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C
C
      SUBROUTINE TRANSPHOCPA(TSSTM,NKMMAX,NLMAX,TSSTJB)
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NKMMAX,NLMAX
      COMPLEX*16 TSSTJB(NKMMAX,NKMMAX),TSSTM(NKMMAX,NKMMAX)
C
C Local variables
C
      INTEGER I,I11,J,J11,KAP(:),MUD(:),NB(:),NR(:),NRL(:),NRM(:),
     &        REL(:,:),SB(:)
      REAL*8 NRS(:)
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE NB,SB,NR,KAP,REL,MUD,NRL,NRM,NRS
      ALLOCATE (NB(NKMMAX),SB(NKMMAX),NR(NKMMAX),KAP(NKMMAX))
      ALLOCATE (REL(NKMMAX,2),MUD(NKMMAX),NRL(NKMMAX),NRM(NKMMAX))
      ALLOCATE (NRS(NKMMAX))
C
C*** End of declarations rewritten by SPAG
C
C
C
C
C
C
      CALL TRANS_HEJB1(NLMAX,NKMMAX,KAP,MUD,NRL,NRM,NRS,NR,REL)
      CALL TRANS_HEJB2(NB,SB,NKMMAX,NLMAX)
C
C     WRITE (12,*) 'CPA PROJECTOR MATRIX AT ENERGY',E
      DO I = 1,NKMMAX
         I11 = REL(NB(I),SB(I))
         DO J = 1,NKMMAX
            J11 = REL(NB(J),SB(J))
            TSSTJB(I11,J11) = TSSTM(I,J)
C           IF (CDABS(TSSTJB(I11,J11)).GT.1.D-10) THEN
C               WRITE (12,100) I,J,I11,J11,TSSTJB(I11,J11)
C           ENDIF
         END DO
      END DO
C
      END
C*==transphow.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE TRANSPHOW(RG,RF,HG,HF,NKMMAX,NLMAX,FF_S1M,FF_G1M,
     &                     FFI_S1M,FFI_G1M,RG1,RF1,FF_S1MZ,FF_G1MZ,
     &                     NZERO)
C
C
      USE MOD_RMESH,ONLY:NRMAX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NKMMAX,NLMAX
      COMPLEX*16 FFI_G1M(NRMAX,NKMMAX,NKMMAX),
     &           FFI_S1M(NRMAX,NKMMAX,NKMMAX),
     &           FF_G1M(NRMAX,NKMMAX,NKMMAX),
     &           FF_G1MZ(NRMAX,NKMMAX,NKMMAX),
     &           FF_S1M(NRMAX,NKMMAX,NKMMAX),
     &           FF_S1MZ(NRMAX,NKMMAX,NKMMAX),HF(NRMAX,NKMMAX,NKMMAX),
     &           HG(NRMAX,NKMMAX,NKMMAX),RF(NRMAX,NKMMAX,NKMMAX),
     &           RF1(NRMAX,NKMMAX,NKMMAX),RG(NRMAX,NKMMAX,NKMMAX),
     &           RG1(NRMAX,NKMMAX,NKMMAX)
      REAL*8 NZERO(NKMMAX,NKMMAX)
C
C Local variables
C
      INTEGER I,I11,J,J11,KAP(:),L,MUD(:),NB(:),NR(:),NRL(:),NRM(:),
     &        REL(:,:),SB(:)
      REAL*8 NRS(:),NZERO2(:,:)
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE NB,SB,NR,KAP,REL,MUD,NRL,NRM,NRS,NZERO2
      ALLOCATE (NB(NKMMAX),SB(NKMMAX),NR(NKMMAX),KAP(NKMMAX))
      ALLOCATE (REL(NKMMAX,2),MUD(NKMMAX),NRL(NKMMAX),NRM(NKMMAX))
      ALLOCATE (NRS(NKMMAX),NZERO2(NKMMAX,NKMMAX))
C
C*** End of declarations rewritten by SPAG
C
C
C
C
C
      NZERO2 = 0.0D0
      CALL TRANS_HEJB1(NLMAX,NKMMAX,KAP,MUD,NRL,NRM,NRS,NR,REL)
      CALL TRANS_HEJB2(NB,SB,NKMMAX,NLMAX)
C
      DO I = 1,NKMMAX
         I11 = REL(NB(I),SB(I))
         DO J = 1,NKMMAX
            J11 = REL(NB(J),SB(J))
            NZERO2(I11,J11) = NZERO(I,J)
            DO L = 1,NRMAX
               FF_S1M(L,I11,J11) = RG(L,I,J)
               FF_G1M(L,I11,J11) = RF(L,I,J)
               FF_S1MZ(L,I11,J11) = RG1(L,I,J)
               FF_G1MZ(L,I11,J11) = RF1(L,I,J)
               FFI_S1M(L,I11,J11) = HG(L,I,J)
               FFI_G1M(L,I11,J11) = HF(L,I,J)
            END DO
         END DO
      END DO
      NZERO = NZERO2
C
      END
