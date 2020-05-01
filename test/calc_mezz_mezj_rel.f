C*==calc_mezz_mezj_rel.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE CALC_MEZZ_MEZJ_REL(IT,ZFLB,ZGLB,JGLB,JFLB,ZGRA,ZFRA,
     &                              JFRA,JGRA,MREG,MIRR,KIRR,KSYMZJ,NME)
C   ********************************************************************
C   *                                                                  *
C   *   calculate the matrix elements IME up to NME                    *
C   *   using the angular matrix elements set up in CALC_AME           *
C   *                                                                  *
C   *  IME                                                             *
C   *   1: IDOS  < LAM |      1      | LAM' >    (l,s)-resolved DOS    *
C   *   2: ISMT  < LAM | sigma(ipol) | LAM' >    spin moment           *
C   *   3: IOMT  < LAM |     l(ipol) | LAM' >    orbital moment        *
C   *   4: IHFF  < LAM |  B_hf(ipol) | LAM' >    hyperfine field       *
C   *   5: ISDM  < LAM |     T(ipol) | LAM' >    spin dipole moment    *
C   *   6: IKDS              (dummy)             kappa-resolved DOS    *
C   *   7: IODN  < LAM |P_dn l(ipol) | LAM' >    orb. mnt. spin down   *
C   *   8: IOUP  < LAM |P_up l(ipol) | LAM' >    orb. mnt. spin up     *
C   *   9: IALF  < LAM | alpha(ipol) | LAM' >    velocity operator     *
C   *                                                                  *
C   *   IPOL= 1,2,3  ==  -1,0,+1  ==  (-),(z),(+)                      *
C   *                                                                  *
C   *   IKM = 2 * l * (j+1/2) + j + mj + 1                             *
C   *                                                                  *
C   *   KIRR = 0: only matrix elements of the type MREG = <ZLB|O|ZRA>  *
C   *   KIRR = 1: also matrix elements of the type MIRR = <ZLB|O|JRA>  *
C   *                                                                  *
C   *   KSYMZJ = 1:  use symmetrized value   (ZJ+JZ)/2                 *
C   *                                                                  *
C   *   *RA:    RHS wave functions for energy and/or atom type A       *
C   *   *LB:    LHS wave functions for energy and/or atom type B       *
C   *                                                                  *
C   ********************************************************************
      USE MOD_RMESH,ONLY:R2DRDI,NRMAX,JRWS
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX,NKM,NCPLWF,AME_G,AME_F,MEZZ,MEZJ
      USE MOD_TYPES,ONLY:IMT,NCPLWFMAX,IKMCPLWF
      IMPLICIT NONE
C*--CALC_MEZZ_MEZJ_REL38
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      LOGICAL CHECK_ME
      PARAMETER (CHECK_ME=.TRUE.)
      REAL*8 RELTOL,THRESHOLD
      PARAMETER (RELTOL=1D-12,THRESHOLD=1D-12)
C
C Dummy arguments
C
      INTEGER IT,KIRR,KSYMZJ,NME
      COMPLEX*16 JFLB(NRMAX,NCPLWFMAX,NKM),JFRA(NRMAX,NCPLWFMAX,NKM),
     &           JGLB(NRMAX,NCPLWFMAX,NKM),JGRA(NRMAX,NCPLWFMAX,NKM),
     &           MIRR(NKMMAX,NKMMAX,3,NMEMAX),
     &           MREG(NKMMAX,NKMMAX,3,NMEMAX),ZFLB(NRMAX,NCPLWFMAX,NKM),
     &           ZFRA(NRMAX,NCPLWFMAX,NKM),ZGLB(NRMAX,NCPLWFMAX,NKM),
     &           ZGRA(NRMAX,NCPLWFMAX,NKM)
C
C Local variables
C
      COMPLEX*16 ADD_F,ADD_G,AIRR,AREG,DIRR,DREG,JFZF(2,2),JGZG(2,2),
     &           ZFJF(2,2),ZFZF(2,2),ZGJG(2,2),ZGZG(2,2)
      REAL*8 FSGN_ME(:),RDIFF
      INTEGER I,IA,IB,IKMA,IKMB,IM,IME,IPOL,IRTOP,J,LAMA,LAMB,NPOL
      CHARACTER*4 TXT_ME(12)
C
C*** End of declarations rewritten by SPAG
C
      DATA (TXT_ME(I),I=1,3)/'CHR ','SPN ','ORB '/
C
      ALLOCATABLE FSGN_ME
C
      ALLOCATE (FSGN_ME(NMEMAX))
C
C-----------------------------------------------------------------------
C       account for beta when calculating spin and orbital moments
C-----------------------------------------------------------------------
      FSGN_ME(:) = 1D0
      IF ( NMEMAX.GT.1 ) THEN
         DO IME = 2,MIN(3,NMEMAX)
            FSGN_ME(IME) = -1D0
         END DO
      END IF
C-----------------------------------------------------------------------
C
      NPOL = 3
C
      IM = IMT(IT)
      IRTOP = JRWS(IM)
C
      MREG(:,:,:,:) = C0
      MIRR(:,:,:,:) = C0
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                                   LAMA
      DO LAMA = 1,NKM
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C                                                                   LAMB
         DO LAMB = 1,NKM
C
C-----------------------------------------------------------------------
C                     check the selection rules
C-----------------------------------------------------------------------
C
            DO IME = 1,NME
               DO IA = 1,NCPLWF(LAMA)
                  IKMA = IKMCPLWF(IA,LAMA)
                  DO IB = 1,NCPLWF(LAMB)
                     IKMB = IKMCPLWF(IB,LAMB)
                     DO IPOL = 1,NPOL
                        IF ( ABS(AME_G(IKMB,IKMA,IPOL,IME)).GT.1D-8 )
     &                       GOTO 20
                     END DO
                  END DO
               END DO
            END DO
C -------------------------------------- all angular matrix elements = 0
            CYCLE
C ---------------------------------- non-0 angular matrix elements found
C
C
C-----------------------------------------------------------------------
C                   calculate radial matrix elements
C-----------------------------------------------------------------------
C
 20         CONTINUE
            CALL CINTABR(ZGLB(1,1,LAMB),ZGRA(1,1,LAMA),ZGZG,
     &                   ZFLB(1,1,LAMB),ZFRA(1,1,LAMA),ZFZF,R2DRDI(1,IM)
     &                   ,NCPLWF(LAMB),NCPLWF(LAMA),IRTOP,NRMAX)
C
C-----------------------------------------------------------------------
C                calculate total matrix elements  MREG
C-----------------------------------------------------------------------
C
            DO IA = 1,NCPLWF(LAMA)
               IKMA = IKMCPLWF(IA,LAMA)
               DO IB = 1,NCPLWF(LAMB)
                  IKMB = IKMCPLWF(IB,LAMB)
C
                  DO IME = 1,NME
                     DO IPOL = 1,NPOL
C
                        ADD_G = AME_G(IKMB,IKMA,IPOL,IME)*ZGZG(IB,IA)
                        ADD_F = AME_F(IKMB,IKMA,IPOL,IME)*ZFZF(IB,IA)
C
                        MREG(LAMB,LAMA,IPOL,IME)
     &                     = MREG(LAMB,LAMA,IPOL,IME) + ADD_G + 
     &                     FSGN_ME(IME)*ADD_F
C
                     END DO
                  END DO
C
               END DO
            END DO
C
C=======================================================================
C                    deal with IRREGULAR terms
C=======================================================================
C
            IF ( LAMB.EQ.LAMA .AND. KIRR.EQ.1 ) THEN
C
C-----------------------------------------------------------------------
C                   calculate radial matrix elements
C-----------------------------------------------------------------------
C
               CALL CINTABR(ZGLB(1,1,LAMB),JGRA(1,1,LAMA),ZGJG,
     &                      ZFLB(1,1,LAMB),JFRA(1,1,LAMA),ZFJF,
     &                      R2DRDI(1,IM),NCPLWF(LAMB),NCPLWF(LAMA),
     &                      IRTOP,NRMAX)
C
               IF ( KSYMZJ.EQ.1 ) THEN
C
                  CALL CINTABR(JGLB(1,1,LAMB),ZGRA(1,1,LAMA),JGZG,
     &                         JFLB(1,1,LAMB),ZFRA(1,1,LAMA),JFZF,
     &                         R2DRDI(1,IM),NCPLWF(LAMB),NCPLWF(LAMA),
     &                         IRTOP,NRMAX)
C
                  DO IA = 1,NCPLWF(LAMA)
                     DO IB = 1,NCPLWF(LAMB)
                        ZGJG(IB,IA) = (ZGJG(IB,IA)-JGZG(IB,IA))/2D0
                        ZFJF(IB,IA) = (ZFJF(IB,IA)-JFZF(IB,IA))/2D0
                     END DO
                  END DO
C
               END IF
C
C-----------------------------------------------------------------------
C                calculate total matrix elements  MIRR
C-----------------------------------------------------------------------
C
               DO IA = 1,NCPLWF(LAMA)
                  IKMA = IKMCPLWF(IA,LAMA)
                  DO IB = 1,NCPLWF(LAMB)
                     IKMB = IKMCPLWF(IB,LAMB)
C
                     DO IME = 1,NME
                        DO IPOL = 1,NPOL
C
                           ADD_G = AME_G(IKMB,IKMA,IPOL,IME)*ZGJG(IB,IA)
                           ADD_F = AME_F(IKMB,IKMA,IPOL,IME)*ZFJF(IB,IA)
C
                           MIRR(LAMB,LAMA,IPOL,IME)
     &                        = MIRR(LAMB,LAMA,IPOL,IME) + ADD_G + 
     &                        FSGN_ME(IME)*ADD_F
C
                        END DO
                     END DO
C
                  END DO
               END DO
C
            END IF
C
         END DO
C                                                                   LAMB
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
      END DO
C                                                                   LAMA
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C              check results against standard matrix elements
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( CHECK_ME ) THEN
         IPOL = 2
         DO I = 1,NKM
            DO J = 1,NKM
               DO IME = 1,3
C
                  AREG = MEZZ(I,J,IT,IME) + MREG(I,J,IPOL,IME)
                  DREG = MEZZ(I,J,IT,IME) - MREG(I,J,IPOL,IME)
                  IF ( ABS(AREG).GT.THRESHOLD ) THEN
                     RDIFF = ABS(DREG/AREG)
                     IF ( RDIFF.GT.RELTOL ) WRITE (6,99001) TXT_ME(IME),
     &                    '  REGULAR',I,J,MEZZ(I,J,IT,IME),
     &                    MREG(I,J,IPOL,IME),DREG,RDIFF
                  END IF
C
                  AIRR = MEZJ(I,J,IT,IME) + MIRR(I,J,IPOL,IME)
                  DIRR = MEZJ(I,J,IT,IME) - MIRR(I,J,IPOL,IME)
                  IF ( ABS(AIRR).GT.THRESHOLD ) THEN
                     RDIFF = ABS(DIRR/AIRR)
                     IF ( RDIFF.GT.RELTOL ) WRITE (6,99001) TXT_ME(IME),
     &                    '  IRREGULAR',I,J,MEZJ(I,J,IT,IME),
     &                    MIRR(I,J,IPOL,IME),DIRR,RDIFF
                  END IF
C
               END DO
            END DO
         END DO
      END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
99001 FORMAT (/,79('*'),/,10X,'conflict for  ',A,2X,A,
     &        '  matrix element',/,2I3,2E17.8,/,6X,2E17.8,2X,2E17.8,/,
     &        6X,2E17.8)
C
      END
