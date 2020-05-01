C*==chifmagind.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHIFMAGIND(RHOCHR,RHOSPN,RHOORB)
C   ********************************************************************
C   *                                                                  *
C   * SUBROUTINE TO CALCULATE THE  CHARGE, SPIN  AND  ORBITAL DENSITY  *
C   *                  WITHIN AN ATOMIC CELL                           *
C   *                                                                  *
C   * 07/02/99 HE                                                      *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:DATSET,LDATSET
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_TYPES,ONLY:IMT,NTMAX,LTXT_T,TXT_T,NT
      USE MOD_RMESH,ONLY:R,NRMAX,JRWS,R2DRDI
      USE MOD_CONSTANTS,ONLY:A0_ANG,PI
      IMPLICIT NONE
C*--CHIFMAGIND17
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NTMAXCHK,NQNEUMAX
      PARAMETER (NTMAXCHK=10,NQNEUMAX=200)
C
C Dummy arguments
C
      REAL*8 RHOCHR(NRMAX,NTMAX),RHOORB(NRMAX,NTMAX),RHOSPN(NRMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 ARG
      REAL*8 BES0,BES2,CHKO(NTMAXCHK),CHKQ(NTMAXCHK),CHKS(NTMAXCHK),
     &       CHRMAX,CHRMIN,DQNEU,DX,FAVR(NQNEUMAX),FORB(NQNEUMAX),
     &       FSPN(NQNEUMAX),ORBMAX,ORBMIN,QNEU,QNEUMAX,RHOR2,RINT(NRMAX)
     &       ,RINTO(NRMAX),RINTS(NRMAX),SCAL,SCALCHR,SCALORB,SPNMAX,
     &       SPNMIN,YMAX,YMIN
      COMPLEX*16 CJLZ
      CHARACTER*80 FILNAM
      INTEGER I,IC,IFIL,IM,IQNEU,IR,IRTOP,IT,LFN,LINESTYLE(10),NC,NQNEU
      SAVE CHKO,CHKQ,CHKS
C
C*** End of declarations rewritten by SPAG
C
      DATA LINESTYLE/1,5,2,3,4,1,1,1,1,1/
C
C ----------------------------- account for spin degeneracy for IREL <=1
      IF ( IREL.LE.1 ) THEN
         DO IT = 1,NT
            IM = IMT(IT)
            IRTOP = JRWS(IM)
            DO IR = 1,IRTOP
               RHOCHR(IR,IT) = RHOCHR(IR,IT)/2.0D0
               RHOSPN(IR,IT) = 0.0D0
               RHOORB(IR,IT) = 0.0D0
            END DO
         END DO
      END IF
C
      IF ( NT.GT.NTMAXCHK ) STOP '<CHIFMAGIND> NT > NTMAXCHK'
      DO IT = 1,NT
         IM = IMT(IT)
         IRTOP = JRWS(IM)
         DO IR = 1,IRTOP
            RINT(IR) = RHOCHR(IR,IT)*R2DRDI(IR,IM)
         END DO
         CALL RRADINT(IM,RINT,CHKQ(IT))
         DO IR = 1,IRTOP
            RINT(IR) = RHOSPN(IR,IT)*R2DRDI(IR,IM)
         END DO
         CALL RRADINT(IM,RINT,CHKS(IT))
         DO IR = 1,IRTOP
            RINT(IR) = RHOORB(IR,IT)*R2DRDI(IR,IM)
         END DO
         CALL RRADINT(IM,RINT,CHKO(IT))
      END DO
C
      WRITE (6,99001)
C
      DO IT = 1,NT
C
         IFIL = 72
C
         WRITE (6,'(///)')
C
C=======================================================================
C     PLOT DENSITIES
C=======================================================================
C
         CALL FILNAMT(DATSET,LDATSET,TXT_T,LTXT_T,IT,NT,'dens.agr',8,
     &                FILNAM,LFN,1,IFIL,'DENS-FILE ',10,NTMAX)
C
         IM = IMT(IT)
         IRTOP = JRWS(IM)
C
         CHRMIN = 0D0
         CHRMAX = 0D0
         DO IR = 1,IRTOP
            RINT(IR) = RHOCHR(IR,IT)*R2DRDI(IR,IM)
            RHOR2 = RHOCHR(IR,IT)*R(IR,IM)**2
            CHRMIN = MIN(CHRMIN,RHOR2)
            CHRMAX = MAX(CHRMAX,RHOR2)
         END DO
C
         SPNMIN = +1D20
         SPNMAX = -1D20
         DO IR = 1,IRTOP
            RINT(IR) = RHOSPN(IR,IT)*R2DRDI(IR,IM)
            RHOR2 = RHOSPN(IR,IT)*R(IR,IM)**2
            SPNMIN = MIN(SPNMIN,RHOR2)
            SPNMAX = MAX(SPNMAX,RHOR2)
         END DO
C
         ORBMIN = +1D20
         ORBMAX = -1D20
         DO IR = 1,IRTOP
            RINT(IR) = RHOORB(IR,IT)*R2DRDI(IR,IM)
            RHOR2 = RHOORB(IR,IT)*R(IR,IM)**2
            ORBMIN = MIN(ORBMIN,RHOR2)
            ORBMAX = MAX(ORBMAX,RHOR2)
         END DO
C
         WRITE (IFIL,99003) '#',DATSET(1:LDATSET)
         WRITE (IFIL,99002) '#',IT,TXT_T(IT)(1:LTXT_T(IT))
C
         SCAL = 2D0
         SCALCHR = 0.0D0
         SCALORB = 1D0
         YMIN = MIN(CHRMIN*SCALCHR,SPNMIN,ORBMIN*SCALORB)
         YMAX = MAX(CHRMAX*SCALCHR,SPNMAX,ORBMAX*SCALORB)
         YMIN = DBLE(INT(YMIN/SCAL)-1)*SCAL
         YMAX = DBLE(INT(YMAX/SCAL)+1)*SCAL
         DX = 0.5D0
C
         NC = 3
         DO IC = 1,NC
            WRITE (IFIL,99005) (IC-1),IC,(IC-1),LINESTYLE(IC),(IC-1),2
         END DO
C
         DO IC = 1,NC
            WRITE (IFIL,99004)
            IF ( IC.EQ.1 ) WRITE (IFIL,99007)
     &                            (R(IR,IM),SCALCHR*RHOCHR(IR,IT)
     &                            *R(IR,IM)**2,IR=1,IRTOP)
            IF ( IC.EQ.2 ) WRITE (IFIL,99007)
     &                            (R(IR,IM),RHOSPN(IR,IT)*R(IR,IM)**2,
     &                            IR=1,IRTOP)
            IF ( IC.EQ.3 ) WRITE (IFIL,99007)
     &                            (R(IR,IM),SCALORB*RHOORB(IR,IT)
     &                            *R(IR,IM)**2,IR=1,IRTOP)
            WRITE (IFIL,99006)
         END DO
C
C=======================================================================
C     CALCULATE AND PLOT SCATTERING AMPLITUDES
C=======================================================================
C
         CALL FILNAMT(DATSET,LDATSET,TXT_T,LTXT_T,IT,NT,'chifmag.agr',
     &                11,FILNAM,LFN,1,IFIL,'f_mag-FILE',10,NTMAX)
C
         WRITE (IFIL,99003) '#',DATSET(1:LDATSET)
         WRITE (IFIL,99002) '#',IT,TXT_T(IT)(1:LTXT_T(IT))
C
         NQNEU = MIN(NQNEUMAX,100)
         QNEUMAX = 4*PI*A0_ANG
         DQNEU = QNEUMAX/DBLE(NQNEU-1)
         DO IQNEU = 1,NQNEU
            QNEU = (IQNEU-1)*DQNEU
            DO IR = 1,IRTOP
               ARG = DCMPLX(QNEU*R(IR,IM),0.0D0)
               BES0 = DREAL(CJLZ(0,ARG))
               BES2 = DREAL(CJLZ(2,ARG))
               RINTS(IR) = BES0*RHOSPN(IR,IT)*R2DRDI(IR,IM)
               RINTO(IR) = (BES0+BES2)*RHOORB(IR,IT)*R2DRDI(IR,IM)
            END DO
            CALL RRADINT(IM,RINTS,FSPN(IQNEU))
            CALL RRADINT(IM,RINTO,FORB(IQNEU))
C
            FAVR(IQNEU) = FSPN(IQNEU) + FORB(IQNEU)
         END DO
         DO IQNEU = NQNEU,1, - 1
            FAVR(IQNEU) = FAVR(IQNEU)/FAVR(1)
            FSPN(IQNEU) = FSPN(IQNEU)/FSPN(1)
            FORB(IQNEU) = FORB(IQNEU)/FORB(1)
         END DO
C
         DX = DQNEU/(4*PI*A0_ANG)
C
         NC = 3
         DO IC = 1,NC
            WRITE (IFIL,99005) (IC-1),IC,(IC-1),LINESTYLE(IC),(IC-1),2
         END DO
C
         DO IC = 1,NC
            WRITE (IFIL,99004)
            IF ( IC.EQ.1 ) WRITE (IFIL,99007)
     &                            ((I-1)*DX,FAVR(I),I=1,NQNEU)
            IF ( IC.EQ.2 ) WRITE (IFIL,99007)
     &                            ((I-1)*DX,FSPN(I),I=1,NQNEU)
            IF ( IC.EQ.3 ) WRITE (IFIL,99007)
     &                            ((I-1)*DX,FORB(I),I=1,NQNEU)
            WRITE (IFIL,99006)
         END DO
C
      END DO
C
      CLOSE (72)
      WRITE (6,*) '--> in <CHIFMAGIND>  --- all done '
99001 FORMAT (' ',//,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*            ******       *     *    **     ****             *'
     &  ,/,10X,
     &  '*            *            **   **   *  *   *    *            *'
     &  ,/,10X,
     &  '*            *            * * * *  *    *  *                 *'
     &  ,/,10X,
     &  '*            *****   ***  *  *  *  ******  *  ***            *'
     &  ,/,10X,
     &  '*            *            *     *  *    *  *    *            *'
     &  ,/,10X,
     &  '*            *            *     *  *    *  *    *            *'
     &  ,/,10X,
     &  '*            *            *     *  *    *   ****             *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
C
99002 FORMAT (A1,9X,'IT=',I2,3X,A2)
99003 FORMAT (A1,9X,'DATASET: ',A)
99004 FORMAT ('@TYPE xy')
99005 FORMAT ('@    s',I1,' color ',I2,/,'@    s',I1,' linestyle ',I2,/,
     &        '@    s',I1,' linewidth ',I2)
99006 FORMAT ('&')
99007 FORMAT (2E15.5)
C
      END
C*==filnamt.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE FILNAMT(DATSET,LDATSET,TXT_T,LTXT_T,IT,NT,EXT,LEXT,
     &                   FILNAM,LFN,KOPEN,IFIL,TEXT,LTEXT,NTMAX)
C   ********************************************************************
C   *                                                                  *
C   *   create a filename for atom type   IT                           *
C   *                                                                  *
C   *   [DATASET_]TXT_T(IT)(1:LTXT_T(IT))[_x(IT)]_EXT                  *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--FILNAMT260
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*80 DATSET,FILNAM
      CHARACTER*(*) EXT,TEXT
      INTEGER IFIL,IT,KOPEN,LDATSET,LEXT,LFN,LTEXT,NT,NTMAX
      INTEGER LTXT_T(NTMAX)
      CHARACTER*8 TXT_T(NTMAX)
C
C Local variables
C
      INTEGER IA,JT,NA
C
C*** End of declarations rewritten by SPAG
C
      IF ( LDATSET.NE.0 ) THEN
         FILNAM = DATSET(1:LDATSET)//TXT_T(IT)(1:LTXT_T(IT))
      ELSE
         FILNAM = TXT_T(IT)(1:LTXT_T(IT))
      END IF
      LFN = LDATSET + LTXT_T(IT)
C
      NA = 0
      DO JT = 1,NT
         IF ( TXT_T(IT).EQ.TXT_T(JT) ) THEN
            NA = NA + 1
            IF ( IT.GE.JT ) IA = NA
         END IF
      END DO
C
      IF ( NA.GT.1 ) THEN
C
         FILNAM = FILNAM(1:LFN)//'_'//CHAR(ICHAR('a')-1+IA)//'_'
         CALL STRING_ADD_N(FILNAM,IT)
C
      END IF
C
      LFN = LEN_TRIM(FILNAM)
      FILNAM = FILNAM(1:LFN)//'_'//EXT(1:LEXT)
      LFN = LFN + 1 + LEXT
C
      IF ( KOPEN.NE.0 ) THEN
C
         OPEN (UNIT=IFIL,FILE=FILNAM)
         WRITE (6,'(10X,2A,I2,A,A)') TEXT(1:LTEXT),':  (',IFIL,') ',
     &                               FILNAM(1:LFN)
C
      END IF
C
      END
