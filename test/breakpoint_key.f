C*==stop_breakpoint.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE STOP_BREAKPOINT(ROUTINE)
C   ********************************************************************
C   *                                                                  *
C   *   setting a BREAKPOINT by setting  BREAK = X  in section CONTROL *
C   *   forces additional output and stop of the program               *
C   *   see CASE list below for the various BREAKPONTs                 *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_MPI,ONLY:MPI,NPROCS,MPI_ID
      USE MOD_CALCMODE,ONLY:BREAKPOINT
      IMPLICIT NONE
C*--STOP_BREAKPOINT14
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*(*) ROUTINE
C
C Local variables
C
      INTEGER IERR,LR
      CHARACTER*74 T
C
C*** End of declarations rewritten by SPAG
C
      IF ( MPI ) WRITE (6,99002) MPI_ID,NPROCS
C
      LR = LEN_TRIM(ROUTINE)
C
      SELECT CASE (BREAKPOINT)
C------------------------------------------------------------------- (1)
C                                                              <SFNDUMP>
      CASE (1)
C
         T = 'shape functions dumped by calling <SFNDUMP> in <SCF>'
C
C------------------------------------------------------------------- (2)
C                                                              <SCFCMNT>
      CASE (2)
C
         T = 'charge moments written to file  QLM_charge_moment.dat'
C
C------------------------------------------------------------------- (3)
C                                                           <SCFMAD_POT>
      CASE (3)
C
         T = 'potential info written to VLM_check_sum.dat and V_madel'
C
C------------------------------------------------------------------- (4)
C                                                               <*SSITE>
      CASE (4)
C
         T = 'single site info dumped'
C
C------------------------------------------------------------------- (5)
C
C      - <SCF>:            write EFERMI, V, B to scratch file
C      - <FPSSITE>:        read rot. angles TETNEW,PHINEW
C                          write rotated single site matrics TSST etc.
C      - <FPSCFNEWPOT>:    read rot. angles TETNEW,PHINEW
C                          rotate CMNT
C                          reread original EFERMI,V,B from scratch file
C                          rotate V, B
C                          write to BREAK_V_rotated.dat
C
      CASE (5)
C
         T = 'potential rotated and dumped to BREAK_V_rotated.dat'
C
C------------------------------------------------------------------- (6)
C
C      - <SCF>:            reread ROTATED EFERMI, V, B
C                          from BREAK_V_rotated.dat
C      - <FPSSITE>:        write single site matrics TSST etc.
C      - <FPSCFNEWPOT>:    print results and stop
C
      CASE (6)
C
         T = 'rotated potential from  BREAK_V_rotated.dat used '
C
C-----------------------------------------------------------------------
C                                                              <>
      CASE (-1)
C
         T = ''
C
C-----------------------------------------------------------------------
      CASE DEFAULT
         T = ' '
      END SELECT
C-----------------------------------------------------------------------
C
      WRITE (6,99001) BREAKPOINT,ROUTINE(1:LR),T
C
      IF ( MPI ) CALL MPI_FINALIZE(IERR)
C
      STOP
C
99001 FORMAT (2(/,1X,79('B')),//,10X,
     &        'test run completed at BREAKPOINT = ',I3,//,10X,
     &        'in subroutine <',A,'>',//,10X,A,//,2(1X,79('B'),/))
99002 FORMAT (2(/,20(' MPI')),/,5X,'MPI process number',I3,'  out of ',
     &        I3,' processes',//,2(20(' MPI'),/))
      END
C*==breakpoint_5.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE BREAKPOINT_5(CMNTMTT,CMNTIST,DROT_QLM)
C   ********************************************************************
C   *                                                                  *
C   *  rotate the charge moments and potential functions B and V       *
C   *  to an equivalent frame                                          *
C   *  specified by the orientation  TETNEW, PHINEW                    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ENERGY,ONLY:EFERMI,EMIN
      USE MOD_ANGMOM,ONLY:L_LM,M_LM
      USE MOD_SYMMETRY,ONLY:IQREPQ
      USE MOD_FILES,ONLY:IFILBREAK,FOUND_SECTION
      USE MOD_RMESH,ONLY:NRMAX,JRCRI
      USE MOD_SITES,ONLY:QMPHI,QMTET,IQBOT,IQTOP,NLMQMAD,NLQMAD,IQAT
      USE MOD_TYPES,ONLY:NLMFPMAX,NTMAX,BNST,VNST,BT,VT,IMT,ITBOT,ITTOP,
     &    TXT_T,LTXT_T
      USE MOD_CALCMODE,ONLY:IREL,BREAKPOINT,SOLVER_FP
      USE MOD_CONSTANTS,ONLY:SQRT_4PI
      IMPLICIT NONE
C*--BREAKPOINT_5140
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='BREAKPOINT_5')
C
C Dummy arguments
C
      REAL*8 CMNTIST(NLMFPMAX,NTMAX),CMNTMTT(NLMFPMAX,NTMAX),
     &       DROT_QLM(NLMQMAD,NLMQMAD,IQBOT:IQTOP)
C
C Local variables
C
      REAL*8 BLMTMP(:),BNST_NEW(:,:,:),CMNTIST_NEW(:,:),CMNTMTT_NEW(:,:)
     &       ,DROT_G_TO_LP(:,:),DROT_QLM_NEW(:,:,:),D_LMPLM,D_WRK(:,:),
     &       PHI_NEW,SUMA,SUMB,TET_NEW,VLMTMP(:),VNST_NEW(:,:,:)
      CHARACTER*40 FMT_QLM_BREAK
      INTEGER IA_ERR,IM,IQ,IT,LM,LMP,M,N
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DROT_G_TO_LP
      ALLOCATABLE CMNTIST_NEW,CMNTMTT_NEW,DROT_QLM_NEW
      ALLOCATABLE D_WRK,BNST_NEW,VNST_NEW
      ALLOCATABLE VLMTMP,BLMTMP
C
      FMT_QLM_BREAK = '(A,I7,2I3,2X,2I3,2X,A,2F20.14)'
C
C-----------------------------------------------------------------------
      IF ( BREAKPOINT.NE.5 )
     &     CALL STOP_MESSAGE(ROUTINE,'BREAKPOINT <> 5')
      IF ( IREL.NE.3 ) CALL STOP_MESSAGE(ROUTINE,'IREL <> 3')
C-----------------------------------------------------------------------
C
      WRITE (6,99002)
C
C-----------------------------------------------------------------------
C
      ALLOCATE (VNST_NEW(NRMAX,NLMFPMAX,NTMAX))
      ALLOCATE (BNST_NEW(NRMAX,NLMFPMAX,NTMAX))
      ALLOCATE (CMNTIST_NEW(NLMFPMAX,NTMAX))
      ALLOCATE (CMNTMTT_NEW(NLMFPMAX,NTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate VNS')
C
C=======================================================================
C              rotation matrices for REAL spherical harmonics
C=======================================================================
C
C------------------------------------------------ transformation g -> l'
C
      CALL INPUT_FIND_SECTION('CONTROL',0)
C
      TET_NEW = 0D0
      PHI_NEW = 0D0
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_SET_REAL('TETNEW',TET_NEW,9999D0,0)
         CALL SECTION_SET_REAL('PHINEW',PHI_NEW,9999D0,0)
      END IF
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      TET_NEW = 0D0
      PHI_NEW = -90D0
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      WRITE (6,99001) ROUTINE,TET_NEW,PHI_NEW
C
C-----------------------------------------------------------------------
C
      ALLOCATE (DROT_G_TO_LP(NLMQMAD,NLMQMAD))
      ALLOCATE (VLMTMP(NRMAX),BLMTMP(NRMAX))
C
      CALL ROTMAT_RYLM(NLQMAD,NLMQMAD,PHI_NEW,TET_NEW,0.0D0,
     &                 DROT_G_TO_LP)
C
C
      ALLOCATE (DROT_QLM_NEW(NLMQMAD,NLMQMAD,IQBOT:IQTOP))
      ALLOCATE (D_WRK(NLMQMAD,NLMQMAD))
C
      M = NLMQMAD
C
      DO IQ = IQBOT,IQTOP
C
C--------------------------- rotate Q-moments from local to global frame
C                 use inverse of DROT_QLM  with  1/DROT_QLM = DROT_QLM^T
C------------------------------------------------- transformation l -> g
C
         D_WRK(1:M,1:M) = TRANSPOSE(DROT_QLM(1:M,1:M,IQ))
C
C------------------------------------------------ transformation l -> l'
C
         DROT_QLM_NEW(1:M,1:M,IQ) = MATMUL(D_WRK(1:M,1:M),DROT_G_TO_LP(1
     &                              :M,1:M))
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         DROT_QLM_NEW(1:M,1:M,IQ) = DROT_G_TO_LP(1:M,1:M)
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END DO
C
C=======================================================================
C
C            charge-moments  CMNTMT  w.r.t. muffin tin sphere
C                            CMNTIS  w.r.t. interstitial regime
C  ---------------------------------------------------------------------
C
C  ---------------------------------------------------------------------
C                  charge moments in ORIGINAL local frame
C  ---------------------------------------------------------------------
C
      WRITE (6,99004) 'ORIGINAL'
C
      DO IT = ITBOT,ITTOP
C
         IQ = IQAT(1,IT)
C
         WRITE (6,99003) 'atom type  IT =',IT,TXT_T(IT)(1:LTXT_T(IT)),
     &                   '      '
C
         DO LM = 1,NLMQMAD
            IF ( ABS(CMNTMTT(LM,IT))+ABS(CMNTIST(LM,IT)).GT.1D-6 )
     &           WRITE (6,FMT=FMT_QLM_BREAK) '      ',LM,L_LM(LM),
     &                  M_LM(LM),NINT(QMTET(IQ)),NINT(QMPHI(IQ)),
     &                  SOLVER_FP,SQRT_4PI*CMNTMTT(LM,IT),
     &                  SQRT_4PI*CMNTIST(LM,IT)
         END DO
      END DO
C
C  ---------------------------------------------------------------------
C                  charge moments in ROTATED local frame
C  ---------------------------------------------------------------------
C
      WRITE (6,99004) 'ROTATED'
C
      DO IT = ITBOT,ITTOP
         WRITE (6,99003) 'atom type  IT =',IT,TXT_T(IT)(1:LTXT_T(IT)),
     &                   '      '
C
         IQ = IQAT(1,IT)
         IF ( IQREPQ(IQ).NE.IQ )
     &         CALL STOP_MESSAGE(ROUTINE,'IQREPQ(IQ) <> IQ')
C
C------------------------------------------------ transformation l -> l'
C
         DO LM = 1,NLMQMAD
            SUMA = 0D0
            SUMB = 0D0
            DO LMP = 1,NLMQMAD
               D_LMPLM = DROT_QLM_NEW(LMP,LM,IQ)
               SUMA = SUMA + CMNTMTT(LMP,IT)*D_LMPLM
               SUMB = SUMB + CMNTIST(LMP,IT)*D_LMPLM
            END DO
            CMNTMTT_NEW(LM,IT) = SUMA
            CMNTIST_NEW(LM,IT) = SUMB
C
            IF ( ABS(CMNTMTT_NEW(LM,IT))+ABS(CMNTIST_NEW(LM,IT))
     &           .GT.1D-6 ) WRITE (6,FMT=FMT_QLM_BREAK) 'BREAK ',LM,
     &                             L_LM(LM),M_LM(LM),NINT(QMTET(IQ)),
     &                             NINT(QMPHI(IQ)),SOLVER_FP,
     &                             SQRT_4PI*CMNTMTT_NEW(LM,IT),
     &                             SQRT_4PI*CMNTIST_NEW(LM,IT)
C
         END DO
C
      END DO
C
C  ---------------------------------------------------------------------
C                      reread potential functions
C  ---------------------------------------------------------------------
C
      REWIND (IFILBREAK)
      READ (IFILBREAK) EMIN,EFERMI,VT,BT,VNST,BNST
      CLOSE (IFILBREAK)
C
C  ---------------------------------------------------------------------
C             potential functions in ROTATED local frame
C  ---------------------------------------------------------------------
C
      DO IT = ITBOT,ITTOP
C
         IQ = IQAT(1,IT)
C
         IF ( IQREPQ(IQ).NE.IQ )
     &         CALL STOP_MESSAGE(ROUTINE,'IQREPQ(IQ) <> IQ')
C
         IM = IMT(IT)
         N = JRCRI(IM)
C
C------------------------------------------------ transformation l -> l'
C
         DO LM = 1,NLMQMAD
            VLMTMP(1:N) = 0.0D0
            BLMTMP(1:N) = 0.0D0
            DO LMP = 1,NLMQMAD
               D_LMPLM = DROT_QLM_NEW(LMP,LM,IQ)
               VLMTMP(1:N) = VLMTMP(1:N) + VNST(1:N,LMP,IT)*D_LMPLM
               BLMTMP(1:N) = BLMTMP(1:N) + BNST(1:N,LMP,IT)*D_LMPLM
            END DO
            VNST_NEW(1:N,LM,IT) = VLMTMP(1:N)
            BNST_NEW(1:N,LM,IT) = BLMTMP(1:N)
C
         END DO
C
      END DO
C
C-----------------------------------------------------------------------
C                   write  ROTATED  potential functions
C-----------------------------------------------------------------------
C
      OPEN (UNIT=IFILBREAK,FILE='BREAK_V_rotated.dat',
     &      FORM='UNFORMATTED')
      WRITE (IFILBREAK) EMIN,EFERMI,VT,BT,VNST_NEW,BNST_NEW
C
C-----------------------------------------------------------------------
C
      CALL STOP_BREAKPOINT(ROUTINE)
C
C-----------------------------------------------------------------------
C
99001 FORMAT (//,10X,A,//,10X,'setting up rotated potential set for',//,
     &        10X,'TETNEW = ',F10.5,5X,'PHINEW = ',F10.5)
99002 FORMAT (///,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*      *****   *****   *****    ***    *    *    ******      *'
     &  ,/,10X,
     &  '*      *    *  *    *  *       *   *   *   *     *           *'
     &  ,/,10X,
     &  '*      *    *  *    *  *      *     *  *  *      *           *'
     &  ,/,10X,
     &  '*      *****   *****   ****   *******  * *       *****       *'
     &  ,/,10X,
     &  '*      *    *  *  *    *      *     *  ** *           *      *'
     &  ,/,10X,
     &  '*      *    *  *   *   *      *     *  *   *     *    *      *'
     &  ,/,10X,
     &  '*      *****   *    *  *****  *     *  *    *     ****       *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99003 FORMAT (/,5X,'charge moments for ',A,I3,3X,A,//,A,21X,'LM  L  M',
     &        10X,'muffin-tin',12X,' interstitial')
99004 FORMAT (//,5X,'charge moments w.r.t. ',A,' local frame ')
      END
C*==breakpoint_5_6_ssite.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE BREAKPOINT_5_6_SSITE(ERYD,TSST,MSST,SSST,MEZZ,MEZJ)
C   ********************************************************************
C   *                                                                  *
C   *  rotate the charge moments and potential functions B and V       *
C   *  to an equivalent frame                                          *
C   *  specified by the orientation  TETNEW, PHINEW                    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,NMEMAX,NK,WKM2
      USE MOD_FILES,ONLY:FOUND_SECTION
      USE MOD_RMESH,ONLY:FULLPOT
      USE MOD_TYPES,ONLY:NTMAX,ITBOT,ITTOP,TXT_T,LTXT_T
      USE MOD_CALCMODE,ONLY:IREL,BREAKPOINT
      IMPLICIT NONE
C*--BREAKPOINT_5_6_SSITE410
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='BREAKPOINT_5_6_SSITE')
      REAL*8 TOL
      PARAMETER (TOL=1D-8)
C
C Dummy arguments
C
      COMPLEX*16 ERYD
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),SSST(NKMMAX,NKMMAX,NTMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 DROT_L_TO_LP(:,:)
      INTEGER IA_ERR,IT,M,N
      LOGICAL INITIALZE
      REAL*8 PHI_NEW,TET_NEW
      SAVE DROT_L_TO_LP
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALZE/.TRUE./
C
      ALLOCATABLE DROT_L_TO_LP
C
C-----------------------------------------------------------------------
      IF ( BREAKPOINT.NE.5 .AND. BREAKPOINT.NE.6 )
     &     CALL STOP_MESSAGE(ROUTINE,'BREAKPOINT <> 5 AND <> 6')
      IF ( .NOT.FULLPOT .AND. IREL.NE.3 )
     &      CALL STOP_MESSAGE(ROUTINE,'.NOT. FULLPOT  and  IREL <> 3')
C-----------------------------------------------------------------------
C
      WRITE (6,99004)
C
CIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
      IF ( INITIALZE .AND. BREAKPOINT.EQ.5 ) THEN
C
C------------------------------------------------ transformation l -> l'
C
         CALL INPUT_FIND_SECTION('CONTROL',0)
C
         TET_NEW = 0D0
         PHI_NEW = 0D0
         IF ( FOUND_SECTION ) THEN
            CALL SECTION_SET_REAL('TETNEW',TET_NEW,9999D0,0)
            CALL SECTION_SET_REAL('PHINEW',PHI_NEW,9999D0,0)
         END IF
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         TET_NEW = 0D0
         PHI_NEW = 90D0
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         WRITE (6,99003) ROUTINE,TET_NEW,PHI_NEW
C
C-----------------------------------------------------------------------
C
         ALLOCATE (DROT_L_TO_LP(NKMMAX,NKMMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate DROT')
C
         CALL ROTMAT(NK,3,PHI_NEW,TET_NEW,0.0D0,DROT_L_TO_LP,NKMMAX)
C
         INITIALZE = .FALSE.
C
      END IF
CIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C
      WRITE (6,99001) ROUTINE,BREAKPOINT,ERYD
C
      N = NKM
      M = NKMMAX
C
      DO IT = ITBOT,ITTOP
C
         WRITE (6,99002) IT,TXT_T(IT)(1:LTXT_T(IT))
C
C--------------------------------------------------- rotate from L to LP
         IF ( BREAKPOINT.EQ.5 ) THEN
C
            CALL CMATSTRUCT('ROT MATRIX              ',DROT_L_TO_LP,N,M,
     &                      3,3,1,TOL,6)
C
            CALL ROTATE(TSST(1,1,IT),'L->G',WKM2,N,DROT_L_TO_LP,M)
C
            CALL CMATSTRUCT('T-MATRIX      (kappa,mu)',WKM2,N,M,3,3,1,
     &                      TOL,6)
C
            CALL ROTATE(MSST(1,1,IT),'L->G',WKM2,N,DROT_L_TO_LP,M)
C
            CALL CMATSTRUCT('M-MATRIX      (kappa,mu)',WKM2,N,M,3,3,1,
     &                      TOL,6)
C
            CALL ROTATE(SSST(1,1,IT),'L->G',WKM2,N,DROT_L_TO_LP,M)
C
            CALL CMATSTRUCT('S-MATRIX      (kappa,mu)',WKM2,N,M,3,3,1,
     &                      TOL,6)
C
            CALL ROTATE(MEZZ(1,1,IT,1),'L->G',WKM2,N,DROT_L_TO_LP,M)
C
            CALL CMATSTRUCT('MEZZ-MATRIX',WKM2,N,M,3,3,1,1.0D-9,6)
C
            CALL ROTATE(MEZJ(1,1,IT,1),'L->G',WKM2,N,DROT_L_TO_LP,M)
C
            CALL CMATSTRUCT('MEZJ-MATRIX',WKM2,N,M,3,3,1,1.0D-9,6)
C
C---------------------------------------------------- no rotation needed
         ELSE
C
            CALL CMATSTRUCT('T-MATRIX      (kappa,mu)',TSST(1,1,IT),N,M,
     &                      3,3,1,TOL,6)
C
            CALL CMATSTRUCT('M-MATRIX      (kappa,mu)',MSST(1,1,IT),N,M,
     &                      3,3,1,TOL,6)
C
            CALL CMATSTRUCT('S-MATRIX      (kappa,mu)',SSST(1,1,IT),N,M,
     &                      3,3,1,TOL,6)
C
            CALL CMATSTRUCT('MEZZ-MATRIX',MEZZ(1,1,IT,1),N,M,3,3,1,
     &                      1.0D-9,6)
C
            CALL CMATSTRUCT('MEZJ-MATRIX',MEZJ(1,1,IT,1),N,M,3,3,1,
     &                      1.0D-9,6)
C
         END IF
C
      END DO
C-----------------------------------------------------------------------
C
99001 FORMAT (//,10X,'running ',A,/,10X,'for BREAKPOINT =',I3,/,10X,
     &        'energy    ERYD =',2F15.10,/)
99002 FORMAT (//,1X,79('T'),//,10X,'single site matrices for type ',i3,
     &        3X,A,/,1X,79('T'),/)
99003 FORMAT (//,10X,A,//,10X,'setting up rotated potential set for',//,
     &        10X,'TETNEW = ',F10.5,5X,'PHINEW = ',F10.5)
99004 FORMAT (///,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*           *****   *****   *****    ***    *    *           *'
     &  ,/,10X,
     &  '*           *    *  *    *  *       *   *   *   *            *'
     &  ,/,10X,
     &  '*           *    *  *    *  *      *     *  *  *             *'
     &  ,/,10X,
     &  '*           *****   *****   ****   *******  * *              *'
     &  ,/,10X,
     &  '*           *    *  *  *    *      *     *  ** *             *'
     &  ,/,10X,
     &  '*           *    *  *   *   *      *     *  *   *            *'
     &  ,/,10X,
     &  '*           *****   *    *  *****  *     *  *    *           *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
      END
