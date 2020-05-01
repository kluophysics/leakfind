C*==tau_std_ktetra.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE TAU_STD_KTETRA(ICPAFLAG,CPACHNG,ERYD,P,IPRINT,ITCPA,
     &                          ICPACONV,PHASK,IE,TSST,MSST,TSSQ,MSSQ,
     &                          TAUQ)
C   ********************************************************************
C   *                                                                  *
C   *    PERFORM THE K-SPACE INTEGRAL USING THE THEDRAHEDRON METHOD    *
C   *                                                                  *
C   *   TK      NELMT non-zero elements of TAU(K)                      *
C   *           tabulated for NKTET K-points                           *
C   *   TW      NELMT non-zero elements of TAU                         *
C   *           determined by tetrahedron integration                  *
C   *           over 1 wedge of BZ                                     *
C   *   KTET    table of tetrahedron k-points in (2*PI/a)              *
C   *   IKCTET  indexlist for k-points belonging to NTETS tetrahedrons *
C   *                                                                  *
C   *   NKTET   number of k-points                                     *
C   *   NTETS   number of tets                                         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:CONC,NTMAX
      USE MOD_LATTICE,ONLY:ALAT
      USE MOD_CPA,ONLY:CPATOL,ITCPAMAX,ICPAALG,NCPA
      USE MOD_SITES,ONLY:NQ,NQMAX,ITOQ,NOQ,ICPA,DROTQ
      USE MOD_ANGMOM,ONLY:NKMMAX,IND0Q,NKMQ,NKKR
      USE MOD_CALCMODE,ONLY:ITEST,LLOYD,KMROT
      USE MOD_SYMMETRY,ONLY:NWEDGE,NELMTMAX,MROTK,IWEDGEROT
      USE MOD_KSPACE,ONLY:NKTET,NGFEP,GFEP,WTAUUV,VTAUUV,UTAUUV,NTAUUV,
     &    QTBZ,JTBZ,ITBZ,NELMT,GBAD,NGBAD,KTET,IKCTET,NTETS
      USE MOD_CONSTANTS,ONLY:PI,C0,C1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='TAU_STD_KTETRA')
C
C Dummy arguments
C
      REAL*8 CPACHNG
      COMPLEX*16 ERYD,P
      INTEGER ICPACONV,ICPAFLAG,IE,IPRINT,ITCPA
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           PHASK(IE),TAUQ(NKMMAX,NKMMAX,NQMAX),
     &           TSSQ(NKMMAX,NKMMAX,NQMAX),TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 ADD,CSUM,DET(:),DETFE(:),DETFEW,DETGT(:),DETGTW,DETK,
     &           DETKFE,DETKGT,EDU,LDET(:),MAUX(:,:),MSSQLU(:,:,:),
     &           PHASE,TAUK(:,:),TBZ(:),TK(:,:),TW(:),W1(:,:)
      REAL*8 BRA(:),CPACORR,CPAERR,CPAERRL,KN2,KVEC(3),KVECL(3),SA,SCAL,
     &       SD,WTET(:)
      INTEGER I,I1,I2,I3,IA_ERR,IBRA(:),ICC,IG,IGFEP,IK,IKM,IQ,IROT,
     &        IROTLAST,IWEDGE,J,N
      LOGICAL RVEC_SAME
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE IBRA,LDET,MAUX,TAUK,TK,TW,BRA,DET,TBZ,W1,MSSQLU
      ALLOCATABLE DETGT,DETFE,WTET
C
      ALLOCATE (DETGT(NKTET),DETFE(NKTET),DET(NKTET))
      ALLOCATE (IBRA(NKTET),LDET(NKTET),BRA(NKTET),WTET(NTETS))
      WTET(1:NTETS) = 1D0
C
      IF ( LLOYD ) ALLOCATE (MSSQLU(NKMMAX,NKMMAX,NQMAX))
      ALLOCATE (TAUK(NKKR,NKKR),TK(NELMTMAX,NKTET))
      ALLOCATE (TW(NELMTMAX),TBZ(NELMTMAX))
      ALLOCATE (MAUX(NKKR,NKKR),W1(NKMMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MAUX')
C
      EDU = ERYD/(2*PI/ALAT)**2
C
C ------------------ calculate energy - dependent terms of str.constants
C
      CALL STRCC(ERYD,.FALSE.)
C
      IROTLAST = 0
C
      CALL CINIT(NKMMAX*NKMMAX*NQ,TAUQ)
C
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      ICPACONV = 0
      CPAERRL = 1.0D+6
      ITCPA = 0
 100  CONTINUE
      ITCPA = ITCPA + 1
C
      CALL CINIT(NELMT,TBZ)
C
      IF ( LLOYD ) THEN
C
         DO IQ = 1,NQ
            N = NKMQ(IQ)
            DO J = 1,N
               CALL ZCOPY(N,MSSQ(1,J,IQ),1,MSSQLU(1,J,IQ),1)
            END DO
            CALL CINVLU(MSSQLU(1,1,IQ),W1,N,NKMMAX)
         END DO
C
         STOP 'LLOYD not yet available in TAU_STD_KTETRA'
      END IF
C
      DETKGT = C0
      DETKFE = C0
      PHASK(IE) = C0
C
C WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
C
      DO IWEDGE = 1,NWEDGE
         IROT = IWEDGEROT(IWEDGE)
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
         DO IK = 1,NKTET
C
            CALL DGEMV('N',3,3,1D0,MROTK(1,1,IROT),3,KTET(1,IK),1,0D0,
     &                 KVEC,1)
C
            IF ( ABS(KVEC(1))+ABS(KVEC(2))+ABS(KVEC(3)).LT.1D-5 )
     &           KVEC(3) = 1D-5
C
C-- if last and current k-points are identical: use previous results ---
            IF ( IWEDGE.GT.1 ) THEN
               CALL DGEMV('N',3,3,1D0,MROTK(1,1,IROTLAST),3,KTET(1,IK),
     &                    1,0D0,KVECL,1)
C
               IF ( RVEC_SAME(3,KVEC,KVECL,1D-5) ) CYCLE
            END IF
C
            CALL STRSET(IK,KVEC,TAUK,MAUX,P)
C
            CALL SETKKR(NQ,NKMQ,IND0Q,MAUX,MSSQ,NQMAX,NKKR,NKMMAX)
C
            CALL CINVLU(MAUX,TAUK,NKKR,NKKR)
C
C--------------------------------------------------------------- pack TK
            I2 = 0
            DO I1 = 1,NELMT
               CSUM = C0
               DO I3 = 1,NTAUUV(I1)
                  I2 = I2 + 1
                  CSUM = CSUM + WTAUUV(I2)*TAUK(UTAUUV(I2),VTAUUV(I2))
               END DO
               TK(I1,IK) = CSUM
            END DO
C-----------------------------------------------------------------------
C
            IF ( IK.EQ.1 ) THEN
               IF ( IWEDGE.EQ.1 ) THEN
                  SCAL = 1.0D0
                  I = 0
                  DO IQ = 1,NQ
                     DO IKM = 1,NKMQ(IQ)
                        I = I + 1
                        SCAL = MAX(SCAL,
     &                         CDABS(MAUX(I,I)/MSSQ(IKM,IKM,IQ)))
                     END DO
                  END DO
                  SCAL = 1.0D0/SCAL
               END IF
            END IF
C
            BRA(IK) = 0D0
            DET(IK) = C1
            I = 0
            DO IQ = 1,NQ
               DO IKM = 1,NKMQ(IQ)
                  I = I + 1
                  ADD = MAUX(I,I)/MSSQ(IKM,IKM,IQ)
                  IF ( ABS(DIMAG(EDU)).LT.1D-8 )
     &                 ADD = DCMPLX(DREAL(ADD),0.0D0)
C
                  SD = SIGN(1D0,DIMAG(DET(IK)))
                  SA = SIGN(1D0,DIMAG(ADD))
                  DET(IK) = DET(IK)*ADD*SCAL
                  BRA(IK) = BRA(IK) + SA*2.5D-1*(SD+SA)
     &                      *(SD-SIGN(1D0,DIMAG(DET(IK))))
               END DO
            END DO
C
            DO IG = 1,NGBAD
               KN2 = (KVEC(1)+GBAD(1,IG))**2 + (KVEC(2)+GBAD(2,IG))
     &               **2 + (KVEC(3)+GBAD(3,IG))**2
C
               ADD = KN2 - EDU
               SD = SIGN(1D0,DIMAG(DET(IK)))
               SA = SIGN(1D0,DIMAG(ADD))
               DET(IK) = DET(IK)*ADD
               BRA(IK) = BRA(IK) + SA*2.5D-1*(SD+SA)
     &                   *(SD-SIGN(1D0,DIMAG(DET(IK))))
            END DO
C
C--------------------------------------------------------- LLOYD FORMULA
            IF ( LLOYD ) THEN
               DETK = C0
               DO IQ = 1,NQ
                  I1 = IND0Q(IQ)
                  N = NKMQ(IQ)
                  DO J = 1,N
                     I = I1 + J
                     DETK = DETK + LOG(MAUX(I,I)/MSSQLU(J,J,IQ))
                  END DO
               END DO
               DETGT(IK) = DETK
C
               DETK = C0
               DO IGFEP = 1,NGFEP
                  KN2 = 0D0
                  DO ICC = 1,3
                     KN2 = KN2 + (KVEC(ICC)+GFEP(ICC,IGFEP))**2
                  END DO
                  DETK = DETK + LOG(KN2-EDU)
               END DO
               DETFE(IK) = DETK
C
            END IF
C-----------------------------------------------------------------------
C
         END DO
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
         CALL TETSUM(WTET,DET,LDET,BRA,IBRA,TW,PHASE,IKCTET,NKTET,NTETS,
     &               TK,DETGT,DETGTW,DETFE,DETFEW,NELMTMAX,NELMT)
C
         DETKGT = DETKGT - DETGTW/DBLE(NWEDGE)
         DETKFE = DETKFE - DETFEW/DBLE(NWEDGE)
         DO I = 1,NELMT
            TBZ(I) = TBZ(I) + TW(I)/DBLE(NWEDGE)
         END DO
C
         IROTLAST = IROT
      END DO
C
C WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
C
C
C------------------------------------------------------------ unpack TBZ
      DO I = 1,NELMT
         IQ = QTBZ(I)
         TAUQ(ITBZ(I),JTBZ(I),IQ) = TBZ(I)
      END DO
C-----------------------------------------------------------------------
C
      IF ( NCPA.GT.0 ) THEN
C
         IF ( ICPAALG.EQ.1 ) THEN
            CALL CPAMILLS(CPAERR,CPACORR,CPACHNG,IPRINT,ICPA,NQ,NKMQ,
     &                    NOQ,ITOQ,CONC,TSST,MSST,TSSQ,MSSQ,TAUQ,NTMAX,
     &                    NQMAX,NKMMAX)
C
         ELSE
            CALL CPANESBET(CPAERR,CPACORR,CPACHNG,IPRINT,ICPA,NQ,NKMQ,
     &                     NOQ,ITOQ,CONC,TSST,MSST,TSSQ,MSSQ,TAUQ,KMROT,
     &                     DROTQ,NTMAX,NQMAX,NKMMAX)
C
         END IF
C
         IF ( IPRINT.EQ.1 ) WRITE (6,99004) CPAERR,CPACORR,CPACHNG
C
         IF ( CPAERR.LE.CPATOL ) THEN
            ICPACONV = 1
            IF ( IPRINT.GE.0 ) WRITE (6,99001) ITCPA,CPAERR,CPACORR,
     &                                CPACHNG
         ELSE IF ( ITCPA.GT.ITCPAMAX ) THEN
            WRITE (6,99002) ITCPA,CPAERR,CPACORR,CPACHNG
            ICPAFLAG = 1
         ELSE IF ( CPAERR.GT.20*CPAERRL ) THEN
            WRITE (6,99003) ITCPA
            WRITE (6,99004) CPAERR,CPACORR,CPACHNG
            ICPAFLAG = 2
         ELSE
            CPAERRL = CPAERR
            GOTO 100
         END IF
C
      END IF
C
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IF ( ITEST.EQ.5 ) THEN
         WRITE (77,'(2f14.5,3x,2f14.5)') DREAL(ERYD),DIMAG(PHASK(IE))
         WRITE (78,'(2f14.5,3x,2f14.5)') DREAL(ERYD),DIMAG(DETKFE)
         WRITE (79,'(2f14.5,3x,2f14.5)') DREAL(ERYD),DIMAG(DETKGT)
      END IF
C         WRITE (100,'(20e14.5)') DREAL(ERYD),phase,
C     &    DIMAG(DETKFE),DIMAG(DETKGT),DIMAG(DETKFE+DETKGT)
      WRITE (100,'(20e14.5)') DREAL(ERYD),DIMAG(2*DETKFE+DETKGT),
     &                        DIMAG(PHASE)
C
      IF ( IPRINT.GE.2 .AND. NCPA.GT.0 )
     &      WRITE (12,'(''E '',2F10.5,'' CPA '',I5,3E12.5)') ERYD,
     &     ICPAFLAG,CPAERR,CPACORR,CPACHNG
C
      DEALLOCATE (IBRA,LDET,MAUX,TAUK,TK,TW,BRA,DET,TBZ,W1)
      DEALLOCATE (DETGT,DETFE)
      IF ( LLOYD ) DEALLOCATE (MSSQLU)
C
C ----------------------------------------------------------------------
C
99001 FORMAT (' CPA converged after',I3,' iterations   ','ERR:',F9.6,
     &        ' CORR:',F9.6,' CHNG:',F9.6)
99002 FORMAT (' CPA-cycle  NOT  converged after ',I3,
     &        ' iterations: ERROR ',F12.8,' CORRECTION ',F15.8,
     &        ' CHANGE ',F15.8,10('!'))
99003 FORMAT (' CPA: ERROR increased by more than ','20*TOL for ITCPA=',
     &        i4,' >>>> iteration stopped ',5X,10('!'))
99004 FORMAT (' CPA: ERROR ',F12.8,'    CORRECTION ',F15.8,
     &        '    CHANGE ',F15.8)
      END
