C*==tetarrsiz.f    processed by SPAG 6.70Rc at 21:40 on 19 Dec 2016
      SUBROUTINE TETARRSIZ(BRAVAIS,NKTET0,NFTET,NKTET,NTETS,BBAS)
C   ********************************************************************
C   *                                                                  *
C   *  specify tetrahedron grid by  NFTET  or  NKTET0                  *
C   *                                                                  *
C   *  NFTET = 0  set NFTET such that NKTET >= NKTET0                  *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *         1 triclinic   primitive      -1     C_i                  *
C   *         2 monoclinic  primitive      2/m    C_2h                 *
C   *         3 monoclinic  base centered  2/m    C_2h                 *
C   *         4 orthorombic primitive      mmm    D_2h                 *
C   *         5 orthorombic base-centered  mmm    D_2h                 *
C   *         6 orthorombic body-centered  mmm    D_2h                 *
C   *         7 orthorombic face-centered  mmm    D_2h                 *
C   *         8 tetragonal  primitive      4/mmm  D_4h                 *
C   *         9 tetragonal  body-centered  4/mmm  D_4h                 *
C   *        10 trigonal    primitive      -3m    D_3d                 *
C   *        11 hexagonal   primitive      6/mmm  D_6h                 *
C   *        12 cubic       primitive      m3m    O_h                  *
C   *        13 cubic       face-centered  m3m    O_h                  *
C   *        14 cubic       body-centered  m3m    O_h                  *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NFCC,NBCC
      PARAMETER (NFCC=8,NBCC=12)
C
C Dummy arguments
C
      INTEGER BRAVAIS,NFTET,NKTET,NKTET0,NTETS
      REAL*8 BBAS(3,3)
C
C Local variables
C
      INTEGER I
      REAL*8 KNBCC(3,NBCC),KNFCC(3,NFCC)
C
C*** End of declarations rewritten by SPAG
C
      DATA KNFCC/1D0,1D0,1D0,1D0, - 1D0,1D0, - 1D0, - 1D0,1D0, - 1D0,
     &     1D0,1D0,1D0,1D0, - 1D0,1D0, - 1D0, - 1D0, - 1D0, - 1D0,
     &     - 1D0, - 1D0,1D0, - 1D0/
C
      DATA KNBCC/1D0,1D0,0D0,0D0, - 1D0,1D0,1D0,0D0,1D0,1D0,0D0, - 1D0,
     &     - 1D0,1D0,0D0, - 0D0, - 1D0,1D0, - 1D0,0D0,1D0, - 1D0,0D0,
     &     - 1D0,0D0,1D0,1D0,0D0,1D0, - 1D0,0D0, - 1D0,1D0,0D0, - 1D0,
     &     - 1D0/
C
C*** End of declarations rewritten by SPAG
C
      IF ( BRAVAIS.EQ.4 ) THEN
C
         STOP
C
C------------------------------------------------------------------- FCC
      ELSE IF ( BRAVAIS.EQ.13 ) THEN
C
         IF ( NFTET.EQ.0 ) THEN
            DO I = 1,100
               NKTET = NINT((14D0*I+27D0*I**2+16D0*I**3)/3D0) + 1
               IF ( NKTET.GT.NKTET0 ) THEN
                  NFTET = I
                  EXIT
               END IF
            END DO
         END IF
C
         NTETS = 32*NFTET**3
         NKTET = NINT((14D0*NFTET+27D0*NFTET**2+16D0*NFTET**3)/3D0) + 1
C
         CALL TETCHECK(BBAS,KNFCC,NFCC)
C
C------------------------------------------------------------------- BCC
      ELSE IF ( BRAVAIS.EQ.14 ) THEN
C
         IF ( NFTET.EQ.0 ) THEN
            DO I = 1,100
               NKTET = NINT((13D0*I+18D0*I**2+8D0*I**3)/3D0) + 1
               IF ( NKTET.GT.NKTET0 ) THEN
                  NFTET = I
                  EXIT
               END IF
            END DO
         END IF
C
         NTETS = 16*NFTET**3
         NKTET = NINT((13D0*NFTET+18D0*NFTET**2+8D0*NFTET**3)/3D0) + 1
C
         CALL TETCHECK(BBAS,KNBCC,NBCC)
C
      ELSE
         WRITE (6,99004) ' ARGUMENT         BRAVAIS   = ',BRAVAIS
         STOP 'in  <TETARRSIZ>   BRAVAIS lattice not implemented '
      END IF
C
C-----------------------------------------------------------------------
      IF ( NFTET.EQ.0 ) THEN
         WRITE (6,99005) BRAVAIS,NKTET0
         STOP 'in  <TETARRSIZ>'
      END IF
C-----------------------------------------------------------------------
      WRITE (6,99001)
      WRITE (6,99003) ' start                    NKTET0: ',NKTET0
      WRITE (6,99003) ' K-vectors                NKTET:  ',NKTET
      WRITE (6,99003) ' tetrahedra               NTETS:  ',NTETS
      WRITE (6,99003) ' TET-grid                 NFTET:  ',NFTET
      WRITE (6,99002)
C-----------------------------------------------------------------------
99001 FORMAT (/,10X,44('='),/,29X,'<TETARRSIZ>',/,10X,44('='),/,10X,
     &        ' set-up of tetrahedron - mesh ')
99002 FORMAT (10X,44('='),/)
99003 FORMAT (10X,A,I6,:,4X,'(',I6,')')
99004 FORMAT (A,I7)
99005 FORMAT (/,1X,79('#'),/,' NFTET could not be set for BRAVAIS =',I3,
     &        3X,'NKTET0  =',I3)
      END
C*==tetcheck.f    processed by SPAG 6.70Rc at 21:40 on 19 Dec 2016
      SUBROUTINE TETCHECK(BBAS,KN,N)
C   ********************************************************************
C   *                                                                  *
C   *   check whether the list of reciprocal mesh points  KN  that     *
C   *   is consistent with the reciprocal lattice assumed for the      *
C   *   TET... soubroutines can be represented as                      *
C   *                                                                  *
C   *             KN = n1*BG1 + n2*BG2 + n3*BG3                        *
C   *                                                                  *
C   *   with ni integers and  BGi the set of basis vectors for the     *
C   *   reciprocal lattice  that are derived from the real space       *
C   *   basis vectors read in or tabulated somewhere else              *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N
      REAL*8 BBAS(3,3),KN(3,N)
C
C Local variables
C
      REAL*8 BB(3,3),BBINV(3,3),BGMAT(3,3),BK(3),CF(3)
      REAL*8 DDOT
      INTEGER I,IERR,J
C
C*** End of declarations rewritten by SPAG
C
      IERR = 0
C
      BGMAT(1:3,1:3) = BBAS(1:3,1:3)
C
      DO J = 1,3
         DO I = 1,3
            BB(I,J) = DDOT(3,BGMAT(1,I),1,BGMAT(1,J),1)
         END DO
      END DO
C
      CALL RINVGJ(BBINV,BB,3,3)
C
      DO I = 1,N
         DO J = 1,3
            BK(J) = DDOT(3,BGMAT(1,J),1,KN(1,I),1)
         END DO
C
         CALL DGEMV('N',3,3,1D0,BBINV,3,BK,1,0D0,CF,1)
C
         DO J = 1,3
            IF ( (CF(J)-NINT(CF(J))).GT.1D-7 ) THEN
               WRITE (6,*) '############# non-integer found ',I,J,CF(J)
               IERR = 1
            END IF
         END DO
      END DO
C
      IF ( IERR.NE.0 ) STOP 'in  <TETCHECK> -- lattice inconsistent'
      END
C*==tetgen.f    processed by SPAG 6.70Rc at 21:40 on 19 Dec 2016
      SUBROUTINE TETGEN(BRAVAIS,BOA,COA,NFTET,KTET,NKTET,IKCTET,NTETS)
C   ********************************************************************
C   *                                                                  *
C   *    THE FOLLOWING ROUTINES GENERATE A SET OF TETRAHEDRA WHICH     *
C   *    FILL A IREDUCIBLE WEDGE OF THE B-Z.                           *
C   *    IKCTET(I,J) (I=1,2,3,4) GIVES THE ADDRESS OF FOUR CORNERS     *
C   *    OF THE J-TH TETRAHEDRON.                                      *
C   *    NTETS IS THE TOTAL NUMBER OF THE TETRAHEDRA.                  *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *         1 triclinic   primitive      -1     C_i                  *
C   *         2 monoclinic  primitive      2/m    C_2h                 *
C   *         3 monoclinic  base centered  2/m    C_2h                 *
C   *         4 orthorombic primitive      mmm    D_2h                 *
C   *         5 orthorombic base-centered  mmm    D_2h                 *
C   *         6 orthorombic body-centered  mmm    D_2h                 *
C   *         7 orthorombic face-centered  mmm    D_2h                 *
C   *         8 tetragonal  primitive      4/mmm  D_4h                 *
C   *         9 tetragonal  body-centered  4/mmm  D_4h                 *
C   *        10 trigonal    primitive      -3m    D_3d                 *
C   *        11 hexagonal   primitive      6/mmm  D_6h                 *
C   *        12 cubic       primitive      m3m    O_h                  *
C   *        13 cubic       face-centered  m3m    O_h                  *
C   *        14 cubic       body-centered  m3m    O_h                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LATTICE,ONLY:BBAS
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER N_EDGEMAX
      PARAMETER (N_EDGEMAX=20)
C
C Dummy arguments
C
      REAL*8 BOA,COA
      INTEGER BRAVAIS,NFTET,NKTET,NTETS
      INTEGER IKCTET(4,NTETS)
      REAL*8 KTET(3,NKTET)
C
C Local variables
C
      CHARACTER*20 FILNAM
      INTEGER I,IA_ERR,IERR,IK,IW,IWK(:,:,:),J,LFILNAM,NW1,NW2,NW3,
     &        N_EDGE
      REAL*8 KA_EDGE(:,:),KE_EDGE(:,:),S,SCLWGT
      CHARACTER*8 LBL_EDGE(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE IWK,LBL_EDGE,KA_EDGE,KE_EDGE
C
      NW1 = 4*NFTET + 1
      NW2 = 3*NFTET + 1
      NW3 = 2*NFTET + 1
      ALLOCATE (IWK(NW1,NW2,NW3),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'ALLOC: TETGEN -> IWK'
C
      ALLOCATE (KA_EDGE(3,N_EDGEMAX),KE_EDGE(3,N_EDGEMAX))
      ALLOCATE (LBL_EDGE(N_EDGEMAX))
      LBL_EDGE(1:N_EDGEMAX) = 'NO LABEL'
C
      IERR = 0
C
      IF ( BRAVAIS.EQ.4 ) THEN
         CALL TETSOR(NFTET,NKTET,KTET,IKCTET,NTETS,IERR,BOA,COA)
      ELSE IF ( BRAVAIS.EQ.8 ) THEN
         CALL TETHEX(NFTET,NKTET,KTET,IKCTET,NTETS,IERR,COA,BRAVAIS)
      ELSE IF ( BRAVAIS.EQ.11 ) THEN
         CALL TETHEX(NFTET,NKTET,KTET,IKCTET,NTETS,IERR,COA,BRAVAIS)
      ELSE IF ( BRAVAIS.EQ.12 ) THEN
         CALL TETSIC(NFTET,NKTET,KTET,IKCTET,NTETS,IERR)
      ELSE IF ( BRAVAIS.EQ.13 ) THEN
         CALL TETFCC(NFTET,NKTET,KTET,IKCTET,NTETS,IWK,NW1,NW2,NW3,IERR)
      ELSE IF ( BRAVAIS.EQ.14 ) THEN
         CALL TETBCC(NFTET,NKTET,KTET,IKCTET,NTETS,IWK,NW1,NW2,NW3,IERR)
      ELSE
         WRITE (6,99008) ' ARGUMENT         BRAVAIS   = ',BRAVAIS
         STOP 'in  <TETGEN>   BRAVAIS lattice not implemented '
      END IF
C
      CALL KDIRTAB(BRAVAIS,10,N_EDGE,LBL_EDGE,KA_EDGE,KE_EDGE,BBAS)
C
      KTET(1,1) = 1D-10
      KTET(2,1) = 1D-10
      KTET(3,1) = 1D-10
C
      IF ( IERR.NE.0 ) THEN
         WRITE (6,99001) ' STOP in  <TETGEN> ----  arrays too small '
         WRITE (6,99001) ' ARGUMENT         BRAVAIS   = ',BRAVAIS
         WRITE (6,99001) ' ARGUMENT         NFTET = ',NFTET
         WRITE (6,99001) ' REQUIRED         NKTET = ',NKTET
         WRITE (6,99001) ' REQUIRED         NTETS = ',NTETS
         STOP
      END IF
C
      IF ( NKTET.GE.100000 ) RETURN
C
C-----------------------------------------------------------------------
C                  create input for the  rasmol  program
C-----------------------------------------------------------------------
C
      FILNAM = 'rasmol_tet_mesh'
      LFILNAM = LEN_TRIM(FILNAM)
      IW = 80
      S = 15D0
      OPEN (IW,FILE=FILNAM(1:LFILNAM)//'.pdb')
C
      WRITE (IW,FMT=99002) 'for '
C
C------------------------------------------ specify corners of unit cell
C
      WRITE (IW,FMT=99005) 1,1,0D0,0D0,0D0,0D0
      DO I = 1,3
         WRITE (IW,FMT=99005) (I+1),(I+1),S*BBAS(1,I),S*BBAS(2,I),
     &                        S*BBAS(3,I),0D0
      END DO
C
      WRITE (IW,FMT=99005) 5,5,S*(BBAS(1,1)+BBAS(1,2)),
     &                     S*(BBAS(2,1)+BBAS(2,2)),
     &                     S*(BBAS(3,1)+BBAS(3,2)),0D0
      WRITE (IW,FMT=99005) 6,6,S*(BBAS(1,1)+BBAS(1,3)),
     &                     S*(BBAS(2,1)+BBAS(2,3)),
     &                     S*(BBAS(3,1)+BBAS(3,3)),0D0
      WRITE (IW,FMT=99005) 7,7,S*(BBAS(1,3)+BBAS(1,2)),
     &                     S*(BBAS(2,3)+BBAS(2,2)),
     &                     S*(BBAS(3,3)+BBAS(3,2)),0D0
      WRITE (IW,FMT=99005) 8,8,S*(BBAS(1,1)+BBAS(1,2)+BBAS(1,3)),
     &                     S*(BBAS(2,1)+BBAS(2,2)+BBAS(2,3)),
     &                     S*(BBAS(3,1)+BBAS(3,2)+BBAS(3,3)),0D0
C
C-------------------------------------------------- specify edges of IBZ
C
      I = 8
      DO J = 1,N_EDGE
         I = I + 1
         WRITE (IW,FMT=99005) I,I,S*KA_EDGE(1:3,J),3D0
         I = I + 1
         WRITE (IW,FMT=99005) I,I,S*KE_EDGE(1:3,J),3D0
      END DO
C
      SCLWGT = 2D0
C
      DO IK = 1,NKTET
         I = I + 1
         WRITE (IW,FMT=99004) I,I,S*KTET(1:3,IK),SCLWGT
      END DO
C
      WRITE (IW,FMT=99007)
      DO I = 9,9 + 2*N_EDGE,2
         WRITE (IW,FMT=99006) I,I + 1
      END DO
      WRITE (IW,'(''END'')')
      CLOSE (IW)
C
C--------------------------------------------------- write RASMOL script
C
      OPEN (IW,FILE=FILNAM(1:LFILNAM)//'.ras')
      WRITE (IW,*) 'load '''//FILNAM(1:LFILNAM)//'.pdb'' '
      WRITE (IW,*) 'background white'
      WRITE (IW,*) 'color temperature'
      WRITE (IW,99001) 'select ',1,'-',8 + 2*N_EDGE
      WRITE (IW,*) 'cpk 0'
      WRITE (IW,99001) 'select ',8 + 2*N_EDGE + 1,'-',
     &                 8 + 2*N_EDGE + NKTET
      WRITE (IW,*) 'cpk 50'
      WRITE (IW,*) 'select all '
      WRITE (IW,*) 'set axes on '
C
      CLOSE (IW)
C
      WRITE (6,99003) FILNAM(1:LFILNAM),FILNAM(1:LFILNAM)
C
      DEALLOCATE (IWK)
C
C-----------------------------------------------------------------------
C
99001 FORMAT (A,I6,A,I6)
99002 FORMAT ('HEADER    k-mesh                 ',A,/,
     &        'SOURCE    SPRKKR - program       ',/,
     &        'AUTHOR    H. Ebert               ',/,
     &        'REMARK    None                   ')
99003 FORMAT (/,10X,'k-mesh stored in rasmol data-file ',A,'.pdb',/,10X,
     &        'view via:   rasmol  -script ',A,'.ras',/)
99004 FORMAT ('ATOM  ',I5,'          ',I5,'    ',3F8.3,'  0.00',F8.3)
99005 FORMAT ('HETATM',I5,'          ',I5,'    ',3F8.3,'  0.00',F8.3)
99006 FORMAT ('CONECT',2I5)
99007 FORMAT ('CONECT    1    2',/,'CONECT    1    3',/,
     &        'CONECT    1    4',/,'CONECT    5    2',/,
     &        'CONECT    5    3',/,'CONECT    5    8',/,
     &        'CONECT    6    2',/,'CONECT    6    4',/,
     &        'CONECT    6    8',/,'CONECT    7    3',/,
     &        'CONECT    7    4',/,'CONECT    7    8')
99008 FORMAT (A,I7)
      END
C*==tetsum.f    processed by SPAG 6.70Rc at 21:40 on 19 Dec 2016
      SUBROUTINE TETSUM(WTET,DET,LDET,BRA,IBRA,TW,PHASE,IKCTET,NKTET,
     &                  NTETS,TK,DETGT,DETGTW,DETFE,DETFEW,NELMTMAX,
     &                  NELMT)
C   ********************************************************************
C   *                                                                  *
C   *  perform the BZ-integral via the thedrahedron method             *
C   *  see: Ph. Lambin + J.P. Vigneron  PRB 29, p.3430 (1984)          *
C   *                                                                  *
C   *  WTET   weight for tetrathedra                                   *
C   *  WK     weight of k-point evalated first and then applied        *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:C0,PI,CONST_2PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 DETFEW,DETGTW,PHASE
      INTEGER NELMT,NELMTMAX,NKTET,NTETS
      REAL*8 BRA(NKTET),WTET(NTETS)
      COMPLEX*16 DET(NKTET),DETFE(NKTET),DETGT(NKTET),LDET(NKTET),
     &           TK(NELMTMAX,NKTET),TW(NELMTMAX)
      INTEGER IBRA(NKTET),IKCTET(4,NTETS)
C
C Local variables
C
      REAL*8 DIV,PHJMP
      COMPLEX*16 F1,F2,F3,F4,LT1,LT2,LT3,LT4,T1,T2,T3,T4,TWOPII,W1,W2,
     &           W3,W4,WK(:)
      INTEGER I,ITET,J1,J2,J3,J4,K
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE WK
C
      ALLOCATE (WK(NKTET))
      WK(1:NKTET) = C0
C
      TWOPII = DCMPLX(0.0D0,CONST_2PI)
C
      DO I = 1,NELMT
         TW(I) = C0
      END DO
      PHASE = C0
      DETGTW = C0
      DETFEW = C0
C
      DO K = 1,NKTET
         LDET(K) = LOG(DET(K)) + TWOPII*BRA(K)
         IBRA(K) = INT(DIMAG(LDET(K))/PI+100.5D0)
         IBRA(K) = INT(DBLE(IBRA(K))-100.0D0)
      END DO
C
C-----------------------------------------------------------------------
C      loop over the tetrahedra to find weight  wk  of each  k-point
C-----------------------------------------------------------------------
C
      DO ITET = 1,NTETS
         J1 = IKCTET(1,ITET)
         J2 = IKCTET(2,ITET)
         J3 = IKCTET(3,ITET)
         J4 = IKCTET(4,ITET)
C
         PHJMP = MAX(IBRA(J1),IBRA(J2),IBRA(J3),IBRA(J4))
     &           - MIN(IBRA(J1),IBRA(J2),IBRA(J3),IBRA(J4))
C
         DIV = 1.D0/MAX(1D0,PHJMP)
C
         IF ( PHJMP.LT.0.01D0 ) THEN
            LT1 = LDET(J1)
            LT2 = LDET(J2)
            LT3 = LDET(J3)
            LT4 = LDET(J4)
C
            T1 = DET(J1)
            T2 = DET(J2)
            T3 = DET(J3)
            T4 = DET(J4)
         ELSE
            LT1 = LDET(J1)*DIV
            LT2 = LDET(J2)*DIV
            LT3 = LDET(J3)*DIV
            LT4 = LDET(J4)*DIV
C
            T1 = EXP(LT1)
            T2 = EXP(LT2)
            T3 = EXP(LT3)
            T4 = EXP(LT4)
         END IF
C
         F1 = T1**2/((T1-T2)*(T1-T3)*(T1-T4))
         F2 = T2**2/((T2-T1)*(T2-T3)*(T2-T4))
         F3 = T3**2/((T3-T1)*(T3-T2)*(T3-T4))
         F4 = T4**2/((T4-T1)*(T4-T2)*(T4-T3))
C
         W1 = (F1+F2*T2/(T2-T1)*(LT2-LT1)+F3*T3/(T3-T1)*(LT3-LT1)
     &        +F4*T4/(T4-T1)*(LT4-LT1))*T1*WTET(ITET)
C
         W2 = (F2+F1*T1/(T1-T2)*(LT1-LT2)+F3*T3/(T3-T2)*(LT3-LT2)
     &        +F4*T4/(T4-T2)*(LT4-LT2))*T2*WTET(ITET)
C
         W3 = (F3+F1*T1/(T1-T3)*(LT1-LT3)+F2*T2/(T2-T3)*(LT2-LT3)
     &        +F4*T4/(T4-T3)*(LT4-LT3))*T3*WTET(ITET)
C
         W4 = (F4+F1*T1/(T1-T4)*(LT1-LT4)+F2*T2/(T2-T4)*(LT2-LT4)
     &        +F3*T3/(T3-T4)*(LT3-LT4))*T4*WTET(ITET)
C
         PHASE = PHASE - (W1*LT1+W2*LT2+W3*LT3+W4*LT4)/DIV
         DETGTW = DETGTW + W1*DETGT(J1) + W2*DETGT(J2) + W3*DETGT(J3)
     &            + W4*DETGT(J4)
         DETFEW = DETFEW + W1*DETFE(J1) + W2*DETFE(J2) + W3*DETFE(J3)
     &            + W4*DETFE(J4)
C
         WK(J1) = WK(J1) + W1
         WK(J2) = WK(J2) + W2
         WK(J3) = WK(J3) + W3
         WK(J4) = WK(J4) + W4
C
      END DO
C
C-----------------------------------------------------------------------
C                  sum over the k-points and apply weight WK
C-----------------------------------------------------------------------
C
      DO K = 1,NKTET
         DO I = 1,NELMT
            TW(I) = TW(I) + TK(I,K)*WK(K)
         END DO
      END DO
C
      DO I = 1,NELMT
         TW(I) = TW(I)/DBLE(NTETS)
      END DO
C
      PHASE = PHASE/PI/DBLE(NTETS)
      DETGTW = DETGTW/DBLE(NTETS)
      DETFEW = DETFEW/DBLE(NTETS)
C
      END
C*==tetsor.f    processed by SPAG 6.70Rc at 21:40 on 19 Dec 2016
      SUBROUTINE TETSOR(NFTET,NKTET,KTET,IKCTET,NTETS,IERR,BOA,COA)
C   ********************************************************************
C   *                                                                  *
C   *                BRAVAIS= 4 orthorombic primitive                  *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 BOA,COA
      INTEGER IERR,NFTET,NKTET,NTETS
      INTEGER IKCTET(4,NTETS)
      REAL*8 KTET(3,NKTET)
C
C Local variables
C
      REAL*8 DKX,DKY,DKZ
      INTEGER IC(4,6),IK(8),IKTET,ITETRA,IX,IY,IZ,JC,JT,NPL,NX,NY,NZ
C
C*** End of declarations rewritten by SPAG
C
      DATA IC/1,2,3,7,1,2,6,7,1,5,6,7,1,5,8,7,1,4,8,7,1,4,3,7/
C
      NX = 2*NFTET + 1
C
      DKX = 0.5D0/DBLE(NX-1)
      NY = MAX(2,1+INT((0.5D0/BOA)/DKX))
      NZ = MAX(2,1+INT((0.5D0/COA)/DKX))
      DKY = (0.5D0/BOA)/DBLE(NY-1)
      DKZ = (0.5D0/COA)/DBLE(NZ-1)
C
      NPL = NX*NY
C
      IKTET = 0
      ITETRA = 0
      DO IZ = 1,NZ
         DO IX = 1,NX
            DO IY = 1,NY
               IKTET = IKTET + 1
               IF ( IKTET.GT.NKTET ) THEN
                  IERR = 1
               ELSE
                  KTET(1,IKTET) = DKX*DBLE(IX-1)
                  KTET(2,IKTET) = DKY*DBLE(IY-1)
                  KTET(3,IKTET) = DKZ*DBLE(IZ-1)
               END IF
C
               IF ( IX.NE.NX .AND. IZ.NE.NZ ) THEN
                  IK(1) = IKTET
                  IK(2) = IKTET + IX
                  IK(3) = IKTET + IX + 1
                  IK(4) = IKTET + 1
                  DO JC = 5,8
                     IK(JC) = IK(JC-4) + NPL
                  END DO
C
                  DO JT = 1,6
                     ITETRA = ITETRA + 1
                     IF ( ITETRA.GT.NTETS ) THEN
                        IERR = 1
                     ELSE
                        DO JC = 1,4
                           IKCTET(JC,ITETRA) = IK(IC(JC,JT))
                        END DO
                     END IF
                  END DO
C
               END IF
            END DO
         END DO
      END DO
C
      IF ( IKTET.NE.NKTET ) IERR = 1
      IF ( ITETRA.NE.NTETS ) IERR = 1
      END
C*==tethex.f    processed by SPAG 6.70Rc at 21:40 on 19 Dec 2016
      SUBROUTINE TETHEX(NFTET,NKTET,KTET,IKCTET,NTETS,IERR,COA,BRAVAIS)
C   ********************************************************************
C   *                                                                  *
C   *                BRAVAIS= 8 tetragonal  primitive                  *
C   *                                                                  *
C   *                BRAVAIS=11 hexagonal   primitive                  *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER BRAVAIS,IERR,NFTET,NKTET,NTETS
      REAL*8 COA
      INTEGER IKCTET(4,NTETS)
      REAL*8 KTET(3,NKTET)
C
C Local variables
C
      REAL*8 DKX,DKY,DKZ
      INTEGER IC(4,6),IK(8),IKTET,ITETRA,IX,IY,IZ,JC,JT,NPL,NTP,NX,NZ
C
C*** End of declarations rewritten by SPAG
C
      DATA IC/1,2,3,7,1,2,6,7,1,5,6,7,1,5,8,7,1,4,8,7,1,4,3,7/
C
      NX = 2*NFTET + 1
C
      IF ( BRAVAIS.EQ.11 ) THEN
         DKX = SQRT(1.0D0/3.0D0)/DBLE(NX-1)
         DKY = (1.0D0/3.0D0)/DBLE(NX-1)
         DKZ = (DKX+DKY)/2.0D0
      ELSE
         DKX = 0.5D0/DBLE(NX-1)
         DKY = DKX
         DKZ = DKX
      END IF
C
      NZ = MAX(2,1+INT((0.5D0/COA)/DKZ))
      DKZ = (0.5D0/COA)/DBLE(NZ-1)
      NPL = NX*(NX+1)/2
C
      IKTET = 0
      ITETRA = 0
      DO IZ = 1,NZ
         DO IX = 1,NX
            DO IY = 1,IX
               IKTET = IKTET + 1
               IF ( IKTET.GT.NKTET ) THEN
                  IERR = 1
               ELSE
                  KTET(1,IKTET) = DKX*DBLE(IX-1)
                  KTET(2,IKTET) = DKY*DBLE(IY-1)
                  KTET(3,IKTET) = DKZ*DBLE(IZ-1)
               END IF
C
               IF ( IX.NE.NX .AND. IZ.NE.NZ ) THEN
                  IK(1) = IKTET
                  IK(2) = IKTET + IX
                  IK(3) = IKTET + IX + 1
                  IK(4) = IKTET + 1
                  DO JC = 5,8
                     IK(JC) = IK(JC-4) + NPL
                  END DO
C
                  IF ( IX.EQ.IY ) THEN
                     NTP = 3
                  ELSE
                     NTP = 6
                  END IF
C
                  DO JT = 1,NTP
                     ITETRA = ITETRA + 1
                     IF ( ITETRA.GT.NTETS ) THEN
                        IERR = 1
                     ELSE
                        DO JC = 1,4
                           IKCTET(JC,ITETRA) = IK(IC(JC,JT))
                        END DO
                     END IF
                  END DO
C
               END IF
            END DO
         END DO
      END DO
C
      IF ( IKTET.NE.NKTET ) IERR = 1
      IF ( ITETRA.NE.NTETS ) IERR = 1
      END
C*==tetsic.f    processed by SPAG 6.70Rc at 21:40 on 19 Dec 2016
      SUBROUTINE TETSIC(NFTET,NKTET,KTET,IKCTET,NTETS,IERR)
C   ********************************************************************
C   *                                                                  *
C   *                BRAVAIS= 12 cubic       primitive                 *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IERR,NFTET,NKTET,NTETS
      INTEGER IKCTET(4,NTETS)
      REAL*8 KTET(3,NKTET)
C
C Local variables
C
      REAL*8 DK
      INTEGER I,IC(4,6,4),IK(8),IKTET,ISW,ITETRA,IX,IY,IZ,JC,JT,ND(4),NN
      INTEGER NL
C
C*** End of declarations rewritten by SPAG
C
      DATA IC/1,2,3,7,1,2,6,7,1,5,6,7,1,5,8,7,1,4,8,7,1,4,3,7,1,3,4,8,1,
     &     3,7,8,1,2,3,7,12*0,1,2,3,7,1,2,6,7,1,5,6,7,12*0,1,2,3,7,
     &     20*0/,ND/6,3,3,1/
C
      NL(I) = I*(I+1)/2
C
      NN = 2*NFTET + 1
      DK = 1D0/DBLE(NN-1)
C
      IKTET = 0
      ITETRA = 0
      DO IZ = 1,NN
         DO IX = IZ,NN
            DO IY = IZ,IX
               IKTET = IKTET + 1
               IF ( IKTET.GT.NKTET ) THEN
                  IERR = 1
               ELSE
                  KTET(1,IKTET) = DK*DBLE(IX-1)
                  KTET(2,IKTET) = DK*DBLE(IY-1)
                  KTET(3,IKTET) = DK*DBLE(IZ-1)
               END IF
C
               IF ( IX.NE.NN .AND. IZ.NE.NN ) THEN
                  IK(1) = IKTET
                  IK(2) = IKTET + IX - IZ + 1
                  IK(3) = IKTET + IX - IZ + 2
                  IK(4) = IKTET + 1
                  IK(5) = IKTET + NL(NN-IZ+1) - NL(IX-IZ) - (IY-IZ+1)
     &                    + NL(IX-IZ-1) + (IY-IZ)
                  IK(6) = IK(5) + IX - IZ
                  IK(7) = IK(6) + 1
                  IK(8) = IK(5) + 1
C
                  ISW = 1
                  IF ( IY.EQ.IZ ) ISW = ISW + 1
                  IF ( IY.EQ.IX ) ISW = ISW + 2
C
C     --- all are simple cases: four corners are tabulated in "ic"
C
                  DO JT = 1,ND(ISW)
                     ITETRA = ITETRA + 1
                     IF ( ITETRA.GT.NTETS ) THEN
                        IERR = 1
                     ELSE
                        DO JC = 1,4
                           IKCTET(JC,ITETRA) = IK(IC(JC,JT,ISW))
                        END DO
                     END IF
                  END DO
C
               END IF
            END DO
         END DO
      END DO
C
      IF ( IKTET.NE.NKTET ) IERR = 1
      IF ( ITETRA.NE.NTETS ) IERR = 1
      END
C*==tetfcc.f    processed by SPAG 6.70Rc at 21:40 on 19 Dec 2016
      SUBROUTINE TETFCC(NFTET,NKTET,KTET,IKCTET,NTETS,IWK,NW1,NW2,NW3,
     &                  IERR)
C   ********************************************************************
C   *                                                                  *
C   *                BRAVAIS=13 cubic       face-centered              *
C   *                                                                  *
C   *  NFTET       1    2    3    4    5     6                         *
C   * ------------------------------------------                       *
C   *  KPOINT     20   89   240  505  916  1505                        *
C   *  TETS       32  256   864 2048 4000                              *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IERR,NFTET,NKTET,NTETS,NW1,NW2,NW3
      INTEGER IKCTET(4,NTETS),IWK(NW1,NW2,NW3)
      REAL*8 KTET(3,NKTET)
C
C Local variables
C
      REAL*8 DK
      INTEGER I,IC(4,6,10),IKTET,IP,ISW,IT(3,8),ITETRA,IX,IXP,IY,IYP,IZ,
     &        IZP,J,K,ND(10),NE
C
C*** End of declarations rewritten by SPAG
C
      DATA IT/0,0,0,1,0,0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,1,1,1,1/,IC/1,2,4,
     &     8,1,2,6,8,1,5,6,8,1,5,7,8,1,3,7,8,1,3,4,8,1,2,3,5,20*0,1,5,6,
     &     7,1,2,3,6,1,3,6,7,4,2,3,6,4,3,6,7,4*0,1,2,4,8,1,2,6,8,1,5,6,
     &     8,12*0,24*0,1,2,4,6,1,4,5,6,16*0,1,2,4,8,1,3,4,8,1,3,7,8,
     &     12*0,24*0,1,2,4,7,1,3,4,7,16*0,1,2,4,8,20*0/,ND/6,1,5,3,0,2,
     &     3,0,2,1/
C
      NE = 4*NFTET
      DK = 1D0/DBLE(NE)
C
      IKTET = 0
      DO IX = 1,NE + 1
         DO IY = 1,IX
            DO IZ = 1,IY
               IF ( IX+IY+IZ.LE.6*NFTET+3 ) THEN
                  IKTET = IKTET + 1
                  IWK(IX,IY,IZ) = IKTET
                  IF ( IKTET.GT.NKTET ) THEN
                     IERR = 1
                  ELSE
                     KTET(1,IKTET) = DK*DBLE(IX-1)
                     KTET(2,IKTET) = DK*DBLE(IY-1)
                     KTET(3,IKTET) = DK*DBLE(IZ-1)
                  END IF
               END IF
            END DO
         END DO
      END DO
      IF ( IKTET.NE.NKTET ) THEN
         WRITE (6,*) 'IKTET, NKTET:',IKTET,NKTET
         STOP 'in <TETFCC>'
      END IF
C
      ITETRA = 0
      DO IX = 1,NE
         DO IY = 1,IX
            DO IZ = 1,IY
               IF ( IX+IY+IZ.LE.6*NFTET+2 ) THEN
                  K = IWK(IX,IY,IZ)
C
C                                   --- switch to the case 1,...,case 10
                  ISW = 1
                  IF ( IX.EQ.IY ) ISW = ISW + 3
                  IF ( IY.EQ.IZ ) ISW = ISW + 6
                  IF ( ABS((IX+IY+IZ-2)*DK-1.5D0).LT.1D-6 ) ISW = ISW + 
     &                 1
                  IF ( ABS((IX+IY+IZ-1)*DK-1.5D0).LT.1D-6 ) ISW = ISW + 
     &                 2
C
                  IF ( ISW.EQ.5 ) THEN
                     ITETRA = ITETRA + 1
                     IF ( ITETRA.GT.NTETS ) THEN
                        IERR = 1
                     ELSE
C                                            --- look for the cube below
                        IF ( IZ.LE.1 ) GOTO 100
                        IKCTET(1,ITETRA) = K
                        IKCTET(2,ITETRA) = IWK(IX,IY,IZ+1)
                        IKCTET(3,ITETRA) = IWK(IX+1,IY,IZ)
                        IKCTET(4,ITETRA) = IWK(IX+1,IY+1,IZ-1)
                        CYCLE
                     END IF
                  END IF
C
                  IF ( ISW.EQ.8 ) THEN
C                                             --- look for the cube left
                     ITETRA = ITETRA + 1
                     IF ( ITETRA.GT.NTETS ) THEN
                        IERR = 1
                     ELSE
                        IF ( IX.LE.1 ) GOTO 100
                        IKCTET(1,ITETRA) = K
                        IKCTET(2,ITETRA) = IWK(IX+1,IY,IZ)
                        IKCTET(3,ITETRA) = IWK(IX,IY+1,IZ)
                        IKCTET(4,ITETRA) = IWK(IX-1,IY+1,IZ+1)
                        CYCLE
                     END IF
                  END IF
C
C                   --- simple cases: four corners are tabulated in "ic"
                  DO I = 1,ND(ISW)
                     ITETRA = ITETRA + 1
                     IF ( ITETRA.GT.NTETS ) THEN
                        IERR = 1
                     ELSE
                        DO J = 1,4
                           IP = IC(J,I,ISW)
                           IXP = IX + IT(1,IP)
                           IYP = IY + IT(2,IP)
                           IZP = IZ + IT(3,IP)
                           IKCTET(J,ITETRA) = IWK(IXP,IYP,IZP)
                        END DO
                     END IF
                  END DO
               END IF
            END DO
         END DO
      END DO
C
      IF ( ITETRA.NE.NTETS ) THEN
         WRITE (6,*) 'ITETRA, NTETS:',ITETRA,NTETS
         STOP 'in <TETFCC>'
      END IF
      RETURN
C
 100  CONTINUE
      WRITE (6,'(A,I4)') ' STOP IN <TETFCC>...FAILS FOR NTETS= ',NTETS
      STOP
      END
C*==tetbcc.f    processed by SPAG 6.70Rc at 21:40 on 19 Dec 2016
      SUBROUTINE TETBCC(NFTET,NKTET,KTET,IKCTET,NTETS,IWK,NW1,NW2,NW3,
     &                  IERR)
C   ********************************************************************
C   *                                                                  *
C   *                BRAVAIS=14 cubic       body-centered              *
C   *                                                                  *
C   *    NFTET       1   2   3    4    5    6                          *
C   * ---------------------------------------                          *
C   *  KPOINT    14  55  140  285  506  819                            *
C   *  TETS      16 128  432 1024 2000                                 *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IERR,NFTET,NKTET,NTETS,NW1,NW2,NW3
      INTEGER IKCTET(4,NTETS),IWK(NW1,NW2,NW3)
      REAL*8 KTET(3,NKTET)
C
C Local variables
C
      REAL*8 DK
      INTEGER I,IC(4,6,6),IKTET,IP,ISW,IT(3,8),ITETRA,IX,IXP,IY,IYP,IZ,
     &        IZP,J,ND(6),NE
C
C*** End of declarations rewritten by SPAG
C
      DATA IT/0,0,0,1,0,0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,1,1,1,1/IC/1,2,4,
     &     8,1,2,6,8,1,5,6,8,1,5,7,8,1,3,7,8,1,3,4,8,1,2,4,8,1,3,4,8,1,
     &     3,7,8,12*0,1,2,4,8,1,2,6,8,1,5,6,8,12*0,1,2,4,8,20*0,1,2,3,7,
     &     1,2,6,7,1,5,6,7,12*0,1,2,3,7,20*0/,ND/6,3,3,1,3,1/
C
      NE = 4*NFTET
      DK = 1D0/DBLE(NE)
C
      IKTET = 0
      DO IX = 1,NE + 1
         DO IY = 1,IX
            DO IZ = 1,IY
               IF ( IX+IY.LE.4*NFTET+2 ) THEN
                  IKTET = IKTET + 1
                  IWK(IX,IY,IZ) = IKTET
                  IF ( IKTET.GT.NKTET ) THEN
                     IERR = 1
                  ELSE
                     KTET(1,IKTET) = DK*DBLE(IX-1)
                     KTET(2,IKTET) = DK*DBLE(IY-1)
                     KTET(3,IKTET) = DK*DBLE(IZ-1)
                  END IF
               END IF
            END DO
         END DO
      END DO
      IF ( IKTET.NE.NKTET ) THEN
         WRITE (6,*) 'IKTET, NKTET:',IKTET,NKTET
         STOP 'in <TETBCC>'
      END IF
C
      ITETRA = 0
      DO IX = 1,NE
         DO IY = 1,IX
            DO IZ = 1,IY
               IF ( IX+IY.LE.4*NFTET+1 ) THEN
C                                    --- switch to the case 1,...,case 6
                  ISW = 1
                  IF ( IX.EQ.IY ) ISW = ISW + 2
                  IF ( IY.EQ.IZ ) ISW = ISW + 1
                  IF ( ABS(IX+IY-1-NE).EQ.0 ) ISW = ISW + 4
C
C           --- all are simple cases: four corners are tabulated in "ic"
C
                  DO I = 1,ND(ISW)
                     ITETRA = ITETRA + 1
                     IF ( ITETRA.GT.NTETS ) THEN
                        IERR = 1
                     ELSE
                        DO J = 1,4
                           IP = IC(J,I,ISW)
                           IXP = IX + IT(1,IP)
                           IYP = IY + IT(2,IP)
                           IZP = IZ + IT(3,IP)
                           IKCTET(J,ITETRA) = IWK(IXP,IYP,IZP)
                        END DO
                     END IF
                  END DO
               END IF
C
            END DO
         END DO
      END DO
C
      IF ( ITETRA.NE.NTETS ) THEN
         WRITE (6,*) 'ITETRA, NTETS:',ITETRA,NTETS
         STOP 'in <TETBCC>'
      END IF
      END
C
