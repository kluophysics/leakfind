C*==scfbroypt2.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SCFBROYPT2(IPRINT,XINP,XOUT,FM,FM1,SM,SM1,UI,VI,AM,BM,
     &                      ALPHA,ITDEPT,ITDEPTMAX,IMIX,IOBROY,IPF,MIT,
     &                      NMAP,NMAPMAX)
C   ********************************************************************
C   *                                                                  *
C   * imix :                                                           *
C   *   3   broyden's            f i r s t  m e t h o d                *
C   *   4   broyden's          s e c o n d  m e t h o d                *
C   *   5   use ibroy=1 and save jacobian accumulated up to itdept-1   *
C   *       to update potential                                        *
C   *   6   use ibroy=2 and save jacobian accumulated up to itdept-1   *
C   *       to update potential                                        *
C       7   anderson's     g e n e r a l i z e d   m e t h o d
C   *                                                                  *
C   * implemented here according to notes of s.b.                      *
C   * broyden iteration scheme following the papers of :               *
C   * srivastava, j. phys. , 17 (1984) , pp l317                       *
C   * c.g. broyden in math.comput., 19 , pp 577, 1965                  *
C   * c.g. broyden in ibid, 21 ,pp 368 ,1967                           *
C   * the method has been generalized to include a metric.the          *
C   * definition of the necessary inner products are similar to the    *
C   * discription given in the notes of m.weinert.the algorithm        *
C   * discribed in the paper srivastava  has been simplified           *
C   * ( see notes of s.b.)                                             *
C   * the files ui,vi are stored on file iobroy=70                     *
C   * broyden's update treats charge and spin on the same footing      *
C   *              s. bluegel , kfa , may 1987                         *
C   * the anderson method (d.g. anderson, j. acm 12, 547 (1964)) has   *
C   * been generalized and reformulated as an improvement of broyden's *
C   * second method. successive linesearch is replaced by successive   *
C   * search on hyperplanes. ( see notes of s.b. )                     *
C   *              s. bluegel , issp , july 1989                       *
C   *------------------------------------------------------------------*
C   * modified by h. ebert   july/august 1994                          *
C   ********************************************************************
      IMPLICIT NONE
C*--SCFBROYPT238
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
C
C Dummy arguments
C
      REAL*8 ALPHA
      INTEGER IMIX,IOBROY,IPF,IPRINT,ITDEPT,ITDEPTMAX,MIT,NMAP,NMAPMAX
      REAL*8 AM(2:ITDEPTMAX-1),BM(2:ITDEPTMAX-1),FM(NMAPMAX),
     &       FM1(NMAPMAX),SM(NMAPMAX),SM1(NMAPMAX),UI(NMAPMAX,2:3),
     &       VI(NMAPMAX,2:3),XINP(NMAPMAX),XOUT(NMAPMAX)
C
C Local variables
C
      REAL*8 CMM,RMIXIV,SMNORM,VMDENO,VMNORM,WIT(2:200)
      REAL*8 DDOT
      INTEGER IBROY,IJ,IT
      SAVE WIT
      EXTERNAL DAXPY,DDOT,DSCAL
C
C*** End of declarations rewritten by SPAG
C
C
C
C*** End of declarations rewritten by SPAG
C
      IBROY = IMIX - 2
      IF ( IBROY.LE.0 .OR. IBROY.GT.5 ) STOP 'ibroyd ?'
C
      IF ( (IBROY.LE.2 .OR. IBROY.EQ.5) .AND. MIT.GT.ITDEPT ) MIT = 1
      IF ( (IBROY.EQ.3 .OR. IBROY.EQ.4) .AND. MIT.GE.ITDEPT )
     &     MIT = ITDEPT
      IF ( IBROY.EQ.3 .AND. MIT.LT.ITDEPT ) IBROY = 1
      IF ( IBROY.EQ.4 .AND. MIT.LT.ITDEPT ) IBROY = 2
C
      IF ( IPRINT.GT.0 ) THEN
         IF ( IBROY.EQ.1 ) WRITE (IPF,FMT=99004) '1st'
         IF ( IBROY.EQ.2 ) WRITE (IPF,FMT=99004) '2nd'
         IF ( IBROY.EQ.3 ) WRITE (IPF,FMT=99005) '1st'
         IF ( IBROY.EQ.4 ) WRITE (IPF,FMT=99005) '2nd'
         IF ( IBROY.EQ.5 ) WRITE (IPF,FMT=99006)
      END IF
C
      RMIXIV = ONE/ALPHA
C
C---->  the comming block is activated only one iteration before
C        broyden iteration scheme is used
C---->  set up of : sm1 = pot(1) ; fm1=fm[1]=f(pot(1)) - pot(1) ;
C
      IF ( MIT.EQ.1 ) THEN
C
         DO IJ = 1,NMAP
            SM1(IJ) = XINP(IJ)
            FM1(IJ) = RMIXIV*(XOUT(IJ)-SM1(IJ))
         END DO
C
         MIT = MIT + 1
C
      ELSE
C
         DO IJ = 1,NMAP
C----> map pot(m) of all mt-spheres into one single vector
C
            SM(IJ) = XINP(IJ)
C
C----> map f[m] = f(m) - pot(m) = f(pot(m)) - pot(m) of all mt-spheres
C      into one single vector
            FM(IJ) = RMIXIV*(XOUT(IJ)-SM(IJ))
C
C----> calculate  sm = pot(m) - pot(m-1)
C----> calculate dfm = f[m] - f[m-1]
C
            SM1(IJ) = SM(IJ) - SM1(IJ)
            FM1(IJ) = FM(IJ) - FM1(IJ)
C
C----> loop to generate u[m] = u(ij,mit)
C
            UI(IJ,3) = ALPHA*FM1(IJ) + SM1(IJ)
         END DO
C
         REWIND IOBROY
         DO IT = 2,MIT - 1
            READ (IOBROY,ERR=50) (UI(IJ,2),IJ=1,NMAP),
     &                           (VI(IJ,2),IJ=1,NMAP)
            AM(IT) = DDOT(NMAP,FM1,1,VI(1,2),1)
            CALL DAXPY(NMAP,-AM(IT),UI(1,2),1,UI(1,3),1)
         END DO
C
C----> print amj = the importance of the history of ui
C
         IF ( IPRINT.GT.0 ) WRITE (IPF,99001) MIT - 1,
     &                             (AM(IT),IT=2,MIT-1)
C
C
         IF ( IBROY.EQ.1 ) THEN
C-------->     b r o y d e n ' s   f i r s t   m e t h o d
C
C
C----> calculate dsmnorm
C
            SMNORM = ZERO
            DO IJ = 1,NMAP
               SMNORM = SMNORM + SM1(IJ)*SM1(IJ)
C
C----> loop to generate x[m] = x(ij,mit)
C
               VI(IJ,3) = ALPHA*SM1(IJ)
            END DO
C
            REWIND IOBROY
            DO IT = 2,MIT - 1
               READ (IOBROY,ERR=50) (UI(IJ,2),IJ=1,NMAP),
     &                              (VI(IJ,2),IJ=1,NMAP)
               BM(IT) = DDOT(NMAP,SM1,1,UI(1,2),1)
               CALL DAXPY(NMAP,-BM(IT),VI(1,2),1,VI(1,3),1)
            END DO
C
C----> complete the evaluation of x[m]
C
            VMDENO = DDOT(NMAP,SM1,1,UI(1,3),1) - SMNORM
            IF ( ABS(VMDENO).LT.1D-70 ) THEN
               STOP 'bry1sn'
C
            ELSE
C
               CALL DSCAL(NMAP,ONE/VMDENO,VI(1,3),1)
C
C----> print bmj = the importance of the history of vi
C
               IF ( IPRINT.GT.0 ) WRITE (IPF,99002) MIT - 1,
     &              (BM(IT),IT=2,MIT-1)
            END IF
C
         ELSE IF ( IBROY.EQ.2 ) THEN
C-------->     b r o y d e n ' s   s e c o n d    m e t h o d
C
C----> calculate x[m]
C
            DO IJ = 1,NMAP
               VI(IJ,3) = FM1(IJ)
            END DO
C
C----> calculate #xm# and normalize x[m]
C
            VMNORM = DDOT(NMAP,VI(1,3),1,FM1,1)
            CALL DSCAL(NMAP,ONE/VMNORM,VI(1,3),1)
C
         ELSE IF ( IBROY.EQ.5 ) THEN
C-------->     g e n e r a l i z e d   a n d e r s o n   m e t h o d
C
C----> calculate v[m]
C
            DO IJ = 1,NMAP
               VI(IJ,3) = FM1(IJ)
            END DO
            REWIND IOBROY
            DO IT = 2,MIT - 1
               READ (IOBROY,ERR=50) (UI(IJ,2),IJ=1,NMAP),
     &                              (VI(IJ,2),IJ=1,NMAP)
Cabc              READ (IOBROY) (UI2(IJ),IJ=1,IMAP), (VI2(IJ),IJ=1,IMAP),WIT
C
               CALL DAXPY(NMAP,-AM(IT)*WIT(IT),VI(1,2),1,VI(1,3),1)
            END DO
C
C----> complete the evaluation of v[m]
C
            VMDENO = DDOT(NMAP,FM1,1,VI(1,3),1)
C
            IF ( ABS(VMDENO).LT.1D-70 ) STOP '<scfbroypt2> VMDENO !'
C
            CALL DSCAL(NMAP,ONE/VMDENO,VI(1,3),1)
C
C----> save wit(mit) for next iteration
C
            WIT(MIT) = VMDENO
C
         END IF
C
C----> update f[m-1] = f[m]  ; pot(m) = pot(m-1)
C
         DO IJ = 1,NMAP
            FM1(IJ) = FM(IJ)
            SM1(IJ) = SM(IJ)
         END DO
C
C----> write on ssd u(ij,mit) and x(ij,mit)
C
         IF ( IBROY.LE.2 .OR. IBROY.EQ.5 ) WRITE (IOBROY,ERR=50)
     &        (UI(IJ,3),IJ=1,NMAP),(VI(IJ,3),IJ=1,NMAP)
C
C----> calculate cmm
C
         IF ( IBROY.LE.2 .OR. IBROY.EQ.5 ) THEN
            CMM = DDOT(NMAP,FM,1,VI(1,3),1)
            IF ( IPRINT.GT.0 ) WRITE (IPF,99003) CMM
C
         ELSE
            CMM = ZERO
         END IF
C
C----> update pot(m+1)
C
         CALL DAXPY(NMAP,ONE-CMM,UI(1,3),1,SM,1)
C
C----> map solution back into each mt-sphere
C
         DO IJ = 1,NMAP
            XOUT(IJ) = SM(IJ)
         END DO
C
         MIT = MIT + 1
         RETURN
C
 50      CONTINUE
         STOP 'iobroy'
C
      END IF
99001 FORMAT (5x,' amj ---> j=2,',i3,:,/,(9x,1p,7D10.2))
99002 FORMAT (5x,' bmj ---> j=2,',i3,:,/,(9x,1p,7D10.2))
99003 FORMAT (5x,' cmm = ',1p,d12.4)
C
99004 FORMAT (' broyden''s ',A3,' method used ')
99005 FORMAT ('jacobian is fixed after broyden''s ',A3,' method used')
99006 FORMAT (' andersons''s method used ')
      END
