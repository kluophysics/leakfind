C*==eos.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      PROGRAM EOS
      USE MOD_LATTICE,ONLY:ABAS
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_SCF,ONLY:ALAT_TAB,EFERMI_TAB,ETOT_TAB,MUEORB_TAB,
     &    MUESPN_TAB,VOL_TAB
      USE MOD_FILES,ONLY:IFILSCLALAT,DATSET,LSYSTEM,SYSTEM,LDATSET
      IMPLICIT NONE
C*--EOS9
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      REAL*8 ALAT_IN(:),EFERMI_IN(:),ETOT_IN(:),MUEORB_IN(:),
     &       MUESPN_IN(:),VOL_IN(:),XX
      INTEGER ISCL,NSCL,NSCLIN
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE ALAT_IN,EFERMI_IN,ETOT_IN
      ALLOCATABLE MUEORB_IN,MUESPN_IN,VOL_IN
C
C#                         HEADER OF DATA FILE
C#  SCF-runs for
C#  DATSET = Ni_SCF-NL4
C#  SYSTEM = Ni
C#  IREL   = 1
C#  ABAS   =     0.00000000    0.50000000    0.50000000
C#               0.50000000    0.00000000    0.50000000
C#               0.50000000    0.50000000    0.00000000
C#  scaling the lattice parameter ALAT =     6.65819984
C#  using   SCL1 = 0.950   SCL2 = 1.050   NSCL =  9
C#
C#     ALAT    EFERMI    TOTDOS    MUESPN    MUEORB ....
C
      READ (5,99002) DATSET,SYSTEM,IREL,ABAS,NSCL
      LDATSET = LEN_TRIM(DATSET)
      LSYSTEM = LEN_TRIM(SYSTEM)
C
      WRITE (6,99002) DATSET(1:LDATSET),SYSTEM(1:LSYSTEM),IREL,ABAS,NSCL
C
      ALLOCATE (ALAT_IN(NSCL),EFERMI_IN(NSCL),ETOT_IN(NSCL))
      ALLOCATE (MUEORB_IN(NSCL),MUESPN_IN(NSCL),VOL_IN(NSCL))
C
      DO ISCL = 1,NSCL
C
         READ (5,*,ERR=100,END=100) ALAT_IN(ISCL),EFERMI_IN(ISCL),XX,
     &                              MUESPN_IN(ISCL),MUEORB_IN(ISCL),
     &                              ETOT_IN(ISCL)
         WRITE (6,99001) ALAT_IN(ISCL),EFERMI_IN(ISCL),MUESPN_IN(ISCL),
     &                   MUEORB_IN(ISCL),ETOT_IN(ISCL)
         NSCLIN = ISCL
      END DO
C
 100  CONTINUE
      NSCL = NSCLIN
C
      ALLOCATE (ALAT_TAB(NSCL),EFERMI_TAB(NSCL),ETOT_TAB(NSCL))
      ALLOCATE (MUEORB_TAB(NSCL),MUESPN_TAB(NSCL),VOL_TAB(NSCL))
C
      ALAT_TAB(1:NSCL) = ALAT_IN(1:NSCL)
      EFERMI_TAB(1:NSCL) = EFERMI_IN(1:NSCL)
      ETOT_TAB(1:NSCL) = ETOT_IN(1:NSCL)
      MUEORB_TAB(1:NSCL) = MUEORB_IN(1:NSCL)
      MUESPN_TAB(1:NSCL) = MUESPN_IN(1:NSCL)
C
      CALL SCFSCALE_ALAT_EOS(NSCL)
C
      STOP 'FERTIG'
99001 FORMAT (12F16.6)
99002 FORMAT (//,12X,A,/,12X,A,/,11X,I2,/,3(12X,3F14.8,/),/,47X,I3,//)
99003 FORMAT ('#',/,'#  SCF-runs for ',/,'#  DATSET = ',A,/,
     &        '#  SYSTEM = ',A,/,'#  ABAS   = ',3F14.8,/,
     &        2('#',11X,3F14.8,/),
     &        '#  scaling the lattice parameter ALAT = ',F14.8,/,
     &        '#  using   SCL1 =',F6.3,'   SCL2 =',F6.3,'   NSCL =',I3,
     &        /,'#',/,
     &        '#     ALAT    EFERMI    TOTDOS    MUESPN    MUEORB',
     &        '                ETOT    SCLNOS      RMSV      RMSB')
      END
C*==xmgrhead.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE XMGRHEAD(DATSET,LDATSET,IDENT1,LIDENT1,IDENT2,LIDENT2,
     &                    FILNAM,LFNMAX,LFN,IFIL,NGRAPH,XMIN,KXMIN,XMAX,
     &                    KXMAX,YMIN1,KYMIN1,YMAX1,KYMAX1,YMIN2,KYMIN2,
     &                    YMAX2,KYMAX2,XTXT,LXTXT,YTXT1,LYTXT1,YTXT2,
     &                    LYTXT2,TITLE,LTITLE,SUBTITLE,LSUBTITLE,
     &                    KNOHEADER)
C **********************************************************************
C *                                                                    *
C *  DATSET    LDATSET     dataset name and length of string           *
C *  IDENT1    LIDENT1     add. identifier and length of string > 0 !! *
C *  IDENT2    LIDENT2     add. identifier and length of string        *
C *  FILNAM    LFNMAX LFN  file name to be created, its max. and       *
C *                        actual length on exit                       *
C *  IFIL                  chanel for xmgr-file                        *
C *  NGRAPH                number of graphs   1  or  2                 *
C *  XMIN      KXMIN       x-start and key to fix (0) or float (1)     *
C *  XMAX      KXMAX       x-end   and key to fix (0) or float (1)     *
C *  YMIN1     KYMIN1      y-start and key to fix (0) or float (1)     *
C *  YMAX1     KYMAX1      y-end   and key ...  for graph 1            *
C *  YMIN2     KYMIN2      y-start and key ...  for graph 2            *
C *  YMAX2     KYMAX2      y-end   and key ...  used if NGRAPH=2       *
C *  XTXT      LXTXT       text for x-axis label                       *
C *  YTXT1     LYTXT1      text for y-axis label for graph 1           *
C *  YTXT2     LYTXT2      text for y-axis label for graph 2           *
C *  TITLE     LTITLE      title and length of string                  *
C *  SUBTITLE  LSUBTITLE   subtitle and length of string               *
C *  KNOHEADER             key to get (false) or suppress (true) title *
C *                                                                    *
C *  NOTE:  only  FILNAM and LFN will be changed on exit               *
C *         apart of IDENT1 any text may have length 0                 *
C *                                                                    *
C *         use '!' instead of '\' for the xmgr escape character       *
C *         in the text variables --- see <XMGRWRITE>, where also      *
C *         the macro list is stored  !UP, !DN,                        *
C *                                                                    *
C **********************************************************************
C
      IMPLICIT NONE
C*--XMGRHEAD129
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*3 CHARSIZE
      PARAMETER (CHARSIZE='1.2')
C
C Dummy arguments
C
      CHARACTER*(*) DATSET,FILNAM,IDENT1,IDENT2,SUBTITLE,TITLE,XTXT,
     &              YTXT1,YTXT2
      INTEGER IFIL,KXMAX,KXMIN,KYMAX1,KYMAX2,KYMIN1,KYMIN2,LDATSET,LFN,
     &        LFNMAX,LIDENT1,LIDENT2,LSUBTITLE,LTITLE,LXTXT,LYTXT1,
     &        LYTXT2,NGRAPH
      LOGICAL KNOHEADER
      REAL*8 XMAX,XMIN,YMAX1,YMAX2,YMIN1,YMIN2
C
C Local variables
C
      REAL*8 DX,DY1,DY2,XA,XB,YA1,YA2,YB1,YB2
      INTEGER IG
C
C*** End of declarations rewritten by SPAG
C
C ----------------------------------------------------------------------
C         create filename   FILNAM(1:LFN)//'.agr'   and open
C ----------------------------------------------------------------------
C
      CALL XMGRFNAM(DATSET,LDATSET,IDENT1,LIDENT1,IDENT2,LIDENT2,'.agr',
     &              4,FILNAM,LFN,LFNMAX)
C
      OPEN (IFIL,FILE=FILNAM(1:LFN))
C
C ----------------------------------------------------------------------
C                       set range  parameters
C ----------------------------------------------------------------------
C
      CALL XMGRTICKS(XMIN,XA,KXMIN,XMAX,XB,KXMAX,0.0D0,DX,1)
      CALL XMGRTICKS(YMIN1,YA1,KYMIN1,YMAX1,YB1,KYMAX1,0.0D0,DY1,1)
      IF ( NGRAPH.EQ.2 ) THEN
         CALL XMGRTICKS(YMIN2,YA2,KYMIN2,YMAX2,YB2,KYMAX2,0.0D0,DY2,1)
         DY1 = 2*DY1
         DY2 = 2*DY2
      END IF
C
C ----------------------------------------------------------------------
C                       write head of xmgr - file
C ----------------------------------------------------------------------
C
      WRITE (IFIL,99008)
      WRITE (IFIL,99009) '@default linewidth 2.5'
      WRITE (IFIL,99009) '@default linestyle 1'
      WRITE (IFIL,99009) '@default color 1'
      WRITE (IFIL,99009) '@default font 0'
      WRITE (IFIL,99009) '@default char size '//CHARSIZE
      WRITE (IFIL,99009) '@default symbol size '//CHARSIZE
      WRITE (IFIL,99009) '@page background fill off'
C
      DO IG = 0,(NGRAPH-1)
         WRITE (IFIL,99002) IG,IG
         WRITE (IFIL,99001)
         WRITE (IFIL,99009) '@    frame linewidth 2.5'
         WRITE (IFIL,99009) '@    xaxis  label font 0'
         WRITE (IFIL,99009) '@    xaxis  label char size '//CHARSIZE
         WRITE (IFIL,99009) '@    xaxis  ticklabel font 0'
         WRITE (IFIL,99009) '@    xaxis  ticklabel char size '//CHARSIZE
         WRITE (IFIL,99009) '@    xaxis  bar linewidth 2.5'
         WRITE (IFIL,99009) '@    xaxis  tick major linewidth 2.5'
         WRITE (IFIL,99009) '@    xaxis  tick minor linewidth 2.5'
         WRITE (IFIL,99009) '@    yaxis  label font 0'
         WRITE (IFIL,99009) '@    yaxis  label char size '//CHARSIZE
         WRITE (IFIL,99009) '@    yaxis  label place spec'
         WRITE (IFIL,99009) '@    yaxis  label place 0.000000, 0.100000'
         WRITE (IFIL,99009) '@    yaxis  ticklabel font 0'
         WRITE (IFIL,99009) '@    yaxis  ticklabel char size '//CHARSIZE
         WRITE (IFIL,99009) '@    yaxis  bar linewidth 2.5'
         WRITE (IFIL,99009) '@    yaxis  tick major linewidth 2.5'
         WRITE (IFIL,99009) '@    yaxis  tick minor linewidth 2.5'
C
C ----------------------------------------------------------- min .. max
         IF ( IG.EQ.0 ) THEN
            IF ( NGRAPH.EQ.2 ) YB1 = YB1 - DY2*0.0001D0
            WRITE (IFIL,99003) 'x',XA,'x',XB
            WRITE (IFIL,99003) 'y',YA1,'y',YB1
         ELSE
            WRITE (IFIL,99007) 'x'
            WRITE (IFIL,99003) 'x',XA,'x',XB
            WRITE (IFIL,99003) 'y',YA2,'y',YB2
         END IF
C
C ------------------------------------------------------------ tick step
         WRITE (IFIL,99004) 'x',DX
         WRITE (IFIL,99005) 'x',DX
         IF ( IG.EQ.0 ) THEN
            WRITE (IFIL,99004) 'y',DY1,'y',DY1
            WRITE (IFIL,99005) 'y',DY1,'y',DY1
         ELSE
            WRITE (IFIL,99004) 'y',DY2,'y',DY2
            WRITE (IFIL,99005) 'y',DY2,'y',DY2
         END IF
C
C ---------------------------------------------------------- view window
         WRITE (IFIL,99006) 'x',0.15,'x',0.85
         IF ( IG.NE.0 ) THEN
            WRITE (IFIL,99006) 'y',0.50,'y',0.85
         ELSE IF ( NGRAPH.NE.2 ) THEN
            WRITE (IFIL,99006) 'y',0.15,'y',0.85
         ELSE
            WRITE (IFIL,99006) 'y',0.15,'y',0.50
         END IF
C
C ------------------------------------------------------------ axis text
         IF ( IG.EQ.0 ) THEN
            CALL XMGRWRITE(IFIL,'x','axis  label ',XTXT,LXTXT,'\\')
            CALL XMGRWRITE(IFIL,'y','axis  label ',YTXT1,LYTXT1,'\\')
         ELSE
            CALL XMGRWRITE(IFIL,'y','axis  label ',YTXT2,LYTXT2,'\\')
         END IF
C
C ---------------------------------------------------------------- title
         IF ( (.NOT.KNOHEADER) .AND. (IG.EQ.(NGRAPH-1)) ) THEN
            CALL XMGRWRITE(IFIL,' ','title',TITLE,LTITLE,'\\')
            CALL XMGRWRITE(IFIL,' ','subtitle',SUBTITLE,LSUBTITLE,'\\')
         END IF
C
      END DO
C
99001 FORMAT ('@  zeroxaxis  bar on')
99002 FORMAT ('@g',i1,' on',/,'@with g',I1)
99003 FORMAT ('@  world ',A,'min ',E18.9,/,'@  world ',A,'max ',E18.9)
99004 FORMAT ('@  ',A,'axis  tick major ',E18.9)
99005 FORMAT ('@  ',A,'axis  tick minor ',E18.9)
99006 FORMAT ('@  view  ',A,'min ',E18.9,/,'@  view  ',A,'max ',E18.9)
99007 FORMAT ('@  ',A,'axis  ticklabel off')
99008 FORMAT ('# Grace project file',/,'#',/,'@version 50005',/,
     &        '@page size 600, 600')
99009 FORMAT (A)
      END
C*==xmgrtable.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE XMGRTABLE(IG,IS,X,Y,W,N,IFIL)
C **********************************************************************
C *                                                                    *
C *   write out the data for graph  IG  and set  IS                    *
C *   W is used as weight parameter for  Y  normally set to  1         *
C *                                                                    *
C **********************************************************************
C
      IMPLICIT NONE
C*--XMGRTABLE293
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IFIL,IG,IS,N
      REAL*8 W
      REAL*8 X(N),Y(N)
C
C Local variables
C
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      IF ( IS.LT.10 ) THEN
         WRITE (IFIL,99001) IG,IS
      ELSE
         WRITE (IFIL,99002) IG,IS
      END IF
C
      DO I = 1,N
         WRITE (IFIL,99004) X(I),Y(I)*W
      END DO
C
      WRITE (IFIL,99003)
C
99001 FORMAT ('@target G',I1,'.S',I1,/,'@type xy')
99002 FORMAT ('@target G',I1,'.S',I2,/,'@type xy')
99003 FORMAT ('&')
99004 FORMAT (2E18.10)
C
      END
C*==xmgrwrite.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE XMGRWRITE(IFIL,T1,T2,TXT,LTXT,ESC)
C **********************************************************************
C *                                                                    *
C *  write the text strings  T1, T2 and TXT  - modify TXT if necessary *
C *                                                                    *
C *  when calling XMGRWRITE  ESC  should be set to '\\'                *
C *  for any compiler / machine this should be converted to            *
C *  have the escape character  '\'  in the xmgr-file at the end       *
C *                                                                    *
C **********************************************************************
C
      IMPLICIT NONE
C*--XMGRWRITE351
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER LSTRMAX
      PARAMETER (LSTRMAX=80)
C
C Dummy arguments
C
      CHARACTER*1 ESC
      INTEGER IFIL,LTXT
      CHARACTER*(*) T1,T2,TXT
C
C Local variables
C
      INTEGER L
      CHARACTER*(LSTRMAX) T
C
C*** End of declarations rewritten by SPAG
C
      IF ( LTXT.LE.0 ) RETURN
C
      T = TXT
      L = LTXT
C
C ----------------------------------------------- reset to original font
      CALL XMGRSUBS(T,L,'!F',2,ESC//'f{}',4)
C
C ------------------------------------ replace macros for spin UP and DN
      CALL XMGRSUBS(T,L,'!UP',3,ESC//'x'//ESC//'c-'//ESC//'C'//ESC//
     &              'f{0}',12)
C
      CALL XMGRSUBS(T,L,'!DN',3,ESC//'x'//ESC//'c/'//ESC//'C'//ESC//
     &              'f{0}',12)
C
C     -------------------------- replace macros for arrow left and right
      CALL XMGRSUBS(T,L,'!AL',3,ESC//'x'//ESC//'c,'//ESC//'C'//ESC//
     &              'f{0}',12)
C
      CALL XMGRSUBS(T,L,'!AR',3,ESC//'x'//ESC//'c.'//ESC//'C'//ESC//
     &              'f{0}',12)
C
C --------------------------- replace '!' by the proper escape character
      CALL XMGRSUBS(T,L,'!',1,ESC,1)
C
      WRITE (IFIL,99001) T1,T2,T(1:L)
99001 FORMAT ('@  ',2A,' "',A,'"')
      END
C*==xmgrsubs.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE XMGRSUBS(T,L,S1,L1,S2,L2)
C **********************************************************************
C *                                                                    *
C *  substitute string  S1  in string  T  by string   S2               *
C *                                                                    *
C **********************************************************************
C
      IMPLICIT NONE
C*--XMGRSUBS424
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER LSTRMAX
      PARAMETER (LSTRMAX=80)
C
C Dummy arguments
C
      INTEGER L,L1,L2
      CHARACTER*(*) S1,S2,T
C
C Local variables
C
      INTEGER I1
      CHARACTER*(LSTRMAX) S
C
C*** End of declarations rewritten by SPAG
C
 100  CONTINUE
      I1 = INDEX(T(1:L),S1(1:L1))
      IF ( I1.EQ.0 ) RETURN
      S = T(1:(I1-1))//S2(1:L2)//T((I1+L1):L)
      T = S
      L = L - L1 + L2
      IF ( L.GT.LSTRMAX ) STOP '<XMGRSUBS>  L > LSTRMAX'
      GOTO 100
C
      END
C*==xmgrpoints.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE XMGRPOINTS(IFIL,NGRAPH,NC1,NC2,LWIDTH,KCOLOR,KSYMBOL)
C **********************************************************************
C *                                                                    *
C *  IFIL                  chanel for xmgr-file                        *
C *  NGRAPH                number of graphs   1  or  2                 *
C *  NC1                   number of curves for graph 1                *
C *  NC2                   number of curves for graph 2                *
C *  LWIDTH                line width                                  *
C *  KCOLOR                key to set color   0 (1) = off (on)         *
C *  KSYMBOL               key to set symbol  0 (1) = off (on)         *
C *                                                                    *
C **********************************************************************
C
      IMPLICIT NONE
C*--XMGRPOINTS484
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IFIL,KCOLOR,KSYMBOL,LWIDTH,NC1,NC2,NGRAPH
C
C Local variables
C
      INTEGER ICOLOR,IG,IS,ISYMBOL,NC
C
C*** End of declarations rewritten by SPAG
C
      DO IG = 0,(NGRAPH-1)
         WRITE (IFIL,99001) IG,IG
         IF ( IG.EQ.0 ) THEN
            NC = NC1
         ELSE
            NC = NC2
         END IF
         DO IS = 0,(NC-1)
C
C ----------------------------------------------------------------------
            IF ( LWIDTH.NE.0 ) THEN
               IF ( IS.LT.10 ) THEN
                  WRITE (IFIL,99002) IS,'symbol linewidth',LWIDTH
               ELSE
                  WRITE (IFIL,99003) IS,'symbol linewidth',LWIDTH
               END IF
            END IF
C
C ----------------------------------------------------------------------
            ICOLOR = 0
            IF ( KCOLOR.GT.0 ) THEN
               ICOLOR = IS + 1
            ELSE IF ( KCOLOR.EQ.-1 ) THEN
               IF ( IG.EQ.0 ) THEN
                  ICOLOR = (IS+1)/2 + 1
               ELSE
                  ICOLOR = IS/2 + 1
               END IF
C
            END IF
C
            IF ( ICOLOR.GT.0 ) THEN
               ICOLOR = MOD((ICOLOR-1),14) + 1
               IF ( IS.LT.10 ) THEN
                  WRITE (IFIL,99002) IS,'symbol color    ',ICOLOR
               ELSE
                  WRITE (IFIL,99003) IS,'symbol color    ',ICOLOR
               END IF
            END IF
C
C ----------------------------------------------------------------------
            ISYMBOL = 0
            IF ( KSYMBOL.GT.0 ) THEN
               ISYMBOL = IS + 1
            ELSE IF ( KSYMBOL.EQ.-1 ) THEN
               ISYMBOL = 2*MOD(IS,2) + 1
            END IF
C
            IF ( ISYMBOL.GT.0 ) THEN
               ISYMBOL = MOD((ISYMBOL-1),10) + 1
               IF ( IS.LT.10 ) THEN
                  WRITE (IFIL,99002) IS,'symbol',ISYMBOL
                  WRITE (IFIL,99002) IS,'symbol size 1.200000'
                  WRITE (IFIL,99002) IS,'symbol pattern 1'
                  WRITE (IFIL,99002) IS,'linestyle    0'
               ELSE
                  WRITE (IFIL,99003) IS,'symbol',ISYMBOL
                  WRITE (IFIL,99003) IS,'symbol size 1.200000'
                  WRITE (IFIL,99003) IS,'symbol pattern 1'
                  WRITE (IFIL,99003) IS,'linestyle    0'
               END IF
            END IF
C ----------------------------------------------------------------------
C
         END DO
C
      END DO
C
99001 FORMAT ('@g',i1,' on',/,'@with g',I1)
99002 FORMAT ('@  s',I1,'   ',A,'   ',I2)
99003 FORMAT ('@  s',I2,'   ',A,'   ',I2)
C
      END
C*==xmgrcurves.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE XMGRCURVES(IFIL,NGRAPH,NC1,NC2,LWIDTH,KCOLOR,KSTYLE)
C **********************************************************************
C *                                                                    *
C *  IFIL                  chanel for xmgr-file                        *
C *  NGRAPH                number of graphs   1  or  2                 *
C *  NC1                   number of curves for graph 1                *
C *  NC2                   number of curves for graph 2                *
C *  LWIDTH                line width                                  *
C *  KCOLOR                key to set color  0 (1) = off (on)          *
C *  KSTYLE                key to set style  0 (1) = off (on)          *
C *                                                                    *
C **********************************************************************
C
      IMPLICIT NONE
C*--XMGRCURVES597
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IFIL,KCOLOR,KSTYLE,LWIDTH,NC1,NC2,NGRAPH
C
C Local variables
C
      INTEGER ICOLOR,IG,IS,ISTYLE,NC
C
C*** End of declarations rewritten by SPAG
C
C --------------------------------------------- exclude white and yellow
C      INTEGER COLORTAB(14),STYLETAB(8)
C      DATA COLORTAB/1,2,4,3,6,7,8,9,10,11,12,13,14,15/
C      DATA STYLETAB/1,2,3,4,5,6,7,8/
C
      DO IG = 0,(NGRAPH-1)
         WRITE (IFIL,99001) IG,IG
         IF ( IG.EQ.0 ) THEN
            NC = NC1
         ELSE
            NC = NC2
         END IF
         DO IS = 0,(NC-1)
C
C ----------------------------------------------------------------------
            IF ( LWIDTH.NE.0 ) THEN
               IF ( IS.LT.10 ) THEN
                  WRITE (IFIL,99002) IS,'linewidth',LWIDTH
               ELSE
                  WRITE (IFIL,99003) IS,'linewidth',LWIDTH
               END IF
            END IF
C
C ----------------------------------------------------------------------
            ICOLOR = 0
            IF ( KCOLOR.GT.0 ) THEN
               ICOLOR = IS + 1
            ELSE IF ( KCOLOR.EQ.-1 ) THEN
               IF ( IG.EQ.0 ) THEN
                  ICOLOR = (IS+1)/2 + 1
               ELSE
                  ICOLOR = IS/2 + 1
               END IF
C
            END IF
C
            IF ( ICOLOR.GT.0 ) THEN
               ICOLOR = MOD((ICOLOR-1),14) + 1
               IF ( IS.LT.10 ) THEN
                  WRITE (IFIL,99002) IS,'color    ',ICOLOR
               ELSE
                  WRITE (IFIL,99003) IS,'color    ',ICOLOR
               END IF
            END IF
C
C ----------------------------------------------------------------------
            ISTYLE = 0
            IF ( KSTYLE.GT.0 ) THEN
               ISTYLE = IS + 1
            ELSE IF ( KSTYLE.EQ.-1 ) THEN
               ISTYLE = 2*MOD(IS,2) + 1
C               IF( IG.EQ.0 ) THEN
C                  ISTYLE = 2*MOD(IS,2) + 1
C               ELSE
C                  ISTYLE = 2*MOD(IS/2,2) + 1
C               END IF
            END IF
C
            IF ( ISTYLE.GT.0 ) THEN
               ISTYLE = MOD((ISTYLE-1),14) + 1
               IF ( IS.LT.10 ) THEN
                  WRITE (IFIL,99002) IS,'linestyle',ISTYLE
               ELSE
                  WRITE (IFIL,99003) IS,'linestyle',ISTYLE
               END IF
            END IF
C ----------------------------------------------------------------------
C
         END DO
C
      END DO
C
99001 FORMAT ('@g',i1,' on',/,'@with g',I1)
99002 FORMAT ('@  s',I1,'   ',A,'   ',I2)
99003 FORMAT ('@  s',I2,'   ',A,'   ',I2)
C
      END
C*==xmgrleg1.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE XMGRLEG1(IFIL,IGRAPH,NC,LEGTAB,X1,Y1)
C **********************************************************************
C *                                                                    *
C *  IFIL                  chanel for xmgr-file                        *
C *  IGRAPH                number of graph    0  or  1                 *
C *  NC                    number of curves for graph 1                *
C *  LEGTAB                NC  legends for the corresponding curves    *
C *  X1 Y1                 POSITION OF LEGEND 1 in WORLD COORDINATES   *
C *                                                                    *
C **********************************************************************
C
      IMPLICIT NONE
C*--XMGRLEG1712
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*3 CHARSIZE
      PARAMETER (CHARSIZE='1.2')
C
C Dummy arguments
C
      INTEGER IFIL,IGRAPH,NC
      REAL*8 X1,Y1
      CHARACTER*(*) LEGTAB(NC)
C
C Local variables
C
      INTEGER IS,LL,LTXT
      CHARACTER*10 STR10
      CHARACTER*80 TXT
C
C*** End of declarations rewritten by SPAG
C
      WRITE (IFIL,99001) IGRAPH,IGRAPH
C
      IF ( NC.GT.1 ) THEN
         WRITE (IFIL,99002) 'on'
         WRITE (IFIL,99002) 'loctype view'
         IF ( X1.LT.1D-6 ) THEN
            WRITE (IFIL,99003) 'x1 ',0.6D0
         ELSE
            WRITE (IFIL,99003) 'x1 ',X1
         END IF
         IF ( IGRAPH.EQ.1 ) THEN
            IF ( Y1.LT.1D-6 ) THEN
               WRITE (IFIL,99003) 'y1 ',0.83D0
            ELSE
               WRITE (IFIL,99003) 'y1 ',Y1
            END IF
         ELSE IF ( Y1.LT.1D-6 ) THEN
            WRITE (IFIL,99003) 'y1 ',0.48D0
         ELSE
            WRITE (IFIL,99003) 'y1 ',Y1
         END IF
         WRITE (IFIL,99002) 'on'
         WRITE (IFIL,99002) 'box color 0'
         WRITE (IFIL,99002) 'box pattern 0'
         WRITE (IFIL,99002) 'box linewidth 0'
         WRITE (IFIL,99002) 'box linestyle 0'
         WRITE (IFIL,99002) 'box fill color 0'
         WRITE (IFIL,99002) 'box fill pattern 0'
         WRITE (IFIL,99002) 'font 0'
         WRITE (IFIL,99002) 'char size '//CHARSIZE
         WRITE (IFIL,99002) 'length 6'
         WRITE (IFIL,99002) 'vgap 1'
         WRITE (IFIL,99002) 'hgap 1'
      END IF
C
      DO IS = 0,(NC-1)
C
         TXT = LEGTAB(IS+1)
         LTXT = LEN_TRIM(TXT)
C
         STR10 = 's'
         CALL STRING_ADD_N(STR10,IS)
         LL = LEN_TRIM(STR10)
C
         CALL XMGRWRITE(IFIL,STR10(1:LL),'  legend ',TXT,LTXT,'\\')
C
      END DO
C
99001 FORMAT ('@g',i1,' on',/,'@with g',I1)
99002 FORMAT ('@    legend ',A)
99003 FORMAT ('@    legend ',A,F6.2)
      END
C*==xmgrlegend.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE XMGRLEGEND(IFIL,NGRAPH,NC1,NC2,LEGTAB1,LEGTAB2)
C **********************************************************************
C *                                                                    *
C *  IFIL                  chanel for xmgr-file                        *
C *  NGRAPH                number of graphs   1  or  2                 *
C *  NC1                   number of curves for graph 1                *
C *  NC2                   number of curves for graph 2                *
C *  LEGTAB1               NC1 legends for the corresponding curves    *
C *  LEGTAB2               NC2 legends for the corresponding curves    *
C *                                                                    *
C **********************************************************************
C
      IMPLICIT NONE
C*--XMGRLEGEND815
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*3 CHARSIZE
      PARAMETER (CHARSIZE='1.2')
C
C Dummy arguments
C
      INTEGER IFIL,NC1,NC2,NGRAPH
      CHARACTER*(*) LEGTAB1(NC1),LEGTAB2(NC2)
C
C Local variables
C
      CHARACTER*6 CIS
      INTEGER IC0,IG,IS,L80,NC
      CHARACTER*80 STR80
C
C*** End of declarations rewritten by SPAG
C
      DO IG = 0,(NGRAPH-1)
C
         WRITE (IFIL,99001) IG,IG
C
         IF ( IG.EQ.0 ) THEN
            NC = NC1
         ELSE
            NC = NC2
         END IF
C
         IF ( NC.GT.1 ) THEN
            WRITE (IFIL,99002) 'on'
            WRITE (IFIL,99002) 'loctype view'
            WRITE (IFIL,99002) 'x1 0.70'
            IF ( IG.EQ.(NGRAPH-1) ) THEN
               WRITE (IFIL,99002) 'y1 0.83'
            ELSE
               WRITE (IFIL,99002) 'y1 0.48'
            END IF
            WRITE (IFIL,99002) 'on'
            WRITE (IFIL,99002) 'box color 0'
            WRITE (IFIL,99002) 'box pattern 0'
            WRITE (IFIL,99002) 'box linewidth 0'
            WRITE (IFIL,99002) 'box linestyle 0'
            WRITE (IFIL,99002) 'box fill color 0'
            WRITE (IFIL,99002) 'box fill pattern 0'
            WRITE (IFIL,99002) 'font 0'
            WRITE (IFIL,99002) 'char size '//CHARSIZE
            WRITE (IFIL,99002) 'length 6'
            WRITE (IFIL,99002) 'vgap 1'
            WRITE (IFIL,99002) 'hgap 1'
         END IF
C
         IC0 = ICHAR('0')
         DO IS = 0,(NC-1)
C
            IF ( IS.LT.10 ) THEN
               CIS = '   s'//CHAR(IC0+IS)//' '
            ELSE IF ( IS.EQ.0 ) THEN
               CIS = '   s'//CHAR(IC0+IS/10)//CHAR(IC0+IS-10*(IS/10))
            END IF
C
            IF ( IG.EQ.0 ) THEN
               STR80 = LEGTAB1(IS+1)
            ELSE
               STR80 = LEGTAB2(IS+1)
            END IF
            L80 = LEN_TRIM(STR80)
C
            CALL XMGRWRITE(IFIL,CIS,' legend ',STR80,L80,'\\')
C
         END DO
C
      END DO
C
99001 FORMAT ('@g',i1,' on',/,'@with g',I1)
99002 FORMAT ('@    legend ',A)
C
      END
C*==xmgrsubscripts.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE XMGRSUBSCRIPTS(TXT,LTXT,LTXTMAX)
C **********************************************************************
C *                                                                    *
C *    change subscripts in  TXT  according to  XMGR syntax            *
C *    works only for concentrations like "A_23B_77" or "A_ .23B_ .77" *
C *    '!' is used as the xmgr escape parameter -                      *
C *    it will be replaced later by   <XMGRWRITE>                      *
C *                                                                    *
C **********************************************************************
C
      IMPLICIT NONE
C*--XMGRSUBSCRIPTS922
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER LTMAX
      PARAMETER (LTMAX=120)
C
C Dummy arguments
C
      INTEGER LTXT,LTXTMAX
      CHARACTER*(*) TXT
C
C Local variables
C
      CHARACTER C
      INTEGER I,I0,I9,IC,ID,K,L
      CHARACTER*(LTMAX) T
C
C*** End of declarations rewritten by SPAG
C
      IF ( LTXTMAX.GT.LTMAX ) STOP '1 in <XMGRSUBSCRIPTS>'
C
      I0 = ICHAR('0')
      I9 = ICHAR('9')
      ID = ICHAR('.')
C
      K = 0
      L = 0
      DO I = 1,LTXT
         C = TXT(I:I)
         IC = ICHAR(C)
         IF ( C.EQ.'_' ) THEN
            K = 1
            T((L+1):(L+2)) = '!s'
            L = L + 2
         ELSE IF ( K.EQ.1 ) THEN
            IF ( C.NE.' ' ) THEN
               IF ( (IC.NE.ID) .AND. ((IC.LT.I0) .OR. (IC.GE.I9)) ) THEN
                  T((L+1):(L+2)) = '!N'
                  L = L + 2
                  K = 0
               END IF
               L = L + 1
               T(L:L) = C
            END IF
         ELSE
            L = L + 1
            T(L:L) = C
         END IF
C
      END DO
C
      IF ( L.GT.LTXTMAX ) STOP '2 in <XMGRSUBSCRIPTS>'
      LTXT = L
      TXT(1:L) = T(1:L)
C
      END
C*==xmgrticks.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE XMGRTICKS(XA0,XA,KA,XB0,XB,KB,DX0,DX,KD)
C   ********************************************************************
C   *                                                                  *
C   *       find a reasonable setting for the plot ticks               *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--XMGRTICKS1004
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 DX,DX0,XA,XA0,XB,XB0
      INTEGER KA,KB,KD
C
C Local variables
C
      REAL*8 DA,DB,DC,DL
C
C*** End of declarations rewritten by SPAG
C
      IF ( KD.EQ.0 ) THEN
         DX = DX0
      ELSE
         DL = LOG10(XB0-XA0)
         DA = INT(DL)
         IF ( DL.GT.0D0 ) DA = DA + 1D0
         DB = 10D0**(DL-DA)
         DC = 10D0**DA
C
         IF ( DB.LT.0.16D0 ) THEN
            DX = 0.02D0*DC
         ELSE IF ( DB.LT.0.3D0 ) THEN
            DX = 0.04D0*DC
         ELSE IF ( DB.LT.0.5D0 ) THEN
            DX = 0.05D0*DC
         ELSE
            DX = 0.1D0*DC
         END IF
      END IF
C
      IF ( KA.EQ.0 ) THEN
         XA = XA0
      ELSE IF ( XA0.GE.0D0 ) THEN
         XA = DBLE(INT(XA0/DX))*DX
      ELSE
         XA = DBLE(INT(XA0/DX)-1)*DX
      END IF
C
      IF ( KB.EQ.0 ) THEN
         XB = XB0
      ELSE
         XB = DBLE(INT(XB0/DX)+1)*DX
      END IF
C
      END
C*==xmgr4head.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE XMGR4HEAD(DATSET,LDATSET,IDENT1,LIDENT1,IDENT2,LIDENT2,
     &                     FILNAM,LFNMAX,LFN,IFIL,NGRAPH,XMIN,KXMIN,
     &                     XMAX,KXMAX,YMIN1,KYMIN1,YMAX1,KYMAX1,YMIN2,
     &                     KYMIN2,YMAX2,KYMAX2,XTXT,LXTXT,YTXT1,LYTXT1,
     &                     YTXT2,LYTXT2,TITLE,LTITLE,SUBTITLE,LSUBTITLE,
     &                     KNOHEADER)
C **********************************************************************
C *                                                                    *
C * dummy REAL*4 version to call the proper REAL*8 version             *
C *                                                                    *
C **********************************************************************
C *                                                                    *
C *  DATSET    LDATSET     dataset name and length of string           *
C *  IDENT1    LIDENT1     add. identifier and length of string > 0 !! *
C *  IDENT2    LIDENT2     add. identifier and length of string        *
C *  FILNAM    LFNMAX LFN  file name to be created, its max. and       *
C *                        actual length on exit                       *
C *  IFIL                  chanel for xmgr-file                        *
C *  NGRAPH                number of graphs   1  or  2                 *
C *  XMIN      KXMIN       x-start and key to fix (0) or float (1)     *
C *  XMAX      KXMAX       x-end   and key to fix (0) or float (1)     *
C *  YMIN1     KYMIN1      y-start and key to fix (0) or float (1)     *
C *  YMAX1     KYMAX1      y-end   and key ...  for graph 1            *
C *  YMIN2     KYMIN2      y-start and key ...  for graph 2            *
C *  YMAX2     KYMAX2      y-end   and key ...  used if NGRAPH=2       *
C *  XTXT      LXTXT       text for x-axis label                       *
C *  YTXT1     LYTXT1      text for y-axis label for graph 1           *
C *  YTXT2     LYTXT2      text for y-axis label for graph 2           *
C *  TITLE     LTITLE      title and length of string                  *
C *  SUBTITLE  LSUBTITLE   subtitle and length of string               *
C *  KNOHEADER             key to get (false) or suppress (true) title *
C *                                                                    *
C *  NOTE:  only  FILNAM and LFN will be changed on exit               *
C *         apart of IDENT1 any text may have length 0                 *
C *                                                                    *
C *         use '!' instead of '\' for the xmgr escape character       *
C *         in the text variables --- see <XMGRWRITE>, where also      *
C *         the macro list is stored  !UP, !DN,                        *
C *                                                                    *
C **********************************************************************
C
      IMPLICIT NONE
C*--XMGR4HEAD1108
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*(*) DATSET,FILNAM,IDENT1,IDENT2,SUBTITLE,TITLE,XTXT,
     &              YTXT1,YTXT2
      INTEGER IFIL,KXMAX,KXMIN,KYMAX1,KYMAX2,KYMIN1,KYMIN2,LDATSET,LFN,
     &        LFNMAX,LIDENT1,LIDENT2,LSUBTITLE,LTITLE,LXTXT,LYTXT1,
     &        LYTXT2,NGRAPH
      LOGICAL KNOHEADER
      REAL*4 XMAX,XMIN,YMAX1,YMAX2,YMIN1,YMIN2
C
C*** End of declarations rewritten by SPAG
C
      CALL XMGRHEAD(DATSET,LDATSET,IDENT1,LIDENT1,IDENT2,LIDENT2,FILNAM,
     &              LFNMAX,LFN,IFIL,NGRAPH,DBLE(XMIN),KXMIN,DBLE(XMAX),
     &              KXMAX,DBLE(YMIN1),KYMIN1,DBLE(YMAX1),KYMAX1,
     &              DBLE(YMIN2),KYMIN2,DBLE(YMAX2),KYMAX2,XTXT,LXTXT,
     &              YTXT1,LYTXT1,YTXT2,LYTXT2,TITLE,LTITLE,SUBTITLE,
     &              LSUBTITLE,KNOHEADER)
C
      END
C*==xmgr4table.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE XMGR4TABLE(IG,IS,X,Y,W,N,IFIL)
C **********************************************************************
C *                                                                    *
C * dummy REAL*4 version to call the proper REAL*8 version             *
C *                                                                    *
C **********************************************************************
C *                                                                    *
C *   write out the data for graph  IG  and set  IS                    *
C *   W is used as weight parameter for  Y  normally set to  1         *
C *                                                                    *
C **********************************************************************
C
      IMPLICIT NONE
C*--XMGR4TABLE1154
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IFIL,IG,IS,N
      REAL*4 W
      REAL*4 X(N),Y(N)
C
C Local variables
C
      INTEGER I
      REAL*8 XR8(N),YR8(N)
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,N
         XR8(I) = DBLE(X(I))
         YR8(I) = DBLE(Y(I))
      END DO
C
      CALL XMGRTABLE(IG,IS,XR8,YR8,DBLE(W),N,IFIL)
C
      END
C*==xmgrfnam.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE XMGRFNAM(DATSET,LDATSET,IDENT1,LIDENT1,IDENT2,LIDENT2,
     &                    FILEXT,LFILEXT,FILNAM,LFN,LFNMAX)
C **********************************************************************
C *                                                                    *
C *  COMPOSE FILENAME FROM DATASET, IDENT1, IDENT2 and FILEXT          *
C *                                                                    *
C *  DATSET    LDATSET     dataset name and length of string           *
C *  IDENT1    LIDENT1     add. identifier and length of string > 0 !! *
C *  IDENT2    LIDENT2     add. identifier and length of string        *
C *  FILEXT    LFILEXT        add extension and length of string       *
C *  FILNAM    LFNMAX LFN  file name to be created, its max. and       *
C *                        actual length on exit                       *
C *                                                                    *
C **********************************************************************
C
      IMPLICIT NONE
C*--XMGRFNAM1207
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER LSTRMAX
      PARAMETER (LSTRMAX=80)
C
C Dummy arguments
C
      CHARACTER*(*) DATSET,FILEXT,FILNAM,IDENT1,IDENT2
      INTEGER LDATSET,LFILEXT,LFN,LFNMAX,LIDENT1,LIDENT2
C
C Local variables
C
      INTEGER LDATSETEFF
      CHARACTER*(LSTRMAX) STR
C
C*** End of declarations rewritten by SPAG
C
      IF ( LFNMAX.GT.LSTRMAX ) STOP '<XMGRFNAM>   LFNMAX > LSTRMAX'
C
      IF ( LIDENT1.LE.0 ) STOP '<XMGRFNAM>   supply IDENT1'
C
      IF ( .FALSE. ) THEN
         WRITE (6,*) 'DATSET=',DATSET(1:LDATSET)
         WRITE (6,*) 'LDATSET=',LDATSET
         WRITE (6,*) 'IDENT1=',IDENT1(1:LIDENT1)
         WRITE (6,*) 'LIDENT1=',LIDENT1
         WRITE (6,*) 'IDENT2=',IDENT2(1:LIDENT2)
         WRITE (6,*) 'LIDENT2=',LIDENT2
         WRITE (6,*) 'FILEXT=',FILEXT(1:LFILEXT)
         WRITE (6,*) 'LFILEXT=',LFILEXT
         WRITE (6,*) 'FILNAM=',FILNAM(1:LFN)
         WRITE (6,*) 'LFN=',LFN
         WRITE (6,*) 'LFNMAX=',LFNMAX
      END IF
C
      LDATSETEFF = LDATSET
      IF ( LDATSETEFF.LE.0 ) THEN
         FILNAM = IDENT1(1:LIDENT1)
         LFN = LIDENT1
      ELSE
         IF ( DATSET(LDATSET:LDATSET).EQ.'_' ) LDATSETEFF = LDATSET - 1
         IF ( (LDATSETEFF-LIDENT1).GT.0 ) THEN
            STR = DATSET((LDATSETEFF-LIDENT1):LDATSETEFF)
            IF ( STR(1:(LIDENT1+1)).EQ.'.'//IDENT1(1:LIDENT1) .OR. 
     &           STR(1:(LIDENT1+1)).EQ.'_'//IDENT1(1:LIDENT1) ) THEN
               FILNAM = DATSET(1:LDATSETEFF)
               LFN = LDATSETEFF
            ELSE
               FILNAM = DATSET(1:LDATSETEFF)//'_'//IDENT1(1:LIDENT1)
               LFN = LDATSETEFF + 1 + LIDENT1
            END IF
         ELSE
            FILNAM = DATSET(1:LDATSETEFF)//'_'//IDENT1(1:LIDENT1)
            LFN = LDATSETEFF + 1 + LIDENT1
         END IF
      END IF
C
      IF ( LIDENT2.GT.0 ) THEN
         FILNAM = FILNAM(1:LFN)//'_'//IDENT2(1:LIDENT2)
         LFN = LFN + 1 + LIDENT2
      END IF
C
      FILNAM = FILNAM(1:LFN)//FILEXT(1:LFILEXT)
      LFN = LFN + LFILEXT
C
      IF ( LFN.GT.LFNMAX ) THEN
         WRITE (6,*) 'LFN=',LFN,' LFNMAX=',LFNMAX
         STOP '<XMGRFNAM>    LFN > LFNMAX'
      END IF
C
      END
C*==xmgrtick.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE XMGRTICK(IFIL,IG,AXIS,SEL,DEL)
C **********************************************************************
C *                                                                    *
C *  IFIL                  chanel for xmgr-file                        *
C *  IG                    number of graph:   0  or  1                 *
C *  AXIS                  x or y                                      *
C *  SEL                   tick type:  major or minor                  *
C *  DEL                   tick step                                   *
C *                                                                    *
C **********************************************************************
C
      IMPLICIT NONE
C*--XMGRTICK1309
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*1 AXIS
      REAL*8 DEL
      INTEGER IFIL,IG
      CHARACTER*5 SEL
C
C*** End of declarations rewritten by SPAG
C
      WRITE (IFIL,99001) IG,AXIS,SEL,DEL
99001 FORMAT ('@with g',I1,/,'@    ',A,'axis  tick ',A,' ',F7.4)
      END
C*==dnrm2.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      REAL*8 FUNCTION DNRM2(N,X,INCX)
      IMPLICIT NONE
C*--DNRM21336
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
C
C Dummy arguments
C
      INTEGER INCX,N
      REAL*8 X(*)
C
C Local variables
C
      REAL*8 ABSXI,NORM,SCALE,SSQ
      INTEGER IX
C
C*** End of declarations rewritten by SPAG
C
C     .. Scalar Arguments ..
C     ..
C     .. Array Arguments ..
C     ..
C
C  Purpose
C  =======
C
C  DNRM2 returns the euclidean norm of a vector via the function
C  name, so that
C
C     DNRM2 := sqrt( x'*x )
C
C
C  -- This version written on 25-October-1982.
C     Modified on 14-October-1993 to inline the call to DLASSQ.
C     Sven Hammarling, Nag Ltd.
C
C
C     .. Parameters ..
C     ..
C     .. Local Scalars ..
C     ..
C     .. Intrinsic Functions ..
C     ..
      IF ( N.LT.1 .OR. INCX.LT.1 ) THEN
         NORM = ZERO
      ELSE IF ( N.EQ.1 ) THEN
         NORM = ABS(X(1))
      ELSE
         SCALE = ZERO
         SSQ = ONE
C        The following loop is equivalent to this call to the LAPACK
C        auxiliary routine:
C        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
C
         DO IX = 1,1 + (N-1)*INCX,INCX
            IF ( X(IX).NE.ZERO ) THEN
               ABSXI = ABS(X(IX))
               IF ( SCALE.LT.ABSXI ) THEN
                  SSQ = ONE + SSQ*(SCALE/ABSXI)**2
                  SCALE = ABSXI
               ELSE
                  SSQ = SSQ + (ABSXI/SCALE)**2
               END IF
            END IF
         END DO
         NORM = SCALE*SQRT(SSQ)
      END IF
C
      DNRM2 = NORM
C
C     End of DNRM2.
C
      END
C*==rvecspro.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE RVECSPRO(A,B,C,V)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the vector-vector operation         V = (A x B) * C   *
C   *                                                                  *
C   *   A,B,C   real vectors of dimension 3                            *
C   *   V       spatial product of A, B and C with sign kept           *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--RVECSPRO1437
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 V
      REAL*8 A(3),B(3),C(3)
C
C Local variables
C
      REAL*8 D(3)
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      CALL RVECXPRO(A,B,D)
C
      V = 0D0
      DO I = 1,3
         V = V + D(I)*C(I)
      END DO
C
      END
C*==rvecxpro.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE RVECXPRO(A,B,C)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the vector-vector operation           C = A x B       *
C   *                                                                  *
C   *   A,B,C   real vectors of dimension 3                            *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--RVECXPRO1482
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 A(3),B(3),C(3)
C
C*** End of declarations rewritten by SPAG
C
      C(1) = A(2)*B(3) - A(3)*B(2)
      C(2) = A(3)*B(1) - A(1)*B(3)
      C(3) = A(1)*B(2) - A(2)*B(1)
C
      END
C*==rvecspat.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE RVECSPAT(A,B,C,V,KEY)
C   ********************************************************************
C   *                                                                  *
C   *   spatial product                          V = A * ( B x C )     *
C   *                                                                  *
C   *   A,B,C   real vectors of dimension 3                            *
C   *                                                                  *
C   *   for  KEY = 1:                            V ->  |V|             *
C   *                                                                  *
C   *                                                                  *
C   *   NOTE:  the input is copied to local arrays to avoid            *
C   *          warnings when checking with ftnchek                     *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--RVECSPAT1522
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER KEY
      REAL*8 V
      REAL*8 A(3),B(3),C(3)
C
C Local variables
C
      REAL*8 AA(3),BB(3),CC(3),DD(3)
      REAL*8 DDOT
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
C
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,3
         AA(I) = A(I)
         BB(I) = B(I)
         CC(I) = C(I)
      END DO
C
      CALL RVECXPRO(BB,CC,DD)
C
      V = DDOT(3,AA,1,DD,1)
C
      IF ( KEY.EQ.1 ) V = ABS(V)
C
      END
C*==string_add_n.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE STRING_ADD_N(STRING,N)
C   ********************************************************************
C   *                                                                  *
C   *  add the integer N to     STRING                                 *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--STRING_ADD_N1577
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N
      CHARACTER*(*) STRING
C
C Local variables
C
      INTEGER I,I0,J,LL,M,NDIG
      REAL*8 X
C
C*** End of declarations rewritten by SPAG
C
      LL = LEN_TRIM(STRING)
C
      M = N
      DO I = 1,10000
         IF ( M.LT.10 ) THEN
            NDIG = I
            EXIT
         END IF
         M = M/10
      END DO
      X = DBLE(N)
      I0 = ICHAR('0')
      DO I = NDIG,1, - 1
         LL = LL + 1
         J = INT(X/10.0D0**(I-1)) - INT(X/10.0D0**I)*10
         STRING(LL:LL) = CHAR(I0+J)
      END DO
      END
C*==ddot.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      REAL*8 FUNCTION DDOT(N,DX,INCX,DY,INCY)
      IMPLICIT NONE
C*--DDOT1625
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER INCX,INCY,N
      REAL*8 DX(*),DY(*)
C
C Local variables
C
      REAL*8 DTEMP
      INTEGER I,IX,IY,M,MP1
C
C*** End of declarations rewritten by SPAG
C
C  Purpose
C  =======
C
C     forms the dot product of two vectors.
C     uses unrolled loops for increments equal to one.
C     jack dongarra, linpack, 3/11/78.
C     modified 12/3/93, array(1) declarations changed to array(*)
C
      DDOT = 0.0D0
      DTEMP = 0.0D0
      IF ( N.LE.0 ) RETURN
      IF ( INCX.EQ.1 .AND. INCY.EQ.1 ) THEN
C
C        code for both increments equal to 1
C
C
C        clean-up loop
C
         M = MOD(N,5)
         IF ( M.NE.0 ) THEN
            DO I = 1,M
               DTEMP = DTEMP + DX(I)*DY(I)
            END DO
            IF ( N.LT.5 ) THEN
               DDOT = DTEMP
               GOTO 99999
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,5
            DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) + DX(I+2)
     &              *DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
         END DO
         DDOT = DTEMP
      ELSE
C
C        code for unequal increments or equal increments
C          not equal to 1
C
         IX = 1
         IY = 1
         IF ( INCX.LT.0 ) IX = (-N+1)*INCX + 1
         IF ( INCY.LT.0 ) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DTEMP = DTEMP + DX(IX)*DY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
         DDOT = DTEMP
         RETURN
      END IF
99999 CONTINUE
      END
