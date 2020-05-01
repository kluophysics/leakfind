C*==symplotuc.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE SYMPLOTUC(RASRAD,RASSCL)
C   ********************************************************************
C   *                                                                  *
C   *        create input for the  rasmol  program                     *
C   *        to plot the crystallographic unit cell                    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:RWS
      USE MOD_FILES,ONLY:DATSET,LDATSET
      USE MOD_TYPES,ONLY:NT,TXT_T,LTXT_T,IMT,NAT
      USE MOD_SITES,ONLY:QBAS,IQAT
      USE MOD_LATTICE,ONLY:ALAT,ABAS
      USE MOD_CONSTANTS,ONLY:A0_ANG
      IMPLICIT NONE
C*--SYMPLOTUC17
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 RASRAD,RASSCL
C
C Local variables
C
      INTEGER I,I0,IA,IM,IQ,IRAD,IRADLIM,IT,IW,LL
      REAL*8 RS,S,SCLIRAD
      CHARACTER*40 STR40
C
C*** End of declarations rewritten by SPAG
C
      IW = 80
C
      S = RASSCL*A0_ANG*ALAT
      RS = RASRAD*RASSCL*A0_ANG*250
      I0 = 16
C
      OPEN (IW,FILE='rasmol_uc.pdb')
C
      WRITE (IW,FMT=99005) 'of '//DATSET(1:(LDATSET-1))
C
C------------------------------------------ specify corners of unit cell
C
      WRITE (IW,FMT=99002) 1,1,0D0,0D0,0D0,0D0
      DO I = 1,3
         WRITE (IW,FMT=99002) (I+1),(I+1),S*ABAS(1,I),S*ABAS(2,I),
     &                        S*ABAS(3,I),0D0
      END DO
C
      WRITE (IW,FMT=99002) 5,5,S*(ABAS(1,1)+ABAS(1,2)),
     &                     S*(ABAS(2,1)+ABAS(2,2)),
     &                     S*(ABAS(3,1)+ABAS(3,2)),0D0
      WRITE (IW,FMT=99002) 6,6,S*(ABAS(1,1)+ABAS(1,3)),
     &                     S*(ABAS(2,1)+ABAS(2,3)),
     &                     S*(ABAS(3,1)+ABAS(3,3)),0D0
      WRITE (IW,FMT=99002) 7,7,S*(ABAS(1,3)+ABAS(1,2)),
     &                     S*(ABAS(2,3)+ABAS(2,2)),
     &                     S*(ABAS(3,3)+ABAS(3,2)),0D0
      WRITE (IW,FMT=99002) 8,8,S*(ABAS(1,1)+ABAS(1,2)+ABAS(1,3)),
     &                     S*(ABAS(2,1)+ABAS(2,2)+ABAS(2,3)),
     &                     S*(ABAS(3,1)+ABAS(3,2)+ABAS(3,3)),0D0
C
C------------------------- specify corners of cube with edge length ALAT
C
      WRITE (IW,FMT=99002) 9,9,S*0D0,S*0D0,S*0D0,3D0
      WRITE (IW,FMT=99002) 10,10,S*1D0,S*0D0,S*0D0,3D0
      WRITE (IW,FMT=99002) 11,11,S*0D0,S*1D0,S*0D0,3D0
      WRITE (IW,FMT=99002) 12,12,S*0D0,S*0D0,S*1D0,3D0
      WRITE (IW,FMT=99002) 13,13,S*1D0,S*1D0,S*0D0,3D0
      WRITE (IW,FMT=99002) 14,14,S*1D0,S*0D0,S*1D0,3D0
      WRITE (IW,FMT=99002) 15,15,S*0D0,S*1D0,S*1D0,3D0
      WRITE (IW,FMT=99002) 16,16,S*1D0,S*1D0,S*1D0,3D0
C
C---------------------------------------------- specify atomic positions
C
      I = I0
      IRAD = 0
      DO IT = 1,NT
         IM = IMT(IT)
         IRAD = MAX(IRAD,NINT(RWS(IM)*RS))
C
         DO IA = 1,NAT(IT)
            I = I + 1
            IQ = IQAT(IA,IT)
C
            WRITE (IW,FMT=99001) I,I,S*QBAS(1,IQ),S*QBAS(2,IQ),
     &                           S*QBAS(3,IQ),DBLE(IT)*3D0/DBLE(NT)
         END DO
      END DO
C
      WRITE (IW,FMT=99007)
      CLOSE (IW)
C
C--------------------------------------------------- write RASMOL script
C
      OPEN (IW,FILE='rasmol_uc.ras')
      WRITE (IW,*) 'load ''rasmol_uc.pdb'' '
      WRITE (IW,*) 'color temperature'
      WRITE (IW,*) 'label true'
      WRITE (IW,*) 'set fontsize 20'
      WRITE (IW,*) 'select 1-16'
      WRITE (IW,*) 'label '' '' '
C
C------------------------ rasmol expects the radii in multiples of 1/250
C--------------------------------------- only values <= 500 are accepted
      IRADLIM = 500
      IF ( IRAD.LE.IRADLIM ) THEN
         SCLIRAD = 1D0
      ELSE
         SCLIRAD = DBLE(IRADLIM)/DBLE(IRAD)
         WRITE (6,99008) SCLIRAD
      END IF
C
      I = I0
      DO IT = 1,NT
         IM = IMT(IT)
         IRAD = MIN(NINT(RWS(IM)*RS*SCLIRAD),IRADLIM)
C
         DO IA = 1,NAT(IT)
            I = I + 1
            IQ = IQAT(IA,IT)
            WRITE (IW,*) 'select ',I
            STR40 = TXT_T(IT)(1:LTXT_T(IT))//' '
            CALL STRING_ADD_N(STR40,IQ)
            LL = LEN_TRIM(STR40)
            STR40 = STR40(1:LL)//':'
            CALL STRING_ADD_N(STR40,IT)
            LL = LEN_TRIM(STR40)
            WRITE (IW,99003) STR40(1:LL)
C
            WRITE (IW,99004) IRAD
         END DO
      END DO
      WRITE (IW,*) 'select all '
      WRITE (IW,*) 'set axes on '
C
      CLOSE (IW)
C
      WRITE (6,99006)
C
      RETURN
C
99001 FORMAT ('ATOM  ',I5,'          ',I5,'    ',3F8.3,'  0.00',F8.3)
99002 FORMAT ('HETATM',I5,'          ',I5,'    ',3F8.3,'  0.00',F8.3)
99003 FORMAT (' label ''',a,''' ')
99004 FORMAT (' cpk ',I5)
99005 FORMAT ('HEADER    crystallographic UNIT CELL ',A,/,
     &        'SOURCE    SPRKKR - program       ',/,
     &        'AUTHOR    H. Ebert               ',/,
     &        'REMARK    None                   ')
99006 FORMAT (/,5X,'unit cell data stored in rasmol data-file ',
     &        ' rasmol_uc.pdb',/,5X,
     &        'view via:   rasmol  -script rasmol_uc.ras')
99007 FORMAT ('CONECT    1    2',/,'CONECT    1    3',/,
     &        'CONECT    1    4',/,'CONECT    5    2',/,
     &        'CONECT    5    3',/,'CONECT    5    8',/,
     &        'CONECT    6    2',/,'CONECT    6    4',/,
     &        'CONECT    6    8',/,'CONECT    7    3',/,
     &        'CONECT    7    4',/,'CONECT    7    8',/,
     &        'CONECT    9   10',/,'CONECT    9   11',/,
     &        'CONECT    9   12',/,'CONECT   13   10',/,
     &        'CONECT   13   11',/,'CONECT   13   16',/,
     &        'CONECT   14   10',/,'CONECT   14   12',/,
     &        'CONECT   14   16',/,'CONECT   15   11',/,
     &        'CONECT   15   12',/,'CONECT   15   16',/,'END')
99008 FORMAT (/,70('*'),/,5X,'warning from <SYMPLOTUC>',/,/,5X,
     &        'RASMOL max. atomic radius exceeded ',/,/,5X,
     &        'atomic radii reduced by ',F10.5,/,70('*'),/)
      END
