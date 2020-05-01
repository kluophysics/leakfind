C*==clxps.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CLXPS(ICALL,CALCINT,GETIRRSOL,IPRINT,TAUT,TSST,MSST,
     &                 MEZZ,MEZJ,GCOR,FCOR,ECOR,SZCOR,KAPCOR,MM05COR,
     &                 NKPCOR,IKMCOR,IZERO,ITXRAY,NCSTMAX,NPOLMAX,
     &                 LCHPMAX)
C   ********************************************************************
C   *                                                                  *
C   *       calculation of core level XPS - spectra                    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:R,NRMAX,JRWS,R2DRDI
      USE MOD_ANGMOM,ONLY:A1_ADA,A2_ADA,NMEMAX,NKMMAX,CGC,NKMQ,NKM,
     &    NCPLWF,WKM1,WKM2
      USE MOD_SITES,ONLY:NQ,IQAT
      USE MOD_CALCMODE,ONLY:ORBPOL
      USE MOD_FILES,ONLY:DATSET,LSYSTEM,SYSTEM,LDATSET,WRTAU,WRTAUMQ,
     &    LRECREAL8,IFILCBWF,IFILCORWF,IFILBUILDBOT,WRBUILDBOT,
     &    FOUND_REAL_ARRAY,STRINP
      USE MOD_ENERGY,ONLY:NETAB,EFERMI,ETAB
      USE MOD_TYPES,ONLY:NTMAX,LCXRAY,NCXRAY,BT,VT,IMT,NT,TXT_T,CTL,
     &    IKMCPLWF,NCPLWFMAX
      USE MOD_CONSTANTS,ONLY:C0,CI,RY_EV,PI
      IMPLICIT NONE
C*--CLXPS25
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CLXPS')
      INTEGER NGEOMAX
      PARAMETER (NGEOMAX=10)
C
C Dummy arguments
C
      LOGICAL CALCINT,GETIRRSOL
      INTEGER ICALL,IPRINT,ITXRAY,LCHPMAX,NCSTMAX,NPOLMAX
      REAL*8 ECOR(NCSTMAX),FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX),
     &       SZCOR(NCSTMAX)
      INTEGER IKMCOR(NCSTMAX,2),IZERO(NCSTMAX),KAPCOR(NCSTMAX),
     &        MM05COR(NCSTMAX),NKPCOR(NCSTMAX)
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),TAUT(NKMMAX,NKMMAX,NTMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      REAL*8 ANGE(2),ANGK(2),ANGQ(2),BCOR(:),BCORS(:),C,DBDR(:),DIR(3),
     &       DVDR(:),EIMAG,EPHOT,EPHOTEV,EVX,EVY,EVZ,FUN(:,:,:,:),JJ,
     &       K0X,K0Y,K0Z,MJ,MS,PHIEPH(NGEOMAX),PHIKEL(NGEOMAX),
     &       PHIQPH(NGEOMAX),RINT(:),TETEPH(NGEOMAX),TETKEL(NGEOMAX),
     &       TETQPH(NGEOMAX),W(:),WZ2,XNORM(2)
      LOGICAL ANGRES,FOUND,MECHECK,STRPOL,USETSS
      CHARACTER*1 CHPOL1(5),SHELL(5)
      CHARACTER*2 CL
      CHARACTER*10 CP,STRE,STRK
      COMPLEX*16 EPS,ERYD,ICY,ME(:),MEA(:,:,:),MEX,MEY,MEZ,P,RSUM(:),
     &           SSST(NKMMAX,NKMMAX,NTMAX),WGTR(:),WMAT(:,:,:),YLM(:),
     &           YLMK,ZF(:,:,:),ZG(:,:,:)
      CHARACTER*80 FILNAM,SPEC
      INTEGER IA_ERR,ICST,IFIL,IGEO,IKM,IKMC,IM,IPOL,IPOLAR,IQ,IR,IRTOP,
     &        IT,ITIMREV,JCST,K,L,LAM,LAMP,LFN,LL,LM,LSP,LSUBSH(0:4),
     &        LTXT_T,M,MMS,MSM05,N,NCST,NGEO,NIN,NLAM,NLMCHPMAX,NPOL,
     &        NTXRSGRP
      CHARACTER*3 MEFORM,SUBSH(0:4),SUBSHP(0:4)
      CHARACTER*40 STRQ
C
C*** End of declarations rewritten by SPAG
C
      DATA SHELL/'K','L','M','N','O'/
      DATA LSUBSH/1,3,3,3,3/
      DATA SUBSH/'1  ','2,3','4,5','6,7','8,9'/
      DATA SUBSHP/'1  ','23 ','45 ','67 ','89 '/
      DATA CHPOL1/'+','-','z','x','y'/
      DATA ERYD/(999999D0,999999D0)/
C
      ALLOCATABLE DBDR,BCOR,DVDR,RINT,WMAT
      ALLOCATABLE W,BCORS,ME,MEA
      ALLOCATABLE FUN,RSUM,YLM,ZG,ZF,WGTR
C
      ALLOCATE (WGTR(NRMAX))
      ALLOCATE (ZG(NRMAX,NCPLWFMAX,NKMMAX),ZF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (FUN(NCSTMAX,3,NPOLMAX,NGEOMAX),BCORS(NTMAX))
      ALLOCATE (DBDR(NRMAX),BCOR(NTMAX),DVDR(NRMAX))
      ALLOCATE (RINT(NRMAX),WMAT(NKMMAX,NKMMAX,NPOLMAX))
      ALLOCATE (W(NPOLMAX),RSUM(NPOLMAX),ME(NPOLMAX))
      ALLOCATE (MEA(NCSTMAX,NKMMAX,NPOLMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MEBIRR')
C
      NLMCHPMAX = (LCHPMAX+1)**2
      ALLOCATE (YLM(NLMCHPMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: YLM')
C
C ======================================================================
C               select MODE to treat time reversal operation
C ======================================================================
      ITIMREV = 1
C
      OPEN (UNIT=IFILCBWF,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=(5+2*(2+2*NRMAX*2))*LRECREAL8)
      OPEN (UNIT=IFILCORWF,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=(5+2*(2+2*NRMAX*2))*LRECREAL8)
C
C ======================================================================
      IF ( WRTAU .OR. ICALL.EQ.2 ) THEN
C ======================================================================
C
         WRITE (6,99001)
         IF ( NPOLMAX.NE.3 ) CALL STOP_MESSAGE(ROUTINE,'NPOLMAX <> 3')
C
         NPOL = 3
C
         CALL INPUT_FIND_SECTION('TASK',0)
C
         CALL SECTION_SET_STRING('ME',MEFORM,'ADA',0)
C
         CALL SECTION_FIND_KEYWORD('MECHECK',MECHECK)
         CALL SECTION_FIND_KEYWORD('LINPOL',STRPOL)
         CALL SECTION_FIND_KEYWORD('USETSS',USETSS)
         CALL SECTION_FIND_KEYWORD('ANGRES',ANGRES)
C
         CALL SECTION_SET_REAL('EPHOT',EPHOTEV,1253.6D0,0)
         EPHOT = EPHOTEV/RY_EV
C
         CALL SECTION_SET_REAL('IME',EIMAG,0.0D0,0)
         IF ( ABS(EIMAG).LT.1D-8 ) EIMAG = 0.01D0
C
         CALL SECTION_SET_INTEGER('IT',ITXRAY,1,0)
         IT = ITXRAY
C
         CALL SECTION_GET_CORE_LEVEL_INFO(CL,NCXRAY(IT),LCXRAY(IT))
C
         IQ = IQAT(1,IT)
         IM = IMT(IT)
         NLAM = NKMQ(IQ)
         C = CTL(IT,1)
         IF ( TXT_T(IT)(2:2).EQ.' ' ) THEN
            LTXT_T = 1
         ELSE
            LTXT_T = 2
         END IF
         IF ( LDATSET.NE.0 ) THEN
            FILNAM = DATSET(1:LDATSET)//TXT_T(IT)(1:LTXT_T)//'.'
         ELSE
            FILNAM = TXT_T(IT)(1:LTXT_T)//'.'
         END IF
         LFN = LDATSET + LTXT_T + 1
C
         FILNAM = FILNAM(1:LFN)//SHELL(NCXRAY(IT))
         LFN = LFN + 1
         SPEC = '  '//SHELL(NCXRAY(IT))
C
         IF ( NCXRAY(IT).NE.1 ) THEN
            SPEC = SPEC(1:3)//SUBSH(LCXRAY(IT))
            LSP = 3 + LSUBSH(LCXRAY(IT))
            FILNAM = FILNAM(1:LFN)//SUBSHP(LCXRAY(IT))
            LFN = LFN + MAX(1,LSUBSH(LCXRAY(IT))-1)
         ELSE
            LSP = 3
         END IF
         REWIND 9
         FILNAM = FILNAM(1:LFN)//'.xps'
         SPEC = SPEC(1:LSP)//' - core level xps spectrum of '//TXT_T(IT)
     &          //' in  '//SYSTEM(1:LSYSTEM)
C
         OPEN (UNIT=7,FILE=FILNAM)
         WRITE (6,'(10X,A,A,/)') 'SPEC-FILE :  ( 7) ',FILNAM
C
         WRITE (6,'(A)') SPEC
         WRITE (6,99002) IT,NCXRAY(IT),LCXRAY(IT)
C
         NCST = 4*LCXRAY(IT) + 2
C
         IF ( ICALL.EQ.2 ) THEN
            CALL WRHEAD(7,FILNAM,'CLXPS     ',NCST)
C
            NTXRSGRP = 1
            WRITE (7,99014) 'NTXRSGRP  ',NTXRSGRP
            WRITE (7,99003) 'SPECTRUM  ',SPEC(1:LSP)
            WRITE (7,99014) 'IT        ',IT
            WRITE (7,99014) 'NCXRAY    ',NCXRAY(IT)
            WRITE (7,99014) 'LCXRAY    ',LCXRAY(IT)
         END IF
C
         IF ( NCST.GT.NCSTMAX ) THEN
            WRITE (6,99017) LCXRAY(IT),NCSTMAX
            STOP
         END IF
C
         CALL CORE(IPRINT,GCOR,FCOR,ECOR,SZCOR,KAPCOR,MM05COR,NKPCOR,
     &             IKMCOR,IZERO,ITXRAY,BCOR,BCORS,NCSTMAX)
C
         IF ( MEFORM.EQ.'GRV' ) THEN
            CALL DVDRSPLINE(VT(1,IT),R(1,IM),DVDR,JRWS(IM))
            CALL DVDRSPLINE(BT(1,IT),R(1,IM),DBDR,JRWS(IM))
         END IF
C
         IF ( ICALL.EQ.1 ) THEN
            IFIL = 6
            WRITE (IFIL,99004) NCST,(NKPCOR(ICST),ICST=1,NCST)
            WRITE (IFIL,99005)
         ELSE
            DO IFIL = 6,7
               WRITE (IFIL,99004) NCST,(NKPCOR(ICST),ICST=1,NCST)
               WRITE (IFIL,99005)
            END DO
         END IF
C
         NETAB(1) = NCST
C
         DO ICST = 1,NCST
C
            ETAB(ICST,1) = ECOR(ICST) + DCMPLX(EPHOT,EIMAG)
C
            IF ( DREAL(ETAB(ICST,1)).LT.EFERMI )
     &           CALL STOP_MESSAGE(ROUTINE,'set higher E_phot')
C
            DO K = 1,NKPCOR(ICST)
               DO N = 1,JRWS(IM)
                  RINT(N) = R2DRDI(N,IM)
     &                      *(GCOR(N,K,ICST)**2+FCOR(N,K,ICST)**2)
               END DO
               CALL RRADINT(IM,RINT,XNORM(K))
            END DO
C
            IF ( ICALL.EQ.1 ) THEN
               IFIL = 6
               WRITE (IFIL,99006) ICST,NCXRAY(IT),LCXRAY(IT),
     &                            KAPCOR(ICST),(2*MM05COR(ICST)+1),
     &                            IKMCOR(ICST,1),XNORM(1),ECOR(ICST),
     &                            ECOR(ICST)*RY_EV,SZCOR(ICST),
     &                            IZERO(ICST)
               IF ( NKPCOR(ICST).EQ.2 ) WRITE (IFIL,99007)
     &              IKMCOR(ICST,2),XNORM(2)
            ELSE
               DO IFIL = 6,7
                  WRITE (IFIL,99006) ICST,NCXRAY(IT),LCXRAY(IT),
     &                               KAPCOR(ICST),(2*MM05COR(ICST)+1),
     &                               IKMCOR(ICST,1),XNORM(1),ECOR(ICST),
     &                               ECOR(ICST)*RY_EV,SZCOR(ICST),
     &                               IZERO(ICST)
                  IF ( NKPCOR(ICST).EQ.2 ) WRITE (IFIL,99007)
     &                 IKMCOR(ICST,2),XNORM(2)
               END DO
            END IF
         END DO
C
C ======================================================================
         IF ( WRTAU .AND. ICALL.EQ.1 ) RETURN
      END IF
C ======================================================================
      WRTAU = .FALSE.
      WRTAUMQ = .FALSE.
C ======================================================================
C
      IQ = IQAT(1,IT)
      IM = IMT(IT)
      IRTOP = JRWS(IM)
      DO IR = 1,IRTOP
         WGTR(IR) = R2DRDI(IR,IM)
      END DO
C
      NLAM = NKMQ(IQ)
      C = CTL(IT,1)
C ------------------------------------ calculate angular matrix elements
C -------------------------------------------- and perform time reversal
      IF ( MEFORM.EQ.'ADA' ) THEN
C
         CALL MEADA_AMEADA_OLD_VERSION
C
         CALL AMETIMREV(A1_ADA,A2_ADA,2,NKMMAX)
C
C -------------------------------------------------- reverse sign for A2
C -------------------------------- to account for T( i f_LAM) = -i f_LAM
C
         A2_ADA(:,:,:) = -A2_ADA(:,:,:)
C
      ELSE
         CALL STOP_MESSAGE(ROUTINE,'MEFORM='//MEFORM//
     &                     ' not implemented')
      END IF
C
      IF ( MECHECK ) THEN
         DO IPOL = 1,3
            CP = '('//CHPOL1(IPOL)//')  '
            CALL RMATSTRUCT('AME  A1'//CP,A1_ADA(1,1,IPOL),NLAM,NKMMAX,
     &                      3,3,0,1D-8,7)
            CALL RMATSTRUCT('AME  A2'//CP,A2_ADA(1,1,IPOL),NLAM,NKMMAX,
     &                      3,3,0,1D-8,7)
C
         END DO
      END IF
      NGEO = 0
      DO IGEO = 1,NGEOMAX
         TETQPH(IGEO) = 0.0D0
         PHIQPH(IGEO) = 0.0D0
         TETEPH(IGEO) = 0.0D0
         PHIEPH(IGEO) = 0.0D0
         TETKEL(IGEO) = 0.0D0
         PHIKEL(IGEO) = 0.0D0
      END DO
      DO IGEO = 0,NGEOMAX
         STRQ = 'ANGQ'
         CALL STRING_ADD_N(STRQ,IGEO)
         LL = LEN_TRIM(STRQ)
         STRE = 'ANGE'//STRQ(5:LL)
         STRK = 'ANGK'//STRQ(5:LL)
         IF ( IGEO.EQ.0 ) LL = 4
         CALL SECTION_SET_REAL_ARRAY(STRQ,ANGQ,NIN,2,0,9999D0,0)
         FOUND = FOUND_REAL_ARRAY
         CALL SECTION_SET_REAL_ARRAY(STRE,ANGE,NIN,2,0,9999D0,0)
         FOUND = FOUND .AND. FOUND_REAL_ARRAY
         CALL SECTION_SET_REAL_ARRAY(STRK,ANGK,NIN,2,0,9999D0,0)
         FOUND = FOUND .AND. FOUND_REAL_ARRAY
         IF ( FOUND ) THEN
            NGEO = MAX(1,IGEO)
            TETQPH(NGEO) = ANGQ(1)
            PHIQPH(NGEO) = ANGQ(2)
            TETEPH(NGEO) = ANGE(1)
            PHIEPH(NGEO) = ANGE(2)
            TETKEL(NGEO) = ANGK(1)
            PHIKEL(NGEO) = ANGK(2)
            IF ( IGEO.EQ.0 ) EXIT
         ELSE IF ( IGEO.GT.0 ) THEN
            EXIT
         END IF
      END DO
C
      IF ( ANGRES .EQV. .FALSE. ) NGEO = 1
C
      IFIL = 6
      WRITE (IFIL,99018)
      WRITE (IFIL,99008) 'TIME REVERSAL MODE ',ITIMREV
      WRITE (IFIL,99016) 'USE t_ss           ',USETSS
      WRITE (IFIL,99016) 'ME-CHECK           ',MECHECK
      WRITE (IFIL,99009) 'ME-FORM            ',MEFORM
      DO IFIL = 6,7
         WRITE (IFIL,99015) 'E_PH (eV) ',EPHOTEV
         WRITE (IFIL,99014) 'NGEO      ',NGEO
C
         DO IGEO = 1,NGEO
            WRITE (IFIL,99014) 'IGEO      ',IGEO
            WRITE (IFIL,99015) 'TET/PHI Q ',TETQPH(IGEO),PHIQPH(IGEO)
            WRITE (IFIL,99015) 'TET/PHI E ',TETEPH(IGEO),PHIEPH(IGEO)
            WRITE (IFIL,99015) 'TET/PHI K ',TETKEL(IGEO),PHIKEL(IGEO)
         END DO
C
      END DO
      IF ( ANGRES .EQV. .FALSE. ) WRITE (6,*)
     &      'CLXPS--ANGLE INTEGRATED CASE'
C
C
      IF ( STRPOL ) THEN
         WRITE (6,*) 'LINEARLY POLARISED LIGHT'
C     For linear polarization results are  stored in arrays ME,W,fun
C     in the first position (corresponding to ipol=1)
      ELSE
         WRITE (6,*) 'CIRCURARLY POLARISED LIGHT'
      END IF
      WRITE (6,99018)
      WZ2 = SQRT(2.0D0)
      DO ICST = 1,NCST
         WRITE (6,'(A,10E12.5)') 'core E ',ECOR(ICST),EPHOT,ETAB(ICST,1)
     &                           ,EPHOT*RY_EV
      END DO
      REWIND (9)
C
C     CORE CORE CORE CORECORE CORE CORE CORECORE CORE CORE CORECORE CORE
C
      CALL READTAU(9,ERYD,0,NCST,TAUT(1,1,1),0,NT,.TRUE.,WKM1,WKM2,0,NQ,
     &             0,NKM,NKMMAX,IPRINT)
C
      DO ICST = 1,NCST
C
         ETAB(ICST,1) = ECOR(ICST) + DCMPLX(EPHOT,EIMAG)
C
         TAUT(:,:,:) = C0
C
         CALL READTAU(9,ERYD,ICST,NCST,TAUT(1,1,IT),IT,NT,.TRUE.,WKM1,
     &                WKM2,0,NQ,1,NKM,NKMMAX,IPRINT)
C
C--------------------------------------------------------- normal   case
C
         CALL SSITE(1,1,IFILCBWF,CALCINT,GETIRRSOL,ERYD,P,IPRINT,NKM,
     &              TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
C------------------------- read in regular conduction band wavefunctions
C
         CALL WAVFUN_READ_REL(IFILCBWF,IT,0,ZG,ZF,ZG,ZF,IRTOP,NCPLWF,
     &                        IKMCPLWF)
C
         IF ( USETSS ) TAUT(:,:,IT) = TSST(:,:,IT)
C
C     --------------------------------------------- PRINT ------ START -
         IF ( IPRINT.GT.0 ) WRITE (6,*) '  ENERGY : ',ERYD
C
         WRITE (6,99010) ICST,ERYD,IQAT(1,IT),IT
         IF ( IPRINT.GE.5 ) THEN
C
            CALL CMATSTRUCT('TAU(KAP,MUE)',TAUT(1,1,IT),NLAM,NKMMAX,3,3,
     &                      0,1D-8,6)
C
            WRITE (6,*) ' '
         END IF
C     --------------------------------------------- PRINT ------- END --
C
         EPS = SQRT((ERYD+C**2)/(2*ERYD+C**2))
C
         IF ( ABS(ETAB(ICST,1)-ERYD).GT.1D-5 )
     &        CALL STOP_MESSAGE(ROUTINE,'E consistency ??')
C
         WRITE (6,99011) ICST,DREAL(ERYD-EPHOT),ERYD
C
         WMAT(:,:,:) = C0
         MEA(:,:,:) = C0
C
         IKMC = IKMCOR(ICST,1)
C
C     ADA=ADA=ADA=ADA=ADA=ADA=ADA=ADA=ADA=ADA=ADA=ADA=ADA=ADA=ADA=ADA
         IF ( MEFORM.EQ.'ADA' ) THEN
C
            DO IKM = 1,NKM
C
               CALL ME_ADA_COR(NCST,NKPCOR,IKMCOR,GCOR,FCOR,
     &                         IKMCPLWF(1,IKM),IKM,IKM,NCPLWF(IKM),ZG,
     &                         ZF,IT,C,NPOL,MEA,WGTR,NCSTMAX,NPOLMAX)
            END DO
C
         ELSE
            WRITE (6,*) ' MEFORM = ',MEFORM
            CALL STOP_MESSAGE(ROUTINE,'MEFORM not implemented')
         END IF
C
C
         IF ( MECHECK ) THEN
C
            WMAT(:,:,:) = C0
C
            DO JCST = 1,NCST
               IKMC = IKMCOR(JCST,1)
               DO IPOL = 1,NPOL
                  DO IKM = 1,NLAM
                     WMAT(IKM,IKMC,IPOL) = MEA(JCST,IKM,IPOL)
                  END DO
               END DO
            END DO
C
            DO IPOL = 1,NPOL
               STRINP = 'ME('//CHPOL1(IPOL)//')  '//MEFORM//'  REG'
               CALL CMATSTRUCT(STRINP,WMAT(1,1,IPOL),NLAM,NKMMAX,3,3,0,
     &                         1D-8,67)
            END DO
            WRITE (67,*) 'CLXPS wrote regular matrix elements'
            WRITE (6,*) ' regular matrix elements written to 67'
         END IF
C     GEO GEO  GEO GEO  GEO GEO  GEO GEO  GEO GEO GEO GEO GEO GEO GEO
C
         DO IGEO = 1,NGEO
C
C
            IF ( ANGRES ) THEN
               K0X = SIN(TETKEL(IGEO)*PI/180D0)
     &               *COS(PHIKEL(IGEO)*PI/180D0)
               K0Y = SIN(TETKEL(IGEO)*PI/180D0)
     &               *SIN(PHIKEL(IGEO)*PI/180D0)
               K0Z = COS(TETKEL(IGEO)*PI/180D0)
               EVX = SIN(TETEPH(IGEO)*PI/180D0)
     &               *COS(PHIEPH(IGEO)*PI/180D0)
               EVY = SIN(TETEPH(IGEO)*PI/180D0)
     &               *SIN(PHIEPH(IGEO)*PI/180D0)
               EVZ = COS(TETEPH(IGEO)*PI/180D0)
C
               DIR(1) = -K0X
               DIR(2) = -K0Y
               DIR(3) = -K0Z
C
               CALL RVECNORM(3,DIR)
C
               CALL CALC_CHPLM(DIR(1),DIR(2),DIR(3),YLM,LCHPMAX,
     &                         NLMCHPMAX)
C
               WRITE (6,'(I5,A,3F10.4)') IGEO,' ->K0 ',K0X,K0Y,K0Z
               WRITE (6,'(I5,A,3F10.4)') IGEO,' ->EV ',EVX,EVY,EVZ
            END IF
C
C
C MSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMS
            DO MSM05 = -1,0
               MS = DBLE(MSM05) + 0.5D0
               IF ( IPRINT.GT.0 ) WRITE (6,*) ' CORE STATE',ICST,
     &              '  MS=',MS,'  GEO=',IGEO
               DO IPOLAR = 1,NPOL
                  RSUM(IPOLAR) = C0
               END DO
C LAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAML
               DO LAM = 1,NLAM
C
                  L = INT(SQRT(DBLE(LAM)/2.0D0)-0.0001D0)
                  JJ = L + SIGN(0.5D0,DBLE(LAM-(2*L*(L+1))-0.0001D0))
                  MJ = DBLE(LAM) - (L*2*(JJ+0.5D0)+JJ+1.0D0)
                  M = NINT(MJ+MS)
                  LM = L*(L+1) + M + 1
                  MMS = NINT(-MS+1.5D0)
C
C
                  IF ( ABS(CGC(LAM,MMS)).GT.1D-8 ) THEN
C
                     IF ( ANGRES ) THEN
                        YLMK = YLM(LM)
                        ICY = CI**L*CGC(LAM,MMS)*DCONJG(YLMK)
                     ELSE
                        ICY = CGC(LAM,MMS)
                     END IF
C
C
C     LAMPLAMP LAMPLAMP LAMPLAMP LAMPLAMP LAMPLAMP LAMPLAMP LAMPLAMP
                     DO LAMP = 1,NLAM
                        IF ( STRPOL ) THEN
C     index 3:  ipol= 1,2,3  ==  (+),(-),(z)
C     M(X) =   [  M(+) + M(-) ] / SQRT(2)
C     M(Y) = I*[ -M(+) + M(-) ] / SQRT(2)
C
                           MEX = (MEA(ICST,LAMP,1)+MEA(ICST,LAMP,2))/WZ2
                           MEY = CI*(-MEA(ICST,LAMP,1)+MEA(ICST,LAMP,2))
     &                           /WZ2
                           MEZ = MEA(ICST,LAMP,3)
C
                           ME(1) = EVX*MEX + EVY*MEY + EVZ*MEZ
                        ELSE
                           DO IPOL = 1,NPOL
                              ME(IPOL) = MEA(ICST,LAMP,IPOL)
                           END DO
                        END IF
                        IF ( STRPOL ) THEN
                           RSUM(1) = RSUM(1) + EPS*ICY*TAUT(LAMP,LAM,IT)
     &                               *ME(1)
                        ELSE
                           DO IPOLAR = 1,NPOL
                              RSUM(IPOLAR) = RSUM(IPOLAR)
     &                           + EPS*ICY*TAUT(LAMP,LAM,IT)*ME(IPOLAR)
                           END DO
                        END IF
C
                        IF ( IPRINT.GT.0 ) THEN
                           IFIL = 80 + MMS
                           DO IPOL = 1,NPOL
                              IF ( ABS(ME(IPOL)).GT.0.000000001D0 )
     &                             WRITE (IFIL,
     &                             '(2I3,'' ME '',4E11.4,2F6.3,8E9.2)')
     &                             LAM,LAMP,ME(1),ME(2),ME(3),ICY,
     &                             TAUT(LAMP,LAM,IT)
                           END DO
                        END IF
C
C     LAMPLAMP LAMPLAMP LAMPLAMP LAMPLAMP LAMPLAMP LAMPLAMP LAMPLAMP
                     END DO
                  END IF
C LAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAMLAML
               END DO
               IF ( STRPOL ) THEN
                  W(1) = DBLE(RSUM(1)*DCONJG(RSUM(1)))
                  WRITE (7,99012) ICST,MS,IGEO,ERYD,W(1)
                  IF ( IPRINT.GT.0 ) WRITE (6,99012) ICST,MS,IGEO,ERYD,
     &                 W(1)
                  IF ( MS.LT.0D0 ) THEN
                     FUN(ICST,1,1,IGEO) = W(1)
                  ELSE
                     FUN(ICST,2,1,IGEO) = W(1)
                     FUN(ICST,3,1,IGEO) = FUN(ICST,1,1,IGEO)
     &                  + FUN(ICST,2,1,IGEO)
                  END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
                  IF ( WRBUILDBOT ) WRITE (IFILBUILDBOT,99019)
     &                 ROUTINE(1:LEN_TRIM(ROUTINE)),ICST,MS,IGEO,ERYD,
     &                 W(1)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
               ELSE
                  DO IPOLAR = 1,NPOL
                     W(IPOLAR) = DBLE(RSUM(IPOLAR)*DCONJG(RSUM(IPOLAR)))
                     IF ( MS.LT.0D0 ) THEN
                        FUN(ICST,1,IPOLAR,IGEO) = W(IPOLAR)
                     ELSE
                        FUN(ICST,2,IPOLAR,IGEO) = W(IPOLAR)
                        FUN(ICST,3,IPOLAR,IGEO)
     &                     = FUN(ICST,1,IPOLAR,IGEO)
     &                     + FUN(ICST,2,IPOLAR,IGEO)
                     END IF
                  END DO
                  WRITE (7,99013) ICST,MS,IGEO,ERYD,W(1),W(2),W(3)
                  IF ( IPRINT.GT.0 ) WRITE (6,99013) ICST,MS,IGEO,ERYD,
     &                 W(1),W(2),W(3)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
                  IF ( WRBUILDBOT ) WRITE (IFILBUILDBOT,99019)
     &                 ROUTINE(1:LEN_TRIM(ROUTINE)),ICST,MS,IGEO,ERYD,
     &                 W(1),W(2),W(3)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
               END IF
C     MSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMSMS
            END DO
C     GEO GEO  GEO GEO GEO GEO GEO GEO GEO GEO GEO GEO GEO GEO GEO GEO
         END DO
C     CORE CORE CORE CORECORE CORE CORE CORECORE CORE CORE CORECORE CORE
      END DO
C
C
      DEALLOCATE (DBDR,BCOR,DVDR,RINT,WMAT)
      DEALLOCATE (W,BCORS,ME,MEA)
      DEALLOCATE (FUN,RSUM,YLM)
C
      STOP
C
C     EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEERROR
C
C=======================================================================
C
99001 FORMAT (' ',//,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*        ****   *              *     *  *****    ****        *'
     &  ,/,10X,
     &  '*       *    *  *               *   *   *    *  *    *       *'
     &  ,/,10X,
     &  '*       *       *                * *    *    *  *            *'
     &  ,/,10X,
     &  '*       *       *        ***      *     *****    ****        *'
     &  ,/,10X,
     &  '*       *       *                * *    *            *       *'
     &  ,/,10X,
     &  '*       *    *  *               *   *   *       *    *       *'
     &  ,/,10X,
     &  '*        ****   ******         *     *  *        ****        *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99002 FORMAT (/,4X,' CORE QUANTUM-NUMBERS  FOR  IT=',I2,':   N=',I2,
     &        '  L=',I2,/)
99003 FORMAT (A10,A)
99004 FORMAT (//,' CORE STATES :',//,' NCST:  ',I4,/,' NKPCOR:',20I4)
99005 FORMAT (/,' ICST  N   L  KAP  MUE  IKM     NORM   ',
     &        '     E(Ry)       E(eV)    <SIGMA_z>  I0')
99006 FORMAT (5I4,'/2',I4,F12.6,F12.4,F12.3,F12.4,I5)
99007 FORMAT (22X,I4,F12.6)
99008 FORMAT (A,I10)
99009 FORMAT (A,A)
99010 FORMAT (I3,' E=',2F6.3,'  IQ=',I2,'  IT=',I2)
99011 FORMAT (' core state ',I2,'  E(initial)=',F8.4,' RYD',
     &        '  E(final)=',2F8.4,' RYD')
99012 FORMAT (' ICST=',I2,'  MS=',F4.1,'  IGEO=',I2,'  E(FIN)=',2F10.5,
     &        ' RY      W=',E14.7)
99013 FORMAT (' ICST=',I2,'  MS=',F4.1,'  IGEO=',I2,'  E(FIN)=',2F10.5,
     &        ' RY     W(+,-,z)=',3(E14.7))
99014 FORMAT (A10,I10)
99015 FORMAT (A10,6F10.5)
99016 FORMAT (A,L3)
99017 FORMAT (//,60('*'),/,10X,' STOP in <CLXPS>',/,10X,'LCXRAY=',I3,
     &        ' too large for NCSTMAX=',I3)
99018 FORMAT (//,60('*'))
99019 FORMAT ('# BUILDBOT: ',A,':  CL-XPS intensity',
     &        '  for ICST, MS, IGEO, ERYD =',I5,F5.2,I5,2F10.6,/,
     &        (1PE22.14))
      END
