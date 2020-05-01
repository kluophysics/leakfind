C*==calcdos_write.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE CALCDOS_WRITE(IEPATH)
C   ********************************************************************
C   *                                                                  *
C   *  write the DOS data to stanadrd files for plotting               *
C   *                                                                  *
C   *  MPI:  collect first the energy resolved DOS data                *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:ITBOT,ITTOP,DOBS_LTEX,DOBS_TEX_GLO,NLT,NTMAX,
     &                   CONC
      USE MOD_ENERGY,ONLY:NETAB,ETAB,NEMAX
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_ANGMOM,ONLY:IDOS,ISMT,IOMT,ISDM,IKDS,NLMAX,NOBSMAX
      USE MOD_FILES,ONLY:IFILDOS,IOTMP,LDATSET,DATSET
      USE MOD_RMESH,ONLY:FULLPOT
      USE MOD_CALCMODE,ONLY:THERMAL_VIBRA_FLUCT
      USE MOD_MPI,ONLY:MPI,MPI_ID
      USE MOD_THERMAL,ONLY:TEMP_LAT,NHAT_MAG
      IMPLICIT NONE
C*--CALCDOS_WRITE21
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IEPATH
C
C Local variables
C
      REAL*8 DDOT
      COMPLEX*16 ERYD
      CHARACTER*80 FILNAM
      INTEGER I,IE,IL,IPOL,IT,LFILDAT,N
      REAL*8 SPF,SPNPOL(:),SWAP(:)
      CHARACTER*5 STR5
c modified by XJQ: total spin-up / spin-down DOS
      real*8 sum_dos_spinup, sum_dos_spindown,
     &       sum_dos_spinup_it, sum_dos_spindown_it
c end-mod-xjq
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE SWAP,SPNPOL
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C                    collect energy resolved DOS data
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      IF ( MPI ) THEN
C
         CALL DRV_MPI_BARRIER
C
         N = 4*NOBSMAX*NLMAX*NTMAX*NEMAX
         ALLOCATE (SWAP(N))
C
         CALL DRV_MPI_REDUCE_R(DOBS_LTEX(0,1,1,1,1),SWAP(1),N)
C
         DEALLOCATE (SWAP)
C
C----------------------------------------------------------- TEMPERATURE
         IF ( THERMAL_VIBRA_FLUCT ) THEN
C
            N = 4*NOBSMAX*NTMAX*NEMAX
            ALLOCATE (SWAP(N))
C
            CALL DRV_MPI_REDUCE_R(DOBS_TEX_GLO(0,1,1,1),SWAP(1),N)
C
            DEALLOCATE (SWAP)
C
         END IF
C----------------------------------------------------------- TEMPERATURE
C
         IF ( MPI_ID.NE.0 ) RETURN
C
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C
C ======================================================================
C ======================================================================
C                              write DOS
C ======================================================================
C ======================================================================
C
      IF ( IREL.GT.1 ) THEN
         SPF = 0.5D0
      ELSE
         SPF = 1.0D0
      END IF
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
      DO IE = 1,NETAB(IEPATH)
C
         ERYD = ETAB(IE,IEPATH)
C
C---------------------------------------------------- (l,s)-resolved DOS
C
         WRITE (IFILDOS,99001) ERYD,(((SPF*(DOBS_LTEX(0,IDOS,IL,IT,IE)+
     &                         DOBS_LTEX(0,ISMT,IL,IT,IE))),IL=1,NLT(IT)
     &                         ),((SPF*(DOBS_LTEX(0,IDOS,IL,IT,IE)-
     &                         DOBS_LTEX(0,ISMT,IL,IT,IE))),IL=1,NLT(IT)
     &                         ),IT=ITBOT,ITTOP)
C
         IF ( FULLPOT ) CYCLE
C
C----------------------------------------------------kappa-resolved DOS
C
         WRITE (13,99002) DREAL(ERYD),
     &                    (((DOBS_LTEX(0,IDOS,IL,IT,IE)),IL=1,NLT(IT)),
     &                    (DOBS_LTEX(0,IDOS,1,IT,IE)),
     &                    (0.5D0*(DOBS_LTEX(0,IDOS,IL,IT,IE)
     &                    -DOBS_LTEX(0,IKDS,IL,IT,IE)),
     &                    0.5D0*(DOBS_LTEX(0,IDOS,IL,IT,IE)
     &                    +DOBS_LTEX(0,IKDS,IL,IT,IE)),IL=2,NLT(IT)),
     &                    IT=ITBOT,ITTOP)
C
C----------------------------------------------- spin dipolar moment T_z
C
         WRITE (14,99002) DREAL(ERYD),
     &                    ((DOBS_LTEX(0,ISMT,IL,IT,IE),IL=1,NLT(IT)),
     &                    (DOBS_LTEX(0,IOMT,IL,IT,IE),IL=1,NLT(IT)),
     &                    (DOBS_LTEX(0,ISDM,IL,IT,IE),IL=1,NLT(IT)),
     &                    IT=ITBOT,ITTOP)
C
      END DO
c modified by XJQ: total spin-up / spin-down DOS
      write(ifildos,'(a)') 'sum of dos on il and it:'
      DO IE = 1,NETAB(IEPATH)
        ERYD = ETAB(IE,IEPATH)
        sum_dos_spinup = 0d0
        sum_dos_spindown = 0d0
        DO IT = ITBOT, ITTOP
          sum_dos_spinup_it = 0d0
          sum_dos_spindown_it = 0d0
          DO IL=1,NLT(IT)
            sum_dos_spinup_it = sum_dos_spinup_it + 
     &                          SPF*(DOBS_LTEX(0,IDOS,IL,IT,IE)+
     &                               DOBS_LTEX(0,ISMT,IL,IT,IE))
            sum_dos_spindown_it = sum_dos_spindown_it + 
     &                            SPF*(DOBS_LTEX(0,IDOS,IL,IT,IE)-
     &                                 DOBS_LTEX(0,ISMT,IL,IT,IE))
          ENDDO
          sum_dos_spinup = sum_dos_spinup + sum_dos_spinup_it * conc(it)
          sum_dos_spindown = sum_dos_spindown +
     &                       sum_dos_spindown_it * conc(it)
        ENDDO
        write(ifildos,'(4e13.6)') eryd, sum_dos_spinup, sum_dos_spindown
      END DO
c end-mod-xjq
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
C----------------------------------------------------------- TEMPERATURE
      IF ( .NOT.THERMAL_VIBRA_FLUCT ) RETURN
C----------------------------------------------------------- TEMPERATURE
C
C
CDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
C                     average spin-projected DOS
CDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
C
      ALLOCATE (SPNPOL(NTMAX))
C
      IF ( TEMP_LAT.GT.0D0 ) THEN
         I = INT(LOG10(TEMP_LAT))
      ELSE
         I = 0
      END IF
      WRITE (STR5,'(I5)') INT(TEMP_LAT)
      FILNAM = DATSET(1:LDATSET)//'_DOS_NHAT_'//STR5((5-I):5)//'K.dat'
      LFILDAT = LDATSET + 10 + (I+1) + 5
      OPEN (IOTMP,FILE=FILNAM(1:LFILDAT))
C
      WRITE (IOTMP,99004) NHAT_MAG,ITBOT,ITTOP,TEMP_LAT
C
      IPOL = 3
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
      DO IE = 1,NETAB(IEPATH)
C
         ERYD = ETAB(IE,IEPATH)
C
C
         DO IT = ITBOT,ITTOP
            SPNPOL(IT) = DDOT(3,DOBS_TEX_GLO(1,ISMT,IT,IE),1,NHAT_MAG,1)
         END DO
C
         WRITE (IOTMP,99003) DREAL(ERYD),
     &                       (((DOBS_TEX_GLO(IPOL,IDOS,IT,IE))
     &                       +SPNPOL(IT))/2D0,
     &                       ((DOBS_TEX_GLO(IPOL,IDOS,IT,IE))-SPNPOL(IT)
     &                       )/2D0,IT=ITBOT,ITTOP)
C
      END DO
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
      CLOSE (IOTMP)
CDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
C
99001 FORMAT (8E10.4,:,/,(10X,7E10.4))
99002 FORMAT (10E12.5)
99003 FORMAT (200E15.6)
99004 FORMAT ('#',/,
     &        '#  spin projected DOS:  (DOS +/- NHAT_MAG*VEC_SIGMA)/2',
     &        /,'#  NHAT_MAG = (',F7.3,',',F7.3,',',F7.3,')',/,
     &        '#  IT =',I4,' ...',I4,/,'#  temperature T =',F6.1,' K',/,
     &        '#')
C
      END
