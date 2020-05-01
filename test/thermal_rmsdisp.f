C*==thermal_rmsdisp.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE THERMAL_RMSDISP(TEMP_LAT,T_DEBYE,MASS,RMSU)
C   ********************************************************************
C   *                                                                  *
C   *  subroutine to calculate the root-mean-square displacement RMSU  *
C   *  of atoms at certain temperature TEMP_LAT within the Debye model *
C   *                                                                  *
C   *               9.D0*HBAR**2         F(T_Debye/T)       1          *
C   *  RMSU**2 = ----------------- * ( --------------  +  ---  )       *
C   *              KB*MASS*T_Debye        T_Debye/T         2          *
C   *                                                                  *
C   *       F(x) - Debye function                                      *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:HBAR_SI,KB_SI,A0_SI,MP_SI
      USE MOD_SITES,ONLY:NQ,ITOQ,NOQ
      USE MOD_TYPES,ONLY:CONC,Z
      USE MOD_FILES,ONLY:IPRINT
      IMPLICIT NONE
C*--THERMAL_RMSDISP20
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 CONST_DF
      PARAMETER (CONST_DF=9.D0*HBAR_SI**2/KB_SI/MP_SI)
C
C Dummy arguments
C
      REAL*8 MASS,RMSU,TEMP_LAT,T_DEBYE
C
C Local variables
C
      REAL*8 COEFF,MASTAB(0:100),RMSUSQ,T_DEBYE_INP,T_DEBYE_TAB(0:100),X
      LOGICAL INITIALIZE
      INTEGER IO,IQ,IT
      REAL*8 THERMAL_DEBYE1
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
      DATA MASTAB/0.0D0,1.008D0,4.003D0,6.941D0,9.012D0,10.810D0,
     &     12.010D0,14.010D0,15.999D0,18.998D0,20.180D0,22.990D0,
     &     24.305D0,26.982D0,28.086D0,30.974D0,32.074D0,35.453D0,
     &     39.948D0,39.090D0,40.080D0,44.956D0,47.900D0,50.942D0,
     &     52.000D0,54.938D0,55.850D0,58.930D0,58.710D0,63.550D0,
     &     65.380D0,69.720D0,72.590D0,74.922D0,78.960D0,79.910D0,
     &     83.800D0,85.470D0,87.620D0,88.910D0,91.220D0,92.910D0,
     &     95.940D0,98.910D0,101.070D0,102.900D0,106.400D0,107.870D0,
     &     112.400D0,114.820D0,118.690D0,121.750D0,127.600D0,126.900D0,
     &     131.300D0,132.910D0,137.340D0,139.910D0,140.120D0,140.910D0,
     &     144.240D0,145.000D0,150.350D0,151.960D0,157.250D0,158.920D0,
     &     162.500D0,164.930D0,167.260D0,168.930D0,173.040D0,174.970D0,
     &     178.490D0,180.950D0,183.850D0,186.200D0,190.200D0,192.220D0,
     &     195.090D0,196.970D0,200.590D0,204.370D0,207.190D0,208.980D0,
     &     210.000D0,0.000D0,0.000D0,0.000D0,0.000D0,0.000D0,0.000D0,
     &     0.000D0,0.000D0,0.000D0,0.000D0,0.000D0,0.000D0,0.000D0,
     &     0.000D0,0.000D0,0.000D0/
      DATA T_DEBYE_TAB/0.0D0,110.0D0,26.0D0,400.0D0,1000.0D0,1250.0D0,
     &     1860.0D0,79.0D0,46.0D0,0.0D0,63.0D0,150.0D0,318.0D0,394.0D0,
     &     625.0D0,0.0D0,0.0D0,0.0D0,85.0D0,100.0D0,230.0D0,359.0D0,
     &     380.0D0,390.0D0,460.0D0,400.0D0,420.0D0,385.0D0,375.0D0,
     &     315.0D0,234.0D0,240.0D0,360.0D0,285.0D0,150.0D0,0.0D0,73.0D0,
     &     56.0D0,147.0D0,256.0D0,250.0D0,275.0D0,380.0D0,0.0D0,382.0D0,
     &     350.0D0,275.0D0,215.0D0,120.0D0,129.0D0,170.0D0,200.0D0,
     &     139.0D0,0.0D0,55.0D0,40.0D0,110.0D0,132.0D0,139.0D0,152.0D0,
     &     157.0D0,0.0D0,166.0D0,107.0D0,176.0D0,188.0D0,186.0D0,
     &     191.0D0,195.0D0,200.0D0,118.0D0,207.0D0,0.0D0,225.0D0,
     &     310.0D0,416.0D0,400.0D0,430.0D0,230.0D0,170.0D0,100.0D0,
     &     0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &     0.0D0,0.0D0,0.D00,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/
C
C-----------------------------------------------------------------------
C
      T_DEBYE_INP = T_DEBYE
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
      IF ( INITIALIZE ) THEN
C------------------------------------- evaluate average mass and T_Debye
         MASS = 0.D0
         T_DEBYE = 0.D0
         X = 0.0D0
         IF ( IPRINT.GT.0 ) WRITE (6,99003)
         DO IQ = 1,NQ
            DO IO = 1,NOQ(IQ)
               IT = ITOQ(IO,IQ)
               IF ( Z(IT).GT.0 ) THEN
                  MASS = MASS + CONC(IT)*MASTAB(Z(IT))
                  T_DEBYE = T_DEBYE + CONC(IT)*T_DEBYE_TAB(Z(IT))
                  X = X + CONC(IT)
               END IF
               IF ( IPRINT.GT.0 ) WRITE (6,99002) IQ,IT,Z(IT),
     &              MASTAB(Z(IT)),T_DEBYE_TAB(Z(IT))
C
            END DO
         END DO
         MASS = MASS/X
         T_DEBYE = T_DEBYE/X
C
         IF ( T_DEBYE_INP.GT.0.1D0 ) T_DEBYE = T_DEBYE_INP
C
         IF ( IPRINT.GT.0 ) WRITE (6,99001) MASS,T_DEBYE
C
         INITIALIZE = .FALSE.
C
         RETURN
      END IF
C=======================================================================
C
      IF ( T_DEBYE_INP.GT.0.1D0 ) T_DEBYE = T_DEBYE_INP
C
      X = T_DEBYE/TEMP_LAT
C
C------------- accounting for zero vibrations --------------------------
      COEFF = (THERMAL_DEBYE1(X)/X + 0.25d0)
C-----------------------------------------------------------------------
C
C      COEFF = THERMAL_DEBYE1(X)/X
C
      RMSUSQ = COEFF*CONST_DF/MASS/T_DEBYE
C
      RMSU = SQRT(RMSUSQ)/A0_SI
C
      WRITE (6,99001) MASS,T_DEBYE,TEMP_LAT,RMSU
C
99001 FORMAT (/,10X,'average mass                  ',F10.3,' [m_p]',/,
     &        10X,'average Debye temperature     ',F10.3,' [K]',/,:,10X,
     &        'temperature                   ',F10.3,' [K]',/,:,10X,
     &        'root-mean-square displacement ',F10.3,' [a.u.]',/)
99002 FORMAT (10X,'IQ=',I4,'  IT=',I4,'   Z(t)=',I4,'   M(t)=',F6.1,
     &        '    T_DEBYE(t)=',F6.1)
99003 FORMAT (/,10X,'evaluating the avarage mass MASS and',
     &        ' estimate for T_Debye',/)
      END
