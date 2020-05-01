C*==aa0001.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      IMPLICIT NONE
C*--AA00013
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      REAL*8 A_GAUSS,EPHOTON_EV,EPHOTON_RY,FWHM,GAUSS,PI,RY_EV,
     &       SIGMA_GAUSS,TIM,TIM_0_GAUSS,TIM_PERIOD,TIM_RANGE
      INTEGER ITIM,NTIM
C
C*** End of declarations rewritten by SPAG
C
      RY_EV = 13.6056981D0
      PI = 3.1415926535897932384626433D0
C
      NTIM = 1000
C
      EPHOTON_EV = 3
      EPHOTON_RY = EPHOTON_EV/RY_EV
C
      WRITE (6,99002) 'E_photon     (ev)  ',EPHOTON_EV
      WRITE (6,99002) 'E_photon     (Hz)  ',EPHOTON_EV*2.41799*10D14
C
      TIM_PERIOD = 2*PI/EPHOTON_RY
      TIM_RANGE = 40*TIM_PERIOD
C
      TIM_0_GAUSS = 10*TIM_PERIOD
      FWHM = 5*TIM_PERIOD
C
      SIGMA_GAUSS = FWHM/(2*SQRT(2*LOG(2D0)))
      A_GAUSS = 1D0/(2*SIGMA_GAUSS**2)
C
      DO ITIM = 1,NTIM
C
         TIM = TIM_RANGE*(ITIM-1)/DBLE(NTIM-1)
C
         GAUSS = EXP(-A_GAUSS*(TIM-TIM_0_GAUSS)**2)
         WRITE (8,99001) TIM,COS(TIM*EPHOTON_RY)*GAUSS
C
      END DO
C
99001 FORMAT (2E20.8)
99002 FORMAT (A,3X,2E20.8)
C
      END
C
