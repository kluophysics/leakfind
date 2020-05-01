C*==potfit.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE POTFIT
      IMPLICIT NONE
C*--POTFIT4
C
C*** Start of declarations rewritten by SPAG
C
C*** End of declarations rewritten by SPAG
C
C   ********************************************************************
C   *                                                                  *
C   *     import potential from other band structure packages          *
C   *                                                                  *
C   *     after reading the potential is interpolated to the           *
C   *     standard SPR-KKR mesh and then written to an auxilary file.  *
C   *     from there the standard subroutines read the data again.     *
C   *                                                                  *
C   ********************************************************************
C   *   subroutine to read in potential, field and grid data           *
C   *   from various sources -- the according format is found          *
C   *   by trial and error                                             *
C   *                                                                  *
C   *   minimum information in the input file is :                     *
C   *                                                                  *
C   *   - TITLE                                                        *
C   *   - NPT(IM)       NUMBER OF R-MESHPOINTS                         *
C   *   - RIN(IR,IM)    R-MESH                                         *
C   *   - VIN(IR,IT,IS) POTENTIAL                                      *
C   *   - RWS(IM)       WIGNER-SEITZ-RADIUS FOR ATOM IT                *
C   *                                                                  *
C   *   the input potential is converted to the internally used        *
C   *   form  V=0.5*(V(UP)+V(DN)) and  B=0.5*(V(UP)-V(DN))             *
C   *   for interpolation RIN may have arbitrary form and VIN is       *
C   *   expected to be the true potential, i.e. not V*R, not V*R**2 .. *
C   *                                                                  *
C   *   the potential is muffin-tin-ized if IMTZ=1                     *
C   *                                                                  *
C   ********************************************************************
C
      WRITE (6,*) '************************************************'
      WRITE (6,*) '************************************************'
      WRITE (6,*) '***      POTFIT NOT AVAILABLE see 6.0.2    *****'
      WRITE (6,*) '************************************************'
      WRITE (6,*) '************************************************'
      END
