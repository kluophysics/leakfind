C*==me_grv_grv.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE ME_GRV_GRV(IFIL_RHSA,ERYDA,IFIL_RHSB,ERYDB,
     &                      ME_CC_BRA_RWF,IT,MZAZB,MZBZA,MIRR_2,MIRR_3,
     &                      MIRR_4,C)
C   ********************************************************************
C   *                                                                  *
C   *    read wave function and the calculate matrix elements          *
C   *                                                                  *
C   *         MZAZB  =   < Z^+(E_a) | H_lam | Z(E_b) >                 *
C   *                                                                  *
C   *         MZBZA  =   < Z^+(E_b) | H_lam | Z(E_a) >                 *
C   *                                                                  *
C   *         MIRR_2 =   < Z^+(E_b) | H_lam * IZAZB | J(E_a) >         *
C   *                  + < Z^+(E_b) | H_lam * IJAZB | Z(E_a) >         *
C   *                                                                  *
C   *         MIRR_3 =   < J^+(E_b) | H_lam * IZAZB | J(E_a) >         *
C   *                  + < Z^+(E_b) | H_lam * IZAJB | Z(E_a) >         *
C   *                                                                  *
C   *         MIRR_4 =   < J^+(E_b) | H_lam * IZAZB | J(E_a) >         *
C   *                  + < Z^+(E_b) | H_lam * IZAJB | Z(E_a) >         *
C   *                                                                  *
C   *     with    H_lam = Grad V * A_lam                               *
C   *                                                                  *
C   *     and the matrix element integral functions                    *
C   *                                                                  *
C   *      IXAYB(r) = < X^+(E_a) | H_lam | Y(E_b) >_r                  *
C   *                                                                  *
C   *  X, Y: = Z (J) regular (irregular) solutions to Dirac equation   *
C   *                                                                  *
C   *    polarisation lam:        ipol = 1,2,3  ==  (+),(-),(z)        *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *  ME_CC_BRA_RWF = .T.   take complex conjugate                    *
C   *                        <BRA| radial wave function                *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:IMT,NCPLWFMAX,NLMFP,NLMFPT,KLMFP,VT,BT,VNST,
     &    BNST,IKMCPLWF_LA,IKMCPLWF_LB,IKMCPLWF_RA,IKMCPLWF_RB
      USE MOD_RMESH,ONLY:JRWS,NRMAX,JRCRI,JRNS1,FULLPOT,R2DRDI
      USE MOD_CONSTANTS,ONLY:CI,C0,Y00,SQRT_2
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,NCPLWF_LA,NCPLWF_LB,NCPLWF_RA,
     &    NCPLWF_RB,NPOL,AME_G,AG_RGNT,NKM_EXT,A_SIGY,NLM_AME_RLM_EXT,
     &    ISMT
      IMPLICIT NONE
C*--ME_GRV_GRV47
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ME_GRV_GRV')
C
C Dummy arguments
C
      REAL*8 C
      COMPLEX*16 ERYDA,ERYDB
      INTEGER IFIL_RHSA,IFIL_RHSB,IT
      LOGICAL ME_CC_BRA_RWF
      COMPLEX*16 MIRR_2(NKMMAX,NKMMAX,3,3),MIRR_3(NKMMAX,NKMMAX,3,3),
     &           MIRR_4(NKMMAX,NKMMAX,3,3),MZAZB(NKMMAX,NKMMAX,3),
     &           MZBZA(NKMMAX,NKMMAX,3)
C
C Local variables
C
      REAL*8 BLM(:,:),GBLM(:,:,:),GVLM(:,:,:),GX,GY,VLM(:,:)
      COMPLEX*16 B_JFRA(:,:,:),B_JFRB(:,:,:),B_JGRA(:,:,:),B_JGRB(:,:,:)
     &           ,B_ZFRA(:,:,:),B_ZFRB(:,:,:),B_ZGRA(:,:,:),
     &           B_ZGRB(:,:,:),GB_JFRA(:,:,:,:),GB_JFRB(:,:,:,:),
     &           GB_JGRA(:,:,:,:),GB_JGRB(:,:,:,:),GB_ZFRA(:,:,:,:),
     &           GB_ZFRB(:,:,:,:),GB_ZGRA(:,:,:,:),GB_ZGRB(:,:,:,:),
     &           GSBLM(:,:,:),GSVLM(:,:,:),GV_JFRA(:,:,:,:),
     &           GV_JFRB(:,:,:,:),GV_JGRA(:,:,:,:),GV_JGRB(:,:,:,:),
     &           GV_ZFRA(:,:,:,:),GV_ZFRB(:,:,:,:),GV_ZGRA(:,:,:,:),
     &           GV_ZGRB(:,:,:,:),IJAJB(:),IJAZB(:),IZAJB(:),IZAZB(:),
     &           JFLA(:,:,:),JFLB(:,:,:),JFRA(:,:,:),JFRB(:,:,:),
     &           JGLA(:,:,:),JGLB(:,:,:),JGRA(:,:,:),JGRB(:,:,:),
     &           OMEGA_AB,OMEGA_BA,PAB_GRV,PBA_GRV,TMC,TMCSQ,
     &           V_JFRA(:,:,:),V_JFRB(:,:,:),V_JGRA(:,:,:),V_JGRB(:,:,:)
     &           ,V_ZFRA(:,:,:),V_ZFRB(:,:,:),V_ZGRA(:,:,:),
     &           V_ZGRB(:,:,:),ZFLA(:,:,:),ZFLB(:,:,:),ZFRA(:,:,:),
     &           ZFRB(:,:,:),ZGLA(:,:,:),ZGLB(:,:,:),ZGRA(:,:,:),
     &           ZGRB(:,:,:)
      INTEGER I,IA_ERR,IFIL_LHSA,IFIL_LHSB,IKM,IM,IPOL,IR,IRTOP,J,JPOL,
     &        K,KGSVLM(:,:),KGVLM(:,:),KGV_ZGR(:,:,:),KVLM(:),LAMA1,
     &        LAMA2,LAMB3,LAMB4,LM,N
      LOGICAL INITIALIZE,K_SELECTION_RULES
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
C
      ALLOCATABLE IZAZB,IZAJB,IJAZB,IJAJB,KGV_ZGR
      ALLOCATABLE JFLA,JGLA,ZGLA,ZFLA,JGRB,JFRB,ZFRB,ZGRB
      ALLOCATABLE JFRA,JGRA,ZGRA,ZFRA,JGLB,JFLB,ZFLB,ZGLB
      ALLOCATABLE VLM,BLM,GVLM,GBLM,GSBLM,GSVLM,KVLM,KGVLM,KGSVLM
C
      ALLOCATABLE V_ZGRA,GV_ZGRA,V_ZFRA,GV_ZFRA
      ALLOCATABLE V_JFRA,GV_JFRA,V_JGRA,GV_JGRA
      ALLOCATABLE B_ZGRA,GB_ZGRA,B_ZFRA,GB_ZFRA
      ALLOCATABLE B_JFRA,GB_JFRA,B_JGRA,GB_JGRA
C
      ALLOCATABLE V_ZGRB,GV_ZGRB,V_ZFRB,GV_ZFRB
      ALLOCATABLE V_JFRB,GV_JFRB,V_JGRB,GV_JGRB
      ALLOCATABLE B_ZGRB,GB_ZGRB,B_ZFRB,GB_ZFRB
      ALLOCATABLE B_JFRB,GB_JFRB,B_JGRB,GB_JGRB
C
      N = NKM
C
      ALLOCATE (ZGLA(NRMAX,NCPLWFMAX,N),ZFLA(NRMAX,NCPLWFMAX,N))
      ALLOCATE (JGLA(NRMAX,NCPLWFMAX,N),JFLA(NRMAX,NCPLWFMAX,N))
      ALLOCATE (ZGRB(NRMAX,NCPLWFMAX,N),ZFRB(NRMAX,NCPLWFMAX,N))
      ALLOCATE (JGRB(NRMAX,NCPLWFMAX,N),JFRB(NRMAX,NCPLWFMAX,N))
C
      ALLOCATE (JFRA(NRMAX,NCPLWFMAX,N),JGRA(NRMAX,NCPLWFMAX,N))
      ALLOCATE (ZGRA(NRMAX,NCPLWFMAX,N),ZFRA(NRMAX,NCPLWFMAX,N))
      ALLOCATE (JGLB(NRMAX,NCPLWFMAX,N),JFLB(NRMAX,NCPLWFMAX,N))
      ALLOCATE (ZFLB(NRMAX,NCPLWFMAX,N),ZGLB(NRMAX,NCPLWFMAX,N))
C
      ALLOCATE (V_ZGRA(NRMAX,NCPLWFMAX,NKM),GV_ZGRA(NRMAX,N,N,3))
      ALLOCATE (V_ZFRA(NRMAX,NCPLWFMAX,NKM),GV_ZFRA(NRMAX,N,N,3))
      ALLOCATE (V_JFRA(NRMAX,NCPLWFMAX,NKM),GV_JFRA(NRMAX,N,N,3))
      ALLOCATE (V_JGRA(NRMAX,NCPLWFMAX,NKM),GV_JGRA(NRMAX,N,N,3))
      ALLOCATE (B_ZGRA(NRMAX,NCPLWFMAX,NKM),GB_ZGRA(NRMAX,N,N,3))
      ALLOCATE (B_ZFRA(NRMAX,NCPLWFMAX,NKM),GB_ZFRA(NRMAX,N,N,3))
      ALLOCATE (B_JFRA(NRMAX,NCPLWFMAX,NKM),GB_JFRA(NRMAX,N,N,3))
      ALLOCATE (B_JGRA(NRMAX,NCPLWFMAX,NKM),GB_JGRA(NRMAX,N,N,3))
C
      ALLOCATE (V_ZGRB(NRMAX,NCPLWFMAX,NKM),GV_ZGRB(NRMAX,N,N,3))
      ALLOCATE (V_ZFRB(NRMAX,NCPLWFMAX,NKM),GV_ZFRB(NRMAX,N,N,3))
      ALLOCATE (V_JFRB(NRMAX,NCPLWFMAX,NKM),GV_JFRB(NRMAX,N,N,3))
      ALLOCATE (V_JGRB(NRMAX,NCPLWFMAX,NKM),GV_JGRB(NRMAX,N,N,3))
      ALLOCATE (B_ZGRB(NRMAX,NCPLWFMAX,NKM),GB_ZGRB(NRMAX,N,N,3))
      ALLOCATE (B_ZFRB(NRMAX,NCPLWFMAX,NKM),GB_ZFRB(NRMAX,N,N,3))
      ALLOCATE (B_JFRB(NRMAX,NCPLWFMAX,NKM),GB_JFRB(NRMAX,N,N,3))
      ALLOCATE (B_JGRB(NRMAX,NCPLWFMAX,NKM),GB_JGRB(NRMAX,N,N,3))
C
      ALLOCATE (KGV_ZGR(NKM,NKM,3))
C
      ALLOCATE (IZAZB(NRMAX))
      ALLOCATE (IZAJB(NRMAX))
      ALLOCATE (IJAZB(NRMAX))
      ALLOCATE (IJAJB(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate IJAJB')
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ( INITIALIZE ) THEN
         CALL AME_INIT('GRV       ',0)
         INITIALIZE = .FALSE.
      END IF
Cc      if( it .gt. 1 ) stop
C      ab1_grv=0d0
C      ab2_grv=0d0
C      ac1_grv=0d0
C      ac2_grv=0d0
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
      IF ( .NOT.ALLOCATED(A_SIGY) ) THEN
C
         ALLOCATE (A_SIGY(NKM_EXT,NKM_EXT,3,NLM_AME_RLM_EXT))
C
         N = NKM_EXT
C
         DO IPOL = 1,3
C
            IF ( IPOL.EQ.1 ) THEN
               JPOL = 2
            ELSE IF ( IPOL.EQ.2 ) THEN
               JPOL = 3
            ELSE
               JPOL = 1
            END IF
C
            DO LM = 1,NLM_AME_RLM_EXT
               A_SIGY(1:N,1:N,JPOL,LM)
     &            = MATMUL(AME_G(1:N,1:N,IPOL,ISMT),AG_RGNT(1:N,1:N,LM))
            END DO
C
         END DO
C
         A_SIGY(1:N,1:N,1,:) = -A_SIGY(1:N,1:N,1,:)
C
      END IF
C
C=======================================================================
C
C         set the chanel number for the LHS wave functions
C
      CALL SET_IFIL_LHS(IFIL_RHSA,IFIL_LHSA)
C
      CALL SET_IFIL_LHS(IFIL_RHSB,IFIL_LHSB)
C
C=======================================================================
C
      IM = IMT(IT)
      IF ( FULLPOT ) THEN
         IRTOP = JRCRI(IM)
      ELSE
         IRTOP = JRWS(IM)
         NLMFPT(IT) = 1
         KLMFP(1,IT) = 1
      END IF
C
      OMEGA_BA = ERYDB - ERYDA
      OMEGA_AB = ERYDA - ERYDB
      PBA_GRV = CI/(OMEGA_BA*(ERYDB+ERYDA+C**2))
      PAB_GRV = CI/(OMEGA_AB*(ERYDB+ERYDA+C**2))
      TMC = 2.0D0*0.5D0*C
      TMCSQ = TMC*C
C
C=======================================================================
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      APPLY GRADIENT TO UNCONVLUTED POTENTIAL (NO SHAPE FUNCTION)
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C------------------------------ supply potential functions and gradients
C
      ALLOCATE (VLM(NRMAX,NLMFP),GVLM(NRMAX,NLMFP,3))
      ALLOCATE (BLM(NRMAX,NLMFP),GBLM(NRMAX,NLMFP,3))
      ALLOCATE (GSVLM(NRMAX,NLMFP,3))
      ALLOCATE (GSBLM(NRMAX,NLMFP,3))
      GVLM(:,:,:) = 0D0
      GBLM(:,:,:) = 0D0
      VLM(:,:) = 0D0
      BLM(:,:) = 0D0
C
      DO IR = 1,IRTOP
         VLM(IR,1) = VT(IR,IT)/Y00
         BLM(IR,1) = BT(IR,IT)/Y00
      END DO
C
      DO LM = 2,NLMFPT(IT)
         IF ( KLMFP(LM,IT).NE.0 ) THEN
            DO IR = JRNS1(IM),IRTOP
               VLM(IR,LM) = VNST(IR,LM,IT)
               BLM(IR,LM) = BNST(IR,LM,IT)
            END DO
         END IF
      END DO
C
      ALLOCATE (KVLM(NLMFP),KGVLM(NLMFP,3),KGSVLM(NLMFP,3))
C
      KVLM(:) = 0
      KVLM(1:NLMFPT(IT)) = KLMFP(1:NLMFPT(IT),IT)
      KGVLM(1:NLMFP,1:3) = 1
C
C=======================================================================
C               gradient  GLM  in real spherical harmonics
C=======================================================================
C
      CALL GRAD_FLM(VLM,KVLM,GVLM,KGVLM,NLMFP,IM)
      WRITE (*,'(A,20I3)') '### KVLM     ',KVLM(1:9)
      WRITE (*,'(A,20I3)') '### KGVLM  1:',KGVLM(1:9,1)
      WRITE (*,'(A,20I3)') '### KGVLM  2:',KGVLM(1:9,2)
      WRITE (*,'(A,20I3)') '### KGVLM  3:',KGVLM(1:9,3)
C
      CALL GRAD_FLM(BLM,KVLM,GBLM,KGVLM,NLMFP,IM)
C
C=======================================================================
C             multiply gradient  GLM  with shape functions
C=======================================================================
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C                 TO BE ADDED FOR FULL POTENTIAL
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C=======================================================================
C     multiply factor  R2DRDI  for integration to ALL potential terms
C=======================================================================
C
      KGSVLM(:,:) = 0
      GSVLM(IRTOP:NRMAX,:,:) = C0
      GSBLM(IRTOP:NRMAX,:,:) = C0
C
      DO LM = 1,NLMFP
C
         IF ( KVLM(LM).NE.0 ) THEN
            DO IR = 1,IRTOP
               VLM(IR,LM) = VLM(IR,LM)*R2DRDI(IR,IM)
               BLM(IR,LM) = BLM(IR,LM)*R2DRDI(IR,IM)
            END DO
         END IF
C
         DO IPOL = 1,3
            IF ( KGVLM(LM,IPOL).NE.0 ) THEN
               DO IR = 1,IRTOP
                  GVLM(IR,LM,IPOL) = GVLM(IR,LM,IPOL)*R2DRDI(IR,IM)
                  GBLM(IR,LM,IPOL) = GBLM(IR,LM,IPOL)*R2DRDI(IR,IM)
               END DO
            END IF
         END DO
C
C???????????????????????????????????????????????????????????????????????
C                index 3:  ipol= 1,2,3  ==  (+),(-),(z)
C                 M(X) =   [  M(+) + M(-) ] / SQRT(2)
C                 M(Y) = I*[ -M(+) + M(-) ] / SQRT(2)
C
C
         IF ( KGVLM(LM,1).NE.0 .OR. KGVLM(LM,2).NE.0 ) THEN
C
            KGSVLM(LM,1) = 1
            KGSVLM(LM,2) = 1
C
            DO IR = 1,IRTOP
C
               GX = GVLM(IR,LM,1)/SQRT_2
               GY = GVLM(IR,LM,2)/SQRT_2
               GSVLM(IR,LM,1) = DCMPLX(GX,+GY)
               GSVLM(IR,LM,2) = DCMPLX(GX,-GY)
C
               GX = GBLM(IR,LM,1)/SQRT_2
               GY = GBLM(IR,LM,2)/SQRT_2
               GSBLM(IR,LM,1) = DCMPLX(GX,+GY)
               GSBLM(IR,LM,2) = DCMPLX(GX,-GY)
C
            END DO
C
         END IF
C
         IF ( KGVLM(LM,3).NE.0 ) THEN
C
            KGSVLM(LM,3) = 1
C
            GSVLM(1:IRTOP,LM,3) = GVLM(1:IRTOP,LM,3)
            GSBLM(1:IRTOP,LM,3) = GBLM(1:IRTOP,LM,3)
C
         END IF
C
      END DO
      WRITE (*,'(A,20I3)') '###-------------------------------spherical'
      WRITE (*,'(A,20I3)') '### KGSVLM 1:',KGSVLM(1:9,1)
      WRITE (*,'(A,20I3)') '### KGSVLM 2:',KGSVLM(1:9,2)
      WRITE (*,'(A,20I3)') '### KGSVLM 3:',KGSVLM(1:9,3)
      IR = 200
      WRITE (*,'(A,20I3)') '###----------------------'
      WRITE (*,'(A,2e20.10)') 'GVLM  4 x',GVLM(IR,4,1)
      WRITE (*,'(A,2e20.10)') 'GVLM  2 y',GVLM(IR,2,2)
      WRITE (*,'(A,2e20.10)') 'GVLM  3 z',GVLM(IR,3,3)
      WRITE (*,'(A,20I3)') '###----------------------'
      WRITE (*,'(A,2e20.10)') 'GSVLM 2 +',GSVLM(IR,2,1)
      WRITE (*,'(A,2e20.10)') 'GSVLM 4 +',GSVLM(IR,4,1)
      WRITE (*,'(A,2e20.10)') 'GSVLM 2 -',GSVLM(IR,2,2)
      WRITE (*,'(A,2e20.10)') 'GSVLM 4 -',GSVLM(IR,4,2)
      WRITE (*,'(A,2e20.10)') 'GSVLM 3 z',GSVLM(IR,3,3)
C=======================================================================
C
      MZAZB(:,:,:) = C0
      MZBZA(:,:,:) = C0
      MIRR_2(:,:,:,:) = C0
      MIRR_3(:,:,:,:) = C0
      MIRR_4(:,:,:,:) = C0
C
C=======================================================================
C
C      evaluate matrix element integral functions of the type
C
C          IXAYB(r) = < X^+(E_a) | H_lam | Y(E_b) >_r
C
C=======================================================================
C
C ----------------------------------------- read in wavefunctions for LA
C
      CALL WAVFUN_READ_REL(IFIL_LHSA,IT,1,ZGLA,ZFLA,JGLA,JFLA,IRTOP,
     &                     NCPLWF_LA,IKMCPLWF_LA)
C
C------------- take complex conjugate of wave function ZGLA if requested
C
      IF ( ME_CC_BRA_RWF ) THEN
         DO IKM = 1,NKM
            ZGLA(:,:,IKM) = DCONJG(ZGLA(:,:,IKM))
            ZFLA(:,:,IKM) = DCONJG(ZFLA(:,:,IKM))
            JGLA(:,:,IKM) = DCONJG(JGLA(:,:,IKM))
            JFLA(:,:,IKM) = DCONJG(JFLA(:,:,IKM))
         END DO
      END IF
C
C ----------------------------------------- read in wavefunctions for RB
C
      CALL WAVFUN_READ_REL(IFIL_RHSB,IT,1,ZGRB,ZFRB,JGRB,JFRB,IRTOP,
     &                     NCPLWF_RB,IKMCPLWF_RB)
C
C--------- convolute the wave functions *RB with the potential functions
C
      CALL MEGRV_CONVOLUTE(IT,VLM,BLM,KVLM,GVLM,GBLM,KGVLM,ZGRB,ZFRB,
     &                     JGRB,JFRB,IKMCPLWF_RB,NCPLWF_RB,V_ZGRB,
     &                     V_ZFRB,V_JGRB,V_JFRB,B_ZGRB,B_ZFRB,B_JGRB,
     &                     B_JFRB,GV_ZGRB,GV_ZFRB,GV_JGRB,GV_JFRB,
     &                     GB_ZGRB,GB_ZFRB,GB_JGRB,GB_JFRB,KGV_ZGR)
C
C=======================================================================
C
C      evaluate the matrix elements of the type
C
C               MZBZA = < Z^+(E_b) | H_lam | Z(E_a) >
C
C=======================================================================
C
C ----------------------------------------- read in wavefunctions for LB
C
      CALL WAVFUN_READ_REL(IFIL_LHSB,IT,1,ZGLB,ZFLB,JGLB,JFLB,IRTOP,
     &                     NCPLWF_LB,IKMCPLWF_LB)
C
C------------- take complex conjugate of wave function ZGRA if requested
C
      IF ( ME_CC_BRA_RWF ) THEN
         DO IKM = 1,NKM
            ZGLB(:,:,IKM) = DCONJG(ZGLB(:,:,IKM))
            ZFLB(:,:,IKM) = DCONJG(ZFLB(:,:,IKM))
            JGLB(:,:,IKM) = DCONJG(JGLB(:,:,IKM))
            JFLB(:,:,IKM) = DCONJG(JFLB(:,:,IKM))
         END DO
      END IF
C
C ----------------------------------------- read in wavefunctions for RA
C
      CALL WAVFUN_READ_REL(IFIL_RHSA,IT,1,ZGRA,ZFRA,JGRA,JFRA,IRTOP,
     &                     NCPLWF_RA,IKMCPLWF_RA)
C
C--------- convolute the wave functions *RA with the potential functions
C
      CALL MEGRV_CONVOLUTE(IT,VLM,BLM,KVLM,GVLM,GBLM,KGVLM,ZGRA,ZFRA,
     &                     JGRA,JFRA,IKMCPLWF_RA,NCPLWF_RA,V_ZGRA,
     &                     V_ZFRA,V_JGRA,V_JFRA,B_ZGRA,B_ZFRA,B_JGRA,
     &                     B_JFRA,GV_ZGRA,GV_ZFRA,GV_JGRA,GV_JFRA,
     &                     GB_ZGRA,GB_ZFRA,GB_JGRA,GB_JFRA,KGV_ZGR)
C
C ------------------------------------------------------- calculate B4A1
C
      CALL MEGRV_B4A1(OMEGA_BA,PBA_GRV,TMC,TMCSQ,IM,IRTOP,ZGLB,ZFLB,
     &                ZGRA,ZFRA,VLM,BLM,KVLM,GSVLM,GSBLM,KGSVLM,V_ZGRA,
     &                V_ZFRA,B_ZGRA,B_ZFRA,MZBZA)
C
C=======================================================================
C
C      evaluate the matrix elements involving irregular solutions
C
C=======================================================================
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C                                                      energy B -- LAMB3
      DO LAMB3 = 1,NKM
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                      energy A -- LAMA2
         DO LAMA2 = 1,NKM
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP JPOL
            DO JPOL = 1,NPOL
C
               CALL MEGRV_A2B3(LAMA2,LAMB3,JPOL,OMEGA_AB,PAB_GRV,TMC,
     &                         TMCSQ,IM,IRTOP,ZFLA,ZGLA,JGLA,JFLA,ZGRB,
     &                         ZFRB,JFRB,JGRB,V_ZGRB,V_ZFRB,B_ZGRB,
     &                         B_ZFRB,V_JFRB,V_JGRB,B_JFRB,B_JGRB,VLM,
     &                         BLM,KVLM,GSVLM,GSBLM,KGSVLM,IZAZB,IZAJB,
     &                         IJAZB,IJAJB,MZAZB,K_SELECTION_RULES)
C
C#######################################################################
C
               IF ( K_SELECTION_RULES ) THEN
C
C***********************************************************************
C***********************************************************************
C     TERM 1    MZBZA(4,1) * TAU(1,2,a) * MZAZB(2,3) * TAU(3,4,b)
C***********************************************************************
C***********************************************************************
C
C***********************************************************************
C***********************************************************************
C     TERM 2    - TAU(3,4,b) * MIRR_2(3,4)                 LAMA1 = LAMA2
C***********************************************************************
C***********************************************************************
C
                  LAMA1 = LAMA2
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C                                                      energy B -- LAMB4
                  DO LAMB4 = 1,NKM
C
                     CALL MEGRV_BXAXW_REL(OMEGA_BA,PBA_GRV,TMC,TMCSQ,IM,
     &                  IRTOP,LAMB4,ZFLB,ZGLB,ZGLB,ZFLB,LAMA1,ZGRA,ZFRA,
     &                  JGRA,JFRA,V_ZGRA,V_ZFRA,B_ZGRA,B_ZFRA,GV_ZGRA,
     &                  GV_ZFRA,GB_ZGRA,GB_ZFRA,V_JGRA,V_JFRA,B_JGRA,
     &                  B_JFRA,GV_JGRA,GV_JFRA,GB_JGRA,GB_JFRA,VLM,BLM,
     &                  KVLM,GSVLM,GSBLM,KGSVLM,IZAZB,IJAZB,KGV_ZGR,
     &                  MIRR_2,LAMB4,LAMB3,JPOL)
C
                  END DO
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB LAMB4
C
C
C***********************************************************************
C***********************************************************************
C     TERM 3    - TAU(1,2,a) * MIRR_3(1,2)                 LAMB3 = LAMB4
C***********************************************************************
C***********************************************************************
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                      energy A -- LAMA1
                  DO LAMA1 = 1,NKM
C
                     CALL MEGRV_BXAXW_REL(OMEGA_BA,PBA_GRV,TMC,TMCSQ,IM,
     &                  IRTOP,LAMB3,ZFLB,ZGLB,JGLB,JFLB,LAMA1,ZGRA,ZFRA,
     &                  ZGRA,ZFRA,V_ZGRA,V_ZFRA,B_ZGRA,B_ZFRA,GV_ZGRA,
     &                  GV_ZFRA,GB_ZGRA,GB_ZFRA,V_ZGRA,V_ZFRA,B_ZGRA,
     &                  B_ZFRA,GV_ZGRA,GV_ZFRA,GB_ZGRA,GB_ZFRA,VLM,BLM,
     &                  KVLM,GSVLM,GSBLM,KGSVLM,IZAZB,IZAJB,KGV_ZGR,
     &                  MIRR_3,LAMA2,LAMA1,JPOL)
C
                  END DO
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA LAMA1
C
C
C***********************************************************************
C***********************************************************************
C     TERM 4    MIRR_4(1,3)           LAMA1 = LAMA2  and   LAMB3 = LAMB4
C***********************************************************************
C***********************************************************************
C
                  LAMA1 = LAMA2
C
                  CALL MEGRV_BXAXW_REL(OMEGA_BA,PBA_GRV,TMC,TMCSQ,IM,
     &                                 IRTOP,LAMB3,ZFLB,ZGLB,JGLB,JFLB,
     &                                 LAMA1,ZGRA,ZFRA,JGRA,JFRA,V_ZGRA,
     &                                 V_ZFRA,B_ZGRA,B_ZFRA,GV_ZGRA,
     &                                 GV_ZFRA,GB_ZGRA,GB_ZFRA,V_JGRA,
     &                                 V_JFRA,B_JGRA,B_JFRA,GV_JGRA,
     &                                 GV_JFRA,GB_JGRA,GB_JFRA,VLM,BLM,
     &                                 KVLM,GSVLM,GSBLM,KGSVLM,IZAZB,
     &                                 IJAJB,KGV_ZGR,MIRR_4,LAMA1,LAMB3,
     &                                 JPOL)
C
C***********************************************************************
C
               END IF
C################################################################ IFLAG1
            END DO
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP JPOL
C
         END DO
C                                                      energy A -- LAMA2
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
      END DO
C                                                      energy B -- LAMB3
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ( IT.LE.-1110 ) THEN
         IM = 200
         DO I = 1,18
            DO J = 1,18
               DO K = 1,3
                  IF ( ABS(MZBZA(I,J,K)).LT.1D-12 ) MZBZA(I,J,K) = 0D0
                  IF ( ABS(MZAZB(I,J,K)).LT.1D-12 ) MZAZB(I,J,K) = 0D0
               END DO
            END DO
         END DO
         WRITE (IM,*) 'MZAZB'
         WRITE (IM,'(2e22.14)') MZAZB
         WRITE (IM,*) 'MZBZA'
         WRITE (IM,'(2e22.14)') MZBZA
         WRITE (IM,*) 'MIRR_2'
         WRITE (IM,'(2e22.14)') MIRR_2
         WRITE (IM,*) 'MIRR_3'
         WRITE (IM,'(2e22.14)') MIRR_3
         WRITE (IM,*) 'MIRR_4'
         WRITE (IM,'(2e22.14)') MIRR_4
         STOP
      END IF
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
      END
C*==megrv_bxaxw_rel.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MEGRV_BXAXW_REL(OMEGA_BA,PBA_GRV,TMC,TMCSQ,IM,IRTOP,
     &                           LAMB,ZFLB,ZGLB,JGLB,JFLB,LAMA,ZGRA,
     &                           ZFRA,JGRA,JFRA,V_ZGRA,V_ZFRA,B_ZGRA,
     &                           B_ZFRA,GV_ZGRA,GV_ZFRA,GB_ZGRA,GB_ZFRA,
     &                           V_JGRA,V_JFRA,B_JGRA,B_JFRA,GV_JGRA,
     &                           GV_JFRA,GB_JGRA,GB_JFRA,VLM,BLM,KVLM,
     &                           GSVLM,GSBLM,KGSVLM,IZAZB,IJAJB,KGV_ZGR,
     &                           MIRR,I,J,JPOL)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:IKMCPLWF_LB,IKMCPLWF_RA,NCPLWFMAX,NLMFP
      USE MOD_RMESH,ONLY:JRWS,NRMAX,JRCRI,FULLPOT
      USE MOD_CONSTANTS,ONLY:C0,CI
      USE MOD_ANGMOM,ONLY:AB1_GRV,AB2_GRV,AC1_GRV,AC2_GRV,NKMMAX,NKM,
     &    IMKM_IKM,NCPLWF_LB,NCPLWF_RA,AME_G,AME_F,NPOL
      IMPLICIT NONE
C*--MEGRV_BXAXW_REL608
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='MEGRV_BXAXW_REL')
      REAL*8 NON0_TOL
      PARAMETER (NON0_TOL=1D-8)
C
C Dummy arguments
C
      INTEGER I,IM,IRTOP,J,JPOL,LAMA,LAMB
      COMPLEX*16 OMEGA_BA,PBA_GRV,TMC,TMCSQ
      REAL*8 BLM(NRMAX,NLMFP),VLM(NRMAX,NLMFP)
      COMPLEX*16 B_JFRA(NRMAX,NCPLWFMAX,NKM),B_JGRA(NRMAX,NCPLWFMAX,NKM)
     &           ,B_ZFRA(NRMAX,NCPLWFMAX,NKM),
     &           B_ZGRA(NRMAX,NCPLWFMAX,NKM),GB_JFRA(NRMAX,NKM,NKM,3),
     &           GB_JGRA(NRMAX,NKM,NKM,3),GB_ZFRA(NRMAX,NKM,NKM,3),
     &           GB_ZGRA(NRMAX,NKM,NKM,3),GSBLM(NRMAX,NLMFP,3),
     &           GSVLM(NRMAX,NLMFP,3),GV_JFRA(NRMAX,NKM,NKM,3),
     &           GV_JGRA(NRMAX,NKM,NKM,3),GV_ZFRA(NRMAX,NKM,NKM,3),
     &           GV_ZGRA(NRMAX,NKM,NKM,3),IJAJB(NRMAX),IZAZB(NRMAX),
     &           JFLB(NRMAX,NCPLWFMAX,NKM),JFRA(NRMAX,NCPLWFMAX,NKM),
     &           JGLB(NRMAX,NCPLWFMAX,NKM),JGRA(NRMAX,NCPLWFMAX,NKM),
     &           MIRR(NKMMAX,NKMMAX,3,3),V_JFRA(NRMAX,NCPLWFMAX,NKM),
     &           V_JGRA(NRMAX,NCPLWFMAX,NKM),V_ZFRA(NRMAX,NCPLWFMAX,NKM)
     &           ,V_ZGRA(NRMAX,NCPLWFMAX,NKM),ZFLB(NRMAX,NCPLWFMAX,NKM),
     &           ZFRA(NRMAX,NCPLWFMAX,NKM),ZGLB(NRMAX,NCPLWFMAX,NKM),
     &           ZGRA(NRMAX,NCPLWFMAX,NKM)
      INTEGER KGSVLM(NLMFP,3),KGV_ZGR(NKM,NKM,3),KVLM(NLMFP)
C
C Local variables
C
      COMPLEX*16 AAB1,AAB2,AABF,AABG,AAV1,AAV2,AAVF,AAVG,FBF,FBG,FV(:),
     &           FVF,FVG,GBF,GBG,GVF,GVG,PBA_OM_TMC,PBA_TMCSQ,VCNTR,
     &           WIRR,WREG
      INTEGER IA,IA_ERR,IB,IDOS,IFLAG2,IKMA,IKMB,IMKMA,IMKMB,IPOL,IR,
     &        ISPN
      REAL*8 RARG
      LOGICAL RNON0
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE FV
C
      RNON0(RARG) = ABS(RARG).GT.NON0_TOL
C
      ALLOCATE (FV(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate FV')
C
      IDOS = 1
      ISPN = 2
C
      PBA_TMCSQ = PBA_GRV*TMCSQ
      PBA_OM_TMC = PBA_GRV*OMEGA_BA*TMC
C
      IF ( FULLPOT ) THEN
         IRTOP = JRCRI(IM)
      ELSE
         IRTOP = JRWS(IM)
      END IF
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
      DO IPOL = 1,NPOL
C
         FV(:) = C0
C
         IFLAG2 = 0
C
C =================================================== ******************
C                                                     *  grad V  and B *
C                                                     ******************
         DO IKMA = 1,NKM
            IMKMA = IMKM_IKM(IKMA)
            IF ( KGV_ZGR(IKMA,LAMA,IPOL).NE.0 ) THEN
C
               DO IB = 1,NCPLWF_LB(LAMB)
                  IKMB = IKMCPLWF_LB(IB,LAMB)
                  IMKMB = IMKM_IKM(IKMB)
C
                  IF ( IKMB.EQ.IKMA .OR. RNON0(AME_G(IKMB,IKMA,2,ISPN))
     &                 ) THEN
C
                     IFLAG2 = 1
C
                     AAVG = PBA_TMCSQ*AME_G(IKMB,IKMA,1,IDOS)
                     AABG = PBA_TMCSQ*AME_G(IKMB,IKMA,2,ISPN)
C
                     AAVF = PBA_TMCSQ*AME_F(IKMB,IKMA,1,IDOS)
                     AABF = PBA_TMCSQ*AME_F(IKMB,IKMA,2,ISPN)
C
                     IF ( (IMKMB.GT.NKM) .OR. (IMKMA.GT.NKM) ) THEN
                        AAVF = 0.0D0
                        AABF = 0.0D0
                     END IF
C
                     DO IR = 1,IRTOP
                        GVG = AAVG*JGLB(IR,IB,LAMB)
     &                        *GV_JGRA(IR,IKMA,LAMA,IPOL)
                        GBG = AABG*JGLB(IR,IB,LAMB)
     &                        *GB_JGRA(IR,IKMA,LAMA,IPOL)
                        FVF = AAVF*JFLB(IR,IB,LAMB)
     &                        *GV_JFRA(IR,IKMA,LAMA,IPOL)
                        FBF = AABF*JFLB(IR,IB,LAMB)
     &                        *GB_JFRA(IR,IKMA,LAMA,IPOL)
C
                        WREG = GVG + GBG + FVF + FBF
C
                        GVG = AAVG*ZGLB(IR,IB,LAMB)
     &                        *GV_ZGRA(IR,IKMA,LAMA,IPOL)
                        GBG = AABG*ZGLB(IR,IB,LAMB)
     &                        *GB_ZGRA(IR,IKMA,LAMA,IPOL)
                        FVF = AAVF*ZFLB(IR,IB,LAMB)
     &                        *GV_ZFRA(IR,IKMA,LAMA,IPOL)
                        FBF = AABF*ZFLB(IR,IB,LAMB)
     &                        *GB_ZFRA(IR,IKMA,LAMA,IPOL)
C
                        WIRR = GVG + GBG + FVF + FBF
C
                        FV(IR) = FV(IR) + WREG*IZAZB(IR)
     &                           + WIRR*IJAJB(IR)
C
                     END DO
C
                  END IF
C
               END DO
            END IF
         END DO
C
C =================================================== ******************
C                                                     *  V/c alfa . a  *
C                                                     ******************
         DO IA = 1,NCPLWF_RA(LAMA)
            IKMA = IKMCPLWF_RA(IA,LAMA)
C
            DO IB = 1,NCPLWF_LB(LAMB)
               IKMB = IKMCPLWF_LB(IB,LAMB)
C
               IF ( RNON0(AB1_GRV(IKMB,IKMA,IPOL)) .OR. 
     &              RNON0(AB2_GRV(IKMB,IKMA,IPOL)) ) THEN
C
                  IFLAG2 = 1
C
                  AAV1 = PBA_OM_TMC*AB1_GRV(IKMB,IKMA,IPOL)
                  AAV2 = PBA_OM_TMC*AB2_GRV(IKMB,IKMA,IPOL)
C
                  DO IR = 1,IRTOP
                     GVF = AAV1*JGLB(IR,IB,LAMB)*V_JFRA(IR,IA,LAMA)
                     FVG = AAV2*JFLB(IR,IB,LAMB)*V_JGRA(IR,IA,LAMA)
                     WREG = GVF - FVG
C
                     GVF = AAV1*ZGLB(IR,IB,LAMB)*V_ZFRA(IR,IA,LAMA)
                     FVG = AAV2*ZFLB(IR,IB,LAMB)*V_ZGRA(IR,IA,LAMA)
                     WIRR = GVF - FVG
C
                     FV(IR) = FV(IR) + WREG*IZAZB(IR) + WIRR*IJAJB(IR)
C
                  END DO
C
               END IF
C
            END DO
         END DO
C
C =================================================== ******************
C                                                     * ibB/c alfa x a *
C                                                     ******************
         DO IA = 1,NCPLWF_RA(LAMA)
            IKMA = IKMCPLWF_RA(IA,LAMA)
C
            DO IB = 1,NCPLWF_LB(LAMB)
               IKMB = IKMCPLWF_LB(IB,LAMB)
C
               IF ( RNON0(AC1_GRV(IKMB,IKMA,IPOL)) .OR. 
     &              RNON0(AC2_GRV(IKMB,IKMA,IPOL)) ) THEN
C
                  IFLAG2 = 1
C
                  AAB1 = PBA_OM_TMC*(-CI)*CI*AC1_GRV(IKMB,IKMA,IPOL)
                  AAB2 = PBA_OM_TMC*(-CI)*CI*AC2_GRV(IKMB,IKMA,IPOL)
C
                  DO IR = 1,IRTOP
                     GBF = AAB1*JGLB(IR,IB,LAMB)*B_JFRA(IR,IA,LAMA)
                     FBG = AAB2*JFLB(IR,IB,LAMB)*B_JGRA(IR,IA,LAMA)
                     WREG = GBF + FBG
C
                     GBF = AAB1*ZGLB(IR,IB,LAMB)*B_ZFRA(IR,IA,LAMA)
                     FBG = AAB2*ZFLB(IR,IB,LAMB)*B_ZGRA(IR,IA,LAMA)
                     WIRR = GBF + FBG
C
                     FV(IR) = FV(IR) + WREG*IZAZB(IR) + WIRR*IJAJB(IR)
C
                  END DO
C
               END IF
C
            END DO
C
         END DO
C
C                                                                IA  IB
C-----------------------------------------------------------------------
C
C#######################################################################
C
         IF ( IFLAG2.EQ.1 ) THEN
C
            CALL CRADINT(IM,FV,VCNTR)
C
C***********************************************************************
C***********************************************************************
C     TERM 1    MZBZA(4,1) * TAU(1,2,a) * MZAZB(2,3) * TAU(3,4,b)
C***********************************************************************
C***********************************************************************
C
            MIRR(I,J,IPOL,JPOL) = MIRR(I,J,IPOL,JPOL) + VCNTR
C
         END IF
C################################################################ IFLAG2
      END DO
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
C
      END
C*==megrv_b4a1.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MEGRV_B4A1(OMEGA_BA,PBA_GRV,TMC,TMCSQ,IM,IRTOP,ZGLB,
     &                      ZFLB,ZGRA,ZFRA,VLM,BLM,KVLM,GSVLM,GSBLM,
     &                      KGSVLM,V_ZGRA,V_ZFRA,B_ZGRA,B_ZFRA,MZBZA)
C   ********************************************************************
C   *                                                                  *
C   *  evaluate the matrix elements of the type                        *
C   *                                                                  *
C   *           MZBZA = < Z^+(E_b) | H_lam | Z(E_a) >                  *
C   *                                                                  *
C   *     with    H_lam   the electric dipole operator                 *
C   *                                                                  *
C   *     in the  grad V  - form                                       *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *  NOTE: MZBZA assumed to be initialized                           *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:CI,C0
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_TYPES,ONLY:IKMCPLWF_LB,IKMCPLWF_RA,NCPLWFMAX,NLMFP
      USE MOD_ANGMOM,ONLY:AB1_GRV,AB2_GRV,AC1_GRV,AC2_GRV,NKMMAX,NKM,
     &    IMKM_IKM,NCPLWF_LB,NCPLWF_RA,AG_RGNT,A_SIGY,NPOL
      IMPLICIT NONE
C*--MEGRV_B4A1874
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='MEGRV_B4A1')
      REAL*8 NON0_TOL
      PARAMETER (NON0_TOL=1D-8)
C
C Dummy arguments
C
      INTEGER IM,IRTOP
      COMPLEX*16 OMEGA_BA,PBA_GRV,TMC,TMCSQ
      REAL*8 BLM(NRMAX,NLMFP),VLM(NRMAX,NLMFP)
      COMPLEX*16 B_ZFRA(NRMAX,NCPLWFMAX,NKM),B_ZGRA(NRMAX,NCPLWFMAX,NKM)
     &           ,GSBLM(NRMAX,NLMFP,3),GSVLM(NRMAX,NLMFP,3),
     &           MZBZA(NKMMAX,NKMMAX,3),V_ZFRA(NRMAX,NCPLWFMAX,NKM),
     &           V_ZGRA(NRMAX,NCPLWFMAX,NKM),ZFLB(NRMAX,NCPLWFMAX,NKM),
     &           ZFRA(NRMAX,NCPLWFMAX,NKM),ZGLB(NRMAX,NCPLWFMAX,NKM),
     &           ZGRA(NRMAX,NCPLWFMAX,NKM)
      INTEGER KGSVLM(NLMFP,3),KVLM(NLMFP)
C
C Local variables
C
      COMPLEX*16 AAB1,AAB2,AABF,AABG,AAV1,AAV2,AAVF,AAVG,CARG,FBF,FBG,
     &           FVF,FVG,FZBZA(:),GBF,GBG,GVF,GVG,PBA_OM_TMC,PBA_TMCSQ,
     &           V2BAR
      LOGICAL CNON0,RNON0
      INTEGER IA1,IA_ERR,IB4,IFLAG1,IKMA1,IKMB4,IMKMA1,IMKMB4,IPOL,IR,
     &        LAMA1,LAMB4,LM
      REAL*8 RARG
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE FZBZA
C
      RNON0(RARG) = ABS(RARG).GT.NON0_TOL
      CNON0(CARG) = ABS(DREAL(CARG)) + ABS(DIMAG(CARG)).GT.NON0_TOL
C
      ALLOCATE (FZBZA(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: FZBZA')
C
      PBA_TMCSQ = PBA_GRV*TMCSQ
      PBA_OM_TMC = PBA_GRV*OMEGA_BA*TMC
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C                                                      energy B -- LAMB4
      DO LAMB4 = 1,NKM
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                      energy A -- LAMA1
         DO LAMA1 = 1,NKM
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
            DO IPOL = 1,NPOL
C
               FZBZA(:) = C0
C
               IFLAG1 = 0
C
C =================================================== ******************
C                                                     *  grad V  and B *
C                                                     ******************
               DO IA1 = 1,NCPLWF_RA(LAMA1)
                  IKMA1 = IKMCPLWF_RA(IA1,LAMA1)
                  IMKMA1 = IMKM_IKM(IKMA1)
C
                  DO IB4 = 1,NCPLWF_LB(LAMB4)
                     IKMB4 = IKMCPLWF_LB(IB4,LAMB4)
                     IMKMB4 = IMKM_IKM(IKMB4)
C
                     DO LM = 1,NLMFP
C
                        IF ( KGSVLM(LM,IPOL).EQ.0 ) CYCLE
C
                        IF ( .NOT.(CNON0(AG_RGNT(IKMB4,IKMA1,LM)) .OR. 
     &                       CNON0(A_SIGY(IKMB4,IKMA1,3,LM))) ) CYCLE
C
                        IFLAG1 = 1
C
                        AAVG = PBA_TMCSQ*AG_RGNT(IKMB4,IKMA1,LM)
                        AABG = PBA_TMCSQ*A_SIGY(IKMB4,IKMA1,3,LM)
C
                        AAVF = PBA_TMCSQ*AG_RGNT(IMKMB4,IMKMA1,LM)
                        AABF = PBA_TMCSQ*A_SIGY(IMKMB4,IMKMA1,3,LM)
C
                        IF ( (IMKMB4.GT.NKM) .OR. (IMKMA1.GT.NKM) ) THEN
                           AAVF = 0.0D0
                           AABF = 0.0D0
                        END IF
C
                        DO IR = 1,IRTOP
                           GVG = AAVG*ZGLB(IR,IB4,LAMB4)
     &                           *ZGRA(IR,IA1,LAMA1)*GSVLM(IR,LM,IPOL)
                           GBG = AABG*ZGLB(IR,IB4,LAMB4)
     &                           *ZGRA(IR,IA1,LAMA1)*GSBLM(IR,LM,IPOL)
                           FVF = AAVF*ZFLB(IR,IB4,LAMB4)
     &                           *ZFRA(IR,IA1,LAMA1)*GSVLM(IR,LM,IPOL)
                           FBF = AABF*ZFLB(IR,IB4,LAMB4)
     &                           *ZFRA(IR,IA1,LAMA1)*GSBLM(IR,LM,IPOL)
C
                           FZBZA(IR) = FZBZA(IR) + GVG + GBG + FVF + FBF
                        END DO
C
                     END DO
                  END DO
               END DO
C
C =================================================== ******************
C                                                     *  V/c alfa . a  *
C                                                     ******************
               DO IA1 = 1,NCPLWF_RA(LAMA1)
                  IKMA1 = IKMCPLWF_RA(IA1,LAMA1)
C
                  DO IB4 = 1,NCPLWF_LB(LAMB4)
                     IKMB4 = IKMCPLWF_LB(IB4,LAMB4)
C
                     IF ( RNON0(AB1_GRV(IKMB4,IKMA1,IPOL)) .OR. 
     &                    RNON0(AB2_GRV(IKMB4,IKMA1,IPOL)) ) THEN
C
                        IFLAG1 = 1
C
                        AAV1 = PBA_OM_TMC*AB1_GRV(IKMB4,IKMA1,IPOL)
                        AAV2 = PBA_OM_TMC*AB2_GRV(IKMB4,IKMA1,IPOL)
C
                        DO IR = 1,IRTOP
                           GVF = AAV1*ZGLB(IR,IB4,LAMB4)
     &                           *V_ZFRA(IR,IA1,LAMA1)
                           FVG = AAV2*ZFLB(IR,IB4,LAMB4)
     &                           *V_ZGRA(IR,IA1,LAMA1)
                           FZBZA(IR) = FZBZA(IR) + GVF - FVG
                        END DO
C
                     END IF
C
                  END DO
               END DO
C
C =================================================== ******************
C                                                     * ibB/c alfa x a *
C                                                     ******************
               DO IA1 = 1,NCPLWF_RA(LAMA1)
                  IKMA1 = IKMCPLWF_RA(IA1,LAMA1)
C
                  DO IB4 = 1,NCPLWF_LB(LAMB4)
                     IKMB4 = IKMCPLWF_LB(IB4,LAMB4)
C
                     IF ( RNON0(AC1_GRV(IKMB4,IKMA1,IPOL)) .OR. 
     &                    RNON0(AC2_GRV(IKMB4,IKMA1,IPOL)) ) THEN
C
                        IFLAG1 = 1
C
                        AAB1 = PBA_OM_TMC*(-CI)
     &                         *CI*AC1_GRV(IKMB4,IKMA1,IPOL)
                        AAB2 = PBA_OM_TMC*(-CI)
     &                         *CI*AC2_GRV(IKMB4,IKMA1,IPOL)
C
                        DO IR = 1,IRTOP
                           GBF = AAB1*ZGLB(IR,IB4,LAMB4)
     &                           *B_ZFRA(IR,IA1,LAMA1)
                           FBG = AAB2*ZFLB(IR,IB4,LAMB4)
     &                           *B_ZGRA(IR,IA1,LAMA1)
                           FZBZA(IR) = FZBZA(IR) + GBF + FBG
                        END DO
C
                     END IF
C
                  END DO
C
               END DO
C                                                                IA1 IB4
C-----------------------------------------------------------------------
C
C#######################################################################
C
               IF ( IFLAG1.EQ.1 ) THEN
C
                  CALL CRADINT(IM,FZBZA,V2BAR)
C
C***********************************************************************
C***********************************************************************
C     TERM 1    MZBZA(4,1) * TAU(1,2,a) * MZAZB(2,3) * TAU(3,4,b)
C***********************************************************************
C***********************************************************************
C
                  MZBZA(LAMB4,LAMA1,IPOL) = MZBZA(LAMB4,LAMA1,IPOL)
     &               + V2BAR
C
               END IF
C################################################################ IFLAG1
            END DO
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
C
         END DO
C                                                      energy A -- LAMA1
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
      END DO
C                                                      energy B -- LAMB4
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
      END
C*==megrv_a2b3.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MEGRV_A2B3(LAMA2,LAMB3,JPOL,OMEGA_AB,PAB_GRV,TMC,TMCSQ,
     &                      IM,IRTOP,ZFLA,ZGLA,JGLA,JFLA,ZGRB,ZFRB,JFRB,
     &                      JGRB,V_ZGRB,V_ZFRB,B_ZGRB,B_ZFRB,V_JFRB,
     &                      V_JGRB,B_JFRB,B_JGRB,VLM,BLM,KVLM,GSVLM,
     &                      GSBLM,KGSVLM,IZAZB,IZAJB,IJAZB,IJAJB,MZAZB,
     &                      K_SELECTION_RULES)
C   ********************************************************************
C   *                                                                  *
C   *  evaluate matrix element integral functions of the type          *
C   *                                                                  *
C   *      IXAYB(r) = < X^+(E_a) | H_lam | Y(E_b) >_r                  *
C   *                                                                  *
C   *      SXAYB : corresponding surface contributions                 *
C   *                                                                  *
C   *  X, Y: = Z (J) regular (irregular) solutions to Dirac equation   *
C   *                                                                  *
C   *     with    H_lam   the electric dipole operator in the form     *
C   *                                                                  *
C   *             H_lam = mc alfa_lam = mc alfa * A_lam                *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:IKMCPLWF_LA,IKMCPLWF_RB,NCPLWFMAX,NLMFP
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_CONSTANTS,ONLY:C0,CI
      USE MOD_ANGMOM,ONLY:AB1_GRV,AB2_GRV,AC1_GRV,AC2_GRV,NKMMAX,NKM,
     &    IMKM_IKM,NCPLWF_LA,NCPLWF_RB,AG_RGNT,A_SIGY
      IMPLICIT NONE
C*--MEGRV_A2B31121
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='MEGRV_A2B3')
      REAL*8 NON0_TOL
      PARAMETER (NON0_TOL=1D-8)
C
C Dummy arguments
C
      INTEGER IM,IRTOP,JPOL,LAMA2,LAMB3
      LOGICAL K_SELECTION_RULES
      COMPLEX*16 OMEGA_AB,PAB_GRV,TMC,TMCSQ
      REAL*8 BLM(NRMAX,NLMFP),VLM(NRMAX,NLMFP)
      COMPLEX*16 B_JFRB(NRMAX,NCPLWFMAX,NKM),B_JGRB(NRMAX,NCPLWFMAX,NKM)
     &           ,B_ZFRB(NRMAX,NCPLWFMAX,NKM),
     &           B_ZGRB(NRMAX,NCPLWFMAX,NKM),GSBLM(NRMAX,NLMFP,3),
     &           GSVLM(NRMAX,NLMFP,3),IJAJB(NRMAX),IJAZB(NRMAX),
     &           IZAJB(NRMAX),IZAZB(NRMAX),JFLA(NRMAX,NCPLWFMAX,NKM),
     &           JFRB(NRMAX,NCPLWFMAX,NKM),JGLA(NRMAX,NCPLWFMAX,NKM),
     &           JGRB(NRMAX,NCPLWFMAX,NKM),MZAZB(NKMMAX,NKMMAX,3),
     &           V_JFRB(NRMAX,NCPLWFMAX,NKM),V_JGRB(NRMAX,NCPLWFMAX,NKM)
     &           ,V_ZFRB(NRMAX,NCPLWFMAX,NKM),
     &           V_ZGRB(NRMAX,NCPLWFMAX,NKM),ZFLA(NRMAX,NCPLWFMAX,NKM),
     &           ZFRB(NRMAX,NCPLWFMAX,NKM),ZGLA(NRMAX,NCPLWFMAX,NKM),
     &           ZGRB(NRMAX,NCPLWFMAX,NKM)
      INTEGER KGSVLM(NLMFP,3),KVLM(NLMFP)
C
C Local variables
C
      COMPLEX*16 AAB1,AAB2,AABF,AABG,AAV1,AAV2,AAVF,AAVG,CARG,FBF,FBG,
     &           FJAJB(:),FJAZB(:),FVF,FVG,FZAJB(:),FZAZB(:),GBF,GBG,
     &           GVF,GVG,PAB_OM_TMC,PAB_TMCSQ
      LOGICAL CNON0,RNON0
      INTEGER IA2,IA_ERR,IB3,IFLAG1,IKMA2,IKMB3,IMKMA2,IMKMB3,IR,LM
      REAL*8 RARG
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE FZAZB,FZAJB,FJAZB,FJAJB
C
      RNON0(RARG) = ABS(RARG).GT.NON0_TOL
      CNON0(CARG) = ABS(DREAL(CARG)) + ABS(DIMAG(CARG)).GT.NON0_TOL
C
      ALLOCATE (FZAZB(NRMAX),FZAJB(NRMAX))
      ALLOCATE (FJAZB(NRMAX),FJAJB(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate FZAZB')
C
      PAB_TMCSQ = PAB_GRV*TMCSQ
      PAB_OM_TMC = PAB_GRV*OMEGA_AB*TMC
C
      K_SELECTION_RULES = .FALSE.
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C                                                      energy B -- LAMB3
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                      energy A -- LAMA2
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP JPOL
C
      FZAZB(:) = C0
      FZAJB(:) = C0
      FJAZB(:) = C0
      FJAJB(:) = C0
C
      IFLAG1 = 0
C
C =================================================== ******************
C                                                     *  grad V  and B *
C                                                     ******************
      DO IB3 = 1,NCPLWF_RB(LAMB3)
         IKMB3 = IKMCPLWF_RB(IB3,LAMB3)
         IMKMB3 = IMKM_IKM(IKMB3)
C
         DO IA2 = 1,NCPLWF_LA(LAMA2)
            IKMA2 = IKMCPLWF_LA(IA2,LAMA2)
            IMKMA2 = IMKM_IKM(IKMA2)
C
            DO LM = 1,NLMFP
C
               IF ( KGSVLM(LM,JPOL).EQ.0 ) CYCLE
C
               IF ( .NOT.(CNON0(AG_RGNT(IKMA2,IKMB3,LM)) .OR. 
     &              CNON0(A_SIGY(IKMA2,IKMB3,3,LM))) ) CYCLE
C
               IFLAG1 = 1
C
               AAVG = PAB_TMCSQ*AG_RGNT(IKMA2,IKMB3,LM)
               AABG = PAB_TMCSQ*A_SIGY(IKMA2,IKMB3,3,LM)
C
               AAVF = PAB_TMCSQ*AG_RGNT(IMKMA2,IMKMB3,LM)
               AABF = PAB_TMCSQ*A_SIGY(IMKMA2,IMKMB3,3,LM)
C
               IF ( (IMKMA2.GT.NKM) .OR. (IMKMB3.GT.NKM) ) THEN
                  AAVF = 0.0D0
                  AABF = 0.0D0
               END IF
C
               DO IR = 1,IRTOP
                  GVG = AAVG*ZGLA(IR,IA2,LAMA2)*ZGRB(IR,IB3,LAMB3)
     &                  *GSVLM(IR,LM,JPOL)
                  GBG = AABG*ZGLA(IR,IA2,LAMA2)*ZGRB(IR,IB3,LAMB3)
     &                  *GSBLM(IR,LM,JPOL)
                  FVF = AAVF*ZFLA(IR,IA2,LAMA2)*ZFRB(IR,IB3,LAMB3)
     &                  *GSVLM(IR,LM,JPOL)
                  FBF = AABF*ZFLA(IR,IA2,LAMA2)*ZFRB(IR,IB3,LAMB3)
     &                  *GSBLM(IR,LM,JPOL)
C
                  FZAZB(IR) = FZAZB(IR) + GVG + GBG + FVF + FBF
C
                  GVG = AAVG*ZGLA(IR,IA2,LAMA2)*JGRB(IR,IB3,LAMB3)
     &                  *GSVLM(IR,LM,JPOL)
                  GBG = AABG*ZGLA(IR,IA2,LAMA2)*JGRB(IR,IB3,LAMB3)
     &                  *GSBLM(IR,LM,JPOL)
                  FVF = AAVF*ZFLA(IR,IA2,LAMA2)*JFRB(IR,IB3,LAMB3)
     &                  *GSVLM(IR,LM,JPOL)
                  FBF = AABF*ZFLA(IR,IA2,LAMA2)*JFRB(IR,IB3,LAMB3)
     &                  *GSBLM(IR,LM,JPOL)
C
                  FZAJB(IR) = FZAJB(IR) + GVG + GBG + FVF + FBF
C
                  GVG = AAVG*JGLA(IR,IA2,LAMA2)*ZGRB(IR,IB3,LAMB3)
     &                  *GSVLM(IR,LM,JPOL)
                  GBG = AABG*JGLA(IR,IA2,LAMA2)*ZGRB(IR,IB3,LAMB3)
     &                  *GSBLM(IR,LM,JPOL)
                  FVF = AAVF*JFLA(IR,IA2,LAMA2)*ZFRB(IR,IB3,LAMB3)
     &                  *GSVLM(IR,LM,JPOL)
                  FBF = AABF*JFLA(IR,IA2,LAMA2)*ZFRB(IR,IB3,LAMB3)
     &                  *GSBLM(IR,LM,JPOL)
C
                  FJAZB(IR) = FJAZB(IR) + GVG + GBG + FVF + FBF
C
                  GVG = AAVG*JGLA(IR,IA2,LAMA2)*JGRB(IR,IB3,LAMB3)
     &                  *GSVLM(IR,LM,JPOL)
                  GBG = AABG*JGLA(IR,IA2,LAMA2)*JGRB(IR,IB3,LAMB3)
     &                  *GSBLM(IR,LM,JPOL)
                  FVF = AAVF*JFLA(IR,IA2,LAMA2)*JFRB(IR,IB3,LAMB3)
     &                  *GSVLM(IR,LM,JPOL)
                  FBF = AABF*JFLA(IR,IA2,LAMA2)*JFRB(IR,IB3,LAMB3)
     &                  *GSBLM(IR,LM,JPOL)
C
                  FJAJB(IR) = FJAJB(IR) + GVG + GBG + FVF + FBF
C
               END DO
C
            END DO
         END DO
      END DO
C
C =================================================== ******************
C                                                     *  V/c alfa . a  *
C                                                     ******************
      DO IB3 = 1,NCPLWF_RB(LAMB3)
         IKMB3 = IKMCPLWF_RB(IB3,LAMB3)
C
         DO IA2 = 1,NCPLWF_LA(LAMA2)
            IKMA2 = IKMCPLWF_LA(IA2,LAMA2)
C
            IF ( RNON0(AB1_GRV(IKMA2,IKMB3,JPOL)) .OR. 
     &           RNON0(AB2_GRV(IKMA2,IKMB3,JPOL)) ) THEN
C
               IFLAG1 = 1
C
               AAV1 = PAB_OM_TMC*AB1_GRV(IKMA2,IKMB3,JPOL)
               AAV2 = PAB_OM_TMC*AB2_GRV(IKMA2,IKMB3,JPOL)
C
               DO IR = 1,IRTOP
                  GVF = AAV1*ZGLA(IR,IA2,LAMA2)*V_ZFRB(IR,IB3,LAMB3)
                  FVG = AAV2*ZFLA(IR,IA2,LAMA2)*V_ZGRB(IR,IB3,LAMB3)
C
                  FZAZB(IR) = FZAZB(IR) + GVF - FVG
C
                  GVF = AAV1*ZGLA(IR,IA2,LAMA2)*V_JFRB(IR,IB3,LAMB3)
                  FVG = AAV2*ZFLA(IR,IA2,LAMA2)*V_JGRB(IR,IB3,LAMB3)
C
                  FZAJB(IR) = FZAJB(IR) + GVF - FVG
C
                  GVF = AAV1*JGLA(IR,IA2,LAMA2)*V_ZFRB(IR,IB3,LAMB3)
                  FVG = AAV2*JFLA(IR,IA2,LAMA2)*V_ZGRB(IR,IB3,LAMB3)
C
                  FJAZB(IR) = FJAZB(IR) + GVF - FVG
C
                  GVF = AAV1*JGLA(IR,IA2,LAMA2)*V_JFRB(IR,IB3,LAMB3)
                  FVG = AAV2*JFLA(IR,IA2,LAMA2)*V_JGRB(IR,IB3,LAMB3)
C
                  FJAJB(IR) = FJAJB(IR) + GVF - FVG
C
               END DO
C
            END IF
C
         END DO
      END DO
C
C =================================================== ******************
C                                                     * ibB/c alfa x a *
C                                                     ******************
      DO IB3 = 1,NCPLWF_RB(LAMB3)
         IKMB3 = IKMCPLWF_RB(IB3,LAMB3)
C
         DO IA2 = 1,NCPLWF_LA(LAMA2)
            IKMA2 = IKMCPLWF_LA(IA2,LAMA2)
C
            IF ( RNON0(AC1_GRV(IKMA2,IKMB3,JPOL)) .OR. 
     &           RNON0(AC2_GRV(IKMA2,IKMB3,JPOL)) ) THEN
C
               IFLAG1 = 1
C
               AAB1 = PAB_OM_TMC*(-CI)*CI*AC1_GRV(IKMA2,IKMB3,JPOL)
               AAB2 = PAB_OM_TMC*(-CI)*CI*AC2_GRV(IKMA2,IKMB3,JPOL)
C
               DO IR = 1,IRTOP
                  GBF = AAB1*ZGLA(IR,IA2,LAMA2)*B_ZFRB(IR,IB3,LAMB3)
                  FBG = AAB2*ZFLA(IR,IA2,LAMA2)*B_ZGRB(IR,IB3,LAMB3)
C
                  FZAZB(IR) = FZAZB(IR) + GBF + FBG
C
                  GBF = AAB1*ZGLA(IR,IA2,LAMA2)*B_JFRB(IR,IB3,LAMB3)
                  FBG = AAB2*ZFLA(IR,IA2,LAMA2)*B_JGRB(IR,IB3,LAMB3)
C
                  FZAJB(IR) = FZAJB(IR) + GBF + FBG
C
                  GBF = AAB1*JGLA(IR,IA2,LAMA2)*B_ZFRB(IR,IB3,LAMB3)
                  FBG = AAB2*JFLA(IR,IA2,LAMA2)*B_ZGRB(IR,IB3,LAMB3)
C
                  FJAZB(IR) = FJAZB(IR) + GBF + FBG
C
                  GBF = AAB1*JGLA(IR,IA2,LAMA2)*B_JFRB(IR,IB3,LAMB3)
                  FBG = AAB2*JFLA(IR,IA2,LAMA2)*B_JGRB(IR,IB3,LAMB3)
C
                  FJAJB(IR) = FJAJB(IR) + GBF + FBG
C
               END DO
C
            END IF
C
         END DO
C
      END DO
C
C                                                                IB3 IA2
C-----------------------------------------------------------------------
C
C#######################################################################
C
      IF ( IFLAG1.EQ.1 ) THEN
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 2-3: INW
         CALL CRADINT_R(IM,FZAZB,IZAZB)
C
         CALL CRADINT_R(IM,FZAJB,IZAJB)
C
         CALL CRADINT_R(IM,FJAZB,IJAZB)
C
         CALL CRADINT_R(IM,FJAJB,IJAJB)
C
C***********************************************************************
C***********************************************************************
C     TERM 1    MZAZB(4,1) * TAU(1,2,a) * MZAZB(2,3) * TAU(3,4,b)
C***********************************************************************
C***********************************************************************
C
         MZAZB(LAMA2,LAMB3,JPOL) = MZAZB(LAMA2,LAMB3,JPOL)
     &                             + IZAZB(IRTOP)
C
         K_SELECTION_RULES = .TRUE.
C
      END IF
C################################################################ IFLAG1
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP JPOL
C                                                      energy B -- LAMB3
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                      energy A -- LAMA2
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
      END
C*==megrv_convolute.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MEGRV_CONVOLUTE(IT,VLM,BLM,KVLM,GVLM,GBLM,KGVLM,ZGR,
     &                           ZFR,JGR,JFR,IKMCPLWF_R,NCPLWF_R,V_ZGR,
     &                           V_ZFR,V_JGR,V_JFR,B_ZGR,B_ZFR,B_JGR,
     &                           B_JFR,GV_ZGR,GV_ZFR,GV_JGR,GV_JFR,
     &                           GB_ZGR,GB_ZFR,GB_JGR,GB_JFR,KGV_ZGR)
C   ********************************************************************
C   *                                                                  *
C   *  convolute the RHS wave functions ZGR,ZFR and JGR,JFR            *
C   *                                                                  *
C   *          XZ_L2[f]  = Sum{L1,L} X_L Z_L1 C(L1,L,L2)               *
C   *                                                                  *
C   *  X_L(r): the potential functions V_L(r) and B_L(r)               *
C   *          the convolution must not change the coupling scheme     *
C   *          i.e. the settings for NCPLWF and IKMCPLWF are unchanged *
C   *                                                                  *
C   *  X_L(r): the gradient of V_L(r) and B_L(r)                       *
C   *          this will lead to a new coupling scheme                 *
C   *          represented by the flag KGV_ZGR                         *
C   *                                1                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,AG_RGNT,AF_RGNT,WKM1
      USE MOD_RMESH,ONLY:NRMAX,JRCRI,FULLPOT,JRWS
      USE MOD_CONSTANTS,ONLY:CI,C0,SQRT_2
      USE MOD_TYPES,ONLY:IMT,NCPLWFMAX,NLMFP
      IMPLICIT NONE
C*--MEGRV_CONVOLUTE1442
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='MEGRV_CONVOLUTE')
C
C Dummy arguments
C
      INTEGER IT
      REAL*8 BLM(NRMAX,NLMFP),GBLM(NRMAX,NLMFP,3),GVLM(NRMAX,NLMFP,3),
     &       VLM(NRMAX,NLMFP)
      COMPLEX*16 B_JFR(NRMAX,NCPLWFMAX,NKM),B_JGR(NRMAX,NCPLWFMAX,NKM),
     &           B_ZFR(NRMAX,NCPLWFMAX,NKM),B_ZGR(NRMAX,NCPLWFMAX,NKM),
     &           GB_JFR(NRMAX,NKM,NKM,3),GB_JGR(NRMAX,NKM,NKM,3),
     &           GB_ZFR(NRMAX,NKM,NKM,3),GB_ZGR(NRMAX,NKM,NKM,3),
     &           GV_JFR(NRMAX,NKM,NKM,3),GV_JGR(NRMAX,NKM,NKM,3),
     &           GV_ZFR(NRMAX,NKM,NKM,3),GV_ZGR(NRMAX,NKM,NKM,3),
     &           JFR(NRMAX,NCPLWFMAX,NKMMAX),JGR(NRMAX,NCPLWFMAX,NKMMAX)
     &           ,V_JFR(NRMAX,NCPLWFMAX,NKM),V_JGR(NRMAX,NCPLWFMAX,NKM),
     &           V_ZFR(NRMAX,NCPLWFMAX,NKM),V_ZGR(NRMAX,NCPLWFMAX,NKM),
     &           ZFR(NRMAX,NCPLWFMAX,NKMMAX),ZGR(NRMAX,NCPLWFMAX,NKMMAX)
      INTEGER IKMCPLWF_R(NCPLWFMAX,NKMMAX),KGVLM(NLMFP,3),
     &        KGV_ZGR(NKM,NKM,3),KVLM(NLMFP),NCPLWF_R(NKMMAX)
C
C Local variables
C
      COMPLEX*16 AME_F,AME_G,A_B,A_GB,A_GV,A_V,GSBLM(:,:),GSVLM(:,:)
      INTEGER I1,I2,IKM1,IKM2,IM,IPOL,IR,IRTOP,KGSVLM(:),LAM,LM
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE GSBLM,GSVLM,KGSVLM
C
      ALLOCATE (GSBLM(NRMAX,3),GSVLM(NRMAX,3))
      ALLOCATE (KGSVLM(3))
C
      IM = IMT(IT)
      IF ( FULLPOT ) THEN
         IRTOP = JRCRI(IM)
      ELSE
         IRTOP = JRWS(IM)
      END IF
C
C=======================================================================
C
      V_ZGR(:,:,:) = C0
      V_ZFR(:,:,:) = C0
      V_JGR(:,:,:) = C0
      V_JFR(:,:,:) = C0
      B_ZGR(:,:,:) = C0
      B_ZFR(:,:,:) = C0
      B_JGR(:,:,:) = C0
      B_JFR(:,:,:) = C0
C
C=======================================================================
C               convolute wave functions with V and B
C             this should not change the coupling scheme
C=======================================================================
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                                    LAM
      DO LAM = 1,NKM
C
C-------------------------------------------------------------------- I2
         DO I2 = 1,NCPLWF_R(LAM)
            IKM2 = IKMCPLWF_R(I2,LAM)
C
C-------------------------------------------------------------------- I1
            DO I1 = 1,NCPLWF_R(LAM)
               IKM1 = IKMCPLWF_R(I1,LAM)
C
C-------------------------------------------------------------------- LM
               DO LM = 1,NLMFP
C
                  AME_G = AG_RGNT(IKM1,IKM2,LM)
                  AME_F = AF_RGNT(IKM1,IKM2,LM)
C
                  IF ( ABS(AME_G-AME_F).GT.1D-8 )
     &                 CALL STOP_MESSAGE(ROUTINE,'AME_G <> AME_F')
C
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
                  IF ( ABS(AME_G).GT.1D-6 .AND. KVLM(LM).NE.0 ) THEN
C
                     DO IR = 1,IRTOP
C
                        A_V = AME_G*VLM(IR,LM)
C
                        V_ZGR(IR,I2,LAM) = V_ZGR(IR,I2,LAM)
     &                     + ZGR(IR,I1,LAM)*A_V
C
                        V_ZFR(IR,I2,LAM) = V_ZFR(IR,I2,LAM)
     &                     + ZFR(IR,I1,LAM)*A_V
C
                        V_JGR(IR,I2,LAM) = V_JGR(IR,I2,LAM)
     &                     + JGR(IR,I1,LAM)*A_V
C
                        V_JFR(IR,I2,LAM) = V_JFR(IR,I2,LAM)
     &                     + JFR(IR,I1,LAM)*A_V
C
                        A_B = AME_G*BLM(IR,LM)
C
                        B_ZGR(IR,I2,LAM) = B_ZGR(IR,I2,LAM)
     &                     + ZGR(IR,I1,LAM)*A_B
C
                        B_ZFR(IR,I2,LAM) = B_ZFR(IR,I2,LAM)
     &                     + ZFR(IR,I1,LAM)*A_B
C
                        B_JGR(IR,I2,LAM) = B_JGR(IR,I2,LAM)
     &                     + JGR(IR,I1,LAM)*A_B
C
                        B_JFR(IR,I2,LAM) = B_JFR(IR,I2,LAM)
     &                     + JFR(IR,I1,LAM)*A_B
C
                     END DO
C
                  END IF
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
               END DO
C-------------------------------------------------------------------- LM
            END DO
C-------------------------------------------------------------------- I1
         END DO
C-------------------------------------------------------------------- I2
      END DO
C                                                                    LAM
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
C
C=======================================================================
C            convolute wave functions with gradient of V and B
C                this WILL change the coupling scheme
C=======================================================================
C
      GV_ZGR(:,:,:,:) = C0
      GV_ZFR(:,:,:,:) = C0
      GV_JGR(:,:,:,:) = C0
      GV_JFR(:,:,:,:) = C0
      GB_ZGR(:,:,:,:) = C0
      GB_ZFR(:,:,:,:) = C0
      GB_JGR(:,:,:,:) = C0
      GB_JFR(:,:,:,:) = C0
C
      KGV_ZGR(:,:,:) = 0
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                                    LAM
      DO LAM = 1,NKM
C
C-------------------------------------------------------------------- I2
         DO IKM2 = 1,NKM
C
C-------------------------------------------------------------------- I1
            DO I1 = 1,NCPLWF_R(LAM)
               IKM1 = IKMCPLWF_R(I1,LAM)
C
C-------------------------------------------------------------------- LM
               DO LM = 1,NLMFP
C
C=======================================================================
C        transfer gradients from CARTESIAN to SPHERICAL coordinates
C=======================================================================
C???????????????????????????????????????????????????????????????????????
C                index 3:  ipol= 1,2,3  ==  (+),(-),(z)
C                 M(X) =   [  M(+) + M(-) ] / SQRT(2)
C                 M(Y) = I*[ -M(+) + M(-) ] / SQRT(2)
C
                  KGSVLM(:) = 0
                  GSVLM(:,:) = C0
                  GSBLM(:,:) = C0
C
                  IF ( KGVLM(LM,1).NE.0 .OR. KGVLM(LM,2).NE.0 ) THEN
C
                     KGSVLM(1) = 1
                     KGSVLM(2) = 1
C
                     DO IR = 1,IRTOP
C
                        GSVLM(IR,1) = (GVLM(IR,LM,1)-CI*GVLM(IR,LM,2))
     &                                /SQRT_2
                        GSVLM(IR,2) = (GVLM(IR,LM,1)+CI*GVLM(IR,LM,2))
     &                                /SQRT_2
C
                        GSBLM(IR,1) = (GBLM(IR,LM,1)-CI*GBLM(IR,LM,2))
     &                                /SQRT_2
                        GSBLM(IR,2) = (GBLM(IR,LM,1)+CI*GBLM(IR,LM,2))
     &                                /SQRT_2
C
                     END DO
C
                  END IF
C
                  IF ( KGVLM(LM,3).NE.0 ) THEN
C
                     KGSVLM(3) = 1
C
                     GSVLM(1:IRTOP,3) = GVLM(1:IRTOP,LM,3)
                     GSBLM(1:IRTOP,3) = GBLM(1:IRTOP,LM,3)
C
                  END IF
C
                  AME_G = AG_RGNT(IKM1,IKM2,LM)
C
                  DO IPOL = 1,3
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
                     IF ( ABS(AME_G).GT.1D-6 .AND. KGSVLM(IPOL).NE.0 )
     &                    THEN
C
                        KGV_ZGR(IKM2,LAM,IPOL) = 1
C
                        DO IR = 1,IRTOP
C
                           A_GV = AME_G*GSVLM(IR,IPOL)
C
                           GV_ZGR(IR,IKM2,LAM,IPOL)
     &                        = GV_ZGR(IR,IKM2,LAM,IPOL)
     &                        + ZGR(IR,I1,LAM)*A_GV
C
                           GV_ZFR(IR,IKM2,LAM,IPOL)
     &                        = GV_ZFR(IR,IKM2,LAM,IPOL)
     &                        + ZFR(IR,I1,LAM)*A_GV
C
                           GV_JGR(IR,IKM2,LAM,IPOL)
     &                        = GV_JGR(IR,IKM2,LAM,IPOL)
     &                        + JGR(IR,I1,LAM)*A_GV
C
                           GV_JFR(IR,IKM2,LAM,IPOL)
     &                        = GV_JFR(IR,IKM2,LAM,IPOL)
     &                        + JFR(IR,I1,LAM)*A_GV
C
                           A_GB = AME_G*GSBLM(IR,IPOL)
C
                           GB_ZGR(IR,IKM2,LAM,IPOL)
     &                        = GB_ZGR(IR,IKM2,LAM,IPOL)
     &                        + ZGR(IR,I1,LAM)*A_GB
C
                           GB_ZFR(IR,IKM2,LAM,IPOL)
     &                        = GB_ZFR(IR,IKM2,LAM,IPOL)
     &                        + ZFR(IR,I1,LAM)*A_GB
C
                           GB_JGR(IR,IKM2,LAM,IPOL)
     &                        = GB_JGR(IR,IKM2,LAM,IPOL)
     &                        + JGR(IR,I1,LAM)*A_GB
C
                           GB_JFR(IR,IKM2,LAM,IPOL)
     &                        = GB_JFR(IR,IKM2,LAM,IPOL)
     &                        + JFR(IR,I1,LAM)*A_GB
C
                        END DO
C
                     END IF
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
                  END DO
C
               END DO
C-------------------------------------------------------------------- LM
            END DO
C-------------------------------------------------------------------- I1
         END DO
C-------------------------------------------------------------------- I2
      END DO
C                                                                    LAM
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
C=======================================================================
      IF ( IPRINT.LT.2 ) RETURN
C=======================================================================
C
      WRITE (6,99001) ROUTINE(1:LEN_TRIM(ROUTINE))
      DO LM = 1,NLMFP
         WRITE (6,99002) LM,KVLM(LM),(KGVLM(LM,IPOL),IPOL=1,3)
      END DO
C
      DO IPOL = 1,3
         WRITE (6,99003) IPOL
         WKM1(1:NKM,1:NKM) = KGV_ZGR(1:NKM,1:NKM,IPOL)
C
         CALL CMATSTRUCT('KGV_ZGR',WKM1,NKM,NKMMAX,3,3,0,1D-8,-6)
      END DO
C
99001 FORMAT (//,1X,79('*'),/,32X,'<',A,'>',/,1X,79('*'),//,10X,
     &        'flags for non-0 convoluted wave functions',//,10X,
     &        'LM  KVLM(LM)  KGVLM(LM,IPOL)',/)
99002 FORMAT (10X,I2,I6,7X,3I3)
99003 FORMAT (//,8X,'KGV_ZGR(IKM2,LAM,IPOL) for IPOL = ',I3,/)
      END
