      subroutine wrmregmirr_scfvb(itemp,ntemp,id_ene,eryd,check_novb)
c
      use mod_constants, only : c0, ci
      use mod_files, only : lrecreal8, reclngwf, reclngwf_sph
      use mod_sites, only : iqbot, iqtop
      use mod_angmom, only : nkm, nkmmax, tsst, msst, ssst, mezz, mezj
      use mod_types, only : nthost, ctl, itbot, ittop
      use mod_energy, only : efermi
      use mod_scfvb_cpa_sigma, only : ifilmreg, ifilmirr, 
     &                                ifilwfvba, ifilwfvbb
      use mod_mpi, only : mpi_id
c
      implicit none
c input
      logical check_novb
      integer itemp, ntemp, id_ene
      complex*16 eryd
c local
      logical lopen
      character*4 cid
      integer ifilwfvba_lhs, ifilwfvbb_lhs, ifilwfvba_sph, ifilwfvbb_sph
      integer iza, izb, it, ntlim
      integer reclng_mreg, reclng_mirr, rec_mreg, rec_mirr
      real*8 c
      complex*16 eryda, erydb, p_dummy
      complex*16 mreg_tab(nkmmax,nkmmax,3),mreg_tba(nkmmax,nkmmax,3),
     &           mirr_2(nkmmax,nkmmax,3,3),mirr_3(nkmmax,nkmmax,3,3),
     &           mirr_4(nkmmax,nkmmax,3,3)
c
      if(mpi_id .ne. 0) return
c
      if(check_novb .and. itemp .ne. 1) return
c
      p_dummy = c0
      ntlim = nthost
c
      write(cid,*) id_ene
      cid=adjustl(cid)
c
      if(check_novb .or. ((.not. check_novb) .and. itemp==1)) then
        ntlim = 1
c
        reclng_mreg = 2*lrecreal8 * nkmmax * nkmmax * 3 * 2
        reclng_mirr = 2*lrecreal8 * nkmmax * nkmmax * 3 * 3 * 3
        open(unit=ifilmreg,file='mreg_e'//trim(cid)//'.dat',
     &       form='unformatted',
     &       access='direct',recl=reclng_mreg,status='replace',err=100)
        open(unit=ifilmirr,file='mirr_e'//trim(cid)//'.dat',
     &       form='unformatted',
     &       access='direct',recl=reclng_mirr,status='replace',err=100)
c
        if(check_novb) then
          print *,'**************************************************' 
          print *,'Write mreg and mirr without scfvb for check'
          print *,'**************************************************'
          write(*,*) 'itbot = ',itbot,', ittop = ',ittop
          write(*,*) 'iqbot = ',iqbot,', iqtop = ',iqtop
          write(*,*)
        else
          return
        endif
      endif
c
c=======================================================================
c                           Loop on energy
c=======================================================================
c
      do iza=1,2
      if(iza==1) then
      eryda = eryd - ci*0d0
      else
      eryda = eryd + ci*0d0
      endif
c
      do izb=1,2
      if(izb==1) then
      erydb = eryd - ci*0d0
      else
      erydb = eryd + ci*0d0
      endif
c
c=======================================================================
c                  Generate wavefunctions Z and J
c=======================================================================
c
        open(unit=ifilwfvba,status='scratch',form='unformatted',
     &       access='direct',recl=reclngwf)
        call set_ifil_lhs(ifilwfvba,ifilwfvba_lhs)
        open (ifilwfvba_lhs,status='scratch',form='unformatted',
     &        access='direct',recl=reclngwf)
        call set_ifil_sph(ifilwfvba,ifilwfvba_sph)
        open (ifilwfvba_sph,status='scratch',form='unformatted',
     &        access='direct',recl=reclngwf_sph)
        open(unit=ifilwfvbb,status='scratch',form='unformatted',
     &       access='direct',recl=reclngwf)
        call set_ifil_lhs(ifilwfvbb,ifilwfvbb_lhs)
        open (ifilwfvbb_lhs,status='scratch',form='unformatted',
     &        access='direct',recl=reclngwf)
        call set_ifil_sph(ifilwfvbb,ifilwfvbb_sph)
        open (ifilwfvbb_sph,status='scratch',form='unformatted',
     &        access='direct',recl=reclngwf_sph)
c
        call runssite(.true.,1,1,ifilwfvba,.true.,eryda,p_dummy,0,
     &                tsst,msst,ssst,mezz,mezj,'none      ')
        call runssite(.true.,1,1,ifilwfvbb,.true.,erydb,p_dummy,0,
     &                tsst,msst,ssst,mezz,mezj,'none      ')
c
        do it=1,ntlim
c
          c = ctl(it,1)
c
          call me_alf_alf(1,nkm,ifilwfvba,eryda,1,nkm,ifilwfvbb,erydb,
     &                    .false.,it,mreg_tab,mreg_tba,
     &                    mirr_2,mirr_3,mirr_4,c,3)
c
c=======================================================================
c              transform and renormalize matrix elements
c=======================================================================
c
          call trans_renorm_mreg_mirr_scfvb(mreg_tab,mreg_tba,
     &      mirr_2,mirr_3,mirr_4,it,iza,izb)
c
c=======================================================================
c                    Write mreg and mirr in files
c=======================================================================
c
          if(check_novb) then
            rec_mreg = ntemp*2*2*nthost + 
     &                 (iza-1)*2 + (izb-1) + it
          else
            rec_mreg = (itemp-1)*2*2*nthost + (iza-1)*2*nthost +
     &                 (izb-1)*nthost + it
          endif
          rec_mirr = rec_mreg
c
          print *,'write mreg and mirr for rec ',rec_mreg
c
          write(ifilmreg,rec=rec_mreg) mreg_tab, mreg_tba
          write(ifilmirr,rec=rec_mirr) mirr_2, mirr_3, mirr_4
c
        enddo
c
        close(ifilwfvba)
        close(ifilwfvbb)
c
      enddo
      enddo
c
      if(itemp==ntemp .and. .not. check_novb) then
        close(ifilmreg)
        close(ifilmirr)
      endif
c
      return
 100  print *,'failed to open filmreg or filmirr'
      stop
      end
      
      subroutine read_mregmirr_scfvb(iza,izb,ivt,
     &  mreg_tab,mreg_tba,mirr_2,mirr_3,mirr_4)
c
      use mod_types,only:nthost
      use mod_angmom,only:nkmmax
      use mod_thermal,only:i_temp_lat,nvibra,nfluct
      use mod_sig,only:nspr
      use mod_scfvb_cpa_sigma,only:lscfvb,ifilmreg,ifilmirr
      use mod_mpi,only:mpi_id
c
      implicit none
c input
      integer iza,izb,ivt
      complex*16 mreg_tab(nkmmax,nkmmax,3),mreg_tba(nkmmax,nkmmax,3),
     &           mirr_2(nkmmax,nkmmax,3,3),mirr_3(nkmmax,nkmmax,3,3),
     &           mirr_4(nkmmax,nkmmax,3,3)
c local
      integer rec_mregmirr
c
      if(.not. lscfvb) return
c
      if(nspr .ne. 1) 
     &  stop 'for scfvb, nspr should be 1 for current version'
c
      if(nfluct .ne. 1)
     &  stop 'nfluct must be 1 for scfvb for current version'
c
      rec_mregmirr = (i_temp_lat-1)*2*2*nthost*nvibra +
     &  (iza-1)*2*nthost*nvibra + (izb-1)*nthost*nvibra + ivt
c
      if(mpi_id==0)
     &  print *,'Read mreg and mirr of rec ',rec_mregmirr
c
      read(ifilmreg,rec=rec_mregmirr) mreg_tab,mreg_tba
      read(ifilmirr,rec=rec_mregmirr) mirr_2,mirr_3,mirr_4
c
      return
      end
      
      subroutine read_mreg_scfvb(iza,izb,ivt,
     &  mreg_tab,mreg_tba)
c
      use mod_types,only:nthost
      use mod_angmom,only:nkmmax
      use mod_thermal,only:i_temp_lat,nvibra,nfluct
      use mod_sig,only:nspr
      use mod_scfvb_cpa_sigma,only:lscfvb,ifilmreg
      use mod_mpi,only:mpi_id
c
      implicit none
c input
      integer iza,izb,ivt
      complex*16 mreg_tab(nkmmax,nkmmax,3),mreg_tba(nkmmax,nkmmax,3)
c local
      integer rec_mreg
c
      if(.not. lscfvb) return
c
      if(mpi_id==0) 
     &  print *,'Read mreg in sigme_aux of ivt ',ivt
c
      if(nspr .ne. 1) 
     &  stop 'for scfvb, nspr should be 1 for current version'
c
      if(nfluct .ne. 1)
     &  stop 'nfluct must be 1 for scfvb for current version'
c
      rec_mreg = (i_temp_lat-1)*2*2*nthost*nvibra +
     &  (iza-1)*2*nthost*nvibra + (izb-1)*nthost*nvibra + ivt
c
      if(mpi_id==0)
     &  print *,'Read mreg of rec ',rec_mreg
c
      read(ifilmreg,rec=rec_mreg) mreg_tab,mreg_tba
c
      return
      end
      
      subroutine checkmregmirr_noscfvb(iza,izb,it,
     &  mreg_tab,mreg_tba,mirr_2,mirr_3,mirr_4)
c
      use mod_types,only:nthost
      use mod_angmom,only:nkmmax
      use mod_thermal,only:n_temp_lat,nvibra
      use mod_scfvb_cpa_sigma,only:ifilmreg,ifilmirr
      use mod_mpi,only:mpi_id
c
      implicit none
c input
      integer iza,izb,it
      complex*16 mreg_tab(nkmmax,nkmmax,3),mreg_tba(nkmmax,nkmmax,3),
     &           mirr_2(nkmmax,nkmmax,3,3),mirr_3(nkmmax,nkmmax,3,3),
     &           mirr_4(nkmmax,nkmmax,3,3)
c local
      integer rec_mregmirr
      integer ikm1, ikm2, mue, nue
      real*8 dr,di
      complex*16 mreg_tab_read(nkmmax,nkmmax,3),
     &           mreg_tba_read(nkmmax,nkmmax,3),
     &           mirr_2_read(nkmmax,nkmmax,3,3),
     &           mirr_3_read(nkmmax,nkmmax,3,3),
     &           mirr_4_read(nkmmax,nkmmax,3,3)
c
      if(mpi_id==0)
     &  print *,'check mreg and mirr without scfvb of ivt ',it
c
      rec_mregmirr = n_temp_lat*2*2*nthost*nvibra +
     &  (iza-1)*2 + izb-1 + it
      read(ifilmreg,rec=rec_mregmirr) mreg_tab_read,mreg_tba_read
      read(ifilmirr,rec=rec_mregmirr) mirr_2_read,mirr_3_read,
     &                                mirr_4_read
c
      do mue=1,3
        do ikm1=1,nkmmax
          do ikm2=1,nkmmax
            dr = abs(real(mreg_tab_read(ikm2,ikm1,mue))-
     &           real(mreg_tab(ikm2,ikm1,mue)))
            di = abs(aimag(mreg_tab_read(ikm2,ikm1,mue))-
     &           aimag(mreg_tab(ikm2,ikm1,mue)))
            if(dr>1e-6 .or. di>1e-6) then
              print *,'mreg_tab differs from that read for rec',
     &rec_mregmirr
              print *,'ikm2 = ',ikm2,'ikm1 = ',ikm1,'mue = ',mue
              print *,'mreg_tab_read = ',mreg_tab_read(ikm2,ikm1,mue)
              print *,'mreg_tab = ',mreg_tab(ikm2,ikm1,mue)
              print *,aimag(mreg_tab_read(ikm2,ikm1,mue))
              print *,aimag(mreg_tab(ikm2,ikm1,mue))
              print *,di
              stop
            endif
            dr = abs(real(mreg_tba_read(ikm2,ikm1,mue))-
     &           real(mreg_tba(ikm2,ikm1,mue)))
            di = abs(aimag(mreg_tba_read(ikm2,ikm1,mue))-
     &           aimag(mreg_tba(ikm2,ikm1,mue)))
            if(dr>1e-6 .or. di>1e-6) then
              print *,'mreg_tba differs from that read for rec',
     &rec_mregmirr
              print *,'ikm2 = ',ikm2,'ikm1 = ',ikm1,'mue = ',mue
              print *,'mreg_tba_read = ',mreg_tba_read(ikm2,ikm1,mue)
              print *,'mreg_tba = ',mreg_tba(ikm2,ikm1,mue)
              stop
            endif
          enddo
        enddo
      enddo
c
      do mue=1,3
        do nue=1,3
          do ikm1=1,nkmmax
            do ikm2=1,nkmmax
              dr = abs(real(mirr_2_read(ikm2,ikm1,nue,mue))-
     &             real(mirr_2(ikm2,ikm1,nue,mue)))
              di = abs(aimag(mirr_2_read(ikm2,ikm1,nue,mue))-
     &             aimag(mirr_2(ikm2,ikm1,nue,mue)))
              if(dr>1e-6 .or. di>1e-6)
     &          stop 'mirr_2 by sigme is different from that read'
              dr = abs(real(mirr_3_read(ikm2,ikm1,nue,mue))-
     &             real(mirr_3(ikm2,ikm1,nue,mue)))
              di = abs(aimag(mirr_3_read(ikm2,ikm1,nue,mue))-
     &             aimag(mirr_3(ikm2,ikm1,nue,mue)))
              if(dr>1e-6 .or. di>1e-6)
     &          stop 'mirr_3 by sigme is different from that read'
              dr = abs(real(mirr_4_read(ikm2,ikm1,nue,mue))-
     &             real(mirr_4(ikm2,ikm1,nue,mue)))
              di = abs(aimag(mirr_4_read(ikm2,ikm1,nue,mue))-
     &             aimag(mirr_4(ikm2,ikm1,nue,mue)))
              if(dr>1e-6 .or. di>1e-6)
     &          stop 'mirr_4 by sigme is different from that read'
            enddo
          enddo
        enddo
      enddo
c
      return 'mreg and mirr without scfvb have no problem'
      end
      
      subroutine trans_renorm_mreg_mirr_scfvb(mreg_tab,mreg_tba,
     &      mirr_2,mirr_3,mirr_4,it,iza,izb)
c
      use mod_constants, only : c0
      use mod_angmom, only : nkmmax, nkm, lmat3
      use mod_sites, only : iqat, nqmax, drotq, mrotq
      use mod_types, only : ntmax
      use mod_calcmode, only : moments_rotated
c
      implicit none
c input
      integer it, iza, izb
      complex*16 mreg_tab(nkmmax,nkmmax,3),mreg_tba(nkmmax,nkmmax,3),
     &           mirr_2(nkmmax,nkmmax,3,3),mirr_3(nkmmax,nkmmax,3,3),
     &           mirr_4(nkmmax,nkmmax,3,3)
c local
      integer mue, nue, iq, m, n
      complex*16 wkm1(nkmmax,nkmmax)
      complex*16 mreg_tab_tmp(nkmmax,nkmmax,3),
     &           mreg_tba_tmp(nkmmax,nkmmax,3),
     &           mirr_2_tmp(nkmmax,nkmmax,3,3),
     &           mirr_3_tmp(nkmmax,nkmmax,3,3),
     &           mirr_4_tmp(nkmmax,nkmmax,3,3)
c
c-----------------------------------------------------------------------
c           transformation from spherical to cartesian coordinates
c              from now on        ipol = 1,2,3  ==  (x),(y),(z)
c-----------------------------------------------------------------------
c
      call cmat_convert_polar(mreg_tab,'S>C')
      call cmat_convert_polar(mreg_tba,'S>C')
      call cmat_convert_polar2(mirr_2,'S>C')
      call cmat_convert_polar2(mirr_3,'S>C')
      call cmat_convert_polar2(mirr_4,'S>C')
c
c-----------------------------------------------------------------------
c  rotation from local to global frame for non-collinear magnetisation
c-----------------------------------------------------------------------
c
      if ( moments_rotated ) then
c
         iq = iqat(1,it)
c
         mreg_tab_tmp(:,:,:)
     &      = mreg_tab(:,:,:)
         mreg_tba_tmp(:,:,:)
     &      = mreg_tba(:,:,:)
c
         mirr_2_tmp(:,:,:,:) = mirr_2(:,:,:,:)
         mirr_3_tmp(:,:,:,:) = mirr_3(:,:,:,:)
         mirr_4_tmp(:,:,:,:) = mirr_4(:,:,:,:)
c
         call me_rotate_reg(iq,nqmax,drotq,mrotq,
     &                      mreg_tab_tmp,mreg_tba_tmp,
     &                      mreg_tab,mreg_tba)
c
         call me_rotate_irr(iq,nqmax,drotq,mrotq,
     &                      mirr_2_tmp,mirr_3_tmp,mirr_4_tmp,
     &                      mirr_2,mirr_3,mirr_4)
c
      end if
c
c----------------------------------------------------------------------
c   in case of energy on real axis multiply by appropriate
c   phase factors (lmat) to obtain mes
c   m(e+,e-), m(e-,e+) and m(e-,e-)
c----------------------------------------------------------------------
c
      wkm1(:,:) = c0
      m = nkmmax
      n = nkm
c
      do mue = 1,3
c
         if ( iza.eq.1 ) then
c
            call cmatmul(n,m,lmat3,mreg_tab(1,1,mue),wkm1)
            mreg_tab(:,:,mue) = wkm1(:,:)
c
            call cmatmul(n,m,mreg_tba(1,1,mue),lmat3,wkm1)
            mreg_tba(:,:,mue) = wkm1(:,:)
c
            do nue = 1,3
c
               call cmatmul(n,m,lmat3,mirr_3(1,1,mue,nue),wkm1)
c
               call cmatmul(n,m,wkm1,lmat3,mirr_3(1,1,mue,nue))
c
            end do
c
         end if
c
         if ( izb.eq.1 ) then
c
            call cmatmul(n,m,mreg_tab(1,1,mue),lmat3,wkm1)
            mreg_tab(:,:,mue) = wkm1(:,:)
c
            call cmatmul(n,m,lmat3,mreg_tba(1,1,mue),wkm1)
            mreg_tba(:,:,mue) = wkm1(:,:)
c
            do nue = 1,3
c
               call cmatmul(n,m,lmat3,mirr_2(1,1,mue,nue),wkm1)
c
               call cmatmul(n,m,wkm1,lmat3,mirr_2(1,1,mue,nue))
c
            end do
c
         end if
c
      end do
c
      return
      end
