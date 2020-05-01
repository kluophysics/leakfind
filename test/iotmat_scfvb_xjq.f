      subroutine wrtmat_scfvb(itemp,ntemp,id_ene,eryd,lnovb)
c
      use mod_constants, only : c0, ci
      use mod_files, only : ifilcbwf, ifilcbwf_lhs, ifilcbwf_sph,
     &                      reclngwf, reclngwf_sph, lrecreal8
      use mod_angmom, only : nkmmax, nkm, tsst, msst, ssst, mezz, mezj
      use mod_types, only : nthost
      use mod_energy, only : efermi
      use mod_scfvb_cpa_sigma, only : ifilscftv, ifilscftv2
      use mod_dmft_ldau, only : dmftsigma, dmftsolver
      use mod_mpi, only : mpi_id
c
      implicit none
c input
      logical lnovb
      integer itemp, ntemp, id_ene
      complex*16 eryd
c local
      character*4 cid
      logical lopen
      integer reclng_tmat,rec_tmat, it, ikm, ie, ne
      real*8 estep
      complex*16 p
c
      inquire(unit=ifilcbwf,opened=lopen)
      if(lopen) close(ifilcbwf)
      inquire(unit=ifilcbwf_lhs,opened=lopen)
      if(lopen) close(ifilcbwf_lhs)
      inquire(unit=ifilcbwf_sph,opened=lopen)
      if(lopen) close(ifilcbwf_sph)
      open(unit=ifilcbwf,status='scratch',form='unformatted',
     &     access='direct',recl=reclngwf)
      call set_ifil_lhs(ifilcbwf,ifilcbwf_lhs)
      open (ifilcbwf_lhs,status='scratch',form='unformatted',
     &      access='direct',recl=reclngwf)
      call set_ifil_sph(ifilcbwf,ifilcbwf_sph)
      open (ifilcbwf_sph,status='scratch',form='unformatted',
     &      access='direct',recl=reclngwf_sph)
c
      p = c0
      call runssite(.true.,0,0,ifilcbwf,.true.,eryd,p,0,
     &              tsst,msst,ssst,mezz,mezj,'none      ')
c
      write(cid,*) id_ene
      cid=adjustl(cid)
c
      if(itemp==1) then
        if(lnovb==.false.) then
          reclng_tmat = 2*lrecreal8 * nkmmax * nkmmax
          open(unit=ifilscftv,file='scftv_e'//trim(cid)//'.dat',
     &         form='unformatted',
     &         access='direct',recl=reclng_tmat,
     &         status='replace',err=100)
          open(unit=ifilscftv2,file='scftv2_e'//trim(cid)//'.dat',
     &         form='formatted',
     &         status='unknown')
        else
          open(unit=ifilscftv2,file='tv_novb_e'//trim(cid)//'.dat',
     &         form='formatted',
     &         status='unknown')
        endif
      endif
c
      do it=1,nthost
        if(lnovb==.false.) then
          rec_tmat = (itemp-1)*nthost + it
          write(ifilscftv,rec=rec_tmat) tsst(:,:,it)
        endif
        write(ifilscftv2,'(2i5)') itemp,it
        do ikm=1,nkmmax
          write(ifilscftv2,'(2e22.14)') tsst(ikm,ikm,it)
        enddo
      enddo ! it
c
      close(ifilcbwf)
      close(ifilcbwf_lhs)
      close(ifilcbwf_sph)
      if(itemp==ntemp) then
        if(lnovb==.false.) then
          close(ifilscftv)
        endif
        close(ifilscftv2)
      endif
c
      print *,'finish writing t-mat'
c
      return
 100  stop 'failed to open filscftv'
      end

      subroutine read_scftv(ivt,tscf)
c
      use mod_files, only : lrecreal8
      use mod_types, only : nthost
      use mod_angmom, only : nkmmax
      use mod_thermal, only : i_temp_lat, n_temp_lat, nvibra, nfluct
      use mod_scfvb_cpa_sigma, only : lscfvb, ifilscftv
      use mod_mpi, only : mpi_id
c
      implicit none
c input
      integer ivt
c inout
      complex*16 tscf(nkmmax,nkmmax)
c local
      logical lfilopen
      integer reclng_tmat, rec_tmat
      integer ikm
      complex*16 d_tl(nkmmax), t0(nkmmax,nkmmax)
c
      if(.not. lscfvb) return
c
      inquire(ifilscftv,opened=lfilopen)
      if(.not. lfilopen) then
        reclng_tmat = 2*lrecreal8 * nkmmax * nkmmax
        open(unit=ifilscftv,file='scftv.dat',form='unformatted',
     &       access='direct',recl=reclng_tmat,status='old')
      endif
c
      if(nfluct .ne. 1)
     &  stop 'nfluct must be 1 for scfvb for current version'
c
      t0(:,:) = tscf(:,:)
c
      if(mpi_id==0) then
        print *,'**************************************************'
        print *,'Read self-consistent t_v of ivt ',ivt
        print *,'**************************************************'
      endif
c
      rec_tmat = (i_temp_lat-1)*nthost*nvibra + ivt
c
      read(ifilscftv,rec=rec_tmat) tscf(:,:)
c
      do ikm = 1, nkmmax
        d_tl(ikm) = tscf(ikm,ikm) - t0(ikm,ikm)
        if(abs(d_tl(ikm)) > 0.00000001) then
          if(mpi_id==0) then
            print *, 'ikm,  ',ikm
            print *, 'tscf: ',tscf(ikm,ikm),', t0: ',t0(ikm,ikm)
          endif
c
c=======================================================================
c If lmax depends on sites or types, code should be modified at least for test case
c=======================================================================
c
        endif
      enddo
c
      return
      end subroutine read_scftv
