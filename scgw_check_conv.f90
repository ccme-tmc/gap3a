      subroutine scgw_check_conv(isc,lconv)
      use bands,   only: bande,nspin,nomaxs,numins,metallic 
      use kpoints, only: nirkp
      use modmpi 
      use selfenergy, only: egk
      use task,    only: eps_sc
      implicit none
      integer,intent(in)  :: isc
      logical,intent(out)  :: lconv
      integer:: isp,nvm,ncm,ierr
      real(8) :: err
      real(8) :: egkold(nirkp,nspin)

      if(.not.allocated(egk)) allocate(egk(nirkp,nspin))

      if(myrank.eq.0) then
        lconv = .false.
        if(isc.gt.0) egkold=egk
        do isp=1,nspin
          nvm = nomaxs(isp)
          ncm = numins(isp)
          if(metallic) then
            egk(:,isp)=bande(nvm,:,isp)
          else
            egk(:,isp)=bande(ncm,:,isp)-bande(nvm,:,isp)
          endif
        enddo

        if(isc.gt.0) then
          err=abs(maxval(egk-egkold))
          write(6,'(a,i5,e12.4)') "#scgw: isc,err=",isc,err
          if(err.lt.eps_sc) lconv=.true.
        endif
      endif !! myrank.eq.0
#ifdef MPI
      call mpi_bcast(lconv,1,mpi_logical,0,mpi_comm_world,ierr)
#endif 
      end subroutine ! sub_check_conv
