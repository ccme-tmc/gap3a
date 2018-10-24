      subroutine scgw_update_bande() 
!
! This sobroutine update bands according to quasi-particle energies
! for states below ibgw and above nbgw, assume same corrections as those
! for the ibgw band and nbgw, respectively
!
      use bands,      only: bande,efermi,eferqp,eferqp0,nspin,ibgw,nbgw,&
     &                      nbmax,eqp,nbandsgw 
      use selfenergy, only: iop_esgw0 
      use kpoints, only: nirkp
      use core,    only: eigcore,ncg
      use modmpi
      implicit none
!     integer,intent(in):: iop 
      integer:: isp, irk,ierr
      real(8):: de_cb,de_vb,sfactor=1.0

      if(myrank.eq.0) then 
        do isp=1,nspin
          de_vb=sum(eqp(ibgw,:,isp)-bande(ibgw,:,isp))/nirkp
          de_cb=sum(eqp(nbgw,:,isp)-bande(nbgw,:,isp))/nirkp
          do irk=1,nirkp
            bande(1:ibgw-1,irk,isp)=bande(1:ibgw-1,irk,isp) + de_vb
            bande(nbgw+1:nbmax,irk,isp)=bande(nbgw+1:nbmax,irk,isp)+de_cb 
            bande(ibgw:nbgw,irk,isp) = eqp(ibgw:nbgw,irk,isp) 
          enddo
        enddo
   
        if(iop_esgw0.eq.1) then 
          bande = bande - eferqp
          efermi = 0.d0
          eferqp0 = eferqp0 + eferqp
          if(ncg.gt.0) eigcore = eigcore + de_vb - eferqp
        endif 

      endif 

#ifdef MPI
      call mpi_bcast(bande,nbmax*nirkp*nspin,mpi_double_precision,0,&
     &                 mpi_comm_world,ierr)

      if(ncg.gt.0) then 
        call mpi_bcast(eigcore,size(eigcore),mpi_double_precision,0,&
     &                 mpi_comm_world,ierr)
      endif 

#endif 

      end subroutine
