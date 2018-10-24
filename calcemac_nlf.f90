!BOP
!
! !ROUTINE: calcemac_nlf
!
! !INTERFACE:
      subroutine calcemac_nlf(emac_a,iom0,nom)
      
! !DESCRIPTION:
!
! This task calculates macroscopic diectric function 
!
! !USES:

      use constants, only: hev
      use freq,      only: nomeg,omega,womeg,nomeg_blk
      use bands,     only: metallic,nomax,numin,nbmaxpol,nspin
      use kpoints,   only: nkp,nirkp
      use mixbasis,  only: matsiz,init_mixbasis,end_mixbasis
      use barcoul,   only: init_barcoul,end_barcoul,barcevtol,end_barcev
      use dielmat,   only: eps,epsw1,epsw2,init_dielmat,end_dielmat,&
     &                     emac,head,c0_head
      use freq,      only: nomeg
      use minmmat,   only: init_minm,end_minm
      use mommat,    only: init_mommat,end_mommat,iop_mommat
      use task,      only: taskname,casename
      use modmpi 
      use liboct_parser
      
! !LOCAL VARIABLES:

      implicit none
      integer, intent(in)   :: iom0 ! the index for starting frequency  
      integer, intent(in)   :: nom  ! the number of freq points 
      complex(8),intent(out):: emac_a(iom0:iom0+nom-1) 
 
      integer :: ierr,ip,comm,iom
      integer :: ikfirst,iklast,iomfirst,iomlast
      integer :: iomcnt(0:nproc_max), iomdsp(0:nproc_max) 
      real(8) :: wpl
      logical :: lprt=.false.
 
      character(20):: sname="calcemac_nlf"

      integer:: nom_p
      integer::rstype
!
!EOP
!BOC            

!     calculate emac without local field effect
#ifdef DEBUG
      lprt = .true.
#endif 
       
#ifdef MPI
      call mpi_set_group(nirkp,1) 
      call mpi_set_range(nproc_row,myrank_row,nirkp,1,ikfirst,iklast)
      call mpi_set_range(nproc_col,myrank_col,nom,iom0,iomfirst, &
     &                     iomlast,iomcnt,iomdsp)
      write(6,*) "myrank_col,iomfirst,iomlast=",myrank_col,iomfirst,iomlast
#else
      iomfirst= iom0
      iomlast = iom0+nom-1
      ikfirst = 1
      iklast  = nirkp
#endif
      nom_p=iomlast-iomfirst+1

      if(lprt) write(6,*) trim(sname)//" - init_dielmat"
      call init_dielmat(0,iomfirst,iomlast)

      !! Calculate the q-dependent integration weights
      write(6,*) "task_emac: bz_calcqdepw"
      call bz_calcqdepw(1)

      call init_mommat(1,nomax,numin,nbmaxpol,nirkp,nspin)

      if(iop_mommat.eq.0) then
        call calcmommat(0,1,nomax,numin,nbmaxpol,0)
      else
        call w2k_readmommat(1,nomax,numin,nbmaxpol)
      endif

      if(lprt) write(6,*) trim(sname)//" - calchead "
      call calchead(ikfirst,iklast,iomfirst,iomlast)

#ifdef MPI
      if(nproc_row.gt.1) then
        write(6,*) sname//":Sum head from different k's"
        call mpi_sum_scalar(0,c0_head,mycomm_row)
        call mpi_sum_array(0,head(iomfirst:iomlast),nom_p,mycomm_row)
      endif
#endif
      if(nproc_col.gt.1) then
#ifdef MPI
        write(6,*) sname//":collect emac data for different freq"
        call MPI_Type_contiguous(1,MPI_DOUBLE_COMPLEX,rstype,ierr)
        call MPI_Type_commit(rstype,ierr)
        call MPI_GatherV(head,nom_p,rstype,emac_a,&
     &                  iomcnt,iomdsp,rstype,0,mycomm_col,ierr)
        call MPI_Type_free(rstype,ierr)
#endif
      else 
        emac_a = head 
      endif 
      call end_dielmat(0) 

      return

      end subroutine calcemac_nlf
!EOC
