!BOP
!
! !ROUTINE: task_mwm
!
! !INTERFACE:
      subroutine task_mwm 
      
! !DESCRIPTION:
!
! This task subroutine calculate M*W*M and save them in files
! !USES:

      use bands,       only: nspin,nbands_c,ibgw,nbgw
      use barcoul,     only: init_barcoul,end_barcoul,barcevtol,        &
     &                       end_barcev,iop_coul_c, lcutoff_in_coul_barc
      use dielmat,     only: init_dielmat,end_dielmat
      use freq,        only: nomeg,nomeg_blk
      use kpoints,     only: nkp,nqp,nirkp
      use mixbasis,    only: matsiz,init_mixbasis,end_mixbasis
      use modmpi      
      use liboct_parser
      
! !LOCAL VARIABLES:

      implicit none

      integer :: ierr      ! error code 
      integer :: iq        ! index for q-points
      integer :: irk       ! index for irreducible k-points 
      integer :: iom       ! index for frequency 
      integer :: iop_minm  ! control the treatment of mwm  
      integer :: isp       ! index for spin 

!     variables used for parallelization 
      integer :: iq_f,iq_l     !! lower and upper bound for iq  
      integer :: irk_f, irk_l  !! lower and upper bound for irk
      integer :: iom_f, iom_l   !! lower and upper bound for iom 
      integer :: ie2_f, ie2_l 

      complex(8),allocatable:: mwm(:,:,:) 
       

      character(len=20):: sname="task_mwm"
      character(len=20):: blk_mwm="mwm"
      real(8) :: tstart,tend        !! variables related to cputime 

! !REVISION HISTORY:
!
! Created 19.07.2008 by Hong Jiang
!      
!EOP
!BOC
      call linmsg(6,'-',sname) 

      call cpu_time(tstart)       
      iop_minm = -1 
#ifdef MPI
      call mpi_set_group(nqp,nirkp) 
      call mpi_set_range(nproc_col,myrank_col,nqp,1,iq_f,iq_l)
      call mpi_set_range(nproc_row,myrank_row,nirkp,1,irk_f,irk_l)
      call mpi_set_range(nproc_3rd,myrank_3rd,nbands_c,1,ie2_f,ie2_l) 
      write(6,*) "myrank_col, iq_f,  iq_l  =",myrank_col,iq_f,iq_l
      write(6,*) "myrank_row, irk_f, irk_l =",myrank_row,irk_f, irk_l
      write(6,*) "myrank_3rd, ie2_f, ie2_l =",myrank_3rd,ie2_f, ie2_l
#else
      write(6,*) "sequential q-loop"
      iq_f  = 1
      iq_l  = nqp
      irk_f = 1
      irk_l = nirkp
      ie2_f = 1 
      ie2_l = nbands_c
#endif

      !! to save memory, split freq points to blocks to save memory
      do iq=iq_f,iq_l
        call init_mixbasis(iq)
        call init_barcoul(iq)
        if (lcutoff_in_coul_barc) then
          call coul_barc_cutoff(iq, iop_coul_c)
        else
          call coul_barc(iq)
        endif
        call coul_setev(iq,barcevtol,iop_coul_c)

        do iom=1,nomeg,nomeg_blk 
          iom_f = iom 
          iom_l = min(iom+nomeg_blk-1,nomeg) 

          allocate(mwm(ie2_f:ie2_l,ibgw:nbgw,iom_f:iom_l)) 

          call init_dielmat(iq,iom_f,iom_l)  !! initialize
          call io_eps('r',iq,iom_f,iom_l,ierr) 

          do isp=1,nspin 
            do irk=irk_f, irk_l
              call calcmwm(iop_minm,isp,iq,irk,ie2_f,ie2_l,           &
     &                     iom_f,iom_l,mwm) 
            enddo 
          enddo
          deallocate(mwm)
          call end_dielmat(iq) 

        enddo !! iom
        call end_barcev
        call end_barcoul(iq)
        call end_mixbasis 
        call flushbuf(6) 
      enddo !! iq
      end subroutine task_mwm
!EOC
