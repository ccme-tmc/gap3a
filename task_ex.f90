!BOP
!
! !ROUTINE: task_ex
!
! !INTERFACE:
      subroutine task_ex 
      
! !DESCRIPTION:
!
! This task subroutine performs "standard" G0W0 and GW0 calculatons
! !USES:

      use bands,       only: nbandsgw,nspin,bande,ibgw,nbgw,eqp,bande0, &
     &                       numin,nomax,nbmax,nbmaxpol,eferqp,efermi
      use barcoul,     only: init_barcoul,end_barcoul,barcevtol,        &
     &                       end_barcev,iop_coul_x,iop_coul_c
      use dielmat,     only: init_dielmat,end_dielmat
      use freq,        only: nomeg
      use kpoints,     only: nkp,nqp,nirkp,kirlist,kirvec
      use mixbasis,    only: init_mixbasis,end_mixbasis
      use mommat,      only: init_mommat,end_mommat
      use selfenergy,  only: init_selfenergy,end_selfenergy,  &
     &                       sigx,sigc,sigsx,sigch,sxcmn,flag_sxc,fn_selfe
      use xcpot,       only: init_xcpot, end_xcpot,vxcnn
      use task,        only: lrestart,casename,nmax_sc
      use modmpi      
      use liboct_parser
      
! !LOCAL VARIABLES:

      implicit none

      integer :: ie        ! index for energy bands 
      integer :: ierr      ! error code 
      integer :: iq        ! index for q-points
      integer :: iq0       ! index for starting k-points, needed for the restart mode 
      integer :: irk       ! index for irreducible k-points 
      integer :: iop_mwm   ! control the treatment of mwm  
      integer :: iop_minm  ! control the treatment of mwm  
      integer :: iop_sxc   ! which approximation to the self-energy
      integer :: iop_eps   ! how to handle dielectric matrix 
      integer :: isp       ! index for spin 
      logical :: lconv
      logical :: lprt=.true.

      character(len=20):: sname="task_ex"
      character(len=20):: blk_gw="ex"
      real(8) :: tstart,tend        !! variables related to cputime 

!     variables used for parallelization
      integer :: iq_f,iq_l     !! lower and upper bound for iq
      integer :: irk_f, irk_l  !! lower and upper bound for irk

! !REVISION HISTORY:
!
! Created 10.08.2011 by Hong Jiang
!      
!EOP
!BOCa

      iop_minm = -1 
      iop_sxc = 1
!
!     allocate the arrays needed for the calculation of the self
!     energy (coulomb matrix, polarization matrix, etc)
!
      call init_selfenergy
      call init_xcpot

#ifdef MPI
      write(6,*) "task_ex: parallel "
      call mpi_set_range(nproc_row,myrank_row,nqp,1,iq_f,iq_l)
      call mpi_set_range(nproc_col,myrank_col,nirkp,1,irk_f,irk_l)
      write(6,*) "myrank_row, iq_f,  iq_l  =",myrank_row,iq_f,iq_l
      write(6,*) "myrank_col, irk_f, irk_l =",myrank_col,irk_f, irk_l
#else
      write(6,*) "task_ex: sequential"
      iq_f  = 1
      iq_l  = nqp
      irk_f = 1
      irk_l = nirkp
#endif
      iq0=iq_f-1 
      if(lrestart) then
        if(lprt) write(6,*) "Restart mode"
        call io_sxcmn('r','d',iq0,iop_sxc,myrank)
        if(iq0.eq.0) then 
          iq0=iq_l 
        elseif(iq0.lt.0) then 
          iq0=iq_f-1 
        endif 
      endif 

!
!     set up vxcnn
!
      if(myrank.eq.0) then 
        call io_vxcmn('r','d',ierr)
        if(ierr.ne.0) then 
          write(6,*) " Fail to read vxc from the vxcnn file"
          write(6,*) " -- calc vxcnn from KS vxc"
          call w2k_calcvxcnn
        endif 
        if(lprt) write(6,*) "task_ex: calc vxcnn done"
      endif 

      do iq=iq0+1,iq_l

        call init_mixbasis(iq)
        call init_barcoul(iq)

        call coul_barc(iq)           !! bare Coulomb matrix 
        call coul_setev(iq,0.d0,iop_coul_x) 

        call calcselfx(iq,-1) 

        call io_sxcmn('w','d',iq,iop_sxc,myrank)
!
!       Deallocate arrays with q-dependent sizes
!
        call end_barcev
        call end_barcoul(iq)
        call end_mixbasis 
        call flushbuf(6) 
      enddo ! iq

#ifdef MPI
      !! when iop_gw==1, the q-loop is skipped so that sumselfe is not needed
      if(nproc_row.gt.1.and.myrank_col.eq.0) then
        write(6,*) "sum selfe from different q-points !"
        call mpi_sum_array(0,sigx,nbandsgw,nirkp,nspin,mycomm_row)
        write(6,*) "sum selfe done !"
      endif
#endif

      if(myrank.eq.0)then
        !! Write the final selfenergies 
        call io_sxcmn('w','d',0,iop_sxc,'')

        if(lprt) write(6,*) "task_ex: calc qp energy  "
        call calceqp(iop_sxc,0,-1)  

        call io_eqp('w',iop_sxc,0,flag_sxc)
        call io_eqp('w',iop_sxc,1,flag_sxc)

        call bandanaly(ibgw,nbgw,nirkp,kirvec,  &
     &         bande0(ibgw:nbgw,:,:),efermi,nspin,'KS')

        call bandanaly(ibgw,nbgw,nirkp,kirvec,eqp,eferqp,nspin,flag_sxc)
      endif

      call end_selfenergy
      call end_xcpot
      return

      end subroutine task_ex
!EOC      
