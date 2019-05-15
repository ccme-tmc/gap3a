!BOP
!
! !ROUTINE: calcescgw
!
! !INTERFACE:
      subroutine calcescgw_cutoff(isxc)
      
! !DESCRIPTION:
!
! This subroutine performs energy-only self-consistent gw calculatons,
!
! !USES:
      use bands,       only: nbandsgw,nspin,bande,ibgw,nbgw,eqp,&
     &                       eferqp,efermi
      use freq,        only: nomeg
      use kpoints,     only: nkp,nqp,nirkp,kirlist,kirvecs
      use minmmat,     only: init_minm,end_minm
      use mixbasis,    only: init_mixbasis,end_mixbasis
      use selfenergy,  only: init_selfenergy,end_selfenergy,  &
     &                       sigx,sigc,sxcmn
      use barcoul,     only: init_barcoul,end_barcoul,barcevtol,end_barcev,&
     &                       iop_coul_x,iop_coul_c
      use dielmat,     only: init_dielmat,end_dielmat
      use xcpot,       only: init_xcpot, end_xcpot
      use task,        only: lrestart,casename,nmax_sc
      use modmpi
      use liboct_parser

      
! !LOCAL VARIABLES:

      implicit none
      integer, intent(in) :: isxc 

      character(len=20):: sname="calcescgw_cutoff"
      integer(4) :: iq        ! (Counter) Runs over q-points
      integer(4) :: ierr
      integer(4) :: isc
      integer    :: iq_f, iq_l 
      logical    :: lconv, lread_eps 

! !REVISION HISTORY:
!
! Created Oct.02.2010 by Hong Jiang
!      
!EOP
!BOCa

!
!     allocate the arrays needed for the calculation of the self
!     energy (coulomb matrix, polarization matrix, etc)
!
      call init_selfenergy(0)
!
!     Calculate the integration weights using the linearized tetrahedron
!     method
!

!
!     Prepare for parallelization
!

      if(lrestart) then 
        lread_eps = .true.
      else
        lread_eps = .false.
      endif 

#ifdef MPI
      call mpi_set_group(nqp,nirkp)
      call mpi_set_range(nproc_col,myrank_col,nqp,1,iq_f,iq_l)
      write(6,*) "myrank_col, iq_f, iq_l   =", myrank_col,  iq_f,  iq_l
#else
      write(6,*) "sequential q-loop"
      iq_f=1
      iq_l=nqp
#endif

!
!     set up vxcnn
!
      if(myrank.eq.0) then
        call init_xcpot
        call w2k_calcvxcnn
      endif
 
!
!     Start self-consistent loop
!
      do isc=0,nmax_sc
        call scgw_check_conv(isc,lconv)
        if(lconv) exit  
!
!       Start loop over q
!
        sigc=0.d0 
        do iq=iq_f,iq_l
          call init_mixbasis(iq)
          call init_barcoul(iq)

          !! exchange part needs to be calculated only at the first
          !iteration 
          if(isc.eq.0) then 
            call coul_barc(iq, iop_coul_x)           !! bare Coulomb matrix 
            call coul_setev(iq,0.d0,iop_coul_x)
            call calcselfx(iq,-1)
            call end_barcev
          endif 

          call coul_barc(iq, iop_coul_c)           !! bare Coulomb matrix 
          call coul_setev(iq,barcevtol,iop_coul_c)

          if(isc.eq.0) then
            call init_minm(iq,0)
            call set_minm(iq,0)
          else
            call init_minm(iq,1)
          endif

          ! Calculate the polarization matrix
          call init_dielmat(iq,1,nomeg)
          call calceps(iq,1,nomeg,0,0,2,lread_eps)

          !Add the corresponding q-term to the correlation term of the
          !self-energy
          call calcselfc(iq,0)
          call end_dielmat(iq)
          call end_minm
          call end_barcev
          call end_barcoul(iq)
          call end_mixbasis
          call flushbuf(6)
        enddo 

#ifdef MPI
        if(nproc_col.gt.1.and.myrank_ra3.eq.0) then
          write(6,*) "sum selfe from different q-points !"
          mycomm=mycomm_col
          if(isc.eq.0) then
            call mpi_sum_array(0,sigx,nbandsgw,nirkp,nspin,mycomm)
          endif 

          if(isxc.eq.0) then !! GW correlation self-energy

            call mpi_sum_array(0,sigc,nomeg,nbandsgw,nirkp,nspin,mycomm)

          elseif(isxc.eq.2) then !! COHSEX correlation self-energy

            call mpi_sum_array(0,sigc,3,nbandsgw,nirkp,nspin,mycomm)

          endif

          write(6,*) "sum selfe done !"
        endif
#endif
        if(myrank.eq.0) call calceqp(isxc,0,isc)
        call scgw_update_bande
        call bandanaly(ibgw,nbgw,nirkp,kirvecs,&
     &                 bande(ibgw:nbgw,:,:),efermi,nspin,'e-GW')

      enddo  !! loop over isc

      if(myrank.eq.0) then
        if(isc.gt.nmax_sc) then
          write(6,*) "calcescgw: WARNING - Not converged!!"
        endif 
      endif 
      call end_selfenergy(0) 

      return
      end subroutine calcescgw
!EOC      
