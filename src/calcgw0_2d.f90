!BOP
!
! !ROUTINE: calcgw0
!
! !INTERFACE:
      subroutine calcgw0_2d(isxc,iref)
      
! !DESCRIPTION:
! This subroutine performs "standard" G0W0 and GW0 calculatons, it can
! be called in two situations, as controlled by iref:
!   iref == 0 -> a standard G0W0 calculation based on some "KS" reference
!   iref == 1 -> G0W0 and GW0 based on some approximate self-consistent GW
! isxc -- determine which approximation to self-energy
!    0 --- GW
!    1 --- exchange-only
!    2 -- - COHSEX
!
! !USES:

      use bands,       only: nbandsgw,nspin,bande,ibgw,nbgw,eqp,bande0, &
     &                       numin,nomax,nbmax,nbmaxpol,eferqp,efermi
      use barcoul,     only: init_barcoul,end_barcoul,barcevtol,        &
     &                       end_barcev,iop_coul_x,iop_coul_c, lcutoff_in_coul_barc
      use dielmat,     only: init_dielmat,end_dielmat
      use freq,        only: nomeg,omega
      use kpoints,     only: nkp,nqp,nirkp,kirvecs
      use mixbasis,    only: init_mixbasis,end_mixbasis,matsiz
      use mommat,      only: init_mommat,end_mommat
      use selfenergy,  only: init_selfenergy,end_selfenergy,sigx,sigc,  & 
     &                       flag_sxc,fn_selfe,isxc 
      use xcpot,       only: init_xcpot,end_xcpot
      use task,        only: lrestart,casename,nmax_sc,time_aniso
      use anisotropy,  only: init_aniso,end_aniso,iop_aniso,aniten,&
                             lmax_q0
      use bzinteg,     only: n_ang_grid
      use modmpi      
      use liboct_parser
      
! !LOCAL VARIABLES:

      implicit none
      integer,intent(in):: isxc, iref

      integer :: iq        ! index for q-points
      integer :: iq0       ! index for starting k-points, needed for the restart mode 
      integer :: iop_minm  ! control the treatment of mwm  
      integer :: ierr      ! error code

!     variables used for parallelization
      integer :: iq_f,iq_l     !! lower and upper bound for iq
      logical :: lread_eps 

      character(len=20):: sname="calcgw0"

! !REVISION HISTORY:
!
! Created 19.07.2008 by Hong Jiang
!      
!EOP
!BOCa

      if(iref.eq.0) then 
        iop_minm = -1 
      else 
        iop_minm = 0 
      endif 

      if(lrestart) then 
        lread_eps = .true.
      else
        lread_eps = .false.
      endif 

      if(isxc.eq.2) then !! COHSEX 
        nomeg=1
        omega(1)=0.d0
      endif

!
!     allocate the arrays needed for the calculation of the self
!     energy (coulomb matrix, polarization matrix, etc)
!
      call init_selfenergy(0) 

#ifdef MPI
      call mpi_set_group(nqp,nirkp) 
      call mpi_set_range(nproc_col,myrank_col,nqp,1,iq_f,iq_l)
      write(6,*) "myrank_col, iq_f, iq_l   =", myrank_col,  iq_f,  iq_l
#else
      write(6,*) "sequential q-loop"
      iq_f=1
      iq_l=nqp
#endif
      !! set up <psi_nk|Vxc^DFT | psi_nk> 
      if(myrank.eq.0) then
        call init_xcpot
        call w2k_calcvxcnn
      endif

      do iq=iq_f,iq_l
        call init_mixbasis(iq)
        call init_barcoul(iq)

        !! bare Coulomb matrix 
        if (lcutoff_in_coul_barc) then
          call coul_barc_cutoff(iq, iop_coul_x)
        else
          call coul_barc(iq)
        endif

        !! exchange self-energies 
        call sub_calc_sigx 

        !! Calculate GW or COHSEX correlation self-energy 
        if(isxc.eq.0.or.isxc.eq.2.or.isxc.eq.3) then
          if (lcutoff_in_coul_barc) then
            call coul_barc_cutoff(iq,iop_coul_c)
          endif
          call sub_calc_sigc
        endif

        call end_barcoul(iq)
        call end_mixbasis 
        call flushbuf(6) 
      enddo ! iq

#ifdef MPI
      if(nproc_col.gt.1.and.myrank_ra3.eq.0) then
        write(6,*) "sum selfe from different q-points !"

        call mpi_sum_array(0,sigx,nbandsgw,nirkp,nspin,mycomm_col)

        if(isxc.eq.0) then !! GW correlation self-energy
          call mpi_sum_array(0,sigc,nomeg,nbandsgw,nirkp,nspin,mycomm_col)
        elseif(isxc.eq.2) then !! COHSEX correlation self-energy
          call mpi_sum_array(0,sigc,3,nbandsgw,nirkp,nspin,mycomm_col)
        elseif(isxc.eq.3) then 
          call mpi_sum_array(0,sigc,1,nbandsgw,nirkp,nspin,mycomm_col)
        endif 
        write(6,*) "sum selfe done !"
      endif
#endif

      if(myrank.eq.0)then
        call io_sxcmn('w','d',0,isxc,1,ierr)
        call calceqp(isxc,0,-1)  

        call io_eqp('w',isxc,0,flag_sxc)
        call io_eqp('w',isxc,1,flag_sxc)

        call bandanaly(ibgw,nbgw,nirkp,kirvecs,  &
     &         bande0(ibgw:nbgw,:,:),efermi,nspin,'KS')

        call bandanaly(ibgw,nbgw,nirkp,kirvecs,eqp,eferqp,nspin,flag_sxc)
      endif
!
!     perform energy-only self-consistent GW0 calculation 
!
      if(isxc.eq.0) call sub_scgw0 

      call end_selfenergy(0)
      if(myrank.eq.0) call end_xcpot
      return

      contains 
      

!-----------------------------------------------------------------------
!          internal subroutine for exchange self-energies              !
!-----------------------------------------------------------------------                            
        subroutine sub_calc_sigx
        integer :: ierr_sx    ! how to handle Selfx matrix 
        if(lrestart) then
          call io_sxcmn('r','d',iq,1,1,ierr_sx)
          if(ierr_sx.eq.0) then
            write(6,*) " Use existing Sx data!"
          endif
        else
          ierr_sx = 1
        endif

        if(ierr_sx.ne.0) then
          call coul_setev(iq,0.d0,iop_coul_x)
          call calcselfx(iq,iop_minm)

          if(myrank_ra3.eq.0) then
            call io_sxcmn('w','d',iq,1,1,ierr)
          endif

          call end_barcev
        endif
        end subroutine 

!-----------------------------------------------------------------------
!          internal subroutine for correlation self-energies           !
!-----------------------------------------------------------------------                            
        subroutine sub_calc_sigc
        integer:: ierr_sc
        real(8) :: time1, time2

        if(lrestart) then
          call io_sxcmn('r','d',iq,isxc,1,ierr_sc)
        else
          ierr_sc = 1
        endif

        if(ierr_sc.eq.0) then
          write(6,*) " Use existing Sc data!"
        else
          call coul_setev(iq,barcevtol,iop_coul_c)
          call init_dielmat(iq,1,nomeg)  !! initialize
          if(iop_aniso.ne.-1)then
            time_aniso = 0.0
            call cpu_time(time1)
            call init_aniso(aniten,iq,matsiz,1,nomeg,lmax_q0,n_ang_grid)
            call cpu_time(time2)
            time_aniso = time_aniso + time2 - time1
          endif

          if (iop_coul_c.eq.2)then
            !call calceps_2d(iq,1,nomeg,0,-1,2,lread_eps)
            call calceps_aniso(iq,1,nomeg,0,-1,2,lread_eps)
          else
            call calceps_aniso(iq,1,nomeg,0,-1,2,lread_eps)
          endif

          if (iop_coul_c.eq.2)then
            !call calcselfc_2d(iq,iop_minm)
            call calcselfc_2d(iq,iop_minm)
          else
            call calcselfc(iq,iop_minm)
          endif

          if(myrank_ra3.eq.0) then
            call io_sxcmn('w','d',iq,isxc,1,ierr)
          endif

          if(iop_aniso.ne.-1)then
            call end_aniso(iq,aniten)
          endif

          call end_dielmat(iq)
          call end_barcev
        endif
        end subroutine 


!----------------------------------------------------------------------
!       Internal subroutine to run energy-only self-consistent GW0      !
!----------------------------------------------------------------------
        subroutine sub_scgw0
        implicit none
        integer(4) :: isc   ! index for the self-consistent iteration
        logical :: lconv

        do isc=0,nmax_sc
          !!  update bande on the root process 
          call scgw_update_bande
          call scgw_check_conv(isc,lconv)
          if(lconv) exit

          !! calculate sigc using stored mwm
          sigc=0.d0
          lrestart = .true.
          do iq=iq_f,iq_l
            if (iop_coul_c.eq.2)then
              !call calcselfc_2d(iq,iop_minm)
              call calcselfc_2d(iq,0)
            else
              call calcselfc(iq,0)
            endif
          enddo

#ifdef MPI
          if(nproc_col.gt.1.and.myrank_row.eq.0.and.myrank_3rd.eq.0) then
            call mpi_sum_array(0,sigc,nbandsgw,nirkp,nomeg,nspin,mycomm_col)
          endif
#endif
          if(myrank.eq.0) call calceqp(0,0,isc)
        enddo  !! loop over isc

        if(myrank.eq.0) then
          if(isc.gt.nmax_sc) then
            write(6,*) "calcscgw0: WARNING - Not converged!!"
          endif
          call io_eqp('w',0,0,"GW0")
          call io_eqp('w',0,1,"GW0")
          call bandanaly(ibgw,nbgw,nirkp,kirvecs,eqp,eferqp,nspin,"GW0")
        endif
        return
        end subroutine 
      
      end subroutine calcgw0_2d
!EOC      
