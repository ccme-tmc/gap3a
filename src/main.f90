!BOP
!
! !ROUTINE: main
!
! !INTERFACE:
      program main

! !DESCRIPTION:
!
! This is the main program unit.  In this version, it includes only "task". 
! All test functionalities have been removed 
! 
! !USES:

      use bzinteg,     only: init_bzinteg
      use bands,       only: nspin
      use barcoul,     only: init_barcoul
      use kpoints,     only: nqp,nirkp
      use task,        only: progname,taskname,time_minm,time_lapack,&
     &                       time_evec,time_eps,time_selfx,time_selfc,&
     &                       time_mwm,time_coul,time_vxc,time_aniso
      use xcpot,       only: init_xcpot, end_xcpot
      use anisotropy,  only: iop_aniso
      use modmpi

!
! !DEFINED PARAMETERS:

      implicit none
!
! !LOCAL VARIABLES:
!
      integer(4) :: info
      integer(4) :: isp
      real(8) :: tstart,tend

! !REVISION HISTORY:
!       
! Created Jan. 2007 by HJ
!
! !TO DO:
!   - Calculation of the dielectric function
!   - Calculation of the screened coulomb potential
!   - and more...
!
! !FILES USED:
!
!EOP
!BOC
!---------------------------------------
!       Initializations
!---------------------------------------
      progname = 'GAP 3.0a'
      call init_mpi()
      call cpu_time(tstart)

      call readingw           ! Read the gw main input file
      call set_struct         !  set struct information
      call set_lapwlo         !* set up LAPW+lo basis 
      if(taskname.eq."nvf") then
        call task_nvf         !* Generate new vector according to GW quasi-particle energy correction
        stop
      endif

      call bz_setkmesh        !* set up k- and q-points grids
      call set_ipw            !* operations related to interstitial plane wave 
      call set_ksref          !* setup up Kohn-Sham reference system including KS energy, xc potential, KS vectors 
      call freq_set_fgrid     !* setup the frequency mesh and integration weights
      call set_mixbasis       !* setup mixed basis set 
      call bz_setkiw          !* normal BZ integration weights
      call init_bzinteg       !* Initialize libbzint

      if(taskname.eq.'chkbz') then 
        write(6,*) "check BZ integration"
        call task_chkbz
        stop
      endif 

      call init_barcoul(0)    !* Initialize bare Coulomb matrix calculations 

      if(taskname.eq.'chkinit') stop 

      select case(trim(taskname))

!      case(id_lapw) 
!        call task_lapw

!      case(id_evec) 
!        call task_evec

!      case(id_mixf)
!        call task_mixf 
!
!  Calculate the macroscopic dielectric function
!
      case("coul")     !! Test barc Coulomb matrix 
        call task_coul
      case("acfd") 
        call task_acfd 
      case("gw_aniso") 
        call task_gw_aniso
      case ("emac")   !! calc macroscopic dielectric function
        call task_emac
      case("ppgw")    !! post-processing GW 
        call task_ppgw
      case("acont")    !! redo analytic continuation  
        call task_acont
      case("eps")
        call task_eps 
      case("gw") 
        call task_gw 
      case("gwb")  ! run GW in a batch mode 
        call task_eps()
        call task_mwm()
        call task_gw()
      case("gwsc") 
        call task_gwsc
      case("gw2wann") 
        call task_gw2wann 
      case("ldau") 
        call task_ldau
      case("mwm") 
        call task_mwm
      case("crpa")
        call task_crpa
      case default
        write(6,*) "main: Task= ",trim(taskname)," is not supported!"
      end select 

      call cpu_time(tend)
      write(6,10) "Total GW Calculations  ",tend-tstart
      write(6,10) "  Minm matrix          ",time_minm 
      write(6,10) "  LAPACK/BLAS routines ",time_lapack
      write(6,10) "  Expand eigen-vectors ",time_evec
      write(6,10) "  Dielectric  matrix   ",time_eps 
      write(6,10) "  LDA/GGA Vxc          ",time_vxc
      write(6,10) "  selfx matrix         ",time_selfx
      write(6,10) "  selfc matrix         ",time_selfc
      write(6,10) "  Bare Coulomb         ",time_coul
      write(6,10) "  M*W*M matrix         ",time_mwm
      !if(iop_aniso.ne.-1) then
      !  write(6,10) "  Anisotropy realted   ",time_aniso
      !endif

 10   format('CPUTIME for    ',A40,f16.2,' seconds')      
      call io_cleanup
      call end_mpi()
      end program 
!EOC
