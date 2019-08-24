!BOP
!
! !ROUTINE: task_eps
!
! !INTERFACE:
      subroutine task_eps 
      
! !DESCRIPTION:
!
! This task subroutine calculate dielectric matrix used for later GW calculations 
! !USES:

      use bands,       only: nbandsgw,nspin,bande,ibgw,nbgw,eqp,bande0, &
     &                       numin,nomax,nbmax,nbmaxpol,eferqp,efermi,metallic
      use barcoul,     only: init_barcoul,end_barcoul,barcevtol,        &
     &                       end_barcev,iop_coul_x,iop_coul_c, &
     &                       lcutoff_in_coul_barc
      use dielmat,     only: init_dielmat,end_dielmat,eps,epsw1,epsw2,head
      use freq,        only: nomeg,nomeg_blk
      use kpoints,     only: nkp,nqp,nirkp
      use mixbasis,    only: init_mixbasis,end_mixbasis,matsiz
      use mommat,      only: init_mommat,end_mommat
      use task,        only: lrestart,casename,nmax_sc
      use anisotropy,  only: aniten,init_aniso,end_aniso,iop_aniso,lmax_q0
      use bzinteg,     only: n_ang_grid
      use modmpi 
      use liboct_parser
      
! !LOCAL VARIABLES:

      implicit none

      integer :: ie        ! index for energy bands 
      integer :: ierr      ! error code 
      integer :: iq        ! index for q-points
      integer :: iq0       ! index for starting k-points, needed for the restart mode 
      integer :: irk       ! index for irreducible k-points 
      integer :: iop_minm   ! control the treatment of mwm  
      integer :: iop_sxc 
      integer :: isp       ! index for spin 
      logical :: lconv

      character(len=20):: sname = "task_eps"
      character(len=20):: blk_gw = "eps"
      real(8) :: tstart,tend        !! variables related to cputime 

      !! variables used for parallelization 
      integer :: iq_f,iq_l
      integer :: iom,nom,iom_f,iom_l 

! !EXTERNAL ROUTINES: 
      !external calceps_new
      external coul_barc

! !INTRINSIC ROUTINES: 
      intrinsic cpu_time      

! !REVISION HISTORY:
!
! Created 19.07.2008 by Hong Jiang
!      
!EOP
!BOCa

      call linmsg(6,'-',sname) 

#ifdef MPI
      write(6,*) "parallel q-loop"

      call mpi_set_group(nqp,nkp) 

      call mpi_set_range(nproc_col,myrank_col,nqp,1,iq_f,iq_l)
      write(6,*) "myrank_col iq_f,iq_l =",myrank_col,iq_f,iq_l
#else
      write(6,*) "sequential q-loop"
      iq_f=1
      iq_l=nqp
#endif
!
! Restart 
!
      iq0=iq_f-1 

      do iom = 1,nomeg,nomeg_blk
        iom_f = iom
        iom_l  = min(iom + nomeg_blk-1,nomeg) 
        nom = iom_l-iom_f+1

        do iq=iq0+1,iq_l

          call init_mixbasis(iq)
          call init_barcoul(iq)

          !! bare Coulomb matrix 
          if (lcutoff_in_coul_barc) then
            call coul_barc_cutoff(iq, iop_coul_c)
          else
            call coul_barc(iq)
          endif

          call coul_setev(iq,barcevtol,iop_coul_c)

          call init_dielmat(iq,iom_f,iom_l)  !! initialize
          if(iop_aniso.ne.-1)then
            call init_aniso(aniten,iq,matsiz,iom_f,iom_l,lmax_q0,n_ang_grid)
          endif


          !Calculate the dielectric matrix
          call calceps(iq,iom_f,iom_l,0,-1,2,.false.)
    
          if(iop_aniso.ne.-1) call end_aniso(iq,aniten)
          call end_dielmat(iq)
          call end_barcev
          call end_barcoul(iq)
          call end_mixbasis 
          call flushbuf(6) 
        enddo ! iq
      enddo ! iom 

      return
      
      end subroutine task_eps
!EOC
