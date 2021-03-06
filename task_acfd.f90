!BOP
!
! !ROUTINE: task_acfd
!
! !INTERFACE:
      subroutine task_acfd
      
! !DESCRIPTION:
!
! This subroutine calculate the ACFD total energy 
!
! !USES:

      use acfd
      use freq,     only: nomeg
      use kpoints,  only: nirkp,idikp,nkp
      use mixbasis, only: init_mixbasis,end_mixbasis
      use barcoul,  only: init_barcoul,end_barcoul,barcevtol,end_barcev,&
                          iop_coul_x, iop_coul_c
      use dielmat,  only: init_dielmat, end_dielmat
      use liboct_parser
      use modmpi
      
      
! !LOCAL VARIABLES:

      implicit none

      integer(4) :: iq,irq        ! (Counter) Runs over q-points
      integer(4) :: ierr 
      integer(4) :: iop_acfd
      integer :: iqfirst, iqlast
      
      real(8) :: tstart         ! Initial CPU-time of the subroutine
      real(8) :: tend         ! Final CPU-time of the subroutine
      real(8) :: time1,time2  ! Initial and final time of each called subroutine
      character(len=10)::blk_acfd="acfd",sname="task_acfd"
      

! !REVISION HISTORY:
!
! Created Sept. 20, 2009
!      
!EOP
!BOC     

      write(6,*) 
      write(6,*) "------------------ task:acfd-------------------------"
      write(6,*) 

!
! Read task-specific input parameters from *.ingw
!     iop_acfd = 0 - just calculate exact exchange (EXX)  energy
!                1 - EXX + RPA correlation energy
!
      ierr = loct_parse_isdef(blk_acfd)
      if(ierr.eq.1) then
        call loct_parse_block_int(blk_acfd,0,0,iop_acfd)
      else
        iop_acfd=1
      endif

#ifdef MPI
      write(6,*) "task_acfd: parallel q-loop"
      call mpi_set_range(nproc_row,myrank_row,nkp,1,iqfirst,iqlast)
      call mpi_set_range(nproc_col,myrank_col,nomeg,1,iom_first,iom_last, &
     &                  iom_cnts,iom_dspl)
      write(6,*) "myrank_row iqfirst,iqlast   =",myrank_row,iqfirst,iqlast
      write(6,*) "myrank_col,iom_first,iom_last=",myrank_col,iom_first,iom_last
#else
      write(6,*) "sequential q-loop"
      iqfirst=1
      iqlast=nkp
      iom_first=1
      iom_last=nomeg
#endif
      call init_acfd(iqfirst,iqlast,iom_first,iom_last) 

      do iq=iqfirst,iqlast   !! Loop over q-points.

        call init_mixbasis(iq) 
        call init_barcoul(iq)

        !! Calculate the q-dependent integration weights
        call bz_calcqdepw(iq)

        !! calculate the bare coulomb matrix and its square root
        call coul_barc(iq)

        if(myrank_col.eq.0) then 
          call coul_setev(iq,0.d0,iop_coul_x)
          call calcexhf(iq)
          call end_barcev
        endif 

        if(iop_acfd.ne.0) then 

          write(6,'(a,f8.4)')" -Use reduced basis,barcevtol =",barcevtol
          call coul_setev(iq,barcevtol,iop_coul_c)
          call init_dielmat(iq,iom_first,iom_last)

          ! Calculate the dielectric matrix
          call calceps(iq,iom_first,iom_last,0,-1,0,.false.) 

          call calcecrpa(iq,iom_first,iom_last) 
          call end_dielmat(iq)
        endif 

        !! Deallocate arrays with q-dependent sizes
        call end_barcev
        call end_barcoul(iq)
        call end_mixbasis
        call flushbuf(6)
      enddo ! iq
!
!    Calculate acfd correlation energy by integrating over q-points  
!
#ifdef MPI
        call mpi_sum_scalar(0,ex_hf,-1)
        call mpi_sum_scalar(0,ec_acfd,-1)
#endif

      if(myrank.eq.0) then
        write(6,*) 
        call boxmsg(6,'-',"ACFD Exc Summary (in Hartree)")
        write(6,'(a,f12.6)') " Exact exchange   =",ex_hf
        write(6,'(a,f12.6)') " ACFD Correlation =",ec_acfd
        write(6,'(a,f12.6)') " Total ACFD Exc   =",ex_hf+ec_acfd
      endif 
      call end_acfd

  101 format(5x,'Data for q-point nr.:',i4,//,5x,'Mixed basis:',/,5x, &
     &           'Number of atomic basis functions:       ',i4,/,5x,   &
     &           'Number of interstitial basis functions: ',i4,/,5x,   &
     &           'Total number of basis functions:        ',i4,/)
      
 1000 format('CPUTIME for ',A20,F16.2,' seconds')
      
      end subroutine task_acfd
!EOC      
