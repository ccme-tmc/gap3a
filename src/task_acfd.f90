!BOP
!
! !ROUTINE: task_acfd
!
! !INTERFACE:
      subroutine task_acfd
      
! !DESCRIPTION:
!
! This subroutine calculate the ACFD total energy 
!  !USES:

      use acfd
      use freq,     only: nomeg
      use kpoints,  only: nirkp,idikp,nkp
      use mixbasis, only: init_mixbasis,end_mixbasis
      use barcoul,  only: init_barcoul,end_barcoul,barcevtol,end_barcev,&
                          iop_coul_x, iop_coul_c
      use dielmat,  only: init_dielmat, end_dielmat
      use liboct_parser
      use modmpi
      use task,     only: casename, fid_outgw
      
      
! !LOCAL VARIABLES:

      implicit none

      integer(4) :: iq,irq        ! (Counter) Runs over q-points
      integer(4) :: ierr 
      integer(4) :: iop_acfd
      integer :: iqfirst, iqlast
      
      real(8) :: tstart         ! Initial CPU-time of the subroutine
      real(8) :: tend         ! Final CPU-time of the subroutine
      real(8) :: time1,time2  ! Initial and final time of each called subroutine
      character(len=10)::blk_acfd="acfd", sname="task_acfd"
      logical :: ldbg=.true.
      

! !REVISION HISTORY:
!
! Created Sept. 20, 2009
!      
!EOP
!BOC     

      write(6,*) 
      write(6,*) "------------------ task:acfd-------------------------"
      write(6,*)

! Read wien2k total and exchange-correlation energy 
      if(myrank.eq.0) then
        open(unit=999,file=trim(casename)//'.exc',status='old',iostat=ierr)
        call errmsg(ierr.ne.0, sname, &
     &        "Energy file not found: "//trim(casename)//'.exc')
        read(999, '(A)') 
        read(999, *) exc_lda
        read(999, *) etot_lda
        ! convert to Hartree
        exc_lda = exc_lda / 2.0D0
        etot_lda = etot_lda / 2.0D0
        close(999)
      endif
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
      write(6,*) "task_acfd: parallel q-loop and w-loop"
      call mpi_set_group(nkp, nomeg)
     !call mpi_set_group(nkp,nirkp)

     ! call mpi_set_range(nproc_row,myrank_row,nkp,1,iqfirst,iqlast)
     ! call mpi_set_range(nproc_col,myrank_col,nomeg,1,iom_first,iom_last, &
     !&                  iom_cnts,iom_dspl)
     ! write(6,*) "nproc_row = ", nproc_row
     ! write(6,*) "myrank_row, iqfirst, iqlast     =", &
     !&    myrank_row,iqfirst,iqlast
     ! write(6,*) "nproc_col = ", nproc_col
     ! write(6,*) "myrank_col, iom_first, iom_last =", &
     !&    myrank_col,iom_first,iom_last
      call mpi_set_range(nproc_col,myrank_col,nkp,1,iqfirst,iqlast)
      call mpi_set_range(nproc_row,myrank_row,nomeg,1,iom_first,iom_last)
     ! call mpi_set_range(1,0,nomeg,1,iom_first,iom_last)
      write(6,*) "nproc_col= ", nproc_col
      write(6,*) "myrank_col, iqfirst, iqlast     =", &
     &    myrank_col,iqfirst,iqlast
      write(6,*) "nproc_row = ", nproc_row
      write(6,*) "myrank_row, iom_first, iom_last =", &
     &    myrank_row,iom_first,iom_last
#else
      write(6,*) "sequential q-loop"
      iqfirst=1
      iqlast=nkp
      iom_first=1
      iom_last=nomeg
#endif
      call init_acfd(iqfirst,iqlast,iom_first,iom_last) 

      do iq=iqfirst,iqlast   !! Loop over q-points.
      !do iq=iqlast,iqfirst,-1   !! Back loop for debug
      !  write(*,*) "iq=",iq
        call init_mixbasis(iq) 
        call init_barcoul(iq)
        !! Calculate the q-dependent integration weights
        call bz_calcqdepw(iq)
        !! calculate the bare coulomb matrix and its square root
        call coul_barc(iq)

        if(myrank_row.eq.0) then 
          call coul_setev(iq, 0.d0, iop_coul_x)
          call calcexhf(iq)
          call end_barcev
        endif 

        if(iop_acfd.ne.0) then 
          write(fid_outgw,'(a,f8.4)')&
     &        " - Use reduced basis, barcevtol = ",barcevtol
          call coul_setev(iq,barcevtol,iop_coul_c)
          call init_dielmat(iq,iom_first,iom_last)
          ! Calculate the dielectric matrix
          call calceps(iq,iom_first,iom_last,0,-1,0,.false.) 
          call calcecrpa(iq,iom_first,iom_last) 
          if(ldbg)then
            write(fid_outgw, *)
            write(fid_outgw, *) "Accumulated ACFD correlation: ", ec_acfd
          endif
          call end_dielmat(iq)
          call end_barcev
        endif 

        !! Deallocate arrays with q-dependent sizes
        call end_barcoul(iq)
        call end_mixbasis
        call flushbuf(fid_outgw)
      enddo ! iq
!
!    Calculate acfd correlation energy by integrating over q-points  
!
#ifdef MPI
      if(ldbg)then
        write(fid_outgw,*)"!!!DEBUG!!!"
        write(fid_outgw,*)"ACFD info in process ",myrank
        write(fid_outgw,*)"Responsible for iq : ",iqfirst," to ",iqlast
        write(fid_outgw,*)"            for iw : ",iom_first," to ",iom_last
        write(fid_outgw,*)"   Exact exchange = ",ex_hf
        write(fid_outgw,*)" ACFD Correlation = ",ec_acfd
      endif
      call mpi_sum_scalar(0,ex_hf,mycomm)
      call mpi_sum_scalar(0,ec_acfd,mycomm)
#endif

      if(myrank.eq.0) then
        call boxmsg(fid_outgw,'-',"ACFD Exc Summary (in Hartree)")
        write(fid_outgw,'(a,f20.8)') "   Exact exchange = ",ex_hf
        write(fid_outgw,'(a,f20.8)') " ACFD Correlation = ",ec_acfd
        write(fid_outgw,'(a,f20.8)') "      LDA/GGA Exc = ",exc_lda
        write(fid_outgw,'(a,f20.8)') "         ACFD Exc = ",ex_hf+ec_acfd
        write(fid_outgw,'(a,f20.8)') "     LDA/GGA Etot = ",etot_lda
        write(fid_outgw,'(a,f20.8)') "        ACFD Etot = ",etot_lda - exc_lda &
     &      + ex_hf + ec_acfd
      endif 
      call end_acfd

  101 format(5x,'Data for q-point nr.:',i4,//,5x,'Mixed basis:',/,5x, &
     &           'Number of atomic basis functions:       ',i4,/,5x,   &
     &           'Number of interstitial basis functions: ',i4,/,5x,   &
     &           'Total number of basis functions:        ',i4,/)
      
 1000 format('CPUTIME for ',A20,F16.2,' seconds')
      
      end subroutine task_acfd
!EOC      
