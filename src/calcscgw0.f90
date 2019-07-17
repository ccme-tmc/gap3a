!BOP
!
! !ROUTINE: calcscgw0
!
! !INTERFACE:
      subroutine calcscgw0 
! !USES:
      use bands,       only: nbandsgw,nspin,bande,ibgw,nbgw,metallic,eferqp,eqp
      use freq,        only: nomeg  
      use kpoints,     only: nirkp,nkp,kirvecs
      use selfenergy,  only: sigc
      use task,        only: nmax_sc,eps_sc,casename,lrestart 
      use modmpi      
! !LOCAL VARIABLES:
      implicit none
      integer(4) :: iq,isc   ! (Counter) Runs over q-points,self-consistent iteration and spin 
      integer(4) :: iprt=1,ierr
      integer :: iqfirst,iqlast 
      real(8) :: err,egk(nirkp,nspin),egkold(nirkp,nspin)
      logical :: lconv
!
! Created 16.09.2008 by Hong Jiang
!      
!EOP
!BOC

#ifdef MPI
      mycomm=mycomm_row
      if(myrank_col.ne.0) return 
      call mpi_set_range(nproc_row,myrank_row,nkp,1,iqfirst,iqlast)
      write(6,*) "myrank_row iqfirst,iqlast   =",myrank_row,iqfirst,iqlast
#else
      iqfirst=1
      iqlast=nkp
#endif

      do isc=0,nmax_sc
        !!  update bande on the root process 
        call scgw_update_bande
        call scgw_check_conv(isc,lconv)
        if(lconv) exit

        !! calculate sigc using stored mwm
        sigc=0.d0
        lrestart = .true.
        do iq=iqfirst,iqlast
          call calcselfc(iq,0)
        enddo

#ifdef MPI
        if(nproc_row.gt.1) then
          call mpi_sum_array(0,sigc,nbandsgw,nirkp,nomeg,nspin,mycomm)
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
      
      end subroutine calcscgw0
!EOC     
