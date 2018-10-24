!BOP
!
! !ROUTINE: calc_sxcmn
!
! !INTERFACE:
      subroutine calc_sxcmn(isxc,isc,isym) 
      
! !DESCRIPTION:
!
! This subroutine calculate the full sxc matrix including the summation
! over q. It can used for self-consistent GW, or for one-shot GW (G0W0)
! calculation considering the contribution of the full Sxc matrix 
!
!   isxc --- indicating which approximation to the selfenergy
!       0/1/2  -- GW / Exchange-only / COHSEX
!   isc  --- count the self-consistent cycle  
!       -1 -- for one-shot GW calculation 
!        0 -- first iteration of SCGW 
!      n>1 -- 
!   isym --- indicating whether Sxc for the full (0) or irreducible (1)
!        BZ are calculated 
! !USES:

      use barcoul,   only: init_barcoul,end_barcoul,barcevtol,&
     &                     end_barcev,iop_coul_c,iop_coul_x
      use bands,     only: nbandsgw,nspin,nomax,numin,nbmaxpol
      use dielmat,   only: init_dielmat,end_dielmat
      use freq,      only: nomeg,nomeg_blk
      use kpoints,   only: nkp,nqp,nirkp
      use minmmat,   only: init_minm, end_minm 
      use mixbasis,  only: init_mixbasis,end_mixbasis
      use selfenergy,only: sigm,sigm_f,qpwf_dmat,iop_dvh,iop_sc,n_sigm
      use task,      only: casename,save_prefix,lrestart
      use modmpi       
      
! !LOCAL VARIABLES:

      implicit none
      integer, intent(in):: isxc     !! which approx. selfenergy  
      integer, intent(in):: isc      !! sc iteration counter 
      integer, intent(in):: isym     !! 0/1 

      integer :: ierr          ! 
      integer :: iq,  iq_f, iq_l    !! the range of q-indices for each "row" process 
      integer :: iom, iom_f, iom_l  !! the range of freq index 
      integer :: isym_eps = 0  !! summation of the full BZ is used when calculating eps 
      integer :: iop_eps = 2   !! the matrix of eps^{-1) - 1 is calculated
      integer :: nktot         !! total number of k-points 
      
      logical :: lprt=.false.  
      logical :: lread_eps=.false.

      character(len=10) :: sname='calc_sxcmn' 
      character(len=5)  :: flag_sxcmn
      
      real(8):: tstart,tend   ! bookkeeping CPU-time 

#ifdef DEBUG
      lprt = .true.
#endif 
      if(lrestart) then
        call io_sxcmn('r','f',0,isxc,isym,ierr)

        if(ierr.eq.0) then 
          write(6,*) "Restarting mode: read Sxcmn from files"
          return 
        else
          call errmsg(.true.,sname,"Fail to read Sxcmn") 
        endif 
      endif

      if(isym.eq.0) then 
        nktot = nkp
      else
        nktot = nirkp
      endif

      !  When iop_sc == 1, a self-consistent GW0 is required,
      !  in which eps is read from the existing files after the first
      !  iteration  
      if(isc.gt.0.and.iop_sc.eq.1) then
        lread_eps = .true.
      else
        lread_eps = .false.
      endif 

      !! calculate the density matrix corresponding to the qp wave functions
      call calcqpwf_dmat

      sigm = 0.d0 
      if(isc.eq.0) sigm_f=0.d0 

#ifdef MPI
      write(6,*) "parallel q-loop"
      call mpi_set_group(nqp,nktot)
      call mpi_set_range(nproc_col,myrank_col,nqp,1,iq_f,iq_l)
      write(6,*) " myrank_col   =", myrank_col
      write(6,*) " iq_f,iq_l   =", iq_f,iq_l
#else
      write(6,*) "sequential q-loop"
      iq_f = 1
      iq_l = nqp
#endif

      do iq=iq_f,iq_l
        call init_mixbasis(iq)
        call init_barcoul(iq)
        call coul_barc(iq)    !! bare Coulomb matrix
        
        !! set v-diag basis for the exchange selfenergy matrix 
        call coul_setev(iq,0.d0,iop_coul_x) 

        !! matrix elements for the Hartree potential 
        if(iq.eq.1.and.iop_dvh.gt.0.and.isc.ge.0)  call calcvhmn(isc)

        call calcselfxmn(iq,isc,isym)           
        call end_barcev

        if(isxc.eq.1) then   !! exchange-only self-energy 
          !!Deallocate arrays with q-dependent sizes
          call end_barcoul(iq)
          call end_mixbasis
          cycle 
        endif 

        !! Calculate the correlation self-energy 
        call coul_setev(iq,barcevtol,iop_coul_c)

        if(isc.eq.-1) then 
          call init_minm(iq,-1) 
        else if(isc.eq.0) then 
          call init_minm(iq,0) 
          call set_minm(iq,0) 
        else 
          call init_minm(iq,1) 
        endif 

        !! Dielectric matrix
        if(nomeg_blk.lt.nomeg) then   ! the freq-loop is split into blocks 
          write(6,*) "The freq-loop is split into blocks"
          do iom = 1, nomeg, nomeg_blk
            iom_f = iom
            iom_l = min(nomeg, iom + nomeg_blk - 1)
            call init_dielmat(iq,iom_f,iom_l)
            call calceps(iq,iom_f,iom_l,0,isc,2,lread_eps)
            call calc_mwm3(iq,iom_f, iom_l, isc,isym)
            call end_dielmat(iq)
          enddo
          !! correlation self-energy matrix 
          call calcselfcmn_mwm(iq,isxc,isc,isym)

        else            !! the freq-loop is not split 
          call init_dielmat(iq,1,nomeg)
          call calceps(iq,1,nomeg,0,isc,2,lread_eps)
          call calcselfcmn(iq,isxc,isc,isym)
          call end_dielmat(iq)
        endif

        call io_sxcmn('w','f',iq,isxc,isym,ierr)

        call end_minm 
        call end_barcev
        call end_barcoul(iq)
        call end_mixbasis 
        call flushbuf(6) 
      enddo 
#ifdef MPI
      if(nproc_col.gt.1.and.myrank_ra3.eq.0) then
        write(6,*) trim(sname)//": sum selfe from different q-points !"
        call mpi_sum_array(0,sigm,nbandsgw,nbandsgw,n_sigm,nktot,nspin, &
     &            mycomm_col)

        if(isc.eq.0) then
          call mpi_sum_array(0,sigm_f,nbandsgw,nbandsgw,n_sigm,nktot,   &
     &          nspin,mycomm_col)
        endif
      endif
#endif
      if(myrank.eq.0)then
        call io_sxcmn('w','f',0,isxc,isym,ierr)
      endif 

      end subroutine calc_sxcmn
!EOC      
