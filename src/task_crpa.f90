!BOP
!
! !ROUTINE: task_crpa
!
! !INTERFACE:
      subroutine task_crpa
      
! !DESCRIPTION:
!
! This subroutine calculate frequency-dependent U from the constrained
! RPA approach 
!
! !USES:
      use constants,only: HeV
      use crpa,     only:iop_crpa,iop_wf,iop_intercell,ncell,rcell,  &
     &                   rcell_cart,init_crpa,end_crpa,nlmorb,vmat,  &
     &                   nsp_crpa,umat,iop_pln_phase 
      use freq,     only: nomeg,nomeg_blk
      use kpoints,  only: nkp,nqp
      use barcoul,  only: init_barcoul,end_barcoul,barcevtol,end_barcev,&
                          iop_coul_c,iop_coul_x
      use dielmat,  only: init_dielmat, end_dielmat,iop_mask_eps
      use minmmat,  only: init_minm, end_minm
      use mixbasis, only: init_mixbasis,end_mixbasis,matsiz
      use liboct_parser
      use modmpi
      use task,     only: casename,spflag,savdir,lrestart

! !LOCAL VARIABLES:
      implicit none

      integer:: ic 
      integer:: iq, iq_f, iq_l
      integer:: iom, iom_f,iom_l
      integer:: iop_minm 
      integer:: lldim 
      integer:: ierr 
      
      real(8):: tini         ! Initial CPU-time of the subroutine
      real(8):: tend         ! Final CPU-time of the subroutine
      real(8):: time1,time2  ! Initial and final time of each called subroutine
 
      character(len=10):: blk_crpa="crpa",sname="task_crpa"
      character(len=20):: blk_cell="CRPA_intercell"
      character(len=80):: fn_wfpln,fn_umat
      character(len=10),external::Int2Str

      logical:: ldbg=.true.


! !REVISION HISTORY:
!
! Created Sept. 20, 2009
!      
!EOP
!BOC     

!
! Read task-specific input parameters 
!
! iop_crpa ---> indicate which cRPA scheme 
!       -1 - calculate the matrix of the bare Coulomb interaction (v)
!        0 - calculate the matrix of the screened Coulomb interaction (W)
!        1 - cRPA U with the mask approach 
!        2 - cRPA U with the weighted mask approach  (as in Sasioglu et al. PRB 83,121101(R)(2011))
!        3 - the projection approach (not implemented)
!        4 - the impurity U approach (not implemented) 
!        5 -- the disentanglement approach  by Miyake et al.(PRB 80, 155134(2009) (not implemented)
!
! iop_wf ---> option about the Wannier function 
!        0 -- use the internally generated WF (not implemented yet!)
!        1 -- use WF generated by dmftproj 
!        2 -- use WF generated by wien2wannier + wannier90 (not implemented yet!)

      iop_minm = -1 

      iop_crpa = 0
      iop_wf   = 1
      iop_intercell = 0
      ierr = loct_parse_isdef(blk_crpa)
      if (ierr.eq.1) then
        call loct_parse_block_int(blk_crpa,0,0,iop_crpa)
        call loct_parse_block_int(blk_crpa,0,1,iop_wf)
        call loct_parse_block_int(blk_crpa,0,2,iop_intercell)
      endif
    
      call sub_set_intercell

      if(iop_crpa.eq.0) then 
        iop_mask_eps = 0
      elseif(iop_crpa.eq.2) then 
        iop_mask_eps = 2
      endif  

     !! Read Wannier projectors from case.pln
      if(iop_wf.eq.1) then   !! read WF projectors from dmftproj.x 
        call crpa_readpln
      elseif(iop_wf.eq.2) then !! read WF from a Wannier90 chk file 
        call crpa_readw2w
      else
        write(6,*) "ERROR: not implemented iop_wf!"
      endif

      call init_crpa(0) 

#ifdef MPI
      write(6,*) "task_crpa: parallel q-loop"
      call mpi_set_group(nqp,nkp)
      call mpi_set_range(nproc_col,myrank_col,nqp,1,iq_f,iq_l)
      write(6,*) "myrank_col, iq_f, iq_l   =", myrank_col,  iq_f,  iq_l
#else
      write(6,*) "sequential q-loop"
      iq_f=1
      iq_l=nkp
#endif
      do iq =iq_f, iq_l  !  iq_f,iq_l   !! Loop over q-points.
        call init_mixbasis(iq) 
        call init_minm(iq,iop_minm) 

        call init_barcoul(iq)
        call bz_calcqdepw(iq)                     !! Calc the q-dependent integration weights
        call coul_barc(iq, iop_coul_c)            !! bare coulomb matrix
        call coul_setev(iq,barcevtol,iop_coul_c)  !! filter eigenvectors of barcoul

        call init_crpa(iq)                        !! initializes mill

        call crpa_calcmill(iq,1,nkp,0,iop_minm)   !! calc M^i_{L,L'}

        call crpa_calcvmat(iq)                    !! calc the bare Coulomb matrix in the WF basis 

        if(iop_crpa.ge.0) then
          do iom = 1, nomeg, nomeg_blk
            iom_f = iom
            iom_l = min(nomeg, iom + nomeg_blk - 1)
            call init_dielmat(iq,iom_f,iom_l)

            if(iop_crpa.le.2) then 
              call calceps(iq,iom_f,iom_l,0,-1,1,lrestart)
            else 
              write(6,*) "CRPA: iop_crpa=",iop_crpa  
            endif
            call crpa_calcumat(iq,iom_f,iom_l)
            call end_dielmat(iq)
          enddo 
        endif 

        !! Deallocate arrays with q-dependent sizes
        call end_crpa(iq) 
        call end_barcev
        call end_barcoul(iq)
        call end_minm
        call end_mixbasis
      enddo ! iqa

#ifdef MPI
      if(nproc_col.gt.1.and.myrank_row.eq.0.and.myrank_3rd.eq.0) then
        write(6,*) "sum wr from different q-points !"
        mycomm=mycomm_col
        lldim=nlmorb**2
        do ic=1,ncell 
          call mpi_sum_array(0,vmat(:,:,:,ic),lldim,lldim,nsp_crpa,mycomm)
          if(iop_crpa.ge.0) then 
            do iom=1,nomeg 
              call mpi_sum_array(0,umat(:,:,:,iom,ic),lldim,lldim,nsp_crpa,mycomm)
            enddo 
          endif 
        enddo 
      endif
#endif
      if(myrank.eq.0) then 
        call crpa_output 
      endif 

      call end_crpa(0) 
  100 format(8f12.6) 

      contains 

        subroutine sub_set_intercell
        integer:: i,nx,ny,nz
        integer:: ierr

        if(iop_intercell.ge.0) then 
          ncell = 1 + iop_intercell
        elseif(iop_intercell.eq.-1) then 
          ncell = 4
        elseif(iop_intercell.eq.-2) then 
          ncell = 10
        else
          write(6,*) "WARNING: Unsupported value for iop_intercell=",&
     &               iop_intercell
          ncell = 1
        endif 

        allocate(rcell(3,ncell),rcell_cart(3,ncell))
        rcell(:,1) = (/0.0,0.0,0.0/)
        rcell_cart(:,1) = (/0.0,0.0,0.0/)
        
        if(iop_intercell.gt.0) then
          ierr = loct_parse_isdef(blk_cell)
          if(ierr.eq.1) then 
            do i=2,ncell
              call loct_parse_block_float(blk_cell,i-2,0,rcell(1,i))
              call loct_parse_block_float(blk_cell,i-2,1,rcell(2,i))
              call loct_parse_block_float(blk_cell,i-2,2,rcell(3,i))
            enddo
          else
            write(6,*) "WARNING: no input block for ",blk_cell
            ncell = 1
          endif 
        elseif(iop_intercell.eq.-1) then
          rcell(:,2) = (/1.0,0.0,0.0/)
          rcell(:,3) = (/0.0,1.0,0.0/)
          rcell(:,4) = (/0.0,0.0,1.0/)
        elseif(iop_intercell.eq.-2) then
          rcell(:,2) = (/1,0,0/)
          rcell(:,3) = (/0,1,0/)
          rcell(:,4) = (/0,0,1/)
          rcell(:,5) = (/2,0,0/)
          rcell(:,6) = (/0,2,0/)
          rcell(:,7) = (/0,0,2/)
          rcell(:,8) = (/1,1,0/)
          rcell(:,9) = (/0,1,1/)
          rcell(:,10)= (/1,0,1/)
        endif 
        do i=1,ncell 
          write(6,100) rcell(1:3,i) 
        enddo 
 100    format("R=",3f6.2) 

        end subroutine 
      end subroutine task_crpa
!EOC      
