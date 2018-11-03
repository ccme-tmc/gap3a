!BOP
!
! !ROUTINE: task_gwsc
!
! !INTERFACE:
      subroutine task_gwsc() 
      
! !DESCRIPTION:
!
! This subroutine performs self-consistent gw calculatons at different levels controled by iop_sc
!
! !USES:

      use barcoul,    only: init_barcoul,end_barcoul,barcevtol,end_barcev
      use bands,      only: nbandsgw,nspin,bande,bande0,eqp,ibgw,nbgw,  &
     &                      eferqp,efermi,eferqp0
      use freq,       only: nomeg,omega  
      use kpoints,    only: nkp,nqp,nirkp,kirvecs
      use mixbasis,   only: init_mixbasis,end_mixbasis
      use selfenergy, only: init_selfenergy,end_selfenergy,       &
     &                      sigm,sigm_f,vxcmn,sxcmn,&
     &                      qpwf_coef,qpwf_dmat,vhmn,vh0mn,  &
     &                      iop_dvh,iop_sxc_sc,iop_sxc,iop_sc,  &
     &                      iop_qsgw,beta_sc
      use task,       only: nmax_sc,eps_sc,casename,savdir,lrestart
      use xcpot,      only: init_xcpot, end_xcpot
      use liboct_parser
      use modmpi 
      
! !LOCAL VARIABLES:

      implicit none

      integer :: ierr         ! 
      integer :: iq           ! index for q-points
      integer :: irk          ! index for k-points in IBZ
      integer :: isc          ! index for self-consistent iteration 
      integer :: isc0 ! index for self-consistent iteration 
      integer :: isp          ! index for spin
      integer :: nomeg_bak    ! the number of freq points read from gw.inp

      integer :: isym_sxc  = 1  
      logical :: lconv

      character(len=10)  :: blk_gwsc='gwsc' 
      character(len=10)  :: sname='task_gwsc' 
      character(len=5)   :: flag_eqp="qsgw"
      character(len=10),external::Int2Str
      
      real(8) :: omega1_bak  ! bakup omega(1) when doing self-consistent COHSEX calculations 


! !REVISION HISTORY:
!
! Created 19.07.2008 by Hong Jiang
!      
!EOP
!BOCa
!
!     Read input specific to task_gwsc
!     iop_sxc: choose the approximation to the self-energy 
!          0 -- GW
!          1 -- Exchange only:  
!          2 -- static COHSEX approximation 
!          
!     iop_sc : control the level of self-consistency
!         -1 :  diagonal, w/o self-consistency, i.e. standard G0W0 and e-GW0 
!         -2 :  energy-only self-consistent GW
!          0 -- one-shot calcultions considering off-diagonal contributions without self-consistency
!          1 -- self-consistent GW with fixed W 
!          2 -- full self-consistency 
!     iop_sxc_sc: self-energy approximation used for self-consistency
!          0 -- full GW self-energies 
!          1 -- exchange-only 
!          2 -- COHSEX 
!     iop_dvh : control whether to update Hartree potential during the self-consistency
!          0 - no update, i.e. assuming the density is fixed 
!          1 - update VH 
!     beta_sc : mixing factor used in self-consistency

!     iop_qsgw : control how to hermitize the correlation selfenergy in self-consistent GW calculations 
!          0 - calculate self-energy at the Fermi energy 
!          1 - offdiagonal elements calculated at the Fermi energy but  
!              diagonal elements calculated at the qp energy (mode B in Ref. [Kotani07]) 
!          2 - FSK04 QSGW (mode A in Ref. [Kotani07])

      ierr = loct_parse_isdef(blk_gwsc)
      if(ierr.ne.1) then 
        iop_sxc=0
        iop_sc=-1

        iop_sxc_sc = 0
        iop_dvh = 0
        beta_sc = 0.5 
        nmax_sc = 10
        eps_sc  = 0.001
      else 
        call loct_parse_block_int(blk_gwsc,0,0,iop_sxc)  
        call loct_parse_block_int(blk_gwsc,0,1,iop_sc)   
        if(iop_sc.gt.0) then 
          call loct_parse_block_int(    blk_gwsc,1,0,iop_sxc_sc)  
          call loct_parse_block_int(    blk_gwsc,1,1,iop_dvh)
          call loct_parse_block_float(  blk_gwsc,1,2,beta_sc)
          call loct_parse_block_int(    blk_gwsc,1,3,nmax_sc)   
          call loct_parse_block_float(  blk_gwsc,1,4,eps_sc)   
        else
          iop_sxc_sc = iop_sxc 
          iop_dvh = 0
          beta_sc = 0.5
          nmax_sc = 10
          eps_sc  = 0.001
        endif 

      endif 

      call errmsg(iop_sxc.lt.0.and.iop_sxc.gt.2,sname, &
     &            "unsupported iop_sxc value!!") 

      if(iop_sxc.eq.2) then !! COHSEX 
        nomeg=1
        omega(1)=0.d0
      endif 

      if(iop_sc.eq.-1) then !! G0W0 and energy-only GW0  
        call calcgw0(iop_sxc,0) 
        return 
      elseif(iop_sc.eq.-2) then !! energy-only self-consistent GW
        call calcescgw(iop_sxc) 
        return 
      endif 
      
!
!     when iop_sxc_sc.eq.2, i.e. using COHSEX for self-consistency, only 
!    \omega=0 needs to be considered, so in this case, reset nomeg=1   
!
      if(iop_sxc_sc.eq.2) then 
        nomeg_bak  = nomeg 
        omega1_bak = omega(1) 
        nomeg=1
        omega(1) = 0.d0
      endif 
!
!     Some initilization 
!
      call init_selfenergy(1)          !! initialize selfenergy

      isc0=0
      if(lrestart) call scgw_set_restart('r',isc0) 
!
!     vxcmn 
!
      if(myrank.eq.0.and.isc0.eq.0)then
        call w2k_calcvxcmn(1)
      endif

!
!     Self-consistent loop
!
      do isc=isc0,nmax_sc
        call scgw_check_conv(isc,lconv) 
        if(lconv) exit 

        !! calculate Sxcmn 
        if(iop_sc.eq.0) then 
          call calc_sxcmn(iop_sxc_sc,-1,isym_sxc)
        else 
          call calc_sxcmn(iop_sxc_sc,isc,isym_sxc)
        endif 

        !! calculate new QP energy
        if(myrank.eq.0)then
          call calceqp(iop_sxc_sc,1,isc) 
          call bandanaly(ibgw,nbgw,nirkp,kirvecs,eqp,eferqp,nspin,"QSGW")
        endif 

        call scgw_update_bande

#ifdef MPI
        call mpi_bcast(qpwf_coef,nbandsgw*nbandsgw*nirkp*nspin, &
     &                 mpi_double_complex,0,mpi_comm_world,ierr)

#endif 
        if(iop_sc.eq.0) exit 
        call scgw_set_restart('w',isc)   
      enddo !! isc 
!
!     Some post-sc analysis and output
!
      if(myrank.eq.0) then 
        if(isc.gt.nmax_sc) then 
          write(6,*) "WARNING: fail to converge in task_gwsc after",&
     &             nmax_sc," sc cycles"
        endif 

        call bandanaly(ibgw,nbgw,nirkp,kirvecs, &
     &               bande0(ibgw:nbgw,:,:),efermi,nspin,"KS")

        call bandanaly(ibgw,nbgw,nirkp,kirvecs,&
     &               bande(ibgw:nbgw,:,:),efermi,nspin,"QSGW")
      endif 

!
!     Perform G0W0 calculation based on the preceding self-consistent results
!
      if(iop_sxc.ne.iop_sxc_sc.and.iop_sc.gt.0) then 

        if(iop_sxc_sc.eq.2) then  
          nomeg=nomeg_bak
          omega(1)=omega1_bak
        endif

        if(myrank.eq.0) then
          do isp=1,nspin 
            do irk=1,nirkp
              call trans_vmn('n',sxcmn(:,:,irk,isp), &
     &                       qpwf_coef(:,:,irk,isp),nbandsgw)
            enddo 
          enddo 
        endif

        call update_vector(qpwf_coef)
        bande0 = bande 

        call calcgw0(iop_sxc,1)
      endif 
!
!     write final QP wavefunctions to a vector file 
!
      bande = bande + eferqp0
      call w2k_writevector(flag_eqp,1)
 
      !! cleanup
      call end_selfenergy(1) 
      return
      end subroutine 

