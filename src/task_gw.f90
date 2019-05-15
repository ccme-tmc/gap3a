!BOP
!
! !ROUTINE: task_gw
!
! !INTERFACE:
      subroutine task_gw 
      
! !DESCRIPTION:
!
! This task subroutine performs "standard" G0W0 and GW0 calculatons with the option that 
! the polarization matrix (eps), M*W*M (mwm) , and exchange selfenergy (selfx) 
! can be calculated in advance and read them from external files. Parallelization is 
! also restructured to allow larger number of processess to be used. These adjustments are
! useful to handle larger systems. 
!
! !USES:

      use bands,       only: nbandsgw,nspin,bande,ibgw,nbgw,eqp,bande0, &
     &                       numin,nomax,nbmax,nbmaxpol,eferqp,efermi
      use barcoul,     only: init_barcoul,end_barcoul,barcevtol,        &
     &                       end_barcev,iop_coul_x,iop_coul_c
      use dielmat,     only: init_dielmat,end_dielmat
      use freq,        only: nomeg
      use kpoints,     only: nkp,nqp,nirkp,kirlist,kirvecs
      use mixbasis,    only: init_mixbasis,end_mixbasis,matsiz
      use mommat,      only: init_mommat,end_mommat
      use selfenergy,  only: init_selfenergy,end_selfenergy,sigx,sigc,  & 
     &                       flag_sxc,fn_selfe,iop_sxc 
      use xcpot,       only: init_xcpot, end_xcpot
      use task,        only: lrestart,casename,nmax_sc
      use modmpi      
      use liboct_parser
      
! !LOCAL VARIABLES:

      implicit none

      integer :: iq0       ! index for starting k-points, needed for the restart mode 
      integer :: iop_mwm   ! control the treatment of mwm  
      integer :: iop_minm  ! control the treatment of mwm  
      integer :: iop_vxc   ! control how to handle vxc matrix 
      integer :: ierr      ! error code 

      character(len=20):: sname="task_gw"
      character(len=20):: blk_gw="gw"

! !REVISION HISTORY:
!
! Created 19.07.2008 by Hong Jiang
!      
!EOP
!BOCa

      call boxmsg(6,'*',"Task:GW") 
!
!     iop_sxc -- determine which approximation to self-energy
!                0 --- GW
!                1 --- exchange-only
!                2 --- static COHSEX
!                3 --- SEX-only 
      iop_minm = -1 
      ierr = loct_parse_isdef(blk_gw)
      if(ierr.ne.1) then 
        iop_sxc = 0 
        iop_vxc = 0 
      else 
        call loct_parse_block_int(blk_gw,0,0,iop_sxc)
        call loct_parse_block_int(blk_gw,0,1,iop_vxc)
      endif 

      call errmsg(iop_sxc.lt.0.and.iop_sxc.gt.3,sname, &
     &      "unsupported iop_sxc value!! ")


      call calcgw0(iop_sxc, 0)

      end subroutine task_gw
!EOC      
