!BOP
!
! !ROUTINE: task_gw2wann
!
! !INTERFACE:
      subroutine task_gw2wann 
      
! !DESCRIPTION:
!
! This task calculate the full matrices for Vxc and Sxc in the Kohn-Sham basis, and
! then transfrom them to the Wannier basis. It can be used for the
! subsequent GW+DMFT calculation or the Wannier interpolation of GW band
! structure. 
!
! !USES:

      use bands,       only: nbandsgw,nspin,bande,ibgw,nbgw,eqp,bande0, &
     &                       numin,nomax,nbmax,nbmaxpol,eferqp,efermi, &
     &                       bande,eferks
      use barcoul,     only: init_barcoul,end_barcoul,barcevtol,        &
     &                       end_barcev,iop_coul_x,iop_coul_c
      use crpa,        only: iop_wf,nlmorb,pln
      use dielmat,     only: init_dielmat,end_dielmat
      use freq,        only: nomeg
      use kpoints,     only: nkp,nqp,nirkp,get_kvec
      use mixbasis,    only: init_mixbasis,end_mixbasis,matsiz
      use mommat,      only: init_mommat,end_mommat
      use selfenergy,  only: init_selfenergy,end_selfenergy,sigx,sigc,  & 
     &                       flag_sxc,fn_selfe,iop_sxc,isym_kbz,hks_wann,&
     &                       sigm,sxc_wann,vxcmn,vxc_wann
   
      use xcpot,       only: init_xcpot, end_xcpot,vxcnn
      use task,        only: lrestart,casename,nmax_sc
      use modmpi      
      use liboct_parser
      
! !LOCAL VARIABLES:

      implicit none

      integer :: iq        ! index for q-points
      integer :: iq0       ! index for starting k-points, needed for the restart mode 
      integer :: iop_mwm   ! control the treatment of mwm  
      integer :: iop_minm  ! control the treatment of mwm  
      integer :: ierr      ! error code 
      integer :: nktot     ! the total number of k-points considered, and it can be nirkp or nkp depending on isym_kbz

!     variables used for parallelization
      integer :: iq_f,iq_l     !! lower and upper bound for iq
      logical:: lread_eps=.false.

      character(len=20):: sname="task_gw2wann"
      character(len=20):: blk_gw2wann="gw2wann"

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
!     isym_bz = 0 - calculate Vxc, Sxc in the full BZ
!             = 1 - calculate Vxc, Sxc in the irreducible BZ

! iop_wf ---> option about the Wannier function 
!        0 -- use the internally generated WF (not implemented yet!)
!        1 -- use WF generated by dmftproj 
!        2 -- use WF generated by wien2wannier + wannier90 (not
!        implemented yet!)


      iop_minm = -1 
      ierr = loct_parse_isdef(blk_gw2wann)
      if(ierr.ne.1) then 
        iop_sxc = 0 
        iop_wf  = 2
        isym_kbz = 0 
      else 
        call loct_parse_block_int(blk_gw2wann,0,0,iop_sxc)
        call loct_parse_block_int(blk_gw2wann,0,1,iop_wf)
        call loct_parse_block_int(blk_gw2wann,0,2,isym_kbz)
      endif 

      call errmsg(iop_sxc.lt.0.and.iop_sxc.gt.2,sname, &
     &      "unsupported iop_sxc value!! ")

      if(isym_kbz.eq.0) then 
        nktot = nkp
      else 
        nktot = nirkp
      endif 

     !! Read Wannier projectors from case.pln
      if(iop_wf.eq.1) then   !! read WF projectors from dmftproj.x 
        call crpa_readpln
      elseif(iop_wf.eq.2) then !! read WF from a Wannier90 chk file 
        call crpa_readw2w
      else
        write(6,*) "ERROR: not implemented iop_wf!"
      endif

      call init_selfenergy(2) 

      if(myrank.eq.0) then
        call init_xcpot
        call w2k_calcvxcmn(isym_kbz) 
      endif

      call calc_sxcmn(iop_sxc,-1,isym_kbz) 

      if(myrank.eq.0)then
        !! transform Vxc_{n,n'} and Sxc_{n,n'} to the matrices in thea Wannier basis 
        call sub_trans_ks2wann
      endif

      call io_sxc_wann('w',isym_kbz,ierr) 

      call end_selfenergy(2) 
      if(myrank.eq.0) call end_xcpot

      return

      contains 

        subroutine sub_trans_ks2wann()
        implicit none 
        integer:: isp,iik,ik,irk,iw,jw,ie
        complex(8),allocatable:: work(:,:)

        allocate(work(nlmorb,nbandsgw)) 

        hks_wann=0.d0
        do isp=1,nspin
          do iik=1,nktot
            call get_kvec(isym_kbz,iik,ik,irk) 
            call sub_ks2wann(sigm(:,:,:,iik,isp), & 
     &       sxc_wann(:,:,:,iik,isp),pln(:,:,iik,isp),work,nbandsgw,nlmorb,&
     &       nomeg+1)
            call sub_ks2wann(vxcmn(:,:,iik,isp), &
     &           vxc_wann(:,:,iik,isp),pln(:,:,iik,isp),work,nbandsgw,nlmorb,1) 
            call sub_ks2wann(vxcmn(:,:,iik,isp), &
     &       vxc_wann(:,:,iik,isp),pln(:,:,iik,isp),work,nbandsgw,nlmorb,1)
            do jw=1,nlmorb
              do iw=1,nlmorb
                do ie=ibgw,nbgw
                  hks_wann(iw,jw,iik,isp) = hks_wann(iw,jw,iik,isp)  &
     &             + pln(iw,ie,iik,isp)*(bande(ie,irk,isp)+eferks)             & 
     &             *conjg(pln(jw,ie,iik,isp))
                enddo
              enddo
            enddo
            call sub_solve_eig(hks_wann(:,:,iik,isp),nlmorb)
          enddo
        enddo
        deallocate(work) 
        end subroutine 
!
!-----------------------------------------------------------------------
!  internal utility subroutine to transform a matrix in the KS basis   !
!  to the  Wannier basis                                               !
!-----------------------------------------------------------------------                            
        subroutine sub_ks2wann(mat_k,mat_w,pln,work,nks,nwf,nt) 
        implicit none
        integer:: nks,nwf,nt
        complex(8):: mat_k(nks,nks,nt),mat_w(nwf,nwf,nt),pln(nwf,nks)
        complex(8):: work(nwf,nks)

        integer:: i
        complex(8):: zone=1.0,czero=0.0  
        do i=1,nt
          call zgemm('n','n',nwf,nks,nks,zone,pln,nwf,mat_k(:,:,i),nks,&
     &             czero,work,nwf) 
          call zgemm('n','c',nwf,nwf,nks,zone,work,nwf,pln,nwf,czero,  & 
     &             mat_w(:,:,i),nwf)
        enddo   
        end subroutine

        subroutine sub_solve_eig(hmat,n)
        integer,intent(in):: n
        complex(8):: hmat(n,n) 
        integer:: i,ierr
        integer:: lwork,rwsize
        real(8),allocatable::rwork(:),ev(:)
        complex(8),allocatable::work(:)

        lwork=2*n
        rwsize=3*n
        allocate(work(lwork),rwork(rwsize),ev(n))
        write(6,*) "Diagonalize H matrix in the Wannier basis"
        call zheev('v','u',n,hmat,n,ev,work,lwork,rwork,ierr)
        call errmsg(ierr.ne.0,sname,'Fail to diag. Hmat by zheev')
        do i=1,n
          write(6,'(i4,f12.6)') i,ev(i) 
        enddo 
        deallocate(work,rwork,ev) 
        end subroutine  

      end subroutine task_gw2wann
!EOC      
