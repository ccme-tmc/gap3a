!BOP
!
! !ROUTINE: calcvhmn
!
! !INTERFACE: 
      subroutine calcvhmn(isc)
      
! !DESCRIPTION:
!
!This subroutine calculate the full exchange self-energy matrix elements
!

! !USES:

      use bands,      only: bande,nbmax,nspin,fspin,    &
     &                      nomaxs,numins,ibgw,nbgw,nbandsgw
      use barcoul,    only: barcev,barcevsq
      use bzinteg,    only: singc2,kwt_bz
      use constants,  only: cone,czero,pi,fourpi,sqrt4pi
      use core,       only: corind, eigcore,ncg
      use kpoints,    only: nirkp,nkp,kqid,weightq,idikp,     &
     &                      kpirind
      use mixbasis,   only: matsiz
      use selfenergy, only: vh0mn,vhmn,qpwf_dmat
      use struk,      only: vi
      use task,       only: time_lapack

      use modmpi

! !INPUT PARAMETERS:
      implicit none
      
      integer(4), intent(in) :: isc  ! indicate the stage of self-consistency, which in the meanwhile have partial 
                                     ! effects on the treatment of Minm   
                                     !    0 -- the first step during the self-consistency, Minm are calculated and saved
                                     !  > 0 -- continue self-consistency, Minm read from scratch space 

! !REVISION HISTORY:
!
! Created Nov. 05, 2009 by H. Jiang 
!
!EOP

!BOC
! !LOCAL VARIABLES:            

      integer(4) :: ie1m,ie1n,ie1,ie2,ie     ! (Counter) Runs over bands.
      integer(4) :: ik,irk  ! index for k-point related to n
      integer(4) :: jk,jrk  ! index for k-point related to m  
      integer(4) :: fid_minm
      integer(4) :: nomx,numn,nbands,nmdim
      integer(4) :: ierr
      integer(4) :: isp     ! Index for spin 
      real(8) ::    time1,time2
      character(len=10)::sname='calcvhmn'
      logical ::    ldbg=.false.

      integer :: mcounts
      complex(8) :: cwt
      complex(8), allocatable:: minm(:,:,:)   ! M_{nm}^i
      complex(8), allocatable:: minn(:)   ! Minn, needed to calculate VH 
      complex(8), allocatable:: mh(:)   ! Minn, needed to calculate VH 

! !INTRINSIC ROUTINES: 
      intrinsic abs
      intrinsic sign
      intrinsic cmplx
      intrinsic conjg

! !EXTERNAL ROUTINES: 
      complex(8), external :: zdotc,zdotu
      external zhemm

      if(ldbg) call linmsg(6,'-','calcvhmn')

      nmdim=(nbgw-ibgw+1)**2
      allocate(minm(matsiz,ibgw:nbgw,ibgw:nbgw), &
     &         minn(matsiz),                     &
     &         mh(matsiz),stat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to allocate minm") 
     
!
!     calculate M_i
! 
      mh=0.d0
      do isp=1,nspin 
        nomx=nomaxs(isp) 

        do ik=1,nkp
          irk=kpirind(ik)
          call cpu_time(time1)
          call readvector(ik,1,isp,0)          !! Read the eigenvector corresponding to the k-point ik
          call expand_evec(ik,1,.true.,isp)  !! Calculate the expansion coeficients of the eigenvectors
          call readvector(ik,2,isp,0)          !! Read the eigenvector corresponding to the k'-point jk
          call expand_evec(ik,2,.true.,isp)  !! Calculate the expansion coeficients of the eigenvectors
          call cpu_time(time2)
       
          cwt=fspin*cone*kwt_bz(ik) 

          call calcminm(ik,1,ibgw,nbgw,ibgw,nbgw,isp,minm)
          call zgemv('n',matsiz,nmdim,cwt,minm,matsiz,     &
     &              qpwf_dmat(:,:,irk,isp),1,cone,mh,1)

        enddo  ! ik
      enddo !! isp

      do isp=1,nspin 
        do irk=1, nirkp
          ik=idikp(irk)  

          call readvector(ik,1,isp,0)          !! Read the eigenvector corresponding to the k-point ik
          call expand_evec(ik,1,.true.,isp)  !! Calculate the expansion coeficients of the eigenvectors
          call readvector(ik,2,isp,0)          !! Read the eigenvector corresponding to the k'-point jk
          call expand_evec(ik,2,.true.,isp)  !! Calculate the expansion coeficients of the eigenvectors

          call calcminm(ik,1,ibgw,nbgw,ibgw,nbgw,isp,minm)
       
          call cpu_time(time1)
          call zgemv('c',matsiz,nmdim,cone,minm,matsiz,mh,1,czero,&
     &              vhmn(:,:,irk,isp),1) 
          call cpu_time(time2)
          time_lapack = time_lapack + time2 - time1

          if(isc.eq.0) then 
            vh0mn(:,:,irk,isp)=vhmn(:,:,irk,isp)
          endif 
        enddo ! irk
      enddo ! isp

      deallocate(minm,mh,minn)

      return
              
      end subroutine calcvhmn
!EOC        
