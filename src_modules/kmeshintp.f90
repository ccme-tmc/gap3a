!--------------------------------------------------------------
!BOP
! !MODULE: kmeshintp
      module kmeshintp
!
! !DESCRIPTION:
!
! Defines the parameters needed for the interpolation of quasi-particle energy between two k-mesh,
! mainly used for 
!           - band structrure plotting 
!           - emac calculations using GW quasi-particle correction 
!
! !PUBLIC VARIABLES:


      integer       :: iop_kip             ! which interpolation method used for k-mesh interpolation 
      character(10) :: eqptag_kip          ! the tag for the eqpH file used for the interpolation 

      integer :: ib0_kip,ib1_kip           ! Number of bands to be interpolated  
      integer :: nvbm(2),ncbm(2)           ! Index for valence band maximum (vbm) and conductance band minimum (cbm) 
      integer :: nsp_kip

      integer :: nbands1
      integer :: nkp1                      ! Number of k-points in the source k-mesh
      integer :: idvk1                     ! Division 
      integer,allocatable:: klist1(:,:)    ! integer coordindates of source k-mesh
      real(8),   allocatable:: kvecs1(:,:)
      real(8),   allocatable:: wk1(:)      ! weight of 2nd k-mesh    

      integer, allocatable :: kk0ind(:,:)

      integer :: nbands2
      integer :: nkp2                      ! Number of k-points in the target k-mesh
      integer :: idvk2                     ! Division 
      integer,allocatable:: klist2(:,:)    ! integer coordindates of target k-mesh
      real(8),   allocatable:: kvecs2(:,:) ! real coordinates 
      real(8),   allocatable:: wk2(:)      ! weight of 2nd k-mesh    

      real(8), allocatable :: eks1(:,:,:)  ! KS enregies on the small k-mesh
      real(8), allocatable :: eqp1(:,:,:)  ! quasiparticle enregies on the small-kmesh

      real(8), allocatable :: eks2(:,:,:)  ! KS enregies on the large k-mesh
      real(8), allocatable :: eqp2(:,:,:)  ! quasiparticle enregies on the large k-mesh

      real(8)::eferks2,eferqp2,eferks1,eferqp1 
      

!EOP
      contains
!BOP
!
! !IROUTINE: init_kmeshintp
!
! !INTERFACE:
      subroutine init_kmeshintp(iop)

! !INPUT PARAMETERS:
        implicit none
        integer,intent(in):: iop

        integer:: ierr
      
!EOP
!BOC
        if(iop.eq.2) then 
          allocate(eks2(nbands2,nkp2,nsp_kip),   &
     &             eqp2(nbands2,nkp2,nsp_kip),   &
     &             klist2(3,nkp2),             &
     &             kvecs2(3,nkp2),             &
     &             wk2(nkp2),                  &
     &             stat=ierr)
        else
          allocate(eks1(ib0_kip:ib1_kip,nkp1,nsp_kip),   &
     &             eqp1(ib0_kip:ib1_kip,nkp1,nsp_kip),   &
     &             klist1(3,nkp1),             &
     &             kvecs1(3,nkp1),             &
     &             wk1(nkp2),                  &
     &             stat=ierr)
        endif 
        call errmsg(ierr.ne.0,"init_kmeshintp","Fail to initialize")
      end subroutine init_kmeshintp

      subroutine end_kmeshintp
        deallocate(eks2,eqp2,klist2,kvecs2,wk2) 
        deallocate(eks1,eqp1,klist1,kvecs1,wk1) 
      end subroutine 

!EOC
      end module kmeshintp
