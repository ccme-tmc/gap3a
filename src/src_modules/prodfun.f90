!-----------------------------------------------------------------------
!BOP
!
! !MODULE: prodfun
      module prodfun

! !PUBLIC VARIABLES:      
      
      integer(4) :: nup                    !number of product radial functions
      integer(4), allocatable :: eles(:,:) ! l's of the product function
      real(8), allocatable :: uprod(:,:)   ! radial product function  
      real(8), allocatable :: usprod(:,:)  ! radial product function (2nd. rel. comp)
      real(8), allocatable :: rp(:)        ! radial mesh points
      real(8), allocatable :: umat(:,:)    ! overlap matrix of the product functions
      contains
      
      subroutine init_prodfun(ip)
      
      implicit none
      
      integer(4), intent(in) :: ip
      
      allocate(rp(ip))
      allocate(eles(2,nup))
      allocate(uprod(ip,nup))
      allocate(usprod(ip,nup))
      allocate(umat(nup,nup))
      
      end subroutine init_prodfun
      
      subroutine end_prodfun
      
      deallocate(rp)
      deallocate(eles)
      deallocate(uprod)
      deallocate(usprod)
      deallocate(umat)
 
      end subroutine end_prodfun      
      
      
      end module prodfun
!EOP      
