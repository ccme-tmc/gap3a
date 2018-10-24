
!---------------------------------------------------------------
!BOP
!
! !MODULE: fouri
      module fouri

! !DESCRIPTION:
!
!This module declares the R-vectors and stars for the Fourier transform
!

! !PUBLIC VARIABLES:

      integer(4) :: nrr                       ! Number of R vectors
      integer(4) :: nst                       ! Number of start
      integer(4) :: nirk                      ! Number of irreducible k-points
      integer(4) :: irdivk                    ! Division of irred. k-points 
      integer(4), allocatable :: irk(:)       ! integer index of irred. k-points 
      integer(4), allocatable :: rindex(:,:) 
      integer(4), allocatable :: rst(:,:)
      integer(4), allocatable :: kir(:,:)
      real(8) :: rmax
      complex(8), allocatable :: ek(:,:)
      complex(8), allocatable :: er(:,:)
      complex(8), allocatable :: ekpl(:,:)

      logical :: setrindex_done = .false.
      
      end module fouri
!EOP      
!---------------------------------------------------------------
