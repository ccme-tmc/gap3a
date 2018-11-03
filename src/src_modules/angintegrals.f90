!--------------------------------------------------------------
!BOP
!
! !MODULE: angintegrals
      module angintegrals
      
! !PUBLIC VARIABLES:

      complex(8), allocatable :: sphar(:,:) !  Value of the spherical harmonic lm (ll=l^2+l+m+1) at the gridpoint  ileb.
      real(8), allocatable :: w_angrid(:) ! weight of the grid point ileb.
      integer(4) :: n_angrid ! number of grip point in the sphere.
      integer(4) :: maxl ! maximum l for which the spherical harmonics  have been calculated.

!
! !DESCRIPTION:
!
! Declares and allocates the variables needed for angular integrals in
! the unit sphere
!
!EOP
      contains
      
      subroutine init_angint(ml)
      
      implicit none
      
      integer(4), intent(in) :: ml
      
      allocate(w_angrid(n_angrid))
      allocate(sphar(n_angrid,ml))
      
      end subroutine init_angint
      
      subroutine end_angint
      
      deallocate(w_angrid)
      deallocate(sphar)
      
      end subroutine end_angint
      end module angintegrals

