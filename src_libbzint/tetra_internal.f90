!BOP
!
! !MODULE: tetra_internal
!
! !INTERFACE:
      module tetra_internal

      integer :: fout = 10
      logical :: ldbg_bzint = .false.
      integer :: nirkp 
      integer :: ntet                  ! Total number of tetrahedra
      integer :: ncore                 ! Maximum number of core states
      integer :: nband                 ! Maximum number of bands.
      integer :: mndg              
      integer,allocatable:: redtet(:)
      integer,pointer:: tetcorn(:,:)   ! Id. numbers of the corner of the tetrahedra
      integer,pointer:: tetweig(:)     ! Weight of each tetrahedron
      integer,pointer:: qweig(:,:)     ! Weight of each q-point
      integer,pointer:: tetln(:)       ! linked tetrahedron for q.
      integer,pointer:: klinkq(:)      ! linked kpoint for q.
      real(8) :: vt                    ! Relative volume of the tetrahedra      
      real(8),pointer :: eband(:,:)    ! Band energies
      real(8),pointer :: ecore(:)      ! core energies
      real(8):: omgga                  ! the frequency to be included
      integer:: sgnfrq                 ! a sign to tell which weight to be calculated
      real(8):: vol_small_tetra
      real(8), allocatable:: x_gauq(:),w_gauq(:)   ! used for Gauss-Legende quadrature 

      end module tetra_internal
!EOC      
