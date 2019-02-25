!--------------------------------------------------------------
!BOP
!
! !MODULE: struk
      module struk

! !PUBLIC VARIABLES:

        integer(4) :: nat    ! number of inequivalent atoms
        integer(4) :: ndf    ! number of atoms (including equivalent ones)
        real(8) :: neltot    ! total number of electrons in the system
        real(8) :: alat(3)   ! lattice constants (a,b,c)
        real(8) :: alpha(3)           ! angles of the unit-cell's unit-vectors
        real(8), allocatable :: zz(:) ! The charge of each nucleus

        integer(4), allocatable :: iatnr(:) ! atom index of (inequivalent) atom iat also indicates cubic and non-cubic symmetries
                                            ! iatnr(iat)] .GT. 0 => cubic 
                                            ! iatnr(iat) .LT. 0  => non-cubic
        integer(4), pointer :: inddf(:) ! index of atoms, inddf(idf) -- iat corresponding to idf 

        integer(4), allocatable :: mult(:) ! number of equivalent atoms for inequivalent atom iat.

        real(8) :: pia(3)              ! reciprocal lattice constants (2pi/a),(2pi/b),(2pi/c)

        real(8), allocatable :: rmt(:)    ! muffin tin radius of atom iat. 
        real(8) :: rmtmin                 ! minimal muffin tin radius 
        real(8),    allocatable :: vmt(:) ! relative muffin tin sphere volume for atom iat.
        real(8),    allocatable :: dh(:)  ! step size of the logarithmic radial mesh for atom iat; set values in readstruct.f90
        integer(4), allocatable :: nrpt(:)! number of radial  (logarithmic) mesh points for atom iat ; set valuces in readstruct.f90
        real(8),    allocatable :: ro(:)  ! first radial mesh point  for atom iat, set in readstruct.f90

        real(8), pointer :: pos(:,:)   ! position vector of atom idf (including equivalent atoms) in the unit-cell
        real(8) :: vol                 ! volume of the direct lattice unit-cell.
        real(8) :: vi                  ! inverse volume of the direct lattice unit-cell.
        real(8) :: br2(3,3)            ! reciprocal bravais matrix: br2(:,i) -> i-th reciprocal vector 
        real(8) :: rbas(3,3)           ! real space lattice vectors rbas(i,:) --> i-th lattice vector 
        real(8) :: gbas(3,3)           ! recipprocal space lattice vectors, == br2(1:3,1:3)/2\pi 
        logical :: ortho      ! .true. for orthogonal lattice, .false. otherwise
        character(len=4) ::lattic ! lattice type
                                  !  -'P'   => primitive latice (cubic,
                                  !            tetragonal, orthorhombic,
                                  !            monoclinic, triclinoc)
                                  !  -'FC'  => face centered
                                  !  -'BC'  => body centered
                                  !  -'HEX' => hexagonal
                                  !  -'CXY' => c-base centered
                                  !  -'CYZ' => a-base centered
                                  !  -'CXZ' => b-base centered
                                  !  -'R'   => rhombohedral

        character(len=10), allocatable :: atomname(:) ! name of atom iat in the unit-cell (e.g. 'Titanium')

        integer(4) :: nsym                     ! number of symmetry operations of space group
        integer(4), allocatable :: imat(:,:,:) !  matrix representation  of (space group) symmetry operation 
        integer(4), allocatable :: izmat(:,:,:)!  symmetry operation matrix in internal coordinate
        real(8), allocatable :: tau(:,:)       ! non-primitive translation  vector for symmetry  operation j

! variables for rotation matrix 
        real(8), pointer :: rotij(:,:,:) ! rotation matrix of atom idf (including equivalent ones) according to the symmetry.
                                         ! For atoms with only one equiv. atom this is always the identity matrix.
        real(8), allocatable :: rotloc(:,:,:) ! local rotation matrix  for atom i (from global to local coordinate  system)



!
! !DESCRIPTION:
!
! Declares and allocates the variables for the lattice structure
!
!EOP
      contains

!BOP
!
! !IROUTINE: init_struk
!
! !INTERFACE:      
        subroutine init_struk(i)

! !USES:        
          use reallocate, only: doreallocate
        
          implicit none
! 
! !INPUT PARAMETERS:        
          
          integer(4) :: i    ! 0 for first allocation, 1 for reallocation
!EOP
!
!BOC          

          if(i.eq.0) then
            ndf=nat*48
            allocate(atomname(nat))
            allocate( iatnr(nat) )
            allocate( mult(nat) )
            allocate( rmt(nat) )
            allocate( vmt(nat) )
            allocate( nrpt(nat),dh(nat),ro(nat))
            allocate( pos(3, ndf) )
            allocate(zz(nat))
            allocate(inddf(ndf))
            allocate( rotloc(3,3,nat) )
            allocate( rotij(3,3,ndf) )

            rotloc(1:3,1:3,1:nat)=0.0d0
            rotij(1:3,1:3,1:ndf)=0.0d0
            iatnr(1:nat)=0
            mult(1:nat)=0
            rmt(1:nat)=0.0d0
            vmt(1:nat)=0.0d0
            pos(1:3,1:ndf)=0.0d0
            inddf=0 
          else
            call doreallocate(pos,3,ndf)
            call doreallocate(inddf,ndf)
            call doreallocate(rotij,3,3,ndf)
          end if
        end subroutine init_struk
!EOC
      end module struk
