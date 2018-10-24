!--------------------------------------------------------------
!BOP
!
! !MODULE: rotmat
      module rotmat

! !PUBLIC VARIABLES:
!
        real(8), pointer :: rotij(:,:,:) ! rotation matrix of atom idf 
!                                          (including equivalent ones) 
!                                          according to the symmetry. 
!                                          For atoms with only one 
!                                          equiv. atom this is always 
!                                          the identity matrix.

        real(8), allocatable :: rotloc(:,:,:) ! local rotation matrix 
!                                               for atom i (from global
!                                               to local coordinate 
!                                               system)

!
! !DESCRIPTION:
!
! Declares and allocates the rotation matrices
!
!EOP
      contains
!BOP
!
! !IROUTINE: init_rotmat
!
! !INTERFACE:      
        subroutine init_rotmat(ndif, nato, i)
        
! !USES:       
          use reallocate, only: doreallocate

          implicit none
          
! !INPUT PARAMETERS:
          
          integer(4) :: ndif ! number of atoms, including equivalent ones
          
          integer(4) :: nato ! number of inequivalent atoms
          
          integer(4) :: i    ! 0 for first allocation, 1 for reallocation
!
!EOP
!
!BOC          
          if(i.eq.0) then
            allocate( rotloc(3,3,nato) )
            allocate( rotij(3,3,ndif) )
            rotloc(1:3,1:3,1:nato)=0.0d0
            rotij(1:3,1:3,1:ndif)=0.0d0
          else
            call doreallocate(rotij, 3,3,ndif)
          end if
        end subroutine init_rotmat
!EOC
      end module rotmat

