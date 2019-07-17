!BOP
!
! !ROUTINE: rotate
!
! !INTERFACE:
      subroutine unrotate(vector,rotmat,rotvec)

! !DESCRIPTION:
!
! Performs a rotation of the vector from the general cartesian coordination
! system into the local one of the jatom-th sphere.
!

! !INPUT PARAMETERS:

      implicit none
      
      real(8), intent(in) ::  vector(3)
      
      real(8), intent(in) ::  rotmat(3,3)

! !OUTPUT PARAMETERS:      
      
      real(8), intent(out) ::  rotvec(3)

! !LOCAL VARIABLES:

      integer(4) :: j, jcoord
      
      real(8)    ::   dotpro
!EOP
!BOC
      do jcoord = 1, 3
        dotpro = 0.0d0
        do j = 1, 3
          dotpro = dotpro + vector(j)*rotmat(j,jcoord)
        enddo  
        rotvec(jcoord) = dotpro
      enddo
      return
      end subroutine unrotate
!EOC      
