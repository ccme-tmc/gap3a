!BOP
!
! !ROUTINE: radmesh
!
! !INTERFACE:
      subroutine radmesh(iat,rmesh)

! !DESCRIPTION:
!
! Generates the radial mesh points for atom iat.      
!
! !USES:

      use struk,      only: nrpt,dh,ro

! !INPUT PARAMETERS:

      implicit none
      
      integer(4), intent(in) :: iat ! Index of the inequivalent atom

! !OUTPUT PARAMETERS:      
      
      real(8), intent(out) :: rmesh(*) ! The radial mesh points
      
! !LOCAL VARIABLES:

      integer(4) :: i   ! Counter
            
      real(8) :: dd     ! Logarithmic step
      
! !REVISION HISTORY:
!
! Created 20th. July 2004 by RGA
!
!EOP
!BOC
     
      dd=dexp(dh(iat))
      
      rmesh(1)= ro(iat)
      
      do i=2,nrpt(iat)
        rmesh(i)=rmesh(i-1)*dd
      enddo
      
      end subroutine radmesh
!EOC              
            
