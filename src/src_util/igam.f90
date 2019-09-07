!BOP
! !ROUTINE: igam
! 
! !INTERFACE: 
        recursive function igam(n) result(tgam)
! !DESCRIPTION:
!
! This function calculates the gamma function for integer values
!
! !INPUT PARAMETERS:        
        implicit none
        integer(4) :: n
! !OUTPUT PARAMETERS:        
        real(8) :: tgam
! !DEFINED PARAMETERS:
! !REVISION HISTORY:
!
!   Created 6, Sept 2019 (MYZ)
!
!EOP
!BOC
        if(n.le.0)then
          write(*,*)"error: igam: n must be positive"
          stop
        elseif(n.le.2)then
          tgam = 1.0d0
        else
          tgam= real(n-1,8)*igam(n-1)
        endif
        end function igam
!EOC
