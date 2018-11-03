!BOP
!
! !Function: veclen - zmy
!
! !INTERFACE:

real(8) function veclen(vec)

! !DESCRIPTION:

! Return the length of 3d-vector vec in the Cartisian coordinate
!
! Use modules

! Variables
    IMPLICIT NONE

! Input and Output
    real(8), dimension(3) :: vec

! !REVISION HISTORY:
!
! Created 29. Mar, 2018 by M.-Y. Zhang 
!
! BOC

    veclen = vec(1)**2 + vec(2)**2 + vec(3)**2
    veclen = sqrt(veclen)
    return

end function veclen
! EOC

