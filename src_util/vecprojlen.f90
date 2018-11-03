!BOP
!
! !Function: vecprojlen - zmy
!
! !INTERFACE:

real(8) function vecprojlen(vec, ref, projtype)

! !DESCRIPTION:

! Return the length of the projection of 3d-vector vec in the Cartisian coordinate against
! a reference vector ref
! if projtype == 'perp', instead of the projection, the distance to the direction of ref
! will be returned
!
! Use modules

! Variables
    IMPLICIT NONE

! Input and Output
    real(8), dimension(3),intent(in) :: vec ! the vector to be projected
    real(8), dimension(3),intent(in) :: ref ! the vector to be projected against
    character(len=4),intent(in) :: projtype
    real(8) :: reflen, vrdot, projlen, vecperp(3)

! External functions
    real(8), external :: veclen
    
! !REVISION HISTORY:
!
! Created 29. Mar. 2018 by M.-Y. Zhang 
!         03. Nov. 2018 by M.-Y. Zhang 
!
! EOP
! BOC
    reflen = veclen(ref)

    vrdot = vec(1)*ref(1) + vec(2)*ref(2) + vec(3)*ref(3)
    projlen = abs(vrdot)/reflen

    if(projtype=="perp") then
        vecperp(1:3) = vec(1:3) - ref(1:3)*projlen
        vecprojlen = veclen(vecperp)
    else
        vecprojlen = projlen
    endif

end function vecprojlen
! EOC

