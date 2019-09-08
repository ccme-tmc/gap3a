!BOP
!
! !Function: vecprojlen - zmy
!
! !INTERFACE:

real(8) function vecprojlen(vec, ref, projtype)

! !DESCRIPTION:

! Return the projection of 3d-vector vec in against a reference vector ref
! if projtype=='perp', the distance to the direction of ref will be returned
! instead of the projection

! Variables
  implicit none

! Input and Output
  real(8),dimension(3),intent(in) :: vec ! the vector to be projected
  real(8),dimension(3),intent(in) :: ref ! the vector to be projected against
  character(len=4),intent(in) :: projtype
  real(8) :: refdir(3), vrddot, vecperp(3)

! External functions
  real(8), external :: veclen
    
! !REVISION HISTORY:
!
!  Created 29. Mar. 2018 by M.-Y. Zhang 
! Modified 03. Nov. 2018 by M.-Y. Zhang 
! Modified 09. Sep. 2019 by M.-Y. Zhang 
!
! EOP
! BOC
    refdir(:) = ref(:)/veclen(ref)
    vrddot = vec(1)*refdir(1) + vec(2)*refdir(2) + vec(3)*refdir(3)

  if(projtype=="perp") then
    vecperp(1:3) = vec(1:3) - refdir(1:3)*vrddot
    vecprojlen = veclen(vecperp)
  elseif(projtype=="para") then
    vecprojlen = abs(vrddot)
  else
    write(*,"(A30,A4)") "vecprojlen: bad projtype = ", projtype
    stop
  endif

  return

end function vecprojlen
! EOC

