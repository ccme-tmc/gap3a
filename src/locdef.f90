!BOP
!
! !ROUTINE: locdef
!
! !INTERFACE:
      subroutine locdef(rbas,gbas,rotloc)
!
! !INPUT PARAMETERS:

      implicit none

      real(8), intent(in) :: rbas(3,3) ! basis vectors of the real lattice

      real(8), intent(in) :: gbas(3,3) ! basis vectors of the reciprocal lattice
! !INPUT/OUTPUT PARAMETERS:

      real(8), intent(inout) :: rotloc(3,3) ! local rotation matrix


! !DESCRIPTION:
!
!    Redefines local rotation matix from struct-file with
!    unitary transformation  u(-1) * S * u  for non-orthogonal lattice
!
! 
! !LOCAL VARIABLES:

      real(8) :: b(3,3)
! !SYSTEM ROUTINES:
 
      intrinsic matmul

!EOP
!BOC

!      write(6,111) (rotloc(1,j),j=1,3)
!      write(6,111) (rotloc(2,j),j=1,3)
!      write(6,111) (rotloc(3,j),j=1,3)
      b=matmul(gbas,rotloc)
!      write(6,111) (b(1,j),j=1,3)
!      write(6,111) (b(2,j),j=1,3)
!      write(6,111) (b(3,j),j=1,3)
      rotloc=matmul(b,rbas)
!      write(6,111) (rotloc(1,j),j=1,3)
!      write(6,111) (rotloc(2,j),j=1,3)
!      write(6,111) (rotloc(3,j),j=1,3)
!111   format(3f10.4)
      return
      end
!EOC
