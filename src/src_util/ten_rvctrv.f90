!BOP
!
! !Function: veclen - zmy
!
! !INTERFACE:

complex(8) function ten_rvctrv(ndim, tensor, vec)

! !DESCRIPTION:

! calculate vec^{T}\cdot tensor\cdot vec, where tensor is complex and vector is real
!
! Use modules

! Variables
    IMPLICIT NONE

! Input and Output
    integer,intent(in) :: ndim
    real(8),intent(in) :: vec(ndim)
    complex(8),intent(in) :: tensor(ndim,ndim)
    complex(8) :: cvec(ndim)

! Local
    integer :: i,j  ! Counter: 1 to ndim

! !REVISION HISTORY:
!
! Created 05. Nov, 2018 by M.-Y. Zhang 
!
! BOC
    cvec(:) = cmplx(vec(:),0.0D0,8)
    ten_rvctrv = cmplx(0.0D0,0.0D0,8)

    do i=1, ndim
        do j=1, ndim
            ten_rvctrv = ten_rvctrv + tensor(i,j) * cvec(i) * cvec(j)
        enddo
    enddo

    return

end function ten_rvctrv
! EOC

