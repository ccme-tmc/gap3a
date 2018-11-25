!------------------------------------------------------------------------
!BOP
!
! !MODULE: constants
      module constants
!
! !DESCRIPTION:
!
! This module declares useful universal constants that can be used throughout the
! program

! !PUBLIC VARIABLES:
        complex(8), parameter :: imag = (0.0d+0,1.0d+0) ! Imaginary unit.
        complex(8), parameter :: czero = 0.0d+0 !  complex zero.
        complex(8), parameter :: cone = (1.0d0,0.0d0)
        real(8), parameter :: pi = 3.14159265358979D0
        real(8), parameter :: twopi = pi*2.0D0
        real(8), parameter :: fourpi=pi*4.0D0
        real(8), parameter :: sqrt4pi=3.544907701811032055D0
        real(8), parameter :: sqrt3=1.73205080757d0
        real(8), parameter :: sqrt2=1.41421356237d0
        real(8), parameter :: pisq= pi*pi
        real(8), parameter :: pisqrt= 1.772453850905516104D0
        real(8), parameter :: HeV = 27.2113961D0
        real(8), parameter :: Bohr2Ans=0.5292083D0
        complex(8), parameter :: ctwopi = cmplx(twopi,0.0D0,8)
        complex(8), parameter :: cpi = cmplx(pi,0.0D0,8)
        
        real(8) :: cfein1 ! 1st. fine structure constant (always equal to 1.)
        real(8) :: cfein2 ! 2nd. fine structure constant 
                          !  = 137.0359895^{-2} for scalar realtivistic calculations
                          !  = 10.^{-22} for non relativistic calculations

        integer(4)::size_cmplx = 16  ! number of bits for each complex number
        integer(4)::size_real = 8    ! number of bits for each real number 
        integer(4)::size_int = 4
!
!EOP
      end module constants

