!BOP
!
! !MODULE: param
      module param
!
! Set the values of parameters that define the maximum allowed size of
! several matrices
!

!
! !PUBLIC DATA MEMBERS:

      integer(4),parameter :: iblock= 128 !.....optimize iblock  for your hardware (32-255)

      integer(4),parameter :: lmax=   13  ! maximum azimutal quantum  number for the muffin-tin  sphere contribution
      integer(4),parameter :: lmmx=  120  ! maxinum number of of collective quantum numbers (l,m) 


      integer(4),parameter :: nkptstart= 100 ! for x-dos set lxdos to 3
      integer(4),parameter :: nmatmax=   5000 
!
! !DESCRIPTION:
!
!EOP

      end module param
