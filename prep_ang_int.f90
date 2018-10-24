!BOP
! !ROUTINE: prep_ang_int
!
! !INTERFACE:
      subroutine prep_ang_int(ml,ngrid)
!
! !DESCRIPTION:
!
! This subroutine sets the grids for the integration of spherical harmonic
! products using Lebedev-Laikov algorithm and calculates the spherical
! harmonics up to ml for those grid points.
! 
! !USES:

      use angintegrals, only: init_angint, n_angrid, sphar, w_angrid
      use lebedev_laikov

! !INPUT PARAMETERS:
!
      implicit none

      integer(4) :: ngrid ! estimated number of gridpoints
!
      integer(4) :: ml    ! maximum l for which the spherical harmonics
!                           are calculated
!
!
! !LOCAL VARIABLES:
      integer(4) :: dimsph ! dimension of the sphar vector (number of 
!                            spherical harmonics) = (ml+1^2)
!
      integer(4) :: ileb   ! index of the grid point
!

      real(8), dimension(3) :: gp     ! stores one gridpoint
!      

      complex(8), allocatable :: sa(:) ! temporary storage for 
!                                                 the spherical harmonics
!
!EOP
!BOC
      dimsph = (ml+1)*(ml+1)
      allocate(sa(dimsph))
      call set_lebedev_laikov_grid(ngrid)

      n_angrid = nleb
      call init_angint(dimsph)
      w_angrid = wleb      
      
      do ileb=1,n_angrid
        gp(1)=xleb(ileb)
        gp(2)=yleb(ileb)
        gp(3)=zleb(ileb)
        call ylm(gp,ml,sa)
        sphar(ileb,1:dimsph)=sa(1:dimsph)
      enddo
      
      deallocate(sa)
      call unset_lebedev_laikov_grid
      return
      end subroutine prep_ang_int
!EOC
