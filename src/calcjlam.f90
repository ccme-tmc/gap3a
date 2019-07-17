!BOP
!
! !ROUTINE: calcjlam
!
! !INTERFACE:
      subroutine calcjlam(iat,blmax,jri,iq)
!
! !DESCRIPTION: 
!
!This subroutine calculates the matrix elements
! $< j_{\lambda}^{|\vec{G}+\vec{q}|}>_{aNL}$ according to equation \ref{jlam}
!
!
! !USES:
!
!
      use constants,  only: cfein1, cfein2
      use mixbasis,   only: umix,usmix,nmix,nmixmax,bigl,maxbigl
      use recipvec,   only: gqleng, ngqlen
      use barcoul,    only: jlam
!
!
! !INPUT PARAMETERS:
!
      implicit none

      integer(4), intent(in) :: iat    ! Atom for which the matrix elements are calculated.
      integer(4), intent(in) :: blmax
      integer(4), intent(in) :: jri    ! Number of radial (logarithmic) mesh points for atom iat
      integer(4), intent(in) :: iq


! !LOCAL VARIABLES:
!
      integer(4) :: igl ! (Counter) runs over G-vector lengths.
      integer(4) :: irm ! (Counter) runs over radial mixed basis functions
      integer(4) :: irp ! (Counter), runs over the radial mesh points.
      integer(4) :: l   ! Angular momentum quantum number of the mixed atomic function irm.

      real(8) :: rint  ! Radial integral of umix_{irm}j_{l}
      real(8) :: x     ! Argument of j_l, =|q+G|.r
      real(8) :: qglen ! Length of q+G

      real(8), dimension(jri) :: a1, b1 ! Temporary arrays for the radial functions
      real(8), dimension(jri) :: abes   ! the radial left wave function (1st. rel. component) multiplied by the corresponding spherical bessel function. 
      real(8), dimension(jri) :: rp     ! Radial mesh points

      real(8), dimension(jri,0:blmax) :: bessl   ! the spherical bessel func of order lambda at the mesh point irp.
      real(8), dimension(blmax+1) :: bestemp,dbt ! the spherical bessel func of order l and its derivative at a given mesh point 
 
 
! !EXTERNAL ROUTINES: 
      external k2cart
      external radmesh
      external rint13
      external sphbes

! !INTRINSIC ROUTINES: 


      intrinsic dsqrt

!
! !REVISION HISTORY:
!
! Created: 17th. March 2004 by MF
! Last Modified: 30th. March 2004 by RGA
!
!EOP
!BOC
!      
!     Generate the radial mesh
!      

      call radmesh(iat,rp)
 
      do igl=1,ngqlen(iq)
        qglen=gqleng(igl,iq)
!
!       Calculate the spherical bessel function at each mesh point
!
        do irp = 1, jri
          x=rp(irp)*qglen
!orig     call sphbes(blmax+1,x,bestemp,dbt)
          call sphbes(blmax,x,bestemp,dbt)
          bessl(irp,0:blmax) = bestemp(1:blmax+1)
        enddo   ! irp

        do irm=1,nmix(iat)
          l = bigl(irm,iat)
          a1(1:jri)=umix(1:jri,irm,iat)
          b1(1:jri)=usmix(1:jri,irm,iat)
!
!         Calculate jlam
!
          do irp = 1, jri
            abes(irp) = bessl(irp,l)*rp(irp)
          enddo          ! irp
!  
!         Integrate the wavefunctions:
!
          call rint13(cfein1,cfein2,abes,abes,a1,b1,rint,iat)
          jlam(irm,igl) = rint
        enddo ! irm
      enddo ! igl
      
      end subroutine calcjlam
!EOC
