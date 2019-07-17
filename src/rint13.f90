!BOP
!
! !ROUTINE: rint13
!
! !INTERFACE:
      subroutine rint13(c1,c2,a,b,x,y,s,jatom)

!
! !DESCRIPTION:
!
!     Perform radial integrals:
!
! $\int\limits_0^{R_{MT}}{ u(r) v(r) + \tfrac{1}{c^2} u^s(r) v^s(r) r^2 dr}$
!
!     with c = 137.037 (274.074) in Hartree (Rydberg) units
!
!     For non-relativistic calculations c = $10^{+11}$ Rydberg units.
!
! Author: D.D.Koelling
!
! !USES:

      use struk, only  : nrpt, dh, ro

! !INPUT PARAMETERS:

      implicit none

      real(8), intent(in) :: c1   ! c1 fine structure constant = 1 
      real(8), intent(in) :: c2   ! c2 fine structure constant 
                                  ! = 0 for non-relativistic calculations 
                                  ! ~ 1/137^2 for scalar-relativistic calculations 
      real(8), intent(in) :: a(*) ! r * u (r) on the radial mesh
      real(8), intent(in) :: b(*) ! r * us(r) on the radial mesh
      real(8), intent(in) :: x(*) ! r * v (r) on the radial mesh
      real(8), intent(in) :: y(*) ! r * vs(r) on the radial mesh
!
! !OUTPUT PARAMETERS:

      real(8), intent(out) :: s    ! the value of the radial integral

!EOP
!
!BOC
!      - jatom  the current type of atom
      integer(4) :: jatom
      integer(4) :: j
      integer(4) :: j1
      integer(4) :: jri
      real(8) ::   d 
      real(8) ::  dx 
      real(8) ::  gtem1 
      real(8) ::  gtem2 
      real(8) ::  p1 
      real(8) ::  p2 
      real(8) ::  r 
      real(8) ::  r1 
      real(8) ::  rnot 
      real(8) ::  z2
      real(8) ::   z4
!       
!        intrinsic functions   
!
      intrinsic          exp, mod
!       
      rnot = ro(jatom)
      dx = dh(jatom)
      jri = nrpt(jatom)
      d = exp(dx)
!
      j = 3 - mod(jri,2)
      j1 = j - 1
      r = rnot*(d**(j-1))
      r1 = r/d
      z4 = 0.0d+0
      z2 = 0.0d+0
!
!        loop ... if (j .ge. jri) exit ... end loop
!
!   10 z4=z4+r*(c1*a(j)*x(j)+c2*b(j)*y(j))
   10 continue
         gtem1 = c1*a(j)*x(j)
         gtem2 = c2*b(j)*y(j)
         z4 = z4 + r*(gtem1+gtem2)
         if (z4 .eq. 0.0d+0) z4 = 2.0d-55
         r = r*d
         j = j + 1
         if (j .ge. jri) goto 20
         z2 = z2 + r*(c1*a(j)*x(j) + c2*b(j)*y(j))
         r = r*d
         j = j + 1
!
!        end loop
!
      goto 10
   20 continue
      p1 = rnot*(c1*a(1)*x(1) + c2*b(1)*y(1))
      p2 = r1*(c1*a(j1)*x(j1) + c2*b(j1)*y(j1))
      s = 2*z2 + 4*z4 + r*(c1*a(j)*x(j) + c2*b(j)*y(j)) + p2
      s = (dx*s+p1)/3.0d+0
      if (j1 .gt. 1) s = s + 0.5d+0*dx*(p1+p2)
!
      return
!
!        end of 'rint13'
!
      end subroutine rint13
!EOC
