!BOP
!
! !ROUTINE: outwin
!
! !INTERFACE:
      subroutine outwin(rel,v,rnot,dh,jri,eh,fl,val,slo,nodes,z)

! !DESCRIPTION:
!
! Integration of the Scalar-relativistic. Schroedinger equation (Rydberg units)
! !USES:
 
      use radwf, only: nrad,a,b
      implicit none

! !INPUT PARAMETERS:

      integer(4) :: jri  ! number of radial grid points
      real(8) :: eh      ! Energy in Hartree 
      real(8) :: fl      ! Angular momentum
      real(8) :: v(nrad) ! rad.sym. Potential in Hartree
      real(8) :: z       ! Nuclear charge
      real(8) :: rnot    ! first radial grid point
      real(8) :: dh      ! log. step length

! !OUTPUT PARAMETERS:

      integer(4) :: nodes ! Number of nodes
      real(8) :: val      ! Value of the wave function at the sphere boundary
      real(8) :: slo      ! Slope of the wave function at the sphere boundary

! !LOCAL VARIABLES:

      integer(4) :: k 
      integer(4) :: iiij

      real(8) :: dg1 
      real(8) :: drdi 
      real(8) :: dg2 
      real(8) :: df1 
      real(8) :: dg3 
      real(8) :: f0 
      real(8) :: sf 
      real(8) :: aa 
      real(8) :: r 
      real(8) :: df2 
      real(8) :: b11 
      real(8) :: det 
      real(8) :: b22 
      real(8) :: rnet(nrad) 
      real(8) :: d(2,3) 
      real(8) :: phi 
      real(8) :: df3
      real(8) :: u 
      real(8) :: y 
      real(8) :: x 
      real(8) :: s 
      real(8) :: e 
      real(8) :: r2 
      real(8) :: r1 
      real(8) :: r3 
      real(8) :: g0
      real(8) :: h83 
      real(8) :: zz 
      real(8) :: c 
      real(8) :: r83sq 
      real(8) :: fllp1
      
      logical    :: rel
!     
!EOP
!
!BOC      
!      dimension        d(2,3),v(nrad),rnet(nrad)
!
!     hartree in ryd
      e=eh*2.d0
!
      do 100 iiij=1,jri
         rnet(iiij)=rnot*(exp(dh*(iiij-1.d0)))
 100  continue
!
      nodes = 0
      zz = z + z
      c = 2.d0 * 137.0359895d0
!
      if(.not.rel) c=1.d+10
!
      fllp1 = fl*(fl + 1.d0)
      r83sq = 64.d0/9.d0
      r1 = 1.d0/9.d0
      r2 = -5.d0*r1
      r3 = 19.d0*r1
      h83 = 8.d0/3.d0
! 
!
      g0 = 1.d0
      if (z .lt. 0.9d0) then
        s = fl+1.d0
        sf = fl
        f0 = fl/c
      else
        aa = zz/c
        s = dsqrt(fllp1 + 1.d0 - aa*aa)
        sf = s
        f0 = g0*(s - 1.d0)/aa
      endif
      do  2  k = 1,3
        r = rnet(k)
        drdi = dh*r
        a(k) = (r**s)*g0
        b(k) = (r**sf)*f0
        d(1,k) = drdi*a(k)*s/r
        d(2,k) = drdi*b(k)*sf/r
    2 continue
!
!
      dg1 = d(1,1)
      dg2 = d(1,2)
      dg3 = d(1,3)
      df1 = d(2,1)
      df2 = d(2,2)
      df3 = d(2,3)
      !write(*, *) dh
      !write(*, *) d
      !write(*, *) v(1:5)
      do  4  k = 4, jri
        r = rnet(k)
        drdi = dh*r
!
!       faktor zwei vor v wegen hartree-rydberg !
!
        phi = (e - 2.d0*v(k)/r)*drdi/c
        u = drdi*c + phi
        x = -drdi/r
        y = -fllp1*x*x/u + phi
        det = r83sq - x*x + u*y
        b11 = a(k-1)*h83 + r1*dg1 + r2*dg2 + r3*dg3
        b22 = b(k-1)*h83 + r1*df1 + r2*df2 + r3*df3
        a(k) = (b11*(h83-x) + b22*u)/det
        b(k) = (b22*(h83+x) - b11*y)/det
        if (a(k)*a(k-1) .lt. 0d0) nodes = nodes + 1
        dg1 = dg2
        dg2 = dg3
        dg3 = u*b(k) - x*a(k)
        df1 = df2
        df2 = df3
        df3 = x*b(k) - y*a(k)
    4 continue
!
!     
      do 120 iiij=1,jri
         b(iiij)=b(iiij)*c/2.d0
 120  continue
!
      val = a(jri)/rnet(jri)
      slo = dg3/(dh*rnet(jri))
      slo = (slo-val)/rnet(jri) 
      return
      end

!EOC
