!BOP
!
! !ROUTINE: calcudot
!
! !INTERFACE:
      subroutine calcudot(iat,zz,l,el,nodel,nodeu,isp)
      
! !DESCRIPTION:
!
! This subroutine calculates the energy derivative of the radial wave
! functions by finite differences.
! It also insures orthogonality between $u_l$ and %\dot{u}_l$
!
! The result is stored in the global variables \verb"ae" and \verb"be"
!
! !USES:
      
      use lapwlo,    only: umt
      use constants, only: cfein1, cfein2
      use radwf,     only: a, b, ae, be, u, us,vr,rel
      use struk,     only: ro, nrpt, dh
      use task,      only: fid_outmb

! !INPUT PARAMETERS:
      
      implicit none

      integer(4), intent(in) :: iat ! inequivalent atom
      integer(4), intent(in) :: isp ! spin
      integer(4), intent(in) :: l   ! azimutal quantum number of udot
      real(8),    intent(in) :: el  ! linearization energy in Hartrees
      real(8),    intent(in) :: zz  ! charge of the nucleus

! !OUTPUT PARAMETERS:
      
      integer(4), intent(out) :: nodel ! number of nodes of u_l(r,E_l- Delta E)
      integer(4), intent(out) :: nodeu ! number of nodes of u_l(r,E_l+ Delta E)

! !LOCAL VARIABLES:
      
      integer(4) :: npt ! Number of radial mesh points
      integer(4) :: irp ! runs over the spherical grid points
    
      real(8) :: cross  ! overlap between u and udot
      real(8) :: delei  ! Divisor factor used to calculate the energy
!                         derivatives by finite difference: 
!                         delei =  1/(4Delta E)
      real(8) :: e1     ! (temporal) Takes sucesively  the values E_l-Delta E 
!                         and E_l+Delta E. Used to calculate udot_l(r,E_l) 
!                         by finite differences.
      real(8) :: fl     ! Angular momentum l
      real(8) :: dx     ! Step size of the logarithmic radial mesh 
      real(8) :: duvb ! (temporal) value of  d u_l(r,E_l-Delta E)/d r at 
!                       the muffin tin radius
      real(8) :: duve ! (temporal) value of  d u_l(r,E_l+Delta E)/d r at 
!                       the muffin tin radius.
!                       (definitive) value of d\dot{u}_l(r,E_l)/dr at the 
!                       muffin tin radius
      real(8) :: ovlp ! square of the norm of u_l; 
!                       |u_l(r,E_1=E_l-Delta E)|^2
      real(8) :: rnot   ! first radial mesh point (for each atom the 
!                         corresponding ro(iat) is stored here).
      real(8) :: trx ! normalization factor of u_l; 
!                       |u_l(r,E_1=E_l-Delta E)|^-1
      real(8) :: try ! orthogonalization factor between u and udot
      real(8) :: uvb ! (temporal) value of u_l(r,E_l+Delta E) at the 
!                       muffin tin radius.
      real(8) :: uve

! !DEFINED PARAMETERS:
      
      real(8), parameter :: dele = 2.0d-3 ! the up and downward 
!                                           energy-shift in Hartrees

 
! !EXTERNAL ROUTINES: 


      external outwin
      external rint13
      

! !INTRINSIC ROUTINES: 

      
      intrinsic sqrt
      
! !REVISION HISTORY:
!
! Created Feb. 5th. 2004
!
!EOP
!BOC


      delei = .250d+0/dele

      fl = dble(l)
      e1= el-dele
      rnot = ro(iat)
      npt = nrpt(iat)
      dx = dh(iat)
!
!     Calculate u_l(r,E_1=E_l-\Delta E):
!
      call outwin(rel,vr(1,iat,isp),rnot,dx,npt,e1,fl,uvb,duvb,nodel,zz)
!
!     Calculate |u_l(r,E_1=E_l-\Delta E)|^2
!
      call rint13(cfein1,cfein2,a,b,a,b,ovlp,iat)
!
!     Calculate the normalization factor:
!
      trx = 1.0d+0/sqrt(ovlp)
!
!     Store the normalized function u_l(r,E_1) in ae, be:
!
      do irp = 1, npt
        ae(irp) = trx*a(irp)
        be(irp) = trx*b(irp)
      enddo
!
!     Normalize the value and slope of the function at the boundaries:
!
      uvb = trx*uvb
      duvb = trx*duvb
      e1 = el + dele
!
!     Calculate u_l(r,E_1=E_l+\Delta E):
!
      call outwin(rel,vr(1,iat,isp),rnot,dx,npt,e1,fl,uve,duve,nodeu,zz)
!
!     Calculate |u_l(r,E_1=E_l+\Delta E)|^2
!
      call rint13(cfein1,cfein2,a,b,a,b,ovlp,iat)
!
!     Calculate the normalization factor:
!
      trx = 1.0d+0/sqrt(ovlp)
!
!     Calculate the value of udot_l(r,E_l) and d udot_l/dr
!     at the sphere boundaries:
!
      uve = delei * (trx * uve - uvb)
      duve = delei * (trx * duve - duvb)

!
!     Calculate the value of udot_l and d udot_l/dr
!     for the whole mesh:
!
      do irp = 1, npt
        ae(irp) = delei * (trx * a(irp) - ae(irp))
        be(irp) = delei * (trx * b(irp) - be(irp))
      enddo
!
!     ---------------------------
!      Insure ortogonalization
!     ---------------------------
!
!     Calculate <u_l|udot_l>:
!
      call rint13(cfein1,cfein2,u(:,l,iat,isp),us(:,l,iat,isp),ae,be,cross,iat)
      try = - cross
      if(try.lt.-0.05) write(fid_outmb,6070) l, try, ovlp
!
!     Set orthogonalized udot_l:
!
      do irp = 1, npt
        ae(irp) = ae(irp) + try * u(irp,l,iat,isp)
        be(irp) = be(irp) + try * us(irp,l,iat,isp)
      enddo
!
!     Store the orthogonalized values at the boundaries in pe and dpe:
!
      umt(3,l,iat,isp) = uve + try * umt(2,l,iat,isp)
      umt(5,l,iat,isp) = duve + try * umt(4,l,iat,isp)
!
!     Calculate |udot_l(r,E_l)|^2 and store it in  pei(j,iat):
!
      call rint13(cfein1,cfein2,ae,be,ae,be,umt(6,l,iat,isp),iat)

 6070 format (10x,'For l=',i2,' correction=',e14.5,' overlap=',e14.5)

      end subroutine calcudot
!EOC      
