!BOP
! 
! !ROUTINE: mb_calcrlamint
!
! !INTERFACE:
      subroutine mb_calcrlamint(iat)

! !DESCRIPTION:
!
! This subroutine calculates the set of integrals:
!
!\begin{equation}
!\verb"rtl(iat,i)"\equiv \left\langle r^bl\right\rangle_{aNL} = %
!\int\limits_0^{R^a_{MT}}(r^a)^{bl+2}\upsilon_{aNL}(r^a)dr^a
!\end{equation}
!
!and 
!
!\begin{equation}
!\verb"rint(iat,N,N',bl)"\equiv \left\langle%
!\begin{array}{c}
!  r_<^{bl} \\
!  r_>^{bl +1} \\
!\end{array}%
!\right\rangle_{aNL,N'bl}
!=\iint\limits_0^{R^a_{MT}}%
!\upsilon_{aNL}(r^a_1)\frac{r_<^bl}{r_>^{bl+1}}\upsilon_{aNL}(r_2^a)%
!(r_1^a)^2dr_1^a(r_2^a)^2dr_2^a
!\end{equation}
!
!that will be needed for the calculation of the bare coulomb potential
!
!
! !USES:

      use constants, only: cfein1, cfein2,pi
      use mixbasis,  only: umix, usmix, nmix, bigl,rtl,rrint
      use struk,     only: vi, rmt,nrpt
      use task,      only: fid_outmb
      
! !INPUT PARAMETERS:

      implicit none

      integer(4), intent(in) :: iat ! number of the corresponding 
!                                     inequivalent atom.
!
! !LOCAL VARIABLES:

      integer(4) :: irm  ! Counter: Runs over radial mixed basis functions.
      integer(4) :: jrm  ! Counter: Runs over radial mixed basis functions.
      integer(4) :: ijrm ! Joint index for (irm,jrm) pairs for compressed storage
      integer(4) :: bl   ! Angular momentum quantum number of the mixed basis functions
      integer(4) :: l1   ! just a counters
      integer(4) :: irp  ! Counter: Runs over radial mesh point.
      integer(4) :: jri !  number of radial mesh points, locally stores nrpt(iat)
      real(8) :: rxov  ! value of the integral
      
      real(8), allocatable :: a1(:), b1(:) ! temporary storage for the mixed radial wavefunction

      real(8), allocatable :: ui(:) ! r u_l(r,E) at mesh points.
      real(8), allocatable :: uj(:) ! r us_l(r,E) at mesh points.
      real(8), allocatable :: rr(:) ! Spherical coordinate of the  mesh point m
      real(8), allocatable :: rrtol(:)


 
! !EXTERNAL ROUTINES: 

      external radmesh
      external rint13
!
! !REVISION HISTORY:
! 
! Created Dic. 2003 by RGA
! Last Modification: July. 20th. 2004 by RGA
!
!EOP
!BOC

!
!     Initializations
!
      jri = nrpt(iat)
      allocate(a1(jri))
      allocate(b1(jri))
      allocate(ui(jri))
      allocate(uj(jri))
      allocate(rr(jri))
      allocate(rrtol(jri))

!
!     Set the radial mesh points
!
      call radmesh(iat,rr)
      rrtol(1:jri)=rr(1:jri)
      bl=0
      ijrm=0

!-----------------------------------------------------------------
!                      Calculate rtl
!----------------------------------------------------------------- 
      write(fid_outmb,*)
      write(fid_outmb,101)
      do irm=1,nmix(iat)      
!
!        if bl has changed, calculate r^(bl*1) for each grid point 
!
        if(bigl(irm,iat).gt.bl)then
          do l1=bl+1,bigl(irm,iat)
            do irp=1,jri
              rrtol(irp)=rrtol(irp)*rr(irp)
            enddo
          enddo
          bl=bigl(irm,iat)
        elseif(bigl(irm,iat).lt.bl)then
          write(6,*)'WARNING!!!!'
          write(6,*)'radial mixed funcs not ordered by increasind bl'
          do l1=bigl(irm,iat),bl
            do irp=1,jri
              rrtol(irp)=rrtol(irp)/rr(irp)
            enddo
          enddo
          bl=bigl(irm,iat)  
        endif
!
!       store the corresponding mixed function in a local array
!
        a1(1:jri)= umix(1:jri,irm,iat)         
        b1(1:jri)=usmix(1:jri,irm,iat)         
        call rint13(cfein1,cfein2,a1,b1,rrtol,rrtol,rxov,iat)
        if(rxov.lt.0.0d0)then
          rtl(irm,iat)=-1.0d0*rxov
          umix(1:jri,irm,iat)=-1.0d0*umix(1:jri,irm,iat)
          usmix(1:jri,irm,iat)=-1.0d0*usmix(1:jri,irm,iat)
        else  
          rtl(irm,iat)=rxov
        endif  

        write(fid_outmb,102)irm,rtl(irm,iat)
      enddo

!-----------------------------------------------------------------
!                      Calculate rrint
!----------------------------------------------------------------- 
      rrtol(1:jri)=rr(1:jri)
      bl=0
      ijrm=0
      write(fid_outmb,20)
      do irm=1,nmix(iat)      
!
!        if bl has changed, calculate r^(bl*1) for each grid point 
!
        if(bigl(irm,iat).gt.bl)then
          do l1=bl+1,bigl(irm,iat)
            do irp=1,jri
              rrtol(irp)=rrtol(irp)*rr(irp)
            enddo
          enddo
          bl=bigl(irm,iat)
        elseif(bigl(irm,iat).lt.bl)then
          write(6,*)'WARNING!!!!'
          write(6,*)'radial mixed funcs not ordered by increasind bl'
          do l1=bigl(irm,iat),bl
            do irp=1,jri
              rrtol(irp)=rrtol(irp)/rr(irp)
            enddo
          enddo
          bl=bigl(irm,iat)  
        endif
!
!       store the corresponding mixed function in a local array
!
        do irp=1,jri 
          a1(irp)= umix(irp,irm,iat)         
          b1(irp)=usmix(irp,irm,iat)
          ui(irp)=cfein1*a1(irp)+cfein2*b1(irp)
        enddo  
        uj(1:jri)=ui(1:jri)
        ijrm=irm+(irm*(irm-1))/2
        call drinteg(iat,rxov)
        rrint(ijrm,iat)=rxov
        write(fid_outmb,21)irm,irm,bigl(irm,iat),bigl(irm,iat),rrint(ijrm,iat)
!
!       now loop over jrm for the double integrals      
!
        do jrm=irm+1,nmix(iat)
          ijrm=irm+(jrm*(jrm-1))/2
          if(bigl(irm,iat).eq.bigl(jrm,iat))then
            do irp=1,jri 
              a1(irp)= umix(irp,jrm,iat)         
              b1(irp)=usmix(irp,jrm,iat)
              uj(irp)=cfein1*a1(irp)+cfein2*b1(irp)
            enddo  
            call drinteg(iat,rxov)
          else
            rxov=0.0d0
          endif
          rrint(ijrm,iat)=rxov
          write(fid_outmb,21)irm,jrm,bigl(irm,iat),bigl(jrm,iat),rrint(ijrm,iat)
        enddo 
      enddo         
   20 format(/,5x,'Double radial integrals',/,12x,'N1',2x,'N2',2x,'L1',  &
     &       2x,'L2',5x,'rrint')
   21 format(10x,4i4,g18.10)
  101 format(5x,'Radial integrals',/,13x,'N',3x,'<r^L|v_L>',9x,         &
     &         '<r^(L+2)|v_L> ')
  102 format(10x,i4,2(1pg18.10,1x))
      deallocate(a1,b1,ui,uj,rr,rrtol)      

      
      CONTAINS

      subroutine drinteg(iat,xint)
!
! The procedure used consists in writing it as:
!
!\begin{equation}
!\begin{aligned}
!&\int\limits_{0}^{R_{MT}}\upsilon_{aNL}(r_1)%
!\Biggl[r_1^{-(L+1)}\int\limits_{0}^{r_1}{r_2^{L}%
!\upsilon_{aN'L}(r_2)\left(r_2\right)^2dr_2}-\\
!&-r_1^{L}\int\limits_{R_{MT}}^{r_1}{{r_2}^{-L-1}%
!\upsilon_{aN'L}(r_2)\left(r_2\right)^2dr_2}\Biggr]
!\left(r_1\right)^2dr_1
!\end{aligned}
!\end{equation}
!
! Using an integration subroutine that return the intermediate values of
! the integral, reduces the calculation to just three integrations
!
! !SEE ALSO: fint.f
!
!

      implicit none
      integer(4), intent(in) :: iat ! Inequivalent atom for which the integral is calculated

! !OUTPUT PARAMETERS:
!
      real(8), intent(out):: xint   ! value of the integral
!
! !DEFINED PARAMETERS:

      integer(4), parameter :: nsym = 6 ! Order of the polynom used to fit
!                                         the function and make the
!                                         analytical integration

!
! !LOCAL VARIABLES:
!
      integer(4) :: irp ! runs over spherical grid points

      real(8), allocatable :: rtlint(:) ! 1-d integral <r^l>, the vector contains the intermediate values, that is rtl(irp)=int_0^{r(irp)}u1 r^l u2 dr
      real(8), allocatable :: irtl1(:)  ! the one dimensionan integral
!                                         <r^(-l-1)>, the vector contains the
!                                         intermediate values, that is
!                                         rtl(irp)=int_RMT^{r(irp)}u1 r^(-l-1)
!                                                  u2 dr
      real(8), allocatable :: rr2(:)    ! the grid points in inverse order
      real(8), allocatable :: uitrl(:)  ! u_NL x r^(L+1)
      real(8), allocatable :: ujtrl(:)  ! u_N'L x r^(L+1)
      real(8), allocatable :: uiorl(:)  ! u_NL / r^L
      real(8), allocatable :: ujorl(:)  ! u_N'L / r^L
      real(8), allocatable :: uij(:)


! !EXTERNAL ROUTINES:

      external fint
!
! !REVISION HISTORY:
!
! Created Jan. 2004
! Last Modified: Feb. 2nd. 2004
!
!EOP
!BOC
!
!     Allocate the needed arrays
!
      allocate(rtlint(jri))
      allocate(irtl1(jri))
      allocate(rr2(jri))
      allocate(uitrl(jri))
      allocate(ujtrl(jri))
      allocate(uiorl(jri))
      allocate(ujorl(jri))
      allocate(uij(jri))

      do irp=1,jri
        rr2(irp)=rr(jri-irp+1)
        uitrl(irp)=ui(irp)*rrtol(irp)
        ujtrl(irp)=uj(irp)*rrtol(irp)
        uiorl(jri-irp+1)=ui(irp)*rr(irp)/rrtol(irp)
        ujorl(irp)=uj(irp)*rr(irp)/rrtol(irp)
      enddo
!
!     Integrate in r_2:
!
      call fint(nsym,jri,rr,uitrl,rtlint)
      call fint(nsym,jri,rr2,uiorl,irtl1)
      do irp=1,jri
        uij(irp)= ujorl(irp)*rtlint(irp)-ujtrl(irp)*irtl1(jri-irp+1)
      enddo
!
!      Integrate in r_1:
!
      call fint(nsym,jri,rr,uij,rtlint)
      xint=rtlint(jri)

      deallocate(rtlint)
      deallocate(irtl1)
      deallocate(rr2)
      deallocate(uitrl)
      deallocate(ujtrl)
      deallocate(uiorl)
      deallocate(ujorl)
      deallocate(uij)

      end subroutine drinteg


      end subroutine mb_calcrlamint
!EOC
