!EOP
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: freq
      module freq

! !PUBLIC VARIABLES:      
      integer :: iop_freq = 3 ! indicating the approach to treat the frequency dependence       
                              ! 0 - imaginary frequency ( plus analytic continuation for GW)  
                              ! 1 - Godby-Needs Generalized Plasmon Model (GPP)  
                              ! 2 - real freq 
                              ! 3 - imaginary freq 

      integer :: nomeg        ! Number of frequencies for convolution 
      integer :: nomeg_blk    ! number for freq points for each block (needed for emac calculations with large nomeg 

      integer :: iop_fgrid   
!     Flag for the type of frequency mesh
!     iop_fgrid = 1  eqdist  -- equidistant frequencies from 0 to omegmax
!               = 2  gaulag  -- Gauss-Laguerre quadrature from 0 to infinity
!               = 3  gauleg2 -- Gauss-Legendre quadrature from 0 to omegamax and from omegmax to infinity
!               = 4  gauleg  -- Gauss-Legendre quadrature from 0 to omegamax
      
      real(8) :: omegmax     ! cutoff on frequency axis
      real(8) :: omegmin     ! starting point on the freq axis
      real(8) :: tol_difmin   
      real(8) :: freq_eta=0.005 ! the imaginary part on the demoninator 

      real(8), allocatable :: omega(:) ! frequency grid
      real(8), allocatable :: womeg(:) ! weights for frequency grid

      contains 

        complex(8) function freqconvl(iom,nomg,om,enk,mwm,omg,womg) 
!!
!! subroutine for the following frequency convolution 
!!            S(i\omega; \epsilon ) 
!!             = 1/\pi \int_0^\infity \frac{ \epsilon - i \omega ) W(i\omega')} 
!!                  { (\epsilon - i \omega)^2 + \omega'^2} d\omega'
!!    
        implicit none
        integer,intent(in) :: iom   !! control how to do convolution 
                                       !!  iom.eq.0 -- direct integraion 
                                       !!  iom.gt.0 -- indicating that om = omg(iom) 
                                       !!              so that special treatment as explained in Developer Guide (H.4) 
                                       !!              is needed 
        integer,intent(in) :: nomg  !! number of frequency points used for the convolution
        real(8),   intent(in) :: om    !! the frequency \omega 
        real(8),   intent(in) :: enk   !! \epsilon 
        complex(8),   intent(in) :: mwm(nomg)  !! W(i\omega')
        real(8),   intent(in) :: omg(nomg)  !! the array for \omega'
        real(8),   intent(in) :: womg(nomg) !! integraion weight 
        integer :: jom
        complex(8) :: ffac ! the frequecy dependent factor in each term
        complex(8) :: comsq
        complex(8) :: ommek   ! i*omega-e_nk
        complex(8) :: sc
        real(8),parameter:: pi=3.14159265358979d0
 
        ommek=cmplx(enk,-1.0d0*om,8)
        sc=0.d0
        if(iom.gt.0) then 
          do jom = 1, nomg    ! integration frequencies
            if(jom.eq.iom) cycle 
            comsq=cmplx(omg(jom)*omg(jom),0.0d0,8)
            ffac =cmplx(womg(jom),0.0d0,8)/(ommek*ommek+comsq)
            sc=sc+ffac*(mwm(jom)-mwm(iom))
          enddo
          sc=ommek*sc/pi+mwm(iom)*sign(0.5d0,enk)
        else 
          do jom = 1, nomg   
            comsq=cmplx(omg(jom)*omg(jom),0.0d0,8)
            ffac =cmplx(womg(jom),0.0d0,8)/(ommek*ommek+comsq)
            sc=sc+ffac*mwm(jom)
          enddo
          sc=ommek*sc/pi
        endif 
        freqconvl=sc        
      end function 
      
      end module freq
