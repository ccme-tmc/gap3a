!
! the main module used as the interface to the libbzint 
! BZ integrations
!
      module bzint 
      implicit none 
      integer:: iop_bzint=0      !! option to perform the standard BZ integration 
                                !!  -1 - the tetrahedron method, weights calculated from all tetrahdra (test purpose only)
                                !!   0 - the tetrahedron method, weights directly calulated from irreducible tetrahedra
                                !!   1 - use the Fermi smearing 

      integer:: iop_bzintq=0    !! option to perform q-dependent BZ intergation (convolution) 
                                !!   0 - standard analytic integration (old) 
                                !!   1 -  numerical integration of weights

      integer:: nomg_ht=1000    !! the number of freq points used for the Hilbert transform (HT) 
                                !! approach to calculate q-BZ weights
      real(8):: omgmax_ht=4.0   !! the upper bound for the HT integration  

      real(8):: esmear = -1.0   !! the energy corresponding to the broadening temperature  in 

      real(8):: eta_freq=0.01   ! this is the value of the small imaginary part that 
                                ! is needed for real frequency calculations 


      integer :: iop_bcor  = 0         ! the option to choose whether to use Bloechl correction 
                                       !  0 -- w/o  Bloechl correciton 
                                       !  1 -- with Bloechl correction  

      integer :: n_gauq = 8            ! number of points used for Gauss  quaduarture 

      real(8) :: tol_unphys_weight = 1.0e+4   ! the tolerance for unphysical weights, 
                                               ! when a weight is found to be larger than this,
                                               ! a warning is printed 
      real(8) :: ztol_vol = 1.e-10              ! tolerance for zero volume
      real(8) :: tol_taylor = 10.0              ! the tolerance for the use of Taylor expansion 
      real(8) :: ztol_sorteq = 1.0*1.e-2

      integer :: nirtet  ! Number of irreducible tetrahedra
      integer :: ntet    ! Number of tetrahedra
      integer, allocatable :: tnodes(:,:)  ! index of the k-points corresponding to the nodes of each tetrahedra for convolution
      integer, allocatable :: tndi(:,:)    ! index of the k-points corresponding to the nodes of each tetrahedra for integration
      integer, allocatable :: link(:,:)    ! index of the tetrahedra linked to by the corresponding q vector
      integer, allocatable :: wtet(:)      ! geometrical weight of each tetrahedron  for integration
      integer, allocatable :: wirtet(:)    ! geometrical weight of each irred. tetrahedron  for integration
      real(8) :: tvol                      ! volume of the tetrahedra relative to the BZ volume
      target tnodes, tndi, link, wtet, wirtet

      contains 

!
! Frequency related factor 
!
        complex(8) function freq_factor(iop,omg,edif)
        implicit none
        integer, intent(in):: iop   ! 2 for refreq and 3 for imfreq
        real(8), intent(in):: omg, edif

        if(iop.eq.2) then
          freq_factor =  1.d0/cmplx(omg - edif,   eta_freq)  &
     &                 - 1.d0/cmplx(omg + edif, - eta_freq)
        elseif(iop.eq.3) then
          freq_factor = cmplx( -2.d0*edif/(omg*omg+edif*edif), 0.d0)
        endif
        end function

!
! Occupation 
!
        real(8) function bzint_smear(iop,e)
        implicit none
        integer, intent(in):: iop   ! 0 for Fermi, 1 for Gauss
        real(8), intent(in):: e 
        real(8) :: wt

        if(iop.eq.0) then

          if(esmear.le.0.d0) then
            if(e.gt.0.d0) then
              wt = 0.d0
            else
              wt = 1.d0
            endif
          else
            if(e/esmear.gt.10.d0) then
              wt = 0.d0
            elseif(e/esmear.lt.-10.d0) then
              wt = 1.d0
            else
              wt = FermiF(e/esmear)
            endif
          endif

        elseif(iop.eq.1) then
          if(abs(e)>10.d0*esmear) then
            wt = 0.d0
          else
            wt = - FermiF_d(e/esmear)/esmear 
          endif
 
        elseif(iop.eq.2) then 
          wt = GaussF(e,esmear) 
        endif
        bzint_smear = wt
        end function

!
!       Fermi function
!
        real(8) function FermiF(x)
        real(8):: x
        FermiF = 1.d0/(1.d0+exp(x))
        end function

!
!       the first derivate of the Fermi function 
!
        real(8) function FermiF_d(x) 
        real(8) :: x
        FermiF_d = -1.d0*exp(x)/(1.d0+exp(x))**2
        end function 
!
!       Gauss function
!
        real(8) function GaussF(x,s)
        real(8),intent(in):: x,s
        real(8):: sr2pi = 2.5066282746310002

        GaussF = 1.d0/(s*sr2pi)*exp(-x*x/(2.d0*s*s))
        end function

      end module 

