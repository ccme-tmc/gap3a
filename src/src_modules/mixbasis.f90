!------------------------------------------------------------------------
!BOP
!
! !MODULE: mixbasis
      module mixbasis
      use recipvec, only: ngq,ngqbarc,indgq
      use reallocate

! !PUBLIC VARIABLES:

      integer :: nspin_mb=1 ! option for constructing mixed basis in the spin-polarized case (nspin=2) 
                              ! 2 -- use both spin up and down radial wave functions  
                              ! 1 -- use only spin up radial wave functions     
      integer :: iopq0     ! option for generating mixbasis at q=0
                              ! 0 --  
                              ! 1 -- 
      logical :: mb_ludot = .false.  ! control whether considering \dot{u}_l when setting up the mixed basis 
      logical :: mb_lcore =.true.   ! control whether considering core states when building the mixed basis 

      integer :: lmbmax=2    ! maximum l of radial functions considered when 
                             ! building the mixed basis set
      integer :: lblmax=0    ! the maximum L for the mixed basis (big l)

      integer,allocatable :: lmbmax_at(:) ! atom-specific lmbmax  
      integer,allocatable :: lblmax_at(:) ! atom-specific lblmax  

      integer :: mbsiz       ! total number of original mixed basis functions
      integer :: matsiz      ! the number of eigenvectors of bare Coulomb matrix used as basis functions
      integer :: matsiz_x    ! the number of eigenvectors of bare Coulomb matrix used as basis functions
      integer :: matsiz_c  ! the number of eigenvectors of bare Coulomb matrix used as basis functions
      integer :: locmatsiz ! Total number of atomic-like mixed  basis functions
      integer :: nmixmax   ! maximum number of mixed basis functions (nl) per atom: nmixmax = max[nmix(iat)]
      integer :: lmixmax   ! maximun number of mixed functions (nlm) per atom
      integer :: lleftmax  ! Maximum L of the left product functions
      integer :: maxbigl   ! Maximum L of the mixed basis functions
      integer :: nlo_mb    ! the maximal number of LO's per l considered in mixed basis 

      real(8) :: kmr       ! maximum G for the mixed basis 
      real(8) :: pwm       ! maximum G for the pw basis (in kmax*kmr units)
      
      integer, allocatable :: lmax_at(:) ! maximum l of occupied orbitals for each atom 
      integer, allocatable :: mbl(:) !  Maximum L of the mixed basis functions for atom iat
      integer, pointer :: bigl(:,:)  ! L of the radial wave function i of atom iat 
      real(8) :: wftol                  ! nonorthogonality tolerance.


      real(8) :: mb_emin = -1.0E10   ! the energy cut-off for mixed basis funcitons, any radial functions  
                                     ! with energy lower than mb_emin are dropped in constructing mixbasis 
      real(8) :: mb_emax = 10.0 
                                     

      real(8), pointer :: rtl(:,:)        ! the matrix elements <r^lambda>_{al,l'}
      real(8), pointer :: rrint(:,:)       ! the matrix elements   <r_<^{lambda}/r_>^{lambda +1}>_{al,l';l_1,l'_1}


      integer, allocatable :: indggq(:)
      integer, allocatable :: nmix(:)         ! number of mixed basis functions for atom iat including only (nl) 
      integer, allocatable :: nmixlm(:)       ! number of mixed basis functions for atom iat including (nlm) 
      integer, allocatable :: locmixind(:,:)  ! index of mixed basis function at each atom 


      integer, dimension(3) :: igm      ! maximum number of G vectors for the theta coeficient expansion of the IPW 

      real(8), allocatable :: cgcoef(:)    ! The Gaunt coefs.
       
      real(8), pointer :: umix(:,:,:)      ! The radial mixed basis functions
      real(8), pointer :: usmix(:,:,:)     ! The radial mixed basis functions (2nd. rel comp)

      complex(8), allocatable :: mxov(:,:) ! the overlap matrix of the spherical product functionsa

      integer :: nrwmax                    ! maximum number of all radial wavefunctions per atom
      real(8),pointer :: s3r(:,:,:,:,:)    ! overlap integrals of three radial wavefunctions 
                                           ! the matrix elements <ll'lambda|lambda'> (Eq.C1.5)

!
! Variables related to the mixed basis functions in the interstitial region 
!
      complex(8), allocatable :: sgi(:,:)    ! transformation matrix from IPW to its orthogonalized form 
      complex(8), allocatable :: ipwint(:)   ! the integral of a plane waves over the interstitial region
      complex(8), allocatable :: mpwipw(:,:) ! the tranformation matrix between mixed basis and plane waves
      complex(8), allocatable :: mpwmix(:,:) ! the tranformation matrix between mixed basis and  plane waves
      complex(8), allocatable :: wi0(:)      


!!private variables
      integer,private :: ierr
     

! !DESCRIPTION:
!
! This module declares and allocates the variabels needed to set up the
! mixed basis, and it's overlap with the LAPW orbitals          
!
!EOP
      contains
      

      subroutine init_mixbasis(iq)
      implicit none 
      integer, intent(in) :: iq

      integer :: ig,igmax

      mbsiz=locmatsiz+ngq(iq)
      igmax=maxval(indgq(:,iq))
      write(6,101)iq,locmatsiz,ngq(iq),mbsiz

      call flushbuf(6) 
      allocate(mpwipw(ngq(iq),ngqbarc(iq)), &
     &         sgi(1:ngq(iq),1:ngq(iq)),    &
     &         wi0(mbsiz),                 &
     &         indggq(1:igmax),             &
     &         stat=ierr)
      if(ierr.ne.0) then
        write(6,*) "init_mpwipw: Fail to allocate mpwipw and sgi"
        stop
      endif
      mpwipw=0.d0
      sgi=0.d0
      wi0 = 0.d0
      indggq=0
      do ig=1,ngqbarc(iq)
        indggq(indgq(ig,iq))=ig
      enddo

  101 format(/,8x,'Mixed basis for q-point nr.:',i4,/,10x, &
     &           'Number of atomic basis functions:       ',i4,/,10x,   &
     &           'Number of interstitial basis functions: ',i4,/,10x,   &
     &           'Total number of basis functions:        ',i4,/)

      endsubroutine 

      subroutine end_mixbasis
        deallocate(mpwipw,sgi,wi0,indggq) 
      end subroutine 


!BOP
! 
! !IROUTINE: init_mixfun
!
! !INTERFACE:
        subroutine init_mixfun(nrp,maxmix,nat)
        
! !INPUT PARAMETERS:        
          implicit none
          integer, intent(in) :: nat   ! number of ineq. atoms 
          integer, intent(in) :: maxmix ! maximum number of radial mixed basis functions
          integer, intent(in) :: nrp    ! number of radial mesh points

!          
! !DESCRIPTION:
!
! Allocates the radial mixed wave functions
!
!EOP
!
!BOC
          allocate(nmix(nat),                 &
     &             mbl(nat),                  &
     &             umix(nrp,maxmix,nat),      &
     &             usmix(nrp,maxmix,nat),     &
     &             bigl(maxmix,nat),          &
     &             stat=ierr) 
          if(ierr.ne.0) then
            write(6,*) "init_mixfun: fail to allocate momory"
            stop
          endif
        end subroutine init_mixfun

        subroutine reinit_mixfun(nrp,maxmix,nat)
!
! Reallocates the radial mixed wave functions
!
          implicit none
          integer, intent(in) :: nat ! number of ineq. atoms
          integer, intent(in) :: maxmix ! maximum number of radial mixed basis functions
          integer, intent(in) :: nrp ! number of radial mesh points

          call doreallocate(umix,nrp,maxmix,nat)
          call doreallocate(usmix,nrp,maxmix,nat)
          call doreallocate(bigl,maxmix,nat)
        end subroutine reinit_mixfun

        subroutine init_mixmat(nat)
!
! Allocates the arrays that store integrals of radial mixed functions
!
          implicit none
          integer, intent(in) :: nat   ! number of inequivalent atoms

! !LOCAL VARIABLES:
          integer :: packdim

          packdim=nmixmax*(nmixmax+1)/2
          allocate(rtl(nmixmax,nat),                   &
     &             rrint(packdim,nat),                &
     &             stat=ierr)
          if(ierr.ne.0) then
            write(6,*) "init_mixmat: fail to allocate momory"
            stop
          endif

          rtl=0.0d+0
          rrint=0.0d+0
        end subroutine init_mixmat

        subroutine end_mixmat() 
          deallocate(rtl,rrint) 
        end subroutine 

        subroutine init_gaunt(n)
          integer, intent(in) :: n
          allocate(cgcoef(n))
        
        end subroutine init_gaunt

        subroutine end_gaunt

          deallocate(cgcoef)

        end subroutine end_gaunt

        subroutine init_s3r() 
        use struk,  only: nat
        use core,   only: ncoremax
        use bands,  only: nspin 
        use lapwlo, only: nt, nLOmax,lomax
        implicit none

        nrwmax = ncoremax + nt*2 + nLOmax*(lomax+1) 
        
        allocate(s3r(nmixmax,nrwmax,nrwmax,nat,nspin), stat=ierr)
        if(ierr.ne.0) then
          write(*,*) "init_mixmat: fail to allocate momory"
          stop
        endif
        s3r = 0.d0
        end subroutine 

        integer function lo_index_s3r(ilo,l,iat)
!
!       this function returns the index of the LO radial wavefunctions
!       in the s3r matrix, as calculated in mb_calcs3r 
!
        use lapwlo, only:nLO_at,lmax
        use core, only: ncore
        implicit none
        integer:: ilo, l, iat
 
        integer :: ir, lp, jlo

        ir = ncore(iat) + (lmax+1)*2
        do lp = 0, l-1
          do jlo = 1, nlo_at(lp,iat)
            ir = ir + 1
          enddo
        enddo
        lo_index_s3r  = ir + ilo 
        return 
        end function 

        subroutine end_s3r()
        deallocate(s3r)
        end subroutine 
      
      end module mixbasis
!EOC 
