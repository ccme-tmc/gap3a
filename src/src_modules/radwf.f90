!--------------------------------------------------------------
!BOP
!
! !MODULE: radwf
      module radwf
!
! !DESCRIPTION:
!
!This module declares and allocates the vectors where the atomic radwf
!functions are stored. 
!
! !PUBLIC VARIABLES:
        implicit none
        
        logical   :: rel                      !! if .true. perform relativistic calculations
        integer(4):: nrad    !=   881 ! maximum number of radwf   mesh points

        real(8),allocatable :: a(:)          ! r u_l(r,E) at mesh poins.
        real(8),allocatable :: ae(:)         ! At a first step, this is r u_l(r,E-Delta E), which is used 
                                             ! to calculate its final value: r udot_l(r,E)
        real(8),allocatable :: b(:)          ! r us_l(r,E) at mesh poins, where s indicates 2nd. relativistic component.
        real(8),allocatable :: be(:)         ! idem as ae for the second relativistic component.

        real(8), allocatable :: u(:,:,:,:)        ! r u_l(r,E) at mesh poins for all l of atom iat
        real(8), allocatable :: udot(:,:,:,:)     ! r udot_l(r,E) at mesh poins for all l's of atom iat.
        real(8), allocatable :: us(:,:,:,:)       ! r us_l(r,E) at mesh poins for all l's of atom iat.
        real(8), allocatable :: usdot(:,:,:,:)    ! idem as udot for the second relativistic component.

        real(8), allocatable :: ulo(:,:,:,:,:)      ! idem as u for local orbitals. 
        real(8), allocatable :: uslo(:,:,:,:,:)     ! idem as us for local orbitals.

!
! radwf core wave functions 
!
        real(8), allocatable :: vr(:,:,:)   ! spherical part (l=0, m=0)  of the total potential r*V  at mesh point j for atom iat.


!
!EOP
      contains
!BOP
        subroutine init_radwf()
        use lapwlo, only: lmax,lomax,nLOmax
        use bands,  only: nspin
        use struk,  only: nat, nrpt

        implicit none

        integer(4) :: ierr

! !DESCRIPTION:
!
! Allocates the matrices declared in the module
!
!EOP
!BOC
        nrad = maxval(nrpt)

        allocate( a(nrad),ae(nrad),b(nrad),be(nrad),            &
     &            u(     nrad,0:lmax,nat,nspin),                &
     &            udot(  nrad,0:lmax,nat,nspin),                &
     &            us(    nrad,0:lmax,nat,nspin),                &
     &            usdot( nrad,0:lmax,nat,nspin),                &
     &            ulo(   nrad,nLOmax,0:lomax,nat,nspin),        &
     &            uslo(  nrad,nLOmax,0:lomax,nat,nspin),        &
     &            stat=ierr) 
        if(ierr.ne.0) then 
          write(6,*) "init_radwf: error in allocate"
          stop
        endif 
 
        a=0.0d0
        ae=0.0d0
        b=0.0d0
        be=0.0d0

        u=0.0d0
        udot=0.0d0
        us=0.0d0
        usdot=0.0d0

        ulo=0.0d0
        uslo=0.0d0
 
        end subroutine init_radwf

        subroutine end_radwf
        deallocate(a,ae,b,be,u,udot,us,usdot,ulo,uslo)
        end subroutine end_radwf
!EOC
      end module radwf

