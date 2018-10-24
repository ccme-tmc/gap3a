!-----------------------------------------------------------------------
!BOP
!
! !MODULE: recipvec
      module recipvec

! !PUBLIC VARIABLES:

      integer(4) :: npw                       ! Total number of pws for integrals 
      integer(4) :: npw2apw                   ! Number of pws for lapw overlaps
      integer(4) :: maxng                     ! Maximum number of G-vectors used in all possible situations 
      integer(4) :: maxngk                    ! Maximum number of LAPW basis functions
      integer(4) :: maxngq                    ! Maximum number of IPW for the mixed basis
      integer(4) :: maxngqbarc                ! Maximum number of IPW for the coulomb matrix
      integer(4), allocatable :: ngk(:)       ! number of LAPW basis functions 
      integer(4), allocatable :: ngkir(:)     ! number of LAPW basis functions   
      integer(4), allocatable :: ngq(:)       ! number of interstitial plane waves for the mixed basis
      integer(4), allocatable :: ngqbarc(:)   ! number of interstitial plane waves for the coulomb matrix
      integer(4), allocatable :: ngqlen(:)    ! number of shells of G vectors
      integer(4) :: maxngqlen                 ! maximum number of shells
      integer(4), allocatable :: indgq(:,:)   ! Index of the G-vector in G+q
      integer(4), allocatable :: indgk(:,:)   ! Index of the G-vector in G+k     
      integer(4), allocatable :: indgkir(:,:) !  Index of the G-vector in G+k      
      integer(4), allocatable :: indgqlen(:,:)! Index of the G-vector in G+k'     
      
      integer(4), allocatable :: gindex(:,:)  ! Integer coordinates of the G vector 
      
      integer(4) :: ngmin, ngmax
      real(8), allocatable :: gqleng(:,:)


      integer(4), allocatable :: ig0(:,:,:)

      real(8), pointer ::  kpg(:,:) ! k+G vector in cartesian coords.

      real(8), pointer ::  kppg(:,:) ! k'+G vector in cartesian coords.

      real(8), pointer ::  rk(:)    ! the length of k+G 


! !DESCRIPTION:
!
!This module declares and the matrix sgi and its size nppw
!
!EOP
      contains
!BOP
!
! !IROUTINE: init_recipvec
!
! !INTERFACE:      
        subroutine init_recipvec

! !USES:
        use kpoints,  only: nkp,nirkp
        implicit none        
!EOP
!BOC
          allocate( rk(1:maxngk) )
          allocate( kpg(1:3,1:maxngk) )
          allocate( kppg(1:3,1:maxngk) )
          allocate(indgk(maxngk,nkp))
          allocate(indgkir(maxngk,nirkp))

          rk=0.0d0
          kpg=0.0d0
          kppg=0.0d0

        end subroutine init_recipvec
!EOC
!BOP
!
! !IROUTINE: end_recipvec
!
! !INTERFACE:        
        subroutine end_recipvec
!EOP
!BOC
        deallocate(rk, kpg, kppg,indgk,indgkir )

        end subroutine end_recipvec

              
      end module recipvec
!----------------------------------------------------------------

