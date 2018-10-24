!--------------------------------------------------------------
!BOP
! !MODULE: core
      module core

! !PUBLIC VARIABLES:
        integer :: iop_core=0 ! option for core treatment 
                              !  iop_core  = 0  --- core states are included in all calculations 
                              !            = 1  --- core states are included in x/c selfenergies but not in polarization energies  
                              !            = 2  --- core states are excluded in all calculations, but kept in 
                              !                     the construction of mixed basis 
                              !            = 3  -- core states are excluded completely
    
        integer :: iopdc=0   ! Option for deep core, as those defined by ecoremin  
                                !  iopdc   = 0  --- dc states are treated as normal core states
                                !          = 1  --- dc states are not included in generating mixed basis  
                                !          = 1  --- dc states are not included in all calculations 
   
        logical :: core_ortho=.false.          ! control whether to force orthogonalization of core states with valence states

        integer :: ncoremax                 ! Max. num of core states  
        integer :: lcoremax                 ! Maximum l for the core states
        integer :: nclmmax                  ! Max. num of core states including lm
        integer :: ncg                      ! Total num of core states including lm
        integer :: ncg_p                    ! the total number of core states considered for polarization 
        integer :: ncg_c                    ! Total num of core states including lm considered for correlation self-energies
        integer :: ncg_x                    ! Total num of core states including lm considered for exchange self-energies 

        integer, allocatable :: corind(:,:) ! indexes of the core states
                                               !  corind(1,:)=ineq. atom
                                               !  corind(2,:)=eq. atom
                                               !  corind(3,:)=core state
                                               !  corind(4,:)=l
                                               !  corind(5,:)=m

        integer,allocatable :: clmind(:,:)  ! index of core states (including l,m)  at each atom  


        integer :: nelfc                  ! the number of frozen core electrons 
        integer :: nelcor                 ! the number of core electrons included in *.corewf  
        real(8)  :: emin                       ! the energy that separate core states and valence states 

        real(8), allocatable :: eigcore(:,:,:)    ! Core states self-energies
        integer, allocatable :: ncore(:)          ! number of core states (nl) of atom iat
        integer, allocatable :: nclm(:)        ! number of core states including (nlm) of atom iat

        character(len=3),allocatable :: symbl(:,:)  ! symbol for core state (1s,2s,1p,1p*...)
        integer,      allocatable :: npqn(:,:)  ! principal quantum number of core state icore
        integer,      allocatable :: kappa(:,:) ! orbital momentum  relativistic quantum number kappa =-s(j+1/2)
        integer,      allocatable :: lcore(:,:) ! angular quantum number 
        integer,      allocatable :: ocmax(:,:) ! maximum occupation of core state icore
!
! radial core wave functions
!
        real(8),allocatable :: ucore(:,:,:,:)
        real(8),allocatable :: uscore(:,:,:,:)


        
! !DESCRIPTION:
!
! Declares and allocates the variables for storing the core wave
! functions and the core state quantum numbers
!
!EOP
!BOC
      contains
        subroutine init_core(nat,nsp)
          use radwf, only: nrad
          implicit none
          integer,intent(in) :: nat  ! number of inequivalent atoms
          integer,intent(in) :: nsp  ! number of spin channel

          integer:: ierr
          external outerr
          
          allocate(ocmax(ncoremax,nat),       &
     &             kappa(ncoremax,nat),       &
     &             lcore(ncoremax,nat),       &
     &             npqn(ncoremax,nat),        &
     &             symbl(ncoremax,nat),       &
     &             eigcore(ncoremax,nat,nsp),    &
     &             ucore(nrad,ncoremax,nat,nsp),        &
     &             uscore(nrad,ncoremax,nat,nsp),       &
     &             stat=ierr)
          if(ierr.ne.0) call outerr("init_core",                        &
     &                              "Error in allocate arrays")
        end subroutine init_core

        subroutine end_core
          deallocate(npqn,kappa,ocmax,symbl,ncore,eigcore,lcore,ucore,uscore)
        end subroutine end_core

!EOC
      end module core

