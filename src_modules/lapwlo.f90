!--------------------------------------------------------------
!BOP
! 
! !MODULE: lapwlo
      module lapwlo
        implicit none 
      
! !PUBLIC VARIABLES:          
        logical :: l_newlo = .false.
        integer :: nlmax  = 11 ! The maxinum of l quantum nubmer, corresponding to lmax in wien2k 
        integer :: nLOmax = 1 ! maximum number of local orbitals (LO) per atom per l 
        integer :: nrwmax = 3 ! maximum number of radial wavefunctions per l per atom 
        integer :: lomax  = 3 ! maximum azimutal quantum  number for local orbitals
        integer :: nloat  = 3 ! the same variable as used in wien2k 
        integer :: nlo_tot    ! total number of used local orbitals, including counting over m-quantum nubmer 

        integer :: lmax     ! maximum azimuthal quantum number in (L)APW+lo basis 
        integer :: nt       ! nt == lmax + 1
        integer :: lnsmax   ! maximum considered azimuthal quantum number l of ul for non-spherical hamilton matrix contributionsa

        real(8) :: rkmax   ! rmt*kmax determines matrix size (convergence)
        real(8) :: kmax    ! plane wave cut-off in lapw 

        logical, allocatable :: lapw(:,:)    ! selects wether the orbital  is a LAPW (.true.) or an APW+lo (.false.)
        real(8), allocatable :: elapw(:,:,:) ! linearization energy for (L)APW(+lo)
        real(8), allocatable :: umt(:,:,:,:) ! collective storage of e,p,pe,dp,dpe,pei 
                                             ! umt(1,:,:) -->   e(:,:)  -- linearization energy E_l for atom iat
                                             ! umt(2,:,:) -->   p(:,:)  -- u_l(r,E)  at r = RMT(iat)
                                             ! umt(3,:,:) -->   pe(:,:) -- ud(r,E) at r = RMT(iat)
                                             ! umt(4,:,:) -->   dp(:,:) -- derivative of u_l(r,E) to r at r=RMT(iat)
                                             ! umt(5,:,:) -->   dpe(:,:) --derivative of ud_l(r,E) to r at r=RMT(iat) 
                                             ! umt(6,:,:) -->   pei(:,:) -- norm of ud_l(r,E) integrated over MT sphere 

        logical, allocatable :: loor(:,:)    ! selects which local orbitals  should be used for atom iat
        integer, allocatable :: lvmax_at(:)  ! the maximum of l for occupied orbitals of each atom 
        integer, allocatable :: nlo(:,:)     ! the number of local orbitals (lo and LO) on for l-th orbital on iat-atom
        integer, allocatable :: nLO_at(:,:)  ! the number of local orbitals (LO only) on for l-th orbital on iat-atom
        real(8), allocatable :: elo(:,:,:,:) ! Linearization energy for local orbitals (including lo and LO)

        real(8), allocatable :: umtlo(:,:,:,:,:) ! collective storage of plo,dplo,pi12lo,pe12lo for the LO functions 
                                               ! umtlo(1,...) -->   plo(:,:)    - ulo_l(r,E_2)  at r = RMT
                                               ! umtlo(2,...) -->   dplo(:,:)   - deriv. of u_l(r,E_2) to r at r=rmt(iat)
                                               ! umtlo(3,...) -->   pi12lo(:,:) - integ. of u_l(r,E_1)*u_l(r,E_2)  
                                               ! umtlo(4,...) -->   pe12lo(:,:) - integ. of ud_l(r,E_1)*u_l(r,E_2)
        
        real(8), allocatable :: abcelo(:,:,:,:,:) ! collective stroage of alo,blo,clo and elo 
                                                ! abcelo(1,:,:,:) -> alo(:,:,:) - local orbital coefficient  A_{lm} for atom iat.
                                                ! abcelo(2,:,:,:) -> blo(:,:,:) - local orbital coefficient  B_{lm} for atom iat.
                                                ! abcelo(3,:,:,:) -> clo(:,:,:) - local orbital coefficient  C_{lm} for atom iat.
                                                ! abcelo(4,:,:,:) -> elo(:,:,:) - local orbital expansion energy E_l for atom iat.
!EOP
      contains

        subroutine init_lapwlo(nat,nsp)
          implicit none
          integer :: nat   ! number of inequivalent atoms
          integer :: nsp   ! spin 
          integer :: ierr

          nlo_tot = 0 
          allocate( lapw(0:lmax, nat),              &
     &              nlo(0:lomax, nat),              &
     &              loor(0:lomax, nat),             &
     &              umt(1:6,0:lmax,nat,nsp),        &
     &              nLO_at(0:lmax,nat),             & 
     &              lvmax_at(nat),                  &
     &              stat=ierr)
          if(ierr.ne.0) then  
            write(6,*) "init_lapwlo: Fail to allocate arrays"
            stop
          endif 

          lapw=.true.
          umt = 0.d0
          nlo=0
          nlo_at = 0 
          loor=.true.
          lvmax_at = 0 
        end subroutine init_lapwlo

        subroutine end_lapwlo
          deallocate(lapw,umt,nlo,loor,lvmax_at)
        end subroutine 

      end module lapwlo

