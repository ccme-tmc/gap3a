!-----------------------------------------------------------------------
!BOP
!
! !MODULE: selfenergy
      module selfenergy
      implicit none 
      
      
! !PUBLIC VARIABLES:
      integer:: iop_sxc=0                   ! the type of selfenergy 
                                           !   0 -- standard GW selfenergy 
                                           !   1 -- exchange-only
                                           !   2 -- static Coulomb hole screened exchange  
                                           !   3 -- PBE0 hybrid functional 

      integer:: iop_sxc_sc=0               ! self-energy approximation used during self-consistency

      integer:: iop_qsgw=0                 ! control how to define QSGW effective Hamiltonian 
                                           !  0 -- FSK mode B 
                                           !  1 -- FSK mode A 

      integer:: iop_gw0 =0               ! control how the GW0 selfconsistency is done 


      integer:: iop_sc=-1                  ! the level of self-consistency (not including energy-only self-consistency)
                                           !  < 0 - G0W0 with diagonal contributions only 
                                           !  0   - G0W0 with off-diagonal contributions no self-consistency 
                                           !  1   - self-consistency with fixed Hartree potential 
                                           !  2   - full self-consistency 

      integer:: iop_dvh=0                  ! control whether considering the change of Hartree potential during the self-consistency
                                           !  0 -- no, i.e. VH is fixed 
                                           !  1 -- yes

      integer:: isym_kbz = 1               ! 0 - calculate Vxc,Sxc etc in the full BZ
                                           ! 1 - calculate Vxc, Sxc etc in the irreducible BZ

      integer:: n_sigm=1
     
! real quantities  
      real(8) :: beta_sc=1.d0 
      real(8), allocatable ::   selfx(:,:,:)  ! exchange self-energy  
      real(8), allocatable ::   selfc(:,:,:)  ! real part of correlation self-energy 
      real(8), allocatable ::   znorm(:,:,:)  ! the normalization factor in correlation selfenergy 
!
! We use sig*** to complex quantities 
!
      complex(8), allocatable ::  sigx(:,:,:)      ! The selfenergy
      complex(8), allocatable ::  sigx_q(:,:,:)    ! q-dependent exchange selfenergy
      complex(8), allocatable ::  sigc(:,:,:,:)    ! The correllation selfenergy, regular part 
      complex(8), allocatable ::  sigc_q(:,:,:,:)  ! q-dependent correllation selfenergy
      target sigc,sigx,sigc_q,sigx_q
      

! sigm* -- full matrix of self-energies
      complex(8), allocatable :: sigm(:,:,:,:,:)   ! full correlation self-energy matrix including off-diagonal matrix 
      complex(8), allocatable :: sigm_f(:,:,:,:,:) ! the frozen part of \Sigma matrix 
                                                   ! shape: (ibgw:nbgw, ibgw:nbgw,0:nomeg, nirkp, nsp)

      complex(8), allocatable :: qpwf_coef(:,:,:,:) ! coefficients of qp functions in terms of KS vectors 
      complex(8), allocatable :: qpwf_dmat(:,:,:,:) ! density matrix of qp wave functions in terms of KS wave functions 

      complex(8), allocatable :: sxcmn(:,:,:,:)    ! Regularized Sxc matrix 
      complex(8), allocatable :: vxcmn(:,:,:,:)    ! full vxc matrix 
      complex(8), allocatable :: vh0mn(:,:,:,:)    ! matrix elements of the original Hartree potential
      complex(8), allocatable :: vhmn(:,:,:,:)     ! matrix elements of the change of Hartree potential
      complex(8), allocatable :: dvmn_old(:,:,:,:) ! matrix elements of the change of Hartree potential

      complex(8), allocatable :: mwm(:,:,:)    ! M*W*M matrix, auxillilary array for Sigma_c 

      integer :: nomeg_mats = 0      ! the number of Matsubara frequencys 
      real(8) :: beta_mats_eV = 40.0 ! 1/kB*T in eV^{-1}  
      real(8) :: beta_mats = 1088.4  ! 1/kB T (the default: T=290) K used to define the  Mutzbura temperature 
                                    ! in units of Ha
      integer :: iop_intpl_mats = 0  ! the option to choose how to interpolation to the Matsbura frequiences   
      real(8),   allocatable :: omega_mats(:)      ! Matsubura frequencies
      complex(8),allocatable :: hks_wann(:,:,:,:)  ! the KS Hamiltonian in the Wannier basis 
      complex(8),allocatable :: vxc_wann(:,:,:,:)  ! Vxc matrix in the Wannier basis
      complex(8),allocatable :: sxc_wann(:,:,:,:,:) ! Sc matrix in the Wannier basis  
      complex(8),allocatable :: sxc_wann_mats(:,:,:,:,:) ! Sc matrix in the Wannier basis in the Matsubara frequencies 

      real(8), allocatable:: egk(:,:)  ! k-dependent band gaps used for checking convergence in scgw 


! Variables related to analytic continuation 
      integer :: iop_es = 0      !! option for Fermi energy shift (FES)
                                !!  0 -- perturbative G0W0 without energy shift
                                !!  1 -- perturbative G0W0 with energy shift
                                !!  2 -- iterative G0W0 with energy shift
                                !!  3 -- iterative G0W0 without energy shift
                                !! -1 -- selfconsitent GW0 
      integer :: iop_esgw0 = 1  !! option to control whether shift the Fermi energy during 
                                !! self-consistent GW0 
      integer :: iop_ac  = 1     !!  option for analytic continuation (AC) of selfenergy
                                !!  0  -- Pade's approximation (PA)
                                !!  1  -- multipole fitting (MPF) 

      integer(4) :: npol_ac        ! number of poles used in the analytical continuation  
      integer(4) :: npar_ac        ! Number of parameters sused in analytical continuuation
      complex(8), allocatable :: sacpar(:,:,:,:)  ! parameters of the functions fitting the selfenergy

      character(len=80):: info_sxc(0:5)  !! information about the currently used selfenergy approximation 

      character(len=20)::  flag_sxc       !! the flag for QP energy file 
      character(len=120):: fn_selfe      !! the name for self-energy file 

      data info_sxc(0:2) / &
     &  "Selfenergy in the GW approximation", &
     &  "Seflenergy in the exchange only approximation", &
     &  "Selfenergy in the static COHSEX approximation"/ 
      
      integer,private::ierr
      real(8),private::memsiz
      
      contains
      
      subroutine init_selfenergy(iop)
!
!  Initialize arrays used related to the selfenergy calculation
!  iop= 0/1/2 
!    0  --  G0W0 and GW0 that needs only diagonal selfenergy 
!    1  --  quasi-particle self-consistent GW 
!    2  --  initialization for GW+DMFT 
!
      use bands,    only: nspin,ibgw,nbgw,nbandsgw
      use freq,     only: nomeg
      use kpoints,  only: nirkp,nkp
      use task,     only: casename 
      use crpa,     only: nbmax_wf,nbmin_wf,nlmorb
      implicit none 
      integer,intent(in):: iop   

      integer:: i, ie,isp,irk,nom,nktot
      real(8):: w0 
      character(20)::sname="init_selfenergy"
      real(8),parameter:: pi=3.14159265358979

      if(iop_sxc_sc.eq.0) then  ! full GW
        nom=nomeg
      else if (iop_sxc_sc.eq.1) then ! exchange-only
        nom=0
      else if (iop_sxc_sc.eq.2) then ! cohsex
        nom=1
      endif
      n_sigm = nom + 1

      if(isym_kbz.eq.0) then
        nktot = nkp
      else
        nktot = nirkp
      endif

      if(iop_sxc.eq.0) then
        flag_sxc="GW"
      elseif(iop_sxc.eq.1) then
        flag_sxc="HF"
      elseif(iop_sxc.eq.2) then
        flag_sxc="CHSX"
      endif
      fn_selfe=trim(casename)//".selfe_"//trim(flag_sxc)

      if(iop.eq.0) then  !! G0W0 and e-only GW0 
        allocate(                                                       &
     &     sigx(ibgw:nbgw,nirkp,nspin),                               &
     &     sigx_q(ibgw:nbgw,nirkp,nspin),                             &
     &     selfx(ibgw:nbgw,nirkp,nspin),                              &
     &     stat=ierr)
        call errmsg(ierr.ne.0,sname,"Fail to init Sx")
        sigx=0.d0
        sigx_q=0.d0
        selfx = 0.d0 

        if(iop_sxc.eq.0) then 
          allocate(                                               &
     &         sigc(nomeg,ibgw:nbgw,nirkp,nspin),          &
     &         sigc_q(nomeg,ibgw:nbgw,nirkp,nspin),          &
     &         znorm(ibgw:nbgw,nirkp,nspin),                    &
     &         selfc(ibgw:nbgw,nirkp,nspin),                   &
     &         sacpar(npar_ac,ibgw:nbgw,nirkp,nspin ),             &
     &         stat=ierr) 
          call errmsg(ierr.ne.0,sname,"Fail to init GW sxc")
          sacpar=0.d0

        elseif(iop_sxc.eq.2) then
          ! perturbative static COHSEX
          ! Reset frequency mesh: only \omega=0 is needed for static COHSEX calculations
          allocate(                                                       &
     &           sigc(3,ibgw:nbgw,nirkp,nspin),               &
     &           sigc_q(3,ibgw:nbgw,nirkp,nspin),               &
     &           selfc(ibgw:nbgw,nirkp,nspin),                          &
     &           stat=ierr)
          call errmsg(ierr.ne.0,sname,"Fail to init CHSX")
        endif 
        sigc_q = 0.d0 
        sigc = 0.d0
        selfc = 0.d0 

      elseif(iop.eq.1) then !! Quasi-particle selfconsistent GW 
        allocate(sigm(ibgw:nbgw,ibgw:nbgw,0:nom,nirkp,nspin),     &
     &           qpwf_coef(ibgw:nbgw,ibgw:nbgw,nirkp,nspin), &
     &           qpwf_dmat(ibgw:nbgw,ibgw:nbgw,nirkp,nspin), &
     &           vxcmn(ibgw:nbgw,ibgw:nbgw,nirkp,nspin),     &
     &           sxcmn(ibgw:nbgw,ibgw:nbgw,nirkp,nspin),     &
     &           sigm_f(ibgw:nbgw,ibgw:nbgw,0:nom,nirkp,nspin), &
     &           dvmn_old(ibgw:nbgw,ibgw:nbgw,nirkp,nspin),  &
     &           egk(nirkp,nspin),                       &
     &         stat=ierr)
        call errmsg(ierr.ne.0,sname,"Fail to init sigmx etc")

        sigm = 0.d0 
        vxcmn = 0.d0 
        sxcmn = 0.d0 
        sigm_f = 0.d0 
        egk = 0.d0 
        dvmn_old = 0.d0 

        qpwf_coef = 0.d0 
        do isp=1,nspin 
          do irk=1,nirkp
            do ie=ibgw,nbgw 
              qpwf_coef(ie,ie,irk,isp) = 1.d0 
            enddo 
          enddo 
        enddo 

        if(iop_dvh.gt.0) then
          allocate(vh0mn(ibgw:nbgw,ibgw:nbgw,nirkp,nspin),         &
     &             vhmn(ibgw:nbgw,ibgw:nbgw,nirkp,nspin),            &
     &             stat=ierr)
          call errmsg(ierr.ne.0,sname,"Fail to init vhmn")
        endif

      elseif(iop.eq.2)  then  !! gw2wann

        ibgw = min(ibgw, nbmin_wf) 
        nbgw = max(nbgw, nbmax_wf) 
        nbandsgw = nbgw - ibgw + 1 

        allocate(  &
     &   sigm(ibgw:nbgw,ibgw:nbgw,0:nomeg,nktot,nspin),     &
     &   vxcmn(ibgw:nbgw,ibgw:nbgw,nktot,nspin),            &
     &   qpwf_coef(ibgw:nbgw,ibgw:nbgw,nirkp,nspin),        &
     &   qpwf_dmat(ibgw:nbgw,ibgw:nbgw,nirkp,nspin),        &
     &   hks_wann(nlmorb,nlmorb,nktot,nspin),               &
     &   vxc_wann(nlmorb,nlmorb,nktot,nspin),               &
     &   sxc_wann(nlmorb,nlmorb,nomeg,nktot,nspin),         &
     &   stat=ierr)
        call errmsg(ierr.ne.0,sname,"Fail to init sigmx etc")
        sigm = 0.d0
        vxcmn = 0.d0
        hks_wann = 0.d0
        vxc_wann = 0.d0
        sxc_wann = 0.d0

        qpwf_coef = 0.d0 
        do isp=1,nspin 
          do irk=1,nirkp
            do ie=ibgw,nbgw 
              qpwf_coef(ie,ie,irk,isp) = 1.d0 
            enddo 
          enddo 
        enddo 

        if(nomeg_mats.gt.0) then 
          allocate(omega_mats(nomeg_mats),                    & 
     &     sxc_wann_mats(nlmorb,nlmorb,nomeg_mats,nktot,nspin),&
     &     stat=ierr)

          call errmsg(ierr.ne.0,sname,"Fail to init sc_wann_mats")
          
          sxc_wann_mats = 0.d0 

          beta_mats = beta_mats_eV*27.211
          w0 =  pi/beta_mats
          do i=1,nomeg_mats
            omega_mats(i) = w0*(2*i-1)
          enddo
        endif

      else 
        write(6,*) trim(sname)//":iop=",iop," is not supported!"
        stop
      endif 
      end subroutine 
!
!     Clean up selfenergy
!
      subroutine end_selfenergy(iop) 
      integer,intent(in):: iop

      if(iop.eq.0) then 

        selectcase(iop_sxc)
        case(0)  
          deallocate(sigx, sigx_q, sigc,sigc_q,sacpar,znorm,selfx,selfc)
        case(1) 
          deallocate(sigx, sigx_q, selfx)
        case(2) 
          deallocate(sigx, sigx_q, sigc,sigc_q,selfx,selfc)
        endselect

      elseif(iop.eq.1) then 
        deallocate(sigm,qpwf_coef,sigm_f,vxcmn,qpwf_dmat,  &
     &           sxcmn,dvmn_old,egk)
        if(iop_dvh.gt.0) deallocate(vh0mn,vhmn)
      
      elseif(iop.eq.2) then 

        deallocate(sigm,qpwf_coef,qpwf_dmat,vxcmn, &
     &             hks_wann,vxc_wann,sxc_wann)  
        if(nomeg_mats.gt.0) deallocate(sxc_wann_mats,omega_mats)

      endif 

      end subroutine 

      end module selfenergy
!EOP
