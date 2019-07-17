!EOP
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bands
      module bands
      
! !PUBLIC VARIABLES:
    
      integer:: iop_metallic=0           ! option for treating metallicity   
                                         !   0 -- metallic if the input band energies have no gap 
                                         !   1 -- do non-metallic calculations even if the input band energies have no gap
                                         !       nomax and numin are set in terms of electron number, in the meanwhile 
                                         !       unoccupied bands are shifted upward by "band_scissor"
      integer:: iop_nvmk = 0             ! test only -- option to control whether to use k-dependent VBM index 
   
      logical :: metallic                ! Indicate whether the system is metallic 
      real(8) :: band_scissor = 0.0   

      integer :: nspin                ! = 1 or 2 for spin-unpolarized and polarized calculations, respectively 
      integer :: fspin                ! = 2  or 1 for spin-unpolarized and polarized calculations, respectively 

      integer :: nbmax               ! the total number of bands 
      integer :: nbmaxpol            ! Number of bands used for polarization calculations 
      integer :: nbmaxsc             ! the number of bands used for correlation selfenergy
      integer :: nomax               ! Highest partially occupied band
      integer :: numin               ! Lowest partially unoccupied band
      integer :: nomaxs(2)           ! highest partially occupied band for each spin 
      integer :: numins(2)           ! Lowest partially unoccupied band for each spin 
 
      integer :: ibgw                 ! index of the first gw band 
      integer :: nbgw                 ! index of the last gw band 

      integer :: nbandsgw            ! Nubmer of bands to for which GW band structure is calculated  
      integer :: nbands_c            ! number of bands considered in the summation for the corr. self-energies 
      integer :: nbands_x            ! number of bands considered in the summation for the exch. self-energies 

      real(8)  :: eminsc,emaxsc          ! the energy range for the calculation of the correlation selfenergy 
      real(8)  :: eminpol,emaxpol        ! Energy range for polarization matrix calculations 
      real(8)  :: emingw,emaxgw          ! The energy range of the bands to which GW corrections are going to be calculated 
      real(8)  :: spinmom                ! spin momentum per unit cell 
     
      integer, allocatable :: nvmk(:,:)    ! index of highest occupied valence band for each k-point
      integer, allocatable :: nv(:)        ! Number of IPW's for each k-point
      
      real(8), allocatable :: bande(:,:,:)  ! the band energies
      real(8), allocatable :: bande0(:,:,:) ! the band energies, used in task_emac to store KS eigen-energies 
      real(8), allocatable :: eqp(:,:,:)    ! the quasiparticle enregies
      real(8), allocatable :: eqp_im(:,:,:) ! imaginary part of the qp energies

      real(8) :: efermi             ! Fermi energy used internally, which will be always zero
      real(8) :: eferks             ! Fermi energy from KS DFT eigen-energies, absolute value 
      real(8) :: eferqp,eferqp0     ! Fermi energy for quasi-particle energies  
      real(8) :: eferhf             ! Fermi energy for HF band energies (perturbatively) 

      contains
      
      subroutine init_bands(nirkp)
        integer:: nirkp

        integer ierr
      
        allocate(nvmk(nirkp,nspin),                  & 
     &           nv(nirkp),                          &
     &           bande(nbmax,nirkp,nspin),           &
     &           bande0(nbmax,nirkp,nspin),          &
     &           stat=ierr) 
        if(ierr.ne.0) then 
          write(6,*) "init_bands: fail to allocate "
          stop
        endif 
        bande=1.0d+4
      
      end subroutine init_bands
      
      subroutine end_bands
      
      deallocate(nvmk,nv,bande,bande0)
      
      end subroutine end_bands

      subroutine bands_set_nvmk()
!
!     set nvmk, the index of the highest occupied valence band for each k 
!     in terms of Fermi energy and bande  
!
      use kpoints, only: nirkp
      implicit none  
      integer:: isp,irk,ie

      write(6,*) " bands_set_nvmk: set k-dependent VBM index"
      do isp=1,nspin 
        do irk=1,nirkp
          nvmk(irk,isp)=1 
          do ie=2,nbmax 
            if(bande(ie,irk,isp)<=efermi) then 
              nvmk(irk,isp) = ie 
            else
              cycle 
            endif 
          enddo 
          write(6,*) 'isp,irk,nvmk=',isp,irk,nvmk(irk,isp) 
        enddo
      enddo 
      end subroutine 

      subroutine bands_getgap(eband,nb,nk,nsp,efer,egap)
      implicit none 
      integer,intent(in) :: nb,nk,nsp
      real(8),intent(in) :: eband(nb,nk,nsp),efer
      real(8),intent(out):: egap   

      integer:: nvm,ncm,ib,ik,isp 
      real(8):: evm,ecm,enk,eg(nsp) 

      do isp=1,nsp
        evm = -1.e10
        ecm = 1.e10
        nvm = 1
        ncm = nb 
        do ik=1,nk 
          do ib=1,nb 
            enk=eband(ib,ik,isp) 
            if(enk<=efer) then 
              if(enk>=evm) then
                evm = enk 
                nvm = ib 
              endif 
            else  
              if(enk<ecm) then 
                ecm = enk  
                ncm = ib 
              endif 
            endif 
          enddo
        enddo
        if(ncm > nvm) then 
          eg(isp) = ecm-evm 
        else 
          eg(isp) = 0.0 
        endif 
      enddo 
      egap = minval(eg)
      end subroutine 

      end module bands
