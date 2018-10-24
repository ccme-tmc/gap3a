!--------------------------------------------------------------
!BOP
!
! !MODULE: xcpot
      module xcpot

! !DESCRIPTION:
!
! Declares the shared variables used for the matrix elements of the exchange-correlation
! potential
!
! !PUBLIC VARIABLES:

      integer,parameter :: nslmax=  100    ! maximum number of (lm) for the potential in the MT region 
      integer :: lxcmax                     ! maximum number of lm combinations
      integer, allocatable :: lxcm(:)       ! maximum number of lm combinations per atom
      integer :: nksxc                      ! maximum number of  stars
      integer :: nkxc                       ! maximum number of k vectors

      integer:: iop_vxc=0                  ! control how to handle KS xc potential 
                                           !  0 -- calculate vxcmn from wien2k output case.r2v
                                           !  1 -- read from the external file 
! variables related to hybrid functionals 
      logical :: lhybrid =.false. 
      real(8),allocatable:: vxc_hyb(:,:,:)  !! the diagonal elements of the HF-correction part of the hybrid functional
!
! Variables related to LDA+U 
!
      logical :: lvorb=.false.
      logical,allocatable     :: lvorb_at(:)       ! whether a particular atom has vorb corrections 
      integer:: natorb                         ! number of atoms to which vorb is added 

      integer,allocatable  :: iatorb(:)         ! index of the atom to which vorb is added  
      integer,allocatable  :: nlorb(:)          ! number of U-orbital at each atom 
      integer,allocatable  :: lorb(:,:)         ! l value of each U-orbital
      complex(8),allocatable  :: vorb(:,:,:,:,: )  ! vorb(-lmaxorb:lmaxorb,-lmaxorb:lmaxorb,1:2,1:natorb,1:nspin)  
      real(8),allocatable     :: uiorb(:,:,:,:)    ! u-integrals 
      complex(8),allocatable  :: dmorb(:,:,:,:)    ! local density matrix

      integer, pointer     :: lmxc(:,:,:)       ! LM combination in the MTS-xc expansion
      integer, allocatable :: ksxc(:,:)         ! integer coodinates of the representative k-vectors of each 
                                                   ! star in the interstitial xc expansion
      real(8),    pointer     :: vxclm(:,:,:,:)    ! vxc-lm MTS coefficients
      complex(8), allocatable :: vxcs(:,:)         ! vxc-star coefficients
      real(8),    allocatable :: uxcu(:,:,:,:,:)
      complex(8), allocatable :: istpw(:,:)        ! Integral over the interstitial of a star and a pw

      real(8), allocatable    :: vorbnn(:,:,:)     ! vorb-diagonal elements

      contains 
        subroutine init_xcpot 
        use bands,      only: nspin,ibgw,nbgw
        use recipvec,   only: npw2apw
        use kpoints,    only: nirkp
        use struk,      only: nat
        implicit none 
        integer:: ia

      !! KS vxc and vorb diagonal elements used for self-energy correction

        allocate(istpw(npw2apw,nksxc),lvorb_at(nat))
        
        if(lvorb) then
          allocate(vorbnn(ibgw:nbgw,nirkp,nspin))
          vorbnn=0.d0
        endif
      
        lvorb_at=.false.
        if(lvorb) then 
          do ia=1,natorb
            lvorb_at(iatorb(ia))=.true.
          enddo
        endif 

        if(lhybrid) then 
          allocate(vxc_hyb(ibgw:nbgw,nirkp,nspin))
          vxc_hyb = 0.d0
        endif 

        end subroutine 

        subroutine end_xcpot
        deallocate(istpw,lvorb_at)
        if(allocated(vorbnn)) deallocate(vorbnn)
        if(allocated(vxc_hyb)) deallocate(vxc_hyb) 
        end subroutine 

      end module xcpot
!EO
          
