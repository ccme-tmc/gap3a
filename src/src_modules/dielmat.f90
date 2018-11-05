!-----------------------------------------------------------------------
!BOP
!
! !MODULE: dielmat
!  a module for the polarization and dielectric function and screened
!  Coulomb interaction 
!
      module dielmat
      use constants,  only: sqrt3
      use task, only: fid_outgw
      implicit none 
      
! !PUBLIC VARIABLES:
      integer :: iop_epsw = 0
      integer :: iop_mask_eps=0          ! option to skip some transitions when calculating eps
                                         ! needed for constrained RPA calculations 
                                         !  0 -- no skipping 
                                         !  1 -- skip in terms of band indices, i.e. bands in a certain range will be skipped
                                         !  2 -- skip in terms of energies, i.e. bands in a certain energy window will be skipped 
      integer :: iop_drude=1             ! control whether to include the Drude term  
 
      real(8):: q0_eps(3)                ! this gives directional components when taking q-> 0 limit , by default it is taken as 
                                         ! (1,1,1)/sqrt(3) 
      real(8):: c0_head                  ! the coefficient of the head term at \omega=0, which is equal to the square of 
                                         ! the plasmon frequency in metallic systems 
      real(8):: omega_plasma=-1.0
      real(8):: eta_head = 0.01

      complex(8), allocatable :: head(:)       ! the head of the dielectric matrix
      complex(8), allocatable :: eps(:,:,:)    ! the dielectric matrix
      target eps 
      complex(8), allocatable :: epsw1(:,:)    ! the vertical   wing of the dielectric matrix  
      complex(8), allocatable :: epsw2(:,:)    ! the horizontal wing of the dielectric matrix
      complex(8), allocatable :: emac(:,:)     ! the macroscopic dielectric function  
                                               !   emac(1,:) -- without local field effect
                                               !   emac(2,:) -- with    local field effect 
      character(len=2) :: bandtype             ! 'KS' or 'GW' energies used for the dielectric function 
!
! Variables related to constrained RPA. They may be defined in a seperate module in the future  
!
      real(8):: wt_excl=0.0                    !! the excluding weight for the transition 
      integer:: noc_excl,nun_excl
      real(8):: occ_win(2),unocc_win(2) 
      integer,allocatable:: ioc_excl(:),iun_excl(:) !! indices of occupied/unoccupied bands to be excluded from eps calculations
      real(8),allocatable:: mask_eps(:,:)     ! a mask used to exclude some transitions when calculating eps
                                               ! this is needed in the constrained RPA, it can be also used when 
                                               ! calculating eps_mac with certain transition neglected, useful 
                                               ! for some analysis 
      integer,private:: ierr
      logical,private:: ldbg=.false.

      contains
      
      subroutine init_dielmat(iq,iomfirst,iomlast)
      use bands,    only: nomaxs,numins,nbmaxpol,nspin
      use core,     only: ncg_p 
      use constants,only: czero,cone
      use kpoints,  only: nirkp
      use mixbasis, only: matsiz
      integer:: iq,iomfirst,iomlast

      integer:: nomx,numn,ie1,ie2
      logical:: lexcl_ie1, lexcl_ie2

      character(120):: sname="init_dielmat"

      nomx = maxval(nomaxs)
      numn = minval(numins)

      c0_head = 0.0
      if(iq.eq.0) then  !! initialize only head-related arrays, in this
                        !! case, the mixed basis is not needed 
        allocate(head(iomfirst:iomlast),&
     &           mask_eps(numn:nbmaxpol,nomx+ncg_p),                       &
     &           stat=ierr)
        if(ierr.ne.0) then
          write(fid_outgw,*) "init_dielmat: Fail to allocate head"
          stop
        endif
        head = czero
        mask_eps = 1.0
        return 
      endif 

      allocate(eps(matsiz,matsiz,iomfirst:iomlast),                     &
     &         mask_eps(numn:nbmaxpol,1:nomx+ncg_p),                    &
     &         stat=ierr)

      if(ierr.ne.0) then
        write(fid_outgw,*) "init_dielmat: Fail to allocate eps"
        stop
      endif
      eps = czero
      mask_eps = 1.0

      if(iq.eq.1) then 
      ! TODO why iq=1 need head and epsw1/2?
        if(ldbg) write(fid_outgw,*) trim(sname)//": allocate head etc"
        allocate(head(iomfirst:iomlast),                                &
     &           epsw1(matsiz,iomfirst:iomlast),                        &
     &           emac(2,iomfirst:iomlast),                              &
     &           epsw2(matsiz,iomfirst:iomlast), &
     &           stat=ierr)
        if(ierr.ne.0) then
          write(6,*) "init_dielmat: Fail to allocate head,epsw1/2,emac"
          stop
        endif
        head=czero
        epsw1=czero
        epsw2=czero
        emac=czero
      endif 
  100 format("init_dielmat: required memory =",f8.2,"MByte")
      end subroutine 

      subroutine end_dielmat(iq)
        integer:: iq

        if(iq.eq.0) then 
          deallocate(head,mask_eps)
        else 
          deallocate(eps,mask_eps)
          if(iq.eq.1) then 
            deallocate(head,epsw1,epsw2,emac)
          endif 
        endif 
      end subroutine 

      end module 
!EOP
