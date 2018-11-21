!EOP                  
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bzinteg
      module bzinteg
      use task,    only: fid_outgw

! !PUBLIC VARIABLES:
      integer:: iop_q0 = 0                      ! control how to handle q=0 singularity
      complex(8), allocatable :: kcw(:,:,:,:,:)   ! Weights for convolutions  when k is the integration  variable
      
!
!     The following variables are used for the standard BZ integration, the weights related with 
!     the multiplicity of the k-points is not included  
!
      real(8), allocatable :: kiw(:,:,:)   ! weights for summation over occupied states,  
                                           !   N_c^{-1} \sum_k \sum_n \theta_{nk} F_{nk}
                                           !  
      real(8), allocatable :: kwfer(:,:,:) ! weights for Fermi surface integrations 
                                           !   N_c^{-1} \sum_k \sum_n F_{nk} \delta(enk-efermi)
  
      real(8), allocatable :: kwt_ibz(:)   ! weight of each k-point in the IBZ corresponding to the following type 
                                           ! BZ integral $ N_c^{-1} \sum_k^{\mathrm{BZ}} f(k) $
      real(8), allocatable :: kwt_bz(:)    

      !! variables related to the treatment of the anisotropy
      integer :: n_ang_grid              ! number of point for the angular grid
      real(8), allocatable :: grid_vec(:,:) ! coordinates of the grid
      real(8), allocatable :: ang_weight(:) ! weight of each point
      real(8), allocatable :: k_max_gama(:) ! the length of each grid point on the unit sphere

      real(8) :: singc1ex ! Coefficients for the 1/q singularity part of exchange self-energy
      real(8) :: singc2ex ! Coefficients for the 1/q^2 singularity part of exchange self-energy
      real(8) :: singc1co ! Coefficients for the 1/q singularity part of exchange self-energy
      real(8) :: singc2co ! Coefficients for the 1/q^2 singularity part of exchange self-energy
      real(8) :: singc, singc1, singc2

      real(8),parameter,private:: pi=3.14159265358979d0
!EOP
      contains

!BOP
!
! !IROUTINE: init_bzinteg
!
! !INTERFACE:
        subroutine init_bzinteg 
! !USES:
        use bands,   only: nomax,numin,nbmax,nbmaxpol,nspin
        use core,    only: ncg
        use freq,    only: nomeg
        use kpoints, only: nkp,nirkp
        implicit none

! !LOCAL VARIABLES:        
        integer :: ierr,iq

!EOP
!BOC
        call linmsg(fid_outgw,'-',"init_bzinteg")
        allocate(kcw(ncg+nomax,numin:nbmaxpol,nkp,nomeg,nspin),       &
     &           stat=ierr) 
        if(ierr.ne.0) then 
          write(6,*) "init_bzinteg: Fail to allocate memory"
          stop
        endif 

        kcw=0.d0

        !write(fid_outgw, "(A10,I2)") "iop_q0:", iop_q0
        if(iop_q0.eq.0) then 
          call set_singc_ex
          call set_singc_co
!        else 
!          if(n_ang_grid.le.0) n_ang_grid = 26
!          allocate(grid_vec(1:3,n_ang_grid), ang_weight(n_ang_grid),    &
!     &            k_max_gama(n_ang_grid),stat=ierr)
!          if(ierr.ne.0) then
!            write(fid_outgw,*) "init_bzinteg: Fail to allocate memory"
!            stop
!          endif
!          grid_vec=0.0d0
!          ang_weight=0.0d0
!          call set_angular_grid
!          call set_singc_1
        endif

        write(fid_outgw,'(a,2f12.6)') "singc1ex, singc2ex =", &
            singc1ex,singc2ex
        write(fid_outgw,'(a,2f12.6)') "singc1co, singc2co =", &
            singc1co,singc2co

        end subroutine init_bzinteg
    
!
! Clean bzinteg 
!
        subroutine end_bzinteg
          deallocate(kcw) 
        end subroutine end_bzinteg

! Calculate coefficient of singular term at q->0 for exchange and
! correlation self-energy
        subroutine set_singc_ex
        use struk,   only: vi
        use kpoints, only: nqp
        integer:: iq
        real(8) :: pi2vi,alfa,intf1,intf2,sumf1,sumf2
        real(8) :: singf1, singf2   !! Auxiliary functions for BZ integrals at singular $\Gamma$ point

        pi2vi=pi*pi*vi
        alfa=(1.0d0/(6.0d0*pi2vi))**(1.0d0/3.0d0)
        singc1=0.d0
        singc2=0.d0
        write(fid_outgw,'(A30,I3)') "Generating Coefs from nqp = ", nqp
        do iq=1,nqp
          call genauxf(iq,alfa,singf1,singf2)
          singc1=singc1 + singf1  
          singc2=singc2 + singf2  
        enddo 

        intf1=1.0d0/(4.0d0*pi2vi*alfa*alfa)
        intf2=1.0d0/(4.0d0*pi2vi*alfa)*sqrt(pi)
        singc1=intf1-singc1
        singc2=intf2-singc2

        singc1ex = singc1
        singc2ex = singc2
        singc=singc2
        endsubroutine 


        subroutine set_singc_co
        use anisotropy, only: iop_aniso
        integer:: iang

        if(iop_aniso.eq.-1)then
            singc1co=singc1ex
            singc2co=singc2ex
        else
            singc1co=singc1ex
            singc2co=singc2ex
            !singc2co=1.0D0
        endif
        !singc=0.0d0
        !do iang=1,n_ang_grid
        !  singc=singc+ang_weight(iang)*k_max_gama(iang)
        !enddo
        end subroutine 


        subroutine set_angular_grid
!
! Sets the angular grid for the integration of the dielectric function around
! the $\Gamma$-point, including anisotropy.
!
        use lebedev_laikov
        use struk,     only: vi
        implicit none
        integer :: iang
        real(8)    :: BZ_vol, prefactor

        BZ_vol= 8.0d0*pi**3*vi
        prefactor = 4.0d0*pi/BZ_vol

        !! Set the angular grid for integration
        call set_lebedev_laikov_grid(n_ang_grid)

        !! Store the grid
        grid_vec(1,1:n_ang_grid)=xleb(1:n_ang_grid)
        grid_vec(2,1:n_ang_grid)=yleb(1:n_ang_grid)
        grid_vec(3,1:n_ang_grid)=zleb(1:n_ang_grid)
        ang_weight(1:n_ang_grid)=wleb(1:n_ang_grid)*prefactor

        call unset_lebedev_laikov_grid

        !! Set K_max for each vector of the grid
        call set_kmax_gama

        end subroutine set_angular_grid


        subroutine set_kmax_gama
! !DESCRIPTION:
!
!! Determines the intersection of the grid vectors (for angular
!integration) with the border of the Zone around $\Gamma$. The Zone border
!is given by the planes equidistant from $\Gamma$ and $\vec{k}_i$ of the BZ
!integration grid, which is closest to $\Gamma$. The equation of tha plane
!is given by $\vec{x}\cdot \vec{k}_i = 0.5 \frac{|\vec{k}_i|^2}{2}$. For
!each grid vector $\vec{g}_j$ we have then: $\verbatim{maxk}_{ji} = 0.5
!\frac{|\vec{k}_i|^2}{2 \hat{g}_j\cdot \vec{k}_i}$.  And
!\verbatim{k_max_gamma}$_i$ is given by the minimum of the positive
!$\verbatim{maxk}_{ji}$.
!
! !USES:

        use kpoints,  only: nkdivs
        use struk,    only: br2
        implicit none

! !LOCAL VARIABLES:
        integer :: i1, i2, i3
        integer :: j, iang
        integer :: isk

        real(8) :: numerator, denominator
        real(8) :: smallk(1:3,1:26)
        real(8) :: k0(1:3)
        real(8) :: maxk
!
! Created Jun. 2008, by RGA
!
!EOP
!BOC
!
        !! determine the k-points closet to gamma
        isk=0
        do i1=-1,1
          do i2=-1,1
            do i3=-1,1
              if(.not.((i1 .eq. 0) .and. (i2 .eq. 0) .and. (i3 .eq. 0)))then
                isk=isk+1
                do j=1,3
                  smallk(j,isk)= dble(i1)*br2(j,1)/dble(nkdivs(1))+         &
     &             dble(i2)*br2(j,2)/dble(nkdivs(2)) + dble(i3)*br2(j,3)   &
     &            /dble(nkdivs(3))
                enddo ! j
              endif
            enddo ! i3
          enddo ! i2
        enddo ! i1

        !! determine k_max_gamma
        do iang=1,n_ang_grid
          k0(1:3)=grid_vec(1:3,iang)
          k_max_gama(iang)=1.0d+3
          do isk=1,26
            denominator = k0(1)*smallk(1,isk)+k0(2)*smallk(2,isk)+k0(3)*  &
     &                  smallk(3,isk)
            if(denominator .gt. 1.0d-10)then
              numerator = 0.5d0 * sum(smallk(1:3,isk)**2)
              maxk = numerator/denominator
              if(maxk .lt. k_max_gama(iang)) k_max_gama(iang) = maxk
            endif
          enddo ! isk
        enddo ! iang
        end subroutine set_kmax_gama

      end module bzinteg
!EOC        
