!EOP                  
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bzinteg
      module bzinteg
      use task,    only: fid_outgw, fid_outdbg

! !PUBLIC VARIABLES:
      integer:: iop_q0 = 0              ! control how to handle q=0 singularity
                                        ! 0 : Use auxiliary functions
                                        ! 1 : anisotropically integrate over supercell Brillouin zone
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
      real(8), allocatable :: qmax_gamma(:) ! the length of each grid point on the unit sphere

      real(8) :: vol_gamma ! the volume in reciprocal space of $Z_{\Gamma}$
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

        ! set the coefficient for the head and wing integration
        ! according to iop_q0
        write(fid_outgw, "(A10,I2)") "iop_q0:", iop_q0

        ! Always generate the grid_vec for the current stage
        ! as the time cost is very small for small n_ang_grid
        ! TODO check memory issue for large n_ang_grid
        if(n_ang_grid.le.0) n_ang_grid = 26
        call set_angular_grid(n_ang_grid)
        write(fid_outgw, "(A20,I5)") "Used n_ang_grid: ", n_ang_grid
        !! Set q_max for each vector of the grid
        call set_qmax_gamma

        if(iop_q0.eq.0) then 
          call set_singc_0
        elseif(iop_q0.eq.1)then
          call set_singc_1
        endif

        singc1ex = 0.0
        singc2ex = singc2
        ! Note: singc1co and singc2co will be recalculated if anisotropy
        ! is switched on !!!
        singc1co = singc1
        singc2co = singc2

        write(fid_outgw,'(a,2f12.6)') "singc1ex, singc2ex =", &
            singc1ex,singc2ex
        write(fid_outgw,'(a,2f12.6)') "singc1co, singc2co =", &
            singc1co,singc2co
        call linmsg(fid_outgw,'-',"done init_bzinteg")

        end subroutine init_bzinteg
    
!
! Clean bzinteg 
!
        subroutine end_bzinteg
          deallocate(kcw) 
          deallocate(grid_vec, ang_weight, qmax_gamma)
        end subroutine end_bzinteg


! Calculate coefficient of singular term at q->0 for exchange and
! correlation self-energy
        subroutine set_singc_0
        use struk,   only: vi
        use kpoints, only: nqp
        integer:: iq
        real(8) :: pi2vi,alfa,intf1,intf2,sumf1,sumf2
        real(8) :: singf1, singf2   !! Auxiliary functions for BZ integrals at singular $\Gamma$ point

        pi2vi=pi*pi*vi
        alfa=(1.0d0/(6.0d0*pi2vi))**(1.0d0/3.0d0)
        singc1=0.d0
        singc2=0.d0
        write(fid_outgw,'(A30,I3)') "Generating F Coefs from nqp = ",nqp
        do iq=1,nqp
          call genauxf(iq,alfa,singf1,singf2)
          write(fid_outgw,'(I4,2F12.8)') iq, singf1, singf2
          singc1=singc1 + singf1  
          singc2=singc2 + singf2  
        enddo 

        intf1=1.0d0/(4.0d0*pi2vi*alfa*alfa)
        intf2=1.0d0/(4.0d0*pi2vi*alfa)*sqrt(pi)
        write(fid_outgw,'(A12,F12.8,A12,F12.8)') " SumTildF1 ",singc1,&
     &                                           " SumTildF2 ",singc2
        write(fid_outgw,'(A12,F12.8,A12,F12.8)') "     SumF1 ",intf1,&
     &                                           "     SumF2 ",intf2
        singc1=intf1-singc1
        singc2=intf2-singc2
        singc=singc2
        end subroutine set_singc_0


        subroutine set_singc_1

        use struk, only: vi
        implicit none
        integer :: iang
        real(8) :: BZ_vol

        real(8),external :: ddot
        
        BZ_vol= 8.0d0*pi**3*vi
        singc1=ddot(n_ang_grid,ang_weight,1,qmax_gamma**2,1)/2.0D0
        singc2=ddot(n_ang_grid,ang_weight,1,qmax_gamma,1)

        singc1=singc1/BZ_vol
        singc2=singc2/BZ_vol
        singc=singc2

        end subroutine set_singc_1 


        subroutine set_angular_grid(ngrid)
!
! Sets the angular grid for the integration of the dielectric function around
! the $\Gamma$-point, including anisotropy.
!
        use lebedev_laikov
        use struk,     only: vi
        implicit none
        integer,intent(inout) :: ngrid
        integer :: ierr
        real(8) :: prefactor

        prefactor = 4.0D0*pi

        !! Set the angular grid for integration
        call set_lebedev_laikov_grid(ngrid)
        !! Store the grid
        ngrid = nleb
        allocate(grid_vec(ngrid,1:3), ang_weight(ngrid),  &
     &          qmax_gamma(ngrid),stat=ierr)
        if(ierr.ne.0) then
          write(fid_outgw,*) "init_bzinteg: Fail to allocate memory"
          stop
        endif
        grid_vec(1:ngrid,1)=xleb(1:ngrid)
        grid_vec(1:ngrid,2)=yleb(1:ngrid)
        grid_vec(1:ngrid,3)=zleb(1:ngrid)
        ang_weight(1:ngrid)=wleb(1:ngrid)*prefactor

        call unset_lebedev_laikov_grid

        end subroutine set_angular_grid


        subroutine set_qmax_gamma
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
        use struk,    only: br2,vi
        implicit none

! !LOCAL VARIABLES:
        integer :: i1, i2, i3
        integer :: j, iang
        integer :: isq

        real(8) :: numerator, denominator
        real(8) :: smallq(1:3,1:26)
        real(8) :: q0(1:3)
        real(8) :: maxq
!
! Created Jun. 2008, by RGA
!
!EOP
!BOC
!
        !! determine the k-points closet to gamma
        isq=0
        do i1=-1,1
          do i2=-1,1
            do i3=-1,1
              if(.not.((i1 .eq. 0) .and. (i2 .eq. 0) .and. (i3 .eq. 0)))then
                isq=isq+1
                do j=1,3
                  smallq(j,isq)= dble(i1)*br2(j,1)/dble(nkdivs(1))+         &
     &             dble(i2)*br2(j,2)/dble(nkdivs(2)) + dble(i3)*br2(j,3)   &
     &            /dble(nkdivs(3))
                enddo ! j
              endif
            enddo ! i3
          enddo ! i2
        enddo ! i1

        vol_gamma = 8.0d0*pi**3*vi / dble(product(nkdivs))
        write(fid_outdbg,*) "#qmax", n_ang_grid
        !! determine q_max_gamma
        do iang=1,n_ang_grid
          q0(1:3)=grid_vec(iang,:)
          qmax_gamma(iang)=1.0d+3
          do isq=1,26
            denominator = q0(1)*smallq(1,isq)+q0(2)*smallq(2,isq)+q0(3)*  &
     &                  smallq(3,isq)
            if(denominator .gt. 1.0d-10)then
              numerator = 0.5d0 * sum(smallq(1:3,isq)**2)
              maxq = numerator/denominator
              if(maxq .lt. qmax_gamma(iang)) qmax_gamma(iang) = maxq
            endif
          enddo ! isq
          write(fid_outdbg,"(I4,4F12.5)") iang, grid_vec(iang,:),&
     &        qmax_gamma(iang)
        enddo ! iang
        end subroutine set_qmax_gamma

      end module bzinteg
!EOC        
