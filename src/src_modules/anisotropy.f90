MODULE ANISOTROPY

    ! this module defines variables and subroutines utilized for dealing
    ! with the anisotropy of dielectric matrix for $q\to0$

    use mixbasis, only: matsiz
    use constants, only: czero,cone,pi
    use task, only: fid_outgw, fid_outdbg

    integer :: iop_aniso = -1  ! control whether to consider the anisotropy of dielectric function around Gamma
                               ! -1 -- use q0_eps, no anisotropy
    complex(8), allocatable :: vec_u_ani(:,:,:)   ! vector u
    complex(8), allocatable :: vec_t_ani(:,:,:)   ! vector t
    complex(8), allocatable :: vec_a_ani(:,:,:)   ! vector a
    complex(8), allocatable :: vec_b_ani(:,:,:)   ! vector b
    complex(8), allocatable :: ten_p_ani(:,:,:)   ! tensor P
    complex(8), allocatable :: ten_a_ani(:,:,:)   ! tensor A
    ! Following quantities should be determined by a Lebedev-Laikov grid
    integer :: nq0 = 0                   ! number of q0 for angular integration
    real(8), allocatable :: q0_sph(:,:)  ! direction of q0, similar to q0_eps
    real(8), allocatable :: wt_q0_sph(:) ! weight of q0_sph
    integer :: lmax_q0 = 4               ! maximum angular momentum to expand eps on q0_sph

    complex(8), allocatable :: head_q0(:,:)       ! the head (nq0, nomega) for q0 
    complex(8), allocatable :: wing1_q0(:,:,:)    ! the vertical wing (nq0, matsiz, nomega) for q0 
    complex(8), allocatable :: wing2_q0(:,:,:)    ! the horizontal wing (nq0, matsiz, nomega) for q0 

    ! smallq define the proximity around Gamma point
    real(8) :: smallq(26,3)
    integer :: smallq_div = 1
    character(len=3) :: smallq_type = "bzl"       ! BZ like ('bzl')or spherical ('sph')
    real(8) :: vol_q0                             ! volume of the defined proximity near Gamma
    real(8),allocatable :: qmax_q0(:)
    real(8),allocatable :: w_q0(:)
    real(8) :: norm_w_q0


    CONTAINS

        SUBROUTINE init_aniso(iomfirst, iomlast)

            use lebedev_laikov

            integer,intent(in) :: iomfirst, iomlast
            integer :: ierr

            call set_lebedev_laikov_grid(nq0)
            nq0 = nleb

            if(iop_aniso.ne.-1)then
                write(fid_outgw,*) "Anisotropy switched on"
                ! TODO avoid body_q0 in future for integration 
                allocate(vec_u_ani(3,matsiz,iomfirst:iomlast),     &
                &        vec_t_ani(3,matsiz,iomfirst:iomlast),     &
                &        vec_a_ani(3,matsiz,iomfirst:iomlast),     &
                &        vec_b_ani(3,matsiz,iomfirst:iomlast),     &
                &        ten_p_ani(3,3,iomfirst:iomlast),          &
                &        ten_a_ani(3,3,iomfirst:iomlast),          &
                &        q0_sph(1:nq0,3),                          &
                &        wt_q0_sph(1:nq0),                         &
                &        head_q0(1:nq0,iomfirst:iomlast),          &
                &        wing1_q0(1:nq0,matsiz,iomfirst:iomlast),  &
                &        wing2_q0(1:nq0,matsiz,iomfirst:iomlast),  &
                &        qmax_q0(1:nq0),                           &
                &        w_q0(1:nq0),                              &
                &        stat=ierr)
                if(ierr.ne.0) then
                    write(fid_outgw,*) " - init_aniso: Fail to allocate aniso"
                    stop
                else
                    write(fid_outgw,*) " - init_aniso: success"
                endif
            endif

            ten_p_ani(:,:,:) = czero
            ten_a_ani(:,:,:) = czero
            vec_u_ani(:,:,:) = czero
            vec_t_ani(:,:,:) = czero
            vec_a_ani(:,:,:) = czero
            vec_b_ani(:,:,:) = czero
            head_q0(:,:) = cone
            wing1_q0(:,:,:) = czero
            wing2_q0(:,:,:) = czero
            q0_sph(:,1) = xleb(:)
            q0_sph(:,2) = yleb(:)
            q0_sph(:,3) = zleb(:)
            wt_q0_sph(:) = wleb(:)
            call unset_lebedev_laikov_grid
            call init_smallq

        END SUBROUTINE init_aniso


        SUBROUTINE end_aniso

            deallocate(vec_u_ani, vec_a_ani, vec_b_ani, vec_t_ani, ten_p_ani, ten_a_ani,&
    &          q0_sph, wt_q0_sph, head_q0, wing1_q0, wing2_q0, qmax_q0, w_q0)
            write(fid_outgw,*) " - end_aniso: success"

        END SUBROUTINE end_aniso


        SUBROUTINE init_smallq
        ! adopted from set_kmax_gama in bzinteg
        use kpoints,  only: nkdivs
        use struk,    only: br2, vi
        implicit none
! !LOCAL VARIABLES:
        integer :: i1, i2, i3
        integer :: j
        integer :: isq,iq0
        real(8) :: qmax_tmp, denominator, numerator
        logical :: ldbg = .true.

        real(8),external :: ddot
!EOP
!BOC
        ! check if smallq_div is legal, i.e. larger than 0
        if (smallq_div.le.0) then
          write(fid_outgw,*) " - init_smallq: illegal smallq_div"
          stop
        endif
        !! determine the q-points closet to gamma
        isq=0
        do i1=-1,1
          do i2=-1,1
            do i3=-1,1
              if(.not.((i1 .eq. 0) .and. (i2 .eq. 0) .and. (i3 .eq. 0)))then
                isq=isq+1
                do j=1,3
                  smallq(isq,j)=dble(i1)*br2(j,1)/dble(nkdivs(1)*smallq_div)+&
     &                          dble(i2)*br2(j,2)/dble(nkdivs(2)*smallq_div)+&
     &                          dble(i3)*br2(j,3)/dble(nkdivs(3)*smallq_div)
                enddo ! j
              endif
            enddo ! i3
          enddo ! i2
        enddo ! i1

        !! determine qmax with initialized by a large value 1000
        qmax_q0 = 1.0D3

        ! calculate qmax along q0_sph in the region defined by smallq
        ! when smallq_div = nkdivs = 1, qmax with q0_sph fills the fisrt BZ
        do iq0=1,nq0
          do isq=1,26
            denominator = sum(q0_sph(iq0,:)*smallq(isq,:))
            if(denominator .gt. 1.0d-10)then
              numerator = 0.5d0 * sum(smallq(isq,1:3)**2)
              qmax_tmp = numerator/denominator
              if(qmax_tmp .lt. qmax_q0(iq0)) qmax_q0(iq0) = qmax_tmp
            endif
          enddo ! isk
        enddo ! iq0
        vol_q0 = vi*(2.0D0*pi)**3/dble(product(nkdivs)*smallq_div**3)

        w_q0(:) = qmax_q0(:)**3/3.0D0/vol_q0
        ! normalization of w_q0(q) in the defined proximity
        norm_w_q0 = 4*pi*ddot(nq0,w_q0,1,wt_q0_sph,1)

        if(ldbg) then
            write(fid_outdbg,"(A7)") "smallq:"
            do isq=1,26
                write(fid_outdbg,"(I3,3F12.6)") isq, smallq(isq,:)
            enddo
            write(fid_outdbg,"(A30, f12.6)") "Volume of smallq (a.u.^-3):", vol_q0
            write(fid_outdbg,"(A15,I4)") "qmax with nq0 ", nq0
            do iq0=1,nq0
                write(fid_outdbg,"(I4,4F12.6)") iq0, q0_sph(iq0,:),qmax_q0(iq0)
            enddo
            write(fid_outdbg,"(A25,F12.7)") "Normalization of w(q) = ", norm_w_q0
        endif

        END SUBROUTINE init_smallq


        SUBROUTINE angint_eps_sph(iomega, head, wv, wh, bodyinv)
        implicit none
        integer,intent(in) :: iomega 
        ! averaged head, vertical and horizontal wing of dielectri matrix inverse
        complex(8),intent(out) :: head 
        complex(8),intent(out) :: wv(matsiz)
        complex(8),intent(out) :: wh(matsiz)
        ! averaged body. The inverse of the body of dielectric matrix
        ! is parsed in
        complex(8),intent(inout) :: bodyinv(matsiz,matsiz) 

        ! local variables
        integer :: lmsq
        complex :: head_tmp
        complex(8),allocatable :: sph_harms(:,:)
        complex(8),allocatable :: h_lm(:)
        complex(8),allocatable :: s_lm(:)
        complex(8),allocatable :: h_w(:)  ! head of invers times w_q0
        integer :: iq0  ! Counter: runs over nq0
        integer :: iom  ! Counter: runs over iom_f:iom_l
        integer :: im   ! Counter: runs over matsiz
        logical :: ldbg=.true.

        external ylm
        complex(8),external :: zdotu

        lmsq = (lmax_q0+1)**2
        allocate(sph_harms(lmsq,nq0), &
     &           h_w(nq0)            &
     &          )

        ! calculate spherical harmonics at q0_sph
        do iq0=1,nq0
            call ylm(q0_sph(iq0,:),lmax_q0,sph_harms(:,iq0))
        enddo

        head_tmp = head
        h_w = head_q0(:,iomega) * w_q0
        head=4.0D0*pi*zdotu(nq0,h_w,1,cmplx(wt_q0_sph(:),0.0D0,8),1)/norm_w_q0
        if(ldbg)then
            write(*,"(A4,2f13.4,A4,2f13.4)") "OH", head_tmp,"AH",head
            write(fid_outgw,"(A6,I3,A2,2e13.4)") "e-100",iomega,"A",head
        endif

        deallocate(sph_harms,h_w)

        END SUBROUTINE angint_eps_sph

END MODULE ANISOTROPY

