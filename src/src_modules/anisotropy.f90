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

    complex(8), allocatable :: head_q0(:)     ! the head (nq0) for q0 at some freq
    complex(8), allocatable :: q0_va(:,:)    ! q0 dot vector a
    complex(8), allocatable :: q0_vb(:,:)    ! q0 dot vector b

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
                &        head_q0(1:nq0),                           &
                &        q0_va(1:nq0,matsiz),                      &
                &        q0_vb(1:nq0,matsiz),                      &
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
            head_q0(:) = cone
            q0_va(:,:) = czero
            q0_vb(:,:) = czero
            q0_sph(:,1) = xleb(:)
            q0_sph(:,2) = yleb(:)
            q0_sph(:,3) = zleb(:)
            wt_q0_sph(:) = wleb(:)
            call unset_lebedev_laikov_grid
            call init_smallq

        END SUBROUTINE init_aniso


        SUBROUTINE end_aniso

            deallocate(vec_u_ani, vec_a_ani, vec_b_ani, vec_t_ani, &
    &                  ten_p_ani, ten_a_ani,&
    &                  q0_sph, wt_q0_sph, head_q0, &
    &                  q0_va, q0_vb, &
    &                  qmax_q0, w_q0)
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
        integer :: lmsq ! (lmax_q0+1)^2
        complex :: head_tmp, bodyinv_tmp
        complex(8),allocatable :: sph_harms(:,:)
        complex(8),allocatable :: h_lm(:)
        complex(8),allocatable :: a_lm(:,:),b_lm(:,:)
        complex(8),allocatable :: h_w(:)  ! head of invers times w_q0
        complex(8),allocatable :: q_aob_q(:) 
        integer :: iq0               ! Counter: runs over nq0
        integer :: iom               ! Counter: runs over iom_f:iom_l
        integer :: im,jm             ! Counter: runs over matsiz
        integer :: ilm,lm1,lm2,lm3   ! Counter: runs over lmsq
        integer :: l1,l2,l3,m1,m2,m3 ! Counters: runs over angular moment
        integer :: i,j               ! Counters: Cartesian axis
        logical :: ldbg=.true.
        complex(8) :: ccoefcoul_q0, ten_aob_tmp(3,3)

        external ylm
        complex(8),external :: zdotu,ten_rvctrv

        lmsq = (lmax_q0+1)**2
        allocate(sph_harms(lmsq,nq0), &
     &           h_w(nq0),            &
     &           q_aob_q(nq0),        &
     &           h_lm(lmsq),          &
     &           a_lm(lmsq,matsiz),   &
     &           b_lm(lmsq,matsiz)    &
     &          )
        ! calculate head_q0 and q0\cdot wing1, q0\cdot wing2
        do iq0=1,nq0
          ccoefcoul_q0=cmplx(4.0D0*pi,0.0D0,8)
          ! head
          head_q0(iq0) = cone / &
     &     (cone+ccoefcoul_q0*ten_rvctrv(3,ten_a_ani(:,:,iomega),q0_sph(iq0,:)))
          ! wings
          do im=1,matsiz
          ! TODO optimize with ZGEMM
            q0_va(iq0,im)=sum(vec_a_ani(:,im,iomega)*cmplx(q0_sph(iq0,:),0.0D0,8))
            q0_vb(iq0,im)=sum(vec_b_ani(:,im,iomega)*cmplx(q0_sph(iq0,:),0.0D0,8))
          enddo
        enddo

        ! calculate spherical harmonics at q0_sph
        do iq0=1,nq0
            call ylm(q0_sph(iq0,:),lmax_q0,sph_harms(:,iq0))
        enddo

        head_tmp = head
        h_w(:) = head_q0(:)*w_q0(:)
        ! TODO if the normalization against w is necessary?
        head=4.0D0*pi*zdotu(nq0,h_w,1,cmplx(wt_q0_sph(:),0.0D0,8),1)/norm_w_q0
        if(ldbg)then
            write(*,"(A4,2f13.4,A4,2f13.4)") "OH", head_tmp,"AH",head
            write(fid_outgw,"(A6,I3,A2,2e13.4)") "e-100",iomega,"A",head
        endif
        ! project head*w, q0va, q0vb on spherical harmonics
        do ilm=1,lmsq
            h_lm(ilm)=zdotu(nq0,h_w*conjg(sph_harms(ilm,:)),1,cmplx(wt_q0_sph(:),0.0D0,8),1) &
     &        /norm_w_q0 * 4.0D0 * pi
            do im=1,matsiz
                a_lm(ilm,im)=zdotu(nq0,q0_va(:,im)*conjg(sph_harms(ilm,:)),1,cmplx(wt_q0_sph(:),0.0D0,8),1)/norm_w_q0
                b_lm(ilm,im)=zdotu(nq0,q0_vb(:,im)*conjg(sph_harms(ilm,:)),1,cmplx(wt_q0_sph(:),0.0D0,8),1)/norm_w_q0
            enddo
        enddo

        do im=1,matsiz
            wv(im)=zdotu(nq0,q0_va(:,im)*h_w(:),1,cmplx(wt_q0_sph(:),0.0D0,8),1)/norm_w_q0
            wh(im)=zdotu(nq0,q0_vb(:,im)*h_w(:),1,cmplx(wt_q0_sph(:),0.0D0,8),1)/norm_w_q0
            do jm=1,matsiz
                bodyinv_tmp = bodyinv(im,jm)
                do i=1,3
                    do j=1,3
                        ten_aob_tmp(i,j) = vec_a_ani(i,im,iomega)*vec_b_ani(j,jm,iomega)
                    enddo
                enddo
                do iq0=1,nq0
                    q_aob_q(iq0) = ten_rvctrv(3, ten_aob_tmp, q0_sph(iq0,:))
                enddo
                bodyinv(im,jm) = bodyinv(im,jm)+ cmplx(4.0D0*pi,0.0D0,8)* &
     &              zdotu(nq0,h_w*q_aob_q,1,cmplx(wt_q0_sph(:),0.0D0,8),1)/norm_w_q0
                if(ldbg) write(fid_outdbg, "(A10,I3,I4,I4,2E13.5)") "diffbody:",iomega,im,jm,bodyinv(im,jm)-bodyinv_tmp
            enddo
        enddo

        ! Use expansion on spherical harmonics. Problem exists!!!
        ! clean wings first
        !wv = czero
        !wh = czero
        !do l1=0,lmax_q0
        !    do m1=-l1,l1
        !        lm1 = 1 + l1*(l1+1) + m1
        !        do l2=0,lmax_q0
        !            if(l2.ne.l1) continue
        !            do m2=-l2,l2
        !                if(m2+m1.ne.0) continue
        !                lm2 = 1 + l2*(l2+1) + m2
        !                do im=1,matsiz
        !                    wv(im) = wv(im) + h_lm(lm1) * a_lm(lm2,im) * (-1)**m2
        !                    wh(im) = wh(im) + h_lm(lm1) * b_lm(lm2,im) * (-1)**m2
        !                enddo
        !            enddo
        !        enddo
        !    enddo
        !enddo
        wv = - wv *sqrt(4.0D0*pi)
        wh = - wh *sqrt(4.0D0*pi)

        deallocate(sph_harms, h_w, a_lm, b_lm, q_aob_q)

        END SUBROUTINE angint_eps_sph

END MODULE ANISOTROPY

