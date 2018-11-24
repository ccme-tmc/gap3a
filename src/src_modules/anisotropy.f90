MODULE ANISOTROPY

    ! this module defines variables and subroutines utilized for dealing
    ! with the anisotropy of dielectric matrix for $q\to0$

    use mixbasis, only: matsiz
    use constants, only: czero,cone,pi
    use task, only: fid_outgw, fid_outdbg, fid_aniso

    implicit none
    integer :: iop_aniso = -1  ! control whether to consider the anisotropy of dielectric function around Gamma
                               ! -1 -- use q0_eps, no anisotropy

    complex(8), allocatable :: vec_u_ani(:,:,:)   ! vector u
    complex(8), allocatable :: vec_t_ani(:,:,:)   ! vector t
    complex(8), allocatable :: vec_a_ani(:,:,:)   ! vector a
    complex(8), allocatable :: vec_b_ani(:,:,:)   ! vector b
    complex(8), allocatable :: ten_p_ani(:,:,:)   ! tensor P
    complex(8), allocatable :: ten_a_ani(:,:,:)   ! tensor A
    ! Following quantities should be determined by a Lebedev-Laikov grid
    !integer :: nq0 = 0                   ! number of q0 for angular integration
    !real(8), allocatable :: q0_sph(:,:)  ! direction of q0, similar to q0_eps
    !real(8), allocatable :: wt_q0_sph(:) ! weight of q0_sph
    integer :: lmax_gamma = 4               ! maximum angular momentum to expand eps on q0
    integer :: lmgsq


    complex(8), allocatable :: head_g(:,:)      ! the head (n_ang_grid,nomega) for q0 at some freq
    complex(8), allocatable :: wv_g(:,:,:)      ! vertical wing for q0
    complex(8), allocatable :: wh_g(:,:,:)      ! horizontal wing for q0
    complex(8), allocatable :: sph_harms_g(:,:) ! Y_lm
    complex(8), allocatable :: h_g_lm(:,:)      ! projection of head_g on Y_lm
    complex(8), allocatable :: qmax_g_lm(:)   ! projection of qmax_gamma on Y_lm

    ! smallq define the proximity around Gamma point
    !real(8) :: smallq(26,3)
    !integer :: smallq_div = 1
    !character(len=3) :: smallq_type = "bzl"  ! BZ like ('bzl')or spherical ('sph')
    !real(8) :: vol_q0                         ! volume of the defined proximity near Gamma
    real(8),allocatable :: w_gamma(:)  ! previous w_q0
    real(8) :: norm_w_gamma ! previous norm_w_q0


    CONTAINS

        SUBROUTINE init_aniso(iq, iomfirst, iomlast)

        use bzinteg, only: grid_vec,n_ang_grid,vol_gamma,qmax_gamma,ang_weight
        implicit none
        
        integer,intent(in) :: iq 
        integer,intent(in) :: iomfirst, iomlast
        integer :: iang  ! Counter over n_ang_grid
        integer :: ilm   ! Counter over lmgsq
        integer :: ierr

        external :: ylm
        complex(8),external :: zdotc


        lmgsq = (lmax_gamma+1)**2
        if(iq.eq.1)then
            if(iop_aniso.ne.-1)then
                write(fid_outgw,*) "Anisotropy switched on"
                ! TODO avoid body_q0 in future for integration 
                allocate(vec_u_ani(3,matsiz,iomfirst:iomlast),      &
                &        vec_t_ani(3,matsiz,iomfirst:iomlast),      &
                &        vec_a_ani(3,matsiz,iomfirst:iomlast),      &
                &        vec_b_ani(3,matsiz,iomfirst:iomlast),      &
                &        ten_p_ani(3,3,iomfirst:iomlast),           &
                &        ten_a_ani(3,3,iomfirst:iomlast),           &
            !    &        q0_sph(1:nq0,3),                          &
            !    &        wt_q0_sph(1:nq0),                         &
                &        head_g(1:n_ang_grid,iomfirst:iomlast),     &
                &        wv_g(1:n_ang_grid,matsiz,iomfirst:iomlast),&
                &        wh_g(1:n_ang_grid,matsiz,iomfirst:iomlast),&
            !    &        qmax_q0(1:nq0),                           &
                &        w_gamma(1:n_ang_grid),                        &
                &        sph_harms_g(lmgsq,n_ang_grid),             &
                &        h_g_lm(lmgsq,iomfirst:iomlast),            &
                &        qmax_g_lm(lmgsq),                          &
                &        stat=ierr)
                if(ierr.ne.0) then
                    write(fid_outgw,*) " - init_aniso: Fail to allocate aniso"
                    stop
                else
                    write(fid_outgw,*) " - init_aniso: success"
                endif
            endif

            ten_p_ani = czero
            ten_a_ani = czero
            vec_u_ani = czero
            vec_t_ani = czero
            vec_a_ani = czero
            vec_b_ani = czero
            head_g = czero
            wv_g = czero
            wh_g = czero
            w_gamma = qmax_gamma**3/vol_gamma/3.0
            norm_w_gamma = sum(w_gamma*ang_weight)

            ! initialize spherical harmonics and projection of qmax_gamma
            do iang=1,n_ang_grid
              call ylm(grid_vec(iang,:),lmax_gamma,sph_harms_g(:,iang))
              ! include the weight of angular grid in sph_harms
              sph_harms_g(:,iang) = sph_harms_g(:,iang)*cmplx(ang_weight(iang),0.0D0,8)
            enddo

            do ilm=1,lmgsq
              qmax_g_lm(ilm)=zdotc(n_ang_grid,sph_harms_g(ilm,:),1,qmax_gamma(:),1)
            enddo

        endif ! iq.eq.1

        END SUBROUTINE init_aniso


        SUBROUTINE end_aniso(iq)

        integer,intent(in) :: iq
        if(iq.eq.1)then
            deallocate(vec_u_ani, vec_a_ani, vec_b_ani, vec_t_ani, &
    &                  ten_p_ani, ten_a_ani, &
    &                  head_g, wh_g, wv_g, w_gamma, h_g_lm, sph_harms_g )
!    &                  q0_sph, wt_q0_sph, &
!    &                  qmax_q0)
            write(fid_outgw,*) " - end_aniso: success"
        endif

        END SUBROUTINE end_aniso


!        SUBROUTINE init_smallq
!        ! adopted from set_kmax_gama in bzinteg
!        use kpoints,  only: nkdivs
!        use struk,    only: br2, vi
!        implicit none
!! !LOCAL VARIABLES:
!        integer :: i1, i2, i3
!        integer :: j
!        integer :: isq,iq0
!        real(8) :: qmax_tmp, denominator, numerator
!        logical :: ldbg = .true.
!
!        real(8),external :: ddot
!!EOP
!!BOC
!        !if (smallq_div.le.0) then
!        !  write(fid_outgw,*) " - init_smallq: illegal smallq_div"
!        !  stop
!        !endif
!        !! determine the q-points closet to gamma
!        isq=0
!        do i1=-1,1
!          do i2=-1,1
!            do i3=-1,1
!              if(.not.((i1 .eq. 0) .and. (i2 .eq. 0) .and. (i3 .eq. 0)))then
!                isq=isq+1
!                do j=1,3
!!                  smallq(isq,j)=dble(i1)*br2(j,1)/dble(nkdivs(1)*smallq_div)+&
!!     &                          dble(i2)*br2(j,2)/dble(nkdivs(2)*smallq_div)+&
!!     &                          dble(i3)*br2(j,3)/dble(nkdivs(3)*smallq_div)
!                  smallq(isq,j)=dble(i1)*br2(j,1)/dble(nkdivs(1))+&
!     &                          dble(i2)*br2(j,2)/dble(nkdivs(2))+&
!     &                          dble(i3)*br2(j,3)/dble(nkdivs(3))
!                enddo ! j
!              endif
!            enddo ! i3
!          enddo ! i2
!        enddo ! i1
!
!        !! determine qmax with initialized by a large value 1000
!        qmax_q0 = 1.0D3
!
!        ! calculate qmax along q0_sph in the region defined by smallq
!        do iq0=1,nq0
!          do isq=1,26
!            denominator = sum(q0_sph(iq0,:)*smallq(isq,:))
!            if(denominator .gt. 1.0d-10)then
!              numerator = 0.5d0 * sum(smallq(isq,1:3)**2)
!              qmax_tmp = numerator/denominator
!              if(qmax_tmp .lt. qmax_q0(iq0)) qmax_q0(iq0) = qmax_tmp
!            endif
!          enddo ! isk
!        enddo ! iq0
!        !vol_q0 = vi*(2.0D0*pi)**3/dble(product(nkdivs)*smallq_div**3)
!        vol_q0 = vi*(2.0D0*pi)**3/dble(product(nkdivs))
!
!        w_gamma(:) = qmax_q0(:)**3/3.0D0/vol_q0
!        ! normalization of w_gamma(q) in the defined proximity
!        norm_w_q0 = ddot(nq0,w_q0,1,wt_q0_sph,1)
!
!        !if(ldbg) then
!        write(fid_aniso,"(A15)") "smallq list:"
!        do isq=1,26
!            write(fid_aniso,"(I3,3F12.6)") isq, smallq(isq,:)
!        enddo
!        write(fid_outgw,"(A40, f12.6)") "Volume of Gamma proximity (a.u.^-3): ", vol_q0
!        write(fid_outgw,"(A40,I4)") "nq0: ", nq0
!        write(fid_aniso,"(A40,I4)") "qmax with nq0: ", nq0
!        do iq0=1,nq0
!            write(fid_aniso,"(I4,4F12.6)") iq0, q0_sph(iq0,:),qmax_q0(iq0)
!        enddo
!        write(fid_outgw,"(A25,F12.7)") "Normalization of w(q) = ", norm_w_q0
!        !endif
!
!        END SUBROUTINE init_smallq


        SUBROUTINE calc_h_w_inv_ang_grid(iom)
        ! calculate the head and wings of the inverse of dielectric matrix
        use bzinteg, only: n_ang_grid,grid_vec

        implicit none
        integer,intent(in) :: iom
        integer :: iang    ! Counter: runs over n_ang_grid
        integer :: im      ! Counter: runs over matsiz
        complex(8) :: ccoefcoul_g

        complex(8),external :: ten_rvctrv

        ccoefcoul_g=cmplx(4.0D0*pi,0.0D0,8)

        ! vertical wing
        call zgemm('n','n',n_ang_grid,matsiz,3,-sqrt(ccoefcoul_g), &
     &      cmplx(grid_vec,0.0D0,8),n_ang_grid, vec_a_ani(:,:,iom),3, &
     &      czero,wv_g(:,:,iom),n_ang_grid)
        ! horizontal wing
        call zgemm('n','n',n_ang_grid,matsiz,3,-sqrt(ccoefcoul_g), &
     &      cmplx(grid_vec,0.0D0,8),n_ang_grid, vec_b_ani(:,:,iom),3, &
     &      czero,wh_g(:,:,iom),n_ang_grid)   

        do iang=1, n_ang_grid
        ! head
          head_g(iang,iom) = cone / &
     &     (cone+ccoefcoul_g*ten_rvctrv(3,ten_a_ani(:,:,iom),grid_vec(iang,:)))
        ! wings
          wv_g(iang,:,iom) = wv_g(iang,:,iom) * head_g(iang,iom)
          wh_g(iang,:,iom) = wh_g(iang,:,iom) * head_g(iang,iom)
        enddo


        END SUBROUTINE calc_h_w_inv_ang_grid


        SUBROUTINE proj_head_on_ylm(iom)
        use bzinteg, only: n_ang_grid

        implicit none
        integer,intent(in) :: iom
        integer :: ilm     ! Counter: runs over lmgsq
        complex(8),external :: zdotc

        do ilm=1,lmgsq
          h_g_lm(ilm,iom)=zdotc(n_ang_grid,sph_harms_g(ilm,:),1,head_g(:,iom),1)
        enddo

        END SUBROUTINE proj_head_on_ylm
        

        SUBROUTINE angint_eps_sph(iom, bodyinv, use_harm)
        ! calculate the anisotropic term in the body of the inverse of dielectric matrix
        ! by direct calculation, or by the use of spherical harmonics
        use bzinteg, only: n_ang_grid, ang_weight

        implicit none
        integer,intent(in) :: iom  ! index of frequency
        complex(8),intent(inout) :: bodyinv(matsiz,matsiz)
        logical :: use_harm        ! flag to use expansion on spherical harmonics
        ! the direct inverse of the body of the dielectric matrix
!
!        ! local variables
        integer :: iang
        complex(8) :: cw
        integer :: im,jm             ! Counter: runs over matsiz
!        complex :: head_tmp, bodyinv_tmp
!        complex(8),allocatable :: a_lm(:,:),b_lm(:,:)
!        complex(8),allocatable :: h_w(:)  ! head of invers, times w_gamma
!        !complex(8),allocatable :: q_aob_q(:) 
!        integer :: iq0               ! Counter: runs over nq0
!        integer :: ilm,lm1,lm2,lm3   ! Counter: runs over lmgsq
!        integer :: l1,l2,l3,m1,m2,m3 ! Counters: runs over angular moment
!        integer :: i,j               ! Counters: Cartesian axis
!        integer :: scheme_head = 2
!        logical :: ldbg=.true.
!        complex(8) :: ccoefcoul_q0, ten_aob_tmp(3,3)
!
        external ylm
        complex(8),external :: zdotu,zdotc

        if(use_harm)then
        else
          ! TODO maybe need optimize
          do iang=1,n_ang_grid
            cw = cmplx(w_gamma(iang)*ang_weight(iang),0.0D0,8)
            do im=1,matsiz
              do jm=1,matsiz
!                bodyinv(im,jm) = bodyinv(im,jm) + &
!     &              cw*wv_g(iang,im,iom)*wh_g(iang,jm,iom)/head_g(iang,iom)
              enddo
            enddo
          enddo
        endif
!
!        allocate(sph_harms(lmgsq,nq0), &
!     &           h_w(nq0),            &
!!    &           q_aob_q(nq0),        &
!     &           qmax_g_lm(lmgsq),          &
!     &           a_lm(lmgsq,matsiz),   &
!     &           b_lm(lmgsq,matsiz)    &
!     &          )
!        ! calculate head_g and q0\cdot wing1, q0\cdot wing2
!        do iq0=1,nq0
!          ccoefcoul_q0=cmplx(4.0D0*pi,0.0D0,8)
!          ! head
!          head_g(iq0,iomega) = cone / &
!     &     (cone+ccoefcoul_q0*ten_rvctrv(3,ten_a_ani(:,:,iomega),q0_sph(iq0,:)))
!          ! wings
!          do im=1,matsiz
!          ! TODO optimize with ZGEMM
!            q0_va(iq0,im)=sum(vec_a_ani(:,im,iomega)*cmplx(q0_sph(iq0,:),0.0D0,8))
!            q0_vb(iq0,im)=sum(vec_b_ani(:,im,iomega)*cmplx(q0_sph(iq0,:),0.0D0,8))
!          enddo
!        enddo
!
!        ! calculate spherical harmonics at q0_sph
!        do iq0=1,nq0
!            call ylm(q0_sph(iq0,:),lmax_gamma,sph_harms(:,iq0))
!            ! check the calculation
!            !if(ldbg)then
!            !    write(*,"(A20,I5,A1,3F10.4,A1)") "Ylm at iq0 = ",iq0, "(", q0_sph(iq0,:), ")"
!            !    do l1=0,lmax_gamma
!            !        do m1=1,2*l1+1
!            !            write(*,"(A3,I2,A3,I2,2F13.4)") "l=",l1,"m=",m1-l1-1,sph_harms(l1**2+m1,iq0)
!            !        enddo
!            !    enddo
!            !endif
!        enddo
!
!        h_w(:) = head_g(:,iomega)*w_q0(:)
!
!        do ilm=1,lmgsq
!        ! project head on spherical harmonics
!            h_g_lm(ilm)=zdotu(nq0,head_g(:,iomega)*conjg(sph_harms(ilm,:)),1,cmplx(wt_q0_sph(:),0.0D0,8),1)!&
!!     &      /norm_w_q0
!            qmax_g_lm(ilm)=zdotu(nq0,conjg(sph_harms(ilm,:)),1,cmplx(qmax_q0(:)*wt_q0_sph(:),0.0D0,8),1)
!        ! project head*w, q0va, q0vb on spherical harmonics
!!            h_g_lm(ilm)=zdotu(nq0,h_w*conjg(sph_harms(ilm,:)),1,cmplx(wt_q0_sph(:),0.0D0,8),1) &
!!     &        /norm_w_q0 * 4.0D0 * pi
!!            do im=1,matsiz
!!                a_lm(ilm,im)=zdotu(nq0,q0_va(:,im)*conjg(sph_harms(ilm,:)),1,cmplx(wt_q0_sph(:),0.0D0,8),1)/norm_w_q0
!!                b_lm(ilm,im)=zdotu(nq0,q0_vb(:,im)*conjg(sph_harms(ilm,:)),1,cmplx(wt_q0_sph(:),0.0D0,8),1)/norm_w_q0
!!            enddo
!        enddo
!
!        write(fid_aniso,"(I8,A26,4E13.4)") iomega, "head_g_eps and Y00 part", head, h_g_lm(1)*sph_harms(1,1)
!        ! check projection of head
!        do iq0=1,nq0
!            head_tmp=sum(sph_harms(:,iq0)*h_g_lm)
!            write(fid_aniso,"(A12,I5,2F13.4,A12,2F13.4)") "head_g",iq0,head_g(iq0,iomega)," sum_Ylm = ",head_tmp
!            write(fid_aniso,"(I8,I9,A43,2E13.4)") iomega, iq0, "expand_Ylm_Diff = ", head_g(iq0,iomega)-head_tmp
!        enddo
!
!        head_tmp = head
!        if(scheme_head.eq.1)then
!            ! Average head in Gamma proximity
!            head=zdotu(nq0,h_w,1,cmplx(wt_q0_sph(:),0.0D0,8),1)/norm_w_q0
!            write(fid_aniso,"(A7,2f13.4,A7,2f13.4)") "Old H", head_tmp,"Ave H",head
!            write(fid_aniso,"(A20,I3,2e13.4)") "Ang. Ave. e-100",iomega,head
!        elseif(scheme_head.eq.2)then
!            ! Use h_00 term only
!            ! according to Eq.(45) in Friedrich, et al PRB 81,125102(2010)
!            head = h_g_lm(1)*sph_harms(1,1)
!            write(fid_aniso,"(A7,2f13.4,A7,2f13.4)") "Old H ", head_tmp,"H00Y00", head
!        elseif(scheme_head.eq.3)then
!            ! according to Eq.(36) in Freysoldt, et al CPC 176,1(2007)
!            head = zdotc(lmgsq, qmax_g_lm, 1, h_g_lm, 1) / vol_q0
!            write(fid_aniso,"(I3,A7,2f13.4,A7,2f13.4)") iomega,"Old H ",head_tmp," SumQH ",head
!        endif
!
!        do im=1,matsiz
!        ! TODO use BLAS-2 routines
!            wv(im)=zdotu(nq0,q0_va(:,im)*h_w(:),1,cmplx(wt_q0_sph(:),0.0D0,8),1)/norm_w_q0
!            wh(im)=zdotu(nq0,q0_vb(:,im)*h_w(:),1,cmplx(wt_q0_sph(:),0.0D0,8),1)/norm_w_q0
!            do jm=1,matsiz
!                bodyinv_tmp = bodyinv(im,jm)
!!                do i=1,3
!!                    do j=1,3
!!                        ten_aob_tmp(i,j) = vec_a_ani(i,im,iomega)*vec_b_ani(j,jm,iomega)
!!                    enddo
!!                enddo
!!                do iq0=1,nq0
!!                    q_aob_q(iq0) = ten_rvctrv(3, ten_aob_tmp, q0_sph(iq0,:))
!!                enddo
!!                bodyinv(im,jm) = bodyinv(im,jm)+ cmplx(4.0D0*pi,0.0D0,8)* &
!!     &              zdotu(nq0,h_w*q_aob_q,1,cmplx(wt_q0_sph(:),0.0D0,8),1)/norm_w_q0
!                ! actually there is no need to calculate q_aob_q, just use q0_va and q0_vb
!                bodyinv(im,jm) = bodyinv(im,jm) + cmplx(4.0D0*pi,0.0D0,8)* &
!     &              zdotu(nq0,h_w*q0_va(:,im)*q0_vb(:,jm),1,cmplx(wt_q0_sph(:),0.0D0,8),1)/norm_w_q0
!                if(ldbg) write(fid_aniso, "(A10,I3,I4,I4,2E13.5)") "diffbody:",iomega,im,jm,bodyinv(im,jm)-bodyinv_tmp
!                ! TODO verify body, and comment the line below
!                bodyinv(im,jm) = bodyinv_tmp
!            enddo
!        enddo
!
!        ! TODO Use expansion on spherical harmonics. Problem exists!!!
!        ! clean wings first
!        !wv = czero
!        !wh = czero
!        !do l1=0,lmax_gamma
!        !    do m1=-l1,l1
!        !        lm1 = 1 + l1*(l1+1) + m1
!        !        do l2=0,lmax_gamma
!        !            if(l2.ne.l1) continue
!        !            do m2=-l2,l2
!        !                if(m2+m1.ne.0) continue
!        !                lm2 = 1 + l2*(l2+1) + m2
!        !                do im=1,matsiz
!        !                    wv(im) = wv(im) + h_g_lm(lm1) * a_lm(lm2,im) * (-1)**m2
!        !                    wh(im) = wh(im) + h_g_lm(lm1) * b_lm(lm2,im) * (-1)**m2
!        !                enddo
!        !            enddo
!        !        enddo
!        !    enddo
!        !enddo
!        wv = - wv *sqrt(4.0D0*pi)
!        wh = - wh *sqrt(4.0D0*pi)
!
!        deallocate(sph_harms, h_w, a_lm, b_lm, qmax_g_lm)
!
        END SUBROUTINE angint_eps_sph
!
        SUBROUTINE angint_invq2_dhead(iom, angint)
        ! calcualte 
        ! \frac{1}{V_{\Gamma}}\int_{V_{\Gamma}}
        ! {\dd{\hat{\mathbf{q}}}\left\lbrace\varepsilon^{-1}_{00}(\mathbf{q}\to0,\textt{iomega})-1\right\rbrace}
        ! $V_{\Gamma}$ is calculated by 
        integer,intent(in) :: iom
        complex(8),intent(out) :: angint
!        complex(8),external :: zdotu
!
!        angint = zdotu(nq0, head_g(:,iomega)-cone, 1, cmplx(qmax_q0(:)*wt_q0_sph(:),0.0D0,8), 1) &
!     &                / cmplx(norm_w_q0 * vol_q0, 0.0D0, 8)
!        
        END SUBROUTINE angint_invq2_dhead

END MODULE ANISOTROPY

