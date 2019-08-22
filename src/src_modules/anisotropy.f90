module anisotropy
!
! !Description:
!
! this module defines variables and subroutines utilized for dealing
! with the anisotropy of dielectric matrix for $q\to0$
!
! !Uses:
  use mixbasis,  only: matsiz
  use constants, only: czero,cone,pi,ctwopi,cpi,sqrt4pi
  use task,      only: fid_outgw, fid_outdbg, fid_aniso
  use bzinteg,   only: grid_vec,vol_gamma,qmax_gamma,ang_weight, set_qmax_q0
!
! Public variables
  implicit none
  integer :: iop_aniso = -1  ! control whether to consider the anisotropy of dielectric function around Gamma
                             ! -1 -- use q0_eps, no anisotropy
                             !  0 -- 3D anisotropy
                             !  1 -- 2D anisotropy
  integer(4) :: n_ang_grid
  integer :: lmax_gamma = 6                   ! maximum angular momentum to expand eps on q0
  integer :: lmgsq                            ! (lmax_gamma+1)^2

  complex(8), allocatable :: vec_u_ani(:,:,:)   ! vector u
  complex(8), allocatable :: vec_t_ani(:,:,:)   ! vector t
  complex(8), allocatable :: vec_a_ani(:,:,:)   ! vector a
  complex(8), allocatable :: vec_b_ani(:,:,:)   ! vector b
  complex(8), allocatable :: ten_p_ani(:,:,:)   ! tensor P
  complex(8), allocatable :: ten_a_ani(:,:,:)   ! tensor A
  complex(8), allocatable :: head_g(:,:)      ! the head (n_ang_grid,nomega) of eps^-1 at q0
  complex(8), allocatable :: wv_g(:,:,:)      ! vertical wing of eps^-1 at q0
  complex(8), allocatable :: wh_g(:,:,:)      ! horizontal wing of eps^-1 at q0
  complex(8), allocatable :: sph_harms_g(:,:) ! Y_lm at q0
  complex(8), allocatable :: h_g_lm(:,:)      ! projection of head_g on Y_lm
  complex(8), allocatable :: qmax_g_lm(:)     ! projection of qmax_gamma on Y_lm

  ! smallq define the proximity around Gamma point
  !real(8) :: smallq(26,3)
  !integer :: smallq_div = 1
  !character(len=3) :: smallq_type = "bzl"  ! BZ like ('bzl')or spherical ('sph')
  !real(8) :: vol_q0                         ! volume of the defined proximity near Gamma
  real(8),allocatable :: w_gamma(:)  ! previous w_q0
  real(8) :: norm_w_gamma ! previous norm_w_q0
!
! derived type for collectively storing data to calculate dielectric anisotropy
!
  type aniso_tensor
    integer(4) :: sys_dim                    ! dimension of the system, 0D, 1D, 2D, 3D
    integer(4) :: lmax                       ! maximum angular momentum to expand eps on q0
    integer(4) :: lmgsq                      ! (lmax+1)^2
    integer(4) :: nang                       ! number of angular grid points on unit sphere
    integer(4) :: msiz                       ! size of dielectric matrix
    real(8) :: vol_q0                        ! volume of the q0 region
    complex(8),allocatable :: vec_u(:,:,:)   ! vector u
    complex(8),allocatable :: vec_t(:,:,:)   ! vector t
    complex(8),allocatable :: vec_a(:,:,:)   ! vector a
    complex(8),allocatable :: vec_b(:,:,:)   ! vector b
    complex(8),allocatable :: ten_p(:,:,:)   ! tensor P
    complex(8),allocatable :: ten_a(:,:,:)   ! tensor A
    real(8),allocatable    :: q0(:,:)        ! Y_lm at q0
    real(8),allocatable    :: w_q0(:)        ! weight
    real(8),allocatable    :: qmax_q0(:)     ! qmax along angle q0
    real(8),allocatable    :: qmax3d3v(:)    ! qmax^3/3vol_q0
    complex(8),allocatable :: head_q0(:,:)   ! the head (n_ang_grid,nomega) of eps^-1 at q0
    complex(8),allocatable :: wv_q0(:,:,:)   ! vertical wing of eps^-1 at q0
    complex(8),allocatable :: wh_q0(:,:,:)   ! horizontal wing of eps^-1 at q0
    complex(8),allocatable :: ylm_q0(:,:)    ! Y_lm at q0
    complex(8),allocatable :: h_ylm_q0(:,:)  ! projection of head_q0 on Y_lm
    complex(8),allocatable :: qmax_ylm_q0(:) ! projection of qmax_q0 on Y_lm
    real(8)                :: norm           ! normalization factor to vol_q0
  end type aniso_tensor

  type(aniso_tensor), pointer, save :: aniten

  type test_at
    integer(4) :: i
    complex(8),allocatable :: c(:,:,:)
  end type test_at
  type(test_at), pointer :: ta

  private :: cdot_over_ang
  private :: cdot_over_lm

contains

subroutine init_test_at(at, i)
  type(test_at), pointer :: at
  integer(4),intent(in)  :: i
  allocate(at)
  at%i = i
  allocate(at%c(i,i,i))
  at%c = czero
end subroutine init_test_at

subroutine end_test_at(at)
  type(test_at), pointer :: at
  deallocate(at%c)
  deallocate(at)
  nullify(at)
end subroutine end_test_at

subroutine init_aniso_oo(at, iq, msiz, iomfirst, iomlast, lmaxq0, nang)
  use lebedev_laikov
  use modmpi,    only: myrank_ra3
  implicit none
  
  type(aniso_tensor),pointer :: at
  integer(4),intent(in) :: iq 
  integer(4),intent(in) :: msiz
  integer(4),intent(in) :: iomfirst, iomlast
  integer(4),intent(in) :: lmaxq0
  integer(4),intent(in) :: nang

! !Local variables
  integer(4) :: ierr
  integer(4) :: iang ! counter over nang
  integer(4) :: ilm  ! counter over (lmaxq0+1)^2
!
! !external routines
  external linmsg, ylm
  complex(8), external :: zdotu, zdotc

  call linmsg(fid_outgw,'-', "init_aniso_oo")
  allocate(at)
  at%sys_dim = 3 - iop_aniso
  at%lmax = lmaxq0
  at%msiz = msiz
  at%lmgsq = (lmaxq0+1)**2
  at%nang = nang

  if(iq.eq.1)then
    write(fid_outgw,*) "Anisotropy switched on"
    ! TODO avoid body_q0 in future for integration 
    call set_lebedev_laikov_grid(at%nang)
    allocate(at%vec_u(3,msiz,iomfirst:iomlast),       &
             at%vec_t(3,msiz,iomfirst:iomlast),       &
             at%vec_a(3,msiz,iomfirst:iomlast),       &
             at%vec_b(3,msiz,iomfirst:iomlast),       &
             at%ten_p(3,3,iomfirst:iomlast),          &
             at%ten_a(3,3,iomfirst:iomlast),          &
             at%q0(at%nang,3),at%w_q0(at%nang),       &
             at%qmax_q0(at%nang),at%qmax3d3v(at%nang),&
             at%head_q0(at%nang,iomfirst:iomlast),    &
             at%wv_q0(at%nang,msiz,iomfirst:iomlast), &
             at%wh_q0(at%nang,msiz,iomfirst:iomlast), &
             at%ylm_q0(at%lmgsq,at%nang),                &
             at%h_ylm_q0(at%lmgsq,iomfirst:iomlast),     &
             at%qmax_ylm_q0(at%lmgsq),                   &
             stat=ierr)
    if(ierr.ne.0) then
      write(fid_outgw,*) " - fail to allocate aniso_tensor object"
      stop
    else
      write(fid_outgw,*) " - success to allocate aniso_tensor"
    endif
    write(fid_outgw,*) " - initializing angular grid"
    write(fid_outgw,*) " - frequency range: ", iomfirst, iomlast
    at%q0(:,1) = xleb(:)
    at%q0(:,2) = yleb(:)
    at%q0(:,3) = zleb(:)
    at%w_q0(:) = wleb(:)*4.0d0*pi
    call unset_lebedev_laikov_grid
    write(fid_outgw,*) " - computing qmax"
    call set_qmax_q0(at%nang, at%q0, at%qmax_q0, at%vol_q0)
    at%ten_p = czero
    at%ten_a = czero
    at%vec_u = czero
    at%vec_t = czero
    at%vec_a = czero
    at%vec_b = czero
    at%head_q0 = czero
    at%wv_q0 = czero
    at%wh_q0 = czero
    at%qmax3d3v = at%qmax_q0**3/at%vol_q0/3.0
    at%norm = sum(at%qmax3d3v*at%w_q0)

    ! initialize spherical harmonics and projection of qmax_gamma
    write(fid_outgw,*) " - computing Y_lm"
    do iang=1,at%nang
      write(fid_aniso,"(A4,I5,A2,4f10.6)") "Ang ",iang,": ", at%q0(iang,:), at%w_q0(iang)
      call ylm(at%q0(iang,:),at%lmax,at%ylm_q0(:,iang))
      ! include the weight of angular grid in sph_harms
      do ilm=1,at%lmgsq
        write(fid_aniso, "(2I6,2f10.6)") iang, ilm, at%ylm_q0(ilm,iang)
        at%ylm_q0(ilm,iang)=at%ylm_q0(ilm,iang)*cmplx(at%w_q0(iang),0.0D0,8)
        enddo
    enddo

    !write(*,"(A25,F12.5)") "Summation of weight/4pi",sum(ang_weight)/4.0D0/pi
    write(fid_outgw,*) " - computing qmax projection on Y_lm"
    do ilm=1,at%lmgsq
      at%qmax_ylm_q0(ilm)=zdotc(at%nang, at%ylm_q0(ilm,:),1,cmplx(at%qmax_q0(:),0.0D0,8),1)
      if(myrank_ra3.eq.0) write(fid_aniso,"(I4,3F13.5)") ilm, at%qmax_ylm_q0(ilm)
    enddo

    ! check completeness of expansion of qmax
    ! TODO fail to output to fid_outdbg for gap3a-mpi
    ! expansion of qmax is not as accurate as eps
    write(fid_outgw,*) " - check completeness of qmax projection on Y_lm (see aniso debug)"
    if(myrank_ra3.eq.0)then
      do iang=1,n_ang_grid
        write(fid_aniso,"(A4,I6,3F13.6)") "Ang ", iang, at%qmax_q0(iang),&
              zdotu(at%lmgsq,at%ylm_q0(:,iang),1,at%qmax_ylm_q0,1)/cmplx(at%w_q0(iang),0.0D0,8)
      enddo
    endif
  endif ! iq.eq.1
  call linmsg(fid_outgw,'-', "end init_aniso_oo")
end subroutine init_aniso_oo

subroutine end_aniso_oo(iq, at)
  integer(4) :: iq
  type(aniso_tensor), pointer :: at
  if(iq.eq.1) then
    deallocate(at%vec_u,at%vec_t,at%vec_a,at%vec_b,at%ten_p,at%ten_a,&
               at%q0,at%w_q0,at%qmax_q0,at%qmax3d3v,at%head_q0,at%wv_q0,&
               at%wh_q0,at%ylm_q0,at%h_ylm_q0,at%qmax_ylm_q0)
  endif
  deallocate(at)
  nullify(at)
end subroutine end_aniso_oo

subroutine init_aniso(iq, iomfirst, iomlast)

  use modmpi,    only: myrank_ra3
  implicit none
  
  integer,intent(in) :: iq 
  integer,intent(in) :: iomfirst, iomlast
  integer :: iang  ! Counter over n_ang_grid
  integer :: ilm   ! Counter over lmgsq
  integer :: ierr

  external :: ylm
  complex(8),external :: zdotc,zdotu,ddot


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
    !  &        q0_sph(1:nq0,3),                          &
    !  &        wt_q0_sph(1:nq0),                         &
      &        head_g(1:n_ang_grid,iomfirst:iomlast),     &
      &        wv_g(1:n_ang_grid,matsiz,iomfirst:iomlast),&
      &        wh_g(1:n_ang_grid,matsiz,iomfirst:iomlast),&
    !  &        qmax_q0(1:nq0),                           &
      &        w_gamma(1:n_ang_grid),                     &
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
      do ilm=1,lmgsq
        sph_harms_g(ilm,iang)=sph_harms_g(ilm,iang)*cmplx(ang_weight(iang),0.0D0,8)
      enddo
    enddo

    !write(*,"(A25,F12.5)") "Summation of weight/4pi",sum(ang_weight)/4.0D0/pi
    do ilm=1,lmgsq
      qmax_g_lm(ilm)=cdot_over_ang('C',n_ang_grid,sph_harms_g(ilm,:),cmplx(qmax_gamma(:),0.0D0,8))
      if(myrank_ra3.eq.0) write(fid_outdbg,"(I4,3F13.5)") ilm, qmax_g_lm(ilm)
    enddo

    ! check completeness of expansion of qmax
    ! TODO fail to output to fid_outdbg for gap3a-mpi
    ! expansion of qmax is not as accurate as eps
    if(myrank_ra3.eq.0)then
      do iang=1,n_ang_grid
        write(fid_outdbg,"(A4,I6,3F13.6)") "Ang ", iang, qmax_gamma(iang),&
   &          cdot_over_lm('N',lmgsq,sph_harms_g(:,iang),qmax_g_lm)/cmplx(ang_weight(iang),0.0D0,8)
      enddo
    endif

  endif ! iq.eq.1

end subroutine init_aniso


subroutine end_aniso(iq)

  integer,intent(in) :: iq
  if(iq.eq.1)then
    deallocate(vec_u_ani, vec_a_ani, vec_b_ani, vec_t_ani, &
               ten_p_ani, ten_a_ani, &
               head_g, wh_g, wv_g, w_gamma, h_g_lm, sph_harms_g )
!                   q0_sph, wt_q0_sph, &
!                   qmax_q0)
    write(fid_outgw,*) " - end_aniso: success"
  endif

end subroutine end_aniso


subroutine calc_h_w_inv_ang_grid(iom)
! calculate the head and wings of the inverse of dielectric matrix

  implicit none
  integer,intent(in) :: iom
  integer :: iang    ! Counter: runs over n_ang_grid
  integer :: im      ! Counter: runs over matsiz
  complex(8) :: ccoefcoul_g

  complex(8),external :: ten_rvctrv
  external :: zgemm

  ccoefcoul_g=cmplx(4.0D0*pi,0.0D0,8)

  ! vertical wing
  call zgemm('n','n',n_ang_grid,matsiz,3,-sqrt(ccoefcoul_g), &
             cmplx(grid_vec,0.0D0,8),n_ang_grid, vec_a_ani(:,:,iom),3, &
             czero,wv_g(:,:,iom),n_ang_grid)
  ! horizontal wing
  call zgemm('n','n',n_ang_grid,matsiz,3,-sqrt(ccoefcoul_g), &
             cmplx(grid_vec,0.0D0,8),n_ang_grid, vec_b_ani(:,:,iom),3, &
             czero,wh_g(:,:,iom),n_ang_grid)   

  do iang=1, n_ang_grid
  ! head
    head_g(iang,iom) = cone / &
     (cone+ccoefcoul_g*ten_rvctrv(3,ten_a_ani(:,:,iom),grid_vec(iang,:)))
  ! wings
    wv_g(iang,:,iom) = wv_g(iang,:,iom) * head_g(iang,iom)
    wh_g(iang,:,iom) = wh_g(iang,:,iom) * head_g(iang,iom)
  enddo

end subroutine calc_h_w_inv_ang_grid


subroutine calc_h_w_inv_ang_grid_oo(at, iom)
! calculate the head and wings of the inverse of dielectric matrix

  implicit none
  type(aniso_tensor), pointer :: at
  integer,intent(in) :: iom
! !Local variables
  integer :: iang    ! Counter: runs over n_ang_grid
  integer :: msiz
  integer :: nang
  integer :: im      ! Counter: runs over msiz
  complex(8) :: ccoefcoul_g

  complex(8),external :: ten_rvctrv
  external :: zgemm

  ccoefcoul_g=cmplx(4.0D0*pi,0.0D0,8)
  msiz = at%msiz
  nang = at%nang

  ! vertical wing
  call zgemm('n','n',nang,msiz,3,-sqrt(ccoefcoul_g), &
             cmplx(at%q0,0.0D0,8),nang, at%vec_a(:,:,iom),3, &
             czero,at%wv_q0(:,:,iom),nang)
  ! horizontal wing
  call zgemm('n','n',nang,msiz,3,-sqrt(ccoefcoul_g), &
             cmplx(at%q0,0.0D0,8),nang, at%vec_b(:,:,iom),3, &
             czero,at%wh_q0(:,:,iom),nang)

  do iang=1, nang
  ! head
    at%head_q0(iang,iom) = cone / &
     (cone+ccoefcoul_g*ten_rvctrv(3,at%ten_a(:,:,iom),at%q0(iang,:)))
  ! wings
    at%wv_q0(iang,:,iom) = at%wv_q0(iang,:,iom) * at%head_q0(iang,iom)
    at%wh_q0(iang,:,iom) = at%wh_q0(iang,:,iom) * at%head_q0(iang,iom)
  enddo

end subroutine calc_h_w_inv_ang_grid_oo


subroutine proj_head_on_ylm(iom)

  implicit none
  integer,intent(in) :: iom
  integer :: ilm     ! Counter: runs over lmgsq
  integer :: iang    ! Counter: runs over n_ang_grid
  complex(8) :: head_g_tmp

  do ilm=1,lmgsq
    h_g_lm(ilm,iom)=cdot_over_ang('C',n_ang_grid,sph_harms_g(ilm,:),head_g(:,iom))
  enddo

  ! check completeness of expansion
  write(fid_outdbg,"(A40,I3)") "Completeness of Ylm expansion of head_g, iom ", iom
  do iang=1,n_ang_grid
    write(fid_outdbg,"(I5,4F14.6)") iang, head_g(iang,iom), &
          cdot_over_lm('N',lmgsq,sph_harms_g(:,iang),h_g_lm(:,iom))/cmplx(ang_weight(iang),0.0D0,8)
  enddo

  ! substract \sqrt{4\pi} for l=0,m=0 term to exclude bare Coulomb part
  h_g_lm(1,iom) = h_g_lm(1,iom) - cmplx(sqrt4pi,0.0D0,8)

end subroutine proj_head_on_ylm


subroutine proj_head_on_ylm_oo(at, iom)

  implicit none
  type(aniso_tensor), pointer :: at
  integer,intent(in) :: iom
  integer :: nlm, nang
  integer :: ilm     ! Counter: runs over lmgsq
  integer :: iang    ! Counter: runs over n_ang_grid
  complex(8) :: head_g_tmp

  nang = at%nang
  nlm = at%lmgsq

  do ilm=1,lmgsq
    at%h_ylm_q0(ilm,iom)=cdot_over_ang('C',nang,at%ylm_q0(ilm,:),at%head_q0(:,iom))
  enddo

  ! check completeness of expansion
  write(fid_outdbg,"(A40,I3)") "Completeness of Ylm expansion of head_g, iom ", iom
  do iang=1,nang
    write(fid_outdbg,"(I5,4F14.6)") iang, at%head_q0(iang,iom), &
          cdot_over_lm('N',nlm,at%ylm_q0(:,iang),at%h_ylm_q0(:,iom))/cmplx(at%w_q0(iang),0.0D0,8)
  enddo

  ! substract \sqrt{4\pi} for l=0,m=0 term to exclude bare Coulomb part
  at%h_ylm_q0(1,iom) = at%h_ylm_q0(1,iom) - cmplx(sqrt4pi,0.0D0,8)

end subroutine proj_head_on_ylm_oo
      


SUBROUTINE angint_eps_sph(iom, bodyinv, use_harm)
! calculate the anisotropic term in the body of the inverse of dielectric matrix
! by direct calculation, or by the use of spherical harmonics

  implicit none
  integer,intent(in) :: iom  ! index of frequency
  complex(8),intent(inout) :: bodyinv(matsiz,matsiz)
  logical :: use_harm        ! flag to use expansion on spherical harmonics
  ! the direct inverse of the body of the dielectric matrix
!
!  ! local variables
  integer :: iang
  complex(8) :: cw
  integer :: im,jm             ! Counter: runs over matsiz
!  complex :: head_tmp, bodyinv_tmp
!  complex(8),allocatable :: a_lm(:,:),b_lm(:,:)
!  complex(8),allocatable :: h_w(:)  ! head of invers, times w_gamma
!  !complex(8),allocatable :: q_aob_q(:) 
!  integer :: iq0               ! Counter: runs over nq0
!  integer :: ilm,lm1,lm2,lm3   ! Counter: runs over lmgsq
!  integer :: l1,l2,l3,m1,m2,m3 ! Counters: runs over angular moment
!  integer :: i,j               ! Counters: Cartesian axis
!  integer :: scheme_head = 2
!  logical :: ldbg=.true.
!  complex(8) :: ccoefcoul_q0, ten_aob_tmp(3,3)
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
!          bodyinv(im,jm) = bodyinv(im,jm) + &
!              cw*wv_g(iang,im,iom)*wh_g(iang,jm,iom)/head_g(iang,iom)
        enddo
      enddo
    enddo
  endif
!
!  h_w(:) = head_g(:,iomega)*w_q0(:)
!
!  do ilm=1,lmgsq
!  ! project head on spherical harmonics
!      h_g_lm(ilm)=zdotu(nq0,head_g(:,iomega)*conjg(sph_harms(ilm,:)),1,cmplx(wt_q0_sph(:),0.0D0,8),1)!&
!!      /norm_w_q0
!      qmax_g_lm(ilm)=zdotu(nq0,conjg(sph_harms(ilm,:)),1,cmplx(qmax_q0(:)*wt_q0_sph(:),0.0D0,8),1)
!  ! project head*w, q0va, q0vb on spherical harmonics
!!      h_g_lm(ilm)=zdotu(nq0,h_w*conjg(sph_harms(ilm,:)),1,cmplx(wt_q0_sph(:),0.0D0,8),1) &
!!        /norm_w_q0 * 4.0D0 * pi
!!      do im=1,matsiz
!!          a_lm(ilm,im)=zdotu(nq0,q0_va(:,im)*conjg(sph_harms(ilm,:)),1,cmplx(wt_q0_sph(:),0.0D0,8),1)/norm_w_q0
!!          b_lm(ilm,im)=zdotu(nq0,q0_vb(:,im)*conjg(sph_harms(ilm,:)),1,cmplx(wt_q0_sph(:),0.0D0,8),1)/norm_w_q0
!!      enddo
!  enddo
!
!  write(fid_aniso,"(I8,A26,4E13.4)") iomega, "head_g_eps and Y00 part", head, h_g_lm(1)*sph_harms(1,1)
!  ! check projection of head
!  do iq0=1,nq0
!      head_tmp=sum(sph_harms(:,iq0)*h_g_lm)
!      write(fid_aniso,"(A12,I5,2F13.4,A12,2F13.4)") "head_g",iq0,head_g(iq0,iomega)," sum_Ylm = ",head_tmp
!      write(fid_aniso,"(I8,I9,A43,2E13.4)") iomega, iq0, "expand_Ylm_Diff = ", head_g(iq0,iomega)-head_tmp
!  enddo
!
!  head_tmp = head
!  if(scheme_head.eq.1)then
!      ! Average head in Gamma proximity
!      head=zdotu(nq0,h_w,1,cmplx(wt_q0_sph(:),0.0D0,8),1)/norm_w_q0
!      write(fid_aniso,"(A7,2f13.4,A7,2f13.4)") "Old H", head_tmp,"Ave H",head
!      write(fid_aniso,"(A20,I3,2e13.4)") "Ang. Ave. e-100",iomega,head
!  elseif(scheme_head.eq.2)then
!      ! Use h_00 term only
!      ! according to Eq.(45) in Friedrich, et al PRB 81,125102(2010)
!      head = h_g_lm(1)*sph_harms(1,1)
!      write(fid_aniso,"(A7,2f13.4,A7,2f13.4)") "Old H ", head_tmp,"H00Y00", head
!  elseif(scheme_head.eq.3)then
!      ! according to Eq.(36) in Freysoldt, et al CPC 176,1(2007)
!      head = zdotc(lmgsq, qmax_g_lm, 1, h_g_lm, 1) / vol_q0
!      write(fid_aniso,"(I3,A7,2f13.4,A7,2f13.4)") iomega,"Old H ",head_tmp," SumQH ",head
!  endif
!
!  do im=1,matsiz
!  ! TODO use BLAS-2 routines
!      wv(im)=zdotu(nq0,q0_va(:,im)*h_w(:),1,cmplx(wt_q0_sph(:),0.0D0,8),1)/norm_w_q0
!      wh(im)=zdotu(nq0,q0_vb(:,im)*h_w(:),1,cmplx(wt_q0_sph(:),0.0D0,8),1)/norm_w_q0
!      do jm=1,matsiz
!          bodyinv_tmp = bodyinv(im,jm)
!!          do i=1,3
!!              do j=1,3
!!                  ten_aob_tmp(i,j) = vec_a_ani(i,im,iomega)*vec_b_ani(j,jm,iomega)
!!              enddo
!!          enddo
!!          do iq0=1,nq0
!!              q_aob_q(iq0) = ten_rvctrv(3, ten_aob_tmp, q0_sph(iq0,:))
!!          enddo
!!          bodyinv(im,jm) = bodyinv(im,jm)+ cmplx(4.0D0*pi,0.0D0,8)* &
!!              zdotu(nq0,h_w*q_aob_q,1,cmplx(wt_q0_sph(:),0.0D0,8),1)/norm_w_q0
!          ! actually there is no need to calculate q_aob_q, just use q0_va and q0_vb
!          bodyinv(im,jm) = bodyinv(im,jm) + cmplx(4.0D0*pi,0.0D0,8)* &
!              zdotu(nq0,h_w*q0_va(:,im)*q0_vb(:,jm),1,cmplx(wt_q0_sph(:),0.0D0,8),1)/norm_w_q0
!          if(ldbg) write(fid_aniso, "(A10,I3,I4,I4,2E13.5)") "diffbody:",iomega,im,jm,bodyinv(im,jm)-bodyinv_tmp
!          ! TODO verify body, and comment the line below
!          bodyinv(im,jm) = bodyinv_tmp
!      enddo
!  enddo
!
!  ! TODO Use expansion on spherical harmonics. Problem exists!!!
!  ! clean wings first
!  !wv = czero
!  !wh = czero
!  !do l1=0,lmax_gamma
!  !    do m1=-l1,l1
!  !        lm1 = 1 + l1*(l1+1) + m1
!  !        do l2=0,lmax_gamma
!  !            if(l2.ne.l1) continue
!  !            do m2=-l2,l2
!  !                if(m2+m1.ne.0) continue
!  !                lm2 = 1 + l2*(l2+1) + m2
!  !                do im=1,matsiz
!  !                    wv(im) = wv(im) + h_g_lm(lm1) * a_lm(lm2,im) * (-1)**m2
!  !                    wh(im) = wh(im) + h_g_lm(lm1) * b_lm(lm2,im) * (-1)**m2
!  !                enddo
!  !            enddo
!  !        enddo
!  !    enddo
!  !enddo
!  wv = - wv *sqrt(4.0D0*pi)
!  wh = - wh *sqrt(4.0D0*pi)
!
!  deallocate(sph_harms, h_w, a_lm, b_lm, qmax_g_lm)
!
end subroutine angint_eps_sph


subroutine angint_eps_sph_oo(at, iom, bodyinv, use_harm)
! calculate the anisotropic term in the body of the inverse of dielectric matrix
! by direct calculation, or by the use of spherical harmonics

  implicit none
  type(aniso_tensor),pointer :: at
  integer,intent(in) :: iom  ! index of frequency
  complex(8),intent(inout) :: bodyinv(at%msiz,at%msiz)
  logical :: use_harm        ! flag to use expansion on spherical harmonics
  ! the direct inverse of the body of the dielectric matrix
!
!  ! local variables
  integer :: iang
  complex(8) :: cw
  integer :: im,jm             ! Counter: runs over matsiz
!  complex :: head_tmp, bodyinv_tmp
!  complex(8),allocatable :: a_lm(:,:),b_lm(:,:)
!  complex(8),allocatable :: h_w(:)  ! head of invers, times w_gamma
!  !complex(8),allocatable :: q_aob_q(:) 
!  integer :: iq0               ! Counter: runs over nq0
!  integer :: ilm,lm1,lm2,lm3   ! Counter: runs over lmgsq
!  integer :: l1,l2,l3,m1,m2,m3 ! Counters: runs over angular moment
!  integer :: i,j               ! Counters: Cartesian axis
!  integer :: scheme_head = 2
!  logical :: ldbg=.true.
!  complex(8) :: ccoefcoul_q0, ten_aob_tmp(3,3)
!
  external ylm
  complex(8),external :: zdotu,zdotc

  if(use_harm)then
  else
    ! TODO maybe need optimize
    do iang=1,n_ang_grid
      cw = cmplx(at%qmax3d3v(iang)*at%w_q0(iang),0.0D0,8)
      do im=1,at%msiz
        do jm=1,at%msiz
!          bodyinv(im,jm) = bodyinv(im,jm) + &
!              cw*wv_g(iang,im,iom)*wh_g(iang,jm,iom)/head_g(iang,iom)
        enddo
      enddo
    enddo
  endif
end subroutine angint_eps_sph_oo
!
subroutine aniso_calc_sing_q0_1_oo(at, iom, minm, singh, singw)
!
  implicit none
  type(aniso_tensor),pointer :: at
  integer,intent(in) :: iom
  complex(8),intent(in) :: minm(at%msiz)
  complex(8),intent(out) :: singh, singw
!
  ! Contribution from head
  ! Two equations below are equivalent, if the expansion of Ylm is complete
  !singh=cdot_over_lm(lmgsq,'C',qmax_g_lm,h_g_lm(:,iom))/ctwopi
  singh=cdot_over_ang('N',at%nang,cmplx(at%qmax_q0*at%w_q0,0.0D0,8),at%head_q0(:,iom)-cone)/ctwopi/cpi

  ! Contribution from wings
  singw=czero
!
end subroutine aniso_calc_sing_q0_1_oo

subroutine aniso_calc_sing_q0_1(iom, minm, singh, singw)
!
  implicit none
  integer,intent(in) :: iom
  complex(8),intent(in) :: minm(matsiz)
  complex(8),intent(out) :: singh, singw
!
  ! Contribution from head
  ! Two equations below are equivalent, if the expansion of Ylm is complete
  !singh=cdot_over_lm(lmgsq,'C',qmax_g_lm,h_g_lm(:,iom))/ctwopi
  singh=cdot_over_ang('N',n_ang_grid,cmplx(qmax_gamma*ang_weight,0.0D0,8),head_g(:,iom)-cone)/ctwopi/cpi

  ! Contribution from wings
  singw=czero
!
end subroutine aniso_calc_sing_q0_1

!
!  SUBROUTINE angint_invq2_dhead(iom, angint)
!  ! calcualte 
!  ! \frac{1}{V_{\Gamma}}\int_{V_{\Gamma}}
!  ! {\dd{\hat{\mathbf{q}}}\left\lbrace\varepsilon^{-1}_{00}(\mathbf{q}\to0,\textt{iomega})-1\right\rbrace}
!  ! $V_{\Gamma}$ is calculated by 
!  integer,intent(in) :: iom
!  complex(8),intent(out) :: angint
!  complex(8),external :: zdotu
!
!  angint = zdotu(nq0, head_g(:,iomega)-cone, 1, cmplx(qmax_q0(:)*wt_q0_sph(:),0.0D0,8), 1) &
!                / cmplx(norm_w_q0 * vol_q0, 0.0D0, 8)
!  
!  END SUBROUTINE angint_invq2_dhead


! Private functions and subroutines
function cdot_over_lm(op,nlm,funca,funcb)
  ! Calculate \sum_{lm}{op(funca_{lm}) * funcb_{lm}}

  implicit none
  character,intent(in) :: op
  integer(4),intent(in) :: nlm
  complex(8),intent(in) :: funca(nlm)
  complex(8),intent(in) :: funcb(nlm)
  complex(8) :: cdot_over_lm

  complex(8),external :: zdotc,zdotu

  if(op.eq.'C'.or.op.eq.'c')then
    cdot_over_lm = zdotc(nlm,funca,1,funcb,1)
  else
    cdot_over_lm = zdotu(nlm,funca,1,funcb,1)
  endif

end function cdot_over_lm

function cdot_over_ang(op,nang,funca,funcb)
  ! Calculate \sum_{iang}{op(funca_{iang}) * funcb_{iang}}
  
  character,intent(in) :: op
  integer(4),intent(in) :: nang
  complex(8),intent(in) :: funca(nang)
  complex(8),intent(in) :: funcb(nang)
  complex(8) :: cdot_over_ang

  complex(8),external :: zdotc,zdotu

  if(op.eq.'C'.or.op.eq.'c')then
    cdot_over_ang = zdotc(nang,funca,1,funcb,1)
  else
    cdot_over_ang = zdotu(nang,funca,1,funcb,1)
  endif

end function cdot_over_ang

end module anisotropy

