MODULE ANISOTROPY

    ! this module defines variables and subroutines utilized for dealing
    ! with the anisotropy of dielectric matrix for $q\to0$

    use mixbasis, only: matsiz
    use task, only: fid_outgw

    integer :: iop_aniso = -1  ! control whether to consider the anisotropy of dielectric function around Gamma
                               ! -1 -- use q0_eps, no anisotropy
    complex(8), pointer :: vec_u_ani(:,:,:)   ! vector u
    complex(8), pointer :: vec_a_ani(:,:,:)   ! vector a
    complex(8), pointer :: ten_p_ani(:,:,:)   ! tensor P
    complex(8), pointer :: ten_a_ani(:,:,:)   ! tensor A
    ! Following quantities should be determined by a Lebedev-Laikov grid
    integer :: nq0 = 0                        ! number of q0 for angular integration
    real(8), pointer :: q0_sph(:,:)           ! similar to q0_eps
    real(8), pointer :: w_q0_sph(:)           ! weight of q0

    complex(8), pointer :: head_q0(:,:)       ! the head (nq0, nomega) for q0 

    CONTAINS

        SUBROUTINE init_aniso(iomfirst, iomlast)

            use lebedev_laikov

            integer,intent(in) :: iomfirst, iomlast
            integer :: ierr

            call set_lebedev_laikov_grid(nq0)
            nq0 = nleb

            if(iop_aniso.ne.-1)then
                write(fid_outgw,*) "Anisotropy switched on"
                allocate(vec_u_ani(3,matsiz,iomfirst:iomlast), &
                &        vec_a_ani(3,matsiz,iomfirst:iomlast), &
                &        ten_p_ani(3,3,iomfirst:iomlast),      &
                &        ten_a_ani(3,3,iomfirst:iomlast),      &
                &        q0_sph(3,nq0),                        &
                &        w_q0_sph(nq0),                        &
                &        head_q0(nq0, iomfirst:iomlast),       &
                &        stat=ierr)
                if(ierr.ne.0) then
                    write(fid_outgw,*) " - init_aniso: Fail to allocate aniso"
                    stop
                else
                    write(fid_outgw,*) " - init_aniso: success"
                endif
            endif

            q0_sph(1,:) = xleb(:)
            q0_sph(2,:) = yleb(:)
            q0_sph(3,:) = zleb(:)
            w_q0_sph(:) = wleb(:)
            call unset_lebedev_laikov_grid

        END SUBROUTINE init_aniso

        SUBROUTINE end_aniso

            deallocate(vec_u_ani, vec_a_ani, ten_p_ani, q0_sph, w_q0_sph, head_q0)
            write(fid_outgw,*) " - end_aniso: success"

        END SUBROUTINE end_aniso

END MODULE ANISOTROPY
      