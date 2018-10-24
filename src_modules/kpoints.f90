!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 
      module kpoints
      
! !PUBLIC VARIABLES:

!      integer :: fflg    ! flag for the q dependent integration.
      integer :: nirkp   ! number of irreducible k-points
      integer :: nkp     ! number of k-points
      integer :: nkp_ene ! number of k-points in case.energy
      integer :: nqp     ! number of q-points
      integer :: idvq    ! Minimum common divisor of qlist
      integer :: idvk    ! Minimum common divisor of klist
      integer :: idvkir  ! Minimum common divisor of kirlist
      
      integer :: nkdivs(3)       ! Number of k-points in each direction in BZ
      integer :: k0shift(4)      ! Shift vector of the k-submesh
      integer, allocatable :: idikp(:)     ! kpoint index of the irred. kpoint
      integer, allocatable :: iksym(:)     ! Index of the symmetry operation relating the k-point of the irred kpoint
      integer, allocatable :: kpirind(:)   ! Index of the irreducible k-point of the given kpoint
      integer, allocatable :: kirlist(:,:) ! List of irreducible k-points in submesh coordinates
      integer, allocatable :: klist(:,:)   ! List of k-points in  submesh coordinates
      integer, allocatable :: wkir(:)      ! Geometrical weight of each irreducible k-point 
      integer, allocatable :: kqid(:,:)    ! kqid(ik,iq) returns the index of k-q 
      integer, allocatable :: qlist(:,:)   ! List of q-points in  submesh coordinates
      integer, allocatable :: weightk(:)   ! Geometrical weight of each  k-point
      integer, allocatable :: weightq(:)   ! Geometrical weight of each  q-point
      integer, allocatable :: g0(:,:)
 
      real(8), allocatable:: kirvecs(:,:)   ! list of irreducible k-points 
      real(8), allocatable:: kvecs(:,:)     ! k-vectors in the full BZ in the units of 
                                           ! reciprocal cell vectors 

      integer :: nirtet  ! Number of irreducible tetrahedra
      integer :: ntet    ! Number of tetrahedra
      integer, allocatable :: tnodes(:,:)  ! index of the k-points corresponding to the nodes of each tetrahedra for convolution
      integer, allocatable :: tndi(:,:)    ! index of the k-points corresponding to the nodes of each tetrahedra for integration
      integer, allocatable :: link(:,:)    ! index of the tetrahedra linked to by the corresponding q vector
      integer, allocatable :: wtet(:)      ! geometrical weight of each tetrahedron  for integration
      integer, allocatable :: wirtet(:)    ! geometrical weight of each irred. tetrahedron  for integration
      real(8) :: tvol                      ! volume of the tetrahedra relative to the BZ volume

      real(8) :: nvel                      ! Number of valence electrons
      real(8) :: nvelgw                    ! Number of valence electrons included in GW bands 

      contains 
      
      function get_kindex(iop,kv) result(ik_out)
      integer:: iop  !! 0/1 -- irreducible/reducible k
      real(8):: kv(3) !! integer k-coordinates
      integer:: ik_out
      real(8):: zero_diff=1.e-5

      integer:: ik

      ik_out = 0
      if(iop.eq.0) then 
        do ik=1,nirkp
          if(maxval(abs( dble(kirlist(:,ik))/idvkir-kv(:)))<zero_diff) then
            ik_out = ik
            exit 
          endif 
        enddo 
      else
        do ik=1,nkp
          if(maxval(abs(dble(klist(:,ik))/idvk-kv(:)))<zero_diff) then
            ik_out = ik
            exit
          endif
        enddo
      endif 
      end function 

      subroutine get_kvec(isym,iik,ik,irk,kvec) 
      integer,intent(in):: isym,iik
      integer,intent(out):: ik,irk
      real(8),intent(out),optional:: kvec(3) 
      if(isym.eq.0) then
        ik=iik
        irk=kpirind(ik)
        if(present(kvec)) kvec=kvecs(:,ik) 
      else
        irk=iik
        ik=idikp(irk)
        if(present(kvec)) kvec=kirvecs(:,irk) 
      endif
      end subroutine

      end module kpoints
