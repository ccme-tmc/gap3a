!BOP
!
! !MODULE: polyhedron

      module polyhedron
      
! !PUBLIC TYPES:

 
      integer(4) :: nnod ! number of nodes defining the polyhedron

      integer(1), dimension(20) :: ntype ! idem as pl, but for internal
!                                          use

      integer(1), allocatable :: pl(:) ! it is an integer number such
!                                        that its binary expression gives
!                                        one for the bits that correspond
!                                        to the intersecting planes

      integer(1), allocatable :: index(:,:) ! it indicates in which order
!                                             the nodes have to be sent to
!                                             genericprism

      real(8) :: ef        ! the Fermi energy
      
      real(8), dimension(4) :: e ! band energies at k
      
      real(8), dimension(4) :: f ! band energies at k-q      
      
      real(8), allocatable  :: nodes(:,:) ! the coordinates of the corners
!                                           of the polyhedron

      real(8), dimension(3,20) :: intnodes ! the coordinates of the 
!                                             intersections of the planes

      end module polyhedron
      
!EOP      
