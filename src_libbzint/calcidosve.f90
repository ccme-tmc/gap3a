!BOP
!
! !ROUTINE: calcidosve
!
! !INTERFACE:
      subroutine calcidosve(nik,nbd,nsp, eband,ntet,tetc,wtet,vt,emin,&
     &                   emax,ne,dens)
!    
! !DESCRIPTION:
!
! This subroutine calculates the integrated density of states between emin and emax
!
     
! !USES:
      use order
      
      implicit none     
      
! !INPUT PARAMETERS:

      integer(4), intent(in) :: nik ! Number of irreducible k-points
      
      integer(4), intent(in) :: nbd ! Maximum number of bands
      
      integer(4), intent(in) :: nsp ! 1 / 2 for spin unpolarized/polarized
      
      real(8), intent(in) :: eband(nbd,nik,nsp) ! Band energies
      
      integer(4), intent(in) :: ntet        ! Number of tetrahedra
      
      integer(4), intent(in) :: tetc(4,*)! id. numbers of the corners
!                                             of the tetrahedra
  
      integer(4), intent(in) :: wtet(*)  ! weight of each tetrahedron
      
      real(8), intent(in)    :: vt         ! the volume of the tetrahedra

      real(8), intent(in)    :: emin        ! minimum energy
      
      real(8), intent(in)    :: emax        ! maximum energy
      
      integer(4), intent(in) :: ne          ! number of energy steps
      
! !OUTPUT PARAMETERS:      
      
      real(8), intent(out)   :: dens(2,*)    ! the energy in the first
!                                               column and the density 
!                                               of states in the second
 
!

! !LOCAL VARIABLES:

      integer(4) :: itet,ine,i,ib,isp
      real(8) :: dele
      real(8), dimension(4) :: ee
      real(8), external :: intdos1t

! !REVISION HISTORY:
!
! Created:  3th. March 2004. by RGA 
!
!EOP
!
!BOC

      dele=(emax-emin)/dble(ne-1)
      do isp=1,nsp
        do ine=1,ne
          dens(1,ine)=emin+dele*dble(ine-1)
          dens(2,ine)=0.0d0
          do itet=1,ntet
            do ib=1,nbd 
              do i=1,4
                ee(i)=eband(ib,tetc(i,itet),isp)
              enddo
              call sort(4,ee)
              dens(2,ine)=dens(2,ine)+wtet(itet)*intdos1t(ee,dens(1,ine),vt)
            enddo
          enddo
        enddo
      enddo
      return
      
      end subroutine calcidosve
      
!EOC      
          
        
      
