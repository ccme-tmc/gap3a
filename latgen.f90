!BOP
!
! !ROUTINE: latgen
!
! !INTERFACE:

      subroutine latgen

! !DESCRIPTION:
!
! LATGEN generates the Bravais matrix BR2(3,3), which transforms
! a reciprocal lattice vector of a special coordinate system (in
! units of $\frac{2\pi}{alat(i)}$) into Cartesian system
!
! JH: each column of BR2 represent a reciprocal lattice vector 
!     gbas differs from BR2 only by a factor of 2 \pi  
!     rbas, on the other hand, is represented by the row-wise mode 
!     i.e. each row of rbas represent a direct space lattice vector 
!
! !USES:

      use constants, only : pi
      use struk,     only : nat,alat, alpha, pia, vi, mult,br2,ortho,    &
     &                      lattic,gbas,rbas,imat,izmat,nsym

! !LOCAL VARIABLES:

      implicit none

      character(10):: sname="latgen"
      logical:: ldbg = .true.
      integer(4) :: i,j
      integer(4) :: ix
      integer(4) :: iy
      integer(4) :: iz

      real(8) :: cosab
      real(8) :: cosac
      real(8) :: rvfac
      real(8) :: sinab
      real(8) :: wurzel
      real(8) :: cosbc
      real(8) :: sinbc
      real(8) :: mtmp(3,3)


! !EXTERNAL ROUTINES: 


      external outerr
      external rotdef
      external rbass


! !INTRINSIC ROUTINES: 


      intrinsic abs
      intrinsic atan
      intrinsic cos
      intrinsic sin
      intrinsic sqrt

! !REVISION HISTORY:
!
!   Last modified Jan. 27th. 2004.
!
!EOP
! ------------------------------------------------
!BOC
!
      call linmsg(6,'-',"Lattice generation")

      ortho = .false.

      pia(1) = 2.0d+0*pi/alat(1)
      pia(2) = 2.0d+0*pi/alat(2)
      pia(3) = 2.0d+0*pi/alat(3)

      select case (lattic(1:1))
      case ('H')       ! hexagonal

        br2(1,1) = 2.0d+0/sqrt(3.0d+0)
        br2(2,1) = 0.0d+0
        br2(3,1) = 0.0d+0

        br2(1,2) = 1.0d+0/sqrt(3.0d+0)
        br2(2,2) = 1.0d+0
        br2(3,2) = 0.0d+0

        br2(1,3) = 0.0d+0
        br2(2,3) = 0.0d+0
        br2(3,3) = 1.0d+0

        rvfac = 2.0d+0/sqrt(3.0d+0)
        ortho = .false.
        call renorm_br2

      case ('S','P')       ! primitive or simple
        ortho=.true.
        do i = 1,3
          if(abs(alpha(i)-pi/2.0d0).gt.0.0001)then
            write(6,1000) i
            ortho = .false.
          endif
        enddo
!
!       includes triclinic, monoclinic, simple orthorhombic, tetragonal
!       and cubic
!
        sinbc=sin(alpha(1))
        cosab=cos(alpha(3))
        cosac=cos(alpha(2))
        cosbc=cos(alpha(1))
        wurzel=sqrt(sinbc**2-cosac**2-cosab**2+2*cosbc*cosac*cosab)

        br2(1,1)= sinbc/wurzel*pia(1)
        br2(2,1)= 0.0d0
        br2(3,1)= 0.0d0

        br2(1,2)= (-cosab+cosbc*cosac)/(sinbc*wurzel)*pia(2)
        br2(2,2)= pia(2)/sinbc
        br2(2,3)= -pia(3)*cosbc/sinbc

        br2(1,3)= (cosbc*cosab-cosac)/(sinbc*wurzel)*pia(3)
        br2(3,2)= 0.0d0
        br2(3,3)= pia(3)

        rvfac= 1.d0/wurzel

      case ('F')       ! face centered

        br2(1,1) = -1.0d+0
        br2(2,1) =  1.0d+0
        br2(3,1) =  1.0d+0

        br2(1,2) =  1.0d+0
        br2(2,2) = -1.0d+0
        br2(3,2) =  1.0d+0

        br2(1,3) =  1.0d+0
        br2(2,3) =  1.0d+0
        br2(3,3) = -1.0d+0

        rvfac = 4.0d+0
        ortho = .true.
        call renorm_br2

      case ('B')       ! body centered

        br2(1,1) = 0.0d+0
        br2(2,1) = 1.0d+0
        br2(3,1) = 1.0d+0

        br2(1,2) = 1.0d+0
        br2(2,2) = 0.0d+0
        br2(3,2) = 1.0d+0

        br2(1,3) = 1.0d+0
        br2(2,3) = 1.0d+0
        br2(3,3) = 0.0d+0

        rvfac = 2.0d+0
        ortho = .true.
        call renorm_br2

      case ('C')       ! base centered

        select case (lattic(2:3))
        case ('XZ')  !xz
          ix=1
          iy=2
          iz=3
        case ('YZ') !yz
          ix=2
          iy=1
          iz=3
        case ('XY') !xy
          ix=1
          iy=3
          iz=2
        case default ! not found
          call outerr(sname,"c_case not found!")
        end select
!
!       orthorombic case
!
        if(abs(alpha(iz)-pi/2.0d0).lt.0.0001) then

          br2(ix,ix) =  1.0d+0
          br2(ix,iy) =  0.0d+0
          br2(ix,iz) =  1.0d+0
          br2(iy,ix) =  0.0d+0
          br2(iy,iy) =  1.0d+0
          br2(iy,iz) =  0.0d+0
          br2(iz,ix) = -1.0d+0
          br2(iz,iy) =  0.0d+0
          br2(iz,iz) =  1.0d+0

          rvfac = 2.0d+0
          ortho = .true.
          call renorm_br2

        else
!
!         monoclinic case
!
          write(6,1000)iz
          sinab=sin(alpha(iz))
          cosab=cos(alpha(iz))
          br2(ix,ix)= pia(ix)*1.0d+0/sinab
          br2(ix,iy)= -pia(iy)*cosab/sinab
          br2(ix,iz)= pia(ix)*1.0d+0/sinab
          br2(iy,ix)= 0.0
          br2(iy,iy)= pia(iy)*1.0d+0
          br2(iy,iz)= 0.0
          br2(iz,ix)=-pia(iz)*1.0d+0
          br2(iz,iy)= 0.0
          br2(iz,iz)= pia(iz)*1.0d+0

          rvfac=2.0/sinab
          ortho=.false.

        endif

      case ('R')       ! rhombohedral

        br2(1,1) =  1.0d+0/sqrt(3.0d+0)
        br2(2,1) = -1.0d+0
        br2(3,1) =  1.0d+0

        br2(1,2) =  1.0d+0/sqrt(3.0d+0)
        br2(2,2) =  1.0d+0
        br2(3,2) =  1.0d+0

        br2(1,3) = -2.0d+0/sqrt(3.0d+0)
        br2(2,3) =  0.0d+0
        br2(3,3) =  1.0d+0

        rvfac = 6.0d+0/sqrt(3.0d+0)
        ortho = .false.
        call renorm_br2

      case default       ! error: wrong lattice, exit 

        call outerr(sname,'wrong lattice.')

      end select 
!
!     define inverse of cellvolume
!
      vi = rvfac/ (alat(1)*alat(2)*alat(3))
!
!     Calculate the basis vectors of the real lattice
!
      do i = 1,3
        do j = 1,3
            gbas(i,j) = br2(i,j)/2.0d+0/pi
        enddo
      enddo      
      call rbass(gbas,rbas)
!
!     define rotation matrices if required
!
      call rotdef

!
! Set symmetry operation matrix in internal coordinates
!
      if(ortho) then 
        call sym2int(rbas,nsym,imat,izmat)
      else    
        izmat = imat
      endif
      
      write(6,101) lattic,ortho
      write(6,*)'Reciprocal lattice basis vectors:'
      do i=1,3
        write(6,103) br2(1:3,i)
      enddo  
      write(6,*)'Real lattice basis vectors:'
      do i=1,3
        write(6,103) rbas(i,1:3)
      enddo  

      if(ldbg) then 
        mtmp=matmul(br2,rbas)
        write(6,*)'br2.rbas='
        do i=1,3
          write(6,103) mtmp(i,1:3)
        enddo

        write(6,105) vi
      endif 
      
      write(6,*)'------------------------------------------------------'
!
!     formats
!
  101 format('Lattice type: ',a4,3x,'ortho =  ',l1,/)
  103 format(10x,3f14.10)
  105 format('inverse of the WZ cell voulme: ',1pg16.8,/)   
 1000 format('alpha(',i2,') not equal 90')

      contains


      subroutine renorm_br2
! !DESCRIPTION: Internal subroutine : Normalize bravais matrix
      implicit none
      integer(2) :: i
      integer(2) :: j
      
      do i = 1, 3
        do j = 1, 3
          br2(j,i)=br2(j,i)*pia(j)   
        enddo
      enddo
      end subroutine renorm_br2
      end subroutine latgen
!EOC
