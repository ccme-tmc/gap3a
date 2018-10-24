!BOP
!
! !MODULE: mommat
      module mommat

! !DESCRIPTION:
!
! Declares the shared variables used for the calculation of the matrix
!elements of the linear moment
!
! For a description of the naming conventions see eqs. \ref{momradintc1}, 
!\ref{momradintc2}, \ref{momradintv1} and \ref{momradintv2}.
!
!
! !PUBLIC VARIABLES:

!
! Arrays for the radial integrals needed to calculate the
!momentum matrix elements        
!
        implicit none 
        integer :: iop_mommat=0 
        logical::lhermit =.true.        
        logical::lrenorm =.false.        

        real(8), allocatable :: iul1ul(  :,:,:) 
        real(8), allocatable :: iul1udl( :,:,:) 
        real(8), allocatable :: iudl1ul( :,:,:) 
        real(8), allocatable :: iudl1udl(:,:,:) 

        real(8), allocatable :: iulul1(  :,:,:) 
        real(8), allocatable :: iuludl1( :,:,:) 
        real(8), allocatable :: iudlul1( :,:,:) 
        real(8), allocatable :: iudludl1(:,:,:) 

        !! LO related 
        real(8), allocatable :: iul1ulol( :,:,:,:) 
        real(8), allocatable :: iudl1ulol(:,:,:,:) 
        real(8), allocatable :: iulol1ul( :,:,:,:) 
        real(8), allocatable :: iulol1udl(:,:,:,:) 
        real(8), allocatable :: iululol1( :,:,:,:) 
        real(8), allocatable :: iulolul1( :,:,:,:) 
        real(8), allocatable :: iuloludl1(:,:,:,:) 
        real(8), allocatable :: iudlulol1(:,:,:,:) 
        real(8), allocatable :: iulol1ulol(:,:,:,:,:) 
        real(8), allocatable :: iulolulol1(:,:,:,:,:) 


        !! core states related 
        real(8), allocatable :: iucl1ul(:,:,:) 
        real(8), allocatable :: iucl1udl(:,:,:) 
        real(8), allocatable :: iucl1ucl(:,:,:,:)
        real(8), allocatable :: iucl1ulol(:,:,:,:)

        real(8), allocatable :: iuclul1(:,:,:) 
        real(8), allocatable :: iucludl1(:,:,:) 
        real(8), allocatable :: iuclucl1(:,:,:,:)
        real(8), allocatable :: iuclulol1(:,:,:,:)

        real(8), allocatable :: iudl1ucl(:,:,:) 
        real(8), allocatable :: iul1ucl(:,:,:) 
        real(8), allocatable :: iulol1ucl(:,:,:,:) 
        real(8), allocatable :: iulucl1(:,:,:) 
        real(8), allocatable :: iudlucl1(:,:,:) 
        real(8), allocatable :: iulolucl1(:,:,:,:)

        complex(8), allocatable :: mmatvv(:,:,:,:,:)
        complex(8), allocatable :: mmatcc(:,:,:,:,:)
        complex(8), allocatable :: mmatcv(:,:,:,:,:)

!EOP        
      contains

!BOP
        subroutine init_mommat(in0,in1,im0,im1,nk,nsp)
        use bands,     only: nomaxs,numins,nbmaxpol
        use constants, only: czero
        use core,      only: ncg,nclmmax
        use lapwlo,    only: nlomax,lmax 
        use struk,     only: nat,ndf
        implicit none
        integer,intent(in):: nk,nsp
        integer:: im0,im1,in0,in1 

        integer(4) ::ierr 

        allocate(mmatvv(3,in0:in1,im0:im1,nk,nsp),                  &
     &           mmatcv(3,ncg,im0:im1,nk,nsp),                      &
     &           mmatcc(3,nclmmax,nclmmax,nat,nsp),                 &
     &           stat=ierr )
        if(ierr.ne.0) then 
          write(6,*) "init_mommat: Fail to allocate memory for mmat"
          stop
        endif
        mmatvv = czero
        mmatcv = czero
        mmatcc = czero
        end subroutine init_mommat

        subroutine end_mommat
          deallocate(mmatvv,mmatcv,mmatcc)
        end subroutine end_mommat
      
        subroutine init_rad_mom
        use bands,   only: nspin
        use lapwlo,  only: lmax,lomax,nLOmax
        use core,    only: ncoremax
        use struk,   only: nat
        implicit none

        allocate(iul1ul(0:lmax-1,nat,nspin)) 
        allocate(iul1udl(0:lmax-1,nat,nspin)) 
        allocate(iul1ulol(nLOmax,0:lomax,nat,nspin)) 

        allocate(iudl1ul(0:lmax-1,nat,nspin)) 
        allocate(iudl1udl(0:lmax-1,nat,nspin))
        allocate(iudl1ulol(nLOmax,0:lomax,nat,nspin)) 

        allocate(iulol1ul(nLOmax,0:lomax-1,nat,nspin)) 
        allocate(iulol1udl(nLOmax,0:lomax-1,nat,nspin)) 
        allocate(iulol1ulol(nLOmax,nLOmax,0:lomax-1,nat,nspin)) 
 
        allocate(iulul1(0:lmax-1,nat,nspin)) 
        allocate(iuludl1(0:lmax-1,nat,nspin)) 
        allocate(iululol1(nLOmax,0:lomax-1,nat,nspin)) 

        allocate(iudlul1(0:lmax-1,nat,nspin)) 
        allocate(iudludl1(0:lmax-1,nat,nspin)) 
        allocate(iudlulol1(nLOmax,0:lomax-1,nat,nspin)) 

        allocate(iulolul1(nLOmax,0:lomax,nat,nspin)) 
        allocate(iuloludl1(nLOmax,0:lomax,nat,nspin)) 
        allocate(iulolulol1(nLOmax,nLOmax,0:lomax-1,nat,nspin)) 
        
        allocate(iucl1ul(  ncoremax,nat,nspin)) 
        allocate(iucl1udl( ncoremax,nat,nspin)) 
        allocate(iucl1ulol(nLOmax,ncoremax,nat,nspin))
        allocate(iucl1ucl( ncoremax,ncoremax,nat,nspin))

        allocate(iuclul1(  ncoremax,nat,nspin)) 
        allocate(iucludl1( ncoremax,nat,nspin)) 
        allocate(iuclulol1(nLOmax,ncoremax,nat,nspin))
        allocate(iuclucl1( ncoremax,ncoremax,nat,nspin))

        allocate(iul1ucl(  ncoremax,nat,nspin)) 
        allocate(iudl1ucl( ncoremax,nat,nspin)) 
        allocate(iulol1ucl(nLOmax,ncoremax,nat,nspin)) 

        allocate(iulucl1(  ncoremax,nat,nspin)) 
        allocate(iudlucl1( ncoremax,nat,nspin)) 
        allocate(iulolucl1(nLOmax,ncoremax,nat,nspin))
      
        iul1ul = 0.0d0 
        iul1udl = 0.0d0 
        iudl1ul = 0.0d0 
        iudl1udl = 0.0d0 
        iulul1 = 0.0d0 
        iuludl1 = 0.0d0 
        iudlul1 = 0.0d0 
        iudludl1 = 0.0d0 
        iul1ulol = 0.0d0 
        iulol1ul = 0.0d0 
        iulol1udl = 0.0d0 
        iudl1ulol = 0.0d0 
        iulol1ulol = 0.0d0 
        iululol1 = 0.0d0 
        iulolul1 = 0.0d0 
        iuloludl1 = 0.0d0 
        iudlulol1 = 0.0d0 
        iulolulol1 = 0.0d0 
        iucl1ul = 0.0d0 
        iul1ucl = 0.0d0 
        iucl1udl = 0.0d0 
        iudl1ucl = 0.0d0 
        iulol1ucl = 0.0d0 
        iuclul1 = 0.0d0 
        iulucl1 = 0.0d0 
        iucludl1 = 0.0d0 
        iudlucl1 = 0.0d0 
        iulolucl1 = 0.0d0
        iuclulol1 = 0.0d0
        iucl1ulol = 0.0d0
        iuclucl1 = 0.0d0
        iucl1ucl = 0.0d0

        end subroutine init_rad_mom

        subroutine end_rad_mom
        deallocate(iul1ul) 
        deallocate(iul1udl) 
        deallocate(iudl1ul) 
        deallocate(iudl1udl) 
        deallocate(iulul1) 
        deallocate(iuludl1) 
        deallocate(iudlul1) 
        deallocate(iudludl1) 
        deallocate(iul1ulol) 
        deallocate(iulol1ul) 
        deallocate(iulol1udl) 
        deallocate(iudl1ulol) 
        deallocate(iulol1ulol) 
        deallocate(iululol1) 
        deallocate(iulolul1) 
        deallocate(iuloludl1) 
        deallocate(iudlulol1) 
        deallocate(iulolulol1) 
      
        end subroutine end_rad_mom

      end module mommat

!EOC          
