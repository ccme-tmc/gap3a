!BOP
!
! !ROUTINE: orthog_corewf
!
! !INTERFACE:
      subroutine orthog_corewf(iat,isp)

! !USES:

      use constants, only: cfein1, cfein2
      use core,      only: ncore,ocmax,symbl,lcore,ucore,uscore
      use lapwlo,    only: loor,lomax
      use radwf,     only: u,us,udot,usdot,ulo,uslo
      use struk,     only: nrpt
!
! !INPUT PARAMETERS:
      
      implicit none 
      
      integer,intent(in) :: iat   ! index of inequivalent atom 
      integer,intent(in) :: isp   ! index for spin 
!
! !DESCRIPTION:
!
!
! Orthogonalizes the core wave-functions to the valence ones.
!
! !LOCAL VARIABLES:

      integer(4) :: icore ! Counter: run over core states 
      integer(4) :: lc    ! Angular momentum quantum number of the core wf
      integer(4) :: npt   ! number of points in the radial mesh

      real(8) :: normcore ! norm of the core wave function
            
      real(8), allocatable :: uc(:),usc(:)      ! local arrays for the core wf
      real(8), allocatable :: uval(:),usval(:)  ! local arrays for the val wf.
      
!
!EOP
!BOC      

      npt=nrpt(iat)
      allocate(uc(npt),usc(npt),uval(npt),usval(npt))

      do icore = 1, ncore(iat) 
        lc=lcore(icore,iat)

        uc(1:npt) = ucore(1:npt,icore,iat,isp)
        usc(1:npt) = uscore(1:npt,icore,iat,isp)
!
! The norm of the core wave function ucore 
!
        call rint13(cfein1,cfein2,uc,usc,uc,usc,normcore,iat)
!
! Orthogonalize to u
!
        uval(1:npt) = u(1:npt,lc,iat,isp)
        usval(1:npt) = us(1:npt,lc,iat,isp)
        call orthogonalize
!
! Orthogonalize to udot
!
        uval(1:npt) = udot(1:npt,lc,iat,isp)
        usval(1:npt) = usdot(1:npt,lc,iat,isp)
        call orthogonalize
!
! Orthogonalize to local orbital
!        
        if(lc.le.lomax)then
          if(loor(lc,iat))then
          
            uval(1:npt) = ulo(1:npt,1,lc,iat,isp)
            usval(1:npt) = uslo(1:npt,1,lc,iat,isp)
            call orthogonalize

          endif ! loor
        endif ! lc .le. lomax    
!
! Restore core wave function
!
        ucore(1:npt,icore,iat,isp) = uc(1:npt)
        uscore(1:npt,icore,iat,isp) = usc(1:npt)
                 
      enddo !icore

      deallocate(uc,usc,uval,usval)  

      contains
      
      subroutine orthogonalize
      
      implicit none

      integer(4) :: irp     ! index of the radial mesh

      real(8) :: overlap_factor 
      real(8) :: renorm_factor
      real(8) :: overlap
      real(8) :: normafter
      real(8) :: normval
!
! The norm of the valence radial function 
!
        call rint13(cfein1,cfein2,uval,usval,uval,usval,normval,iat)
!
! Overlap between valence and core
!        
        call rint13(cfein1,cfein2,uval,usval,uc,usc,overlap,iat)
!
! Orthogonalization (Gramm-Schmith like)
!       
        overlap_factor = overlap/normval 
        do irp = 1, npt
          uc(irp) = uc(irp) - overlap_factor * uval(irp)
          usc(irp) = usc(irp) - overlap_factor * usval(irp)
        enddo ! irp
!
! Norm of the orthogonalized core function
!
        call rint13(cfein1,cfein2,uc,usc,uc,usc,normafter,iat)
!
! Renormalization
!        
        renorm_factor = sqrt(normcore/normafter)
        do irp = 1, npt
          uc(irp) = renorm_factor * uc(irp)
          usc(irp) = renorm_factor * usc(irp)
        enddo  
      end subroutine orthogonalize  
      
      end subroutine orthog_corewf
!EOC  
