!-----------------------------------------------------------------------
!BOP
!
! !MODULE: eigenvec
      module eigenvec
      use task,     only: iop_scratch
      implicit none

! !PUBLIC VARIABLES:      
      logical:: ltranvec=.false.          ! control whether the vectors are transformed after reading from disk
      integer:: ik_cur_evec, jk_cur_evec  ! index for current k and k-q
      complex(8), allocatable :: zzk(:,:) ! Eigenvector coefficients at k
      complex(8), allocatable :: zzq(:,:) ! Eigenvector coefficients at k-q
      integer, allocatable :: kzz(:,:,:)            ! reciprocal lattice vectors  
      complex(8), allocatable :: zzkall(:,:,:,:)   ! zzk incore 

      complex(8), allocatable :: alm(:,:,:) !coefficient alm(k+G) for atom i (including equivalent atoms)
      complex(8), allocatable :: blm(:,:,:) !coefficient alm(k+G) for atom i (including equivalent atoms)
      complex(8), allocatable :: clm(:,:,:,:) !coefficient alm(k+G) for atom i (including equivalent atoms)
      complex(8), allocatable :: alfa(:,:,:)
      complex(8), allocatable :: beta(:,:,:)
      complex(8), allocatable :: gama(:,:,:,:)
      complex(8), allocatable :: alfp(:,:,:)
      complex(8), allocatable :: betp(:,:,:)
      complex(8), allocatable :: gamp(:,:,:,:)

      real(8), allocatable :: almr(:,:,:) ! coefficient alm(kn) for atom i (including equivalent atoms) real part
      real(8), allocatable :: almi(:,:,:) ! coefficient alm(kn) for atom i (including equivalent atoms) imaginary part
      real(8), allocatable :: blmr(:,:,:) ! coefficient blm(kn) for atom i (including equivalent atoms) real part
      real(8), allocatable :: blmi(:,:,:) ! coefficient blm(kn) for atom i (including equivalent atoms) imaginary part
      real(8), allocatable :: clmr(:,:,:) ! coefficient clm(kn) for atom i (including equivalent atoms) real part
      real(8), allocatable :: clmi(:,:,:) ! coefficient clm(kn) for atom i (including equivalent atoms)  imaginary part

      logical:: lcmplx       !! indicate whether the WIEN2K vectors are complex or not  
      logical:: lsymvector   !! indicate whether the WIEN2K vectors are symmetric or not 

      character(len=120)::vfname !! Direct access vector file name 
      integer::vfunit,vfrecl   !! Direct access vector file unit and record length

      contains
!EOP      
!BOP
!
! !IROUTINE: init_eigenvec
!
! !INTERFACE:      
        subroutine init_eigenvec

! !USES:
        use bands,     only: nbmax,nspin
        use kpoints,   only: nirkp,nkp
        use recipvec,  only: maxngk
        use lapwlo,    only: nt,nLOmax,lomax
        use constants, only: czero
        use struk,    only: ndf

        implicit none     
        integer :: ierr,nk    
         
!EOP
!BOC
        call linmsg(6,'-',"Init_eigenvec") 
        write(6,*) "Memory used for  eigenvecs: "
        write(6,*) "  maxngk =",maxngk 
        write(6,*) "  nbmax  =",nbmax
        write(6,*) "  Memory(KB)",nint(1.0e-3*maxngk*nbmax*64*2) 
        ik_cur_evec = 0 
        jk_cur_evec = 0 

        allocate(zzk(1:maxngk,1:nbmax),       &
     &           zzq(1:maxngk,1:nbmax),       &
     &           kzz(1:3,maxngk,nirkp),       &
     &           alfa(nbmax,nt*nt,ndf),       &
     &           beta(nbmax,nt*nt,ndf),       &
     &           alfp(nbmax,nt*nt,ndf),       &
     &           betp(nbmax,nt*nt,ndf),       &
                 stat=ierr) 
        if(ierr.ne.0) then 
          write(6,*) "ERROR init_eigenvec: fail to allocate memory"
          stop
        endif 
        if(nlomax.gt.0) then 
          allocate(gama(nbmax,nLOmax,(lomax+1)**2,ndf),       &
     &             gamp(nbmax,nLOmax,(lomax+1)**2,ndf),       &
                 stat=ierr)
        endif 
 
        kzz=0
        zzk=czero
        zzq=czero

        nk=nkp
        if(lsymvector) nk=nirkp
        if(iop_scratch.eq.0) then
          allocate(zzkall(maxngk,nbmax,nk,nspin),stat=ierr)
          if(ierr.ne.0) then
            write(6,*) "ERROR init_eigenvec: fail to allocate zzkal"
            write(6,'(a,f8.3)') " Required memory (MB)",1.d0*maxngk*nbmax*nk &
     &                *16/2.d0**20
            stop
          endif
          write(6,'(a,f8.3)') "Memory used for zzk (MB)",1.d0*maxngk*nbmax*nk &
     &                *16/2.d0**20

        endif
        end subroutine init_eigenvec

        subroutine end_eigenvec
          deallocate(zzk,zzq,kzz,alfa,beta,alfp,betp)
          if(allocated(gama)) deallocate(gama,gamp)
          if(iop_scratch.eq.0) deallocate(zzkall)
          ik_cur_evec = 0 
          jk_cur_evec = 0 
        end subroutine end_eigenvec
      
      end module eigenvec
!EOP      
