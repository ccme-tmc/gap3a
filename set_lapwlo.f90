!BOP
!
! !ROUTINE: set_lapwlo
!
! !INTERFACE:
      subroutine set_lapwlo

! !DESCRIPTION:
!
! This subroutine set up informaton related to the LAPW basis set  
!
! !USES:

      use bands,      only: nspin 
      use lapwlo,     only: lmax,elapw,elo,lomax,nlmax,nlomax,nloat,lapw,  &
     &                      umt,abcelo,umtlo,nlo_at,elapw,nlo,l_newlo
      use radwf,      only: init_radwf
      use struk,      only: nat
      use task,       only: casename,iop_dftpkg,spflag
!
! !LOCAL VARIABLES:
!
      implicit none
     
      integer :: fid 
      integer :: ierr 
      integer :: iat,ieq,idf  ! Index for inequivalent/equvalent/all atoms
      integer :: irm,imix,im  ! Index mixed basis functions.
      integer :: isp
      integer :: ilo,ilo0     ! (Counter) runs over kind of local  orbitals (nlo)
      integer :: jri          ! Number of radial mesh points (for each atom  the corresponding nrpt(iat) is stored here).
      integer :: l,m          ! (Counter) runs over angular momentum l
      integer :: lms   
      integer :: maxnt    ! Maximum l for gaunt coefficients
      logical :: lprt =.false.
      character(20):: sname="set_lapwlo"
      
!

!
! !EXTERNAL ROUTINES: 
!
!
! !REVISION HISTORY:
! Created  Dec. 10, 2009 by JH 
!
!EOP
!BOC
!
!

!
!     READ (L)APW+(lo) information including linearization energy
!
      call w2k_readin1

      allocate(elapw(0:lmax,nat,nspin),                 &
     &         abcelo(1:4,0:nloat-1,0:lomax,nat,nspin),  &
     &         elo(0:lomax,0:nloat-1,nat,nspin),           &
     &         umtlo(1:4,1:nlomax,0:lomax,nat,nspin)    )
      do isp=1,nspin
        fid=999
        open(unit=fid,file=trim(casename)//".energy"//spflag(isp),&
     &         action='read',iostat=ierr)
        call errmsg(ierr.ne.0,sname,"Fail to open case.energy")
        do iat = 1,nat
          if(l_newlo) then 
            read(fid,'(1000f12.5)') elapw(0:lmax,iat,isp)
            read(fid,'(1000f12.5)') ((elo(l,ilo,iat,isp),l=0,lomax),ilo=0,nloat-1) 
          else 
            read(fid,'(100f9.5)') elapw(0:lmax,iat,isp)
            read(fid,'(100f9.5)') ((elo(l,ilo,iat,isp),l=0,lomax),ilo=0,nloat-1) 
          endif 
        enddo
        close(fid)
      enddo

      do isp=1,nspin
        ! set linearization energy
        do iat=1,nat
          do l=0,lmax
            if(lapw(l,iat)) then
              umt(1,l,iat,isp)=elapw(l,iat,isp)
            else
              umt(1,l,iat,isp)=elapw(l,iat,isp)-200.d0
            endif
            write(6,201) iat,l,umt(1,l,iat,isp)
            if(l.le.lomax) then 
              abcelo(4,:,l,iat,isp) = elo(l,:,iat,isp)
            endif 
          enddo
        enddo
      enddo
 201  format("E_l     for iat=",i4," l=",i4,f10.5)
 202  format("E_l(LO) for iat=",i4," l=",i4,10f10.5)
 

!
!     Calculate radial wave functions for LAPW basis 
!
      
      call init_radwf

      if(iop_dftpkg.eq.0) then 
        call w2k_readvsp
        do isp=1,nspin
          do iat = 1, nat
            call calcrwf(iat,isp)       ! val/cond states:  calculate  the radial functions u, udot and ulo
          enddo
        enddo
      endif 

!
!     Calculate the cofficients A,B,C of the local orbitals:
!
      do isp=1,nspin
        do iat = 1, nat
          do l = 0, lomax

            if(lapw(l,iat)) then  !! LAPW+LO
              ilo0 = 1 
            else
              ilo0 = 0      !! APW+lo+LO
            endif 

            do ilo = ilo0, nlo_at(l,iat)  
              call calcloabc(iat,l,ilo,lapw(l,iat),isp)
            enddo
          enddo
        enddo ! iat
      enddo ! isp


  101 format(' Max. nr. of MT-sphere wavefunctions per atom ',i6,/,     &
     &       ' Total  nr. of MT-sphere wavefunctions        ',i6)

      return
    
      contains 

      subroutine sub_get_nloat
      integer    i,l,k,ios,ii
      real*8, allocatable :: elo(:,:)

      nloat=1000
      allocate (elo(0:lomax,1:nloat))
      ii=1
      nloat=0
      open(99,file=trim(casename)//".energy")  
      DO i=1,nat
        read(99,*)
        do
          write(6,*) "ii=",ii,"lomax=",lomax
          read(99,'(f12.5)',iostat=ios) elo(0:lomax,1:ii)
          if (ios.ne.0) exit
          ii=ii+1
          backspace(99)
        enddo
        nloat=ii
      ENDDO

      deallocate(elo)
      close(99) 
      end subroutine 

      end subroutine set_lapwlo
!EOC
