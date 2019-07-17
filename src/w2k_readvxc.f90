!BOP
!
! !ROUTINE: readvxc
!
! !INTERFACE:
      subroutine w2k_readvxc
      
! !DESCRIPTION:
!
! This subroutine reads the exchange-correlation potential, all xc
! potential data are supposed to be stored in case.vxc, including those
! of vorb in the case of LDA+U calculations 
!
!{\bf WIEN2k interface}
!
! !USES:

      use bands,      only: nspin,ibgw,nbgw 
      use constants,  only: pi,czero
      use kpoints,    only: nkp,kpirind,idikp 
      use xcpot,      only: nslmax,ksxc,lmxc,lxcm,lxcmax,nksxc,vxclm,   &
     &                      vxcs,vorb,natorb,nlorb,lorb,lvorb,lvorb_at, &
     &                      iatorb,lhybrid,vxc_hyb
      use reallocate, only: doreallocate
      use struk,      only: nat,nrpt
      use task,       only: casename  
    
      
! !LOCAL VARIABLES:
      
      implicit none

      integer(4) :: fid      
      integer(4) :: irp ! (Counter), runs over the radial mesh points.
      integer(4) :: isp     ! index for spin 
      integer(4) :: j
      integer(4) :: iat,ia     ! (Counter) Runs over inequivalent atoms
      integer(4) :: il,l,m,m1
      integer(4) :: ikxc    ! (Counter) Runs over k-vectors
      integer(4) :: npt     ! Number of radial mesh points (for each atom the corresponding nrpt(iat) is stored here).
      integer(4) :: nptmax   ! maximum value of nrpt
      integer(4) :: lxc     ! (Counter) l index of the xc component.

      integer :: itmp
      integer :: ierr   ! error index 
      integer :: nb,ie,ik,irk
      real(8) :: rtmp
      character(len=20):: sname="w2k_readvxc"

      real(8), allocatable :: vxclmryd(:)
      complex(8) :: vxcsryd

! !REVISION HISTORY:
!
! Created: 26.07.05 by RGA

!EOP
!BOC
!
!
      call linmsg(6,'-',"READVXC")
      fid=999
      open(unit=fid,file=trim(casename)//".vxc",action='read',iostat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to open case.vxc")

!     Allocate needed arrays
      nptmax=maxval(nrpt)
      allocate(lxcm(nat))
      allocate(lmxc(2,nslmax,nat))
      allocate(vxclm(nptmax,nslmax,nat,nspin))
      vxclm=0.0d0

      allocate(vxclmryd(nptmax))
!
!     Read the MT-Sphere part
!
      isp=1
      lxcmax=0
      read(fid,10)
      do iat=1, nat      
        npt=nrpt(iat)
        read(fid,11) lxcm(iat)
        if(lxcm(iat).gt.nslmax) call outerr(sname,'nslmax too small')
        if(lxcm(iat).gt.lxcmax) lxcmax=lxcm(iat)
        do lxc=1, lxcm(iat)
          read(fid,12,end=1) lmxc(1,lxc,iat),lmxc(2,lxc,iat)
    1     read(fid,13) (vxclmryd(irp),irp=1,npt)
          read(fid,14)
          vxclm(1:npt,lxc,iat,isp)=0.5d0*vxclmryd(1:npt)
        enddo ! lxc
        read(fid,15)

      enddo ! iat
!
!     Read the interstitial part
!
      read(fid,16) nksxc
      allocate(ksxc(1:3,1:nksxc))
      allocate(vxcs(1:nksxc,nspin) )
      do ikxc=1, nksxc
        read(fid,17,end=2)(ksxc(j,ikxc),j=1,3),vxcsryd
        vxcs(ikxc,isp)=0.5*vxcsryd
    2 enddo ! ikxc

!
! Read vxc for spin down if necessary
!
      if(nspin.eq.2) then 
        isp=2
        read(fid,10)
        do iat=1, nat
          npt=nrpt(iat)
          read(fid,11) itmp
          do lxc=1, lxcm(iat) 
            read(fid,12,end=21) itmp,itmp 
   21       read(fid,13) (vxclmryd(irp),irp=1,npt)
            read(fid,14)
            vxclm(1:npt,lxc,iat,isp)=0.5d0*vxclmryd(1:npt)
          enddo ! lxc
          read(fid,15)
        enddo ! iat
!
!     Read the interstitial part
!
        read(fid,16) itmp
        do ikxc=1, nksxc
          read(fid,17,end=22) itmp,itmp,itmp,vxcsryd
          vxcs(ikxc,isp)=0.5*vxcsryd
   22   enddo ! ikxc
      endif 

      if(lhybrid) then 
        read(fid,*) 
        read(fid,*) 
        do isp=1,nspin  
          do ik=1,nkp 
            read(fid,*) rtmp,rtmp,rtmp,nb
            irk = kpirind(ik)
            do ie=1,nb
              read(fid,*) itmp,rtmp 
              if(ik.eq. idikp(irk).and.ie.ge.ibgw.and.ie.le.nbgw) &
     &         vxc_hyb(ie,irk,isp) = rtmp*0.5d0
            enddo
          enddo
        enddo
      endif 
      close(fid)
!
!     reallocate the arrays to their actual size
!
      call doreallocate(lmxc,2,lxcmax,nat)
      call doreallocate(vxclm,nptmax,lxcmax,nat,nspin)
!
! Read vorb if necessary 
!

      if(lvorb) then 
        open(fid,file=trim(casename)//".vorb",action='read',iostat=ierr)
        call errmsg(ierr.ne.0,sname,"Fail to open case.vorb")
        vorb=czero      

        do isp=1,nspin 
          read(fid,*)      
          do ia=1,natorb 
            read(fid,*) 
            do il=1,nlorb(ia) 
              read(fid,*) 
              l=lorb(il,ia) 
              do m=-l,l
                do m1=-l,l
                  read(fid,103) vorb(m,m1,il,ia,isp)
                enddo
              enddo  
            enddo ! il 
          enddo ! ia
        enddo ! isp
        vorb = vorb * 0.5d0    !! change the unit from Ry. to Hartree  
        close(fid)
      endif  ! if_ldau

!
! Print some diagonasis information 
!
      write(6,*) "### Number of lm component for each atom "
      write(6,'(a5,a10)') " iat ","lxcm(iat)"  
      do iat=1,nat
        write(6,'(i5,i10)') iat,lxcm(iat) 
      enddo 
      write(6,*) "# Number of Institial plane wave compoents=",nksxc
    
   10 format(//)
   11 format(/,15x,i3,//)
   12 format(15x,i3,5x,i2,/)
   13 format(3x,4e19.12)
   14 format(/)
   15 format(///)
   16 format(//,13x,i6)    
   17 format(3x,3i5,2e19.12)
  103 format(2f18.11)
   
      return
      
      end subroutine w2k_readvxc
!EOC   
