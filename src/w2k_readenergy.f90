!BOP
!
! !ROUTINE: w2k_readenergy
!
! !INTERFACE:
      subroutine w2k_readenergy(iop) 
      
! !DESCRIPTION:
!
! This subroutine reads the DFT eigenvalues for every k-point
!
! {\bf WIEN2k interface}
!      
! !USES:
      
      use bands,    only: nspin,bande,init_bands,nv,nbmax 
      use constants,only: hev
      use eigenvec, only: lsymvector
      use kmeshintp,only: kvecs2,eks2,nkp2,wk2
      use kpoints,  only: kpirind,nirkp,nkp,nirtet,tvol,tndi,wirtet,nvel,&
     &                    klist,idvk,kirlist,idvkir,idikp,nkp_ene
      use lapwlo,   only: nlo_tot
      use recipvec, only: ngk, ngkir, maxngk
      use struk,    only: nat 
      use task,     only: casename,spflag

! !LOCAL VARIABLES:

      implicit none
      integer,intent(in)::iop  !! == 1, read the band energy for a general k-mesh and store them in arrays defined
                               !!       module kmeshintp
                               !! == 0,  read the band energy for a regular k-mesh 
                               !! == -1, just read maxngk,nbmax and nirkp

      integer,parameter::nkp_max=10000  !! maximal number of k-points that is possible (?!
      
      integer :: iat      ! Indexes inequivalent atoms
      integer :: ik,irk   ! Indexes irreducible k-points
      integer :: ib       ! Indexes bandsa
      integer :: isp      ! Index for spin 
      integer :: i
      integer :: io
      integer :: iu
      integer :: ic
      integer :: fid
      integer :: l
      integer :: nwf,nbk
      real(8) :: kvec(3),w
      real(8) :: en
      
      integer :: ierr 
      character(10) :: kname
      character(10), parameter :: sname = 'w2k_readenergy'
!EOP
!BOC
!

!----------------------------------------------------------------------!
      call linmsg(6,'-',sname)
      fid=999
      open(unit=fid,file=trim(casename)//".energy"//spflag(1),&
     &     action='read',iostat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to open case.energy")


!   extract the following information from case.energy
!      nkp_ene -- the number of k-points 
!      maxngk  -- the maximal number of plane waves used for the
!                 expansion of KS vectors 
!      nbmax   -- number of bands for each k, in wien2k case.energy, 
!                 the number of bands is k-dependent, but to simplify 
!                 to the code, this k-dependence is neglected by 
!                 neglecting additional bands for some k-points
      
      do iat = 1,nat
        read(fid,*)
        read(fid,*)
      enddo

      nkp_ene=0
      maxngk = 0
      nbmax = 100000
      do ik=1,nkp_max
        read(fid,'(67x,2i6,5x)',iostat=ierr) nwf,nbk
        if(ierr.ne.0) exit
        nkp_ene=nkp_ene+1
        if(nbk.lt.nbmax)   nbmax=nbk
        if(nwf.gt.maxngk) maxngk=nwf
        do ib=1,nbk
          read(fid,*)
        enddo
      enddo
      close(fid)
      if(iop.eq.-1) return 

      if(iop.eq.1) then
        nirkp=nkp_ene
        nkp2=nirkp
        allocate(kvecs2(3,nkp2),wk2(nkp2))
      endif
!
!     allocate the vector for the number of bands of each k-point
!     and the array for the band energies
!     Set them by default to 1000
!
      call init_bands(nirkp) 
      allocate(ngkir(nirkp))

      do isp=1,nspin 
        open(unit=fid,file=trim(casename)//".energy"//spflag(isp),&
     &     action='read',iostat=ierr)
        call errmsg(ierr.ne.0,sname,"Fail to open case.energy")

        do iat = 1,nat
          read(fid,*) 
          read(fid,*) 
        enddo
      
        if(lsymvector.or.iop.eq.1) then 
          do ik = 1, nirkp
            read(fid,100,err=907) kvec(1:3),kname,nwf,nbk,w

            if(iop.eq.1) then
               kvecs2(1:3,ik)=kvec(1:3)
               wk2(ik)=w
            endif

            if(isp.eq.1) then 
              ngkir(ik)=nwf
              nv(ik)=nwf - nlo_tot
            endif 

            do ib=1,nbk
              read(fid,*,err=908) i,en
              if(ib.le.nbmax) bande(ib,ik,isp)=en 
            enddo
          enddo
        else 
          do ik = 1, nkp
            read(fid,100,err=907) kvec(1:3),kname,nwf,nbk,w
            irk = kpirind(ik) 
            if(ik.eq. idikp(irk)) then 
              if(isp.eq.1) then 
                ngkir(irk)=nwf
                nv(irk)=nwf - nlo_tot
              endif 
              do ib=1,nbk
                read(fid,*,err=908) i,en
                if(ib.le.nbmax) bande(ib,irk,isp) = en 
              enddo
            else 
              do ib=1,nbk
                read(fid,*,err=908)  
              enddo 
            endif 
          enddo
        endif 
        close(fid)
      enddo ! isp

      !! Convert band energy values to the unit of Hartree unit
      bande=bande*0.5d+0
      return

  100 format(3e19.12,a10,2i6,f5.1)
  907 call outerr(sname,'error reading energy file, code 907')
  908 call outerr(sname,'error reading energy file, code 908')
      end subroutine w2k_readenergy
      
!EOC      
