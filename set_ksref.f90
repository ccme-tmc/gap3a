!BOP
!
! !ROUTINE: set_ksref
!
! !INTERFACE:
      subroutine set_ksref 
      
! !DESCRIPTION:
!
! This subroutine set up Kohn-Sham reference system, including core
! states, KS energies and orbitals, KS exchange-correlation potential 
!
! !USES:
      use bands,    only: nspin,bande,bande0,nomax,numin,nomaxs,numins, &
     &                    nbmax,nbmaxpol,nbmaxsc,emaxpol,emaxsc, &
     &                    efermi,eferks,eferqp,eferqp0,  &
     &                    ibgw,nbgw,nbandsgw,emingw,emaxgw,eqp,eqp_im,  &
     &                    metallic,iop_metallic,spinmom,band_scissor,   &
     &                    bands_getgap,nbands_x,nbands_c 
  
      use core,     only: ncg,eigcore,core_ortho,ncg_x,ncg_c,ncg_p,     &
     &                    iop_core,lcore,ncore
      use constants,only: hev
      use eigenvec, only: init_eigenvec
      use kpoints,  only: nvelgw,nirkp,nvel
      use lapwlo,   only: lmax,nlomax,lvmax_at
      use mommat,   only: init_rad_mom
      use struk,    only: nat
      use task,     only: iop_dftpkg
      use xcpot,    only: uxcu,lxcmax,uiorb,lvorb,natorb 

! !LOCAL VARIABLES:

      implicit none
      
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
      integer :: nwf,nkp_ene,ncbm
      integer    :: ierr 
      real(8) :: kvec(3),w,kvec0(3) 
      real(8) :: eryd,egap
      
      character(10) :: kname
      character(2):: spflag(2)
      character(10), parameter :: sname = 'set_ksref'
!EOP
!BOC
!

!
!     read core states and calculate integrals needed for momentum matrix elements 
!     related to core states
!
      if(iop_dftpkg.eq.0) then 
        call w2k_readcore 

!        !! correct lvmax_at 
!        do iat=1,nat
!          lvmax_at(iat)=max(lvmax_at(iat),maxval(lcore(1:ncore(iat),iat)))
!          write(6,'(a,2i5)') "lvmax for atom",iat,lvmax_at(iat) 
!        enddo 
          
      else 
        write(6,*) "ERROR: unsupported option for dftpkg"
        stop
      endif

      if(core_ortho) then 
        do isp=1,nspin 
          do iat=1,nat 
            call orthog_corewf(iat,isp)
          enddo
        enddo
      endif
!
!    radial integral related with momentum matrix 
!
      call init_rad_mom
      do isp=1,nspin
        do iat = 1, nat
          call momradintc(iat,isp)
          call momradintv(iat,isp)
        enddo
      enddo
!
!     Read KS band energies  and set some band parameters
!
      if(iop_dftpkg.eq.0) then
        call w2k_readenergy(0)
      endif 
      
      !! internal subroutine to set up some band parameters
      call sub_setbandpar   

!
!     Setup KS eigenvectors 
!
      call init_eigenvec     
      if(iop_dftpkg.eq.0) then 
        call w2k_genvectorf
      endif 
!
!     read xc potential data and calculate integrals related to the
!     calculation of vxc matrix elements 
!
      if(iop_dftpkg.eq.0) then  

        call w2k_readvxc

        allocate(uxcu((nlomax+2)**2,(lmax+1)**2,lxcmax,nat,nspin))  
        uxcu=0.d0

        if(lvorb) then 
          allocate(uiorb((nlomax+2)*(nlomax+1),2,natorb,nspin))
          uiorb = 0.d0 
        endif 

        do isp=1,nspin
          do iat = 1, nat
            call w2k_vxccub(iat,isp)
            call w2k_calcuxcu(iat,isp)
          enddo ! iat
        enddo ! isp
      endif 
  
      return

      contains 

!**********************************************************************!
!     internal subroutine --  set up some band parameters              !
!**********************************************************************!
      subroutine sub_setbandpar
      integer :: nbmax_old 

      !! set nbmaxpol, the number of bands used for polmat calculations 
      call sub_set_nbmax(emaxpol,nbmaxpol) 
      write(6,301) "Nr. of bands for polarization(nbmaxpol):",nbmaxpol

      !! set the number of bands used for selfc if emaxsc is set 
      call sub_set_nbmax(emaxsc,nbmaxsc)
      write(6,301) "Nr. of bands for Sc(nbmaxpol):",nbmaxpol
!
!     Set the number of core states considered in the summation over
!     states
!
      ncg_x = ncg 
      if(iop_core.eq.0) then 
        ncg_c = ncg 
        ncg_p = ncg 
      elseif(iop_core.eq.1) then 
        ncg_p = 0
        ncg_c = ncg 
      else 
        ncg_c = 0 
        ncg_p = 0 
      endif 

      write(6,'(a,i5)') " Nr. core states(ncg):            ",ncg
      write(6,'(a,i5)') " Nr. core states for selfx(ncg_x):",ncg_x
      write(6,'(a,i5)') " Nr. core states for selfc(ncg_c):",ncg_c
      write(6,'(a,i5)') " Nr. core states for eps  (ncg_p):",ncg_p

      !! calculate the fermi energy and shift all band energies with
      !! respect to it  
      if(efermi>1.e2) then 
        call fermi(nirkp,nbmax,bande,nvel,nspin,efermi,egap,'KS') 
      else 
        write(6,'(a)') "Use the Fermi energy from case.ingw"
        call bands_getgap(bande,nbmax,nirkp,nspin,efermi,egap) 
      endif 

      if(egap.le.0.d0) then 
        metallic=.true.
      else 
        metallic=.false.
      endif 

      eferks = efermi 
      eferqp0 = efermi

      bande = bande - efermi 

      if(ncg.gt.0) then 
        eigcore=eigcore - efermi 
      endif 

      efermi = 0.d0 
      bande0 = bande


! Determine some band parameters 
!      nomax --- index for highest occupied band  (VBM)
!      numin --- index for lowest unoccupied band (CBM)
!      nbgw --- the number of band states for which GW self-energy are calculated  
      do isp=1,nspin 
        nomaxs(isp)=0
        numins(isp)=1000 
        do ik = 1, nirkp
          io=0
          iu=1000
          do ib=1,nbmax
            if(bande(ib,ik,isp).lt.efermi)then
              if(ib.gt.io)io=ib
            else
              if(ib.lt.iu)iu=ib
            endif    
          enddo
          if(io.gt.nomaxs(isp))nomaxs(isp)=io
          if(iu.lt.numins(isp))numins(isp)=iu
        enddo    
      enddo 
      if(nspin.eq.1) then 
        nomaxs(2)=nomaxs(1)
        numins(2)=numins(1)
      endif 
!
!    Check whether nomax is consistent with nvel 
!
      if(metallic.and.iop_metallic.eq.1) then 
        write(6,*) "WARNING: trying to do an insulating calculation &
     &with metallic initials bands!!!"
        metallic = .false.
        write(6,'(a)') "  reset nomax in terms of electron numbers"
        if(nspin.eq.1) then 
          nomaxs=nint(nvel/2)
        else
          nomaxs(1)=nint(nvel+spinmom)/2
          nomaxs(2)=nint(nvel-spinmom)/2
        endif
        numins=nomaxs+1
        egap = minval(bande(numins(1),:,1))-maxval(bande(nomaxs(1),:,1))

        write(6,'(a,f8.3)') "  val. and cond. bands overlap:",egap
        if(band_scissor<=0.0.and.egap < 0.0 ) band_scissor = abs(egap)+0.01
        if(band_scissor>0.0) then 
          write(6,'(a,f6.3,a)') "  unoccupies states are shifted by",&
     &      band_scissor, " Hartree" 
        endif 
      endif 

      write(6,301)'Highest occupied band:', nomaxs
      write(6,301)'Lowest unoccupied band:', numins 

      nomax=maxval(nomaxs(1:nspin)) 
      numin=minval(numins(1:nspin)) 

!
!     Set the total number of bands considered in the summation over
!     states for the calculations the exchange (x) and correlation
!     self-energies
!
      nbands_x = nomax + ncg_x
      nbands_c = nbmaxsc + ncg_c 

      write(6,301) "Number of bands considered in Sx(nbands_x):",nbands_x
      write(6,301) "Number of bands considered in Sc(nbands_c):",nbands_c


      call errmsg(nomax.lt.1,    sname,"nomax < 1")
      call errmsg(nomax.ge.nbmax,sname,"nomax >= nbmax")
      call errmsg(numin.lt.1,    sname,"numin<1")
      call errmsg(numin.ge.nbmax,sname,"numin>=nbmax")

!     Determine the range of GW bands ibgw..nbgw and the number of valence electrons included in gw bands
      if(ibgw.le.0) then 
        ibgw=1
        nbgw=nbmax 
        do ib=2,nbmax-1 
          if(   (     emingw.gt.maxval(bande(ib-1,:,:))                 &
     &         .and.emingw.le.minval(bande(ib  ,:,:)))                  &
     &      .or.(     emingw.le.maxval(bande(ib  ,:,:))                 &
     &         .and.emingw.ge.minval(bande(ib  ,:,:)))) ibgw=ib 

          if(   (     emaxgw.lt.minval(bande(ib+1,:,:))                 &
     &         .and.emaxgw.ge.maxval(bande(ib  ,:,:)))                  &
     &      .or.(     emaxgw.le.maxval(bande(ib  ,:,:))                 &
     &         .and.emaxgw.ge.minval(bande(ib  ,:,:)))) nbgw=ib
        enddo 

        if(ibgw.ge.numin ) then 
          write(6,*) "set_ksref: WARNING - range of gw bands!!"
          write(6,*)  " ibgw,numin=",ibgw,numin
          write(6,*)  " set ibgw to 1 "
          ibgw = 1
        endif 

        if(nbgw.le.nomax) then
          write(6,*) "set_ksref: WARNING - range of gw bands!!"
          write(6,*)  " nbgw,nomax=",nbgw,nomax
          write(6,*)  " set nbgw to nbmax"
          nbgw=nbmax 
        endif
      else 
        nbgw=min(nbgw,nbmax)
      endif 

      allocate( eqp(ibgw:nbgw,nirkp,nspin), &
     &          eqp_im(ibgw:nbgw,nirkp,nspin))
      
      nbandsgw=nbgw-ibgw+1
      nvelgw=nvel-2.d0*(ibgw-1)         

      write(6,301)'Nr. of bands (nbmax):               ', nbmax
      write(6,301)'Nr. of bands used in P (nmaxpol):   ', nbmaxpol
      write(6,301)'Nr. of gw bands (nbandsgw):         ', nbandsgw
      write(6,301)'Range of GW bands (ibgw,nbgw):      ',ibgw,nbgw
      write(6,301)"Nr. val. electrons (nvel):          ",int(nvel)
      write(6,301)"Nr. val. electrons in GW (nvelgw):  ",int(nvelgw)
  301 format(a,2i4)

      end subroutine 

      subroutine sub_set_nbmax(emx,nmx)
      real(8):: emx
      integer:: nmx
      integer :: isp,ik,ib

      if(emx .le. 0.0d0 ) then
        nmx = nbmax
      else
        nmx= 0
        do isp=1,nspin
          do ik=1,nirkp
            do ib=1,nbmax
              if( bande(ib,ik,isp).le.emx.and.ib.gt.nmx) then
                 nmx = ib
              endif
            enddo
          enddo
        enddo
      endif
      end subroutine 

      end subroutine set_ksref
      
!EOC      
