      subroutine bandanaly(ib,nb,nkp,kvecs,eband,efermi,nsp,title)
      implicit none
      integer, intent(in) :: ib,nb,nkp,nsp
      real(8),intent(in)::eband(ib:nb,nkp,nsp),efermi,kvecs(3,nkp)
      character(len=*),intent(in):: title


      integer :: i,iat,ik,nc,nv,isp,ikv,ikc     ! Indexes inequivalent atoms
      integer :: fid,info
      integer :: nat,nwf,nbk
      integer :: numin(nsp),nomax(nsp),ikvm(nsp),ikcm(nsp)
      real(8) :: kvec(3),w,kvec0(3)
      real(8) :: evbm,ecbm,egap(3),ebmin,ebmax
      real(8) :: bande(ib:nb,nkp,nsp)
      logical :: lmetal

      real(8),parameter::hev=27.2113961

      call boxmsg(6,'-',trim(title)//" Band Analysis")
      write(6,'(a,2i5)'  ) "  Range of bands considered: ",ib,nb
      write(6,'(a,f10.4)') "  EFermi= ",efermi
      bande=eband
!
! Find Valence Band Maximum (VBM) and Cond. Band Minimum (CBM)
!
      if(maxval(bande).lt.efermi .or. minval(bande).gt.efermi ) then 
        write(6,*) "WARNING from bandanaly: "
        write(6,*) " - Fermi energy outside the energy range of bande!" 
        write(6,'(a,f10.4)') "  Fermi Energy(eV):   ",efermi*hev
        write(6,'(a,f10.4)') " - Maximal energy:    ",maxval(bande)*heV
        write(6,'(a,f10.4)') " - Minimal energy:    ",minval(bande)*heV
        return 
      endif 

      nomax = 1
      numin = nb
      ikvm  = 1
      ikcm  = 1
      lmetal = .false.
      do isp=1,nsp
        do ik = 1, nkp
          nv=0
          nc=nb
          do i=ib,nb
            if(bande(i,ik,isp).lt.efermi)then
              if(i.gt.nv) nv = i 
            else
              if(i.lt.nc) nc = i  
            endif
          enddo
          if(nv.gt.nomax(isp)) nomax(isp) = nv 
          if(nc.lt.numin(isp)) numin(isp) = nc 

          nv = nomax(isp) 
          nc = numin(isp) 
          ikv=ikvm(isp) 
          ikc=ikcm(isp) 
          if(bande(nv,ik,isp).gt.bande(nv,ikv,isp)) ikvm(isp) = ik
          if(bande(nc,ik,isp).lt.bande(nc,ikc,isp)) ikcm(isp) = ik

        enddo
        if(nomax(isp).lt.1) then
          write(6,*) "ERROR: nomax < 1  !!!"
          write(6,*) "  --- Check the Fermi energy"
          stop
        endif

        if(nomax(isp).ge.nb) then
          write(6,*) "ERROR: nomax is larger than the number of available&
     & bands    !!!"
          write(6,*) "  --- Check the Fermi energy"
          stop
        endif

        if(numin(isp).lt.1) then
          write(6,*) "ERROR: numin <1  !!!"
          write(6,*) "  --- Check the Fermi energy"
          stop
        endif

        if(numin(isp).ge.nb) then
          write(6,*) "ERROR: numin >= nb  !!!"
          write(6,*) "  --- Check the Fermi energy"
          stop
        endif

        if(nomax(isp).ge.numin(isp)) then
          write(6,*) "Valence and Conductance bands overlap: metallic!"
          lmetal=.true.
        endif
      enddo

!
!     set evbm
!
      if(lmetal) then
        evbm = efermi 
      else 
        if(nsp.eq.1) then 
          evbm = bande(nomax(1),ikvm(1),1)
        else 
          evbm = max(bande(nomax(1),ikvm(1),1),bande(nomax(2),ikvm(2),2))
        endif 
      endif 
      bande = (bande - evbm)*hev

      do isp=1,nsp
        nv = nomax(isp)
        nc = numin(isp)
        ikv = ikvm(isp) 
        ikc = ikcm(isp) 
        if(nsp.eq.2) write(6,*) "Band analysis for spin",isp    
        write(6,101) "Band index for VBM and CBM=",nv,nc

        egap(1)= bande(nc,ikc,isp)-bande(nv,ikv,isp) 
        egap(2)=minval(bande(nc:nb,ikv,isp))-bande(nv,ikv,isp)
        egap(3)= bande(nc,ikc,isp)-maxval(bande(ib:nv,ikc,isp))

        if(ikvm(isp).eq.ikcm(isp)) then ! direct gap 
          write(6,112) trim(title),egap(1)
          write(6,113) kvecs(:,ikv),ikv
        else 
          write(6,114) trim(title),egap(1:3)
          write(6,115) kvecs(:,ikv),ikv,kvecs(:,ikc),ikc
        endif

        write(6,*) "Range of each band with respect to VBM (eV):"
        write(6,'(a5,3a12)') 'n ','Bottom','Top','Width'
        do i=ib,nb
          ebmin=minval(bande(i,1:nkp,isp))
          ebmax=maxval(bande(i,1:nkp,isp))
          write(6,'(i5,3F12.3)') i,ebmin,ebmax,ebmax-ebmin
        enddo
      enddo
      call linmsg(6,'-','')

  100 format(3e19.12,a10,2i6,f5.1)
  101 format(a,2i4)
 112  format(':BandGap_',a,' =  ',f12.3,' eV')
 113  format("  Direct gap at k=  ",3f8.3,' ik=',i5)
 114  format(':BandGap_',a,' =  ',3f12.3,' eV')
 115  format("  Indirect gap, k(VBM)=",3f8.3,' ik=',i5,/,&
     &       "                k(CBM)=",3f8.3,' ik=',i5)
      end subroutine 
