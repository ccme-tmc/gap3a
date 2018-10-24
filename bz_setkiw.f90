!BOP
!
! !ROUTINE: bz_setkiw
!
! !INTERFACE:
      subroutine bz_setkiw()
      
! !DESCRIPTION:
!
! This subroutine calculates the weights for BZ integrations.
!
!
! !USES:
!
      use bands,     only: nbmax, bande,efermi,nspin,&
     &                     numins, nomaxs, nomax,metallic,numin
      use bzinteg,   only: kiw,kwfer,kwt_bz,kwt_ibz
      use bzint,     only: iop_bzint,bzint_smear,ntet,link,tvol,tnodes,&
     &                     wtet
      use kpoints,   only: nkp,nirkp,kpirind,wkir,idikp
      use struk,     only: nat
      use task,      only: fid_outkpt
      
     
! !LOCAL VARIABLES:

      implicit none

      integer :: iat ! (Counter) Runs over inequivalent atoms
      integer :: ic  ! (Counter) Runs over core states
      integer :: ik,irk,ie
      integer :: isp,nomx,numn
      integer :: ncbm,nkir(nirkp) 
      real(8), allocatable :: enk(:,:) ! Band for which the weights are being calculated
      real(8), allocatable :: wnk(:,:) 
      real(8):: wt
      real(8):: zero_occ = 1.e-4

      logical:: ldbg=.true.
      integer:: ierr
      integer:: fout
      
!
!EOP
!BOC

      allocate(kiw(nbmax,nirkp,nspin),                                &
     &         kwfer(nbmax,nirkp,nspin),                              &
     &         kwt_bz(nkp),                                           &
     &         kwt_ibz(nirkp),                                        &
     &         stat=ierr)
      if(ierr.ne.0) then
        write(6,*) "bz_setkiw: Fail to allocate memory"
        stop
      endif

      kiw=0.0d0
      kwfer=0.0d0
      
      do isp=1,nspin 
        ncbm = numins(isp)

        if(iop_bzint.eq.0) then !! perform BZ integration over tetrahedra in the IBZ
          allocate(enk(nbmax,nirkp))
          do irk=1,nirkp
            enk(1:nbmax,irk)=bande(1:nbmax,irk,isp)
          enddo
          call tetiw(nirkp,nbmax,enk,efermi,kiw(:,:,isp))

          if(metallic) then
            call tetiwsurf(nirkp,nbmax,enk,efermi,kwfer(:,:,isp))
          endif
          deallocate(enk) 

        elseif(iop_bzint.eq.1) then !! Set the weights by the smearing
          allocate(enk(nbmax,nirkp))
          do irk=1,nirkp
            enk(1:nbmax,irk)=bande(1:nbmax,irk,isp)
          enddo
          do irk=1,nirkp
            wt = (1.d0*wkir(irk))/nkp
            do ie=1,nbmax
              kiw(ie,irk,isp)   = wt*bzint_smear(0,enk(ie,irk))
              if(metallic) kwfer(ie,irk,isp) = wt*bzint_smear(1,enk(ie,irk))
            enddo
          enddo
          deallocate(enk) 
!
!         adjust nomaxs and numins 
!
          nomx = nomaxs(isp)
          numn = numins(isp)
          do irk = 1, nirkp
            do ie=1,nbmax
              wt = kiw(ie,irk,isp)*nkp 
              if(wt.gt.zero_occ.and.ie.gt.nomx) nomx = ie 
              if(wt.lt.1.d0-zero_occ .and. ie.lt.numn) numn=ie
            enddo
          enddo
          if(nomx.gt.nomaxs(isp)) then 
            write(6,*) " nomax in terms of the occupation:",nomx
            nomaxs(isp) = nomx
          endif 
          if(numn.lt.numins(isp)) then 
            write(6,*) " numin in terms of the occupation:",numn 
            numins(isp) = numn 
          endif 

        else
          write(6,*) "ERROR: unsupported option iop_bzint=",iop_bzint
          stop
        endif

        do irk=1,nirkp
          kiw(:,irk,isp)=kiw(:,irk,isp)/wkir(irk)  
          if(metallic) kwfer(:,irk,isp) = kwfer(:,irk,isp)/wkir(irk) 
        enddo
      enddo

      nomax = maxval(nomaxs) 
      numin = minval(numins)

      !! set up the band-independent k-weight in the IBZ and BZ
      do irk=1,nirkp
        kwt_ibz(irk)=kiw(1,irk,1)*wkir(irk) 
      enddo
      do ik=1,nkp
        irk=kpirind(ik)
        kwt_bz(ik)=kiw(1,irk,1) 
      enddo

      fout=fid_outkpt 
      if(ldbg) then
        write(fout,100) "ie","irk","wkir","enk","kiw","kwfer"
        do isp=1,nspin 
          do irk=1,nirkp
            write(fout,*) 
            do ie = max(1,nomax-5), min(nomax+5,nbmax)
              write(fout,101) ie,irk,wkir(irk),bande(ie,irk,isp),&
     &           kiw(ie,irk,isp),kwfer(ie,irk,isp)
            enddo
          enddo
        enddo
      endif 

 100  format("#BZ integration weights ",/,"#",a5,2a6,3a12)
 101  format(3i6,10f12.6)

      end subroutine bz_setkiw
!EOC
         
