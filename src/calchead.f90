!BOP
!
! !ROUTINE: calchead
!
! !INTERFACE:
      subroutine calchead(ikfirst,iklast,iomfirst,iomlast)

!
! !DESCRIPTION: 
!
!This subroutine calculate the dielectric matrix at the $\Gamma$ point.      
!

! !USES:

      use bands,       only: bande,nbmaxpol,fspin,    & 
     &                       nomaxs,numins,nspin,metallic
      use bzinteg,     only: kcw,kwfer
      use constants,   only: czero, pi,cone,hev
      use core,        only: corind, eigcore, ncg_p, iop_core 
      use freq,        only: omega,iop_freq
      use kpoints,     only: idikp, nirkp, nkp, wkir,kpirind
      use mommat,      only: mmatvv,mmatcv
      use dielmat,     only: head,c0_head,q0_eps,mask_eps,iop_drude,&
     &                       omega_plasma,eta_head 
      use struk,       only: vi
      use task,        only: fid_outgw, fid_outdbg
      use anisotropy,  only: ten_p_ani, vec_u_ani, iop_aniso
      use mixbasis,    only: matsiz

      implicit none
      integer,intent(in):: ikfirst,iklast
      integer,intent(in):: iomfirst,iomlast
      
! !LOCAL VARIABLES:
      logical:: ldbg=.true.
      character(10):: sname="calchead"

      integer :: icg  ! index for core states of all atoms 
      integer :: ic   ! index for core states on one atom 
      integer :: iat  ! index for inequivalent atoms
      integer :: ie,ie1,ie2   ! Counter, runs over occupied and unoccupied bands,respectively 
      integer :: ie12 ! Joint index for (ie1,ie2) pairs for compressed storage
      integer :: ik,irk,iik ! Counter, runs over kpoints in the full and irreducible BZ
      integer :: iom       ! Counter, runs over frequencies
      integer :: i,j,k ! Counter: Cartisian directions
      integer :: ie12max
      integer :: isp 
      integer :: nvbm,ncbm
      integer :: kwt
     
      real(8) :: om      ! temporay variable for frequency value  
      real(8) :: wpl2    
      real(8) :: cpivi   ! 4*pi/vcell
      real(8) :: edif    ! energy differences
      real(8) :: edsq    ! edif^2
      real(8) :: coef,pnmkq2
      complex(8):: pnmkq 

      complex(8), allocatable :: p_ani_im(:,:,:),u_ani_im(:,:,:)
      complex(8), allocatable :: termcv(:),termvv(:),vwe(:),vwc(:)

!
! !EXTERNAL ROUTINES: 
!
      complex(8), external :: zdotu, ten_rvctrv
      real(8), external :: veclen

!
! !INTRINSIC ROUTINES: 
!
      intrinsic abs
!
! !REVISION HISTORY:
!
! Created 11.02.05 by RGA
! Last modified Aug. 13, 2010 by JH
!
!EOP
!BOC

!TODO: 
!  -- Add and check the treatment of metallic systems 
!

      if(ldbg) call linmsg(6,'-',sname)
      do isp=1,nspin 
        nvbm = nomaxs(isp)
        ncbm = numins(isp) 

        allocate(termvv(1:nbmaxpol*nvbm))
        allocate(vwe(1:nbmaxpol*nvbm))
        allocate(termcv(ncbm:nbmaxpol))
        allocate(p_ani_im(3,3,ncbm:nbmaxpol))
        allocate(u_ani_im(3,matsiz,ncbm:nbmaxpol))
        allocate(vwc(ncbm:nbmaxpol))

        do irk=ikfirst,iklast

          if(irk.eq.1.and.isp.eq.1) then 
            head(iomfirst:iomlast) = cone 
          endif 

          ik = idikp(irk)
          kwt = wkir(irk)

          call crpa_setmask(irk,irk,isp)

          coef=4.0d0*pi*vi*kwt*fspin

          if(ldbg) write(fid_outgw,*) " - coef=",coef
!
!         the part related to metallic systems
!
          if(metallic) then 
            if(ldbg) write(fid_outgw,*) " - intraband transition"
            if(ldbg) write(fid_outgw,*) " - ncbm,nvbm=",ncbm,nvbm 
            do ie=ncbm,nvbm 
              pnmkq2=sum(abs(mmatvv(1:3,ie,ie,irk,isp))**2)/3.d0
              c0_head=c0_head+coef*kwfer(ie,irk,isp)*pnmkq2*mask_eps(ie,ie+ncg_p)
            ! TODO: anisotropy for intraband transition
            enddo 
          endif 
!
!         core-valence
!
          if(iop_core.eq.0)then

            if(ldbg) write(fid_outgw,*) " - core states"

            do icg=1,ncg_p
              iat=corind(1,icg)
              ic=corind(3,icg)
              do ie2=ncbm,nbmaxpol
                edif=bande(ie2,irk,isp)-eigcore(ic,iat,isp)
                edsq=edif*edif

                if(iop_aniso.ne.-1) then
                  !write(*,*) "use ten_p for ani"
                  do i=1,3
                    do j=1,3
                      p_ani_im(i,j,ie2) = mmatcv(j,icg,ie2,irk,isp) * &
     &                          conjg(mmatcv(i,icg,ie2,irk,isp)) / &
     &                          cmplx(edsq,0.0D0,8)
                    enddo
                  enddo
                else
                !! two treatments are equivalent with q0_eps=(1,1,1)/\sqrt{3}
                !! old treatment concerning q -> 0 : averaging over three directions
     !            termcv(ie2)=(1.d0/(3.0d0*edsq))* &
     !&                           sum(abs(mmatcv(:,icg,ie2,irk,isp))**2)
                 !! new treatment : choosing a particular direction 
                  pnmkq=sum(mmatcv(1:3,icg,ie2,irk,isp)*q0_eps(1:3))
                  termcv(ie2)=abs(pnmkq)**2/edsq/veclen(q0_eps)**2
                endif
              enddo ! ie2
              !write(*,*) "before iom with iomfisrt = ", iomfirst
              !write(*,*) "                iomlast  = ", iomlast

              do iom=iomfirst,iomlast
                vwc(ncbm:nbmaxpol)=kcw(icg,ncbm:nbmaxpol,ik,iom,isp)
                if(iop_aniso.ne.-1) then
                  do i=1,3
                    do j=1,3
                      ten_p_ani(i,j,iom) = &
     &                 zdotu(nbmaxpol-ncbm+1,p_ani_im(i,j,:),1,vwc(:),1)
                    enddo
                  enddo
                  ten_p_ani(:,:,iom) = cmplx(coef,0.0D0,8)*ten_p_ani(:,:,iom)
                  !write(*,*) "out ie2 with iom = ", iom
                  head(iom)=head(iom)-ten_rvctrv(3,ten_p_ani(:,:,iom),q0_eps)
                else
                  head(iom)=head(iom)-coef*zdotu(nbmaxpol-ncbm+1,vwc,1,termcv,1)
                endif
              enddo ! iom 
              !write(*,*) "after iom"
            enddo ! icg
          endif ! iop_core .eq.  0
!
!     Valence-valence
!
          if(ldbg) write(fid_outgw,*) " - normal band states"
          ie12=0
          do ie2=ncbm,nbmaxpol
            do ie1=1,nvbm
              ie12=ie12+1
              edif=bande(ie2,irk,isp)-bande(ie1,irk,isp)
              edsq=edif*edif
         
!              call wrnmsg(edsq.lt.1.e-10,sname,"degenerate CB and VB")
!              call wrnmsg(edif.lt.0.0,sname,"Wrong ordered VB and CB")

              if(edsq.lt.1.0d-20)then
                termvv(ie12) = czero  
                if(ldbg) then 
                  write(fid_outgw,*) "#Check: nearly degenerate bands: kcw=", & 
     &              kcw(ncg_p+ie1,ie2,ik,1,isp)
                endif 
 
              else 

!               pnmkq=sum(mmatvv(1:3,ie1,ie2)*q0_eps(1:3))
!               termvv(ie12)=fspin*abs(pnmkq)**2/edsq

                !! old 
                termvv(ie12) = mask_eps(ie2,ie1+ncg_p)/(3.0d0*edsq)        &
     &             *sum(abs(mmatvv(:,ie1,ie2,irk,isp))**2)
              endif

            enddo ! ie2
          enddo ! ie1
          ie12max=ie12

          do iom=iomfirst,iomlast

            ie12=0
            do ie2=ncbm,nbmaxpol
              do ie1=1,nvbm
                ie12=ie12+1
                vwe(ie12) = kcw(ncg_p+ie1,ie2,ik,iom,isp)
              enddo 
            enddo 

            head(iom)=head(iom)-coef*zdotu(ie12max,termvv,1,vwe,1)
          enddo
        enddo ! irk 

        deallocate(termcv,termvv,vwc,vwe,p_ani_im,u_ani_im)
      enddo ! isp 

      !! add the plasmon contributions
      write(fid_outdbg, *) "metallic : ", metallic
      write(fid_outdbg, *) "iop_drude : ", iop_drude
      if(metallic.and.iop_drude.eq.1) then
        write(fid_outgw,*) " Intraband contribution!"
        write(fid_outgw,'(a,f12.4)')" Calc. plasmon freq. (eV):",sqrt(c0_head)*hev

        if(omega_plasma.gt.0.0) then 
          wpl2 = omega_plasma**2
          write(fid_outgw,'(a,f12.6)') " Using input plasmon freq. (eV)",omega_plasma*hev
        else
          wpl2 = c0_head 
        endif 

        do iom=iomfirst,iomlast
          om = omega(iom)
          if(abs(om).lt.1.e-20) then
            write(fid_outgw,*)"  WARNING - om = 0.0 is included"
            write(fid_outgw,*)"  -- set it to 1.e-20!"
            om = 1.e-20
          endif
!old 
!          if(iop_freq.eq.2) then
!            head(iom) = head(iom) - cmplx(wpl2/om**2,0.d0)
!          else
!            head(iom) = head(iom) + cmplx(wpl2/om**2,0.d0)
!          endif
! new 
           if(iop_freq.eq.3) then ! imaginary freq.
             head(iom) = head(iom)+wpl2*(1.d0/(om+eta_head)**2)
           elseif(iop_freq.eq.2) then ! real freq.  
             head(iom) = head(iom)+wpl2*cmplx(-1.d0/(om**2+eta_head**2), &
     &             eta_head/(om*(om**2+eta_head**2)))
           endif
           
        enddo 
      endif

      end subroutine calchead
!EOC      
