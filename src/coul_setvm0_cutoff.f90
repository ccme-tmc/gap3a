!BOP
!
! !ROUTINE: coul_setvm0
!
! !INTERFACE:
      subroutine coul_setvm0_cutoff(iq, icutoff)

! !DESCRIPTION:
!
! This subroutine calculates the matrix of the bare coulomb potential according to equations
! \ref{coulsphdef},\ref{deficouloc}, \ref{coulijdefin} and
! \ref{ipw-coul-sphere}.
!
! !USES:

      use barcoul
      use constants,  only: cone, czero, imag, pi,fourpi,twopi
      use kpoints,    only: idvq, qlist
      use mixbasis,   only: lmbmax,bigl,maxbigl,mbl,nmix,nmixmax,rtl,   &
     &                      rrint,locmatsiz,mbsiz,mpwipw,wi0,          &
     &                      mbsiz
      use recipvec,   only: gindex, gqleng, indgq, indgqlen, ngqbarc,   &
     &                      ngq, ngqlen,maxngqlen
      use struk,      only: mult,ndf,nrpt,vi,pos,nat,rotij,rotloc,alat,rbas
      use task,       only: fid_outdbg
!
! !INPUT PARAMETERS: 

      implicit none

      integer, intent(in) :: iq  ! index of the q-point
      integer, intent(in) :: icutoff ! the cut-off option. See barcoul

!
! !LOCAL VARIABLES:

      integer :: iat       ! Indexes the atoms (inequivalent) for the columns of barc
      integer :: idf       ! Indexes all the atoms  for the columns of barc
      integer :: ierr
      integer :: ieq       ! Indexes the equivalent atoms for the columns of barc
      integer :: igl
      integer :: iipw      ! Indexes the basis function in the  interstitial region (columns of barc; im(ipw) = locmatsiz + iipw)
      integer :: ijdf      ! joint idex for (idf,jdf) pairs for compressed storage
      integer :: ijrm      ! joint index for (irm,jrm) pairs for  compressed storage
      integer :: im      ! Indexes the mixed wave function for the columns of barc (barc(im,jm))
      integer :: ipw       ! Indexes the plane waves 
      integer :: irm       ! Indexes the mixed basis functions of atom iat
      integer :: jat       ! Indexes the atoms (inequivalent) for the files of barc
      integer :: jdf       ! Indexes all the atoms for the files of barc
      integer :: jeq       ! Indexes the equivalent atoms for the files of barc
      integer :: jm      ! Indexes the mixed wave function for the files of barc (barc(im,jm)).
      integer :: jrm       ! Indexes the mixed basis functions  of atom jat
      integer :: ipw0     ! Starting value of the index for G vectors. =2 if iq=1 (if q=0, discard G=0)
      integer :: l1        ! Angular momentum quantum number of the  atomic function irm
      integer :: l2        ! Angular momentum quantum number of the atomic function jrm
      integer :: lm34      ! joint index for the spherical harmonic Y(l1+l2,m3+m4)
      integer :: m1        ! z-component angular momentum quantum number of the atomic function irm 
      integer :: m2        ! z-component angular momentum quantum number of the atomic function jrm
      integer :: m3,m4     ! z-component angular momentum quantum  number of the atomic functions irm,jrm  after rotation
      integer :: ylsize    ! Size of the array for spherical harmonics
      integer, dimension(3) :: iqvec ! q-vector integer coordinates  also used for the G-vector
      
      integer :: kmix   ! Counter: Runs over mixed basis functions.
      integer :: immax

      real(8) :: gpr          ! Scalar product G.r
      real(8) :: prefac       ! multipicative prefactor
      real(8) :: qglen,qg1len ! |qg1|
      real(8) :: qgxy        ! (G+q)_xy
      real(8) :: qgz         ! (G+q)_z
      real(8) :: tg           ! value of \tilde{g}
      real(8) :: minu
      real(8) :: vqg
      
      real(8), dimension(3) :: gvec ! G-vector of the reciprocal  lattice.
      real(8), dimension(3) :: qg1  ! The sum q+G_1
      real(8), dimension(3) :: qg2  ! qg1 rotated to the internal coordinates of the corresponding equivalent atom
      real(8), dimension(3) :: qg3  ! qg2 rotated to general coordinates
      real(8), dimension(3) :: qvec ! q-point for which the bare  coulomb potential is calculated.
      real(8), allocatable :: rtlij(:,:,:,:)
      complex(8) :: exev      ! exact eigenvalues of the coulomb matrix
      complex(8) :: dj3       ! the rotation matrix d*(l1,m3,m1)
      complex(8) :: dj4       ! the rotation matrix d(l2,m4,m2)
      complex(8) :: expg      ! e^(iG.r)       
      complex(8) :: stc,sum34      

      complex(8), allocatable :: sph(:,:)  ! spherical harmonics for all G's
      complex(8), allocatable :: mat1(:,:),mat2(:,:),vtemp(:,:)

!! debug 

      logical :: lprt =.true.
      logical :: ldbg =.true.
      character(len=20)::sname='coul_setvm0_cutoff'
  
 
! !EXTERNAL ROUTINES: 

      
      real(8), external :: getcgcoef
      real(8), external :: vecprojlen
      complex(8), external :: getdjmm
      
      external calcjlam
      external k2cart
      external rotate
      external ylm
      external zgemm
!
! !REVISION HISTORY:
! 
! Created 16th. March 2004 by RGA
! modified 31. March 2005 by RGA
! modified 08. Sept  2019 by MYZ
!
!EOP
!BOC
     
      if(lprt) call linmsg(6,'-','coul_setvm0_cutoff (2D)')
      
      if(lprt) write(6,*) "coul_setvm0_cutoff: allocate space"     
      ylsize=(maxbigl+1)*(maxbigl+1)
      allocate(vtemp(1:ngqbarc(iq),1:locmatsiz),   &
     &         rtlij(nmixmax,nat,nmixmax,nat),     &
     &         sph(ylsize,ngqbarc(iq)),            &
     &         mat2(1:ngq(iq),1:locmatsiz),        &
     &         jlam(nmixmax,maxngqlen),            &
     &         stat=ierr ) 
      call errmsg(ierr.ne.0,sname,"Fail to allocate local arrays")
      vmat=0.d0 

      rtlij(:,:,:,:)=0.0d0
      prefac=4.d0*pi*sqrt(vi)

      !! Set ipw0: Starting value of the index for G vectors. =2 if iq=1 
      !! (if q=0, discard G=0)
      ipw0=1
      if(iq.eq.1)  ipw0=2

      !! Calculate the matrix sigma for the structure constants
      if(lprt) write(6,*) "Calculate the matrix sigma, cutoff", icutoff
      !call coul_strcnst(iq, 4*(lmbmax+1))
      call coul_strcnst_cutoff(iq, 4*(lmbmax+1), icutoff)

      !! calculate the cartesian coordinates of the q-point
      iqvec(1:3)=qlist(1:3,iq)
      call k2cart(iqvec,idvq,qvec)

      !! For q=0 calculate the corrections of the structure constants
      if(lprt) write(6,*)" For q=0 calc the corr. of struc. const"
      ! TODO modify coul_barcq0_cutoff
      if(iq.eq.1) call coul_barcq0_cutoff(icutoff)

      !! Calculate all the products rtl*rtl
      if(lprt) write(6,*) " Calculate all the products rtl*rtl"
      do iat = 1, nat
        do irm = 1, nmix(iat)
          do jat = 1, nat
            do jrm = 1, nmix(jat)
              rtlij(jrm,jat,irm,iat)=rtl(irm,iat)*rtl(jrm,jat)
            enddo ! jrm
          enddo ! irm
        enddo ! jat
      enddo ! iat

      if(lprt) write(6,*) " Loop over inequivalent atoms:"
      idf=0
      im=0
      do iat = 1, nat    !* Loop over inequivalent atoms:
      
        if(lprt) write(6,*) " calculate the matrix elements jlam"  
        call calcjlam(iat,mbl(iat),nrpt(iat),iq)

        if(lprt) write(6,*) " Loop over equivalent atoms"
        do ieq = 1, mult(iat)    !* Loop over equivalent atoms:
          idf=idf+1
          ijdf=(idf*(idf+1))/2
!
!         Calculate Y_lm(q+G) for all G
!
          do iipw=1,ngqbarc(iq)
            iqvec(1:3)=gindex(:,indgq(iipw,iq))
            call k2cart(iqvec,1,gvec)
            qg1(1:3) = qvec(1:3) + gvec(1:3)
            call rotate(qg1,rotij(1:3,1:3,idf),qg2)
            call rotate(qg2,rotloc(1:3,1:3,iat),qg3)
            call ylm(qg3,maxbigl,sph(1:ylsize,iipw))
          enddo

          do irm = 1, nmix(iat)  !* Loop over mixed functions:
            l1=bigl(irm,iat)
            do m1=-l1,l1
              im = im + 1
              jm=0  
              jdf=0    
              do jat=1,nat
                do jeq=1,mult(jat)
                  jdf=jdf+1
                  if(jdf.ge.idf)then
                    ijdf=idf+(jdf*(jdf-1))/2
                  else
                    ijdf=jdf+(idf*(idf-1))/2
                  endif 
                     
                  do jrm = 1, nmix(jat)
                    l2=bigl(jrm,jat)
                    do m2=-l2,l2
                      jm=jm+1
                      if((iq.ne.1).or.(l1.ne.0).or.(l2.ne.0)) then
                        sum34=czero
                        do m3=-l1,l1
                          do m4=-l2,l2
                            tg=gettildeg(l1,l2,m3,m4)
                            dj3=conjg(getdjmm(idf,l1,-m3,m1))
                            dj4=getdjmm(jdf,l2,m4,m2)
                            if(jdf.ge.idf)then
                              lm34=(l1+l2)*(l2+l1+1)+m4+m3+1
                              minu=(-1.0d0)**m3
                              stc=tg*sgm(lm34,ijdf)
                            else
                              lm34=(l1+l2)*(l2+l1+1)-m4-m3+1
                              minu=(-1.0d0)**(m4+l1+l2)
                              stc=tg*conjg(sgm(lm34,ijdf))
                            endif
                            sum34=sum34+minu*stc*dj3*dj4
                          enddo ! m4
                        enddo ! m3  
                        vmat(im,jm)=rtlij(jrm,jat,irm,iat)*sum34

                        if((m1.eq.m2).and.(l1.eq.l2).and.(idf.eq.jdf))then
                          if(jrm.ge.irm)then
                            ijrm=irm+(jrm*(jrm-1))/2
                          else
                            ijrm=jrm+(irm*(irm-1))/2
                          endif  
                          vmat(im,jm)=vmat(im,jm)+cmplx(fourpi &
     &                       *rrint(ijrm,iat)/dble(2*l1+1),0.0d0,8)
                        endif
                      endif ! if((iq.ne.1).or.(l1.ne.0).or.(l2.ne.0)) 
                      ! q=l1=l2=0 already calculated in coul_barcq0
                    enddo ! m2
                  enddo ! jrm
                enddo ! jeq
              enddo ! jat  
!
!             Calculation of the matrix element between an atomic 
!             mixed function and an IPW
!
              if(ipw0.eq.2) vtemp(1,im) = czero
              do ipw = ipw0, ngqbarc(iq)
                iqvec(1:3)=gindex(:,indgq(ipw,iq))
                call k2cart(iqvec,1,gvec)
                gpr=dble(iqvec(1))*pos(1,idf)+ &
     &              dble(iqvec(2))*pos(2,idf)+ &
     &              dble(iqvec(3))*pos(3,idf)
                expg=cmplx(cos(2.0d0*pi*gpr),-sin(2.0d0*pi*gpr),8)
                igl=indgqlen(ipw,iq)
                qglen=gqleng(igl,iq)
                qg1(1:3) = qvec(1:3) + gvec(1:3)
                qgxy=vecprojlen(qg1,rbas(axis_cut_coul,:),'perp')
                qgz=vecprojlen(qg1,rbas(axis_cut_coul,:),'para')
                vqg=fourpi/(qglen*qglen) * &
                    (1.0d0-exp(-qgxy*zcut_coul)*cos(qgz*zcut_coul))
                !if(ldbg) then
                !  write(*,"(7F15.6)") qg1(1:3),qgxy,qgz,fourpi/(qglen*qglen),vqg
                !endif
                vtemp(ipw,im)=cmplx(prefac*jlam(irm,igl)*vqg,0.0d0,8)   &
     &                *sph(l1*(l1+1)+m1+1,ipw)*((-imag)**l1)*expg
                !TODO why minus imag?
              enddo ! ipw       
            enddo ! m1  
          enddo ! irm
        enddo ! ieq
      enddo ! iat

      call zgemm('n','n',ngq(iq),locmatsiz,ngqbarc(iq),cone,mpwipw,     &
     &           ngq(iq),vtemp,ngqbarc(iq),czero,mat2,ngq(iq))
      do iipw=1,ngq(iq)
        do im=1,locmatsiz
          vmat(locmatsiz+iipw,im)=mat2(iipw,im)
          vmat(im,locmatsiz+iipw)=conjg(mat2(iipw,im))
        enddo
      enddo    
      deallocate(mat2,jlam,rtlij,sph,vtemp)
!
!     Calculation of the matrix elements between two IPW's
!
      allocate(mat1(1:ngq(iq),1:ngqbarc(iq)), &
     &         mat2(1:ngq(iq),1:ngq(iq)), &
     &         stat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to allocate mat1,mat2")

      if(lprt) write(6,*) "Calc matrix elements between two IPW's"

      mat1(:,1)=czero 
      do ipw = ipw0, ngqbarc(iq)     
        iqvec(1:3)=gindex(:,indgq(ipw,iq))
        call k2cart(iqvec,1,gvec)
        qg1(1:3) = qvec(1:3) + gvec(1:3)
        qgxy=vecprojlen(qg1,rbas(axis_cut_coul,:),'perp')
        qgz=vecprojlen(qg1,rbas(axis_cut_coul,:),'para')
        igl=indgqlen(ipw,iq)
        qglen=gqleng(igl,iq)
        vqg=fourpi/(qglen*qglen) * &
            (1.0d0-exp(-qgxy*zcut_coul)*cos(qgz*zcut_coul))
        do iipw=1,ngq(iq)
          mat1(iipw,ipw)=mpwipw(iipw,ipw)*cmplx(vqg,0.0d0,8)
        enddo  
      enddo
      call zgemm('n','c',ngq(iq),ngq(iq),ngqbarc(iq),cone,mat1,ngq(iq), &
     &           mpwipw,ngq(iq), czero,mat2,ngq(iq))
      do ipw=1,ngq(iq)
        do iipw=1,ngq(iq)
          vmat(locmatsiz+ipw,locmatsiz+iipw)=mat2(ipw,iipw)
        enddo  
      enddo
      deallocate(mat1,mat2)

      end subroutine coul_setvm0_cutoff
