!OP
!
! !ROUTINE: calcminm
!
! !INTERFACE:
      subroutine calcminm(ik,iq,nstart,nend,mstart,mend,isp,minm)

! !DESCRIPTION:
!
!This subroutine calculates the matrix elements $M^i_{nm}(\vec{k},\vec{q})$ 
!of equation \ref{Mdef} for n, m both belongings to LAPW states 
!
! !USES:
      use bands,       only: nv
      use barcoul,     only: barcvm,im_g0,iop_coul,barcevsq
      use constants,   only: czero, cone, pi
      use core,        only: ncore
      use eigenvec,    only: alfa, beta, gama, alfp, betp, gamp,     &
     &                       zzq, zzk
      use kpoints,     only: klist, qlist, idvk, idvq, kpirind,kqid
      use lapwlo,      only: lmax,lomax,nlo_at
      use mixbasis,    only: nmix,bigl,indggq,locmatsiz,matsiz,mpwipw,  &
     &                       mbsiz,s3r,lo_index_s3r 
      use recipvec,    only: ig0, ngk, indgk, gindex, ngmax,    &
     &                       ngmin, ngq, indgq, ngqbarc, maxngk
      use struk,       only: mult, vi, pos,nat
      use task,        only: time_lapack,time_minm

! !INPUT PARAMETERS:

      implicit none

      integer, intent(in) :: ik  ! the order index of kvec
      integer, intent(in) :: iq  ! index for q-point 
      integer, intent(in) :: nstart,nend,mstart,mend  !! ranges of n and m, respectively
      integer, intent(in) :: isp     ! Index for spin 
      complex(8), intent(out):: minm(matsiz,mstart:mend,nstart:nend)

! !LOCAL VARIABLES:

      integer :: jk
      integer :: bl      ! Angular momentum quantum number of the mixed atomic functions (big l)
      integer :: bm      ! z-component angular momentum quantum number of the mixed atomic function (big m)
      integer :: iat     ! index for inequivalent atoms.
      integer :: icg     ! index over all core states (including a,l,m) 
      integer :: idf     ! index over all atoms, including equivalent ones
      integer :: ie1,ie2 ! index over the eigenvalues in  the corresponding k-point.
      integer :: ieq     ! index for atoms over the equivalent atoms of each type.
      integer :: imix    ! index for mixed basis over all MT-sphere mixed basis functions (of all the atoms)
      integer :: irm     ! index for the radial mixed basis functions of each atom.a
      integer :: iv1,iv2 ! index for the (L)APW basis  functions (excluding core states).
      integer :: l1      ! l-quantum number of the (L)APW eigenfunction at k
      integer :: l2      ! l-quantum number of the (L)APW eigenfunction at k' = k - q
      integer :: m1      ! m-quantum number of the (L)APW eigenfunction at k 
      integer :: m2      ! m-quantum number of the (L)APW eigenfunction at k' = k - q
      integer :: l2min   ! Lowest allowed value of lambda = | bl - l2 |
      integer :: l2max   ! Highest allowed value of lambda = bl + l2
      integer :: lm1     ! Counter: Runs over MTS (L)APW basis  functions (l1,m1) 
      integer :: lm2     ! Counter: Runs over MTS (L)APW basis   functions l2,m2)
      integer :: iipw    ! Counter: Runs over IPW's
      integer :: la1,lb1,lc1,la2,lb2,lc2,i,jpw
      integer :: ilo1,ilo2
      integer :: maxig
      integer :: iz1,iz2
      real(8) :: sumzz

      integer, dimension(3) :: ikvec,igc ! Indexes of G_1+G'-G
      integer ::ierr
      integer :: nv1, nv2
      integer :: ndim,mdim,nmdim
      
      real(8) :: sqvi,x, arg
      real(8) :: qvec(3),kqvec1(3),kqvec2(3)
      real(8) :: time1,time2,tstart,tend

      complex(8) :: angint ! Angular integral calculated by fourylm
      complex(8) :: suml12,sumlm2
      complex(8) :: phs,phs1,phs2
      complex(8) :: sm
      complex(8) :: na,nb,nc

! The following are terms in Eq. \ref{dbracketorg}
      complex(8) :: aa  !  A^{a*}_{l1,m1}(k-q+K')*A^a_{\l2\m2}(k+G)<bl,l2,l1|l2> 
      complex(8) :: ab  !  A^{a*}_{l1,m1}(k-q+K')*B^a_{l2,m2}(k+G)<bl,l2,l1|l2> 
      complex(8) :: ac  !  A^{a*}_{l1,m1}(k-q+K')*C^a_{l2,m2}(k+G)<bl,l2,l1|l2>
      complex(8) :: ba  !  B^{a*}_{l1,m1}(k-q+K')* A^a_{l2,m2}(k+G)<bl,l2,l1|l2> 
      complex(8) :: bb  !  B^{a*}_{l1,m1}(k-q+K')*B^a_{l2,m2}(k+G)<bl,l2,l1|l2> 
      complex(8) :: bc  !  B^{a*}_{l1,m1}(k-q+K')* C^a_{l2,m2}(k+G)<bl,l2,l1|l2> 
      complex(8) :: ca  !  C^{a*}_{l1,m1}(k-q+K')*A^a_{l2,m2}(k+G)<bl,l2,l1|l2> 
      complex(8) :: cb  !  C^{a*}_{l1,m1}(k-q+K')*B^a_{l2,m2}(k+G)<bl,l2,l1|l2> 
      complex(8) :: cc  !  C^{a*}_{l1,m1}(k-q+K')*C^a_{l2,m2}(k+G)<bl,l2,l1|l2> 

      logical :: lprt=.false.
      logical :: notf

      integer,allocatable:: jipw(:,:)    ! Counter: Runs over IPW's
      complex(8),allocatable:: tmat(:,:),tmat2(:,:),mnn(:,:),mmat(:,:,:)

      character(8), parameter :: sname = 'calcminm'
!
! !EXTERNAL ROUTINES: 
      real(8), external :: getcgcoef


! !INTRINSIC ROUTINES: 


      intrinsic abs
      intrinsic cpu_time
      intrinsic min
      intrinsic nint

     
! !REVISION HISTORY:
! 
! Created  23th. Feb. 2004 by RGA
! Last modified Jan. 16, 2007 by JH

! !TO DO:
!
! - Check complex and real cases
!
!EOP
!
!BOC

!
! Expand eigenvectors 
!
      if(lprt) call linmsg(6,'-','calcminm') 
      call cpu_time(tstart)

      ndim=nend - nstart + 1
      mdim=mend - mstart + 1
      nmdim=ndim*mdim

      if(lprt) call linmsg(6,'-','calcminm')
      allocate(mmat(mbsiz,mstart:mend,nstart:nend),stat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to allocate mmat")
      mmat=czero

      jk=kqid(ik,iq)
      do i=1,3
        qvec(i)=dble(qlist(i,iq))/dble(idvq)
        x=dble(klist(i,ik)-klist(i,jk))/dble(idvk)-qvec(i)
        igc(i)=nint(x)
      enddo  

      idf=0
      imix = 0
      ! TODO Any documentation to this large loop?
      do iat = 1, nat           !! Loop over atoms
        do ieq = 1, mult(iat)   !! Loop over equivalent atoms
          idf = idf + 1
          arg=pos(1,idf)*qvec(1)+pos(2,idf)*qvec(2)+pos(3,idf)*qvec(3)
          phs=cmplx(cos(2.0d0*pi*arg),-sin(2.0d0*pi*arg),8)

          do irm = 1, nmix(iat)   !!  Loop over mixed functions:
            bl=bigl(irm,iat)
              
            do bm=-bl,bl
              imix = imix + 1

              suml12=czero
              do l1=0,lmax           !! loop over l1:
                l2min=abs(bl-l1)
                l2max=min(bl+l1,lmax)
                la1 = l1 + ncore(iat) + 1
                lb1 = la1 + lmax+1

                do l2=l2min,l2max           !! loop over l2:
                  la2 = l2 + ncore(iat) + 1 
                  lb2 = la2 + (lmax+1)

                  do m1=-l1,l1                !! loop over m1:
                    lm1 = l1*l1+l1+m1+1
                    m2= -bm + m1

                    if(abs(m2).gt.l2) cycle 
                    lm2 = l2*l2+l2+m2+1
                    angint = getcgcoef(l2,bl,l1,m2,bm)        !! Calculate the angular integral:

                    if(abs(angint).lt.1.0d-8) cycle 

                    do ie1=nstart,nend                   !! Loop over eigenfunctions at k
                      aa = s3r(irm,la1,la2,iat,isp)*alfa(ie1,lm1,idf)
                      ba = s3r(irm,lb1,la2,iat,isp)*beta(ie1,lm1,idf)

                      ab = s3r(irm,la1,lb2,iat,isp)*alfa(ie1,lm1,idf)
                      bb = s3r(irm,lb1,lb2,iat,isp)*beta(ie1,lm1,idf)
                    
                      na = aa + ba
                      nb = ab + bb

                      do ilo1=1,nlo_at(l1,iat) 
                        lc1 = lo_index_s3r(ilo1,l1,iat)       
                        na = na + s3r(irm,lc1,la2,iat,isp)*gama(ie1,ilo1,lm1,idf)
                        nb = nb + s3r(irm,lc1,lb2,iat,isp)*gama(ie1,ilo1,lm1,idf)
                      enddo 

                      do ie2=mstart,mend
                        sm = na*alfp(ie2,lm2,idf)+nb*betp(ie2,lm2,idf)

                        do ilo2=1,nLO_at(l2,iat)
                          lc2 = lo_index_s3r(ilo2,l2,iat)
                          
                          ac=s3r(irm,la1,lc2,iat,isp)*alfa(ie1,lm1,idf)
                          bc=s3r(irm,lb1,lc2,iat,isp)*beta(ie1,lm1,idf)
                          
                          nc = ac + bc 
                          do ilo1=1,nLO_at(l1,iat) 
                            lc1 = lo_index_s3r(ilo1,l1,iat)       
                            nc=nc+s3r(irm,lc1,lc2,iat,isp)*gama(ie1,ilo1,lm1,idf)
                          enddo     
                          sm = sm + nc*gamp(ie2,ilo2,lm2,idf)
                        enddo

                        mmat(imix,ie2,ie1)=mmat(imix,ie2,ie1) &
     &                   + phs*angint*sm
                      enddo ! ie2
                    enddo !ie1

                  enddo ! m1     
                enddo ! l2
              enddo ! l1
            enddo ! bm
          enddo ! irm
        enddo ! ieq
      enddo ! iat

!
!     --------------------
!     Interstitial region:
!     --------------------
!     iv1: Loop over G:
!     iv2: loop over G'
!
      sqvi=sqrt(vi)
      maxig=size(indggq)
      nv1 = nv(kpirind(ik))
      nv2 = nv(kpirind(jk))
      allocate(jipw(1:nv1,1:nv2),                                       &
     &         tmat(1:nv2,1:nv1),                                       &
     &         tmat2(1:nv2,ndim),                                       &
     &         mnn(mdim,ndim),                                          &
     &         stat=ierr) 
      call errmsg0(ierr,sname,"Fail to allocate work arrays")

      !! setup jipw 
      jipw(:,:)=0
      do iv2=1, nv2 
        do iv1=1, nv1
          ikvec(1:3)= gindex(:,indgk(iv1,ik))-gindex(:,indgk(iv2,jk))   &
     &            +igc(1:3)
          if((minval(ikvec).ge.ngmin).and.(maxval(ikvec).le.ngmax)) then
            jpw=ig0(ikvec(1),ikvec(2),ikvec(3))
            if(jpw.gt.maxig)then
              jipw(iv1,iv2)=0    
!              write(6,'(a,4i6)' )'Warning: jipw=0 in calcminm',iv1,iv2, &
!     &                           jpw,ngqbarc(iq)
            else
              jipw(iv1,iv2)=indggq(jpw)
            endif    
          endif  
        enddo ! iv2
      enddo ! iv1

      do iipw = 1, ngq(iq)
        !! setup tmat = mpwipw
        do iv1=1, nv1 
          do iv2=1, nv2 
            if(jipw(iv1,iv2).gt.0)then
              tmat(iv2,iv1)=mpwipw(iipw,jipw(iv1,iv2))
            else 
              tmat(iv2,iv1)=czero
            endif    
          enddo ! iv2
        enddo ! iv1  

        call cpu_time(time1)
        call zgemm('n','n',nv2,ndim,nv1,cone,tmat,nv2,  &
     &           zzk(:,nstart:nend),maxngk,czero,tmat2,nv2 )
        call zgemm('t','n',mdim,ndim,nv2,cone,zzq(:,mstart:mend), &
     &           maxngk,tmat2,nv2,czero,mnn,mdim)
        call cpu_time(time2)
        time_lapack=time_lapack+time2-time1

        do ie1 = nstart, nend 
          do ie2= mstart,mend  
            mmat(iipw+locmatsiz,ie2,ie1)=mnn(ie2-mstart+1,ie1-nstart+1)*sqvi
          enddo 
        enddo 

      enddo ! iipw
      deallocate(jipw,tmat,tmat2,mnn)

!
! transform to the eigenvectors of the coulomb matrix.
! 
      call cpu_time(time1)
      call zgemm('c','n',matsiz,nmdim,mbsiz,cone,barcvm,mbsiz,mmat, &
     &           mbsiz,czero,minm,matsiz)
      call cpu_time(time2)
      time_lapack=time_lapack+time2-time1
      deallocate(mmat) 

      if(iop_coul.gt.-1.and.iq.eq.1) then 
        do ie1=nstart,nend 
          do ie2=mstart,mend
            if(ie1.eq.ie2) then 
              minm(im_g0,ie2,ie1) = sqrt(vi)*barcevsq(im_g0) 
            else
              minm(im_g0,ie2,ie1) = 0.d0
            endif 
          enddo
        enddo
      endif

      call cpu_time(tend)
      time_minm=time_minm+tend-tstart

      end subroutine calcminm
!EOC      


