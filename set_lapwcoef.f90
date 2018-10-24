!BOP
!
! !ROUTINE: set_lapwcoef

! !INTERFACE:
      subroutine set_lapwcoef(ik,iop,ifrot,isp)
      
! !DESCRIPTION:
!
! This subroutine calculates the coefficients Alm, Blm and Clm of the
!augmentation for LAPW wave functions, according to the formulas:
!
!\begin{subequations}
!\begin{align}
!A^{\alpha}_{lm}(\vec{k}+\vec{G})=&\frac{4\pi}{\sqrt{\Omega}}R_{MT_{\alpha}}^2 %
!   i^l Y^*_{lm}(\hat{k+G})e^{i(\vec{k}+\vec{G})\cdot\vec{r}_{\alpha}}a_l \\    
!B^{\alpha}_{lm}(\vec{k}+\vec{G})=&\frac{4\pi}{\sqrt{\Omega}}R_{MT_{\alpha}}^2 %
!   i^l Y^*_{lm}(\hat{k+G})e^{i(\vec{k}+\vec{G})\cdot\vec{r}_{\alpha}}b_l \\    
!C^{\alpha}_{lm}(\vec{k}+\vec{G})=&\frac{4\pi}{\sqrt{\Omega}}R_{MT_{\alpha}}^2 %
!   i^l Y^*_{lm}(\hat{k+G})e^{i(\vec{k}+\vec{G})\cdot\vec{r}_{\alpha}}c_l     
!\end{align}
!\end{subequations}
!
!
! !USES:

      use bands,     only: nv
      use constants, only: pi, imag,czero
      use eigenvec,  only: alm,blm,clm
      use kpoints,   only: idvk, klist, kpirind
      use lapwlo,    only: nt,lapw,nlo_tot,abcelo,umt,lmax,lomax,nlo_at,nlomax
      use struk,     only: mult, pos, rmt, vi, nat,rotij,rotloc
      use recipvec,  only: gindex, indgk, rk
      use task,      only: casename
      

! !INPUT PARAMETERS:

      implicit none

      integer(4), intent(in) :: ik    ! index of the k-point
      integer(4), intent(in) :: iop   !  =1 --- alfa, 2 --- alfap
      logical,    intent(in) :: ifrot ! true if local rotations are taken  into account 
      integer(4), intent(in) :: isp   ! index for spin 
      
! !LOCAL VARIABLES:
      
      real(8), allocatable :: al(:,:,:)   ! coefficient al(kn) for atom i
      real(8), allocatable :: bl(:,:,:)   ! coefficient bl(kn) for atom i
      complex(8), allocatable ::  yl(:,:) ! temporary storage of   spherical harmonics for local orbitals

      integer(4) ::   n 
      integer(4) :: iat 
      integer(4) :: l,m 
      integer(4) :: idf 
      integer(4) :: ieq 
      integer(4) :: irk
      integer(4) :: i 
      integer(4) :: ilm 
      integer(4) :: nvk,ngk
      integer(4) :: jlo,ilo0
      integer(4) :: i1 
      integer(4) :: i2 
      integer(4) :: igvec(3)
      
      real(8) :: rkn
      real(8) :: rmt2        !! square of RMT(iat)
      real(8) :: dfj(0:lmax) 
      real(8) :: fj(0:lmax)
      real(8) :: rhoatm 
      real(8) :: kvec(3)
      real(8) :: arg 
      real(8) :: u0, ud0, du0,dud0 
      real(8) :: vec(3) 
      real(8) :: rotv1(3) 
      real(8) :: rotv2(3)
      real(8) :: alo,blo,clo,comr,comi
      real(8) :: imsgn
      
      complex(8) :: cil
      complex(8), allocatable ::  phs(:)
      real(8), allocatable ::  kpg(:,:)

      logical :: ldbg=.false.

 
! !EXTERNAL ROUTINES: 

      
      external lohns
      external rotate
      external sphbes
      external ylm

!
! !REVISION HISTORY:
! 
! Last modified 1st. Dec 2004 by RGA
!
!EOP
!BOC

! set the sign of the imagniary part

      if(ldbg) write(6,*) "------------ start set_lapwcoef --------"
      if(iop.eq.1) then 
        imsgn=1.d0
      else 
        imsgn=-1.d0
      endif 

      irk = kpirind(ik)
      nvk = nv(irk)
      ngk = nvk+nlo_tot

      allocate(al(nvk,0:nt,nat))
      allocate(bl(nvk,0:nt,nat))
      allocate(yl((nt+1)*(nt+1),ngk))
      allocate(phs(ngk))
      allocate( kpg(1:3,1:ngk))
      al=0.d0
      bl=0.d0
      yl=0.d0
      phs=0.d0

      call sub_setkpg
!
!     precompute al(kn) and bl(kn)  
!
      if(ldbg) write(6,*) "- calculate al and bl -"
      do n = 1, nvk
        rkn = rk(n)
        do iat = 1, nat
          rmt2 = rmt(iat)**2
 
          call sphbes(nt-1,rmt(iat)*rkn,fj,dfj)
          do l = 0,lmax 
            u0   = umt(2,l,iat,isp)
            ud0  = umt(3,l,iat,isp)
            du0  = umt(4,l,iat,isp)
            dud0 = umt(5,l,iat,isp)

            if(lapw(l,iat)) then
              al(n,l,iat) = rkn*dfj(l)*ud0 - fj(l)*dud0 
              bl(n,l,iat) = fj(l)*du0 - rkn*dfj(l)*u0
            else
              al(n,l,iat) = fj(l)/u0/rmt2
              bl(n,l,iat) = 0.d0
            endif
          enddo
        enddo
      enddo

!
! Calculate all required spherical harmonics
!
      do i=1,3
        kvec(i)=dble(klist(i,ik))/dble(idvk)
      enddo 

      idf = 0
      do iat = 1, nat
        rhoatm = 4.0d+0*pi*sqrt(vi)*rmt(iat)*rmt(iat)
        do ieq = 1, mult(iat)
          idf = idf + 1

          if(ldbg) write(6,*) "iat,ieq=",iat,ieq   
          if(ldbg) write(6,*) "- calculate spherical harmonics -"

          do i = 1, ngk
            igvec(1:3)=gindex(:,indgk(i,ik))
            arg = pos(1,idf)*(dble(igvec(1)) + kvec(1))   &  ! x
     &          + pos(2,idf)*(dble(igvec(2)) + kvec(2))   &  ! y
     &          + pos(3,idf)*(dble(igvec(3)) + kvec(3))      ! z
            phs(i) = cmplx(dcos(2.0d+0*pi*arg),dsin(2.0d+0*pi*arg),8)

            vec(1:3) = kpg(1:3,i)
            if(ifrot)then
              call rotate (vec,rotij(1:3,1:3,idf),rotv1)
              call rotate (rotv1,rotloc(1:3,1:3,iat),rotv2)
              call ylm (rotv2,nt,yl(1,i))
            else  
              call ylm (vec,nt,yl(1,i))
            endif
          enddo ! i
!
! calcalate alm,blm,clm
!
          if(ldbg) write(6,*) "- calculate alm,blm and clm -"

          cil = (1.0d+0,0.0d+0)
          ilm = 0

          if(ldbg) then 
            write(99,*) "A_lm, B_lm, C_lm on the atom",idf
            write(99,100) "i","ilm","A_lm","B_lm","C_lm"
          endif 

          do l = 0,lmax
            do m = -l, l

              ilm = ilm + 1

              do i = 1, nvk
                comr=rhoatm*dble(cil*phs(i)*conjg(yl(ilm,i)))
                comi=rhoatm*aimag(cil*phs(i)*conjg(yl(ilm,i)))*imsgn
                alm(i,ilm,idf) = al(i,l,iat)*cmplx(comr,comi,8)
                blm(i,ilm,idf) = bl(i,l,iat)*cmplx(comr,comi,8)
              enddo ! i

              do i = nvk+1,nvk+nlo_tot
                alm(i,ilm,idf) = czero
                blm(i,ilm,idf) = czero
                if(nlomax.gt.0.and.l.le.lomax) clm(i,:,ilm,idf) = czero
              enddo ! i 

              if( l .gt. lomax) cycle 

              if(lapw(l,iat)) then  !! LAPW+LO
                ilo0 = 1
              else                  !! APW+lo+LO
                ilo0 = 0
              endif 

              do jlo= ilo0,nlo_at(l,iat)
                call lohns (iat,mult,i1,i2,l,jlo-ilo0+1)
                alo=abcelo(1,jlo,l,iat,isp)
                blo=abcelo(2,jlo,l,iat,isp)
                clo=abcelo(3,jlo,l,iat,isp)
                
                do i = nvk+i1,nvk+i2  
                  comr=rhoatm*dble(cil*phs(i)*conjg(yl(ilm,i)))
                  comi=rhoatm*aimag(cil*phs(i)*conjg(yl(ilm,i)))*imsgn
                  alm(i,ilm,idf)=alo*cmplx(comr,comi,8)
                  blm(i,ilm,idf)=blo*cmplx(comr,comi,8)

                  if(jlo.ge.1) then 
                    clm(i,jlo,ilm,idf)=clo*cmplx(comr,comi,8)
                  endif 
                  if(ldbg) then 
                    if(jlo.eq.0) then 
                      write(99,101) i,ilm,alm(i,ilm,idf),blm(i,ilm,idf)
                    else
                      write(99,101) i,ilm,alm(i,ilm,idf),blm(i,ilm,idf),&
     &                 clm(i,jlo,ilm,idf) 
                    endif
                  endif  
 
                enddo ! i
              enddo ! jlo
            enddo ! m

            cil = cil*imag
          enddo ! l

        enddo ! ieq 
      enddo !iat

 100  format(a6,a6,a24,a24,a24) 
 101  format(i6,i6,f12.6,f12.6,f12.6,f12.6,f12.6,f12.6,f12.6)

      deallocate(al,bl,yl,phs,kpg)
      return

      contains 

      subroutine sub_setkpg()
      implicit none
      
! !LOCAL VARIABLES:
      integer(4) :: ig,i
      integer(4), dimension(3) :: ikvec
      real(8), dimension(3) :: kvec
      real(8), dimension(3) :: gvec
      external k2cart 
!
!BOC
      ikvec(1:3)=klist(1:3,ik)
      call k2cart(ikvec,idvk,kvec)

      kpg(1:3,1:ngk)=0.0d0
      rk(1:ngk)=0.0d0
      do ig=1,ngk
        ikvec(1:3)=gindex(:,indgk(ig,ik))
        call k2cart(ikvec,1,gvec)
        kpg(:,ig) = kvec + gvec
        rk(ig) = sqrt( sum(kpg(1:3,ig)**2) )
      enddo
      
      end subroutine 

      end subroutine 
!EOC
