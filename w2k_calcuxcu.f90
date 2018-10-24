!BOP
! 
! !ROUTINE: w2k_calcuxcu
!
! !INTERFACE:
      subroutine w2k_calcuxcu(iat,isp)

! !DESCRIPTION:
!
! This subroutine calculates the set of integrals:
!
! \begin{equation}
!\langle LM\tilde{\lambda} | \tilde{\lambda'}\rangle_a\equiv%
!\int^{R^a_{MT}}_{0}{\upsilon_{aLM}(r^a)%
!\tilde{u}^*_{\lambda}(r^a,E_{\lambda})\tilde{u}_{\lambda'}(r^a,E_{\lambda'})\left(r^a\right)^2dr^a}
!\end{equation}
!
! where $\tilde{u}$ can be $u$, $\dot{u}$ or $u^{lo}$
! when lvorb is true, uiorb(:,:,:,:) are also calculated  
!
! !USES:

      use constants, only: cfein1, cfein2, pi
      use core,      only: ncore,ucore,uscore
      use lapwlo,    only: lmax,lomax,nlo_at 
      use xcpot,     only: lxcm, lmxc, vxclm, uxcu, &
     &                     natorb,iatorb,nlorb,lorb,uiorb,lvorb
      use radwf,     only: a,u,ulo,udot,b,us,uslo,usdot
      use struk,     only: atomname,nrpt
      use task,      only: fid_outdbg
      
! !INPUT PARAMETERS:
    
      implicit none
      
      integer, intent(in) :: iat 
      integer, intent(in) :: isp


! !LOCAL VARIABLES:
!

      integer :: bl
      integer :: l,l1  ! Angular momentum quantum number of the atomic function iul
      integer :: l12 ! Index for packed storage of (l1,l2)
      integer :: l2  ! Angular momentum quantum number of the atomic function iur
      integer :: lxc  ! Indexes the radial mixed function
      integer :: lim1,lim2
      integer :: ir1,ir2 ! 1 for u, 2 for udot, 3 for ulo
      integer :: ir12   ! Index for packed storage of (ir,jr)
      integer :: irp     ! Indexes the radial mesh point.
      integer :: npt     ! number of radial mesh points, locally stores nrpt(iat) a
      integer :: ia,il
      integer :: nr1, nr2  !! number of radial functions for a particular l
      integer :: ilo,jlo,ij
      logical:: lvorb_iat
!
      real(8), allocatable :: a1(:) ! r u_l(r,E) at mesh points.
      real(8), allocatable :: b1(:) ! r us_l(r,E) at mesh points.
      real(8), allocatable :: a2(:) ! r u_l(r,E) at mesh points.
      real(8), allocatable :: b2(:) ! r us_l(r,E) at mesh points.
      real(8), allocatable :: a3(:) ! r u_l(r,E) at mesh points.
      real(8), allocatable :: b3(:) ! r us_l(r,E) at mesh points.
      real(8), allocatable :: a4(:) ! r u_l(r,E) at mesh points.
      real(8), allocatable :: b4(:) ! r us_l(r,E) at mesh points.
      real(8) :: rint               ! value of the integral.
      real(8), allocatable :: rp(:) ! Spherical coordinate of the mesh point irp.

 
!
! !EXTERNAL ROUTINES: 
!


      external rint13
      external radmesh

 
! 
! !INTRINSIC ROUTINES: 
!


      intrinsic dexp
      intrinsic iabs
      intrinsic isign

! 
! !REVISION HISTORY:
! 
! Created Dic 2003
! Last Modified: Feb. 2nd. 2004
!
!EOP
!BOC

     
      npt = nrpt(iat)
      allocate(a1(npt),a2(npt),a3(npt),b1(npt),b2(npt),b3(npt),a4(npt),  &
     &         b4(npt),rp(npt))

      call radmesh(iat,rp)

      do lxc=1,lxcm(iat)
        bl=abs(lmxc(1,lxc,iat))

        do l2=0,lmax   ! right radial function 

          if(l2.gt.lomax) then 
            nr2 = 2 
          else 
            nr2 = nLO_at(l2,iat) + 2
          endif 

          do l1=0,lmax

            l12 = l1+1 + l2*(lmax+1) 

            if(.not.( abs(bl-l2).le.l1  .and. l1.le.(bl+l2) ) .and.   &
     &         .not.( iabs(bl-l1).le.l2 .and. l2.le.(bl+l1) ) ) cycle 

            if(l1.gt.lomax) then 
              nr1 = 2 
            else 
              nr1 = nLO_at(l1,iat) + 2
            endif

            ir12 = 0 

            do ir2 = 1, nr2 

              if(ir2.eq.1) then 
                a2(:) = u( 1:npt,l2,iat,isp)
                b2(:) = us(1:npt,l2,iat,isp)
              elseif(ir2.eq.2) then 
                a2(:)=udot( 1:npt,l2,iat,isp)
                b2(:)=usdot(1:npt,l2,iat,isp)
              else
                a2(:)= ulo( 1:npt,ir2-2,l2,iat,isp)
                b2(:)= uslo(1:npt,ir2-2,l2,iat,isp)
              endif 

              do irp=1,npt
                a3(irp)=a2(irp)*vxclm(irp,lxc,iat,isp)/rp(irp)
                b3(irp)=b2(irp)*vxclm(irp,lxc,iat,isp)/rp(irp)
              enddo
        
              do ir1 = 1, nr1
     
                if(ir1.eq.1) then 
                  a1(:)=u( 1:npt,l1,iat,isp)
                  b1(:)=us(1:npt,l1,iat,isp)
                elseif(ir1.eq.2) then 
                  a1(:)=udot( 1:npt,l1,iat,isp)
                  b1(:)=usdot(1:npt,l1,iat,isp)
                else
                  a1(:)=ulo( 1:npt,ir1-2,l1,iat,isp)
                  b1(:)=uslo(1:npt,ir1-2,l1,iat,isp)
                endif 

                a1 = a1/rp
                b1 = b1/rp 
                    
                call rint13(cfein1,cfein2,a1,b1,a3,b3,rint,iat)

                ir12 = (ir2-1)*nr1 + ir1
                uxcu(ir12,l12,lxc,iat,isp) = rint
!                ir12 = (ir1-1)*nr2 + ir2
!                uxcu(ir12,l12,lxc,iat,isp) = rint
              enddo ! jr  
            enddo ! l2
          enddo  ! ir
        enddo ! l1  
      enddo ! lxc

!
! Calculate integrals needed for LDA+U based GW calculations 
!
      lvorb_iat=.false.
      if(lvorb) then 
        do ia=1,natorb
          if( iatorb(ia) .eq. iat )  then
            lvorb_iat=.true.
            exit
          endif
        enddo
      endif 
 
      if(lvorb_iat) then 

#ifdef DEBUG        
        write(fid_outdbg,*)
        write(fid_outdbg,*) '#uiorb for atom ',iat
        write(fid_outdbg,*)
#endif

        do il=1,nlorb(ia)  
          l=lorb(il,ia) 
          a1(1:npt) = u(1:npt,l,iat,isp)
          b1(1:npt) = us(1:npt,l,iat,isp) 
          a2(1:npt) = udot(1:npt,l,iat,isp)
          b2(1:npt) = usdot(1:npt,l,iat,isp)

! <u_l u_l>
          call rint13(cfein1,cfein2,a1,b1,a1,b1,rint,iat)
          uiorb(1,il,ia,isp) = rint 

! <u_l udot_l > 

          call rint13(cfein1,cfein2,a1,b1,a2,b2,rint,iat)
          uiorb(2,il,ia,isp) = rint

! <udot_l udot_l >

          call rint13(cfein1,cfein2,a2,b2,a2,b2,rint,iat)
          uiorb(3,il,ia,isp) = rint

          ij = 3
          do jlo=1,nlo_at(l,iat) 
            a3(1:npt) = ulo(1:npt,jlo,l,iat,isp)
            b3(1:npt) = uslo(1:npt,jlo,l,iat,isp)

! <u_l ulo_l> 
            ij = ij + 1
            call rint13(cfein1,cfein2,a1,b1,a3,b3,rint,iat)
            uiorb(ij,il,ia,isp) = rint

! <udot_l ulo_l>
            ij = ij + 1
            call rint13(cfein1,cfein2,a2,b2,a3,b3,rint,iat)
            uiorb(ij,il,ia,isp) = rint

            do ilo=1,jlo 
              a4(1:npt) = ulo(1:npt,ilo,l,iat,isp)
              b4(1:npt) = uslo(1:npt,ilo,l,iat,isp)

! <u2_l u2_l>
              ij = ij + 1
              call rint13(cfein1,cfein2,a4,b4,a3,b3,rint,iat)
              uiorb(ij,il,ia,isp) = rint
            enddo
          enddo 

#ifdef DEBUG        
          write(fid_outdbg,'(6e16.8)') uiorb(1:6,il,ia,isp)  
#endif 
        enddo     
      endif 

      deallocate(a1,a2,a3,a4,b1,b2,b3,b4,rp)
!   10 format(4i4,d16.9)
  100 format(/,5x,'Integrals <u_(l1)|vxc_(LM)|u_(l2)> for atom ',a10)
  200 format(4a5,6a18) 
  201 format(4i5,6e18.10) 

    2 format(10x,3i4,'     ',i4,' dot ',1pe19.11)
    3 format(10x,3i4,'     ',i4,' lo  ',1pe19.11)
    4 format(10x,3i4,' dot ',i4,' dot ',1pe19.11)
    5 format(10x,3i4,' dot ',i4,' lo  ',1pe19.11)
    6 format(10x,3i4,' lo  ',i4,' lo  ',1pe19.11)

      end subroutine w2k_calcuxcu
!EOC

