!BOP
! 
! !ROUTINE: mb_calcs3r
!
! !INTERFACE:
      subroutine mb_calcs3r(iat,isp)

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
!
! !USES:

      use struk,     only: atomname,nrpt
      use constants, only: cfein1, cfein2
      use core,      only: ncore,lcore,ucore,uscore
      use lapwlo,    only: lmax,lomax,nlo_at,nt
      use mixbasis,  only: nmix,umix,usmix,bigl,s3r,nrwmax
      use radwf,     only: a,u,ulo,udot,b,us,uslo,usdot 
      use task,      only: fid_outmb
      
! !INPUT PARAMETERS:
    
      implicit none
      
      integer, intent(in) :: iat ! index for inequivalent atoms
      integer, intent(in) :: isp ! index for spin 

! !LOCAL VARIABLES:
!

      integer :: l1         !  l for atomic function iul (left in the product)
      integer :: l2         !  l for atomic function iur (right in the product)
      integer :: lb         !  bigl, l for the radial mixed basis function 
      integer :: irm        ! Indexes the radial mixed function
      integer :: ir,ir1,ir2 ! Indexes the radial function
      integer :: it1,it2
      integer :: ilo        ! index for LO 
      integer :: idot,jdot  ! 0 for u, 1 for udot
      integer :: irp        ! Indexes the radial mesh point.
      integer :: npt        ! number of radial mesh points,locally stores nrpt(iat)
      integer :: nrw        ! number of radial wavefunctions 

      real(8), allocatable :: a1(:) ! r u_l(r,E) at mesh points.
      real(8), allocatable :: b1(:) ! r us_l(r,E) at mesh points.
      real(8), allocatable :: a2(:) ! r u_l(r,E) at mesh points.
      real(8), allocatable :: b2(:) ! r us_l(r,E) at mesh points.
      real(8), allocatable :: a3(:) ! r u_l(r,E) at mesh points.
      real(8), allocatable :: b3(:) ! r us_l(r,E) at mesh points.
      real(8), allocatable :: rp(:) ! Spherical coordinate of the mesh point irp.
      real(8), allocatable :: rwf_all(:,:,:) ! radial wavefunctions
      integer, allocatable :: lrw_all(:)     ! angular momentum of radial wavefunctions 
      real(8) :: rint               ! value of the integral.
      character*5::orb_info(4) 
 
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
!
!     Initializations     
! 

      orb_info(1) = ' core'
      orb_info(2) = '     '
      orb_info(3) = ' dot '
      orb_info(4) = ' lo  ' 

      npt = nrpt(iat)
      allocate(a1(npt),a2(npt),a3(npt))
      allocate(b1(npt),b2(npt),b3(npt))
      allocate(rp(npt))
      allocate(rwf_all(npt,2,nrwmax),lrw_all(nrwmax))

      !! Calculate the mesh of radial points
      call radmesh(iat,rp)

      !* first collect all radial wavefunctions 
      do ir=1,ncore(iat)  ! core states
        lrw_all(ir) = lcore(ir,iat)
        rwf_all(1:npt,1,ir) = ucore(1:npt,ir,iat,isp)
        rwf_all(1:npt,2,ir) = uscore(1:npt,ir,iat,isp)
      enddo
      ir = ncore(iat) 

      do idot=0,1
        do l1=0,lmax
          ir = ir + 1 
          lrw_all(ir) = l1
          if(idot.eq.0) then 
            rwf_all(1:npt,1,ir)=u(1:npt,l1,iat,isp)
            rwf_all(1:npt,2,ir)=us(1:npt,l1,iat,isp)
          else
            rwf_all(1:npt,1,ir)=udot(1:npt,l1,iat,isp)
            rwf_all(1:npt,2,ir)=usdot(1:npt,l1,iat,isp)
          endif 
        enddo
      enddo
      do l1=0,lomax
        do ilo=1,nLO_at(l1,iat)
          ir = ir + 1
          lrw_all(ir) = l1
          rwf_all(1:npt,1,ir) = ulo(1:npt,ilo,l1,iat,isp)
          rwf_all(1:npt,2,ir) = uslo(1:npt,ilo,l1,iat,isp)
        enddo
      enddo 

      nrw = ir

      s3r(:,:,:,iat,isp)=0.0d0

      write(fid_outmb,100) atomname(iat)
      do irm=1,nmix(iat)

        a(1:npt)=umix(1:npt,irm,iat)
        b(1:npt)=usmix(1:npt,irm,iat)
        lb = bigl(irm,iat)

        do ir2=1,nrw

          l2 = lrw_all(ir2) 
          a2(1:npt)=rwf_all(1:npt,1,ir2)
          b2(1:npt)=rwf_all(1:npt,2,ir2)

          if(ir2.le.ncore(iat)) then  
            it2 = 1    !! core 
          else if (ir2.gt.ncore(iat).and.ir2.le.ncore(iat)+nt) then  
            it2 = 2    !! u
          else if (ir2.gt.ncore(iat)+nt.and.ir2.le.ncore(iat)+nt*2) then
            it2 = 3    !! udot
          else 
            it2 = 4    !! ulo
          endif 
        
          do ir1=1,ir2
            l1=lrw_all(ir1) 

            if(ir1.le.ncore(iat)) then  
              it1 = 1    !! core 
            else if (ir1.gt.ncore(iat).and.ir1.le.ncore(iat)+nt) then  
              it1 = 2    !! u
            else if (ir1.gt.ncore(iat)+nt.and.ir1.le.ncore(iat)+nt*2) then
              it1 = 3    !! udot
            else
              it1 = 4    !! ulo
            endif
 
            if((iabs(l1-l2).le.lb).and.(l1+l2.ge.lb))then
              a1(1:npt)=rwf_all(1:npt,1,ir1)
              b1(1:npt)=rwf_all(1:npt,2,ir1)

              do irp=1,npt
                a3(irp)=a1(irp)*a2(irp)/rp(irp)
                b3(irp)=b1(irp)*b2(irp)/rp(irp)
              enddo 

              call rint13(cfein1,cfein2,a,b,a3,b3,rint,iat)
              s3r(irm,ir1,ir2,iat,isp)=rint
              s3r(irm,ir2,ir1,iat,isp)=rint

              write(fid_outmb,101) irm,lb,l1,orb_info(it2),l2,orb_info(it1),rint

            endif  
          enddo ! ir1
        enddo  ! ir2
      enddo ! irm

      deallocate(rp,a1,a2,a3,b1,b2,b3,rwf_all,lrw_all) 
!   10 format(4i4,d16.9)
  100 format(/,5x,'Integrals <v_(NL)u_(l1)|u_(l2)> for atom ',a10,/,    &
     &       13x,'N',3x,'L',2x,'l1',1x,'u_',4x,'l2',1x,'u_',8x,     &
     &       '<v u|u>')
  101 format(10x,3i4,a5,i4,a5,1pe19.11)

   11 format(10x,3i4,' core',i4,' core',1pe19.11)
   12 format(10x,3i4,' core',i4,'     ',1pe19.11)
   13 format(10x,3i4,' core',i4,' dot ',1pe19.11)
   14 format(10x,3i4,' core',i4,' lo  ',1pe19.11)
   21 format(10x,3i4,'     ',i4,' core',1pe19.11)
   22 format(10x,3i4,'     ',i4,'     ',1pe19.11)
   23 format(10x,3i4,'     ',i4,' dot ',1pe19.11)
   24 format(10x,3i4,'     ',i4,' lo  ',1pe19.11)
   31 format(10x,3i4,' dot ',i4,' core',1pe19.11)
   32 format(10x,3i4,' dot ',i4,'     ',1pe19.11)
   33 format(10x,3i4,' dot ',i4,' dot ',1pe19.11)
   34 format(10x,3i4,' dot ',i4,' lo  ',1pe19.11)
   41 format(10x,3i4,' lo  ',i4,' core',1pe19.11)
   42 format(10x,3i4,' lo  ',i4,'     ',1pe19.11)
   43 format(10x,3i4,' lo  ',i4,' dot ',1pe19.11)
   44 format(10x,3i4,' lo  ',i4,' lo  ',1pe19.11)
      end subroutine mb_calcs3r
!EOC

