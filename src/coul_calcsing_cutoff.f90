!BOP
!
! !ROUTINE: coul_calcsing_cutoff
!
! !INTERFACE:
      subroutine coul_calcsing_cutoff
!
! !DESCRIPTION: 
!
!This subroutine calculates the matrix elements
! $< j_{\lambda}^{|\vec{G}+\vec{q}|}>_{aNL}$ according to equation \ref{jlam}
!
!
! !USES:
!
!
      use constants,  only: cfein1, cfein2, czero, pi
      use lapwlo,     only: kmax
      use mixbasis,   only: kmr,pwm,umix,usmix,nmix,nmixmax,bigl,maxbigl
      use q0barc,     only: nglen, sing, glen0, phase
      use struk,      only: pia,ndf,mult,nrpt,pos,br2,ortho,nat,rbas
      use barcoul,    only: axis_cut_coul,zcut_coul
!
!
! !INPUT PARAMETERS:
!
      implicit none

      


! !LOCAL VARIABLES:
!
      character*20:: sname="coul_calcsing_cutoff"
      integer(4) :: i
      integer(4) :: iat
      integer(4) :: idf
      integer(4) :: ieq    
      integer(4) :: igl
      integer(4) :: ig1
      integer(4) :: ig2
      integer(4) :: ig3
      integer(4) :: ipw
      integer(4) :: irm ! (Counter) runs over radial mixed basis functions
      integer(4) :: irp ! (Counter), runs over the radial mesh points.
      integer(4) :: jat
      integer(4) :: jdf
      integer(4) :: jeq
      integer(4) :: jgl
      integer(4) :: jri    ! Number of radial (logarithmic)  mesh points for atom iat
      integer(4) :: l   ! Angular momentum quantum number of the   mixed atomic function irm.
      integer(4) :: ng
      integer(4) :: ng1
      integer(4) :: ng2
      integer(4) :: ng3

      integer(4), dimension(3) :: igvec ! Integer coordinates of the G-vector
      integer(4), allocatable :: gind(:,:),gind4(:)
      
      real(8) :: rint  ! Radial integral of umix_{irm}j_{l}
      real(8) :: x     ! Argument of j_l, =|q+G|.r
      real(8) :: qglen ! Length of q+G
      real(8) :: gmax
      real(8) :: kk
      real(8) :: q
      real(8) :: glprev
      real(8) :: gpr
      real(8) :: gxy   ! projection of vector G on xOy plane
      real(8) :: gz    ! projection of vector G on z axis

      real(8), dimension(3) :: gvec ! Coordinates of the vector G
      real(8), dimension(3) :: dpos
      real(8), allocatable :: a1(:), b1(:) ! Temporary arrays for the radial functions
      real(8), allocatable :: sinf(:)    
      real(8), allocatable :: rp(:)     ! Radial mesh points
      real(8), allocatable :: glen(:)   ! Temporary reciprocal latt vecs
      
      complex(8) :: expg

      integer(4):: iprt,ierr 
 
 
! !EXTERNAL ROUTINES: 


      external k2cart
      external radmesh
      external rint13
      real(8),external :: vecprojlen

! !INTRINSIC ROUTINES: 

      intrinsic dsqrt

!
! !REVISION HISTORY:
!
!  Created: 17th. March 2004 by MF
! Modified: 30th. March 2004 by RGA
! Last Modified: 09th. Sept 2019 by MYZ
!
!EOP
!BOC
!      
      iprt = 0


      gmax=kmr*kmax*1.0d+1
      ng1=idint(gmax*pia(1))+1
      ng2=idint(gmax*pia(2))+1
      ng3=idint(gmax*pia(3))+1
      ng = (2*ng1+1)*(2*ng2+1)*(2*ng3+1)

      allocate(glen(ng),gind(3,ng),gind4(ng),stat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to allocate glen etc.")
      
!      write(*,13) maxgxc,maxgcoul,gmax
      ipw=0
      do ig1=-ng1,ng1
        igvec(1)=ig1
        do ig2=-ng2,ng2
          igvec(2)=ig2
          do ig3=-ng3,ng3
            igvec(3)=ig3
            gvec=matmul(br2,dble(igvec))
            kk=sqrt(sum(gvec(1:3)**2))
            if(kk.lt.gmax)then
              ipw = ipw + 1
              glen(ipw)=kk
              if(ortho)then
                do i=1,3
                  q=gvec(i)/pia(i)
                  gind(i,ipw)=nint(q+sign(1.0d-2,q))
                enddo
              else  
                gind(1:3,ipw) = igvec(1:3)
              endif  
            endif
          enddo
        enddo
      enddo
      ng = ipw
!     sort by increasing length using shell algorithm
      call shelsort(ng,gind(:,1:ng),glen(1:ng))

      glprev=-1.0d0
      igl=0
      do ipw=1,ng
        if(abs(glen(ipw)-glprev).gt.1.0d-10)then
          igl=igl+1
          glprev=glen(ipw)
        endif
        gind4(ipw)=igl
      enddo
      nglen=igl

      allocate(glen0(nglen))
      glen0=0.d0
      jgl=0
      do ipw=1,ng
        if(gind4(ipw).ne.jgl)then
          jgl=gind4(ipw)
          glen0(jgl)=glen(ipw)
        endif
      enddo    

      deallocate(glen)
!     
      allocate(phase(ndf,ndf,nglen),stat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to allocate phase")

      idf=0
      phase=czero
      do iat = 1, nat
        do ieq = 1, mult(iat)
          idf=idf+1
          jdf=0    
          do jat=1,nat
            do jeq=1,mult(jat)
              jdf=jdf+1
              dpos(1:3)= pos(1:3,idf)-pos(1:3,jdf)
              do ipw = 2, ng
                igl=gind4(ipw)
                igvec(1:3)=gind(1:3,ipw)
                gpr=dble(igvec(1))*dpos(1)+dble(igvec(2))*dpos(2)+ &
     &              dble(igvec(3)*dpos(3))
                call k2cart(igvec,1,gvec)
                expg=cmplx(cos(2.0d0*pi*gpr),-sin(2.0d0*pi*gpr),8)
                gxy=vecprojlen(gvec,rbas(axis_cut_coul,:),'perp')
                gz=vecprojlen(gvec,rbas(axis_cut_coul,:),'para')
                phase(idf,jdf,igl)=phase(idf,jdf,igl)+expg * &
                  cmplx(1.0d0-exp(-gxy*zcut_coul)*cos(gz*zcut_coul),0.0d0,8)
              enddo ! ipw
            enddo ! jeq
          enddo ! jat
        enddo ! ieq
      enddo ! iat

      deallocate(gind,gind4)
      allocate(sing(nat,nmixmax,nglen),stat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to allocate sing")
 
      sing=0.0d0
      do iat=1,nat
        jri=nrpt(iat)
        allocate(rp(jri))          
        allocate(a1(jri))            
        allocate(b1(jri))
        allocate(sinf(jri))          
        call radmesh(iat,rp)
        do igl=2,nglen
          qglen=glen0(igl)
!
!         Calculate the spherical bessel function at each mesh point
!
          do irp = 1, jri
            x=rp(irp)*qglen
            sinf(irp) = sin(x)
          enddo   ! irp
          do irm=1,nmix(iat)
            l=bigl(irm,iat)
            if(l.eq.0)then
              a1(1:jri)= umix(1:jri,irm,iat)
              b1(1:jri)=usmix(1:jri,irm,iat)
!             Integrate the wavefunctions:
              call rint13(cfein1,cfein2,sinf,sinf,a1,b1,rint,iat)
              sing(iat,irm,igl) = rint
            endif  
          enddo ! irm
        enddo ! igl
        deallocate(rp)
        deallocate(a1)
        deallocate(b1)
        deallocate(sinf)
      enddo !iat  
      
      end subroutine coul_calcsing_cutoff
!EOC
