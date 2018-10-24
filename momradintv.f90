!BOP
!
! !ROUTINE: momradintv
!
! !INTERFACE:
      subroutine momradintv(iat,isp)
      
! !DESCRIPTION:
!      
! This subroutine calculates the set of integrals of radial functions
! needed for the calculation of the momentum matrix elements (Valence
!states only). 
! 
! The naming convention for the corresponding integrals is:
! 
! \begin{subequations}\label{momradintv1}
!\begin{align}
! \texttt{iul1ul(iat,l)}  = &\int\limits_0^{R^a_{MT}}{ ( u_{l+1}(r) u'_l(r)r^2dr - l u_{l+1}(r)u_l(r)/r ) r^2 dr }
! \texttt{iul1udl(iat,l)} = &\int\limits_0^{R^a_{MT}}{ ( u_{l+1}(r)\dot{u}'_l(r)r^2dr}
!     - l \int\limits_0^{R^a_{MT}}{u_{l+1}(r)\dot{u}_l(r)rdr}\\
! \texttt{iudl1ul(iat,l)}=&\int\limits_0^{R^a_{MT}}{ \dot{u}_{l+1}(r)u'_l(r)r^2dr}
!     - l \int\limits_0^{R^a_{MT}}{\dot{u}_{l+1}(r)u_l(r)rdr}\\
! \texttt{iudl1udl(iat,l)}=&\int\limits_0^{R^a_{MT}}{ \dot{u}_{l+1}(r)\dot{u}'_l(r)r^2dr} 
!     - l \int\limits_0^{R^a_{MT}}{\dot{u}_{l+1}(r)\dot{u}_l(r)rdr}\\
! \texttt{iul1ulol(iat,l)}=&\int\limits_0^{R^a_{MT}}{ u_{l+1}(r)u'^{lo}_l(r)r^2dr}
!     - l \int\limits_0^{R^a_{MT}}{u_{l+1}(r)u^{lo}_l(r)rdr}\\
! \texttt{iulol1ul(iat,l)}=&\int\limits_0^{R^a_{MT}}{ u^{lo}_{l+1}(r) u'_l(r)r^2dr}
!     - l \int\limits_0^{R^a_{MT}}{u^{lo}_{l+1}(r)u_l(r)rdr}\\
! \texttt{iulol1udl(iat,l)}=&\int\limits_0^{R^a_{MT}}{ u^{lo}_{l+1}(r)\dot{u}'_l(r)r^2dr}
!     - l \int\limits_0^{R^a_{MT}}{u^{lo}_{l+1}(r)\dot{u}_l(r)rdr}\\
! \texttt{iudl1ulol(iat,l)}=&\int\limits_0^{R^a_{MT}}{ \dot{u}_{l+1}(r)u'^{lo}_l(r)r^2dr}
!     - l \int\limits_0^{R^a_{MT}}{\dot{u}_{l+1}(r)u^{lo}_l(r)rdr}\\
! \texttt{iulol1ulol(iat,l)} = & \int\limits_0^{R^a_{MT}}{ u^{lo}_{l+1}(r)u'^{lo}_l(r)r^2dr}
!     - l \int\limits_0^{R^a_{MT}}{u^{lo}_{l+1}(r)u^{lo}_l(r)rdr}
! \end{align}
! \end{subequations}
!
! \begin{subequations}\label{momradintv2}
!\begin{align}
! \texttt{iulul1(iat,l)}= &\int\limits_0^{R^a_{MT}}{ ( u_{l}(r)u'_{l+1}(r) + (l+2)/r u_{l}(r)u_{l+1}(r)) r^2dr}\\
! \texttt{iuludl1(iat,l)}=&\int\limits_0^{R^a_{MT}}{ ( u_{l}(r)\dot{u}'_{l+1}(r) + (l+2)/r u_{l}(r)\dot{u}_{l+1}(r) ) r^2 dr}\\
! \texttt{iudlul1(iat,l)}=&\int\limits_0^{R^a_{MT}}{ ( \dot{u}_{l}(r)u'_{l+1}(r) + (l+2)/r \dot{u}_{l}(r) u_{l+1}(r)) r^2 dr}\\
! \texttt{iudludl1(iat,l)}=&\int\limits_0^{R^a_{MT}}{\dot{u}_{l}(r)\dot{u}'_{l+1}(r)%
!r^2dr}+(l+2)\int\limits_0^{R^a_{MT}}{\dot{u}_{l}(r)\dot{u}_{l+1}(r)rdr}\\
! \texttt{iululol1(iat,l)}=&\int\limits_0^{R^a_{MT}}{u_{l}(r)u'^{lo}_{l+1}(r)r^2dr}+
! (l+2)\int\limits_0^{R^a_{MT}}{u_{l}(r)u^{lo}_{l+1}(r)rdr}\\
! \texttt{iulolul1(iat,l)}=&\int\limits_0^{R^a_{MT}}{u^{lo}_{l}(r)u'_{l+1}(r)r^2dr}+
! (l+2)\int\limits_0^{R^a_{MT}}{u^{lo}_{l}(r)u_{l+1}(r)rdr}\\
! \texttt{iuloludl1(iat,l)}=&\int\limits_0^{R^a_{MT}}{u^{lo}_{l}(r)\dot{u}'_{l+1}(r)r^2dr}+
! (l+2)\int\limits_0^{R^a_{MT}}{u^{lo}_{l}(r)\dot{u}_{l+1}(r)rdr}\\
! \texttt{iudlulol1(iat,l)}=&\int\limits_0^{R^a_{MT}}{\dot{u}_{l}(r)u'^{lo}_{l+1}(r)r^2dr}+
! (l+2)\int\limits_0^{R^a_{MT}}{\dot{u}_{l}(r)u^{lo}_{l+1}(r)rdr}\\
!\texttt{iulolulol1(iat,l)}=&\int\limits_0^{R^a_{MT}}{u^{lo}_{l}(r)u'^{lo}_{l+1}(r)r^2dr}+
! (l+2)\int\limits_0^{R^a_{MT}}{u^{lo}_{l}(r)u^{lo}_{l+1}(r)rdr}
! \end{align}
! \end{subequations}

! !USES:

      use lapwlo,  only: nt,lomax,nlo_at
      use mommat,  only: iul1ul, iul1udl, iudl1ul, iudl1udl, iulul1,     &
     &                   iuludl1, iudlul1, iudludl1, iul1ulol, iulol1ul, &
     &                   iulol1udl, iudl1ulol, iulol1ulol, iululol1,     &
     &                   iulolul1, iuloludl1, iudlulol1, iulolulol1
      use radwf,   only: u, ulo, udot, us, uslo, usdot
      use struk,   only: nat,atomname,nrpt
      use task,    only: fid_outmom

! !INPUT PARAMETERS:
      implicit none
      
      integer, intent(in) :: iat  ! Index for inequivalent atoms
      integer, intent(in) :: isp  ! Index for spin 
      
! !LOCAL VARIABLES:
  
      integer :: l    ! Angular momentum quantum number
      integer :: nr  ! number of radial mesh points
      integer :: ilo, jlo ! index for LO orbitals 
      real(8) :: fl      ! angular momentum quantum number
      real(8) :: iuup    ! \int(u_l du_l'/dr r^2 dr)
      real(8) :: iuu    ! \int(u_l u_l' r dr)
      
      real(8), allocatable :: rr(:) ! Coordinate of the radial mesh
      real(8), allocatable :: ul(:) ! u_l and u_{l+1}
      real(8), allocatable :: upl(:) ! radial derivative of the radial function du_l/dr
      real(8), allocatable :: uor(:) ! radial function over r: u_l/r
      real(8), allocatable :: udl(:) ! energy derivative of the radial function udot_l      
      real(8), allocatable :: udpl(:) ! radial derivative of udot: d udot_{l}/dr
      real(8), allocatable :: udor(:) ! energy derivative of the radial function over r : udot_l/r
      real(8), allocatable :: ulol(:) ! local orbital radial function ulo_l      
      real(8), allocatable :: ulopl(:) ! radial derivative of ulo: d ulo_{l}/dr
      real(8), allocatable :: uloor(:) ! local orbital the radial function over r : udot_l/r
            
! !DEFINED PARAMETERS:
 
      real(8), parameter :: two = 2.0d+0
      real(8), parameter :: one = 1.0d+0
      real(8), parameter :: cf2 = 1.0d-22
      logical:: ldbg = .false. 
!
! !EXTERNAL ROUTINES: 
      external rint13
      external radmesh
      
 
! !REVISION HISTORY:
!
! Created: 16. August 2004 by RGA
!
!EOP
!BOC
!
!     Allocate necesary arrays
!     
      nr=nrpt(iat)
      allocate(rr(nr))
      allocate(ul(nr))
      allocate(upl(nr))
      allocate(uor(nr))
      allocate(udl(nr))
      allocate(udpl(nr))
      allocate(udor(nr))
      allocate(ulol(nr))
      allocate(ulopl(nr))
      allocate(uloor(nr))
!
!     calculate the radial mesh points
!
      call radmesh(iat,rr)
!      
!     loop over l      
!      
      write(fid_outmom,*)
      write(fid_outmom,11)atomname(iat)

      do l=0,nt-2
        fl=dble(l)
!
!      store the needed shared arrays in the local ones
!        
        ul(1:nr)=u(1:nr,l+1,iat,isp)
        upl(1:nr)=two*(us(1:nr,l,iat,isp))
        uor(1:nr)=u(1:nr,l,iat,isp)/rr(1:nr)
        udl(1:nr)=udot(1:nr,l+1,iat,isp)
        udpl(1:nr)=two*(usdot(1:nr,l,iat,isp))
        udor(1:nr)=udot(1:nr,l,iat,isp)/rr(1:nr)
!        
!       iul1ul = int(u_{l+1} du_l/dr r^2 dr)-l*int(u_{l+1} u_l r dr)
!        
        call rint13(one,cf2,ul,ul,upl,upl,iuup,iat)
        call rint13(one,cf2,ul,ul,uor,uor,iuu,iat)
        iul1ul(l,iat,isp) = iuup - fl*iuu

!        
!       iul1udl = int(u_{l+1} d udot_l/dr r^2 dr)-l*int(u_{l+1} udot_l r dr)
!        
        call rint13(one,cf2,ul,ul,udpl,udpl,iuup,iat)
        call rint13(one,cf2,ul,ul,udor,udor,iuu,iat)
        iul1udl(l,iat,isp)=iuup-fl*iuu
!        
!       iudl1ul = int(udot_{l+1} du_l/dr r^2 dr)-l*int(udot_{l+1} u_l r dr)
!        
        call rint13(one,cf2,udl,udl,upl,upl,iuup,iat)
        call rint13(one,cf2,udl,udl,uor,uor,iuu,iat)
        iudl1ul(l,iat,isp)=iuup-fl*iuu
!        
!       iudl1udl = int(udot_{l+1} d udot_l/dr r^2 dr)-l*int(udot_{l+1} udot_l r dr)
        call rint13(one,cf2,udl,udl,udpl,udpl,iuup,iat)
        call rint13(one,cf2,udl,udl,udor,udor,iuu,iat)
        iudl1udl(l,iat,isp)=iuup-fl*iuu

        do ilo=1,nLO_at(l+1,iat)  
          ulol(1:nr)=ulo(1:nr,ilo,l+1,iat,isp)
!            
!         iulol1ul = int(ulo_{l+1} du_l/dr r^2 dr) - l*int(ulo_{l+1} u_l r dr)
!           
          call rint13(one,cf2,ulol,ulol,upl,upl,iuup,iat)
          call rint13(one,cf2,ulol,ulol,uor,uor,iuu,iat)
          iulol1ul(ilo,l,iat,isp)=iuup-fl*iuu
!          
!         iulol1ul = int(ulo_{l+1} d udot_l/dr r^2 dr)-l*int(ulo_{l+1} udot_l r dr)
!
          call rint13(one,cf2,ulol,ulol,udpl,udpl,iuup,iat)
          call rint13(one,cf2,ulol,ulol,udor,udor,iuu,iat)
          iulol1udl(ilo,l,iat,isp)=iuup-fl*iuu
        enddo 

!        
!------------------------------------------------------        
!       And now for local orbitals, only for l<= lomax
!------------------------------------------------------        
        do jlo=1,nLO_at(l,iat)
          ulopl(1:nr)=two*(uslo(1:nr,jlo,l,iat,isp))
          uloor(1:nr)=ulo(1:nr,jlo,l,iat,isp)/rr(1:nr)
            
          !! iul1ulol = int(u_{l+1} d ulo_l/dr r^2 dr) - l*int(u_{l+1} ulo_l r dr)
          call rint13(one,cf2,ul,ul,ulopl,ulopl,iuup,iat)
          call rint13(one,cf2,ul,ul,uloor,uloor,iuu,iat)
          iul1ulol(jlo,l,iat,isp)=iuup - fl*iuu

          !! iudl1ulol = int(udot_{l+1} d ulo_l/dr r^2 dr) - l*int(udot_{l+1} ulo_l r dr)
          call rint13(one,cf2,udl,udl,ulopl,ulopl,iuup,iat)
          call rint13(one,cf2,udl,udl,uloor,uloor,iuu,iat)
          iudl1ulol(jlo,l,iat,isp)=iuup-fl*iuu

          do ilo=1,nLO_at(l+1,iat) 
            ulol(1:nr)=ulo(1:nr,ilo,l+1,iat,isp)
!          
!           iulol1ulol = int(ulo_{l+1} d ulo_l/dr r^2 dr)-l*int(ulo_{l+1} ulo_l r dr)
!          
            call rint13(one,cf2,ulol,ulol,ulopl,ulopl,iuup,iat)
            call rint13(one,cf2,ulol,ulol,uloor,uloor,iuu,iat)
            iulol1ulol(ilo,jlo,l,iat,isp)=iuup-fl*iuu

          enddo 
        enddo 
                 
!=========================================================================
!                 end the calculation of iul1ul like integrals           !
!=========================================================================
        
!=========================================================================
!               begin the calculation of iulul1 like integrals           !
!=========================================================================
!
!       store the needed shared arrays in the local ones
!        
        ul(1:nr)=u(1:nr,l,iat,isp)
        upl(1:nr)=two*(us(1:nr,l+1,iat,isp))
        uor(1:nr)=u(1:nr,l+1,iat,isp)/rr(1:nr)
        udl(1:nr)=udot(1:nr,l,iat,isp)
        udpl(1:nr)=two*(usdot(1:nr,l+1,iat,isp))
        udor(1:nr)=udot(1:nr,l+1,iat,isp)/rr(1:nr)
!        
!       iulul1 = int(u_l} du_{l+1}/dr r^2 dr)+(l+2)*int(u_l} u_{l+1} r dr)
!        
        call rint13(one,cf2,ul,ul,upl,upl,iuup,iat)
        call rint13(one,cf2,ul,ul,uor,uor,iuu,iat)
        iulul1(l,iat,isp) = iuup + (fl+two)*iuu

!        
!       iuludl1 = int(u_l d udot_{l+1}/dr r^2 dr)+(l+2)*int(u_l udot_{l+1} r dr)
!        
        call rint13(one,cf2,ul,ul,udpl,udpl,iuup,iat)
        call rint13(one,cf2,ul,ul,udor,udor,iuu,iat)
        iuludl1(l,iat,isp) = iuup + (fl+two)*iuu
        
!        
!       iudlul1 = int(udot_l du_{l+1}/dr r^2 dr)+ (l+2)*int(udot_{l+1}u_{l+1} r dr)
!
        call rint13(one,cf2,udl,udl,upl,upl,iuup,iat)
        call rint13(one,cf2,udl,udl,uor,uor,iuu,iat)
        iudlul1(l,iat,isp) = iuup + (fl+two)*iuu
        
!        
!       iudludl1 = int(udot_l d udot_{l+1}/dr r^2 dr)+(l+2)*int(udot_l udot_{l+1} r dr)
!
        call rint13(one,cf2,udl,udl,udpl,udpl,iuup,iat)
        call rint13(one,cf2,udl,udl,udor,udor,iuu,iat)
        iudludl1(l,iat,isp) = iuup + (fl+two)*iuu

        if(ldbg) write(fid_outmom,10) l, iul1ul(l,iat,isp),&
     &   iul1udl(l,iat,isp),iudl1ul(l,iat,isp),iudl1udl(l,iat,isp),&
     &   iulul1(l,iat,isp), iuludl1(l,iat,isp),iudlul1(l,iat,isp), &
     &   iudludl1(l,iat,isp) 

        do ilo=1,nlo_at(l,iat) 
          ulol(1:nr)=ulo(1:nr,ilo,l,iat,isp)
!           
!         iulolul1 = int(ulo_l du_{l+1}/dr r^2 dr)+(l+2)*int(ulo_{l+1} u_l r dr)
!            
          call rint13(one,cf2,ulol,ulol,upl,upl,iuup,iat)
          call rint13(one,cf2,ulol,ulol,uor,uor,iuu,iat)
          iulolul1(ilo,l,iat,isp)=iuup+(fl+two)*iuu

!            
!         iulolul1 = int(ulo_l d udot_{l+1}/dr r^2 dr)+(l+2)*int(ulo_{l+1} udot_l r dr)
!            
          call rint13(one,cf2,ulol,ulol,udpl,udpl,iuup,iat)
          call rint13(one,cf2,ulol,ulol,udor,udor,iuu,iat)
          iuloludl1(ilo,l,iat,isp)=iuup+(fl+two)*iuu
        enddo 
        
!        
!------------------------------------------------------        
!       And now for local orbitals, l< lomax
!------------------------------------------------------        
        do jlo=1,nLO_at(l+1,iat) 
          ulopl(1:nr)=two*(uslo(1:nr,jlo,l+1,iat,isp))
          uloor(1:nr)=ulo(1:nr,jlo,l+1,iat,isp)/rr(1:nr)
!            
!         iululol1 = int(u_l d ulo_{l+1}/dr r^2 dr)+(l+2)*int(u_l ulo_{l+1} r dr)
!
          call rint13(one,cf2,ul,ul,ulopl,ulopl,iuup,iat)
          call rint13(one,cf2,ul,ul,uloor,uloor,iuu,iat)
          iululol1(jlo,l,iat,isp) = iuup+(fl+two)*iuu
!            
!         iudlulol1 = int(udot_l d ulo_{l+1}/dr r^2 dr)+(l+2)*int(udot_l ulo_{l+1} r dr)
!            
          call rint13(one,cf2,udl,udl,ulopl,ulopl,iuup,iat)
          call rint13(one,cf2,udl,udl,uloor,uloor,iuu,iat)
          iudlulol1(jlo,l,iat,isp) = iuup+(fl+two)*iuu

          do ilo=1,nLO_at(l,iat) 
            ulol(1:nr)=ulo(1:nr,ilo,l,iat,isp)
!          
!           iulol1ulol = int(ulo_l d ulo_{l+1}/dr r^2 dr)-(l+2)*int(ulo_l ulo_{l+1} r dr)
!           
            call rint13(one,cf2,ulol,ulol,ulopl,ulopl,iuup,iat)
            call rint13(one,cf2,ulol,ulol,uloor,uloor,iuu,iat)
            iulolulol1(ilo,jlo,l,iat,isp) = iuup+(fl+two)*iuu
          enddo 
        enddo

!=========================================================================
!                 end the calculation of iulul1 like integrals           !
!=========================================================================
      enddo ! l
!
!     Deallocate local arrays
!
      deallocate(rr)
      deallocate(ul)
      deallocate(upl)
      deallocate(uor)
      deallocate(udl)
      deallocate(udpl)
      deallocate(udor)
      deallocate(ulol)
      deallocate(ulopl)
      deallocate(uloor)

   10 format(i3,8(1pg15.7))
   11 format(5x,'radial integrals for the momentum matrix elements',   &
     &      ' of atom ',a10,/,  &
     &       2x,'l',4x,'iul1ul',8x,'iul1udl',8x,'iudl1ul',8x,'iudl1udl'&
     &      ,8x,'iulul1',8x,'iuludl1',8x,'iudlul1',8x,'iudludl1')  

      end subroutine momradintv
!EOC      
