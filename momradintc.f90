!BOP
!
! !ROUTINE: momradintc
!
! !INTERFACE:
      subroutine momradintc(iat,isp)
      
! !DESCRIPTION:
!      
! This subroutine calculates the set of integrals of radial functions
! needed for the calculation of the momentum matrix elements (Only those
! including core states). 
! 
! The naming convention for the corresponding integrals is:
! 
! \begin{subequations}\label{momradintc1}
!\begin{align}
! \texttt{iucl1ul(iat,ic)}=&\int\limits_0^{R^a_{MT}}{%
! u^{\mathrm{core}}_l(r^a)u'_{l-1}(r^a)(r^a)^2dr}-
! (l-1)\int\limits_0^{R^a_{MT}}{u^{\mathrm{core}}_l(r^a)u_{l-1}(r^a)r^adr}\\
! \texttt{iul1ucl(iat,ic)}=&\int\limits_0^{R^a_{MT}}{%
! u_{l+1}(r^a)u'^{\mathrm{core}}_l(r^a)(r^a)^2dr}-
! l\int\limits_0^{R^a_{MT}}{u_{l+1}(r^a)u^{\mathrm{core}}_l(r^a)r^adr}\\
! \texttt{iucl1udl(iat,ic)}=&\int\limits_0^{R^a_{MT}}{%
! u^{\mathrm{core}}_l(r^a)\dot{u}'_{l-1}(r^a)(r^a)^2dr}-
! (l-1)\int\limits_0^{R^a_{MT}}{u^{\mathrm{core}}_l(r^a)\dot{u}_{l-1}(r^a)r^adr}\\
! \texttt{iudl1ucl(iat,ic)}=&\int\limits_0^{R^a_{MT}}{%
! \dot{u}_{l+1}(r^a)u'^{\mathrm{core}}_l(r^a)(r^a)^2dr}-
! l\int\limits_0^{R^a_{MT}}{\dot{u}_{l+1}(r^a)u^{\mathrm{core}}_l(r^a)r^adr}\\
! \texttt{iucl1ulol(iat,ic)}=&\int\limits_0^{R^a_{MT}}{%
! u^{\mathrm{core}}_l(r^a)u'^{lo}_{l-1}(r^a)(r^a)^2dr}-
! (l-1)\int\limits_0^{R^a_{MT}}{u^{\mathrm{core}}_{l}(r^a)u^{lo}_{l-1}(r^a)r^adr}\\
! \texttt{iulol1ucl(iat,ic)}=&\int\limits_0^{R^a_{MT}}{%
! u^{lo}_{l+1}(r^a)u'^{\mathrm{core}}_l(r^a)(r^a)^2dr}-
! l\int\limits_0^{R^a_{MT}}{u^{lo}_{l+1}(r^a)u^{\mathrm{core}}_l(r^a)r^adr}\\
! \texttt{iucl1ucl(iat,ic,jc)}=&\int\limits_0^{R^a_{MT}}{%
! u^{\mathrm{core}}_{l+1}(r^a)u'^{\mathrm{core}}_l(r^a)(r^a)^2dr}-
! l\int\limits_0^{R^a_{MT}}{u^{\mathrm{core}}_{l+1}(r^a)u^{\mathrm{core}}_l(r^a)r^adr}
! \end{align}
! \end{subequations}
!
! \begin{subequations}\label{momradintc2}
!\begin{align}
! \texttt{iuclul1(iat,ic)}=&\int\limits_0^{R^a_{MT}}{%
! u^{\mathrm{core}}_{l}(r^a)u'_{l+1}(r^a)(r^a)^2dr}+
! (l+2)\int\limits_0^{R^a_{MT}}{u^{\mathrm{core}}_{l}(r^a)u_{l+1}(r^a)r^adr}\\
! \texttt{iulucl1(iat,ic)}=&\int\limits_0^{R^a_{MT}}{%
! u_{l-1}(r^a)u'^{\mathrm{core}}_{l}(r^a)(r^a)^2dr}+
! (l+1)\int\limits_0^{R^a_{MT}}{u_{l-1}(r^a)u^{\mathrm{core}}_{l}(r^a)r^adr}\\
! \texttt{iucludl1(iat,ic)}=&\int\limits_0^{R^a_{MT}}{%
! u^{\mathrm{core}}_{l}(r^a)\dot{u}'_{l+1}(r^a)(r^a)^2dr}+
! (l+2)\int\limits_0^{R^a_{MT}}{u^{\mathrm{core}}_{l}(r^a)\dot{u}_{l+1}(r^a)r^adr}\\
! \texttt{iudlucl1(iat,ic)}=&\int\limits_0^{R^a_{MT}}{%
! \dot{u}_{l-1}(r^a)u'^{\mathrm{core}}_{l}(r^a)(r^a)^2dr}+
! (l+1)\int\limits_0^{R^a_{MT}}{\dot{u}_{l-1}(r^a)u^{\mathrm{core}}_{l}(r^a)r^adr}\\
! \texttt{iuclulol1(iat,ic)}=&\int\limits_0^{R^a_{MT}}{%
! u^{\mathrm{core}}_{l}(r^a)u'^{lo}_{l+1}(r^a)(r^a)^2dr}+
! (l+2)\int\limits_0^{R^a_{MT}}{u^{\mathrm{core}}_{l}(r^a)u^{lo}_{l+1}(r^a)r^adr}\\
! \texttt{iulolucl1(iat,ic)}=&\int\limits_0^{R^a_{MT}}{%
! u^{lo}_{l-1}(r^a)u'^{\mathrm{core}}_{l}(r^a)(r^a)^2dr}+
! (l+1)\int\limits_0^{R^a_{MT}}{u^{lo}_{l-1}(r^a)u^{\mathrm{core}}_{l}(r^a)r^adr}\\
! \texttt{iuclucl1(iat,ic,jc)}=&\int\limits_0^{R^a_{MT}}{%
! u^{\mathrm{core}}_{l}(r^a)u'^{\mathrm{core}}_{l+1}(r^a)(r^a)^2dr}+
! (l+2)\int\limits_0^{R^a_{MT}}{u^{\mathrm{core}}_{l}(r^a)u^{\mathrm{core}}_{l+1}(r^a)r^adr}
! \end{align}
! \end{subequations}

! !USES:

      use core,    only: ncore,lcore,ucore,uscore
      use lapwlo,  only: nt,lomax,nlo_at
      use mommat, only: iul1ucl, iudl1ucl, iulol1ucl, iucl1ul,       &
     &                  iucl1udl, iucl1ulol, iuclul1, iucludl1,      &
     &                  iuclulol1, iulucl1, iudlucl1, iulolucl1,     &
     &                  iucl1ucl, iuclucl1
      use radwf,   only: u, ulo, udot, us, uslo, usdot
      use struk,   only: nat,nrpt
      

! !INPUT PARAMETERS:
      implicit none
      
      integer, intent(in) :: iat  ! index inequivalent atoms
      integer, intent(in) :: isp  ! index for spin 
      
! !LOCAL VARIABLES:
  
      integer :: ic,jc   ! Runs over core states
      integer :: l,l2    ! Angular momentum quantum number
      integer :: npt     ! number of radial mesh points

      integer :: ilo,jlo ! index for LO 
      integer :: m,m0,npo2
      real(8) :: fl,fl2      ! angular momentum quantum number
      real(8) :: iuup    ! \int(u_l du_l'/dr r^2 dr)
      real(8) :: iuu    ! \int(u_l u_l' r dr)
      
      real(8), dimension(4) :: c
      real(8), allocatable :: rr(:) ! Coordinate of the radial mesh
      real(8), allocatable :: ucl(:) ! core radial function u_l
      real(8), allocatable :: ucpl(:) ! radial derivative of the core radial function du_l/dr
      real(8), allocatable :: ucor(:) ! core radial function over r: u_l/r
      real(8), allocatable :: ul(:) ! radial function u_l
      real(8), allocatable :: upl(:) ! radial derivative of the radial function du_l/dr
      real(8), allocatable :: uor(:) ! radial function over r: u_l/r
      real(8), allocatable :: udl(:) ! energy derivative of the radial function udot_l      
      real(8), allocatable :: udpl(:) ! radial derivative of udot: d udot_{l}/dr
      real(8), allocatable :: udor(:) ! energy derivative of the radial function over r : udot_l/r
      real(8), allocatable :: ulol(:) ! local orbital radial function ulo_l      
      real(8), allocatable :: ulopl(:) ! radial derivative of ulo: d ulo_{l}/dr
      real(8), allocatable :: uloor(:) ! local orbital the radial function over r : udot_l/r
      
! !DEFINED PARAMETERS:

      integer, parameter :: np = 4
      real(8), parameter :: two = 2.0d+0
      real(8), parameter :: one = 1.0d+0
      real(8), parameter :: cf2 = 1.0d-22

! !EXTERNAL ROUTINES: 

      real(8), external :: polynom
      external rint13
      external radmesh
      

! !INTRINSIC ROUTINES: 
 
      intrinsic iabs
      intrinsic isign
! !REVISION HISTORY:
!
! Created: 16. August 2004 by RGA
! Last modfied: Oct. 22, 2013 by JH
!
!EOP
!BOC
!
!     Allocate necesary arrays
!     
      npt=nrpt(iat)
      allocate(rr(npt))
      allocate(ucl(npt))
      allocate(ucpl(npt))
      allocate(ucor(npt))
      allocate(ul(npt))
      allocate(upl(npt))
      allocate(uor(npt))
      allocate(udl(npt))
      allocate(udpl(npt))
      allocate(udor(npt))
      allocate(ulol(npt))
      allocate(ulopl(npt))
      allocate(uloor(npt))
!
!     calculate the radial mesh points
!
      call radmesh(iat,rr)
!      
!     loop over l      
!      
      npo2=np/2
      do ic=1,ncore(iat)
        l = lcore(ic,iat) 
        fl= dble(l)

!=========================================================================
!               begin the calculation of iul1ul like integrals           !
!=========================================================================
!
!      store the needed shared arrays in the local ones
!        
        ucl(1:npt) = ucore(1:npt,ic,iat,isp)
        ucor(1:npt)= ucore(1:npt,ic,iat,isp)/rr(1:npt)
        ul(1:npt)  = u(1:npt,l+1,iat,isp)
        udl(1:npt) = udot(1:npt,l+1,iat,isp)
!
!       calculate the radial derivative of udore: ucpl = r*ducore_l/dr        
!
        ucpl(1)=0.0d0
        do m=2,npt
          if(m.le.npo2)then
            m0=1
          elseif(m.gt.npt-npo2)then
            m0=npt-np+1
          else
            m0=m-npo2
          endif
          ucpl(m) = polynom(1,np,rr(m0),ucl(m0),c,rr(m)) - ucor(m)
        enddo !m
!        
!       iul1ul = int(u_{l+1} ducore_l/dr r^2 dr)-l*int(u_{l+1} ucore_l r dr)
!        
        call rint13(one,cf2,ul,ul,ucpl,ucpl,iuup,iat)
        call rint13(one,cf2,ul,ul,ucor,ucor,iuu,iat)
        iul1ucl(ic,iat,isp)=iuup-fl*iuu

!        
!       iudl1ucl = int(udot_{l+1} ducore_l/dr r^2 dr)-l*int(udot_{l+1} ucore_l r dr)
!        
        call rint13(one,cf2,udl,udl,ucpl,ucpl,iuup,iat)
        call rint13(one,cf2,udl,udl,ucor,ucor,iuu,iat)
        iudl1ucl(ic,iat,isp)=iuup-fl*iuu
!        
!------------------------------------------------------        
!       And now for local orbitals, l< lomax
!------------------------------------------------------        
        do ilo=1,nLO_at(l+1,iat) 
          ulol(1:npt)=ulo(1:npt,ilo,l+1,iat,isp)
!        
!         iulol1ucl = int(ulo_{l+1} ducore_l/dr r^2 dr) - l*int(ulo_{l+1} ucore_l r dr)
!
          call rint13(one,cf2,ulol,ulol,ucpl,ucpl,iuup,iat)
          call rint13(one,cf2,ulol,ulol,ucor,ucor,iuu,iat)
          iulol1ucl(ilo,ic,iat,isp)=iuup-fl*iuu
        enddo 

        if(l.gt.0)then
          upl(1:npt) = two*(us(1:npt,l-1,iat,isp))
          uor(1:npt) = u(1:npt,l-1,iat,isp)/rr(1:npt)
          udpl(1:npt)= two*(usdot(1:npt,l-1,iat,isp))
          udor(1:npt)= udot(1:npt,l-1,iat,isp)/rr(1:npt)
!        
!         iucl1ul = int( ucore_l du_{l-1}/dr r^2 dr)-(l-1)*int(ucore_l u_{l-1} r dr)
!        
          call rint13(one,cf2,ucl,ucl,upl,upl,iuup,iat)
          call rint13(one,cf2,ucl,ucl,uor,uor,iuu,iat)
          iucl1ul(ic,iat,isp) = iuup-(fl-one)*iuu
!        
!         iul1udl = int(ucore_l d udot_{l-1}/dr r^2 dr)-(l-1)*int(ucore_l udot_{l-1} r dr)
!        
          call rint13(one,cf2,ucl,ucl,udpl,udpl,iuup,iat)
          call rint13(one,cf2,ucl,ucl,udor,udor,iuu,iat)
          iucl1udl(ic,iat,isp) = iuup-(fl-one)*iuu
        
!------------------------------------------------------        
!       And now for local orbitals, l< lomax
!------------------------------------------------------        
          do jlo=1,nLO_at(l-1,iat) 
            ulopl(1:npt)=two*(uslo(1:npt,jlo,l-1,iat,isp))
            uloor(1:npt)=ulo(1:npt,jlo,l-1,iat,isp)/rr(1:npt)
!            
!           iucl1ulol = int(ucore_l d ulo_{l-1}/dr r^2 dr)-(l-1)*int(ucore_l ulo_{l-1} r dr)
!            
            call rint13(one,cf2,ucl,ucl,ulopl,ulopl,iuup,iat)
            call rint13(one,cf2,ucl,ucl,uloor,uloor,iuu,iat)
            iucl1ulol(jlo,ic,iat,isp)=iuup-(fl-one)*iuu
          enddo 
        endif ! l > 0
          
!=========================================================================
!                 end the calculation of iul1ul like integrals           !
!=========================================================================
        
!=========================================================================
!               begin the calculation of iulul1 like integrals           !
!=========================================================================
!
!      store the needed shared arrays in the local ones
!        
        upl(1:npt)=two*(us(1:npt,l+1,iat,isp))
        uor(1:npt)=u(1:npt,l+1,iat,isp)/rr(1:npt)
        udpl(1:npt)=two*(usdot(1:npt,l+1,iat,isp))
        udor(1:npt)=udot(1:npt,l+1,iat,isp)/rr(1:npt)
        call rint13(one,cf2,ucl,ucl,upl,upl,iuup,iat)
        call rint13(one,cf2,ucl,ucl,uor,uor,iuu,iat)
!        
!       iulul1 = int(ucore_l} du_{l+1}/dr r^2 dr)+(l+2)*int(ucore_l} u_{l+1} r dr)
!        
        iuclul1(ic,iat,isp)=iuup+(fl+two)*iuu
!        
!       iuludl1 = int(ucore_l d udot_{l+1}/dr r^2 dr)+ (l+2)*int(ucore_l udot_{l+1} r dr)
!        
        call rint13(one,cf2,ucl,ucl,udpl,udpl,iuup,iat)
        call rint13(one,cf2,ucl,ucl,udor,udor,iuu,iat)
        iucludl1(ic,iat,isp)=iuup+(fl+two)*iuu
!        
!------------------------------------------------------        
!       And now for local orbitals, l< lomax
!------------------------------------------------------        
        do ilo=1, nLO_at(l+1,iat) 
          ulopl(1:npt)=two*(uslo(1:npt,ilo,l+1,iat,isp))
          uloor(1:npt)=ulo(1:npt,ilo,l+1,iat,isp)/rr(1:npt)
!            
!           iuclulol1 = int(ucore_l d ulo_{l+1}/dr r^2 dr)+(l+2)*int(ucore_l ulo_{l+1} r dr)
!            
          call rint13(one,cf2,ucl,ucl,ulopl,ulopl,iuup,iat)
          call rint13(one,cf2,ucl,ucl,uloor,uloor,iuu,iat)
          iuclulol1(ilo,ic,iat,isp)=iuup+(fl+two)*iuu
        enddo ! loor
        
        if(l.gt.0)then
          ul(1:npt)=u(1:npt,l-1,iat,isp)
          udl(1:npt)=udot(1:npt,l-1,iat,isp)
!        
!         iudlucl1 = int(u_{l-1} ducore_l/dr r^2 dr)+(l+1)*int(u_{l-1} ucore_l r dr)
!        
          call rint13(one,cf2,ul,ul,ucpl,ucpl,iuup,iat)
          call rint13(one,cf2,ul,ul,ucor,ucor,iuu,iat)
          iulucl1(ic,iat,isp)=iuup+(fl+one)*iuu
        
!        
!         iudlucl1 = int(udot_{l-1} ducore_l/dr r^2 dr)+ (l+1)*int(udot_{l-1} ucore_l r dr)
!        
          call rint13(one,cf2,udl,udl,ucpl,ucpl,iuup,iat)
          call rint13(one,cf2,udl,udl,ucor,ucor,iuu,iat)
          iudlucl1(ic,iat,isp)=iuup+(fl+one)*iuu
        
!        
!------------------------------------------------------        
!       And now for local orbitals, l< lomax
!------------------------------------------------------        
          do ilo=1,nLO_at(l-1,iat) 
            ulol(1:npt)=ulo(1:npt,ilo,l-1,iat,isp)
!           
!           iulolucl1 = int(ulo_{l-1} ducore_l/dr r^2 dr)+(l+1)*int(ulo_{l-1} ucore_l r dr)
!            
            call rint13(one,cf2,ulol,ulol,ucpl,ucpl,iuup,iat)
            call rint13(one,cf2,ulol,ulol,ucor,ucor,iuu,iat)
            iulolucl1(ilo,ic,iat,isp)=iuup+(fl+one)*iuu
          enddo 
        endif ! l > 0

!---------------------------------------------------
!       ANd now the core core elements
!---------------------------------------------------        

        do jc=1,ncore(iat)
          l2= lcore(jc,iat) 
          fl2=dble(l2)
          if(l2.eq.l-1)then
            ul(1:npt)=ucore(1:npt,jc,iat,isp)
            uor(1:npt)=ucore(1:npt,jc,iat,isp)/rr(1:npt)
!
!           calculate the radial derivative of ucore(jc)
!            
            upl(1)=0.0d0
            do m=2,npt
              if(m.le.npo2)then
                m0=1
              elseif(m.gt.npt-npo2)then
                m0=npt-np+1
              else
                m0=m-npo2
              endif
              upl(m)=polynom(1,np,rr(m0),ul(m0),c,rr(m))-uor(m)
            enddo !m
!        
!           iucl1ul = int(ucore_l ducore_{l2}/dr r^2 dr)-(l2)*int(ucore_l ucore_{l2} r dr)
!        
            call rint13(one,cf2,ucl,ucl,upl,upl,iuup,iat)
            call rint13(one,cf2,ucl,ucl,uor,uor,iuu,iat)
            iucl1ucl(ic,jc,iat,isp)=iuup-fl2*iuu

          elseif(l2.eq.l+1)then  
            ul(1:npt)=ucore(1:npt,jc,iat,isp)
            uor(1:npt)=ucore(1:npt,jc,iat,isp)/rr(1:npt)
!
!           calculate the radial derivative of ucore(jc)
!            
            upl(1)=0.0d0
            do m=2,npt
              if(m.le.npo2)then
                m0=1
              elseif(m.gt.npt-npo2)then
                m0=npt-np+1
              else
                m0=m-npo2
              endif
              upl(m)=polynom(1,np,rr(m0),ul(m0),c,rr(m))-uor(m)
            enddo !m
!        
!           iudlucl1 = int(u_{l-1} ducore_l/dr r^2 dr)+(l+1)*int(u_{l-1} ucore_l r dr)
            call rint13(one,cf2,ucl,ucl,upl,upl,iuup,iat)
            call rint13(one,cf2,ucl,ucl,uor,uor,iuu,iat)
            iuclucl1(ic,jc,iat,isp)=iuup+(fl+two)*iuu
          endif ! l2
        enddo ! jc  
      enddo ! ic
!
!     Deallocate local arrays
!
      deallocate(rr)
      deallocate(ucl)
      deallocate(ucpl)
      deallocate(ucor)
      deallocate(ul)
      deallocate(upl)
      deallocate(uor)
      deallocate(udl)
      deallocate(udpl)
      deallocate(udor)
      deallocate(ulol)
      deallocate(ulopl)
      deallocate(uloor)

      end subroutine momradintc
!EOC      
