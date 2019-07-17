!BOP
!
! !ROUTINE: calcrwf
!
! !INTERFACE: 
      subroutine calcrwf(iat,isp)
      
! !DESCRIPTION:
!
!This subroutine calculates all the radial wave functions $u_l(r,E_l)$, 
!$\dot{u}_l(r,E_l)$ and $u_{lo}(r,E_l)$ for atom \verb"iat".
!
! !USES:
      
      use constants, only: cfein1,cfein2
      use lapwlo,    only: lmax,lomax,nt,nLO_at,umt,umtlo,abcelo,lapw
      use radwf,     only: a,u,us,b,udot,ae,usdot,be,ulo,uslo
      use struk,     only: zz,atomname,nrpt
      use task,      only: fid_outmb
      use modmpi,    only: myrank
      
! !INPUT PARAMETERS:
      
      implicit none
      
      integer(4), intent(in) :: iat  ! Index of the inequivalent atom
      integer(4), intent(in) :: isp  ! Index of spin
      
! !LOCAL VARIABLES:

      integer(4) :: npt   ! Number of radial mesh points  (for each atom the corresponding nrpt(iat)  is stored here).
      integer(4) :: l     ! angular momentum quantum number
      integer(4) :: ilo   !  index for LOs 
      integer(4) :: nLO_a ! the number of LO per atom per l
      integer(4) :: irp   ! counter of the radial mesh points
      integer(4) :: nodel ! number of nodes of u_l(r,E_l- Delta E)
      integer(4) :: nodes ! number of nodes of u_l(r,E_l)
      integer(4) :: nodeu ! number of nodes of u_l(r,E_l+ Delta E)
      real(8)    :: el    ! expansion energy for the radial wavefunction
      

 
! !EXTERNAL ROUTINES: 

      
      external calcu
      external calcudot
      external radmesh
      external rint13

!
! !REVISION HISTORY:
!
! Created 9. July 2004 by RGA
!
!EOP
!BOC
      npt=nrpt(iat)  

      if(myrank.eq.0) then 
         write(fid_outmb,*) "###Calc radfunc for atom ",iat
         write(fid_outmb,*) "-- (L)APW states --"
      endif 

      do l = 0,lmax 
        el = 0.5d+0*umt(1,l,iat,isp)
!
!       Calculate the function at el (u_l(r,E_l)):
!
        call calcu(iat,zz(iat),l,el,umt(2,l,iat,isp),umt(4,l,iat,isp),  &
     &             nodes,isp)

        u(1:npt,l,iat,isp)  = a(1:npt)
        us(1:npt,l,iat,isp) = b(1:npt)
!
!       Calculate the energy derivative of the radial function
!
        call calcudot(iat,zz(iat),l,el,nodel,nodeu,isp)

        udot(1:npt,l,iat,isp) = ae(1:npt) 
        usdot(1:npt,l,iat,isp)= be(1:npt) 

        if(myrank.eq.0) write(fid_outmb,6030) l, umt(2,l,iat,isp),     &
     &                    umt(4,l,iat,isp),                             &
     &                    umt(3,l,iat,isp), umt(5,l,iat,isp),           &
     &                    umt(6,l,iat,isp), nodel, nodes, nodeu
      enddo 


!
!     Calculate the radial function for local orbitals.
!
      if(myrank.eq.0) write(fid_outmb,*) "-- LO states --"

      do l = 0, lomax
        do ilo= 1, nLO_at(l,iat) 
          el = 0.50D+0*abcelo(4, ilo, l,iat,isp)
          call calcu(iat,zz(iat),l,el,umtlo(1,ilo,l,iat,isp),               &
     &               umtlo(2,ilo,l,iat,isp),nodes,isp)

!
!         Calculate pi12lo(l,iat)=<u_l(r,E_l)|u_l(r,E_{lo})>:
!
          call rint13(cfein1,cfein2,u(1,l,iat,isp),us(1,l,iat,isp),a,b, &
     &          umtlo(3,ilo,l,iat,isp),iat)
!
!         Calculate pe12lo(l,iat)=<\dot{u}_l(r,E_l)|u_l(r,E_{lo})>:
!
          call rint13(cfein1,cfein2,udot(1,l,iat,isp),                  &
     &               usdot(1,l,iat,isp),a,b,umtlo(4,ilo,l,iat,isp),iat)
!
!         Store the radial wave function :
!
          do irp = 1, npt
            ulo( irp,ilo,l,iat,isp) = a(irp)
            uslo(irp,ilo,l,iat,isp) = b(irp)
          enddo

          if(myrank.eq.0) write(fid_outmb,6031) ilo, l, el,             &
     &                   umtlo(1:4,ilo,l,iat,isp),nodel, nodes, nodeu
        enddo
      enddo 

 6020 format (/,10x,'Potential parameters for atom ',a10,/,11x,'l',5x,  &
     &  'u(r)',10x, 'u''(r)',9x,'du/de',8x,'du''/de',6x,'norm-u''')
 6021 format (/,10x,'Local orbital potential parameters for atom',a10,/,&
     &  11x,'l',5x,'u(r)',10x,  'u''(r)',5x,'norm u1u2',8x,'norm ue1u2')
 6030 format (10x,i2,5e14.6,5x,3i2)
 6031 format (10x,2i2,5e14.6,5x,3i4)

      
      end subroutine calcrwf
!EOC       
        
