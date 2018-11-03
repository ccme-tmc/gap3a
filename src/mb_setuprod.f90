!BOP
!
! !ROUTINE: mb_setuprod
!
! !INTERFACE:
       subroutine mb_setuprod(iat,npt)
       
! !DESCRIPTION:
!
! This subroutine calculates the radial product functions for atom iat,
! and their overlap matrix
!
!
! !USES:
      use bands,     only: nspin 
      use constants, only: cfein1,cfein2
      use core,      only: ncore,lcore,eigcore,iop_core,ucore,uscore
      use lapwlo,    only: lomax,nLO_at,lapw,elo
      use mixbasis,  only: lmbmax_at,nspin_mb,mb_ludot,mb_emin,mb_emax
      use prodfun,   only: eles, nup, rp, umat, uprod,usprod
      use radwf,     only: a,b,u,ulo,us,uslo,udot,usdot
      use task,      only: fid_outmb

! !INPUT PARAMETERS:

      implicit none

      integer, intent(in) :: iat ! Atom for which the product functions are calculated
      integer, intent(in) :: npt ! Number of radial grid points, 
                                 ! if equal 0, then just count nup without doing anything else
 
! !LOCAL VARIABLES:

      integer :: iup,jup   !* Index for radial function (u) product 
      integer :: irf,jrf       !* Index for radial mesh 
      integer :: l1,l2      !* Index for angular momentum quantum number 
      integer :: isp       !* Index for spin 
      integer :: ierr 
      integer :: nrf   ! the total number of radial functions 

      real(8), dimension(1:npt) :: c,d,u1,u2,us1,us2
      real(8),allocatable:: u_all(:,:),us_all(:,:)
      integer,allocatable:: l_all(:)  
      logical ::ldbg = .false.
      character(len=20):: sname="mb_setuprod"
      
      
      external rint13

! !REVISION HISTORY:
!
! Created 18. May 2004 by RGA
!
!
!EOP
!BOC            

!* --------------------------------------------------------------------!
!* first collect all radial functions to be considered to setup        !
!* mixed-basis                                                         !
!----------------------------------------------------------------------!

      ! count the number of radial functions to be considered       
      call sub_count_nrf(nrf) 
      if(npt.eq.0) then 
        nup = nrf*(nrf+1)/2*nspin_mb
        write(fid_outmb,*) "mb_setuprod: # of uprod for atom",iat,nup
        return 
      endif
      
      !* Collect all radial functions to be considered  
      allocate(u_all(npt,nrf),us_all(npt,nrf),l_all(nrf),stat=ierr) 
      call errmsg(ierr.ne.0,sname,"Fail to allocate u_all ...!")
     
      do isp=1,min(nspin_mb,nspin)
        call sub_collect_u(isp) 

        iup = 0 
        do jrf=1, nrf
          l2 = l_all(jrf) 
          u2 = u_all(:,jrf) 
          us2 = us_all(:,jrf) 

          do irf = 1, jrf 
            l1 = l_all(irf) 
            u1 = u_all(:,irf) 
            us1 = us_all(:,irf) 

            iup = iup + 1

            eles(1,iup) = l1 
            eles(2,iup) = l2  
            uprod(:,iup) = u1*u2/rp
            usprod(:,iup) = us1*us2/rp 
          enddo 
        enddo
      enddo   
      deallocate(u_all,us_all,l_all) 

!
!     calculate the overlap matrix of product functions      
!
      do iup=1,nup
        a(1:npt)=uprod(1:npt,iup)
        b(1:npt)=usprod(1:npt,iup)
        do jup=iup,nup
          c(1:npt)=uprod(1:npt,jup)
          d(1:npt)=usprod(1:npt,jup)
          call rint13(cfein1,cfein2,a,b,c,d,umat(iup,jup),iat)
        enddo  
      enddo 


      contains 
 
!
!     count the total number of radial functions to be considered 
!
      subroutine sub_count_nrf(nrf)
      integer, intent(out):: nrf
      integer:: ic,l,nLO
      
      ! count the total number 
      nrf = 0

      ! core states 
      do ic=1,ncore(iat)
        if(eigcore(ic,iat,1).gt.mb_emin) nrf = nrf+1
      enddo

      ! normal LAPW radial functions 
      if(mb_ludot) then
        nrf = nrf + (lmbmax_at(iat)+1)*2
      else
        nrf = nrf + (lmbmax_at(iat)+1)
      endif

      ! LO basis 
      do l=0,min(lomax,lmbmax_at(iat)) 
        call sub_set_nLO(l,nLO,1)
        nrf = nrf + nLO
      enddo

      write(6,100) iat,nrf 
 100  format("# rad. func. to be considered on the atom ",i5,":",i5)  

      endsubroutine sub_count_nrf

!
!     Set the number of LOs for each l 
!     
      subroutine sub_set_nLO(l,nLO_a,isp)
      integer :: l,nLO_a,isp
      integer :: ilo
      
      if(l.gt.lomax) then
        nLO_a = 0
      else
        nlo_a = 0
        do ilo=1,nlo_at(l,iat) 
          if(elo(l,ilo,iat,isp).lt.mb_emax) then 
            nlo_a = ilo 
          else
            exit 
          endif 
        enddo 
      endif 
      end subroutine 

!
!     Collect all radial functions 
!
      subroutine sub_collect_u(isp) 
      integer, intent(in):: isp 
      integer:: irf,ic,l,nlo,ilo

      irf = 0 
      ! core states 
      do ic=1,ncore(iat)
        if(eigcore(ic,iat,1).gt.mb_emin) then 
          irf = irf + 1 
          l_all(irf) = lcore(ic,iat)
          u_all(:,irf) = ucore(1:npt,ic,iat,isp)
          us_all(:,irf) = uscore(1:npt,ic,iat,isp)
        endif 
      enddo 

      do l=0,lmbmax_at(iat)
        irf = irf + 1 
        l_all(irf) = l 
        u_all(:,irf) = u(1:npt,l,iat,isp) 
        us_all(:,irf) = us(1:npt,l,iat,isp) 

        if(mb_ludot) then
          irf = irf + 1
          l_all(irf) = l
          u_all(:,irf) = udot(1:npt,l,iat,isp)
          us_all(:,irf) = usdot(1:npt,l,iat,isp)
        endif 

        call sub_set_nLO(l,nLO,isp)
        do ilo=1,nLO
          irf = irf + 1
          l_all(irf) = l
          u_all(:,irf) = ulo(1:npt,ilo,l,iat,isp)
          us_all(:,irf) = uslo(1:npt,ilo,l,iat,isp)
        enddo
      enddo  ! l

      end subroutine 

     
      end subroutine mb_setuprod
!EOC
