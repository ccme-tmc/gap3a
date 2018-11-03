!BOP
!
! !ROUTINE: coul_setvm1
!
! !INTERFACE:
      subroutine coul_setvm1(iop,iq)
      
! !DESCRIPTION:
!
! This subroutine calculates the bare coulomb potential matrix by
!expanding it in plane waves and then transforming to the mixed basis.
!(For test purposes only)      
!
! !USES:
     
      use barcoul,    only: vmat,iop_coul,rcut_coul, &
      &                     zcut_coul,acut_coul,bcut_coul
      use constants,  only: cone, czero, fourpi,twopi
      use kpoints,    only: qlist, idvq
      use recipvec,   only: indgq, indgqlen, gqleng, ngqbarc
      use mixbasis,   only: mbsiz, locmatsiz,mpwmix

!
! !INPUT PARAMETERS: 

      implicit none
      integer, intent(in):: iop ! 0 -- calculate v matrix 
                                ! 1 -- calculate v^{1/2} matrix 
      integer, intent(in) :: iq ! index of the q-point

!
! !LOCAL VARIABLES:
      
      integer :: i 
      integer :: ipw
      integer :: ipin
      real(8) :: ev,qgl,ks_tf2
      real(8) :: pow
      complex(8), allocatable :: tmat(:,:)
      
! !EXTERNAL ROUTINES: 
 
      external zgemm
    
!EOP
!BOC
   
      allocate(tmat(mbsiz,ngqbarc(iq)))
      allocate(mpwmix(mbsiz,ngqbarc(iq)))
      if(iq.eq.1) then
        call coul_wmix0
      endif
      call coul_mpwmix(iq)

      vmat(:,:)=czero

      if(iop_coul.eq.4) ks_tf2 = (1.d0/rcut_coul)**2

      if(iop.eq.0) then 
        pow = 1.0
      else
        pow = 0.5
      endif 

      ipin=1
      if(iq.eq.1.and.iop_coul.eq.-1)then
        ipin=2
        tmat(1:mbsiz,1)=czero
      endif

      do ipw=ipin,ngqbarc(iq)
        qgl = gqleng(indgqlen(ipw,iq),iq)
        if(iop_coul.eq.-1) then 
          ev = fourpi/(qgl*qgl)
        elseif(iop_coul.eq.0) then 
          if(qgl.lt.1.0d-10) then 
            ev = twopi*rcut_coul**2
          else 
            ev = fourpi/(qgl*qgl)*(1.d0-cos(qgl*rcut_coul))
          endif 
        ! TODO iop_coul.eq.1
        ! TODO iop_coul.eq.2
        elseif(iop_coul.eq.3) then 
          ev = fourpi/(qgl*qgl+ks_tf2)
        endif 
 
        ev=ev**pow
        tmat(:,ipw) = mpwmix(:,ipw)*ev 
      enddo 
   
      call zgemm('n','c',mbsiz,mbsiz,ngqbarc(iq),cone,tmat,mbsiz,    &
     &            mpwmix,mbsiz,czero,vmat,mbsiz) 


      deallocate(tmat,mpwmix)
      
      
      end subroutine coul_setvm1
!EOC            
     
               
