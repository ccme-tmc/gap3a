!BOP
!
! !ROUTINE: coul_mpwipw
!
! !INTERFACE:
      subroutine coul_mpwipw(iq)
      
! !DESCRIPTION:
!
! This subroutine calculates the matrix elements between PW's and
! orthogonalized IPW's.
!
! !USES:

      use constants,  only: cone, czero
      use mixbasis,   only: igm,ipwint,sgi,mpwipw
      use recipvec,   only: npw,gindex, ig0, indgq, ngqbarc, ngq
      use task,       only: fid_outdbg
!
! !INPUT PARAMETERS:
      
      implicit none
      integer, intent(in) :: iq

! !LOCAL VARIABLES:

      integer :: info
      integer :: ipw    ! (Counter): runs over plane waves
      integer :: jipw   ! (Coutner): runs over IPW's
      integer :: iprt=0

      integer, dimension(3) :: ig ! integer coodintates of G-G'
      
      complex(8), allocatable :: tmat(:,:)
      character*20:: sname= "coul_mpwipw"

 
! !EXTERNAL ROUTINES: 


      external zgemm
!
! !REVISION HISTORY:
! 
! Created: 5th. Jan. 2005 by RGA

!EOP
!BOC

!
!     Calculate the interstitial mixed basis functions
!
      if(iprt.gt.0) write(6,*) "coul_mpwipw: diagipw"
      call diagipw(iq)

      allocate(tmat(ngq(iq),ngqbarc(iq)),stat=info)
      call errmsg(info.ne.0,sname,"Fail to allocate tmat")
!
!     Calculate the integral between pw's and IPW's
!

      if(iprt.gt.0) write(6,*) "coul_mpwipw: set tmat"
      do ipw=1,ngqbarc(iq)
        do jipw=1,ngq(iq)
          ig(1:3)=gindex(:,indgq(ipw,iq))-gindex(:,indgq(jipw,iq))
          if(abs(ig(1)).gt.igm(1))write(*,*)'danger!! ig(1)=',ig(1)
          if(abs(ig(2)).gt.igm(2))write(*,*)'danger!! ig(2)=',ig(2)
          if(abs(ig(3)).gt.igm(3))write(*,*)'danger!! ig(3)=',ig(3)
          tmat(jipw,ipw)=ipwint(ig0(ig(1),ig(2),ig(3)))
        enddo ! jipw  
      enddo ! ipw  

      if(iprt.gt.0) write(6,*) "coul_mpwipw: call zgemm"
      call zgemm('t','n',ngq(iq),ngqbarc(iq),ngq(iq),cone,sgi,ngq(iq),  & 
     &           tmat,ngq(iq),czero,mpwipw,ngq(iq))

#ifdef DEBUG
      write(fid_outdbg,*) "### mpwipw for iq=",iq 
      write(fid_outdbg,*)   
      write(fid_outdbg,'(2a5,2a15,2a24)') 'jipw','ipw','gindex(ipw)',   &
     &                            'gindex(jipw)','ipwint',              &
     &                            'mpwipw(jipw,ipw)'

      do ipw=1,ngqbarc(iq),ngqbarc(iq)/10
        do jipw=1,ngq(iq),ngq(iq)/10
          ig(1:3)=gindex(:,indgq(ipw,iq))-gindex(:,indgq(jipw,iq))
          write(fid_outdbg,'(8i5,4f12.6)') jipw,ipw,                    &
     &             gindex(:,indgq(ipw,iq)),gindex(:,indgq(jipw,iq)),    &
     &            ipwint(ig0(ig(1),ig(2),ig(3))),mpwipw(jipw,ipw)
        enddo 
      enddo 

      write(fid_outdbg,*) "### sgi ###"
      do jipw=1,ngq(iq),ngq(iq)/10
        do ipw=1,ngq(iq),ngq(iq)/10
          write(fid_outdbg,'(2i5,4f12.6)') ipw,jipw,sgi(ipw,jipw)
        enddo 
      enddo  
#endif 
!
!     Transform to the mixed basis (orthogonalized IPW's)
      deallocate(tmat)
      return
 
      end subroutine coul_mpwipw
!EOC
