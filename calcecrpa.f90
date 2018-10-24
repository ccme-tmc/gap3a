!BOP
!
! !ROUTINE: calcecrpa
!
! !INTERFACE:
      subroutine calcecrpa(iq,iomfirst,iomlast) 

!
! !DESCRIPTION: 
!
!  This subroutine calculate the \omega- and q-decompsed RPA correlation energy    
!  E_c^{RPA} ( i\omega,q), the final total RPA correlation energy will be 
!          Ec= N_c^{-1} \sum_\vecwq \frac{1}{2\pi} \int d\omega E_c^{RPA}(i\omega,\vecwq)
!
! !USES:

      use acfd,        only: ec_acfd
      use constants,   only: cone, czero,twopi
      use freq,        only: omega,womeg
      use bzinteg,     only: kwt_bz
      use mixbasis,    only: matsiz
      use dielmat,     only: eps,head,epsw1,epsw2 

! !LOCAL VARIABLES:

      implicit none
      integer(4),intent(in) :: iq                 !! index for the irreducible of q-point 
      integer(4),intent(in) :: iomfirst,iomlast   !! the range of freq mesh points to be calculated 
      
      integer(4) :: info
      integer(4) :: iom  ! Counter, runs over frequencies
      integer(4) :: im   ! index for mixed basis 
      integer(4), allocatable :: ipiv(:)
      integer(4) :: msiz
      integer :: lwork
      
      
      real(8) :: tstart,tend
      real(8) :: ecwq
     
      complex(8) :: tr,det                        
      complex(8), allocatable :: work(:)
      complex(8), pointer :: eps0(:,:)
      logical:: ldbg=.true.
      character(len=10)::sname="calcecrpa"
      
       
! 
! !EXTERNAL ROUTINES: 
!

  
      complex(8), external :: zdotu

      external zgetrf

       
!
! !INTRINSIC ROUTINES: 
!



! !REVISION HISTORY:
!
! Created 31.01.2007 by JH
!
!EOP
!BOC

      msiz=matsiz
      if(iq.eq.1) then 
        msiz=matsiz+1
        allocate(eps0(msiz,msiz),stat=info)
        call errmsg(info.ne.0,sname,"Fail to allocate eps0")
      endif 

      lwork=msiz*64
      allocate(ipiv(msiz),work(lwork),stat=info)
      call errmsg(info.ne.0,sname,"Fail to allocate ipiv, work")

      do iom=iomfirst,iomlast
        if(iq.eq.1) then 
          eps0(1,1)=head(iom)
          eps0(2:msiz,1)=epsw1(1:matsiz,iom)
          do im=1,matsiz
            eps0(im+1:msiz,im+1)=eps(im:matsiz,im,iom)
          enddo 
        else 
          eps0 => eps(:,:,iom)
        endif 
!
!       Calculate the trace of $ 1 - \epsilon(q,iu)  $
!
        tr=czero
        do im=1,msiz
          tr = tr + (cone - eps0(im,im))
        enddo
!
!       Calculate the determinant of $ \epsilon(q,iu) $
!
        call zhetrf('l',msiz,eps0,msiz,ipiv,work,lwork,info)

        call errmsg0(info,sname,"Fail to call zgetrf")

        det = cone
        do im=1,msiz
          det = det*eps0(im,im) 
        enddo
        ecwq=(log(abs(det))+real(tr))/twopi

        if(ldbg) write(6,100) omega(iom),womeg(iom),real(tr),real(det),ecwq
        ec_acfd=ec_acfd+ecwq*womeg(iom)*kwt_bz(iq)
      enddo 

      deallocate(ipiv,work)
      if(iq.eq.1) then 
        deallocate(eps0)
      else 
        nullify(eps0)
      endif 

 100  format("omega=",f8.3,1x,"weight=",f8.3,1x,  &
     &       "Tr(1-\epsilon)=",f8.3,1x,&
     &       "Det(\epsilon)= ",f12.3,1x,&
     &       "Ec(iom,q)=",f12.4)
     
      end subroutine calcecrpa
!EOC      
