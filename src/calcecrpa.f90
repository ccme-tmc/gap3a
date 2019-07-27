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
      use task,        only: fid_outgw, fid_outdbg

! !LOCAL VARIABLES:

      implicit none
      integer(4),intent(in) :: iq                 !! index for the irreducible of q-point 
      integer(4),intent(in) :: iomfirst,iomlast   !! the range of freq mesh points to be calculated 
      
      integer(4) :: info
      integer(4) :: gap
      integer(4) :: iom  ! Counter, runs over frequencies
      integer(4) :: im   ! index for mixed basis 
      integer(4), allocatable :: ipiv(:)
      integer(4) :: msiz
      integer :: lwork
      
      
      real(8) :: tstart,tend
      real(8) :: ecwq0
      real(8) :: tr,det                        
      real(8), allocatable :: rwork(:), eps_eigen(:)
    
      complex(8) :: temp
      complex(8), allocatable :: work(:)
      !complex(8), pointer :: eps0(:,:)
      complex(8), allocatable :: eps0(:,:)
      logical :: ldbg=.true.
      logical :: done=.false.
      character(len=10)::sname="calcecrpa"
      
       
! 
! !EXTERNAL ROUTINES: 
!

  
      complex(8), external :: zdotu

      external zhetrf, zheev

       
!
! !INTRINSIC ROUTINES: 
!



! !REVISION HISTORY:
!
! Created 31.01.2007 by JH
!
!EOP
!BOC

      ! for iq!=1, the head is set to one such that it has no effect
      ! on the RPA correlation
      msiz=matsiz+1
      allocate(eps0(msiz,msiz),stat=info)
      call errmsg(info.ne.0,sname,"Fail to allocate eps0")

      !lwork=msiz*64 ! for zhetrf NB
      lwork=msiz*65 ! for zheev: NB+1
      allocate(ipiv(msiz),work(lwork),eps_eigen(msiz),rwork(3*msiz),stat=info)
      call errmsg(info.ne.0,sname,"Fail to allocate lapack work arrays")

      !if(ldbg) then
      !  write(*,*) "#msiz = ", msiz, "; matsiz = ", matsiz
      !end if
      do iom=iomfirst,iomlast
        !if(ldbg) then
        !  write(*, "(A10,I5)") "#iom=", iom
        !endif 
        eps0 = czero
        do im=1,msiz
          eps0(im,im)=cone
        enddo 

        if(iq.eq.1) then
          !if(ldbg) then
          !  write(*,*) "#eps(1,1) = ", eps(1,1,iom)
          !endif
          eps0(1,1)=head(iom)
          eps0(2:msiz,1)=epsw1(1:matsiz,iom)
          eps0(1,2:msiz)=epsw2(1:matsiz,iom)
        endif

        do im=1,matsiz
          eps0(im+1:msiz,im+1)=eps(im:matsiz,im,iom)
          !eps0(2:msiz,im+1)=eps(1:matsiz,im,iom)
        enddo 
        !tr=0.0D0
        !do im=1,msiz
        !  !tr = tr + (cone - eps_eigen(im))
        !  !tr = tr + (cone - eps0(im,im))
        !  tr = tr + (1.0D0 - real(eps0(im,im),8))
        !enddo
!
!       diagonalize the dielectric matrix
!
        !call zhetrf('l',msiz,eps0,msiz,ipiv,work,lwork,info)
        call zheev('N','L',msiz,eps0,msiz,eps_eigen,work,lwork,rwork,info)
        ! Intel 2018.0 with -O3 may confront bug with zheev. Use 2017.1 or 2018.1 instead
        !call zgetrf(msiz,msiz,eps0,msiz,ipiv,info)
        call errmsg0(info,sname,"Fail to call zhetrf")
        !! sort the eigenvalues by its real value
        !gap=msiz/2
        !do while(gap.ge.1)
        !  done=.false.
        !  do while(.not.done)
        !    done=.true.
        !    do im=1,msiz-gap
        !      if (real(eps_eigen(im)).lt.real(eps_eigen(im+gap))) then
        !        temp=eps_eigen(im)
        !        eps_eigen(im)=eps_eigen(im+gap)
        !        eps_eigen(im+gap)=temp
        !        done=.false.
        !      endif
        !    enddo
        !  enddo
        !  gap=gap/2
        !enddo

        if(ldbg) then
          write(fid_outdbg,1000) womeg(iom), kwt_bz(iq)
          do im=1,msiz
            write(fid_outdbg,"(2I5,es26.17)") iq, im, eps_eigen(im) - 1.0D0
          end do
        end if
!
!       Calculate the trace of $ 1 - \epsilon(q,iu)  $
!
        tr=0.0D0
        do im=1,msiz
          !tr = tr + (cone - eps_eigen(im))
          !tr = tr + (cone - eps0(im,im))
          !tr = tr + (1.0D0 - real(eps0(im,im),8))
          tr = tr + (1.0D0 - eps_eigen(im))
        enddo
!
!       Calculate the determinant of $ \epsilon(q,iu) $
!
        !det = cone
        det = 1.0D0
        do im=1,msiz
          det = det*eps_eigen(im)
          !det = det*eps0(im,im)
        enddo
        !ecwq0=(log(abs(det))+real(tr))/twopi
        ecwq0=(log(det)+tr)/twopi

        if(ldbg)then
          !write(6,100) omega(iom),womeg(iom),real(tr),real(det),ecwq0
          write(fid_outgw,100) omega(iom),womeg(iom),tr,det,ecwq0
          !write(6,101) aimag(tr),aimag(det)
        endif
        ec_acfd=ec_acfd+ecwq0*womeg(iom)*kwt_bz(iq)
      enddo 

      deallocate(ipiv,work,eps_eigen,rwork)
      deallocate(eps0)
      !if(iq.eq.1) then 
      !  deallocate(eps0)
      !else 
      !  nullify(eps0)
      !endif 

 100  format("omega=",f8.3,1x,"weight=",f8.3,1x,  &
     &       "Tr(1-\epsilon)=",f15.5,1x,&
     &       "Det(\epsilon)= ",f15.5,1x,&
     &       "Ec(iom,q)=",f12.7)
 101  format(46x,f15.5,16x,f15.5)
1000  format("#Frequency weight = ", es18.10, " kpoint weight = ", es12.6)
     
     
      end subroutine calcecrpa
!EOC      
