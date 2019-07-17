!BOP
!
! !ROUTINE: coul_barc
!
! !INTERFACE:
      subroutine coul_barc(iq)

! !DESCRIPTION:
!
! This subroutine calculates the matrix of the bare coulomb potential according to equations
! \ref{coulsphdef},\ref{deficouloc}, \ref{coulijdefin} and
! \ref{ipw-coul-sphere}.
!
! !USES:

      use barcoul,    only: vmat,ev,iop_coulvm
      use mixbasis,   only: mbsiz,locmatsiz
      use recipvec,   only: ngq
      use task,       only: fid_outdbg,time_coul
!
! !INPUT PARAMETERS: 

      implicit none

      integer, intent(in) :: iq  ! index of the q-point

!
! !LOCAL VARIABLES:
      integer :: im,jm,iipw,ipw
      integer :: lwork,rwsize
      real(8),    allocatable :: rwork(:)   ! workspace array for the diagonalization subroutine
      complex(8), allocatable :: work(:)    ! workspace array for the diagonalization subroutine

!! debug 
      integer :: ierr 
      logical :: lprt =.false.
      character(len=20)::sname='coul_barc'
      real(8) :: tstart,tend        !! variables related to cputime
  
 
! !EXTERNAL ROUTINES: 

      
      real(8), external :: getcgcoef
      complex(8), external :: getdjmm
      
      external calcjlam
      external k2cart
      external rotate
      external ylm
      external zgemm
!
! !REVISION HISTORY:
! 
! Created 16th. March 2004 by RGA
! Last modified 31. March 2005 by RGA
!
!EOP
!BOC
      call cpu_time(tstart)
      if(lprt) call linmsg(6,'-','coul_barc')
!
!     calculates the matrix between PW's and orthogonalized IPW's. 
!
      call coul_mpwipw(iq)  
      
!
!     calculate bare coulomb potential matrix via plane waves expansion 
!
      if(iop_coulvm.eq.1) then 
        call coul_setvm1(0,iq)
      else
        call coul_setvm0(iq) 
      endif 

      if(lprt) then 
        write(fid_outdbg,*) "### barc from MT for iq= ",iq
        write(fid_outdbg,'(2a5,a24)') 'i','j ','barc(i,j)'
        write(fid_outdbg,*) "MT-MT"
        do jm=1,locmatsiz,locmatsiz/10
          do im=1,locmatsiz,locmatsiz/10
            write(fid_outdbg,'(2i5,10f12.6)') im,jm,vmat(im,jm)
          enddo
        enddo

        write(fid_outdbg,*) "MT-IS"
        do jm=1,locmatsiz,locmatsiz/10
          do iipw=1,ngq(iq),ngq(iq)/10
            im=iipw+locmatsiz
            write(fid_outdbg,'(2i5,10f12.6)') im,jm,vmat(im,jm)
          enddo
        enddo

        write(fid_outdbg,*) "IS-IS"
        do ipw=1,ngq(iq),ngq(iq)/10
          do iipw=1,ngq(iq),ngq(iq)/10
            im=iipw+locmatsiz
            jm=ipw +locmatsiz
            write(fid_outdbg,'(2i5,10f12.6)') im,jm,vmat(im,jm)
          enddo
        enddo
      endif 
!
!     Diagonalize the bare coulomb matrix
!
      if(lprt) write(6,*) "calculate square root of barc "
      lwork=2*mbsiz
      rwsize=3*mbsiz
      allocate(work(lwork),rwork(rwsize),stat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to allocate workspace !!! ")

      call zheev('v','u',mbsiz,vmat,mbsiz,ev,work,lwork,rwork,ierr)
      call errmsg(ierr.ne.0,sname,'Fail to diag. barc by zheev')
      deallocate(work,rwork)

!#ifdef DEBUG
       write(fid_outdbg,*) "### barc eigenvalues ### ", mbsiz
       do im=1,mbsiz
         write(fid_outdbg,'(i5,e16.6)') im,ev(im)
       enddo
!#endif
      call cpu_time(tend)
      time_coul = time_coul + tend - tstart
      end subroutine  
