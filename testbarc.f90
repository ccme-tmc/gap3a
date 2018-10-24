!BOP
!
! !ROUTINE: testbarc
!
! !INTERFACE:
      subroutine testbarc(iq)

! !DESCRIPTION:
!
! This subroutine calculates the matrix of the sqare root of the bare coulomb
! potential. The resulting matrix is stored back in barc.
!
! !USES:

      use constants,  only: czero,pi,cone
      use barcoul,    only: vmat
      use kpoints,    only: qlist,idvq
      use mixbasis,   only: nmixmax,matsiz,locmatsiz
      use recipvec,   only: gindex, indgq
      use struk,      only: vi
      use task,       only: fid_outdbg
!
! !INPUT PARAMETERS: 

      implicit none
      
      integer(4), intent(in) :: iq ! index of the q-point
!
! !LOCAL VARIABLES:

      integer(4) :: i,j  ! just some counters

      integer(4) :: info   ! output status of the diagonalization subroutine
      integer(4) :: lwork  ! size of the workspace array work
      integer(4) :: rwsize ! size of the workspace array rwork
      integer(4) :: gap
      
      integer(4), dimension(3) :: iqvec,igvec
      
      real(8) :: qg1len, kk, t1,t2, x
      real(8), dimension(3) :: qg1, qvec, gvec
      real(8), allocatable :: rwork(:)     ! workspace array for the diagonalization subroutine
      real(8), allocatable :: ev(:),ev2(:) ! eigenvalues of barc
      real(8), allocatable :: exev(:)      ! exact eigenvalues of barc

      complex(8), allocatable :: vtemp(:,:) ! temporary storage for the eigenvectors of barc
      complex(8), allocatable :: vsq(:,:) ! the square root of barc
      complex(8), allocatable :: work(:)  ! workspace array for  the diagonalization subroutine
      complex(8), allocatable :: barc(:,:),sqbarc(:,:)
    
      character(len=67) :: errmsg !the error message
      logical :: done

! !DEFINED PARAMETERS:
 
        character(len=8) , parameter :: sname = 'testbarc'        
 
! !EXTERNAL ROUTINES: 

      external coul_barc
      external outerr
      external zheev
!
! !INTRINSIC ROUTINES: 
!
      intrinsic abs
      intrinsic conjg
      intrinsic sqrt
      intrinsic cpu_time
!
! !REVISION HISTORY:
! 
! Created 4th. Jan. 2005 by RGA
! Modified 14.08.2007 by JH
!
!EOP
!BOC
!
!     Set up the workspace for the diagonalization subroutine
!
      allocate(vtemp(matsiz,matsiz))
      allocate(barc(matsiz,matsiz))
      allocate(sqbarc(matsiz,matsiz))
      allocate(ev(matsiz))
      allocate(ev2(matsiz))
      allocate(exev(matsiz))
      lwork=2*matsiz
      rwsize=3*matsiz
      allocate(work(lwork))
      allocate(rwork(rwsize))
      allocate(vsq(matsiz,matsiz))

!
! check whether sqbarc \times sqbarc equal to barc
!
      call zgemm('n','n',matsiz,matsiz,matsiz,cone,sqbarc,matsiz,  &
     &              sqbarc,matsiz,czero,vsq,matsiz)   
      write(6,*) " Difference between sqbarc**2 and barc"
      write(6,'(2e12.5)') maxval(abs(vsq-barc)),  &
     &               sum(abs(vsq-barc))/(matsiz**2) 
      
!
!     Diagonalize the bare coulomb matrix
!
      vtemp(:,:)=barc(1:matsiz,1:matsiz)
      call zheev('v','u',matsiz,vtemp,matsiz,ev,work,lwork,rwork,info)
      call errmsg0(info,sname,"calling zheev")
      
      write(fid_outdbg,*) "### barc direct from mb ###" 
      do i=1,matsiz,matsiz/10
        write(fid_outdbg,*)'eigenvector:',i,'eval =',ev(i)
        do j=1,matsiz,matsiz/10
           write(fid_outdbg,11)j,vtemp(j,i)
        enddo
        write(fid_outdbg,*)
      enddo     

      vsq(:,:)=barc(:,:)
      call coul_setvm1(0,iq)
      vtemp(:,:)=vmat(1:matsiz,1:matsiz)
      call zheev('v','u',matsiz,vtemp,matsiz,ev2,work,lwork,rwork,info)
      call errmsg0(info,sname,"calling zheev (2)")

      write(6,*) " Difference between two different barc"
      write(6,'(2e12.5)') maxval(abs(vsq-barc)),  &
     &               sum(abs(vsq-barc))/(matsiz**2) 

      write(fid_outdbg,*) "### barc from pw and transform to mb ###" 
      do i=1,matsiz,matsiz/10
        write(fid_outdbg,*)'eigenvector:',i,'eval =',ev2(i)
        do j=1,matsiz,matsiz/10
           write(fid_outdbg,11)j,vtemp(j,i)
        enddo
        write(fid_outdbg,*)
      enddo
!
!       calculate the exact eigenvalues 
!
      if(iq.eq.1)then
        do i=1,matsiz
!
!           Calculate the G vector          
! 
          igvec(1:3)=gindex(:,indgq(matsiz-i+2,iq))
          call k2cart(igvec,1,gvec)
          qg1len = sum(gvec**2)
          exev(i)=4.0d0*pi/qg1len
        enddo    
      else  
        iqvec(1:3)=qlist(1:3,iq)
        call k2cart(iqvec,idvq,qvec)
        do i=1,matsiz

!           Calculate the G vector          
          igvec(1:3)=gindex(:,indgq(matsiz-i+1,iq))
          call k2cart(igvec,1,gvec)

!           Calculate q+G:
          qg1(1:3) = qvec(1:3) + gvec(1:3)

!           Calculate |q+G|^2:

          qg1len = sum(qg1**2) 
          exev(i)=4.0d0*pi/qg1len
        enddo    

! sorting 
        gap=matsiz/2
        do while(gap.ge.1)
          done=.false.
          do while(.not.done)
            done=.true.
            do i=1,matsiz-gap
              if(exev(i).gt.exev(i+gap))then
                kk=exev(i)
                exev(i)=exev(i+gap)
                exev(i+gap)=kk
                done=.false.
              endif
            enddo
          enddo
          gap=gap/2
        enddo
      endif 

      write(fid_outdbg,*) "### Comparison of barc eigenvalues ###"
      do i=1,matsiz,matsiz/20
        write(fid_outdbg,12)i,ev(i),ev2(i),exev(i)
      enddo  

      deallocate(work)
      deallocate(rwork)
      deallocate(vtemp)
      deallocate(ev)
      deallocate(ev2)
      deallocate(exev)
      deallocate(vsq)
    
   11 format(i4,' (',f18.10,',',f18.10,')')
   12 format(i4,3g18.10)
   13 format(2i4,3f18.10)
   14 format(2i4,2('(',g18.10,',',g18.10,')'))
   15 format(i4,'(',g18.10,',',g18.10,')')
   16 format(2i4,3f18.10,' err')
 1000 format('CPUTIME for ',a20,f16.2,' seconds')  

      return
         
      end subroutine testbarc
!EOC            
            
            
        
      
