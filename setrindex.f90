!BOP
!
! !ROUTINE: setrindex
!
! !INTERFACE:
      subroutine setrindex

! !DESCRIPTION:
!
!Sets the indexes of the real space lattice vectors used for  Fourier interpolation 
! ordered by increasing length
!
! !USES:

      use constants,  only: pi
      use fouri,      only: rindex,nrr,rst,nst,rmax
      use struk,      only: ortho,rbas,nsym,izmat
      

! !LOCAL VARIABLES:

      implicit none     
      
      integer(4) :: gap       
      integer(4) :: i         ! (Counter): runs over coordinates
      integer(4) :: ippw      ! (Counter): runs over plane waves
      integer(4) :: ir1       ! (Counter): run over x-coord of G
      integer(4) :: ir2       ! (Counter): run over y-coord of G 
      integer(4) :: ir3       ! (Counter): run over z-coord of G 
      integer(4) :: jppw      ! (Counter): runs over plane waves
      integer(4) :: nr1       ! Max. ir1
      integer(4) :: nr2       ! Max. ir2
      integer(4) :: nr3       ! Max. ir3
      integer(4) :: nr       ! Maximum number of plane waves in intipw
      integer(4) :: nrmin,nrmax,ist,j,isym

      integer(4), dimension(3) :: irvec ! Integer coordinates of the G-vector
      integer(4), allocatable :: invrindex(:,:,:)
      
      real(8) :: rr

      real(8), dimension(3) :: rvec     ! Cartesian coordinates of the R-vector
      integer(4), allocatable :: rind(:,:) ! Temporary storage for rindex
      real(8), allocatable :: rlen(:)   ! Temporary storage for all the  qpg's
      
      logical :: done
      logical :: lprt= .false.


      integer :: ierr
      character(len=10)::sname="setrindex"

! !EXTERNAL ROUTINES: 


      external r2cart

! !REVISION HISTORY:
!
! Created May 2004 by RGA
! Last Modified 9. Nov. 2005 by RGA
!
!EOP
!BOC

#ifdef DEBUG
      lprt=.true.
#endif

      do i=1,3
        rvec(i)=sqrt(sum(rbas(i,:)**2)) 
      enddo
      nr1=nint(rmax/rvec(1))*2
      nr2=nint(rmax/rvec(2))*2
      nr3=nint(rmax/rvec(3))*2
      nr = (2*nr1+1)*(2*nr2+1)*(2*nr3+1)

      allocate(rlen(nr),rind(3,nr),stat=ierr)
      call errmsg(ierr.ne.0,sname,"error in allocation") 

      if(lprt) then 
         write(6,*) '--------- R vectors generation ----------'
         write(6,*)'  Parameters'
         write(6,12) rmax,nr1,nr2,nr3,nr
      endif 
      
      ippw=0
      do ir1=-nr1,nr1
        irvec(1)=ir1
        do ir2=-nr2,nr2
          irvec(2)=ir2
          do ir3=-nr3,nr3
            irvec(3)=ir3
!
!           Transform irvec to cartesian coordinates
!
            do i=1,3
              rvec(i)=sum(dble(irvec(1:3))*rbas(1:3,i))
            enddo
            rr=sqrt(sum( rvec(:)**2) )
            if(rr.le.rmax )then
              ippw = ippw + 1
              rlen(ippw)=rr
              rind(1:3,ippw) = irvec(1:3)
            endif
          enddo
        enddo
      enddo
      nrr = ippw
      
!
!     sort by increasing length using shell algorithm
!
      call shelsort(nrr,rind,rlen) 

      allocate(rindex(1:3,1:nrr))
      allocate(rst(2,1:nrr))
      rst(1:2,1:nrr)=0
      rindex(1:3,1:nrr)=rind(1:3,1:nrr)

!
!    generate the inverse of rindex      
!
      nrmax=maxval(rindex)+1
      nrmin=minval(rindex)-1
      allocate(invrindex(nrmin:nrmax,nrmin:nrmax,nrmin:nrmax))
      invrindex(nrmin:nrmax,nrmin:nrmax,nrmin:nrmax)=0
      do ippw=1,nrr
        ir1=rindex(1,ippw)
        ir2=rindex(2,ippw)
        ir3=rindex(3,ippw)
        invrindex(ir1,ir2,ir3)=ippw
      enddo  
      rst(1,1)=1
      rst(2,1)=nsym
      ist=1
      do ippw=2,nrr
        if(rst(1,ippw).eq.0)then
          ist=ist+1
          rst(1,ippw)=ist
          do isym=1,nsym
            do i=1,3
              irvec(i)=0
              do j=1,3
                irvec(i)=irvec(i)+izmat(j,i,isym)*rindex(j,ippw)
              enddo ! j
            enddo ! i

            jppw=invrindex(irvec(1),irvec(2),irvec(3))
            if(jppw.le.0)then
              write(6,111) ippw,rindex(1:3,ippw),isym,irvec(1:3)
            else  
              rst(1,jppw)=ist
              rst(2,jppw)=rst(2,jppw)+1
            endif  
          enddo ! isym
        endif
      enddo ! ippw
      nst=ist      

      if(lprt) write(6,*) " nrr=",nrr
      if(lprt) write(6,*) " nst=",nst

      deallocate(invrindex,rlen,rind)

   10 format(i6,3i4,4f12.6,2i4)
  111 format(2(i6,3i4)) 
   12 format('rmax = ',f15.10,' nr1 =',i4,' nr3 =',i4,' nr3 =',i4, &
     &       ' nr =',i8)
 6110 format (2(3i2,/),3i2)
 6111 format(4x,i4)
 
      end subroutine setrindex
!EOC
