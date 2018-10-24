!BOP
!
! !SUBROUTINE: sorteq
!
! !INTERFACE:
      subroutine sorteq(v,idx,sigeq)
      
! !DESCRIPTION:
!
! This subroutine sorts the values of the vector v (integer) in an
! order that the first one is the one with the biggest possibility of 
! being equal to the others, then the second one, third one and the 
! forth one. Optionally, it also returns an integer vector that gives 
! the old idxes in vector v according to the new ones
!
! !INPUT PARAMETERS:
      use tetra_internal, only: ztol_sorteq
      implicit none
            
      real(8), intent(inout) :: v(4)
      integer(4), intent(out) :: idx(4)  ! the relation between the sorted and unsorted arrays
      integer(4), intent(out) :: sigeq    

! !LOCAL VARIABLES:

      integer(4) :: i,j,itmp  
      integer(4) :: sig(5)     ! sig(1:4) number of items in the array equal to
                               !  this specific item. For example, if all are 
                               !  different sig(1:4)=1. sig(5) is the sum 
      real(8)    :: vtmp
      
! !SYSTEM ROUTINES:

      intrinsic abs
            
! !REVISION HISTORY:
!
! Created 03rd. Nov. 2004 by XZL
!
!EOP
!BOC
      
      sig=0 
      do i=1,4
        do j=1,4
          if(abs(v(i)-v(j)).le.ztol_sorteq) sig(i)=sig(i)+1
        enddo
      enddo
      
      sig(5)=sig(1)+sig(2)+sig(3)+sig(4)

      do i=1,4
        idx(i)=i
      enddo

      select case (sig(5))
      case (16)
        continue
      case (10)                  ! three of them are equal
        do i=1,3                  ! we want to make it as a=b=c while d not
          if(sig(i).eq.1) then
            vtmp=v(i)                            ! we make sig(4)=1
            v(i)=v(4)
            v(4)=vtmp
            itmp=idx(i)
            idx(i)=idx(4)
            idx(4)=itmp
            itmp=sig(4)
            sig(4)=sig(i)
            sig(i)=itmp
          endif
        enddo
      case (8)                   ! doubly equal but not all
        do i=3,4                                 ! make the first two equal
          if(abs(v(i)-v(1)).lt.(abs(v(2)-v(1))))then
            vtmp=v(i)
            v(i)=v(2)
            v(2)=vtmp
            itmp=idx(i)
            idx(i)=idx(2)
            idx(2)=itmp
            itmp=sig(i)
            sig(i)=sig(2)
            sig(2)=itmp
          endif
        enddo
      case (6)                   ! only two of them are equal
        j=1
        do i=1,4                         ! make the first one with sig(1)=2
          if((sig(i).eq.2)) then
            vtmp=v(j)
            v(j)=v(i)
            v(i)=vtmp
            itmp=idx(j)
            idx(j)=idx(i)
            idx(i)=itmp
            itmp=sig(j)
            sig(j)=sig(i)
            sig(i)=itmp
            j=j+1
           endif
        enddo
      case (4)      ! all different, nothing to do
        continue
      case default   ! None of the above... ERROR
         write(6,*)'ERROR in sorteq: case not found'
         write(6,*)'sig(i)=', sig
         write(6,*)'v = ',v
         stop "ERROR in sorteq"
      end select
      sigeq = sig(5)
      
      end subroutine sorteq
!EOC
 
