!BOP
! 
! !ROUTINE: covsrt
!
! !INTERFACE:
      subroutine covsrt(covar,npc,ma,ia,mfit)

! !DESCRIPTION:
!
! Expand in storage the covariance matrix \texttt{covar}, so as to take
! into account parameters that are being held fixed. (For the latter, return
! zero covariances.)
!
! !INPUT PARAMETERS:
      
      implicit none
      
      integer(4), intent(in) :: ma 
      integer(4), intent(in) :: mfit 
      integer(4), intent(in) :: npc 
      integer(4), intent(in) :: ia(1:ma)

! !INPUT/OUTPUT PARAMETERS:      

      complex(8), intent(inout) :: covar(1:npc,1:npc)

! !LOCAL VARIABLES:
      
      integer(4) :: i
      integer(4) :: j
      integer(4) :: k
      
      complex(8)    :: swap

! !REVISION HISTORY:
!
! Original subroutine: covsrt.for (c) copr. 1986-92 numerical recipes
! software &669i..
! Last modified: 7th. Jul. 2005 by RGA
!
!EOP
!BOC
      do i=mfit+1,ma
        do j=1,i
          covar(i,j)=0.
          covar(j,i)=0.
        enddo ! i
      enddo ! j
      k=mfit
      do j=ma,1,-1
        if(ia(j).ne.0)then
          do i=1,ma
            swap=covar(i,k)
            covar(i,k)=covar(i,j)
            covar(i,j)=swap
          enddo ! i
          do i=1,ma
            swap=covar(k,i)
            covar(k,i)=covar(j,i)
            covar(j,i)=swap
          enddo ! i
          k=k-1
        endif
      enddo ! j
      return
      end subroutine covsrt
!EOC


