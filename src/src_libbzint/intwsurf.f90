!BOP
!
! !ROUTINE: intwsurf
!
! !INTERFACE:
       subroutine intwsurf(efer,iweight)
 
! !DESCRIPTION:
!
!   This subroutine calculates the surface integration weight of each k-point for
!   each band

! !USES:
       
       use order,          only: sort
       use tetra_internal, only: ntet,eband,tetweig,tetcorn,nirkp, &
     &                           nband,vt,fout
       implicit none     

! !INPUT PARAMETERS:
       real(8), intent(in)  :: efer                 ! fermi energy

! !OUTPUT PARAMETERS:
       real(8), intent(out) :: iweight(nband,nirkp) ! the weight of each k-point for each band

! !LOCAL VARIABLES:
       integer(4) :: it,i,ib,kin
       integer(4), dimension(4) :: ik
       real(8) :: term
       real(8), dimension(4) :: ee
       real(8), dimension(4) :: w1t

! !EXTERNAL ROUTINES:
       external  intweight1t
! 
!EOP
!BOC
      iweight = 0.d0 

      do it=1,ntet
        do ib=1,nband 
          do i=1,4
            ee(i)=eband(ib,tetcorn(i,it))
          enddo
          call sort(4,ee,ik)
          w1t(1:4)=0.0d0
          if((ee(1).lt.efer).and.(efer.lt.ee(4))) then
            call ksurf(ee,efer,w1t)     
            do i=1,4
              term=w1t(i)*tetweig(it)
              kin=tetcorn(ik(i),it)
              iweight(ib,kin)=iweight(ib,kin)+term*6.0d0*vt
            enddo
          endif
        enddo
      enddo

      return
      
      end subroutine intwsurf
      
!EOC
