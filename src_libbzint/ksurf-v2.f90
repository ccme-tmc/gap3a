!BOP
!
! !ROUTINE: ksurf
!
! !INTERFACE: 
      subroutine ksurf(ei,esurf,weight)
!     
! !DESCRIPTION:
!
! This subroutine calculates the integration over an energy surface in the
! k-mesh, the area of the Fermi surface inside the tetrahedron is calculated
! which is the essential item to decide the wtmp of the integration over
! each vertex. 
!

! !USES:
 
! !INPUT PARAMETERS:
      use  order, only: sort
      implicit none

      real(8), intent(in) :: ei(4)   ! the band energies or energy difference at k
      real(8), intent(in) :: esurf   ! the energy or energy differnce surface

! !OUTPUT PARAMETERS:
      real(8), intent(out) :: weight(4) ! the wtmp at each corner
       
! !LOCAL VARIABLES:
     
      integer :: i, j
      integer :: sort_ind(4),indx ! Number of nodes below esurf
            
      real(8) :: d21, d31, d41, d32, d42, d43  ! energy differences       
      real(8) :: ds1, ds2, ds3, ds4 ! esurf-eo(i) 
      real(8) :: wt         ! total wtmp
      real(8) :: wtmp(4),eo(4),zero_eps=1.e-6
     

! !REVISION HISTORY:
! 
! Created 17th. Jan 2005 by XZL
! Last version: May 2007 by RGA
!EOP
!BOC
      wtmp(1:4)=0.0d0
      eo = ei
      call sort(4,eo,sort_ind) 
!
!     check how the Fermi energy cross the tetrahedron
!
      if(esurf.lt.eo(1))then
        indx=0
      else
        indx=1
        do i=2,4
          if(eo(i).le.esurf)indx=indx+1
        enddo  
      endif

!
!     calculate the surface integration in different cases
!
      select case(indx)
        case(0,4) ! no intersection  
          continue
        case(1) ! Only eo(1)< esurf, one triangle

          ds1=esurf-eo(1) 
          d21=eo(2)-eo(1) 
          d31=eo(3)-eo(1) 
          d41=eo(4)-eo(1) 
          
          !! total wtmp
          wt= ds1*ds1/(2.d0*d21*d31*d41)
          wtmp(2)=wt*ds1/(3.d0*d21)
          wtmp(3)=wt*ds1/(3.d0*d31) 
          wtmp(4)=wt*ds1/(3.d0*d41)
          wtmp(1)=wt-wtmp(2)-wtmp(3)-wtmp(4)

        case(2) ! eo(1)<eo(2)<esurf, two triangles
          ds1=esurf-eo(1) 
          ds2=esurf-eo(2) 
          ds3=eo(3)-esurf 
          ds4=eo(4)-esurf 
          d21=eo(2)-eo(1) 
          d31=eo(3)-eo(1) 
          d41=eo(4)-eo(1) 
          d32=eo(3)-eo(2) 
          d42=eo(4)-eo(2) 
          d43=eo(4)-eo(3) 
          
! total wtmp          
          wt=(ds1*ds4*d32+ds2*ds3*d41)/(2.0*d31*d41*d32*d42)

! wtmp(1)
          if(d43.gt.zero_eps) then 
            wtmp(1) = ds4**3/(6.0*d42*d43*d41*d41) &
     &               -ds3**3/(6.0*d32*d43*d31*d31) 
            wtmp(2) = ds4**3/(6.0*d41*d43*d42*d42) &
     &               -ds3**3/(6.0*d31*d43*d32*d32) 
          else 
            wtmp(1) = (3*d31*d32*ds3**2 - d31*ds3**3 - 2*d32*ds3**3) &
     &         /(6.*d31**3*d32**2) &
     &        + ((3*d31**2*d32**2*ds3 - 3*d31**2*d32*ds3**2  &
     &        - 6*d31*d32**2*ds3**2 + d31**2*ds3**3 + 2*d31*d32*ds3**3 &
     &        + 3*d32**2*ds3**3)*d43)/(6.*d31**4*d32**3)

            wtmp(2) = (3*d32*d31*ds3**2 - d32*ds3**3 - 2*d31*ds3**3) &
     &         /(6.*d32**3*d31**2) &
     &        + ((3*d32**2*d31**2*ds3 - 3*d32**2*d31*ds3**2  &
     &        - 6*d32*d31**2*ds3**2 + d32**2*ds3**3 + 2*d32*d31*ds3**3 &
     &        + 3*d31**2*ds3**3)*d43)/(6.*d32**4*d32**3)

          endif

! wtmp(3)

          if(d21.gt.zero_eps) then 
            wtmp(3)=  ds1**3/(6.0*d21*d41*d31*d31) &
     &            - ds2**3/(6.0*d21*d42*d32*d32) 

            wtmp(4) = ds1**3/(6.0*d21*d31*d41*d41) &
     &            - ds2**3/(6.0*d21*d32*d42*d42)
          else
            wtmp(3) =  (3*d32*d42*ds2**2 - d32*ds2**3 &
     &               - 2*d42*ds2**3)/(6.*d32**3*d42**2) & 
     &        +  (d21*(3*d32**2*d42**2*ds2 - 3*d32**2*d42*ds2**2 &
     &        - 6*d32*d42**2*ds2**2 + d32**2*ds2**3 + 2*d32*d42*ds2**3 &
     &        + 3*d42**2*ds2**3))/(6.*d32**4*d42**3) 
           
            wtmp(4) =  (3*d42*d32*ds2**2 - d42*ds2**3 &
     &               - 2*d32*ds2**3)/(6.*d42**3*d32**2) &
     &        +  (d21*(3*d42**2*d32**2*ds2 - 3*d42**2*d32*ds2**2 &
     &        - 6*d42*d32**2*ds2**2 + d42**2*ds2**3 + 2*d42*d32*ds2**3 &
     &        + 3*d32**2*ds2**3))/(6.*d42**4*d32**3)
          endif 

!! consistence check 
          if(abs(wt-sum(wtmp)).gt.1.e-4) then 
            write(6,*) "WARNING: inconsistent total wtmp!"
          endif 

        case(3) ! ee(1), ee(2) and ee(3) < efer, one triangle.

          ds4=eo(4)-esurf 
          d41=eo(4)-eo(1) 
          d42=eo(4)-eo(2) 
          d43=eo(4)-eo(3) 
          
          !! total wtmp
          wt= ds4*ds4/(2.d0*d41*d42*d43)
          wtmp(1)=wt*ds4/(3.d0*d41)
          wtmp(2)=wt*ds4/(3.d0*d42)
          wtmp(3)=wt*ds4/(3.d0*d43)
          wtmp(4)=wt-wtmp(1)-wtmp(2)-wtmp(3)

      end select
      do i=1,4
        weight(sort_ind(i))=wtmp(i)
      enddo

      end subroutine ksurf

!EOC        
