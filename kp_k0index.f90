!BOP
!
! !ROUTINE: kp_k0index
!
! !INTERFACE:
      subroutine kp_k0index(k,k0list,nk0,k0ind)
!
! !DESCRIPTION:
!
! Assigns to each k-point in klist2 the closest k-point in klist
!
! !USES:
      use struk,  only: br2
      implicit none
      integer,intent(in) :: nk0
      real(8),intent(in) :: k(3)
      real(8),intent(in) :: k0list(3,nk0)
      integer,intent(out):: k0ind(4)
  
! !LOCAL VARIABLES:
      integer(4) :: i1,i2,i3  
      integer(4) :: ig(3)
      integer(4) :: ik0
      integer(4) :: ik0vec(3)
      integer(4) :: tmp_index
      
      real(8)    :: gvec(3)
      real(8)    :: kvec(3)
      real(8)    :: k0(3),k0vec(3)
      real(8)    :: difk(3)
      real(8)    :: length_difk
      real(8)    :: minlength

!EOP
!BOC
      call kfrac2cart(k,kvec)
      minlength=1.0d+5
      do ik0=1, nk0
        k0=k0list(1:3,ik0)
        call kfrac2cart(k0,k0vec)
        do i1=-1,1
          do i2=-1,1
            do i3=-1,1
              gvec=dble(i1)*br2(:,1)+dble(i2)*br2(:,2)+dble(i3)*br2(:,3)
              difk(1:3)=kvec(1:3)-k0vec(1:3)-gvec(1:3)
              length_difk=sqrt(sum(difk(:)**2)) 
              if(length_difk.lt.minlength)then
                tmp_index=ik0
                minlength=length_difk
                ig(1)=i1
                ig(2)=i2
                ig(3)=i3
              endif
            enddo
          enddo
        enddo    
      enddo
              
      k0ind(1)=tmp_index
      k0ind(2:4)=ig(1:3)
      k0=k0list(1:3,tmp_index)
      call kfrac2cart(k0,k0vec)
      write(6,100) kvec,k0vec,minlength

  100 format("k=",3f8.4,"k_0=",3f8.4,"\Delta k=",f8.4) 
      end subroutine kp_k0index
!EOC       
            
     
