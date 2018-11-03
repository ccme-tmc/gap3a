!BOP
!
! !ROUTINE: setkk0index
!
! !INTERFACE:
      subroutine setkk0index
!
! !DESCRIPTION:
!
! Assigns to each k-point in klist2 the closest k-point in klist
!
! !USES:

      use kmeshintp
      use kpoints
      use struk
      use task,   only: fid_outdbg
      
! !LOCAL VARIABLES:
      implicit none
  
      integer(4) :: i1,i2,i3  
      integer(4) :: ig(3)
      integer(4) :: ik0
      integer(4) :: ik2
      integer(4) :: ik2vec(3)
      integer(4) :: ik0vec(3)
      integer(4) :: tmp_index
      
      real(8)    :: gvec(3)
      real(8)    :: k2vec(3)
      real(8)    :: k0vec(3)
      real(8)    :: difk(3)
      real(8)    :: length_difk
      real(8)    :: minlength

!EOP
!BOC
      allocate(kk0ind(4,nkp2))
      do ik2=1, nkp2
        ik2vec=klist2(1:3,ik2)
        call k2cart(ik2vec,idvk2,k2vec)
        minlength=1.0d+5
        do ik0=1, nkp
          ik0vec=klist(1:3,ik0)
          call k2cart(ik0vec,idvk,k0vec)
          do i1=-1,1
            do i2=-1,1
              do i3=-1,1
                gvec(1:3)=dble(i1)*br2(1:3,1)+dble(i2)*br2(1:3,2)+      &
     &                    dble(i3)*br2(1:3,3)
                difk(1:3)=k2vec(1:3)-k0vec(1:3)-gvec(1:3)
                length_difk=difk(1)*difk(1)+difk(2)*difk(2)+difk(3)*difk(3)
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
              
        kk0ind(1,ik2)=tmp_index
        kk0ind(2:4,ik2)=ig(1:3)
        ik0vec=klist(1:3,tmp_index)
        write(fid_outdbg,'(2i6,3i4,2i6,3i4,2x,3i3)')ik2,ik2vec,idvk2,tmp_index, &
     &           ik0vec,idvk,ig
        call k2cart(ik0vec,idvk,k0vec)
        write(fid_outdbg,'(3f12.6,4x,3f12.6,2x,f18.10)')k2vec,k0vec,minlength
              
      enddo
      close(999)
      end subroutine setkk0index
!EOC       
            
     
