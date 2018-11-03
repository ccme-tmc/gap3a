!BOP
! !ROUTINE: bz_setiksym
!
! !INTERFACE:
       subroutine bz_setiksym

! !DESCRIPTION:
!
! This subroutine set the index of the symmetry operation relating each
! k-point with the corresponding irreducible one.
!
! !USES:
      
      use constants, only: pi
      use kpoints,   only: nirkp, nkp, klist, kirlist, kpirind, iksym,  &
     &                     idvk, idikp, g0
      use struk,     only: alat,br2,ortho,lattic,nsym,imat,tau,izmat,   &
     &                     rbas,gbas


! !LOCAL VARIABLES:

      implicit none
      
      integer(4) :: i
      integer(4) :: j
      integer(4) :: ikp
      integer(4) :: isim
      integer(4) :: ikvec(3)
      integer(4) :: irkvec(3)
      integer(4) :: jkvec(3)
      integer(4) :: dg(3)
      integer(4) :: idivg0
      
      real(8) :: dtk
      real(8) :: rkvec(3)

      character(20)  :: sname="bz_setiksym"
      character(120) :: msg

      logical :: found

!! !REVISION HISTORY:
!       
! Created Nov. 2006 by RGA
!
!EOP
!BOC
!
      if(allocated(iksym))deallocate(iksym)
      allocate(iksym(nkp))
      if(allocated(idikp))deallocate(idikp)
      allocate(idikp(nirkp))
      if(allocated(g0))deallocate(g0)
      allocate(g0(1:3,nkp))
      g0=0
      do ikp=1,nkp
        ikvec(1:3)=klist(1:3,ikp)
        irkvec(1:3)=kirlist(1:3,kpirind(ikp))
        isim=1
        found=(ikvec(1).eq.irkvec(1)).and.(ikvec(2).eq.irkvec(2)).and.  &
     &          (ikvec(3).eq.irkvec(3))  
        if(found)then
          idikp(kpirind(ikp))=ikp
          iksym(ikp)=1          
        else   
          do while ((isim.le.nsym).and.(.not.found))
            do i=1,3
              jkvec(i)=mod(izmat(i,1,isim)*irkvec(1)+izmat(i,2,isim)*irkvec(2)&
     &                +izmat(i,3,isim)*irkvec(3),idvk)
              dg(i)=(1-isign(1,jkvec(i)))/2
              jkvec(i)=jkvec(i)+dg(i)*idvk
            enddo
            found=(ikvec(1).eq.jkvec(1)).and.(ikvec(2).eq.jkvec(2)).and.&
     &           (ikvec(3).eq.jkvec(3))  
            if(found)then
              iksym(ikp)=isim
              g0(1:3,ikp)=dg(1:3)
            else  
              isim=isim+1
            endif  
          enddo  ! while
          if(.not.found)then
            write(msg,*) 'iksym not found for kpoint nr.',ikp
            call outerr(sname,msg)
          endif  
        endif
        rkvec(1:3)=dble(klist(1:3,ikp))/dble(idvk)
        dtk=rkvec(1)*tau(1,iksym(ikp))+rkvec(2)*tau(2,iksym(ikp))+&
            rkvec(3)*tau(3,iksym(ikp))
      enddo ! ikp
      dg(1:3)=1
      if(ortho.or.(lattic(1:3).eq.'CXZ'))then
        idivg0 = 1
        call cartezian(nkp,alat,br2,g0,idivg0) 
      endif        

      return
      
      end subroutine bz_setiksym
!EOC          
