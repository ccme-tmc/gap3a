!BOP
!
! !ROUTINE: trans_vector
!
! !INTERFACE:
      subroutine trans_vector(iop,ib0,ib1,tmat)

! !DESCRIPTION:
!
! This subroutine makes a unitary transform to the wave functions in a given range to generate 
! a set of new wave functions. This is used for self-consistent GW  
!
! !USES:

      use eigenvec,    only: zzk,zzq 
      use recipvec,    only: maxngk

      implicit none 
! !INPUT PARAMETERS:
      integer,intent(in)     :: iop      !! 1/2 to choose zzk/zzq 
      integer,intent(in)     :: ib0, ib1 !! only bands in ib..nb are transformed 
      complex(8),intent(in)  :: tmat(ib0:ib1,ib0:ib1)     !! transformation matrix (ib:nb, ib:nb) 
!EOP
!BOC
! !LOCAL VARIABLES:
 
      integer:: ierr ! error code 
      integer:: nbs   ! total band numbers 
      complex(8):: cone=1.d0, czero=0.d0
      complex(8), allocatable:: ztmp(:,:),ctmp(:,:) 
      character(len=20):: sname="trans_vector"
!
! !REVISION HISTORY:
!
! Created: Jan. 22, 2010 by H. Jiang
!
      nbs=ib1-ib0+1
      allocate(ztmp(maxngk,ib0:ib1),stat=ierr)
      call errmsg0(ierr,sname,"Fail to allocate ztmp")
      if(iop.eq.1) then 
        ztmp=zzk(:,ib0:ib1)
        call zgemm('n','n',maxngk,nbs,nbs,cone,ztmp,maxngk,&
     &             tmat,nbs,czero,zzk(:,ib0:ib1),maxngk)
      else 
        allocate(ctmp(ib0:ib1,ib0:ib1),stat=ierr)
        call errmsg0(ierr,sname,"Fail to allocate ctmp")
        ctmp=conjg(tmat)
        ztmp=zzq(:,ib0:ib1) 
        call zgemm('n','n',maxngk,nbs,nbs,cone,ztmp,maxngk,&
     &             ctmp,nbs,czero,zzq(:,ib0:ib1),maxngk)

        deallocate(ctmp)
      endif 
      deallocate(ztmp)
      end subroutine 
!EOC


