!BOP
! 
! !ROUTINE: kip_qpeintp

!
! !INTERFACE:
      subroutine kip_qpeintp 

! !DESCRIPTION:

! This subroutine interpolate the qp energy correction from a sparse 
! mesh to a dense mesh and add it to bande

! !USES:

      use constants, only: hev, pi
      use struk,     only: nat,rbas,alat,ortho
      use kmeshintp, only: nbands1,eks1,eqp1,kvecs1,nkp1,     &
     &                     nbands2,eks2,eqp2,kvecs2,nkp2,     &
     &                     nvbm,ncbm,ib0_kip,ib1_kip,nsp_kip,      &
     &                     eferks2,eferqp2

!
! !LOCAL VARIABLES:
      implicit none

      integer(4) :: ib,ib0,ib1,nbintp 
      integer(4) :: ik
      integer(4) :: isp
      real(8):: deval,decond
      complex(8),allocatable:: dek1(:,:),dek2(:,:)
      logical::  lprt=.true.
      character(13), parameter :: sname = 'kip_qpeintp'
!EOP
!BOC
      if(lprt) write(6,*) " - readeqp"
      ib0=ib0_kip
      ib1=min(ib1_kip,nbands2)
 
      nbintp=ib1-ib0+1
      write(6,'(a,3i5)') ' ib0,ib1,nbintp=',ib0,ib1,nbintp
      write(6,'(a,3i5)') ' nkp1,nkp2',nkp1,nkp2

      allocate(dek1(nkp1,ib0:ib1),dek2(nkp2,ib0:ib1))
!
! Transform klist band to internal coordinates 
!
      if(ortho)then
        call cart2int(nkp1,rbas,alat,kvecs1)
        call cart2int(nkp2,rbas,alat,kvecs2)
      endif

      do isp=1,nsp_kip 
        do ik=1,nkp1
          do ib=ib0,ib1
            dek1(ik,ib) =cmplx(eqp1(ib,ik,isp)-eks1(ib,ik,isp),0.d0,8) 
          enddo 
        enddo 

        call fourintp(dek1,nkp1,kvecs1,dek2,nkp2,kvecs2,nbintp)

        do ik=1,nkp2
          do ib=ib0,ib1
            eqp2(ib,ik,isp)=eks2(ib,ik,isp)+real(dek2(ik,ib))
          enddo 
        enddo 

        deval=sum(real(dek2(:,ib0)))/nkp2
        decond=sum(real(dek2(:,ib1)))/nkp2
        write(6,100) " shift VBs not calculated in GW by",deval*hev
        write(6,100) " shift CBs not calculated in GW by",decond*hev

        eqp2(1:ib0-1,:,isp)= eks2(1:ib0-1,:,isp)+deval          
        eqp2(ib1+1:nbands2,:,isp) = eks2(ib1+1:nbands2,:,isp)+decond 
      enddo

      eferks2=(maxval(eks2(nvbm(1),:,1))+minval(eks2(ncbm(1),:,1)))/2.d0
      eferqp2=(maxval(eqp2(nvbm(1),:,1))+minval(eqp2(ncbm(1),:,1)))/2.d0
      call bandanaly(ib0,ib1,nkp2,kvecs2,eks2(ib0:ib1,:,:),&
     &               eferks2,nsp_kip,"KS")

      call boxmsg(6,'-',"QSGW(interpolated) Bands ")
      call bandanaly(ib0,ib1,nkp2,kvecs2,eqp2(ib0:ib1,:,:), &
     &               eferqp2,nsp_kip,"QSGW(intpl)")

      deallocate(dek1,dek2) 
      return
 100  format(a,f10.3,' eV')
      end subroutine 

!EOC
