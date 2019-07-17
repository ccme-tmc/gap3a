!BOP
!
! !ROUTINE: setkpg
!
! !INTERFACE:
      subroutine setkpg(ik,isq)
      
! !USES:

      use bands,    only: nv
      use kpoints,  only: idvk, klist, kpirind
      use lapwlo,   only: nlo_tot
      use recipvec, only: gindex, indgk, kpg, kppg, rk, npw
      
      implicit none
      
! !INPUT PARAMETERS: 
      
      integer(4), intent(in) :: ik ! Index of the k-point
      integer(4), intent(in) :: isq ! =1 for k, =2 for k'
      
! !LOCAL VARIABLES:
      integer(4) :: ig,i,irk,nvmax
      integer(4), dimension(3) :: ikvec
      real(8), dimension(3) :: kvec
      real(8), dimension(3) :: gvec
      external k2cart 
! !REVISION HISTORY:
!
! Created: 25th. March 2004 by RGA
!
!EOP
!BOC
      irk=kpirind(ik)
      ikvec(1:3)=klist(1:3,ik)

      call k2cart(ikvec,idvk,kvec)
      nvmax=maxval(nv)

      select case (isq)
        case(1)
          kpg(1:3,1:nvmax)=0.0d0
          rk(1:nvmax)=0.0d0

          do ig=1,nv(irk)+nlo_tot
             ikvec(1:3)=gindex(:,indgk(ig,ik))
             call k2cart(ikvec,1,gvec)
             kpg(:,ig) = kvec + gvec
             rk(ig) = sqrt( sum(kpg(1:3,ig)**2) )
          enddo

        case(2)
          kppg(1:3,1:nvmax)=0.0d0
          rk(1:nvmax)=0.0d0
          do ig=1,nv(irk)+nlo_tot
            ikvec(1:3)=gindex(:,indgk(ig,ik))
            call k2cart(ikvec,1,gvec)
            kppg(:,ig) = kvec + gvec
            rk(ig) = sqrt( sum(kppg(1:3,ig)**2) )
          enddo

      end select  
!   10 format(i4,4f12.6)      
      
      end subroutine setkpg
!EOC
