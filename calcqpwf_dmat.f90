!BOP
!
! !ROUTINE: calcqpwf_dmat
!
! !INTERFACE:
      subroutine calcqpwf_dmat

! !DESCRIPTION:
!       This subroutine calculates the density matrix of QP wavefunctions 
!     using the KS wave functions as the basis. Only states in the range of 
!     ibgw..nbgw are taken into account. This quantity is used for the calculation 
!     if selfxmn, dvhmn in the self-consistent GW. To treat metallic systems 
!     the weight arising from BZ integration is absorbed in qpwf_dmat 
! !USES:
!
      use bands,      only: bande,nbandsgw,nspin,ibgw,nbgw,eferqp,nomaxs
      use bzinteg,    only: kiw
      use kpoints,    only: nkp,nirkp,wkir,idikp
      use selfenergy, only: qpwf_coef,qpwf_dmat
      implicit none 

      integer :: ie1,ie2,ie  !! index for bands 
      integer :: ik, irk !! index for k-points in BZ and IBZ, respectively 
      integer :: isp !! index for spin
      real(8) :: fnk
      complex(8):: tr

      logical :: lprt=.false.

#ifdef DEBUG
      lprt = .true.
#endif 

      qpwf_dmat=0.d0
      do isp=1,nspin
        do irk=1,nirkp
          ik=idikp(irk)
          do ie2=ibgw,nbgw 
            do ie1=ibgw,nbgw 
              do ie=ibgw,nomaxs(isp) 
                fnk = kiw(ie,irk,isp)*nkp 
                qpwf_dmat(ie1,ie2,irk,isp) = qpwf_dmat(ie1,ie2,irk,isp) &
     &             + fnk*qpwf_coef(ie1,ie,irk,isp)                    &
     &            *conjg(qpwf_coef(ie2,ie,irk,isp))
              enddo  ! ie
            enddo ! ie1
          enddo !! ie2

          if(lprt) then 
            write(6,*) "dmat for (irk,isp)=",irk,isp
            tr=0.d0
            do ie2=ibgw,nbgw 
              write(6,100) qpwf_dmat(:,ie2,irk,isp)
              tr=tr+qpwf_dmat(ie2,ie2,irk,isp)
            enddo 
            write(6,*) "Trace=",tr
            write(6,*)
          endif 
 100  format(100f6.2)
        enddo !! irk 
      enddo !! isp  
        
      end subroutine 
