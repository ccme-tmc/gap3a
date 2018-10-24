!BOP
!
! !ROUTINE: genvectorf
!
! !INTERFACE:
      subroutine writminm(ik)

! !DESCRIPTION:
!
! This subroutine write minm matrix to direct-access files 
!
! !USES:


      use bands,     only: nomax,numin,nbmaxpol,nbandsgw,     &
     &                     nbmax
      use kpoints,   only: kpirind,idikp
      use minmmat
      
      implicit none 
      integer(4),intent(in) :: ik

      integer(4) :: ikp,irkp


      if(lmicm) write(cmunit,rec=ik) micm
      if(lminc) write(ncunit,rec=ik) minc
      if(lmicc) write(ccunit,rec=ik) micc

      if(lirkp ) then 
        write(nmunit,rec=ik) minm
      else 
        ikp= ik
        irkp = kpirind(ikp)

        if(.not.lminmir) then
          write(nmunit,rec=ikp) minm
        else
          write(nmunit,rec=ikp)  minm(1:matsiz,1:nomax,            &
     &                             numin:nbmaxpol)
          if(idikp(irkp).eq.ikp) then  
            write(nmirunit,rec=irkp) minm
          endif 
        endif 
      endif 
      end subroutine 
