!BOP
!
! !ROUTINE: set_minm
!
! !INTERFACE:
      subroutine set_minm(iq,iop)

! !DESCRIPTION:
!
! This subroutine calculate all Minm elements and write them into the
! minm file.
! Input: 
!    iq         - index for q-point 
!    iop   - control how to transform minm
!       == 0  --  calculate Minm from KS vectors and save to the file 
!       == 1  --  update Minm in the file by transforming "m" by qpwf_coef
!       == 2  --  update Minm in the file by transforming both "m" and "n"
! !USES:

      use bands,       only: nbmax,nbmaxpol,nspin,    & 
     &                       nomaxs,numins,ibgw,nbgw
      use constants,   only: cone, czero
      use core,        only: ncg
      use kpoints,     only: nkp,kqid,kpirind,idikp
      use minmmat,     only: io_minm
      use mixbasis,    only: matsiz
      use selfenergy,  only: qpwf_coef

! !INPUT PARAMETERS:
      
      implicit none
      integer, intent(in) :: iq    ! Indices for k, q and spin  
      integer, intent(in) :: iop    

! !LOCAL VARIABLES:
      integer:: ik,irk  ! index for k  in BZ and IBZ, respectively 
      integer:: jk,jrk  ! index for k' in BZ and IBZ, respectively 
      integer:: ierr
      integer:: nvbm,ncbm
      integer:: isp
      logical:: ldbg=.false.

      real(8) :: time1,time2,tstart,tend
      character(len=15)::sname='set_minm'
      complex(8),allocatable :: minm(:,:,:) 

! !REVISION HISTORY:
!
! Created: Aug. 4, 2010 by Jiang, H
!EOP
!BOC
      if(ldbg) call linmsg(6,'-',sname)

      do isp=1,nspin 
        nvbm=nomaxs(isp) 
        ncbm=numins(isp)

        do ik=1,nkp
          irk=kpirind(ik)
          jk=kqid(ik,iq)
          jrk=kpirind(jk)

          if(ldbg) write(6,*) " - expand vectors"
          if(iop.eq.0) then 
            call readvector(ik,1,isp,iop)          !! Read the eigenvector corresponding to the k-point ik
            call expand_evec(ik,1,.true.,isp)  !! Calculate the expansion coeficients of the eigenvectors
            call readvector(jk,2,isp,iop)          !! Read the eigenvector corresponding to the k'-point jk
            call expand_evec(jk,2,.true.,isp)  !! Calculate the expansion coeficients of the eigenvectors
          endif 

          if(idikp(irk).eq.ik) then   ! for k in IBZ
            call sub_setminm('cm',1,ncg,ncbm,nbmaxpol)
            call sub_setminm('nc',ibgw,nbgw,1,ncg)
            call sub_setminm('nm',ibgw,nbgw,1,nbmax)
            if(ibgw.gt.1) then 
              call sub_setminm('nm',1,ibgw-1,ncbm,nbmaxpol)
            endif 
          else   !! for k not in IBZ
            call sub_setminm('cm',1,ncg,ncbm,nbmaxpol)
            call sub_setminm('nm',1,nvbm,ncbm,nbmaxpol)
          endif 
        enddo !ik
      enddo ! isp

      contains 
        subroutine sub_setminm(nmflag,nst,nend,mst,mend)
        character(2):: nmflag
        integer:: nst,nend,mst,mend

        complex(8),allocatable::minm(:,:,:)

        allocate(minm(matsiz,nst:nend,mst:mend))
        if(iop.eq.0) then 
          if(nmflag.eq.'nm') then
            call calcminm(ik,iq,nst,nend,mst,mend,isp,minm)
          elseif(nmflag.eq.'nc') then
            call calcminc(ik,iq,nst,nend,mst,mend,isp,minm)
          elseif(nmflag.eq.'cm') then
            call calcmicm(ik,iq,nst,nend,mst,mend,isp,minm)
          endif
        else 
          call io_minm('r',nmflag,minm,nst,nend,mst,mend,ik,isp)
          call trans_minm(iop,nmflag,minm,nst,nend,mst,mend,ik,iq,isp) 
        endif 
        call io_minm('w',nmflag,minm,nst,nend,mst,mend,ik,isp)

        deallocate(minm) 
 100    format("All Mi",a2," for n/c=[",i4,",",i4,'], and m/c=[',i4,',',&
     &     i4,'] are already contained in the Minm file')
 101    format("All Mi",a2," for n/c=[",i4,",",i4,'], and m/c=[',i4,',',&
     &     i4,'] are new and have to be calculated from the scratch')
 102    format("Mi",a2," for n/c= ",i4,",m/c=",i4,                      &
     &     ' is new and has to be calculated from the scratch')
        endsubroutine

      end subroutine set_minm
!EOC

