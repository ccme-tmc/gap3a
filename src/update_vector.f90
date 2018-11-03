      subroutine update_vector(coef)
!
! This subroutine updates the direct access vector file during the Quasi-particle 
! Self-consistent GW (QSGW) 
!
      use bands,       only: nbmax,nbandsgw,ibgw,nbgw,nspin
      use constants,   only: cone, czero
      use eigenvec,    only: lsymvector,vfunit,zzkall,zzk
      use kpoints,     only: kpirind,nirkp,nkp,idikp
      use recipvec,    only: maxngk, gindex, indgk,ngkir
      use task,        only: iop_scratch
      implicit none 
      complex(8), intent(in):: coef(ibgw:nbgw,ibgw:nbgw,nirkp,nspin)

      integer :: ib 
      integer :: ierr
      integer :: ik,irk,iik,irec
      integer :: isp
      integer :: ng,nb, nk
      complex(8),allocatable:: zztmp(:,:) 
      character(20):: sname="update_vector"

      allocate(zztmp(maxngk,nbmax),stat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to allocate zztmp")

      do isp=1,nspin 
        if(lsymvector) then 
          nk=nirkp
        else 
          nk=nkp
        endif 

        do iik=1,nk

          if(lsymvector) then 
            irk=iik
            ik=idikp(irk)
            irec = irk + (nspin-1)*nirkp
          else
            ik=iik
            irk=kpirind(ik)
            irec = ik + (nspin-1)*nkp
          endif

          if(iop_scratch.gt.0) then 
            read(vfunit,rec=irec,iostat=ierr) zztmp
            call errmsg(ierr.ne.0,sname,"Fail to read vectord file")
          else
            zztmp=zzkall(:,:,iik,isp)
          endif

          zzk=zztmp
          ng=maxngk
          nb=nbandsgw
          call zgemm('n','n',ng,nb,nb,cone,zztmp(:,ibgw:nbgw),maxngk,& 
     &              coef(:,:,irk,isp),nb,czero,zzk(:,ibgw:nbgw), &
     &              maxngk)

          if(iop_scratch.gt.0) then 
            write(vfunit,rec=irec,iostat=ierr) zzk
            call errmsg(ierr.ne.0,sname,"Fail to write vectord file")
          else 
            zzkall(:,:,iik,isp)=zzk
          endif 
        enddo !! iik
      enddo !! isp
      deallocate(zztmp)

      end subroutine

