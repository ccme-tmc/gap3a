!BOP
!
! !ROUTINE: readvector        
!
! !INTERFACE:
      subroutine readvector(ik,iop,isp,iop_trans)

! !DESCRIPTION:
!
! This subroutine reads the eigenvector coefficient of the LAPW
!wavefunctions, adapted for the use in the self-consistent GW calculations      
!
! {\bf WIEN2k interface}
!
! !USES:

        use bands,       only: nbmax,nbandsgw,ibgw,nbgw
        use constants,   only: czero, pi
        use eigenvec,    only: zzk,zzq,lsymvector,vfunit,vfname,zzkall
        use kpoints,     only: kpirind,iksym,klist,idvk,g0,kirlist,nkp,nirkp
        use struk,       only: tau,imat,izmat,ortho
        use recipvec,    only: maxngk, gindex, indgk,ngkir,ig0,indgkir
        use selfenergy,  only: qpwf_coef
        use task,        only: iop_scratch,time_evec

        implicit none

! !INPUT PARAMETERS:
        
        integer, intent(in) :: ik    ! index of the corresponding k-point
        integer, intent(in) :: iop   ! indicates whether an eigenvector (iop==1) or its c.c. (iop==2)  should be read 
        integer, intent(in) :: isp   ! indicates from which spin file to read
        integer, intent(in):: iop_trans ! control whether the vectors are transformed (iop_trans > 0)
        
!EOP
!BOC      
! !LOCAL VARIABLES:
      
        integer:: ngs
        integer:: irk,irec
        integer:: indgs(maxngk)
        integer:: isym
        integer:: ib,ig,igt,indg,indgt
        integer:: i,ierr
        integer:: ngk
        integer:: igvec(3),rmat(3,3),igtvec(3)
        
        real(8):: eval(nbmax)
        real(8):: lentau
        real(8):: arg
        real(8):: gkvec(3),tvec(3)
        real(8):: tstart,tend
        
        complex(8) :: fact
        complex(8),allocatable:: zztmp(:,:) , phase(:)
        
! !DEFINED PARAMETERS:
 
        character(len=10) , parameter :: sname = 'readvector'        
                                     
!
! !REVISION HISTORY:
!
! Last mofified: 25th. March 2004 by RGA
!

      call cpu_time(tstart)
      allocate(zztmp(maxngk,nbmax),phase(maxngk),stat=ierr) 
      call errmsg(ierr.ne.0,sname,"Fail to allocate zztemp")

      irk=kpirind(ik)
      ngk=ngkir(irk)

      if(.not.lsymvector)then  !!  no translation related to the symmetry operation
        irec = ik+(isp-1)*nkp
        if(iop_scratch.gt.0) then
          read(vfunit,rec=irec,iostat=ierr) zztmp
          call errmsg(ierr.ne.0,sname,"Fail to read vectord file")
        else
          zztmp=zzkall(:,:,ik,isp)
        endif

        if(iop.eq.1)then
          do ib=1,nbmax 
            zzk(:,ib) = zztmp(:,ib) 
          enddo 
        else
          do ib=1,nbmax
            zzq(:,ib) = conjg(zztmp(:,ib))
          enddo
        endif

      else !! translation related to the symmetry operation is needed
        irec = irk + (isp-1)*nirkp
        if(iop_scratch.gt.0) then
          read(vfunit,rec=irec,iostat=ierr) zztmp
          call errmsg(ierr.ne.0,sname,"Fail to read vectord file")
        else
          zztmp=zzkall(:,:,irk,isp)
        endif
        isym=iksym(ik)
        tvec=tau(:,isym)
        if(ortho) then 
          rmat=imat(:,:,isym)
        else 
          rmat=izmat(:,:,isym)
        endif 

        do ig=1,ngk 
          indg=indgkir(ig,irk)     ! the G index corresponding to ig
          igvec(:)=gindex(:,indg)
!          igvec=matmul(rmat,igvec)  !-g0(:,ik)
          gkvec(:)=dble(igvec)+dble(klist(:,ik))/dble(idvk)
          arg=sum(gkvec*tvec)
          phase(ig)=cmplx(cos(2.0d0*pi*arg),-sin(2.0d0*pi*arg),8)
        enddo
        
        do ib=1,nbmax 
          if(iop.eq.1) then 
            zzk(1:ngk,ib)=phase(1:ngk)*zztmp(1:ngk,ib)
          else 
            zzq(1:ngk,ib)=conjg(phase(1:ngk)*zztmp(1:ngk,ib))
          endif 
        enddo
      endif !! lsymvector
  
      deallocate(zztmp,phase) 

      if(iop_trans.gt.0) then 
        call trans_vector(iop,ibgw,nbgw,qpwf_coef(:,:,irk,isp))
      endif 
      call cpu_time(tend)
      time_evec = time_evec + tend - tstart
      return
  
      end subroutine readvector
!EOC      
      
      
