!BOP
!
! !ROUTINE: w2k_genvectorf        
!
! !INTERFACE:
      subroutine w2k_genvectorf

! !DESCRIPTION:
!
! This subroutine reads the eigenvector coefficient of the LAPW
!wavefunctions      
!
! {\bf WIEN2k interface}
!
! !USES:

        use bands,       only: nbmax,nv,bande0,nspin,eferks
        use constants,   only: czero
        use eigenvec,    only: kzz,zzk,lsymvector,vfunit,vfname,vfrecl, &
     &                         lcmplx,zzkall
        use kpoints,     only: nirkp,nkp,kpirind,idikp,iksym,g0
        use recipvec,    only: indgkir,indgk,ig0,maxngk,ngkir,npw,gindex
        use struk,       only: nat,ortho,imat,izmat
        use task,        only: casename,spflag,iop_scratch 
        use modmpi,      only: myrank

        implicit none

! !LOCAL VARIABLES:
      
        integer:: i,iat,ik,irk,ib,ig,isp,isym,irec
        integer:: nvk,nbk,nwf
        integer:: ierr,itmp
        integer:: nrecl
        integer:: fid
        integer:: symat(3,3),igvir(3),igvec(3) 
        character(10),external::int2str
 
        real(8)    :: kvec(3),weight, eigval
        real(8)    :: emist

        integer,    allocatable:: kzz0(:,:)
        real(8),    allocatable:: zzk_r(:)
        complex(8), allocatable:: zzk_c(:)

        character(len=3)  ::    ipgr
        character(len=10) ::    kname   
        
        logical:: ldbg=.false.
        character( 20):: sname= 'w2k_genvectorf'        
        character(200):: msg 
                                     
!
! !REVISION HISTORY:
!
! Last mofified: 25th. March 2004 by RGA
!
!EOP
!BOC    

      call linmsg(6,'-',sname)

      if(iop_scratch.gt.0) then
        write(6,*) "  Prepare direct access vector file ",trim(vfname)  
        vfrecl=16*nbmax*maxngk
 
        if(lsymvector) then 
          nrecl=nirkp
        else 
          nrecl=nkp
        endif 
        write(6,*) "  Record size (KB): ",nint(1.D-3*vfrecl*4)
        write(6,*) "  Total  size (KB): ",nint(1.0D-3*nrecl*vfrecl*4)*nspin
        open(vfunit,file=vfname,action='readwrite',form='unformatted',  &
     &      access='direct',recl=vfrecl,iostat=ierr)
        call errmsg0(ierr,sname,"Fail to open "// trim(vfname))
      endif

      if(.not.lsymvector) then 
        allocate(kzz0(1:3,maxngk))
        kzz0=0
      endif 

      if(.not.lcmplx) then 
        allocate(zzk_r(maxngk))
        zzk_r = 0.d0
      else
        allocate(zzk_c(maxngk))
        zzk_c = 0.d0
      endif 

      fid=999
      do isp=1,nspin 
        open(unit=fid,file=trim(casename)//'.vector'//spflag(isp),&
     &       action='read',form='unformatted',iostat=ierr)
        call errmsg(ierr.ne.0,sname,"Fail to open case.vector")

        !! skip the head that contains the linearization eneriges
        do iat=1,nat
          read(fid) emist
          read(fid) emist
        enddo
!
! For eigenvector using symmetry
!
        if(lsymvector) then 
          do irk=1,nirkp  
            read(fid,iostat=ierr) kvec(1:3),kname,nwf,nbk,weight
            if (ierr.ne.0) then 
              write(6,*) "!!! ERROR when reading case.vector: kvec ..."
              goto 900
            endif 
           
            read(fid,iostat=ierr) (kzz(1:3,i,irk),i=1,nwf)
            if (ierr.ne.0) then
              write(6,*) "!!! ERROR when reading case.vector: kzz ..."
              goto 900
            endif

            do i=1,nwf
              indgkir(i,irk)=ig0(kzz(1,i,irk),kzz(2,i,irk),kzz(3,i,irk))
            enddo

            zzk=czero
            do ib=1,nbk
              read(fid)  itmp,eigval  

              if(lcmplx) then 
                read(fid,iostat=ierr) (zzk_c(i),i=1,nwf) 
              else 
                read(fid,iostat=ierr) (zzk_r(i),i=1,nwf)
              endif 

              if (ierr.ne.0) then
                write(6,*) "!!! ERROR when reading case.vector: zzk ..."
                goto 900
              endif

              if(ib.gt.nbmax) cycle 
              
              eigval=eigval*0.5d0-eferks
              if( abs(bande0(ib,irk,isp)-eigval).gt.1.d-5) then
                 write(6,*) " !! WARNING !! "
                 write(6,*) " - Inconsistent vector and energy files"
                 write(6,'(2f20.10)') eigval, bande0(ib,irk,isp)
              endif
              if(lcmplx) then 
                zzk(1:nwf,ib) = zzk_c(1:nwf)
              else 
                zzk(1:nwf,ib)=cmplx(zzk_r(1:nwf),0.0d0,8)
              endif 

            enddo ! ib

            if(iop_scratch.gt.0) then
              irec = irk + (isp-1)*nirkp
              if(iop_scratch.eq.1.or.(iop_scratch.eq.2.and.myrank.eq.0))&
     &         write(vfunit,rec=irec) zzk

            else
              zzkall(:,:,irk,isp)=zzk
            endif

          enddo
 
          !! set the indgk array from indgkir by applying the symetries 
          do ik=1,nkp
            irk=kpirind(ik)
            if(idikp(irk).eq.ik)then
              indgk(:,ik)=indgkir(:,irk)
            else
              isym=iksym(ik)
              if(ortho)then
                symat(1:3,1:3)=imat(1:3,1:3,isym)
              else
                symat(1:3,1:3)=izmat(1:3,1:3,isym)
              endif

              do ig=1,ngkir(irk)
                igvir(1:3)=gindex(:,indgkir(ig,irk))
                igvec=matmul(symat,igvir)-g0(:,ik)
                indgk(ig,ik)=ig0(igvec(1),igvec(2),igvec(3))
              enddo ! ig
            endif
          enddo ! ik

        else 
!
! For eigenvectors without using any symmetry
!
          do ik=1,nkp
            irk=kpirind(ik)
            read(fid,iostat=ierr) kvec(1:3),kname,nwf,nbk,weight

            if (ierr.ne.0) then
              write(6,*) "!!! ERROR when reading case.vector: kvec ..."
              goto 900
            endif


            if(ngkir(irk).ne.nwf ) then
              write(msg,100)  "inconsistency occurs!!", &
     &              "  ik,irk=",ik,irk,                 &
     &              "  ngkir(irk),nwf=",ngkir(irk),nwf
              call outerr(sname,msg)
            endif
 
            read(fid,iostat=ierr) (kzz0(1:3,i),i=1,nwf)
            if (ierr.ne.0) then
              write(6,*) "!!! ERROR when reading case.vector: kzz"
              goto 900
            endif

            do i=1,nwf
              indgk(i,ik)=ig0(kzz0(1,i),kzz0(2,i),kzz0(3,i))
              if(indgk(i,ik).gt.npw) then 
                write(msg,*) "indgk(i,ik) > npw! --- i,ik,indgk=", &
     &                      i,ik,indgk
                call outerr(sname,msg)
              endif 
            enddo

            if(idikp(irk).eq.ik) then 
              do i=1,nwf
                kzz(1:3,i,irk)=kzz0(1:3,i) 
                indgkir(i,irk)=indgk(i,ik)
              enddo
              if(ldbg) write(6,*) "w2k_genvectorf: irk=",irk
            endif 

            zzk=czero
            do ib=1,nbk
              read(fid) itmp,eigval
              if(lcmplx) then
                read(fid,iostat=ierr) (zzk_c(i),i=1,nwf)
              else
                read(fid,iostat=ierr) (zzk_r(i),i=1,nwf)
              endif

              if (ierr.ne.0) then
                write(6,*) "!!! ERROR when reading case.vector: zzk"
                goto 900
              endif

              if(ib.gt.nbmax) cycle 

              eigval=eigval*0.5d0-eferks
              if( abs(bande0(ib,irk,isp)-eigval).gt.1.d-5) then 
                 write(6,*) " !! WARNING !! "
                 write(6,*) " -- Inconsistent vector and energy files"
                 write(6,'(2f20.10)') eigval, bande0(ib,irk,isp) 
              endif 

              if(lcmplx) then 
                zzk(1:nwf,ib) = zzk_c(1:nwf)
              else 
                zzk(1:nwf,ib) = cmplx(zzk_r(1:nwf),0.0d0,8)
              endif
            enddo ! ib 

            if(iop_scratch.gt.0) then
              irec = ik + (isp-1)*nkp
              if(iop_scratch.eq.1.or.(iop_scratch.eq.2.and.myrank.eq.0))&
     &         write(vfunit,rec=irec) zzk
            else
              zzkall(:,:,ik,isp)=zzk
            endif
          enddo
        endif !! lsymvector 

        close(fid)
      enddo  ! loop over isp

      if(.not. lsymvector) deallocate(kzz0)
      if(.not. lcmplx ) deallocate(zzk_r) 

      return
  900 call outerr(sname,"error reading file case.vector")
  100 format(a,/,a10,2i5,/,a20,2i5)
  
      end subroutine w2k_genvectorf
!EOC      
      
      
