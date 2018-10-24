!BOP
! 
! !ROUTINE: task_kpinterp

!
! !INTERFACE:
      subroutine task_kpinterp 

! !DESCRIPTION:

! This subroutine calculate gw qp energies at a dense mesh using kp-interpolation 
 

! !USES:

      use bands,     only: nspin,ibgw,nbgw,nbandsgw,eferks,numin, &
     &                     metallic,nomax
      use constants, only: hev, pi
      use struk,     only: nat,br2,rbas,alat,ortho
      use kmeshintp, only: nbandsintp,nvbm,ncbm,xcoorband,kk0ind,       &
     &                     nbands2,eks1,eqp1,klist1,nkp1,idvk1,         &
     &                     eks2,eqp2,klist2,nkp2,idvk2,dek1,dek2
      use kpoints,   only: klist, kpirind, nkp, idvk
      use mommat,    only: mmatvv, init_mommat, end_mommat
      use struk,     only: br2
      use task,      only: casename 
      use liboct_parser

!
! !LOCAL VARIABLES:

      implicit none

      integer(4) :: ib,ib0,ib1,jb
      integer(4) :: ik,ikvbm(1),ikcbm(1),ikstart,ikend 
      integer(4) :: ik2, irk
      integer(4) :: isp
      integer(4) :: nvm,ncm
      integer(4) :: fid
      integer(4) :: info 
      integer(4) :: ik0vec(3), ik2vec(3),ig(3)
      
      real(8) :: ksvbm(nspin),kscbm   !! VBM and CBM from KS band structure 
      real(8) :: qpvbm(nspin),qpcbm   !! VBM and CBM from QP band structure 
      real(8) :: global_ksvbm, global_qpvbm
      real(8) :: ksgap,qpgap
      real(8) :: efryd
      real(8) :: k0vec(3), k2vec(3),gvec(3)
      real(8), allocatable :: eval0(:),eval1(:)
      
      complex(8), allocatable :: cmat(:,:)
      complex(8), allocatable :: pmat(:,:,:)
      
      
      character(13) :: sname = 'task_kpinterp',blk_band='band'
      character(80) :: fname 
      character(2)  :: spflag
   
! !REVISION HISTORY:

! created on July 1st, 2004 by XZL
! Last Modified 06.03.2006 by RGA

!EOP

!BOC

      info = loct_parse_isdef(blk_band)
      if(info.eq.1) then
         call loct_parse_block_float(blk_band,0,0,efryd)
         call loct_parse_block_int(blk_band,0,1,ikstart)
         call loct_parse_block_int(blk_band,0,2,ikend)
      else
         efryd=0.0d0
         ikstart=1
         ikend=0
      endif 
      eferks=efryd*0.5d0
!
!     Read the quasiparticle energies
!
      call linmsg(6,'-',"BAND STRUCTURE")
      write(6,*) "band: readeqp"
      call readeqp 
!
!    Read DFT energy on the larger k-mesh 
!
      write(6,*) "band:  readeband"
      call readeband(ikstart,ikend) 
      eks2(:,:,:)=eks2(:,:,:)-eferks

! The index for the highest band to be interpolated 
      ib0=ibgw 
      ib1=min(nbands2,nbgw) 
      nbandsintp=ib1-ib0+1
      write(6,'(a,3i5)') ' ib0,ib1,nbandsintp=',ib0,ib1,nbandsintp 
      nomax=ib1
      numin=ib0
      
!!
!! Transform klist band to internal coordinates 
!!
!      if(ortho)then
!        call cart2int0(nkp1,rbas,alat,klist1,idvk1)
!        call cart2int0(nkp2,rbas,alat,klist2,idvk2)
!      endif
!
!     Calculate the momentum matrix elements
!
      call calcmommat(0,ib0,ib1,ib0,ib1,1)
      call setkk0index
             
      write(6,*)
      write(6,*)' band: calc bandstruct using k.p perturbation theory'
      write(6,*)
      
      

      if(nspin.eq.1) then 
        spflag=''
      else 
        spflag='up'
      endif 
      fid=999

      call boxmsg(6,'-',"BAND STRUCTURE SUMMARY")

      do isp=1,nspin 
        call init_mommat(isp,ibkp,nbkp,ibkp,nbkp)

        if(nspin.eq.2) write(6,*) "BANDSTRUCT for Spin ",isp
        if(isp.eq.2)  spflag='dn' 
        nvm=nvbm(isp)
        ncm=ncbm(isp)
        allocate(cmat(nbandsintp,nbandsintp))
        allocate(eval0(nbandsintp))
        allocate(eval1(nbandsintp))
        allocate(pmat(1:3,nbandsintp,nbandsintp))

        do ik2=1,nkp2
          ik=kk0ind(1,ik2)
          ig(1:3)=kk0ind(2:4,ik2)
          gvec(1:3)=dble(ig(1))*br2(1:3,1)+dble(ig(2))*br2(1:3,2)+      &
     &                    dble(ig(3))*br2(1:3,3)
          irk=kpirind(ik)
          eval0(1:nbandsintp)=eqp1(ib0:ib1,irk,isp)
          do ib=ib0,ib1
            do jb=ib0,ib1
              pmat(1:3,jb-ib0+1,ib-ib0+1)=mmatvv(1:3,jb,ib,ik,isp)
            enddo
          enddo    
          write(6,*)'kpseceqn for ik=',ik2
          call flushbuf(6)
          ik2vec=klist2(1:3,ik2)
          call k2cart(ik2vec,idvk2,k2vec)
          ik0vec=klist(1:3,ik)
!          k0vec(1:3)=k0vec(1:3)
          call k2cart(ik0vec,idvk,k0vec)
          k0vec(1:3)=k0vec(1:3)+gvec(1:3)
          call kpseceqn(1,nbandsintp,nbandsintp,k0vec,eval0,pmat,  &
     &                  k2vec,eval1,cmat)     
          eqp2(ib0:ib1,ik2,isp)=eval1(1:nbandsintp)
        enddo 
        deallocate(cmat,eval0,eval1,pmat)

!
! Bands that are not calculated in GW are shifted by a constant for valence and couductance band, respectively
!
!        eqp2(1:ib0-1,:,isp)= eks2(1:ib0-1,:,isp) & 
!     &               +sum(real(dek2(:,ib0:nvm)))/(nkp2*(nvm-ib0+1)) 

!        eqp2(ib1+1:nbands2,:,isp)= eks2(ib1+1:nbands2,:,isp) & 
!     &               +sum(real(dek2(:,ncm:ib1)))/(nkp2*(ib1-ncm+1)) 


!
! take the VBM as the zero point and convert to the unit of eV
!
!
! DFT band structure 
!
 
        write(6,*) " Band of VBM: ",nvm
        write(6,*) " Band of CBM: ",ncm
        write(6,'(a,f10.4)') " KS VBM (Ry.)=",maxval(eks2(nvm,1:nkp2,isp))*2.d0
        write(6,'(a,f10.4)') " KS CBM (Ry.)=",minval(eks2(ncm,1:nkp2,isp))*2.d0
        if(nvm.ge.ncm)then
          if(nspin.eq.1)then
            write(6,*)'The system is metallic'
          else
            write(6,*)'Spin ',spflag,' is metallic'
          endif
        endif    
        eks2(:,:,isp)=eks2(:,:,isp)*hev
        write(6,*) "  KS band structure:"
        ksvbm(isp)=maxval(eks2(nvm,1:nkp2,isp))
        kscbm=minval(eks2(ncm,1:nkp2,isp))
        ikvbm=maxloc(eks2(nvm,1:nkp2,isp))
        ikcbm=minloc(eks2(ncm,1:nkp2,isp))

        ksgap=kscbm - ksvbm(isp)
! This is only valid IF there is a gap and the system is not spin
! polarized     
!        eks2(:,:,isp)=eks2(:,:,isp)-ksvbm
          
        write(6,112)  ksgap
        write(6,113)  ikvbm,ikcbm

        if(ikvbm(1).ne.ikcbm(1)) then 
          write(6,*) "  --- indirect gap:"
          write(6,114) eks2(ncm,ikvbm,isp)-eks2(nvm,ikvbm,isp)
          write(6,115) eks2(ncm,ikcbm,isp)-eks2(nvm,ikcbm,isp)
        endif 


        write(6,*)
        write(6,*) "  GW band structure:"
        eqp2(:,:,isp)=eqp2(:,:,isp)*hev
        qpvbm(isp)=maxval(eqp2(nvm,1:nkp2,isp))
        qpcbm=minval(eqp2(ncm,1:nkp2,isp))
        ikvbm=maxloc(eqp2(nvm,1:nkp2,isp))
        ikcbm=minloc(eqp2(ncm,1:nkp2,isp))

        qpgap=qpcbm-qpvbm(isp)
! This is onlz valid IF there is a gap, otherwise the alignment is wrong!!!!        
!        eqp2(:,:,isp)=eqp2(:,:,isp)-qpvbm

        write(6,112) qpgap
        write(6,113)  ikvbm,ikcbm

        if(ikvbm(1).ne.ikcbm(1)) then 
          write(6,*) "  --- indirect gap:"
          write(6,114) eqp2(ncm,ikvbm,isp)-eqp2(nvm,ikvbm,isp)
          write(6,115) eqp2(ncm,ikcbm,isp)-eqp2(nvm,ikcbm,isp)
        endif

!
! Print out the range of each band 
!
        write(6,*) " Range of each band:"
        write(6,'(a5,6a12)') '  n  ','KS Min','KS Maxval','KS width', &
     &      'GW Min',' GW Max ',' GW width ' 
        do ib=ib0,ib1
          write(6,'(i5,6F12.6)') ib,                               &
     &     minval(eks2(ib,1:nkp2,isp)),maxval(eks2(ib,1:nkp2,isp)),  & 
     &     maxval(eks2(ib,1:nkp2,isp))-minval(eks2(ib,1:nkp2,isp)),  &
     &     minval(eqp2(ib,1:nkp2,isp)),maxval(eqp2(ib,1:nkp2,isp)),  &
     &     maxval(eqp2(ib,1:nkp2,isp))-minval(eqp2(ib,1:nkp2,isp))
        enddo 
      enddo  ! isp

      if(nspin.eq.1) then 
        spflag=''
      else 
        spflag='up'
      endif 

      do isp=1,nspin
        if(isp.eq.2)  spflag='dn' 
        if(.not.metallic)then
          global_ksvbm=maxval(ksvbm)         
          global_qpvbm=maxval(qpvbm)         
          eks2(:,:,:)=eks2(:,:,:)-global_ksvbm
          eqp2(:,:,:)=eqp2(:,:,:)-global_qpvbm
        endif  
!
!     Write the bandstructure to disk
!
        write(6,*) "task_band: write bandstructure"

        fname=trim(casename)//'.bandeks'//trim(spflag) 
        open(unit=fid,file=fname,action='write',iostat=info) 
        call errmsg(info.ne.0,sname,"Fail to open file "//trim(fname))

        fname=trim(casename)//'.bandeqp'//trim(spflag)
        open(unit=fid+1,file=fname,action='write',iostat=info)
        call errmsg(info.ne.0,sname,"Fail to open file "//trim(fname))

        do ib=ib0,ib1    
          do ik=1,nkp2
            write(fid,111)  xcoorband(ik),eks2(ib,ik,isp)  
            write(fid+1,111)xcoorband(ik),eqp2(ib,ik,isp)  
          enddo
          write(fid,*)
          write(fid+1,*)
        enddo
        call flushbuf(fid)
        call flushbuf(fid+1)
        close(fid) 
        close(fid+1) 
      enddo  ! isp
      call end_mommat
!      deallocate(dek1,dek2)

 100  format(3e19.12,a10,2i6,f5.1)
 111  format(2f16.8)
 112  format(10x,' Fundamental Band Gap = ',f12.4)
 113  format(10X," k-point index for VBM and CBM: ",2i5)
 114  format(10x," Gap at VBM",f12.4)
 115  format(10x," Gap at CBM",f12.4)

      return


      end subroutine task_kpinterp

!EOC
