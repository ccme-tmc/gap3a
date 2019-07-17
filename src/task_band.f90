!BOP
! 
! !ROUTINE: task_band

!
! !INTERFACE:
      subroutine task_band 

! !DESCRIPTION:

! This subroutine calculate gw band structure by interpolating the qp energy corrections 

! !USES:

      use constants, only: hev, pi
      use struk,     only: nat,br2,rbas,ortho
      use kmeshintp, only: nvbm,ncbm,                        &
     &                     nbands2,eks1,eqp1,nkp1,idvk1,         &
     &                     eks2,eqp2,kvecs1,kvecs2,nkp2,dek1,dek2
      use task,      only: casename 
      use liboct_parser

!
! !LOCAL VARIABLES:

      implicit none

      integer(4) :: i 
      integer(4) :: ii 
      integer(4) :: iat 
      integer(4) :: ib,ib0,ib1 
      integer(4) :: ik,ikvbm(1),ikcbm(1),ikstart,ikend 
      integer(4) :: isp
      integer(4) :: j 
      integer(4) :: kpid1 
      integer(4) :: m 
      integer(4) :: n 
      integer(4) :: l 
      integer(4) :: nbi 
      integer(4) :: nmbd 
      integer(4) :: nwf
      integer(4) :: nvm,ncm
      logical::    lprt=.false.
      integer(4) :: fid
      integer(4) :: info 

      integer(4),  dimension(3) :: onek
      
      integer(4), allocatable :: kkp(:,:,:)
      integer(4), allocatable :: ktmp(:,:)
      
      real(8) :: kx 
      real(8) :: ky 
      real(8) :: kz 
      real(8) :: w 
      real(8) :: ksvbm,kscbm   !! VBM and CBM from KS band structure 
      real(8) :: qpvbm,qpcbm   !! VBM and CBM from QP band structure 
      real(8) :: resu 
      real(8) :: x 
      real(8) :: x0
      real(8) :: ksgap,qpgap
      
      real(8), dimension(3) :: coordinas
      real(8), dimension(8) :: nod0

      real(8), allocatable :: delta(:,:,:)

      character(10) :: kname

      character(13), parameter :: sname = 'task_band',blk_band='band'
      character(80) :: fname 
      character(2)  :: spflag
   
! !REVISION HISTORY:

! created on July 1st, 2004 by XZL
! Last Modified 06.03.2006 by RGA

!EOP

!BOC

#ifdef DEBUG
      lprt = .true.
#endif 


      info = loct_parse_isdef(blk_band)
      if(info.eq.1) then
         call loct_parse_block_int(blk_band,0,0,ikstart)
         call loct_parse_block_int(blk_band,0,1,ikend)
      else
         ikstart=1
         ikend=0
      endif 
      
!
!     Read the quasiparticle energies
!
      call linmsg(6,'-',"BAND STRUCTURE")
      if(lprt) write(6,*) "band: readeqp"
      call readeqp 
!
!    Read DFT energy on the larger k-mesh 
!
      if(lprt) write(6,*) "band:  readeband"
      call readeband(ikstart,ikend) 

! The index for the highest band to be interpolated 
      ib0=ibgw 
      ib1=min(nbands2,nbgw) 
      nbandsintp=ib1-ib0+1
      write(6,'(a,3i5)') ' ib0,ib1,nbandsintp=',ib0,ib1,nbandsintp 
!
! Transform klist band to internal coordinates 
!
      if(ortho)then
        call cart2int(nkp1,rbas,alat,kvecs1)
        call cart2int(nkp2,rbas,alat,kvecs2)
      endif
!
!     Calculate the correction to the energies
!
      allocate(dek1(nkp1,ib0:ib1),dek2(nkp2,ib0:ib1))

      write(6,*)
      write(6,*)' band: calc bandstruct using Fourier interpolation'
      write(6,*)

      if(nspin.eq.1) then 
        spflag=''
      else 
        spflag='up'
      endif 
      fid=999

      call boxmsg(6,'-',"BAND STRUCTURE SUMMARY")

      do isp=1,nspin 
        if(nspin.eq.2) write(6,*) "BANDSTRUCT for Spin ",isp
        if(isp.eq.2)  spflag='dn' 
        nvm=nvbm(isp)
        ncm=ncbm(isp)

        do ik=1,nkp1
          dek1(ik,:) =cmplx(eqp1(ib0:ib1,ik,isp)-eks1(ib0:ib1,ik,isp),0.d0,8) 
        enddo 

        call fourintp(dek1,nkp1,kvecs1,dek2,nkp2,kvecs2,nbandsintp)

        do ib=ib0,ib1
          do ik=1,nkp2
            eqp2(ib,ik,isp)=eks2(ib,ik,isp)+real(dek2(ik,ib))
          enddo 
        enddo
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

        eks2(:,:,isp)=eks2(:,:,isp)*hev
        write(6,*) "  KS band structure:"
        ksvbm=maxval(eks2(nvm,1:nkp2,isp))
        kscbm=minval(eks2(ncm,1:nkp2,isp))
        ikvbm=maxloc(eks2(nvm,1:nkp2,isp))
        ikcbm=minloc(eks2(ncm,1:nkp2,isp))

        ksgap=kscbm - ksvbm
        eks2(:,:,isp)=eks2(:,:,isp)-ksvbm
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
        qpvbm=maxval(eqp2(nvm,1:nkp2,isp))
        qpcbm=minval(eqp2(ncm,1:nkp2,isp))
        ikvbm=maxloc(eqp2(nvm,1:nkp2,isp))
        ikcbm=minloc(eqp2(ncm,1:nkp2,isp))

        qpgap=qpcbm-qpvbm
        eqp2(:,:,isp)=eqp2(:,:,isp)-qpvbm

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

!
!     Write the bandstructure to disk
!
        if(lprt) write(6,*) "task_band: write bandstructure"

        fname=trim(casename)//'.bandeks'//trim(spflag) 
        open(unit=fid,file=fname,action='write',iostat=info) 
        call errmsg(info.ne.0,sname,"Fail to open "//trim(fname))

        fname=trim(casename)//'.bandeqp'//trim(spflag)
        open(unit=fid+1,file=fname,action='write',iostat=info)
        call errmsg(info.ne.0,sname,"Fail to open "//trim(fname))

        do ib=ib0,ib1    
          do ik=1,nkp2
            write(fid,111)  ik,eks2(ib,ik,isp)  
            write(fid+1,111)ik,eqp2(ib,ik,isp)  
          enddo
          write(fid,*)
          write(fid+1,*)
        enddo
        close(fid) 
        close(fid+1) 
      enddo  ! isp

      deallocate(dek1,dek2)

 100  format(3e19.12,a10,2i6,f5.1)
 111  format(i5,10f16.8)
 112  format(10x,' Fundamental Band Gap = ',f12.4)
 113  format(10X," k-point index for VBM and CBM: ",2i5)
 114  format(10x," Gap at VBM",f12.4)
 115  format(10x," Gap at CBM",f12.4)

      return


      end subroutine 

!EOC
