!BOP
!
! !ROUTINE: io_eqp
!
! !INTERFACE:
      subroutine io_eqp(rw,isxc,iform,flag)

! !DESCRIPTION:
! 
! This subroutine writes the qp energies to file
!
! !USES:
!
      use bands,     only: bande0,ibgw,nbgw,nspin,nomaxs,numins,eqp, &
     &                     metallic,eferks,eferqp,eferhf,eqp_im
      use constants, only: hev
      use freq,      only: nomeg
      use kpoints,   only: idvkir, kirlist, nirkp, idikp
      use recipvec,  only: ngk
      use selfenergy, only: sigx,sigc,sacpar,npar_ac,iop_es,iop_ac,      &
     &                      selfc,selfx,znorm,vxcnn
      use task,      only: casename
      
       
! !LOCAL VARIABLES:

      implicit none
      character :: rw             !! 'r' for reading and 'w' for writing 
      integer,intent(in):: isxc   !! 0/1/2/3 (GW/HF/COHSEX/SEX) 
      integer,intent(in):: iform  !! format to write the eqp 
                                   !! 0 -- *.eqpH
                                   !! 1 -- *.eqpeV
                                   !! 2 -- *.eqpeV, but with only KS and QP energies  
      character(len=*),intent(in)::flag !! flag to the qp energy file to differentiate different
                                        !! QP energies. e.g."_g0w0","gw0","ex" 
      integer:: fid    
      integer(4) :: isp
      integer(4) :: ie   !(Counter) Runs over bands
      integer:: ierr
      integer(4) :: ik  !(Counter) Runs over k-points
      integer(4) :: irk  !(Counter) Runs over irreducible k-points
      integer(4), dimension(3) :: ikvec
      integer(4) :: nv,ncm,nv1,nv2,ikvm(1)
      
      real(8) :: enk,ehf,eks,egw,degw,dehf,sxnk,scnk,vxcnk,znk,degwvm
      real(8) :: evmks,evmgw,evmhf
      character(80)::fname
      character(10)::sname="io_eqp"

! !REVISION HISTORY:
!
! Created 16.08.05 by RGA
!
!EOP
!BOC

      fid=999

!
!     set file name 
!
      if(iform.eq.0) then
        if(flag.eq.'') then
          fname=trim(casename)//".eqpH"
        else
          fname=trim(casename)//".eqpH_"//trim(flag)
        endif
      else 
        if(flag.eq.'') then
          fname=trim(casename)//".eqpeV"
        else
          fname=trim(casename)//".eqpeV_"//trim(flag)
        endif
      endif 
!
!     Write qp energy
!
      if(rw.eq.'w'.or.rw.eq.'W') then 
        open(unit=fid,file=fname,action='write',iostat=ierr)
        call errmsg(ierr.ne.0,sname,"Fail to open "//trim(fname))

        if(iform.eq.0) then 
          write(fid,11) ibgw,nbgw, nirkp, nspin,eferks 
          do isp=1,nspin 
            do irk=1,nirkp
              ik=idikp(irk)
              ikvec(1:3)=kirlist(1:3,irk)
              write(fid,12) irk,ik,ikvec,idvkir
              do ie=ibgw,nbgw
                eks=bande0(ie,irk,isp)
                egw=eqp(ie,irk,isp)
                write(fid,13)  ie, eks, egw
              enddo !! ie
              write(fid,*)
            enddo !! irk  
          enddo !! isp   

        elseif(iform.eq.1) then 
          write(fid,*) "# Quasiparticle Energy with iop_sxc=",isxc,&
     &               "(0:gw, 1:ex, 2:cohsex)"  
          write(fid,*) "# Fermi energy is set as zero" 
          write(fid,*) "# Number of spin channels (nspin) =", nspin
          write(fid,*) "# Number of k points (nirkp)     =",  nirkp
          write(fid,*) "# Band range (ibgw,nbgw)         =",  ibgw,nbgw

          if(isxc.eq.0) then 
            write(fid,*) "# Parameters used:"
            write(fid,*) "#  Analytic contiuation (iop_ac)  =", iop_ac
            write(fid,*) "#  Fermi level shift    (iop_es)  =", iop_es
            write(fid,*) "#  Nr.freq points  (nomeg)        =", nomeg
            write(fid,*) "#  Number of AC poles (npar_ac/2) =", npar_ac/2
          endif 

          if(metallic) then 
            evmks = 0.0
            evmgw = eferqp
            evmhf = eferhf
            degwvm = 0.d0  
          else 
            if(nspin.eq.1) then 
              nv=nomaxs(1)
              evmks = maxval(bande0(nv,:,1))
              evmgw=maxval(eqp(nv,:,1))
              evmhf=maxval(bande0(nv,:,1)+selfx(nv,:,1)-vxcnn(nv,:,1))

              !! get the GW correction to VBM 
              ikvm  = maxloc(bande0(nv,:,1))
              degwvm = eqp(nv,ikvm(1),1) - bande0(nv,ikvm(1),1)
            else 
              nv1=nomaxs(1)
              nv2=nomaxs(2)
              evmks=max(maxval(bande0(nv1,:,1)),maxval(bande0(nv2,:,2)))
              evmgw=max(maxval(eqp(nv1,:,1)),   maxval(eqp(nv2,:,2)))
              evmhf=max(maxval(bande0(nv1,:,1)+selfx(nv1,:,1)-vxcnn(nv1,:,1)),&
     &                  maxval(bande0(nv2,:,1)+selfx(nv2,:,1)-vxcnn(nv2,:,1)))

              ikvm  = maxloc(bande0(nv1,:,1))
              degwvm = eqp(nv1,ikvm(1),1) - bande0(nv1,ikvm(1),1)
              if(maxval(bande0(nv1,:,1)) < maxval(bande0(nv2,:,2))) then 
                ikvm  = maxloc(bande0(nv2,:,2))
                degwvm = eqp(nv2,ikvm(1),2) - bande0(nv2,ikvm(1),2)
              endif 
            endif 
          endif 
          write(6,'(a,f10.4)') ":DeltaE_QP(VBM) [eV] = ",degwvm*hev

          evmks = evmks*hev
          evmgw = evmgw*hev
          evmhf = evmhf*hev
          
          do isp=1,nspin
            if(nspin.eq.2) write(fid,*) '# ispin=',isp
            do irk=1,nirkp
              ik=idikp(irk)
              ikvec(1:3)=kirlist(1:3,irk)
              write(fid,21) irk,ik,ikvec,idvkir

              if(isxc.eq.0) then   !! GW output 
                write(fid,20)
              elseif(isxc.eq.1) then  !! EX-only output 
                write(fid,30) 
              elseif(isxc.eq.2.or.isxc.eq.3) then  !! COHSEX output 
                write(fid,40)
              endif 

              do ie=ibgw,nbgw
                enk = bande0(ie,irk,isp)*hev
                sxnk = selfx(ie,irk,isp)*hev
                vxcnk = vxcnn(ie,irk,isp)*hev
                dehf = sxnk - vxcnk 
                eks = enk - evmks
                ehf = enk+sxnk-vxcnk-evmhf 

                if(isxc.eq.0) then   
                  scnk = selfc(ie,irk,isp)*hev
                  znk  = znorm(ie,irk,isp)
                  degw = (eqp(ie,irk,isp)-bande0(ie,irk,isp))*hev
                  egw  = eqp(ie,irk,isp)*hev-evmgw
                  write(fid,22) irk,ie,eks,egw,ehf,sxnk,scnk,vxcnk, &
     &                        degw,dehf,znk,eqp_im(ie,irk,isp)

                elseif(isxc.eq.1) then 

                  write(fid,32) irk,ie,eks,ehf,sxnk,vxcnk,dehf

                elseif(isxc.eq.2.or.isxc.eq.3) then 
                  scnk = selfc(ie,irk,isp)*hev
                  degw = (eqp(ie,irk,isp)-bande0(ie,irk,isp))*hev
                  egw  = eqp(ie,irk,isp)*hev-evmgw
                  write(fid,42) irk,ie,eks,egw,ehf,sxnk,scnk,vxcnk,&
     &                        degw,dehf
                endif 
              enddo
            enddo
            write(fid,*)
          enddo

        elseif(iform.eq.2) then 
          write(fid,*) "# QSGW QP Energy with iop_sxc=",isxc,&
     &               "(0:gw, 1:ex, 2:cohsex)"  
          write(fid,*) "# Fermi energy is set as zero" 
          write(fid,*) "# Number of spin channels (nspin) =", nspin
          write(fid,*) "# Number of k points (nirkp)     =",  nirkp
          write(fid,*) "# Band range (ibgw,nbgw)         =",  ibgw,nbgw

          if(isxc.eq.0) then 
            write(fid,*) "# Parameters used:"
            write(fid,*) "#  Analytic contiuation (iop_ac)  =", iop_ac
            write(fid,*) "#  Fermi level shift    (iop_es)  =", iop_es
            write(fid,*) "#  Nr.freq points  (nomeg)        =", nomeg
            write(fid,*) "#  Number of AC poles (npar_ac/2) =", npar_ac/2
          endif 

          if(metallic) then 
            evmks = 0.0
            evmgw = eferqp
          else 
            if(nspin.eq.1) then 
              nv=nomaxs(1)
              evmks = maxval(bande0(nv,:,1))
              evmgw=maxval(eqp(nv,:,1))

              !! get the GW correction to VBM 
              ikvm  = maxloc(bande0(nv,:,1))
              degwvm = eqp(nv,ikvm(1),1) - bande0(nv,ikvm(1),1)
            else 
              nv1=nomaxs(1)
              nv2=nomaxs(2)
              evmks=max(maxval(bande0(nv1,:,1)),maxval(bande0(nv2,:,2)))
              evmgw=max(maxval(eqp(nv1,:,1)),   maxval(eqp(nv2,:,2)))

              ikvm  = maxloc(bande0(nv1,:,1))
              degwvm = eqp(nv1,ikvm(1),1) - bande0(nv1,ikvm(1),1)
              if(maxval(bande0(nv1,:,1)) < maxval(bande0(nv2,:,2))) then 
                ikvm  = maxloc(bande0(nv2,:,2))
                degwvm = eqp(nv2,ikvm(1),2) - bande0(nv2,ikvm(1),2)
              endif 
            endif 
            write(6,'(a,f10.4)') ":DeltaE_QP(VBM) [eV] = ",degwvm*hev
          endif 

          evmks = evmks*hev
          evmgw = evmgw*hev
          do isp=1,nspin
            if(nspin.eq.2) write(fid,*) '# ispin=',isp
            do irk=1,nirkp
              ik=idikp(irk)
              ikvec(1:3)=kirlist(1:3,irk)
              write(fid,21) irk,ik,ikvec,idvkir

              write(fid,60)
              do ie=ibgw,nbgw
                enk = bande0(ie,irk,isp)*hev
                eks = enk - evmks
                degw=(eqp(ie,irk,isp)-bande0(ie,irk,isp))*hev
                egw = eqp(ie,irk,isp)*hev-evmgw
                write(fid,61) irk,ie,eks,egw,degw
              enddo
            enddo
            write(fid,*)
          enddo
        endif 
        close(fid)
      endif 
!
!     Read QP energies
!
          
  11  format(4i6,e20.12)
  12  format(6i6)
  13  format(i4,2e20.12) 
  20  format('# irk',3x,'ie',6x,'E_KS',6x,'E_GW',6x,'E_HF', &
     &       8x,'Sx',8x,'Sc',7x,'Vxc',5x,'DE_GW',5x,'DE_HF',&
     &       7x,'Znk',6x,'ImSc')
  21  format('# irk=',i4,'ik=',i4,'  ikvec=',3i4,'  idvk=',i4)
  22  format(2i5,10f10.3)

  30  format('# irk',3x,'ie',6x,'E_KS',6x,'E_HF',8x,'Sx', &
     &        7x,'Vxc',5x,'DE_HF')
  32  format(2i5,5f10.3)

  40  format('# irk',3x,'ie',6x,'E_KS',4x,'E_CHSX',6x,'E_HF',8x,'Sx', &
     &        8x,'Sc',7x,'Vxc',3x,'DE_CHSX',5x,'DE_HF')
  42  format(2i5,10f10.3) 


  60  format('# irk',3x,'ie',6x,'E_KS',6x,'E_GW',6x,'DE_GW')
  61  format(2i5,3f10.3) 

      return 
      end subroutine io_eqp
!EOC      
