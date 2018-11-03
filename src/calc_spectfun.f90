!BOP
!
! !ROUTINE: task_ppgw
!
! !INTERFACE:
       subroutine calc_spectfun(iop,n_min,n_max,om_max,tag) 
       use bands,      only: bande,nspin,eferqp,ibgw,nbgw,efermi,eqp,   &
     &                       nbandsgw,nbmax,fspin 
       use bzinteg,    only: kwt_ibz
       use bzint,      only: bzint_smear 
       use constants,  only: hev
       use freq,       only: omega,nomeg 
       use kpoints,    only: nirkp
       use selfenergy, only: selfx,sigc,vxcnn
       use task,       only: casename
       use liboct_parser
       implicit none
       
! !DESCRIPTION:
!
! This subroutine perform the analytic continuation of the selfenergy to
! real frequencies and calculates the quasi-particle energies. The
! selfenergy and exchange correlation potential are read from file, thus, a
! previous run of the GW cycle is needed.        
!
      integer,intent(in):: iop
      integer,intent(in):: n_min,n_max
      character(len=*)::tag 
      
      integer :: ierr,fid
      integer :: isp,ie,ip,irk,iom
      integer :: nom 
      real(8) :: om_min, om_max 
      real(8) :: enk,delta,omg,rtmp 
      real(8) :: exnk,om_step,efer
      real(8) :: sumdos_qp,sumdos_ks,sumspf_ef,sumspf
      real(8),allocatable :: oms(:),spf(:), spf_n(:),dos_qp(:),dos_ks(:)
      character(20)::sname="calc_spectfun"
      character(80)::fname 


! !EXTERNAL ROUTINES: 
      character(10),external::int2str


!EOP
!BOC
      !! set the real frequency grid 
      nom = 2000 
      allocate(oms(nom),spf(nom),spf_n(nom),dos_qp(nom),dos_ks(nom),&
     &   stat=ierr) 
      call errmsg(ierr.ne.0,sname,"fail to allocate oms,...") 
          
      om_step = (om_max - om_min)/(nom-1) 
      do iom=1,nom
        oms(iom) = om_min + (iom-1.0)*om_step 
      enddo 

      if(tag.eq.'') then 
        fname = trim(casename)
      else
        fname = trim(casename)//'-'//trim(tag) 
      endif 

      if(iop.eq.0) then 
        fid = 0 
      else
        fid = 999
        open(fid,file=trim(fname)//"-Ank.dat",action='write')
      endif

      spf = 0.d0 
      do isp=1,nspin 
        do irk=1,nirkp 
          do ie = n_min, n_max 
            exnk = bande(ie,irk,isp) - vxcnn(ie,irk,isp) + selfx(ie,irk,isp)

            if(fid.gt.0) then 
              write(fid,100) ie,irk, kwt_ibz(irk),bande(ie,irk,isp)
            endif 

            call calc_spectfun_nk(spf_n,oms+eferqp,nom,exnk,sigc(:,ie,irk,isp),&
     &         omega,nomeg,fid) 
            spf = spf + spf_n*kwt_ibz(irk)*fspin
          enddo
        enddo
      enddo 
      if(fid.gt.0) then 
        close(fid) 
      endif 

      ! calculate DOS directly from QP and KS energies 
      dos_qp = 0.d0 
      dos_ks = 0.d0 
      do isp=1,nspin
        do irk=1,nirkp
          do ie = n_min, n_max
            do iom=1,nom 
              dos_qp(iom) = dos_qp(iom)+fspin*kwt_ibz(irk) &
     &         *bzint_smear(2,oms(iom)-eqp(ie,irk,isp)+eferqp)
              dos_ks(iom) = dos_ks(iom)+fspin*kwt_ibz(irk) &
     &         *bzint_smear(2,oms(iom)-bande(ie,irk,isp))
            enddo 
          enddo 
        enddo
      enddo 

!     Determine the Fermi energy in terms of the spectral function 

      open(fid,file=trim(fname)//"-specfun.dat",action='write') 

      sumspf_ef = 0.d0
      sumspf = 0.d0
      sumdos_qp = 0.d0 
      sumdos_ks = 0.d0 
      do iom=1,nom 
        sumspf = sumspf + om_step * spf(iom)
        if(oms(iom).le.0.d0) then 
          sumspf_ef = sumspf_ef + om_step * spf(iom)
          sumdos_ks = sumdos_ks + om_step * dos_ks(iom) 
          sumdos_qp = sumdos_qp + om_step * dos_qp(iom)
        endif

        write(fid,120) oms(iom)*HeV,spf(iom)/HeV,dos_qp(iom)/HeV, &
     &   dos_ks(iom)/HeV 
      enddo 
      write(6,110) "Integrating DOS_QP up to E_Fermi",sumdos_qp 
      write(6,110) "Integrating DOS_KS up to E_Fermi",sumdos_ks
      write(6,110) "Integrating spectr. func. up to E_Fermi",sumspf_ef
      write(6,110) "Integrating spectr. func. to infinity ",sumspf
      close(fid)  
      deallocate(oms,spf,spf_n) 

      return
 100  format(/,"# Spectral function for n=",i4," irk=",i4," k-weight=",&
     & f10.6," Enk=",f12.6) 
 110  format(a40,F10.4) 
 120  format(4F12.6) 
      
      end subroutine calc_spectfun 
!EOC      

