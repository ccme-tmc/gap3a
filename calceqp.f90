!BOP
!
! !ROUTINE: calceqp
!
! !INTERFACE:
      subroutine calceqp(isxc,iod,isc_gw)

! !DESCRIPTION:
! 
! Given the matrix elements $\langle
! \Psi_{n\vec{k}}|\Sigma(\vec{k},\omega)|\Psi_{n\vec{k}}\rangle$, $\langle
! \Psi_{n\vec{k}}|V^{xc}|\Psi_{n\vec{k}}\rangle$ and
! $\varepsilon^{DFT}_{n\vec{k}}$. this subroutine calculates the
! quasi-particle energies $\varepsilon^{qp}_{n\vec{k}}$
! 

!
! !USES:
!
      use bands,      only:bande,bande0,nbandsgw,nspin,ibgw,nbgw,efermi,&
     &                     eferqp,eqp,eqp_im,nomaxs,numins,metallic,&
     &                     eferhf
      use kpoints,    only: nirkp,nvelgw
      use constants,  only: hev,czero, cone
      use selfenergy, only: sigx,sigc,selfx,vxcnn,selfc,znorm, &
     &                      iop_es,iop_ac,npar_ac,sacpar, &
     &                      sigm,qpwf_coef,vhmn,vh0mn,iop_dvh,   &
     &                      vxcmn,iop_qsgw,beta_sc,sxcmn,   &
     &                      dvmn_old,iop_gw0,isym_kbz
      use freq,       only: omega,nomeg
      use task,       only: casename,nmax_sc,eps_sc,fid_outqp
      implicit none 
      integer, intent(in) :: isxc    ! which exchange-correlation selfenergy
                                     !  0 -- GW, 1 -- HF, 2 -- COHSEX
      integer, intent(in) :: iod     ! whether or not to consider off-diaginal contributions 
                                     !  0 -- w/o offdigonal contributions (fixed wave functions)
                                     !  1 -- with offdiagonal contributions 
      integer, intent(in) :: isc_gw  ! stage of self-consistency in GW
                                     !   -1 for non-self-consistent

! !LOCAL VARIABLES:
      integer:: fout 
      real(8):: egap
      logical:: ldbg=.false.
      character(20):: sname="calceqp"
      logical:: lwarn = .false.

!  ! internal subroutines 
!    sub_eqp_gw
!    sub_eqp_ex
!    sub_eqp_chsx
!    sub_eqpod_static
!    sub_eqpod_gw

! !EXTERNAL ROUTINES: 
      external fermi

!
!EOP
!BOC

#ifdef DEBGU
      ldbg = .true.
#endif 
      write(6,*) trim(sname)//": calc QP energies"
!      call linmsg(6,'-',trim(sname))

      fout = fid_outqp

      if(iod.eq.0) then    !! consider only diagonal contributions, without wave functions self-consistency 
        if(isc_gw .ge. 0) iop_es = -1 
        select case(isxc) 
        case(0)       ! full GW 
          call sub_eqp_gw 
        case(1)       ! exchange only 
          call sub_eqp_ex
        case(2)       ! COHSEX 
          call sub_eqp_chsx 
        case(3)       ! PBE0 (not implemented yet) 
          call sub_eqp_sx
        case(4)      ! PBE0 (not implemented yet) 
          call sub_eqp_pbe0
        endselect  
      else  
        !! calculate QP energies with the contribution of off-diagonal contributions 
        call sub_eqpod_static   
      endif 

      call fermi(nirkp,nbandsgw,eqp,nvelgw,nspin,eferqp,egap,'QP')

      !* Print out the band gap at the HF level as a by-product
      if(iod.eq.0) then 
        call fermi(nirkp,nbandsgw,bande0(ibgw:nbgw,:,:)+selfx-vxcnn, &
     &           nvelgw,nspin,eferhf,egap,'HF')
      endif 

      contains 
!
!     calculate exchange-only quasi-particle energies 
!
        subroutine sub_eqp_ex
        implicit none 
        integer:: isp,irk,ie
        real(8):: delta
        call linmsg(6,'-','calceqp:exchange-only') 
        do isp=1,nspin
          do irk=1,nirkp
            do ie=ibgw,nbgw
              selfx(ie,irk,isp)=real(sigx(ie,irk,isp))
              delta=selfx(ie,irk,isp)-vxcnn(ie,irk,isp)
              eqp(ie,irk,isp)=bande(ie,irk,isp)+delta
            enddo ! ie
          enddo ! irk
        enddo ! isp
        end subroutine 

!
!     PBE0
!
        subroutine sub_eqp_pbe0
        implicit none
        integer:: isp,irk,ie
        real(8):: delta

        call linmsg(6,'-','calceqp:pbe0')
        do isp=1,nspin
          do irk=1,nirkp
            do ie=ibgw,nbgw
              selfx(ie,irk,isp)=real(sigx(ie,irk,isp))
              delta=0.25*(selfx(ie,irk,isp)-vxcnn(ie,irk,isp))
              eqp(ie,irk,isp)=bande(ie,irk,isp) + delta
            enddo ! ie
          enddo ! irk
        enddo ! isp
        end subroutine
!
!    COHSEX  
!
        subroutine sub_eqp_chsx
        implicit none 
        integer:: isp,irk,ie
        real(8):: delta

        do isp=1,nspin
          do irk=1,nirkp
            do ie=ibgw,nbgw
              selfx(ie,irk,isp)=real(sigx(ie,irk,isp))
              selfc(ie,irk,isp)=real(sigc(1,ie,irk,isp))

              sigc(2,ie,irk,isp)=sigc(2,ie,irk,isp)+sigx(ie,irk,isp) 
              delta=selfx(ie,irk,isp)+selfc(ie,irk,isp)-vxcnn(ie,irk,isp)
              eqp(ie,irk,isp)=bande(ie,irk,isp)+delta
            enddo ! ie
          enddo ! irk
        enddo ! isp
        end subroutine

!
!    SEX-only in COHSEX  
!
        subroutine sub_eqp_sx
        implicit none
        integer:: isp,irk,ie
        real(8):: delta

        do isp=1,nspin
          do irk=1,nirkp
            do ie=ibgw,nbgw
              selfx(ie,irk,isp)=real(sigx(ie,irk,isp))
              selfc(ie,irk,isp)=real(sigc(1,ie,irk,isp))

              delta=selfx(ie,irk,isp)+selfc(ie,irk,isp)-vxcnn(ie,irk,isp)
              eqp(ie,irk,isp)=bande(ie,irk,isp)+delta
            enddo ! ie
          enddo ! irk
        enddo ! isp
        end subroutine


!
!  Calculate GW quasi-particle energies for given diagonal self-energy matrix elements
!
        subroutine sub_eqp_gw
!       The critical parameters used in this subroutine is iop_es, 
!       which control how Fermi energy shift is treated
!       iop_es =    
!         -1 -- selfconsitent GW0  
!          0 -- perturbative G0W0 without energy shift
!          1 -- perturbative G0W0 with energy shift
!          2 -- iterative G0W0    with energy shift
!
        implicit none
        integer(4) :: ie,isp,irk,isc   !(Counter) index for band, spin,irreducible k-points, iteration
        integer(4) :: ierr
        integer(4) :: nvm,ncm    ! band index for valence band maximum and conduction band minimum 
        real(8) :: delta   ! absolute value of the difference between qp energies of succesive iterations
        real(8) :: egap,egap0,egap_old
        real(8) :: es
        real(8) :: enk,vxcnk,znk,sxnk,scnk,snk,enk0
        real(8) :: efshift
        real(8) :: zvm,zcm,dvm,dcm,evmks,evmhf,evmgw 
        complex(8) :: sig,ein
        complex(8) :: dsig

        real(8):: omg_ac(nomeg) ! Parameters of the selfenergya
        complex(8) ::sc_ac(nomeg)

        lwarn = .true.
        if(iop_es.ge.0) then 
          write(6,*) "# Parameters used:"
          write(6,*) "#  Analytic contiuation (iop_ac)  =", iop_ac
          write(6,*) "#  Fermi level shift    (iop_es)  =", iop_es
          write(6,*) "#  Nr.freq points  (nomeg)        =", nomeg
          write(6,*) "#  Number of AC poles (npar_ac/2) =", npar_ac/2
        endif 

        nvm=nomaxs(1)
        ncm=numins(1)

        eqp(ibgw:nbgw,1:nirkp,1:nspin)=bande(ibgw:nbgw,1:nirkp,1:nspin)
        eferqp = 0.0
        es=0.d0
        egap_old=100.0
        ierr=0
        do isc=0,nmax_sc 
          do isp=1,nspin      ! Loop over spin
            do irk=1,nirkp
              do ie=ibgw,nbgw
                enk0=bande0(ie,irk,isp)    ! original KS band energies 
                vxcnk=vxcnn(ie,irk,isp)   
                sxnk=real(sigx(ie,irk,isp))
                enk=bande(ie,irk,isp)

                !! enk is the energy at which the self-energy is calculated 

                if(iop_es.eq.-1) then  !! used for GW0 
                  ein=cmplx(enk,0.d0)

                else if(iop_es.eq.2)  then 

                  ein=cmplx(eqp(ie,irk,isp)-es,0.d0)

                elseif(iop_es.eq.3) then

                  ein=cmplx(eqp(ie,irk,isp),0.d0)
            
                else 
                  ein=cmplx(enk,0.d0) 
                endif 
             
                if(real(ein).ge.0.d0) then 
                  omg_ac = omega 
                  sc_ac = sigc(1:nomeg,ie,irk,isp)
                else 
                  omg_ac = - omega
                  sc_ac = conjg(sigc(1:nomeg,ie,irk,isp))
                endif 
          
                !! Calculate the new quasi-particle energy 
                call calcacfreq(0,iop_ac,nomeg,omg_ac,sc_ac,npar_ac, &
     &                  sacpar(:,ie,irk,isp),ein,sig,dsig)
                znk=1.0d0/(1.0d0-real(dsig))
                scnk=real(sig)
                snk=scnk+sxnk

                selfx(ie,irk,isp)=sxnk
                selfc(ie,irk,isp)=scnk
                znorm(ie,irk,isp)=znk

                if(znk.gt.1.d0 .or.znk.lt.0.d0) then
                  write(6,*) "WARNING: nonphysical Znk found!"
                  write(6,*) " -- check case.outqp for details "
                  write(fout,100) irk,ie,enk*hev,dsig,znk
                  znk=1.0
                  dsig = 0.d0
                endif

                select case(iop_es) 
                case (-1) !! one step in selfconsitent GW0
                  if(iop_gw0.eq.1) then ! direct iteration 
                    delta = snk-vxcnk + enk0- enk
                  elseif(iop_gw0.eq.2) then  ! linear expansion w.r.t. the KS energy
                    delta = znk*(snk-vxcnk) + enk0- enk 
                  else                       ! linear expansion w.r.t. the previous QP energy
                    delta = znk*(snk-vxcnk + enk0- enk)
                  endif
                case(0)
                  delta = znk*(snk-vxcnk)
                case(1) 
                  delta = znk*(snk-vxcnk-es*dsig)
!                  delta = znk*(snk-vxcnk)+(1.d0-znk)*es
                case (2) 
                  if(isc.eq.0 ) then
                    delta=znk*(snk-vxcnk)
                  else
                    delta=snk-vxcnk
                  endif
                case (3) 
                  delta = snk-vxcnk 
                endselect 
                eqp(ie,irk,isp) = enk + delta
                eqp_im(ie,irk,isp) = dimag(sig)*znk
              enddo ! ie
            enddo ! irk
          enddo ! isp

          !! iop_es <= 0 ==> selfconsitent GW0, no loop is needed 
          if( iop_es.le.0) exit 

          call fermi(nirkp,nbandsgw,eqp,nvelgw,nspin,eferqp,egap,'none')


          if(isc.eq.0) then 
            egap0=egap
            write(6,8) 
          endif 

          if(egap.gt.0.d0) then 
            write(6,9 ) isc,eferqp,egap*hev
          else 
            write(6,10) isc,eferqp,-egap
          endif 

          if(isc.ne.0.and..not.metallic) then 
            if(abs(egap-egap0).gt.0.2d0) then   
              write(6,*) "calceqp: WARNING --- Band gap deviates from &
     &initial GW gap more than 0.2 H"
              write(6,*) "   Most probably the iteration is not stable!"
              write(6,*) "   iteration will exit here"
              ierr=1
              exit 
            endif 
          endif 

          if( (metallic.and.abs(es-eferqp).le.eps_sc) &
     &      .or. (.not.metallic.and.abs(egap-egap_old).le.eps_sc)) exit 

          es=eferqp
          egap_old=egap
        enddo  !  do it

        if(isc.gt.nmax_sc) ierr=1 

        if(ierr.ne.0) then 
          write(6,*) "calceqp: WARNING --- Fail to converge!!! "
        endif 

   2    format(2i4,4f14.6)
   3    format(2x,'ikp ie',7x,'E_lda',9x,'S_x',11x,'V_xc',12x,'E_qp')
       
   8    format( ' #iter',7x,"Ef_QP",10x,"es",9x,"Eg/eV",2x,'Eg(k=0)/eV')
   9    format( 'isc=',i6,' insulating, E_F=',g12.4,' Eg/eV   =',g12.3)
  10    format( "isc=",i6,' metallic,   E_F=',g12.4,' DOS(E_F)=',g12.3)
 100    format('WARNING: nonphysical Znk',/,'irk=',i4,'  ie=',i4, &
     &  '  enk=',f8.3,' eV',/,' dsig=',2g12.4, '  Znk=',f8.3,/, &
     &  ' - reset to 1.0') 
        endsubroutine 

        subroutine sub_eqpod_static
!
!    calculate quasi-particle energies with the off-diagonal contributions, 
!    i.e. considering full self-energy matrix, within the static approximation. 
!    In this case, full self-energy matrix is stored in sigm. This can include 
!     - bare exchange (HF)
!     - hybrid functionals 
!     - COHSEX  
!     - full GW 
!
        implicit none
        integer:: isp,irk,iem,ien
        integer:: n,ierr,lwork

        real(8):: delta
        real(8),allocatable:: rwork(:),eig(:),eig0(:) 
        complex(8),allocatable:: work(:) 
        complex(8),allocatable:: hmat(:,:),dvmn(:,:),sxcmn(:,:) 

        n=nbgw-ibgw+1   !! the size of matrix is n*n
        lwork=2*n
        allocate(sxcmn(ibgw:nbgw,ibgw:nbgw),& 
     &           hmat(ibgw:nbgw,ibgw:nbgw), &  
     &           eig(ibgw:nbgw),            &
     &           dvmn(ibgw:nbgw,ibgw:nbgw),  &
     &           work(2*n),                 &
     &           rwork(3*n),stat=ierr)
        call errmsg(ierr.ne.0,sname,"Fail to allocat hmat etc.")

        do isp=1,nspin
          do irk=1,nirkp

            if(ldbg) write(6,*) "irk=",irk,"isp=",isp 

            !! first calculate exchange-only quasi-particle energies as
            ! a by-product of the full GW calculation
          
            call sub_herm_sxc(sxcmn,irk,isp,isxc)
            dvmn=sxcmn - vxcmn(:,:,irk,isp)

            if(iop_dvh.gt.0.and.isc_gw.gt.0) then 
              dvmn = dvmn + vhmn(:,:,irk,isp)-vh0mn(:,:,irk,isp)
            endif 

            !! set up the effective Hamiltonian matrix 
            hmat=0.d0
            do ien=ibgw,nbgw
              hmat(ien,ien)=bande0(ien,irk,isp) 
            enddo
            if(isc_gw.lt.5) then 
              hmat=hmat+dvmn
            else 
              hmat=hmat+beta_sc*dvmn   &   
     &               +(1.d0-beta_sc)*dvmn_old(:,:,irk,isp)
            endif 
            dvmn_old(:,:,irk,isp) = dvmn

            !! solve the eigen-problem 
            call zheev('v','u',n,hmat,n,eig,work,lwork,rwork,ierr)
            call errmsg0(ierr,sname,"fail to call zheev")

            eqp(ibgw:nbgw,irk,isp)=eig(ibgw:nbgw)
            qpwf_coef(ibgw:nbgw,ibgw:nbgw,irk,isp)=hmat(:,:)

          enddo ! irk
        enddo ! isp
        deallocate(hmat,dvmn,eig,work,rwork)

  100   format(6x,100i8)      !! index lies
  101   format(i6,100f8.3)    !! matrix elements
  102   format(i6,f10.4,4x,100f6.2)  !! solution 
  103   format(6x,100f10.4)   !! diagonal 
  104   format(20x,100i6)      !! index lies

        end subroutine

!
! This subroutine calculates an effective static and Hermitian Sxc matrix assuming 
! that the Sc_mn(i \omega) have been calculated 
! 
        subroutine sub_herm_sxc(sxc,irk,isp,isxc)
        implicit none 
        integer, intent(in):: irk, isp !! index for irreducible k-points and spin 
        complex(8), intent(out) :: sxc(ibgw:nbgw,ibgw:nbgw)
        integer, intent(in):: isxc 

        integer:: ie1, ie2, iom
        real(8):: enk,emk 
        complex(8):: ein,snn,dsig,smn,snm
        complex(8):: apar(npar_ac),sc_ac(nomeg)  
        complex(8),allocatable:: tmat(:,:),sc(:,:)
        real(8) :: omg_ac(nomeg) 

        allocate(sc(ibgw:nbgw,ibgw:nbgw))
        allocate(tmat(ibgw:nbgw,ibgw:nbgw))

        !! initialize Sxc 
        sxc=0.d0 
!
!       Bare exchange
!
        do ie2=ibgw,nbgw
          sxc(ie2,ie2)=real(sigm(ie2,ie2,0,irk,isp))
          do ie1=ibgw,ie2-1
            sxc(ie1,ie2)=0.5d0*( sigm(ie1,ie2,0,irk,isp) &
     &                          +conjg(sigm(ie2,ie1,0,irk,isp)) )
            sxc(ie2,ie1)=conjg(sxc(ie1,ie2))
          enddo
        enddo

        if(isxc.eq.1) return  

!
!       COHSEX correlation selfenergy
!
        if(isxc.eq.2) then 
          do ie2=ibgw,nbgw 
            sxc(ie2,ie2)=sxc(ie2,ie2)+real(sigm(ie2,ie2,1,irk,isp))
            do ie1=ibgw,ie2-1 
              smn = 0.5d0*(sigm(ie1,ie2,1,irk,isp) &
     &              + conjg(sigm(ie2,ie1,1,irk,isp)))
              sxc(ie1,ie2)=sxc(ie1,ie2)+smn 
              sxc(ie2,ie1)=sxc(ie2,ie1)+conjg(smn)
            enddo
          enddo

          return 

        endif 

!
!       GW correlation selfenergy in the FvSK scheme
!
        tmat=qpwf_coef(:,:,irk,isp)

        do iom=1,nomeg 
          call trans_vmn('g',sigm(:,:,iom,irk,isp),tmat,nbandsgw)
        enddo 

        do ie2=ibgw,nbgw 
          enk=bande(ie2,irk,isp)
          if(enk.ge. 0.d0) then  
            omg_ac=omega
            sc_ac=sigm(ie2,ie2,1:nomeg,irk,isp) 
          else 
            omg_ac=-omega 
            sc_ac=conjg(sigm(ie2,ie2,1:nomeg,irk,isp))
          endif 
        
          if(iop_qsgw.eq.0) then 
            ein=0.d0
          else
            ein=cmplx(enk,0.d0)
          endif 

          call calcacfreq(0,iop_ac,nomeg,omg_ac,sc_ac,npar_ac, &
     &                apar,ein,snn,dsig)    
          sc(ie2,ie2)=real(snn)

          do ie1=ibgw,ie2-1 
            emk=bande(ie1,irk,isp)

            if(enk.ge.0.d0) then 
              omg_ac=omega
              sc_ac=sigm(ie1,ie2,1:nomeg,irk,isp)
            else 
              omg_ac= - omega
              sc_ac=conjg(sigm(ie2,ie1,1:nomeg,irk,isp))
            endif 

            if(iop_qsgw.le.1) then 
              ein=0.d0 
            else 
              ein=cmplx(enk,0.d0) 
            endif 
              
            call calcacfreq(0,iop_ac,nomeg,omg_ac,sc_ac,npar_ac, &
     &                apar,ein,smn,dsig)

            if(emk.ge.0.d0) then
              omg_ac=omega
              sc_ac=sigm(ie2,ie1,1:nomeg,irk,isp)
            else
              omg_ac= - omega
              sc_ac=conjg(sigm(ie1,ie2,1:nomeg,irk,isp))
            endif

            if(iop_qsgw.eq.2) then 
              ein=cmplx(emk,0.d0)
            endif

            call calcacfreq(0,iop_ac,nomeg,omg_ac,sc_ac,npar_ac, &
     &                apar,ein,snm,dsig)
            smn=0.5d0*(smn+conjg(snm))
            sc(ie1,ie2)=smn
            sc(ie2,ie1)=conjg(smn)
          enddo
        enddo

        call trans_vmn('c',sc,tmat,nbandsgw)
        do ie2=ibgw,nbgw 
          sc(ie2,ie2)=real(sc(ie2,ie2))
          do ie1=ibgw,ie2-1
            sc(ie1,ie2)=0.5d0*(sc(ie1,ie2)+conjg(sc(ie2,ie1)))
            sc(ie2,ie1)=conjg(sc(ie1,ie2))
          enddo
        enddo

        sxc=sxc+sc
        deallocate(sc)
        deallocate(tmat)

        end subroutine  

      end subroutine calceqp

!EOC
          
