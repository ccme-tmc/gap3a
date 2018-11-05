!BOP
!
! !ROUTINE: bz_calcqdepw
!
! !INTERFACE:
      subroutine bz_calcqdepw(iq)
      
! !DESCRIPTION:
!
! This subroutine calculates the weights for q dependent BZ integrations.
! \textbf{needs checking!!}
!
!
! !USES:
!
      use bands,       only: bande,efermi,nbmaxpol,&
     &                       nspin,numins,nomaxs
      use bzint,       only: iop_bzintq,omgmax_ht,nomg_ht,eta_freq,&
     &                       bzint_smear,freq_factor, &
     &                       ntet,link,tvol,tnodes,wtet
      use bzinteg,     only: kcw,kwt_bz 
      use core,        only: eigcore,ncg_p,corind
      use freq,        only: nomeg, omega, iop_freq
      use kpoints,     only: nkp,kqid,kpirind,idikp
      use struk,       only: nat
      use task,        only: casename,savdir 
!
! !INPUT PARAMETERS:

      implicit none
      integer,intent(in) :: iq  !ID. number of the q-point
      
! !LOCAL VARIABLES:
      integer :: fout
      integer :: ib      ! counter, run over bands
      integer :: icg,ic  ! counter, run over core states
      integer :: isp     ! Index for spin 
      integer :: jb  ! counter, run over bands
      integer :: iom ! counter, run over frequencies
      integer :: ik,iat,irk,jk,jrk
      integer :: ncbm,nvbm,nmax
      
      integer, allocatable :: lt(:) ! index of the q-linked tetrahedron 
      integer, allocatable :: lkq(:) ! index of the q-linked kpoint 
      
      real(8) :: emaxb ! maximum energy of the second band
      real(8) :: edif, edsq,omsq

      
      real(8), allocatable :: enk(:,:) ! Pair of bands for which the weights are being calculated
      real(8), allocatable :: cwpar(:,:,:,:)   ! weight between the two bands
      character(len=10),external:: int2str
 
      logical :: ldbg=.false. 
      character(len=20):: sname="bz_calcqdepw"
      character(len=80):: fn_kcw
      real(8) :: tstart,tend 

!
!EOP
!BOC
!

!
!     Set lt to the linkd tetrahedron index
!

      call cpu_time(tstart) 

      allocate(lt(ntet),lkq(nkp))
      lt(1:ntet)=link(1:ntet,iq)  
      lkq(1:nkp)=kqid(1:nkp,iq)

      nmax = ncg_p + nbmaxpol 
      allocate(enk(nmax,nkp),cwpar(nmax,nmax,nkp,2))

      do isp=1,nspin 
        ncbm = numins(isp)
        nvbm = nomaxs(isp) 

        !! prepare a temporal array to store band eneriges including the
        !! core states 
        do icg=1,ncg_p
          iat=corind(1,icg)
          ic=corind(3,icg)
          enk(icg,1:nkp)=eigcore(ic,iat,isp)
        enddo 

        do ik=1,nkp
          irk=kpirind(ik)
          enk(ncg_p+1:nmax,ik) = bande(1:nbmaxpol,irk,isp) 
        enddo 
        
        if(iop_bzintq.le.0) then  
          !! calculate q-BZ weights using the generalized tetrahedron method  
          call sub_bzintq_0

        elseif(iop_bzintq.eq.1) then  !! using the direct summation  

          call sub_bzintq_1
        else 
          call errmsg(.true.,sname,"unsupported option for iop_bzintq")
        endif  
      enddo ! isp

      if(ldbg) then 
        fout=999
        fn_kcw=trim(savdir)//trim(casename)//".kcw-q"//trim(int2str(iq))
        open(fout,file=fn_kcw,action='write')
        write(fout,'(4g16.10)') kcw 
      endif 

      deallocate(lt,lkq,enk,cwpar)

      call cpu_time(tend) 
      if(iq.eq.1) call write_cputime(tend-tstart,sname) 

      contains 

        subroutine sub_bzintq_0
!
!       Internal subroutine to calculate q-BZ weights analytically 
!       (iop_bzintq=0) or numerically (iop_bzintq=1) 
!
        implicit none 
        integer:: io,iom 
        real(8):: coef,omg2,eta
        complex(8):: comg
        real(8),allocatable:: w2_ht(:,:,:,:),omg_ht(:) 
        real,parameter:: pi = 3.14159265358979d0

        if(iop_bzintq.eq.0 .or. iop_bzintq.eq.-1 ) then 
          do iom=1,nomeg
            call tetcw(nkp,ntet,nmax,wtet,enk,tnodes,lt,lkq,tvol, &
     &           efermi,omega(iom),iop_freq,cwpar(:,:,:,1))

            if(iop_freq.eq.2) then  !! calculate the imaginary part for real freq 
              call tetcw(nkp,ntet,nmax,wtet,enk,tnodes,lt,lkq,tvol,  &
     &           efermi,omega(iom),4,cwpar(:,:,:,2))
            endif

            if(iop_freq.eq.2) then   !! real freq
              kcw(1:ncg_p+nvbm,ncbm:nbmaxpol,1:nkp,iom,isp)= &
     &          cmplx(cwpar(1:ncg_p+nvbm,ncg_p+ncbm:nmax,1:nkp,1),&
                      cwpar(1:ncg_p+nvbm,ncg_p+ncbm:nmax,1:nkp,2))
            else        !! imag freq 
              kcw(1:ncg_p+nvbm,ncbm:nbmaxpol,1:nkp,iom,isp)= & 
     &          cmplx(cwpar(1:ncg_p+nvbm,ncg_p+ncbm:nmax,1:nkp,1),0.0)
            endif
          enddo ! iom

        else if (iop_bzintq.eq.-2) then 

          ! calculate q-BZ weights using the Hilbert transform (HT) approach
          eta=eta_freq 
          if(iop_freq.eq.3) eta = 0.d0 

          allocate( w2_ht(1:ncg_p+nvbm, ncbm:nbmaxpol, nkp, nomg_ht),&
     &            omg_ht(nomg_ht) )
          coef = omgmax_ht/nomg_ht*(2.0/pi) 
          do io=1,nomg_ht 
            omg_ht(io)=io*omgmax_ht/nomg_ht 
            call tetcw(nkp,ntet,nmax,wtet,enk,tnodes,lt,lkq,tvol, &
     &           efermi,omg_ht(io),4,cwpar(:,:,:,1))
            w2_ht(:,:,:,io) = cwpar(1:ncg_p+nvbm,ncg_p+ncbm:nmax,1:nkp,1)
          enddo 

          do iom=1,nomeg 
            omg2 = omega(iom)**2
            if(iop_freq.eq.3) omg2 = - omg2

            kcw(:,:,:,iom,isp) = 0.d0  
            do io=1, nomg_ht 
              comg = cmplx(omg_ht(io),-eta) 
              kcw(:,:,:,iom,isp) = kcw(:,:,:,iom,isp) &
     &          + w2_ht(:,:,:,io)*(coef*comg/(comg**2-omg2))
            enddo
          enddo   
          deallocate(w2_ht,omg_ht)
        endif 
        end subroutine 

        subroutine sub_bzintq_1
!
!       Internal subroutine to calculate q-BZ weights using the simple summation 
!       with finite-temperature smearing and life-time broadening 
!
        implicit none 
        real(8):: e_i, e_j, fnm,omg
        complex(8):: kw

        do iom=1,nomeg 

          omg = omega(iom) 
          do ik=1,nkp
            jk=kqid(ik,iq)

            do jb = ncbm,nbmaxpol
              do ib= 1,nvbm+ncg_p
                e_i = enk(ib,ik)-efermi
                e_j = enk(jb+ncg_p,jk)-efermi
                edif = e_j - e_i
                if( abs(edif) < 1.d-10 ) then 
                  kw = 0.d0 
                else
                  fnm = bzint_smear(0,e_i)*(1.0-bzint_smear(0,e_j))
                  kw = fnm*freq_factor(iop_freq,omg,edif)/nkp
                endif 

                kcw(ib,jb,ik,iom,isp)=kw

              enddo
            enddo
          enddo
        enddo
        end subroutine 

      end subroutine bz_calcqdepw
!EOC          
         
