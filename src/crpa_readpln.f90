!BOP
!
! !Routine : crpa_readpln
! !INTERFACE:
      subroutine crpa_readpln()  

!
!! DESCRIPTION :

!  This subroutine reads the outcRPA file as obtained from the construction of Wannier functions 
!

!! USES:
      use constants, only: pi
      use crpa,    only: nlorb_wf, nbmax_wf,nbmin_wf,nlmorb,pln,nbk_wf,& 
     &                   info_orb,nbands_wf,iop_pln_phase,wf_centers,&
     &                   wf_centers_new, wf_centers_lm
      use bands,   only: nspin
      use kpoints, only: nirkp,nkp,get_kindex,kpirind
      use eigenvec,only: lsymvector 
      use struk,   only: pos 
      use task,    only: casename 
      implicit none 

!! LOCAL VARIABLES:
      integer:: nlorb_in, nk_in
      integer:: i,iorb,itmp,iat,isrt,l,nmorb,ityp,ikv(3),idv,ik_in,ik
      integer:: nk_wf
      integer:: fid=999 
      integer:: isp,ie,ib,nbk_bot,nbk_top
      integer:: ilm,im 
 
      integer:: ierr 
      real(8):: rtmp,p_re,p_im
      real(8):: eup, elow 

      real(8):: kvec(1:3),dR(3) 
      real(8):: arg
      complex(8):: phs

      character(len=10):: orb_typ
      character(len=120):: fn_pln

      character(20):: sname="crpa_readpln"
      logical:: ldbg=.false.

      call linmsg(6,'-',sname)

      fn_pln=trim(casename)//".pln"
      open(fid,file=trim(fn_pln),action='read',iostat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to open the outwf file " &
     &       //trim(fn_pln))

      read(fid,*) elow, eup   !! line 1: 
      read(fid,*) nk_in       !! line 2:  

      write(6,202) "The energy window for WF projectors:",elow,eup
      write(6,101) "The number of k-points",nk_in 

      nk_wf=nkp

      if(nk_in.ne.nk_wf) then 
        write(6,*) sname//"- ERROR: inconsistent nkp !!!"
        write(6,*) "internal nkp=",nkp
        write(6,*) "nkp from pln file=",nk_in 
        write(6,*) "lsymvector=",lsymvector
        stop
      endif 

      read(fid,*)               !! line 3 
      read(fid,*) nlorb_in      !! line 4 
      write(6,101) " Nr. of correlated orbitals (in terms of l) to be considered:",nlorb_in
      if(nlorb_wf.eq.0) then 
        nlorb_wf = nlorb_in
      else 
        if(nlorb_wf.ne.nlorb_in) then 
          write(6,*) "WARNING: inconsistent nlorb_wf"
          nlorb_wf = nlorb_in
          deallocate(wf_centers_new)
        endif 
      endif 

      allocate(info_orb(4,nlorb_wf)) 
      read(fid,*)               !! line 5 -> comment 
      read(fid,*)               !! line 6 (with the multiplicity of each sort of atoms) 
      read(fid,*)               !! line 7 -> comment 
      nlmorb = 0 
      do iorb=1,nlorb_wf           !! lines 8 - 7+nlorb: 
        read(fid,*) iat,isrt,itmp,l,nmorb,orb_typ
        info_orb(1,iorb) = iat
        info_orb(2,iorb) = l
        info_orb(3,iorb) = nmorb
        if(orb_typ.eq.'complex') then
          ityp = 0                    ! spherical complex
        elseif(orb_typ.eq.'cubic') then 
          if(l.eq.2.and.nmorb.eq.5) then 
            ityp = 1                    ! cubic harmonic 
          elseif(l.eq.2.and.nmorb.eq.2) then 
            ityp = 2                  ! d-eg 
          elseif(l.eq.2.and.nmorb.eq.3) then 
            ityp = 3                  ! d-t2g
          else
            write(6,*) "unknown type of wannier functions"
            ityp = 99
          endif 
        endif 
          
        info_orb(4,iorb) = ityp
        nlmorb = nlmorb + nmorb
      enddo !! iorb

      !! set the centers of Wannier functions
      allocate(wf_centers(3,nlorb_wf)) 
      do iorb=1,nlorb_wf
        iat = info_orb(1,iorb)
        wf_centers(:,iorb) = pos(:,iat)
      enddo

      write(6,*) "Centers of Wannier functions"
      write(6,300) "n","center_x","center_y","center_z"
      do iorb=1,nlorb_wf
        write(6,301) iorb,wf_centers(1:3,iorb)
      enddo

      if(iop_pln_phase.eq.1.and.allocated(wf_centers_new)) then 
        write(6,*) "Shifted centers of Wannier functions"
        write(6,300) "n","center_x","center_y","center_z"
        do iorb=1,nlorb_wf
          write(6,301) iorb,wf_centers_new(1:3,iorb)
        enddo
      endif

 300  format(a5,3a10)
 301  format(i5,3f10.3)

      read(fid,*)
      read(fid,*) nbmin_wf,nbmax_wf
      nbands_wf = nbmax_wf - nbmin_wf + 1
      write(6,102) "Lower/upper bound for WF bands:",nbmin_wf,nbmax_wf
      allocate( pln(1:nlmorb,nbmin_wf:nbmax_wf,1:nk_wf,1:nspin), &
    &           wf_centers_lm(3,nlmorb), & 
    &           nbk_wf(2,nk_wf,nspin), &
    &          stat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to allocate pln")

      pln = 0.d0
 
      read(fid,*)
      read(fid,*)
      do isp=1,nspin 
        do ik_in=1,nk_wf
          read(fid,*) ikv(1:3),idv,rtmp,nbk_bot,nbk_top
          ik = get_kindex(1,dble(ikv)/idv)
          call errmsg(ik.eq.0.or.ik.ne.ik_in,sname,"ERROR: inconsistent k-mesh!")

          kvec(1:3) = dble(ikv)/idv

          nbk_wf(1,ik,isp) = nbk_bot
          nbk_wf(2,ik,isp) = nbk_top

          do ib=nbk_bot,nbk_top

            ilm=0
            do iorb=1,nlorb_wf

              if(iop_pln_phase.eq.1.and.allocated(wf_centers_new)) then
                dR(1:3) = wf_centers_new(:,iorb)-wf_centers(:,iorb)
                arg = sum(dR*kvec(:))
                phs = cmplx(cos(2.0d0*pi*arg),-sin(2.0d0*pi*arg),8)
              else
                phs = 1.0
              endif

              do im=1,info_orb(3,iorb)
                ilm = ilm+1
                read(fid,*) itmp,isrt,itmp,itmp,p_re,p_im
                pln(ilm,ib,ik,isp)=cmplx(p_re,p_im,8)*phs
              enddo !!im
            enddo !!iorb 
          enddo !!ib
        enddo !!ik
      enddo !!isp
      close(fid)

      !! set the centers of each (lm)-resoluted Wannier function
      ilm = 0 
      do iorb=1,nlorb_wf

        if(iop_pln_phase.eq.1.and.allocated(wf_centers_new)) then
           wf_centers(:,iorb) = wf_centers_new(:,iorb)
        endif 

        do im=1,info_orb(3,iorb) 
          ilm = ilm + 1
          wf_centers_lm(:,ilm) = wf_centers(:,iorb) 
        enddo 

      enddo 

      !! check the unitarity of the pln matrix 
      call crpa_chk_unitary(nk_wf,nspin) 

      if(ldbg) then 
        write(6,*) 
        write(6,100) "WF projectors:"
        do ik=1,nk_wf
          write(6,*)
          write(6,101) "Wannier projections for ik=",ik
          do ib=nbmin_wf,nbmax_wf
            write(6,'(i5,100f10.4)') ib,pln(:,ib,ik,1)
          enddo
        enddo
      endif 

  100 format(a) 
  101 format(a,i5) 
  102 format(a,2i5) 
  105 format(a,5i5)
  202 format(a,2f8.3)
  203 format(a," n=",i5," ik=",i5,f8.3)
      end subroutine crpa_readpln
