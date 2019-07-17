!BOP
!
! !Routine : crpa_readw2w
! !INTERFACE:
      subroutine crpa_readw2w()  

!
!! DESCRIPTION :

!  This subroutine reads the Wannier function data from the output of 
!  wannier90, as interfaced to wien2k by wien2wannier (w2w) file as 
!  obtained from the construction of Wannier functions 
!

!! USES:
      use bands,   only: nspin
      use constants, only:pi,Bohr2Ans
      use crpa,    only: nbmax_wf,nbmin_wf,nlmorb,pln,nbk_wf,nlorb_wf,&
     &                   info_wf,info_orb,nbands_wf,iop_pln_phase,    &
     &                   wf_centers_lm,wf_centers,wf_centers_new,     &
     &                   tol_shift_centers   
      use kpoints, only: nirkp,nkp,get_kindex,kpirind,kvecs 
      use lapack,  only: la_invsym
      use eigenvec,only: lsymvector 
      use task,    only: casename 
      use struk,   only: inddf,pos,alat,ortho
      implicit none 

!! LOCAL VARIABLES:
      integer:: i,j,k,ik,ial,iw,il,ia,l,im,m,ilm,iorb
      integer:: nk_wf,nbas
      integer:: fid1=998,fid2=999 
      integer:: isp,ie,ib,nbk_bot,nbk_top

      integer:: nbands,nbands_excl,nntot
 
      integer:: ierr,itmp
      real(8):: rtmp,x1,x2

      real(8):: kvec(1:3)
      real(8):: arg
      complex(8):: phs

      real(8):: real_latt(3,3) 
      real(8):: recip_latt(3,3) 
      real(8),allocatable:: wf_spreads(:) 

      character(len=120)::stmp
      character(len=120)::fn_chk,fn_inwf,fn_chkdn

      complex(8):: ctmp
      complex(8),allocatable:: umat(:,:,:),umat_dis(:,:,:) 
      
      logical:: ldbg=.false.
      logical:: ldisentangled,ltmp
      character(20):: sname="crpa_readw2w"

      call boxmsg(6,'-',sname)

      if(nspin.eq.1) then 
        fn_chk=trim(casename)//".chk"
      else
        fn_chk=trim(casename)//".chkup"
        fn_chkdn = trim(casename)//".chkdn"
      endif 

      fn_inwf=trim(casename)//".inwf"
      call LinMsg(6,'*',"start crpa_readw2w")
      call LinMSG(6,'*',"Wannier projection from wannier90 output")

      !*
      !* first read information from case.w2win 
      !* 
      open(fid1,file=fn_inwf,action='read',iostat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to open "//trim(fn_inwf))

      open(fid2,file=fn_chk,form='unformatted',action='read',iostat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to open "//trim(fn_chk))

      read(fid1,*)           
      read(fid1,*) nbmin_wf,nbmax_wf 
      write(6,102) "Upper and lower bound ",nbmin_wf,nbmax_wf
      nbands_wf = nbmax_wf - nbmin_wf + 1  
      
      read(fid1,*) itmp, nlmorb
      write(6,101) "Number of WF's:",nlmorb

      read(fid2) stmp(1:33)  !! header 
      read(fid2) nbands
      read(fid2) nbands_excl
      read(fid2) (itmp,i=1,nbands_excl) 
      read(fid2) ((real_latt(i,j),i=1,3),j=1,3) 
      read(fid2) ((recip_latt(i,j),i=1,3),j=1,3)
 
      write(6,*) "Real Lattice Vectors:"
      do j=1,3
        write(6,'(3F12.4)') real_latt(1:3,j) 
      enddo 
      write(6,*) "Reciprocal Lattice Vectors:"
      do j=1,3
        write(6,'(3F12.4)') recip_latt(1:3,j) 
      enddo 

      read(fid2) nk_wf
      read(fid2) (itmp,i=1,3) 
      read(fid2) (rtmp,i=1,nk_wf*3) 
      read(fid2) nntot
      read(fid2) itmp
      read(fid2) stmp(1:20)

      write(6,101) "Number of k-points:",nk_wf
      write(6,101) "Number of bands:",nbands

      allocate(pln(1:nlmorb,nbmin_wf:nbmax_wf,1:nk_wf,1:nspin), &
     &         nbk_wf(2,nk_wf,nspin), &
     &         info_wf(4,1:nlmorb),      &
     &         wf_centers_lm(3,nlmorb), & 
     &         wf_spreads(nlmorb),   & 
     &         stat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to allocate pln etc.")

      isp=1
      pln = 0.d0 
      nbk_wf(1,:,:) = nbmin_wf
      nbk_wf(2,:,:) = nbmax_wf

      !! try to extract the information that defines the correlated
      !! orbitals in terms of "l" only 
      !! -- determine nlorb_wf by counting how many different ("ia","l") in
      !!    the WN projectors. Assumes that WF belonging to the same group
      !!    always comes together 
      nlorb_wf = 1

      im = 1
      ial = 1
      do iw=1,nlmorb
        read(fid1,*) nbas 
        do i=1,nbas 
          read(fid1,*) ia,l,m,x1,x2
        enddo  
        if(iw.gt.1) then 
          if(l.ne.info_wf(1,iw-1).or.ia.ne.info_wf(3,iw-1))then
            nlorb_wf = nlorb_wf + 1
            im = 1 
          else 
            im = im + 1 
          endif 
        endif 
        info_wf(1,iw) = l
        info_wf(2,iw) = im
        info_wf(3,iw) = ia
        info_wf(4,iw) = inddf(ia)
        write(6,105) "iw,atmlmorb",iw,info_wf(:,iw)
      enddo 
      write(6,101) "The number of WF groups:",nlorb_wf
    
      allocate(info_orb(4,nlorb_wf),wf_centers(3,nlorb_wf))

      info_orb(1,1) = info_wf(3,1)
      info_orb(2,1) = info_wf(1,1)
      info_orb(3,1) = 1

      il=1
      do iw=2,nlmorb
        if(info_wf(1,iw).ne.info_wf(1,iw-1).or. &
     &     info_wf(3,iw).ne.info_wf(3,iw-1)) then
          il = il + 1
          info_orb(1,il) = info_wf(3,iw)
          info_orb(2,il) = info_wf(1,iw)
          info_orb(3,il) = 1
        else
          info_orb(3,il) = info_orb(3,il) + 1
        endif 
      enddo

      !! check the consistency
      allocate(umat(nlmorb,nlmorb,nk_wf)) 

      read(fid2) ldisentangled
      if(ldisentangled) then 
        allocate(umat_dis(nbands,nlmorb,nk_wf))
        read(fid2) rtmp
        read(fid2) (ltmp,i=1,nbands*nk_wf)
        read(fid2) (itmp,i=1,nk_wf) 
        read(fid2) (((umat_dis(ib,iw,ik),ib=1,nbands),iw=1,nlmorb),ik=1,nk_wf)
      endif 
      
      read(fid2) (((umat(ib,iw,ik),ib=1,nlmorb),iw=1,nlmorb),ik=1,nk_wf)
      read(fid2) ((((ctmp,i=1,nlmorb),j=1,nlmorb),k=1,nntot),l=1,nk_wf)
      read(fid2) ((wf_centers_lm(i,iw),i=1,3),iw=1,nlmorb)
      read(fid2) (wf_spreads(iw),iw=1,nlmorb)

      do ik=1,nk_wf 
        if(ldisentangled) then 
          call zgemm('c','c',nlmorb,nbands,nlmorb,1.d0,umat(:,:,ik), &
     &     nlmorb,umat_dis(:,:,ik),nbands,0.d0,pln(:,:,ik,isp),nlmorb)
        else
          pln(:,:,ik,isp)=conjg(transpose(umat(:,:,ik)))
        endif 
      enddo
      close(fid1) 
      close(fid2)
      call sub_shift_centers(isp) 

      if(nspin.eq.2) then 
        isp = 2
        write(6,*) "Read spin-down Wannier projectors"
        open(fid2,file=fn_chkdn,form='unformatted',action='read',iostat=ierr)
        call errmsg(ierr.ne.0,sname,"Fail to open "//trim(fn_chk))
        read(fid2) stmp(1:33)  !! header 
        read(fid2) nbands
        read(fid2) nbands_excl
        read(fid2) (itmp,i=1,nbands_excl)
        read(fid2) (rtmp,i=1,9)
        read(fid2) (rtmp,i=1,9)
        read(fid2) nk_wf
        read(fid2) (itmp,i=1,3)
        read(fid2) (rtmp,i=1,nk_wf*3)
        read(fid2) nntot
        read(fid2) itmp
        read(fid2) stmp(1:20)
        read(fid2) ldisentangled
        if(ldisentangled) then
          read(fid2) rtmp
          read(fid2) (ltmp,i=1,nbands*nk_wf)
          read(fid2) (itmp,i=1,nk_wf)
          read(fid2) (((umat_dis(ib,iw,ik),ib=1,nbands),iw=1,nlmorb),ik=1,nk_wf)
        endif

        read(fid2) (((umat(ib,iw,ik),ib=1,nlmorb),iw=1,nlmorb),ik=1,nk_wf)
        read(fid2) ((((ctmp,i=1,nlmorb),j=1,nlmorb),k=1,nntot),l=1,nk_wf)
        read(fid2) ((wf_centers_lm(i,iw),i=1,3),iw=1,nlmorb)
        read(fid2) (wf_spreads(iw),iw=1,nlmorb)

        write(6,*)
        do ik=1,nk_wf
          if(ldisentangled) then
            call zgemm('c','c',nlmorb,nbands,nlmorb,1.d0,umat(:,:,ik), &
     &       nlmorb,umat_dis(:,:,ik),nbands,0.d0,pln(:,:,ik,isp),nlmorb)
          else
            pln(:,:,ik,isp)=conjg(transpose(umat(:,:,ik)))
          endif
        enddo
        close(fid2) 
        call sub_shift_centers(isp) 
      endif 

     !! check the unitarity of the pln matrix 
      call crpa_chk_unitary(nk_wf,nspin)

      if(ldbg) then
        do isp=1,nspin 
          if(nspin.eq.2) then 
            write(6,101) "Wannier projectors for ispin = ",isp
          else
            write(6,100) "WF projectors:"
          endif 

          do ik=1,nk_wf 
            write(6,*) 
            write(6,101) "ik=",ik
            do ib=nbmin_wf,nbmax_wf
              write(6,400) ib,pln(:,ib,ik,isp)
            enddo 
          enddo 
        enddo 
      endif 

      deallocate(umat)
      if(ldisentangled) deallocate(umat_dis)
      call LinMsg(6,'*',"end of crpa_readw2w")
 100  format(a)
 101  format(A,i5) 
 102  format(A,2i5)
 105  format(A,5i5)
 203  format(a," n=",i5," ik=",i5,f8.3)
 400  format(i5,100f10.6) 

      contains 

        subroutine sub_shift_centers(ispin)
        integer,intent(in):: ispin
        integer:: iat
        integer:: iorb
        integer:: iw 

        write(6,*) "Tol. to shift Wannier centers:",tol_shift_centers
        write(6,*) "Original Centers of Wannier functions"
        write(6,300) "n","center_x","center_y","center_z","spread"
        do iw=1,nlmorb
          write(6,301) iw,wf_centers_lm(1:3,iw),wf_spreads(iw)
        enddo
        call sub_convert_centers
        write(6,*) " --> in WIEN2k internal coordinates"
        do iw=1,nlmorb
          write(6,303) iw,wf_centers_lm(1:3,iw)
        enddo

        if(iop_pln_phase.eq.1) then  !! shift the Wannier centers 
          do iorb=1,nlorb_wf 
            write(6,*) iorb
            if(allocated(wf_centers_new)) then 
              wf_centers(1:3,iorb) = wf_centers_new(1:3,iorb) 
            else
              iat = info_orb(1,iorb)
              wf_centers(:,iorb) = pos(:,iat)
            endif
          enddo

          do ik=1,nk_wf
            do ib=nbmin_wf,nbmax_wf
              ilm = 0
              do iorb=1,nlorb_wf
                !! 
                do im=1,info_orb(3,iorb)
                  ilm = ilm+1
                  call sub_set_phase(phs,ilm,iorb,ik)
                  pln(ilm,ib,ik,isp) = pln(ilm,ib,ik,isp)*phs
                enddo
              enddo !!iorb
            enddo !!ib
          enddo !!ik

          ilm = 0  
          do iorb=1,nlorb_wf
            do im=1,info_orb(3,iorb)
              ilm = ilm+1
              wf_centers_lm(:,ilm) = wf_centers(:,iorb) 
            enddo
          enddo

          write(6,*) "Shifted Wannier Centers (in WIEN2k format)"
          do iw=1,nlmorb
            write(6,303) iw,wf_centers_lm(1:3,iw)
          enddo

        else
          ilm=0 
          do iorb=1,nlorb_wf 
            do im=1,info_orb(3,iorb)
              ilm = ilm+1
              if(im.eq.1) then
                wf_centers(:,iorb)=wf_centers_lm(:,ilm) 
              else
                if(sum(abs(wf_centers_lm(:,ilm)-wf_centers(:,iorb)))/3 &
     &             .gt.1.e-4) then 
                  write(6,*) "WARNING: abnormal Wann. centers detected!"
                endif 
              endif 
            enddo
          enddo
        endif

 300    format(a5,4a10)
 301    format(i5,4f10.3)
 303    format(i5,3f10.3)
        end subroutine 

        subroutine sub_set_phase(phs,ilm,iorb,ik) 
        complex(8),intent(out):: phs
        integer,intent(in):: ilm, iorb, ik
        real(8):: kvec(3), R0(3), arg

        kvec(1:3) = kvecs(:,ik)
        if(iop_pln_phase.eq.1) then
          R0(1:3) = wf_centers(:,iorb)-wf_centers_lm(:,ilm)
          if(ik.eq.1) write(6,500) ilm,R0(1:3) 

          arg = sum(R0*kvec(:))
          phs = cmplx(cos(2.0d0*pi*arg),sin(2.0d0*pi*arg),8)
        else
          phs = 1.0
        endif
  500   format("Wannier center shift for ilm=",i3,3F6.2) 

        end subroutine



        subroutine sub_convert_centers
!
!       Convert the Wannier centers to wien2k-consistent internal
!       coordinates 
!
        real(8):: a(3,3), r(3) 
        real(8):: x,x_new
        integer:: iw, i

        wf_centers_lm = wf_centers_lm/Bohr2Ans

        if(ortho) then !! if the structure is orthogonal including F or B symmetry
          do iw=1,nlmorb
            r(1:3) = wf_centers_lm(1:3,iw)/alat(1:3) 
            wf_centers_lm(1:3,iw) = r(1:3) 
          enddo 
        else   
          a = real_latt
          call la_invsym(a,3)
          do iw=1,nlmorb
            r = matmul(a,wf_centers_lm(:,iw)) 
            wf_centers_lm(:,iw) = r
          enddo 
        endif

        write(6,600) "ilm","Original","Rounded" 

        do iw=1,nlmorb
          r = wf_centers_lm(1:3,iw)
          do i=1,3 
            x = wf_centers_lm(i,iw) 
            call sub_round_up(x,x_new) 
            wf_centers_lm(i,iw) = x_new
          enddo
          write(6,601) iw,r(1:3),wf_centers_lm(1:3,iw) 
        enddo 
 600    format(a5,a18,4x,a18)
 601    format(i5,3f6.3,4x,3f6.3) 
        end subroutine 

        subroutine sub_round_up(x,x_new) 
        real(8),intent(in)::x
        real(8),intent(out)::x_new
        integer:: i,n_x0
        real(8):: x0,x0_list(20) 
        data x0_list(1:13) /0.0, 1.0, 0.5, 0.25, 0.75, &      !! # 5
     &                      0.333333, 0.666667, 0.166667, 0.833333, & !+ 4
     &                      0.125, 0.375, 0.625, 0.875/             ! + 4 
        n_x0 = 13
        x_new = x  
        do i=1,n_x0
          x0 = x0_list(i) 
          if(abs ( abs(x) - x0) < tol_shift_centers ) then 
            x_new = dsign(x0, x) 
            exit 
          endif

        enddo 
        end subroutine 
      end subroutine crpa_readw2w
