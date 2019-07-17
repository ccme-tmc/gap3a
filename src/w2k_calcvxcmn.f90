!BOP
!
! !ROUTINE: w2k_calcvxcmn
!
! !INTERFACE:
      subroutine w2k_calcvxcmn(isym) 

! !DESCRIPTION:
!
!This subroutine calculates the matrix elements $v^{xc}_{nn}(\vec{k})$ 
!of equation \ref{vxc-10}, only for valence states.
!   Internal subroutines: 
!      sub_calcvxcmn_mt
!      sub_calcvxcmn_is
!      sub_calcvorbmn 
!
! !USES:
      use xcpot,       only:uxcu,lvorb,lorb,dmorb,lxcm,lmxc,istpw, &
     &                  vorbnn,vxc_hyb,natorb,iatorb,nlorb,vorb,ksxc,&
     &                  vxclm,lhybrid,vxcs,nksxc,uiorb
      use bands,       only: nv,ibgw,nbgw,nbandsgw,nspin
      use bzinteg,     only: kiw 
      use constants,   only: cone, czero, imag
      use eigenvec,    only: alfa, beta, gama, zzk
      use kpoints,     only: nirkp, nkp, idikp,wkir,kpirind,get_kvec
      use lapwlo,      only: lmax,lomax,nlo_at 
      use recipvec,    only: npw,gindex,ig0,indgkir,indgk, maxngk
      use selfenergy,  only: vxcmn
      use struk,       only: iatnr, mult,nat
      use task,        only: casename,time_lapack,time_vxc,lrestart, &
     &                       savdir
      implicit none
      integer, intent(in):: isym !  0/1 -- calculate sigx in the full/irreducible BZ

! !LOCAL VARIABLES:
      integer:: ik,irk      ! the index of full and irred. k-points
      integer:: iik, nktot
      integer:: iat,idf,ieq 
      integer:: ierr
      integer:: isp     ! index for spin 
      integer:: ien,iem ! Counter: run over the eigenvalues in the corresponding k-point.
      integer:: fout=999
      integer:: ikvec(3) ! Indexes of G_1+G'-G
      real(8):: kvec(3)
      real(8):: tstart,tend,time1,time2
      complex(8),allocatable:: vxc_mt(:,:),vxc_is(:,:),vxc_orb(:,:) 

      complex(8),external:: zdotc
      real(8),external:: getcgcoef

      character(20) :: sname='w2k_calcvxcmn'
      character(120):: msg,fn
      logical:: ldbg = .true.
!

! !REVISION HISTORY:
! 
! Created  Nov. 5, 2009 by H. Jiang
!EOP
!
!BOC

!
!  Calculate the integral between stars and planwaves
!
      call cpu_time(tstart)

      if(lrestart) then 
        call io_vxcmn('r','f',isym,ierr)        !! read vxcmn from the file
        if(ierr.eq.0) then
          write(6,*) "Restart mode: read vxcmn from files"
          return
        endif
      endif

      call intstipw

      if(isym.eq.0) then
        nktot = nkp
      else
        nktot = nirkp
      endif

      allocate(vxc_mt(ibgw:nbgw,ibgw:nbgw), &
     &         vxc_is(ibgw:nbgw,ibgw:nbgw), &
     &         vxc_orb(ibgw:nbgw,ibgw:nbgw), &
     &         stat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to allocate vxc_mt/is!") 

      if(ldbg) then 
        fn=trim(savdir)//trim(casename)//".outvxc"
        open(fout,file=fn,action='write') 
      endif 

      do isp=1,nspin 

        do iik=1, nktot                      !* Loop over the k-points:
          vxc_mt = 0.d0
          vxc_is = 0.d0
          vxc_orb = 0.d0

          ! set ik, irk, jk and jrk 
          if(isym.eq.0) then
            ik=iik
            irk=kpirind(ik)
          else
            irk=iik
            ik=idikp(irk)
          endif

          call readvector(ik,1,isp,0)        !* Read the eigenvector corresponding to the k-point ikp
          call expand_evec(ik,1,.true.,isp)  !* Calculate the expansion coeficients of the eigenvectors

          !! Contributiosn from Muffin-Tin sphere        
          idf=0
          do iat = 1, nat                       !*  Loop over inequivalent atoms:
            do ieq = 1, mult(iat)               !* Loop over equivalent atoms:
              idf = idf + 1
              call  sub_calcvxcmn_mt(iat,idf)

              if(lvorb)  then !! contribution from LDA+U correction
                call sub_calcvorbmn(iat,idf)
              endif 

            enddo ! ieq
          enddo ! iat
         
          !!Interstitial region:
          call sub_calcvxcmn_is
          
          !! sum up different contributions 
          vxcmn(:,:,iik,isp) = vxc_mt + vxc_is + vxc_orb

          !! Debugging mode: output different contributions 
          if(ldbg) then 
            call get_kvec(isym,iik,ik,irk,kvec)
            write(fout,10) iik,kvec(1:3)
            do ien=ibgw,nbgw
              do iem=ibgw,nbgw 
                write(fout,11) iem,ien,vxc_mt(iem,ien),vxc_is(iem,ien),&
     &           vxcmn(iem,ien,iik,isp) 
              enddo
            enddo 
          endif 

          !! force Hermiticity
          do ien=ibgw,nbgw 
            vxcmn(ien,ien,iik,isp)=real(vxcmn(ien,ien,iik,isp))
            do iem=ibgw,ien-1
              vxcmn(iem,ien,iik,isp)=0.5d0*( vxcmn(iem,ien,iik,isp) & 
     &           + conjg(vxcmn(ien,iem,iik,isp)))
              vxcmn(ien,iem,iik,isp)=conjg(vxcmn(iem,ien,iik,isp))
            enddo 
          enddo    
          
        enddo ! iik
      enddo  ! isp 
      if(ldbg) then 
        close(fout) 
      endif 

      deallocate(vxc_mt,vxc_is,vxc_orb) 
      deallocate(vxclm,vxcs,lmxc,lxcm,ksxc,uxcu) 

      call io_vxcmn('w','f',isym,ierr)       !! write vxcmn to the file

      call cpu_time(tend)
      time_vxc = time_vxc + tend - tstart 
          
    1 format(3e16.6,2i6)
    2 format(2i6,2e24.12)     
   10 format(/,i6,3e16.6,"  # ik,kvec(1:3)")
   11 format(2i6,6f12.6)

      contains 

        subroutine sub_calcvxcmn_mt(iat,idf)
        implicit none 
        integer,intent(in):: iat, idf 

        integer :: bl      ! Angular quantum number of the   mixed atomic functions (big l)
        integer :: bm      ! Magnetic quantum number of the mixed atomic function (big m)
        integer :: ir1, ir2, nr1, nr2, ir12 
        integer :: l1      ! l-quantum number of the (L)APW eigenfunction  at k
        integer :: l2      ! l-quantum number of the (L)APW eigenfunction  at k' = k - q
        integer :: l12
        integer :: lxc     ! (Counter) l index of the xc component
        integer :: m1      ! m-quantum number of the (L)APW eigenfunction  at k
        integer :: m2      ! m-quantum number of the (L)APW eigenfunction   at k' = k - q
        integer :: l1m1    ! Counter: Runs over MTS (L)APW basis  functions (l1,m1)
        integer :: l2m2    ! Counter: Runs over MTS (L)APW basis functions l2,m2)
        integer :: isgn    ! index for the sign of bm
        integer :: ir,jr   ! index for the radial functons with same l
        integer :: ijr     ! accumulative index for (ir,jr) as used in uxcu
        complex(8) :: angint   ! Angular integral calculated by fourylm
        complex(8) :: io1
        complex(8) :: af,caf

        do l2=0,lmax                     !*    loop over l2:
          nr2 = 2 
          if (l2.le.lomax) nr2 = 2 + nLO_at(l2,iat) 

          do l1=0,lmax                      !* Loop over l1
            l12 = l1+1 + l2*(lmax+1)

            nr1 = 2 
            if (l1.le.lomax) nr1 = 2 + nLO_at(l1,iat) 

            do m1=-l1,l1                  !* loop over m1
              l1m1 = l1*l1+l1+m1+1

              do lxc = 1, lxcm(iat)             !*  Loop over vxc components:
                bl=abs(lmxc(1,lxc,iat))
                bm=lmxc(2,lxc,iat)
                io1=cone
                if(lmxc(1,lxc,iat).lt.0) io1=imag

                if(.not.(abs(bl-l2).le.l1 .and.l1.le.bl+l2)) cycle 
                do isgn = 1, -1, -2
                  bm = bm*isgn
                  m2 = m1-bm

                  if(abs(m2).gt.l2) cycle 
                  
                  l2m2 = l2*l2 + l2 + m2 + 1

                  ! Calculate the angular integral:
                  if(iatnr(iat).gt.0)then
                    select case (bl)
                    case(3,7,9)
                      if(bm.lt.0)then
                        angint=-imag*getcgcoef(bl,l2,l1,bm,m2)
                      else  
                        angint=imag*getcgcoef(bl,l2,l1,bm,m2)
                      endif  
                    case default  
                      angint=getcgcoef(bl,l2,l1,bm,m2)
                    end select     
                  else
                    if(bm.lt.0)then
                      angint=io1*getcgcoef(bl,l2,l1,bm,m2)/sqrt(2.0d0)
                    elseif(bm.eq.0)then
                      angint=getcgcoef(bl,l2,l1,bm,m2)
                    else
                      angint=io1*getcgcoef(bl,l2,l1,bm,m2)        &
     &                 *(-1.0)**bm*sign(1,lmxc(1,lxc,iat))/sqrt(2.0d0)
                    endif
                  endif 
      
                  if(abs(angint).lt.1.0d-12) cycle 

                  do ien=ibgw,nbgw
                    do iem=ibgw,nbgw 

                      ir12 = 0 
                      do ir2 = 1, nr2  
                        if(ir2.eq.1) then 
                          af = alfa(ien,l2m2,idf)
                        elseif(ir2.eq.2) then 
                          af = beta(ien,l2m2,idf)
                        else 
                          af = gama(ien,ir2-2,l2m2,idf)
                        endif 

                        do ir1 =1, nr1 
                          ir12 = ir12 + 1

                          if(ir1.eq.1) then  
                            caf=conjg(alfa(iem,l1m1,idf))
                          elseif(ir1.eq.2) then 
                            caf=conjg(beta(iem,l1m1,idf))
                          else
                            caf=conjg(gama(iem,ir1-2,l1m1,idf))
                          endif 

                          vxc_mt(iem,ien)=vxc_mt(iem,ien) &
     &                     + caf*af*uxcu(ir12,l12,lxc,iat,isp)*angint
                        enddo ! ir
                      enddo ! jr
                    enddo !ien 
                  enddo ! iem 

                  if(bm.eq.0) exit 

                enddo ! isgn
              enddo ! lxc
            enddo ! m1     
          enddo ! l1
        enddo ! l2
        endsubroutine sub_calcvxcmn_mt    

        subroutine sub_calcvxcmn_is
        implicit none 

        integer :: i,iv1,iv2 ! Counter: run over the (L)APW basis functions (excluding core states).
        integer :: ig1,ig2
        integer :: nvk 
        integer :: ikxc    ! Counter: Runs over IPW's
        integer :: ia,il,l
        complex(8) :: sumiv2   ! intermediate sum (over G')
        integer, allocatable :: jipw(:,:)
        complex(8), allocatable :: inti(:,:),tmat1(:,:)
        complex(8), allocatable :: tvec1(:),tvec2(:)

        nvk=nv(irk) 
        allocate(jipw(nvk,nvk),inti(nvk,nvk),tmat1(nvk,ibgw:nbgw),       &
     &         tvec1(nvk),tvec2(nvk))
  
        do iv2=1, nvk
          if(isym.eq.1) then 
            ig2 = indgkir(iv2,irk)
          else 
            ig2 = indgk(iv2,ik)
          endif 

          do iv1=1, nvk

            if(isym.eq.1) then
              ig1 = indgkir(iv1,irk)
            else
              ig1 = indgk(iv1,ik)
            endif

            do i=1,3
              ikvec(i) = gindex(i,ig1)-gindex(i,ig2)
            enddo 

            jipw(iv1,iv2) = ig0(ikvec(1),ikvec(2),ikvec(3))
          enddo ! iv1
        enddo ! iv2   
       
        do ikxc = 1, nksxc  !! loop over the star 

          do iv2=1, nvk   !* loop over G'
            do iv1=1, nvk   !* loop over G 
              inti(iv1,iv2) = istpw(jipw(iv1,iv2),ikxc)
            enddo ! iv2
          enddo ! iv1

          call cpu_time(time1) 
          call zgemm('n','n',nvk,nbandsgw,nvk,vxcs(ikxc,isp),inti,nvk,&
     &           zzk(:,ibgw:nbgw),maxngk,czero,tmat1,nvk)
          call cpu_time(time2) 
          time_lapack = time_lapack + time2 - time1

          do ien=ibgw,nbgw 
            tvec2(1:nvk)=tmat1(1:nvk,ien)
            do iem=ibgw,nbgw
              tvec1(1:nvk) = zzk(1:nvk,iem)
              vxc_is(iem,ien)=vxc_is(iem,ien)+zdotc(nvk,tvec1,1,tvec2,1)
            enddo 
          enddo ! ie
        enddo ! ikxc
        deallocate(jipw,inti,tmat1,tvec1,tvec2) 

        end subroutine sub_calcvxcmn_is


!
!This subroutine calculates the matrix elements $v^{orb}_{nn}(\vec{k})$ 
!
        subroutine sub_calcvorbmn(iat,idf)
        implicit none 
        integer,intent(in)::iat,idf

! !LOCAL VARIABLES:
        integer(4) :: ik      ! the order index of kvec
        integer(4) :: ien,iem ! Counter: run over the eigenvalues in the corresponding k-point.
        integer(4) :: m1      ! m-quantum number of the (L)APW eigenfunction  at k
        integer(4) :: m2      ! m-quantum number of the (L)APW eigenfunction   at k' = k - q
        integer(4) :: l1m1    ! Counter: Runs over MTS (L)APW basis  functions (l1,m1)
        integer(4) :: l2m2    ! Counter: Runs over MTS (L)APW basis functions l2,m2)
        integer(4) :: ia,il,l,itmp
      
        real(8) :: udud,uu2,udu2
        complex(8) :: af,bt,gm,caf,cbt,cgm
        complex(8) :: nimm(-3:3,-3:3,ibgw:nbgw)
        complex(8) :: trvd,sumvv,vmm  
!EOP
!
!BOC
        do ia=1,natorb
          if(iatorb(ia) .ne.iat) cycle   

#ifdef DEBUG
          write(9,*) 
          write(9,*) '#vorb for atom ',iat
          write(9,*) 
#endif

          do il=1,nlorb(ia) 
            l=lorb(il,ia) 
          
            udud=uiorb(3,il,ia,isp)
            uu2= uiorb(4,il,ia,isp)
            udu2=uiorb(5,il,ia,isp)

            do m2=-l,l
              do m1=-l,l 
                l1m1= l*l+l+m1+1
                l2m2= l*l+l+m2+1
                vmm=vorb(m1,m2,il,ia,isp)

                do ien=ibgw,nbgw 
                  af=alfa(ien,l1m1,idf)
                  bt=beta(ien,l1m1,idf)
                  gm=gama(ien,1,l1m1,idf)

                  caf=conjg(alfa(ien,l2m2,idf))
                  cbt=conjg(beta(ien,l2m2,idf))
                  cgm=conjg(gama(ien,1,l2m2,idf))

                  nimm(m1,m2,ien)= af*caf  + bt*cbt*udud + gm*cgm       &
     &                         +(af*cgm+gm*caf)*uu2                     &
     &                         +(bt*cgm+gm*cbt)*udu2 
                
                  dmorb(m1,m2,il,ia) = dmorb(m1,m2,il,ia)                 &
     &               + wkir(irk)*kiw(ien,irk,isp)*nimm(m1,m2,ien)

                  do iem=ibgw,nbgw 

                    caf=conjg(alfa(iem,l2m2,idf))
                    cbt=conjg(beta(iem,l2m2,idf))
                    cgm=conjg(gama(iem,1,l2m2,idf))
                    sumvv= vmm*(af*caf + bt*cbt*udud + gm*cgm             &
     &                 +(af*cgm+gm*caf)*uu2+(bt*cgm+gm*cbt)*udu2 )
                    vxc_orb(iem,ien)=vxc_orb(iem,ien) + sumvv
                  enddo ! ien
                enddo ! iem
              enddo   ! m1
            enddo ! m2

#ifdef DEBUG
            write(9,'(a,6e16.6)') 'uiorb',uiorb(1:6,il,ia,isp) 
            do ien=ibgw,nbgw 
              write(9,*) 
              write(9,*) '--- nimm for l=',l,'and i=',ien
              write(9,*) 
              do m2=-l,l
                do m1=-l,l
                  write(9,'(2E16.8)') nimm(m1,m2,ien)
                enddo
              enddo
            enddo 
#endif
          enddo  ! il
        enddo ! ia
        end subroutine sub_calcvorbmn

      end subroutine w2k_calcvxcmn

