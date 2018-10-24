!BOP
!
! !ROUTINE: w2k_calcvxcnn
!
! !INTERFACE:
      subroutine w2k_calcvxcnn()

! !DESCRIPTION:
!
!This subroutine calculates the matrix elements $v^{xc}_{nn}(\vec{k})$ 
!of equation \ref{vxc-10}, only for valence states.
!
!
! !USES:

      use xcpot,       only: uxcu,lvorb,lorb,dmorb,lxcm,lmxc,istpw, &
     &                  vorbnn,vxc_hyb,natorb,iatorb,nlorb,vorb,ksxc,&
     &                  vxclm,lhybrid,vxcs,nksxc
      use bands,       only: nv, ibgw,nbgw,nspin
      use constants,   only: cone, czero, imag
      use eigenvec,    only: alfa, beta, gama, zzk
      use kpoints,     only: idvkir, kirlist, nirkp, idikp
      use lapwlo,      only: lmax,lomax,nlo_at
      use recipvec,    only: npw,gindex, ig0, indgkir,ngmax,ngmin, maxngk
      use selfenergy,  only: vxcnn
      use struk,       only: iatnr, mult,nat
      use task,        only: casename,time_lapack,fid_outdbg,time_vxc, &
     &                       lrestart
      implicit none

! !LOCAL VARIABLES:
      integer :: irk  ! the order index of the irred. k-point
      integer :: ikp  ! the order index of kvec

      integer :: i
      integer :: isp  ! index for spin 

      integer :: bl      ! Angular quantum number of the   mixed atomic functions (big l)
      integer :: bm      ! Magnetic quantum number of the mixed atomic function (big m)
      integer :: iat     ! Counter: runs over inequivalent atoms.
      integer :: idf     ! Counter: runs over atoms, including  equivalent ones
      integer :: ie     ! Counter: run over the eigenvalues in the corresponding k-point.
      integer :: ieq     ! Counter: runs over the equivalent atoms of each type.
      integer :: iv1,iv2 ! Counter: run over the (L)APW basis functions (excluding core states).
      integer, allocatable :: jipw(:,:)
      integer :: l1      ! l-quantum number of the (L)APW eigenfunction  at k
      integer :: l2      ! l-quantum number of the (L)APW eigenfunction  at k' = k - q
      integer :: l12
      integer :: ir1,ir2,ir12,nr1,nr2
      integer :: lxc     ! (Counter) l index of the xc component
      integer :: m1      ! m-quantum number of the (L)APW eigenfunction  at k
      integer :: m2      ! m-quantum number of the (L)APW eigenfunction   at k' = k - q
      integer :: nvk 
      integer :: l1m1    ! Counter: Runs over MTS (L)APW basis  functions (l1,m1)
      integer :: l2m2    ! Counter: Runs over MTS (L)APW basis functions l2,m2)
      integer :: ikxc    ! Counter: Runs over IPW's
      integer :: isgn 
      integer, dimension(3) :: ik ! Indexes of G_1+G'-G
      integer :: ia,il,l,itmp
      integer :: ierr
      character(len=20) :: sname="w2k_calcvxcnn"
      character(len=200):: msg
      
      real(8) :: kvec(3)
      real(8) :: time1,time2,tstart,tend

      complex(8) :: angint ! Angular integral calculated by fourylm
      complex(8), allocatable :: sumiv2(:)   ! intermediate sum (over G')
      complex(8) :: io1
      complex(8) :: af,caf
      complex(8), allocatable :: inti(:,:),tmat1(:,:)
      complex(8), allocatable :: tvec1(:),tvec2(:)
      complex(8), allocatable :: vxcmt(:,:),vxci(:,:)
      complex(8) :: trvd  
 
!
! !EXTERNAL ROUTINES: 

      real(8), external :: getcgcoef
      complex(8), external :: zdotc
      
      external expand_evec
      external intstipw
      external readvector
      external zgemm
     
!EOP
!
!BOC


!
!  Calculate the integral between stars and planwaves
!
      call cpu_time(tstart)
      if(lrestart) then 
        call io_vxcmn('r','d',1,ierr)
        if(ierr.eq.0) then
          write(6,*) "Restart mode: read eps from files"
          return
        endif
      endif
  
      call intstipw

      allocate(vxcmt(ibgw:nbgw,nirkp))
      allocate(vxci(ibgw:nbgw,nirkp))
      allocate(sumiv2(ibgw:nbgw))

      do isp=1,nspin 

        vxcmt(:,:)=czero
        vxci(:,:)=czero
        if(lvorb) dmorb=0.d0
        do irk=1, nirkp                      !* Loop over the k-points:
          ikp=idikp(irk) 

          call readvector(ikp,1,isp,0)          !* Read the eigenvector corresponding to the k-point ikp
          call expand_evec(ikp,1,.true.,isp)  !* Calculate the expansion coeficients of the eigenvectors

!-------------------------------------------------------------
!      Muffin-Tin Spheres        
!-------------------------------------------------------------
          idf=0
          do iat = 1, nat                       !*  Loop over inequivalent atoms:
            do ieq = 1, mult(iat)               !* Loop over equivalent atoms:
              idf = idf + 1

              do l2=0,lmax                     !*    loop over l2:
              do l1=0,lmax                      !* Loop over l1

                !! nr1,nr2 : the number of radial functions 
                if(l1.gt.lomax) then 
                  nr1 = 2
                else
                  nr1 = nLO_at(l1,iat) + 2
                endif 
                
                if(l2.gt.lomax) then 
                  nr2 = 2
                else 
                  nr2 = nLO_at(l2,iat) + 2 
                endif                 
                 
                l12 = l1+1 + l2*(lmax+1)

                do m1=-l1,l1                  !* loop over m1
                  l1m1 = l1*l1+l1+m1+1

                  do lxc = 1, lxcm(iat)             !*  Loop over vxc components:
                    bl = abs(lmxc(1,lxc,iat))
                    bm = lmxc(2,lxc,iat)

                    io1 = cone

                    if(lmxc(1,lxc,iat).lt.0) io1=imag

                    if( .not.( (l1.ge.abs(bl-l2)).and.(l1.le.bl+l2))) cycle 


                    do isgn = 1, -1, -2 
                      
                      bm = bm*isgn

                      m2 = m1 - bm

                      if(abs(m2).gt.l2) cycle 
                      
                      l2m2 = l2*l2 + l2 + m2 + 1

                      !* Calculate the angular integral:
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
     &                     *(-1.0)**bm*sign(1,lmxc(1,lxc,iat))/sqrt(2.0d0)
                        endif
                      endif 
      
                      if(abs(angint).lt.1.0d-12) cycle 

                      do ie=ibgw,nbgw

                        ir12 = 0 

                        do ir2 = 1,nr2  
                          if(ir2.eq.1) then 
                            af = alfa(ie,l2m2,idf)
                          elseif(ir2.eq.2) then 
                            af = beta(ie,l2m2,idf)
                          else
                            af = gama(ie,ir2-2,l2m2,idf)
                          endif

                          do ir1 = 1, nr1
                            ir12 = (ir2-1)*nr1 + ir1

                            if(ir1.eq.1) then 
                              caf = conjg(alfa(ie,l1m1,idf))
                            elseif(ir1.eq.2) then
                              caf = conjg(beta(ie,l1m1,idf))
                            else
                              caf = conjg(gama(ie,ir1-2,l1m1,idf))
                            endif 

                            vxcmt(ie,irk) = vxcmt(ie,irk) &
     &                       + caf*af*uxcu(ir12,l12,lxc,iat,isp)*angint
                          enddo
                        enddo
                      enddo ! ie 
                     
                      if(bm.eq.0) exit  

                    enddo ! isign

                  enddo ! lxc
                enddo ! m1     
              enddo ! l1
              enddo ! i2

!
!  Calculate contribution from LDA+U correction if necessary
!
              if(lvorb)  then 
                call calcvorbnn(iat,idf,irk,isp)
              endif 
 
            enddo ! ieq
          enddo ! iat

!     --------------------
!     Interstitial region:
!     --------------------

          nvk=nv(irk) 
          allocate(jipw(nvk,nvk),inti(nvk,nvk),tmat1(nvk,ibgw:nbgw),       &
     &           tvec1(nvk),tvec2(nvk))
          do iv1=1, nvk
            if(indgkir(iv1,irk).gt.npw) then
              write(msg,*) "indgkir(iv1,irk) > npw", &
     &            "iv1,irk,indgkir=",iv1,irk,indgkir(iv1,irk)
              call outerr(sname,msg)
            endif 
               
            do iv2=1, nvk
              if(indgkir(iv2,irk).gt.npw) then
                write(msg,*)  "indgkir(iv2,irk) > npw",&
     &             "  iv2,irk,indgkir=",iv2,irk,indgkir(iv2,irk)
                call outerr(sname,msg)
              endif

              do i=1,3
                ik(i)=gindex(i,indgkir(iv1,irk))-gindex(i,indgkir(iv2,irk))
              enddo 
              if((minval(ik).ge.ngmin).and.(maxval(ik).le.ngmax)) then
                jipw(iv1,iv2)=ig0(ik(1),ik(2),ik(3))
              else
                jipw(iv1,iv2)=0
                write(6,'(a,i6,2x,2i4,2x,3i4)')                         &
     &       'WARNING: jipw=0 in w2k_calcvxcnn',jipw(iv1,iv2),ngmin,ngmax,ik
              endif    
            enddo ! iv2
          enddo ! iv1    
!
!     Loop over the mixed basis functions:
!
          do ikxc = 1, nksxc
            do iv1=1, nvk   !* loop over G 
              do iv2=1, nvk   !* loop over G'
                if(jipw(iv1,iv2).gt.0)then
                  inti(iv1,iv2)=istpw(jipw(iv1,iv2),ikxc)
                else  
                  inti(iv1,iv2)=czero
                endif   
              enddo ! iv2
            enddo ! iv1

            call cpu_time(time1) 
            call zgemm('n','n',nvk,nbgw-ibgw+1,nvk,vxcs(ikxc,isp),inti,nvk, &
     &             zzk(:,ibgw:nbgw),maxngk,czero,tmat1,nvk)
            call cpu_time(time2) 
            time_lapack=time_lapack+time2-time1

            do ie=ibgw,nbgw
              tvec1(1:nvk)=zzk(1:nvk,ie)
              tvec2(1:nvk)=tmat1(1:nvk,ie)
              sumiv2(ie)=zdotc(nvk,tvec1,1,tvec2,1)
              vxci(ie,irk)=vxci(ie,irk)+sumiv2(ie)
            enddo ! ie
          enddo ! ikxc
          deallocate(jipw,inti,tmat1,tvec1,tvec2) 

        enddo ! irk   

        do irk=1,nirkp
          do ie=ibgw,nbgw
            vxcnn(ie,irk,isp)=real(vxcmt(ie,irk)+vxci(ie,irk))
          enddo
        enddo
 
        if(lhybrid) then 
          do ie=ibgw,nbgw 
            do irk=1,nirkp 
              vxcnn(ie,irk,isp)=vxcnn(ie,irk,isp) + vxc_hyb(ie,irk,isp)
            enddo
          enddo
        endif 

        if(lvorb) then 
          do irk=1,nirkp
            do ie=ibgw,nbgw
              vxcnn(ie,irk,isp)=vxcnn(ie,irk,isp)+vorbnn(ie,irk,isp)
            enddo
          enddo

         ! Some tests for vorb 
          write(fid_outdbg,*) "--- test vorb ---"

          do ia=1,natorb
            iat=iatorb(ia)
            do il=1,nlorb(ia) 
              l=lorb(il,ia) 

              write(fid_outdbg,*) 
              write(fid_outdbg,*) '--- dm for l=',l,'iat=',iat
              write(fid_outdbg,*) 
              trvd=czero
              do m2=-l,l
                do m1=-l,l
                  trvd = trvd+vorb(m2,m1,il,ia,isp)*dmorb(m1,m2,il,ia) 
                  write(fid_outdbg,'(2e16.8)') dmorb(m1,m2,il,ia) 
                enddo 
              enddo  
              write(fid_outdbg,*) 
              write(fid_outdbg,'(a,2e16.8)') "Tr(rho.V)=",trvd
            enddo
          enddo 
          write(fid_outdbg,*)
        endif  ! lvorb

      enddo ! isp   

      call io_vxcmn('w','d',1,ierr) !! write vxcnn to the hard disk 
!
! deallocate all allocated arrays that are not used any more once vxcnn is calculated 
!
      deallocate(vxclm,vxcs,lmxc,lxcm,ksxc,uxcu) 
      deallocate(vxcmt,vxci,sumiv2)
      call cpu_time(tend)
      time_vxc = time_vxc + tend-tstart
          
    1 format(3e16.6,3i6)
    2 format(i6,2e24.12)     

      end subroutine w2k_calcvxcnn

