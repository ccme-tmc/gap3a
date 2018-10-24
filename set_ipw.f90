!BOP
!
! !ROUTINE: set_ipw
!
! !INTERFACE:
      subroutine set_ipw

! !DESCRIPTION:
!   This driver subroutine contains all initialization operations related to interstitial plane wave (IPW) 
!      1. Sets the indexes of the reciprocal lattice vectors used for the IPW expansion of 
!         the mixed basis and stores them in memory, ordered by increasing length
!
! !USES:

      use lapwlo,     only: kmax
      use mixbasis,   only: igm,ipwint,kmr,pwm
      use kpoints,    only: nqp,idvq, qlist
      use recipvec  
      use struk,      only: pia,br2,ortho,nat,nrpt
      use task,       only: casename,iop_dftpkg,fid_outdbg
      
! !INPUT PARAMETERS:

      implicit none     
      

! !LOCAL VARIABLES:
      logical:: ldbg = .false.
      character(20):: sname="set_ipw"
      integer(4) :: i         ! (Counter): runs over coordinates
      integer(4) :: ierr      ! error control
      integer(4) :: ig,igq
      integer(4) :: ippw      ! (Counter): runs over plane waves
      integer(4) :: ig1       ! (Counter): run over x-coord of G
      integer(4) :: ig2       ! (Counter): run over y-coord of G 
      integer(4) :: ig3       ! (Counter): run over z-coord of G 
      integer(4) :: igl
      integer(4) :: iq
      integer(4) :: jgl
      integer(4) :: ng1       ! Max. ig1
      integer(4) :: ng2       ! Max. ig2
      integer(4) :: ng3       ! Max. ig3
      integer(4) :: kxcmax

      integer(4), dimension(3) :: igvec,iqvec ! Integer coordinates of the G-vector

      real(8) :: kk           ! Length of qpg
      real(8) :: q
      real(8) :: gmax
      real(8) :: maxgxc
      real(8) :: maxgcoul
      real(8) :: maxlenmb
      real(8) :: maxlencoul
      real(8) :: glprev,gqlen,bleng

      real(8), dimension(3) :: gpq      ! q+G in cartesian coordinates
      real(8), dimension(3) :: gvec,qvec     ! Cartesian coordinates of the
!                                         G-vector

      integer(4), allocatable :: gind(:,:) ! Temporary storage for gindex
      real(8), allocatable ::    glen(:)   ! Temporary storage for all the
      real(8), allocatable ::    lentemp(:,:) 
      
! !EXTERNAL ROUTINES: 


      external k2cart

! !REVISION HISTORY:
!
! Created May 2004 by RGA
! Last Modified 9. Nov. 2005 by RGA
!
!EOP
!BOC

      call linmsg(6,"-",sname)

!
!     read maxngk from the energy file
!
      if(iop_dftpkg.eq.0) then
        call w2k_readenergy(-1)  
      endif 

      call init_recipvec
!
! get the kxcmax, the cutoff of IPWs used in the expansion of vxc
!
      call sub_get_kxcmax

!
! set igm(1:3) -- ! maximum number of G vectors for the theta coeficient
! expansion of the IPW
!

      bleng = sqrt(br2(1,1)*br2(1,1)+br2(2,1)*br2(2,1)+&
     &  br2(3,1)*br2(3,1))
      ng1 = idint(kmr*kmax*pwm/bleng)+1
      bleng = sqrt(br2(1,2)*br2(1,2)+br2(2,2)*br2(2,2)+&
     &  br2(3,2)*br2(3,2))
      ng2 = idint(kmr*kmax*pwm/bleng)+1
      bleng = sqrt(br2(1,3)*br2(1,3)+br2(2,3)*br2(2,3)+&
     &  br2(3,3)*br2(3,3))
      ng3 = idint(kmr*kmax*pwm/bleng)+1
     
      if(ldbg) then  
        write(6,*) "  kxcmax=",kxcmax 
        write(6,*) "  ng1,ng2,ng3=",ng1,ng2,ng3
      endif 

      igm(1)=max(4*ng1,2*kxcmax+ng1)
      igm(2)=max(4*ng2,2*kxcmax+ng2)
      igm(3)=max(4*ng3,2*kxcmax+ng3)

!----------------------------------------------------------------------!
!               Set gindex: index of G-vectors in integer coordiate    !
!----------------------------------------------------------------------!
      ng1=igm(1)/2
      ng2=igm(2)/2
      ng3=igm(3)/2
      npw = (2*ng1+1)*(2*ng2+1)*(2*ng3+1)
      allocate(glen(npw),gind(3,npw))
      maxgxc=2.0d0*kmax+4.0d0
      maxgcoul=kmr*kmax*(pwm+2.0d0)
      gmax=max(maxgxc,maxgcoul)
      ippw=0
      do ig1=-ng1,ng1
        do ig2=-ng2,ng2
          do ig3=-ng3,ng3

            igvec(1)=ig1
            igvec(2)=ig2
            igvec(3)=ig3
            do i=1,3
              gvec(i)=(dble(igvec(1))*br2(i,1)+dble(igvec(2))*          &
     &                 br2(i,2)+dble(igvec(3))*br2(i,3))
            enddo 

            kk=sqrt(sum(gvec(1:3)**2))
            if(kk.lt.2.0d0*kmax+4.0d0) npw2apw=npw2apw+1

            if( kk.lt.gmax )then
              ippw = ippw + 1
              glen(ippw)=kk
              if(ortho)then
                do i=1,3
                  q=gvec(i)/pia(i)
                  gind(i,ippw)=nint(q+sign(1.0d-2,q))
                enddo
              else  
                gind(1:3,ippw) = igvec(1:3)
              endif  
            endif
          enddo
        enddo
      enddo
      npw = ippw

      ! sort by increasing length using shell algorithm
      call shelsort(npw,gind,glen(1:npw))

      allocate(gindex(1:3,1:npw),stat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to allocate gindex,gindl")

      do ippw=1,npw
        gindex(:,ippw)=gind(:,ippw)
      enddo

     
!----------------------------------------------------------------------!
!     generate the inverse of gindex                                   !
!----------------------------------------------------------------------!
      ngmax=maxval(gindex)
      ngmin=minval(gindex)
      allocate(ig0(ngmin:ngmax,ngmin:ngmax,ngmin:ngmax))
      ig0(ngmin:ngmax,ngmin:ngmax,ngmin:ngmax)=0
      do ippw=1,npw
        ig1=gindex(1,ippw)
        ig2=gindex(2,ippw)
        ig3=gindex(3,ippw)
        ig0(ig1,ig2,ig3)=ippw
      enddo  

!---------------------------------------------------------------------!
!                 calculate the integral of IPW                       !
!---------------------------------------------------------------------!

      allocate(ipwint(1:npw))
      do ippw=1,npw
        igvec(1:3)=gindex(:,ippw)
        call int1ipw(igvec,ipwint(ippw))
      enddo

!---------------------------------------------------------------------!
!               set ngq(nqp) and ngqbarc*nqp)                         !
!-------------------------------------------------------------------- !

      allocate(ngq(nqp),ngqbarc(nqp))
      maxlenmb = kmr*kmax
      maxlencoul = pwm * maxlenmb
      do iq=1, nqp
        ngq(iq)=0
        ngqbarc(iq)=0
        iqvec(1:3) = qlist(1:3,iq)
        call k2cart(iqvec,idvq,qvec)
        do ig = 1, npw
          igvec(1:3)= gindex(:,ig)
          call k2cart(igvec,1,gvec)
          gpq = gvec + qvec
          gqlen=sqrt(gpq(1)*gpq(1)+gpq(2)*gpq(2)+gpq(3)*gpq(3))
          if(gqlen.lt.maxlenmb) ngq(iq)=ngq(iq)+1
          if(gqlen.lt.maxlencoul) ngqbarc(iq)=ngqbarc(iq)+1
        enddo ! ig
      enddo ! iq

      maxngq=maxval(ngq)
      maxngqbarc=maxval(ngqbarc)

      if(maxngq.eq.npw) write(6,*)'WARNING !! maxngq = npw !!!'
      if(maxngqbarc.eq.npw) write(*,*)'WARNING!! maxngqbarc = npw !!!'


!----------------------------------------------------------------------!
!             set indgq                                                !
!----------------------------------------------------------------------!

      allocate(indgq(maxngqbarc,nqp),indgqlen(maxngqbarc,nqp),ngqlen(nqp))

      allocate(lentemp(maxngqbarc,nqp))
      indgq=0
      lentemp=0.0d0

      do iq=1, nqp
        iqvec(1:3) = qlist(1:3,iq)
        call k2cart(iqvec,idvq,qvec)

        igq=0
        do ig = 1, npw
          igvec(1:3)= gindex(:,ig)
          call k2cart(igvec,1,gvec)
          gpq = gvec + qvec
          gqlen=sqrt(gpq(1)*gpq(1)+gpq(2)*gpq(2)+gpq(3)*gpq(3))

          if(gqlen.lt.maxlencoul)then
            igq=igq+1
            gind(:,igq)=igvec(1:3)
            glen(igq)=gqlen
          endif
        enddo ! ig

        call errmsg(igq.ne.ngqbarc(iq),sname,"igq != ngqbarc")
        call shelsort(ngqbarc(iq),gind(:,1:ngqbarc(iq)),              &
     &                glen(1:ngqbarc(iq)))
        glprev=-1.0d0
        igl=0
        do ig=1,ngqbarc(iq)
          indgq(ig,iq)=ig0(gind(1,ig),gind(2,ig),gind(3,ig))
          if(abs(glen(ig)-glprev).gt.1.0d-6)then
            igl=igl+1
            glprev=glen(ig)
            lentemp(igl,iq)=glen(ig)
          endif
          indgqlen(ig,iq)=igl
        enddo
        ngqlen(iq)=igl

      enddo ! iq
      maxngqlen=maxval(ngqlen)
      allocate(gqleng(maxngqlen,nqp))
      gqleng(1:maxngqlen,1:nqp)=lentemp(1:maxngqlen,1:nqp)

      deallocate(lentemp,gind,glen)

      if(ldbg) then 
        write(6,*) " write detailed info into outdbg"
        write(fid_outdbg,*) 
        write(fid_outdbg,*) "### gindex ###"
        do ippw=1,npw,npw/20
          write(fid_outdbg,'(4i5)') ippw,gindex(:,ippw)
        enddo  

        write(fid_outdbg,*) 
        write(fid_outdbg,*) "### indgw ###"
        do iq=1,nqp
          write(fid_outdbg,*) 
          write(fid_outdbg,*) "iq= ",iq
          do ig=1,ngqbarc(iq),ngqbarc(iq)/20
            write(fid_outdbg,'(3i5)') ig,indgq(ig,iq),indgqlen(ig,iq)
          enddo
        enddo 
      endif 

      write(6,200) "Nr. of IPW(npw):",npw
      write(6,200) "Nr. of pws for lapw overlaps (npw2apw):",npw2apw
      write(6,200) "Max nr. of IPW for APW (maxngk):",maxngk
      write(6,200) "Max nr. of IPW for Mixbasis (maxngq) :",maxngq
      write(6,200) "Max nr. of IPW for bare Coulomb (maxngqbarc) :",maxngqbarc

      return


   10 format(i6,4f12.6)
   12 format('kmax = ',f15.10,'kmr = ',f15.10,'pwm = ',f15.10)
   13 format('maxgxc = ',f15.10,'maxgcoul = ',f15.10,'gmax = ',f15.10)
  101 format('|G_max| = ',f15.10,' Nr. of planewaves =',i6)
  200 format(a50,i5) 

      contains 

        subroutine sub_get_kxcmax 
        implicit none 
        integer:: fid
        integer:: j,iat, npt, nlm,lxc,lxcmax,lm(2),irp,nksxc,ikxc
        integer,allocatable::ksxc(:,:)
        real(8):: tmp

        fid=999
        open(unit=fid,file=trim(casename)//".vxc",action='read',iostat=ierr)
        call errmsg(ierr.ne.0,sname,"Fail to open case.vxc")

        read(fid,10)
        do iat=1, nat
          npt=nrpt(iat)
          read(fid,11) lxcmax
          do lxc=1, lxcmax
            read(fid,12) lm(1:2)
            read(fid,13) (tmp,irp=1,npt)
            read(fid,14)
          enddo ! lxc
          read(fid,15)
        enddo
        read(fid,16) nksxc
        allocate(ksxc(1:3,1:nksxc))
        do ikxc=1, nksxc
          read(fid,17)(ksxc(j,ikxc),j=1,3),tmp
        enddo ! ikxc
        close(fid)
        kxcmax=maxval(abs(ksxc(1:3,1:nksxc)))
   10   format(//)
   11   format(/,15x,i3,//)
   12   format(15x,i3,5x,i2,/)
   13   format(3x,4e19.12)
   14   format(/)
   15   format(///)
   16   format(//,13x,i6)
   17   format(3x,3i5,2e19.12)
        end subroutine 
 
      end subroutine set_ipw
!EOC
