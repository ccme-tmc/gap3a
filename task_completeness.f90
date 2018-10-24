!BOP
!
! !ROUTINE: task_completeness
!
! !INTERFACE:
      subroutine task_completeness
      
! !DESCRIPTION:
!
! This subroutine tests the completeness of the mixed basis
!      
! !USES:

      use bands,        only: nomax, numin
      use kpoints,      only: klist, qlist, nkp, nqp
      use lapwlo,       only: nt
      use liboct_parser
      use mixbasis,     only: init_mixbasis, end_mixbasis, mbsiz
      use minmmat,      only: init_minm, end_minm, minm
      use rspevec,      only: ibmin, ibmax, ibmin2, ibmax2, iik, jjk
      use task,         only: taskname
!      
! !LOCAL VARIABLES:
      
      implicit none
      
      integer(4) :: ib1, ib2, iq, iq1
      integer(4) :: ik1(3),ik2(3),iqvec(3),ngr
      integer(4) :: ierr
      
      real(8) ::  tt(2)
      real(8) :: relermt, releri
      real(8) :: maxrelermt, maxreleri
      real(8) :: minrelermt, minreleri
            
      character(len=17) :: fnat
      character(len=16) :: fnint
      character(len=18) :: fntot
      
      logical :: pflag
      
! !INTRINSIC ROUTINES:
      
      intrinsic cpu_time     

! !REVISION HISTORY:
!     
! Created 16. Sept 2005 by RGA
!
!EOP
!BOC
      call linmsg(6,'-',"task:comp")
      ierr=loct_parse_isdef(taskname)
      if(ierr.ne.0) then
        call loct_parse_block_int(taskname,0,0,iik)
        call loct_parse_block_int(taskname,0,1,jjk)
        call loct_parse_block_int(taskname,0,2,ibmin)
        call loct_parse_block_int(taskname,0,3,ibmax)
        call loct_parse_block_int(taskname,0,4,ibmin2)
        call loct_parse_block_int(taskname,0,5,ibmax2)
      else
        iik=1
        jjk=1
        ibmin=nomax
        ibmax=numin
        ibmin2=nomax
        ibmax2=numin
      endif
      fnat='intevppat-'//achar(iik+48)//'-'//achar(jjk+48)//'.out'
      fnint='intevppi-'//achar(iik+48)//'-'//achar(jjk+48)//'.out'
      fntot='intevpptot-'//achar(iik+48)//'-'//achar(jjk+48)//'.out'
      pflag=.false.
      open(unit=71,file=fnat,status='old',err=901)
      open(unit=72,file=fnint,status='old',err=902)
      open(unit=73,file=fntot,status='old',err=903)
      goto 1
  901   pflag=.true.
        open(unit=71,file=fnat,status='unknown')
  902   if(.not.pflag)pflag=.true.
        open(unit=72,file=fnint,status='unknown')
  903   if(.not.pflag)pflag=.true.
        open(unit=73,file=fntot,status='unknown')
    1 continue
      fnat='intevpmat-'//achar(iik+48)//'-'//achar(jjk+48)//'.out'
      fnint='intevpmi-'//achar(iik+48)//'-'//achar(jjk+48)//'.out'
      fntot='intevpmtot-'//achar(iik+48)//'-'//achar(jjk+48)//'.out'
      open(unit=74,file=fnat,status='unknown')
      open(unit=75,file=fnint,status='unknown')
      open(unit=76,file=fntot,status='unknown')
      ik1(1:3)=klist(1:3,iik)
      ik2(1:3)=klist(1:3,jjk)
      iqvec(1:3)=ik2(1:3)-ik1(1:3)
      iq=1
      do iq1=1,nqp
        if((iqvec(1).eq.qlist(1,iq1)).and.(iqvec(2).eq.qlist(2,iq1))    &
     &     .and.(iqvec(3).eq.qlist(3,iq1)))iq=iq1
      enddo
        call init_mixbasis(iq)
        call init_minm
      call cpu_time(tt(1))
      call coul_mpwipw(iq)
      call cpu_time(tt(2))
      write(98,102)tt(2)-tt(1)
      call flushbuf(98)
      allocate(minm(1:mbsiz,ibmin:ibmax,ibmin2:ibmax2))
      call cpu_time(tt(1))
      call calcminm(iik,iq,ibmin,ibmax,ibmin2,ibmax2,1,.false.)
      call cpu_time(tt(2))
      write(98,103)tt(2)-tt(1)
      call flushbuf(98)
      ngr=16*nt*nt/3
      call prep_ang_int(nt,ngr)
     
      if(pflag)then
        call cpu_time(tt(1))
        do ib1=ibmin,ibmax
          do ib2=ibmin2,ibmax2
            call intevecpp(iik,jjk,ib1,ib2)
          enddo  
        enddo  
        call cpu_time(tt(2))
        write(98,104)tt(2)-tt(1)
        call flushbuf(98)
      endif  
      rewind(71)
      rewind(72)
      rewind(73)
      call cpu_time(tt(1))
      maxreleri=0.0d0
      maxrelermt=0.0d0
      minreleri=1.0d+3
      minrelermt=1.0d+3
      do ib1=ibmin,ibmax
        do ib2=ibmin2,ibmax2
          call intevecpm(jjk,ib1,ib2,relermt,releri)
          if(abs(releri).gt.maxreleri)maxreleri=abs(releri)
          if(abs(relermt).gt.maxrelermt)maxrelermt=abs(relermt)
          if(abs(releri).lt.minreleri)minreleri=abs(releri)
          if(abs(relermt).lt.minrelermt)minrelermt=abs(relermt)
        enddo  
      enddo 
      write(6,*)'mt. min/max error',minrelermt,maxrelermt 
      write(6,*)'i. min/max error',minreleri,maxreleri 
      call cpu_time(tt(2))
      write(98,105)tt(2)-tt(1)
      call flushbuf(98)
      deallocate(minm)
      call end_mixbasis
      call end_minm
      close(71)
      close(72)
      close(73)
      close(74)
      close(75)
      close(76)
  101 format(' diagsgi    -->',f14.5,' seconds')
  102 format(' coul_mpwipw -->',f14.5,' seconds')
  103 format(' calcminm     -->',f14.5,' seconds')
  104 format(' integral prods -->',f14.5,' seconds')
  105 format(' integral mix   -->',f14.5,' seconds')
      return
      
      end subroutine task_completeness            
!EOC
