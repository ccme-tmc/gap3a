!BOP
!
! !ROUTINE: w2k_readmommat
!
! !INTERFACE:
      subroutine w2k_readmommat(nst,nend,mst,mend)
      
! !DESCRIPTION:
!
! This subroutine calculates the momentum matrix elements of the (L)APW+lo
! eigenfunctions      
!
! !USES:
      
      use bands,       only: nspin
      use mommat,      only: mmatvv
      use kpoints,     only: nirkp
      use task,        only: fid_outmom,casename
      
! !INPUT PARAMETERS:

!EOP  
!BOC  
      implicit none
      integer, intent(in) :: nst,nend,mst,mend   !! the range of 1st band index 
         
! !LOCAL VARIABLES:

      integer :: ierr 
      integer :: ik,irk              ! index of the k point.
      integer :: ie1,ie2,nb1,nb2     ! Counter: run over bands.
      integer :: isp         ! index for spin 
      integer :: nmin,nmax   ! k-dependent number of bands
      
      integer :: fout,fin  

      complex(8) :: p12(3),p21(3)
      real(8) :: renorm
      logical :: ldbg=.true.
      character(20):: sname="w2k_readmommat"
      character(80):: fin_name  

! !EXTERNAL ROUTINES: 
      external expand_evec
      external readvector
      external calcmmatvv
      external calcmmatcv
 
      if(ldbg) call linmsg(6,'-',trim(sname))  

      fin = 999
      fout = fid_outmom

      mmatvv = 0.d0

      do isp=1,nspin 

        !! set mommat file name 
        if(nspin.eq.1) then 
          fin_name = trim(casename)//".mommat"
        else 
          if(isp.eq.1) then 
            fin_name = trim(casename)//".mommatup"
          else
            fin_name = trim(casename)//".mommatdn"
          endif 
        endif 
        open(unit=fin,file=trim(fin_name),action='read',iostat=ierr) 
        call errmsg(ierr.ne.0,sname,"Fail to open "//trim(fin_name))
        read(fin,*) 

        do irk=1,nirkp
           
          read(fin,*) 
          read(fin,101) nmin,nmax 
          read(fin,*) 
          
          do nb1=nmin,nmax
            do nb2=nb1,nmax
              read(fin,102) ie1,ie2,p12(1:3)
              if ( (ie1.ge.nst.and.ie1.le.nend).and. &
     &             (ie2.ge.mst.and.ie2.le.mend) ) then 
                 mmatvv(:,ie1,ie2,irk,isp) = p12
              elseif( (ie2.ge.nst.and.ie2.le.nend).and. &
     &                (ie1.ge.mst.and.ie1.le.mend) ) then 
                 mmatvv(:,ie2,ie1,irk,isp) = conjg(p12)
              endif 
            enddo
          enddo 

          !! write out mommat 
          if(ldbg) then  
            write(fout,*) "#momentum matrix for ik=",ik,"irk=",irk
            write(fout,'(a1,a5,a6,3a12)')'#',"ie1","ie2","px","py","pz"
            do ie1=nst,nend
              do ie2=mst,mend
                p12 = mmatvv(:,ie1,ie2,irk,isp)
                write(fout,'(2i6,4f12.6)')ie1,ie2,real(p12), &
     &            real(sum(p12*conjg(p12))/3.d0)
              enddo ! ie2
            enddo 
          endif 

        enddo ! irk
        close(fin)
      enddo ! isp
 101  format(27x,2i5)
 102  format(3X,2I4,6E13.6)
      end subroutine w2k_readmommat
!EOC
            
