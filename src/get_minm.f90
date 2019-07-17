      subroutine get_minm(iop,nmflag,minm,nst,nend,mst,mend,ik,iq,isp)
!
!     Read/Write Minm 
!      Main input:
!        iop - control how to get Minm 
!            -1 -> calculate on the fly
!            -2 -> calculate on the fly and transform "m" by qpwf_coef afterwards
!             0 -> read directly from the file 
!             1 -> read from the file and transform "m" by qpwf_coef
!             2 -> read from the file and transform both "n" and "m"
!
      use eigenvec, only: ik_cur_evec, jk_cur_evec
      use kpoints,  only: kqid
      use minmmat,  only: io_minm,check_minm  
      use mixbasis, only: matsiz
      implicit none

      integer, intent(in) :: iop                     !! control how to get Minm 
      character(2), intent(in):: nmflag
      complex(8), intent(inout):: minm(matsiz,mst:mend,nst:nend)
      integer, intent(in) :: nst,nend,mst,mend
      integer, intent(in):: ik,iq,isp

      character(len=10):: sname="get_minm"
      integer:: jk  ! index for k  in BZ and IBZ, respectively 
      integer:: m_stat
     
      !! check whether Minm is availale in the Minm file 
      m_stat=0
      if(iop.ge.0) then 
        call check_minm(nmflag,nst,nend,mst,mend,ik,isp,m_stat)
      endif 

      if(iop.lt.0.or.m_stat.ne.1) then !! calculate on the fly 
        jk=kqid(ik,iq)
        if(ik_cur_evec.ne.ik) then 
          call readvector(ik,1,isp,iop)      !! Read the eigenvector corresponding to the k-point ik
          call expand_evec(ik,1,.true.,isp)  !! Calculate the expansion coeficients of the eigenvectors
          ik_cur_evec = ik
        endif 
        if(jk_cur_evec.ne.jk) then
          call readvector(jk,2,isp,iop)      !! Read the eigenvector corresponding to the k'-point jk
          call expand_evec(jk,2,.true.,isp)  !! Calculate expansion coeficients of the eigenvectors
          jk_cur_evec = jk 
        endif

        if(nmflag.eq.'nm') then 
          call calcminm(ik,iq,nst,nend,mst,mend,isp,minm)
        elseif(nmflag.eq.'nc') then
          call calcminc(ik,iq,nst,nend,mst,mend,isp,minm)
        elseif(nmflag.eq.'cm') then
          call calcmicm(ik,iq,nst,nend,mst,mend,isp,minm)
        endif 
        if(iop.eq.-2) call trans_minm(1,nmflag,minm,nst,nend,mst,mend,  &
     &   ik,iq,isp)
      else !! read from the file and make some transform if required
        call io_minm('r',nmflag,minm,nst,nend,mst,mend,ik,isp)
      endif 
      if(iop.gt.0) call trans_minm(iop,nmflag,minm,nst,nend,mst,mend,ik,iq,isp)
      end subroutine
