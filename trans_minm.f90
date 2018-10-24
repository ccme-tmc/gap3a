      subroutine trans_minm(iop,nmflag,minm,nst,nend,mst,mend,ik,iq,isp)
!
!     Transform Minm by qpwf_coef 
!     Main inputs:
!       iop -- control how to transform  Minm(i,m,n) 
!         0 -> do nothing 
!         1 -> transform only "m" 
!         2 -> transform both "n" and "m" 
!       nmflag -- "nm","nc" or "cm" 
!
      use bands,      only: ibgw, nbgw
      use constants,  only: cone, czero
      use kpoints,    only: kqid,kpirind,idikp
      use mixbasis,   only: matsiz 
      use selfenergy, only: qpwf_coef
      use task,       only: time_lapack
      implicit none 
      integer, intent(in):: iop
      integer, intent(in):: ik,iq,isp
      character(2), intent(in):: nmflag
      integer, intent(in):: nst,nend,mst,mend
      complex(8), intent(inout):: minm(matsiz,mst:mend,nst:nend)

      character(20):: sname = "trans_minm"
      integer:: ierr 
      integer:: irk,jrk,jk 
      integer:: mdim,ndim,nmdim
      integer:: m1,n1,m2,n2,inm1,inm2
      real(8):: tstart,tend
      complex(8),allocatable::mmat(:,:,:),cmat(:,:)

      call cpu_time(tstart)
      if(iop.le.0) return 

      irk=kpirind(ik)
      jk=kqid(ik,iq)
      jrk=kpirind(jk)
!
!     Transform both "n" and "m"
!
      ndim = nend - nst + 1
      mdim = mend - mst + 1
      nmdim = mdim*ndim 

      if(iop.eq.2.and.nmflag.eq.'nm') then 
        allocate(cmat(nmdim,nmdim),mmat(matsiz,mst:mend,nst:nend),&
     &           stat=ierr)
        call errmsg(ierr.ne.0,sname,"fail to allocate cmat,mmat")

        inm2=0
        do n2 = nst, nend
          do m2 = mst,mend
            inm2 = inm2 + 1
            inm1=0
            do n1 = nst, nend
              do m1 = mst, mend
                inm1 = inm1 + 1
                cmat(inm1,inm2)=f_trans(n1,n2,irk) & 
     &                        *conjg(f_trans(m1,m2,jrk))
              enddo
            enddo
          enddo
        enddo
        mmat = minm
        call zgemm('n','n',matsiz,nmdim,nmdim,cone,mmat,matsiz, &
     &             cmat,nmdim,czero,minm,matsiz)
        deallocate(mmat,cmat)
!
!     Transform only "m"
!
      elseif(    (iop.eq.1.and.nmflag.ne.'nc') &
     &       .or.(iop.eq.2.and.nmflag.eq.'cm')) then
        allocate(cmat(mst:mend,mst:mend),&
     &           mmat(matsiz,mst:mend,1),stat=ierr)
        call errmsg(ierr.ne.0,sname,"fail to allocate cmat,mmat")
        do m2=mst,mend
          do m1=mst,mend 
            cmat(m1,m2) = conjg(f_trans(m1,m2,jrk))
          enddo
        enddo 

        do n1 = nst, nend  
          mmat(:,:,1) = minm(:,mst:mend,n1)  
          call zgemm('n','n',matsiz,mdim,mdim,cone,mmat,matsiz, & 
     &             cmat,mdim,czero,minm(:,:,n1),matsiz)
        enddo 

        deallocate(mmat,cmat)
!
!     Transform only "n"
!
      elseif(iop.eq.2.and.nmflag.eq.'nc') then
        allocate(cmat(nst:nend,nst:nend),&
     &           mmat(matsiz,mst:mend,nst:nend),stat=ierr)
        call errmsg(ierr.ne.0,sname,"fail to allocate cmat,mmat")
        do n2=nst,nend 
          do n1=nst,nend 
            cmat(n1,n2) = f_trans(n1,n2,irk)
          enddo
        enddo  
        mmat=minm 
        call zgemm('n','n',matsiz*mdim,ndim,ndim,cone,mmat,matsiz*mdim,& 
     &              cmat,ndim,czero,minm,matsiz*mdim)
        deallocate(mmat,cmat)
      endif
      call cpu_time(tend)
      time_lapack = time_lapack + tend - tstart

      contains 
         complex(8) function f_trans(n,m,ikir)
         integer:: n, m, ikir

         if(n.ge.ibgw.and.n.le.nbgw.and.m.ge.ibgw.and.m.le.nbgw) then 
           f_trans = qpwf_coef(n,m,ikir,isp) 
         else 
           if(n.eq.m) then 
             f_trans = cone 
           else
             f_trans = czero 
           endif
         endif  
         end function 

      end subroutine 

