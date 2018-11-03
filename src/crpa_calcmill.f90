      subroutine crpa_calcmill(iq,ikst,ikend,isym,iop_minm)
      use crpa,      only: nlmorb,mill,nbk_wf,pln,wf_centers_lm,&
     &                     iop_pln_phase
      use bands,     only: nspin
      use constants, only: cone,pi
      use kpoints,   only: nirkp,idikp,nkp,kqid,kpirind,wkir,klist,& 
     &                     qlist,idvk 
      use bzinteg,   only: kwt_bz,kwt_ibz
      use mixbasis,  only: matsiz

      implicit none
      integer,intent(in):: iq,ikst,ikend
      integer,intent(in):: isym
      integer,intent(in):: iop_minm

      integer:: isp,ib0,ib1,irk,ik,jk,jrk,iik
      integer:: nst,nend,mst,mend,nmdim,inm,ie1,ie2
      integer:: ill,lldim,ilm1,ilm2
      integer:: ierr
      real(8):: arg
      real(8):: kvec(3), kqvec(3) 
      complex(8):: coef, phs 
      complex(8),allocatable:: minm(:,:,:),tmat(:,:)

      logical:: ldbg=.true.
      character(len=20):: sname="crpa_calcmill"

      mill = 0.0 
      do isp=1,nspin 
        do iik=ikst,ikend
          ! set ik, irk, jk and jrk
          if(isym.eq.0) then
            ik=iik
            irk=kpirind(ik)
            coef= kwt_bz(ik)
          else
            irk=iik
            ik=idikp(irk)
            coef= kwt_ibz(irk)
          endif
          jk=kqid(ik,iq)
          jrk=kpirind(jk)

          kvec = dble(klist(:,ik))/idvk 
          !kqvec= dble(klist(:,ik)-qlist(:,iq))/idvk 
          kqvec = dble(klist(:,jk))/idvk 

          nst  = nbk_wf(1,ik,isp)
          nend = nbk_wf(2,ik,isp)
          mst  = nbk_wf(1,jk,isp)
          mend = nbk_wf(2,jk,isp)
          nmdim= (nend-nst+1)*(mend-mst+1)
          lldim= nlmorb*nlmorb
          allocate(minm(matsiz,mst:mend,nst:nend), &
     &             tmat(nmdim,lldim), stat=ierr)
          call errmsg(ierr.ne.0,sname,"Fail to allocate minm and tmat")

          call get_minm(iop_minm,'nm',minm,nst,nend,mst,mend,ik,iq,isp)
 
          ill=0
          do ilm2=1,nlmorb
            do ilm1=1,nlmorb
              ill = ill+1
            
              if(iop_pln_phase.eq.2) then 
                arg = -sum(wf_centers_lm(:,ilm1)*kvec(:))  &
     &                +sum(wf_centers_lm(:,ilm2)*kqvec(:))   
                phs = cmplx(cos(2.0d0*pi*arg),sin(2.0d0*pi*arg),8)
                !if(ldbg) then 
                !  write(6,*) "crpa_calcmill: phase_factor=",phs 
                !endif 
              else 
                phs = 1.0 
              endif 
              
              inm=0
              do ie1=nst,nend 
                do ie2=mst,mend
                  inm = inm+1
                  tmat(inm,ill) = conjg(pln(ilm1,ie1,ik,isp)) &
     &                       *pln(ilm2,ie2,jk,isp)*phs 

                enddo
              enddo
            enddo
          enddo

          call zgemm('n','n',matsiz,lldim,nmdim,coef,minm,matsiz,&
     &            tmat,nmdim,cone,mill(:,:,isp),matsiz)
          deallocate(minm,tmat)
        enddo ! loop over iik
      enddo ! loop over ispin
      
      end subroutine
