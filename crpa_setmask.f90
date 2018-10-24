      subroutine crpa_setmask(irk,jrk,isp)
      use bands,    only: bande,nomaxs,numins,nbmaxpol,nspin
      use dielmat,  only: iop_mask_eps, mask_eps, noc_excl,nun_excl, &
     &                    iun_excl,ioc_excl,occ_win, unocc_win, wt_excl
      use crpa,     only: weight_corr, nbmin_wf, nbmax_wf
      use core,     only: ncg_p 
      implicit none 
      integer:: irk,jrk,isp
      integer:: nomx,numn,ie1,ie2,i
      logical:: lexcl_ie1, lexcl_ie2
      real(8):: e1, e2, wt

      !! imask_eps is used to skip some transitions when calculating eps 
      !! the skipped transitions can be defined by ioc_excl(:) and iun_excl(:)
      !! but it can also be defined in terms of some energy window 
      mask_eps = 1.0

      if(iop_mask_eps.eq.0) return 

      nomx = nomaxs(isp)
      numn = numins(isp)
      if (iop_mask_eps.eq.1) then !! use the band index to define the mask window
        if(noc_excl.gt.0.and.nun_excl.gt.0) then 
          do ie1 = 1, nomx
            lexcl_ie1 = .false.
            do i = 1, noc_excl
              if (ie1 == ioc_excl(i)) then
                lexcl_ie1 = .true.
                exit
              endif
            enddo
            if (.not.lexcl_ie1) cycle

            do ie2 = numn, nbmaxpol
              lexcl_ie2 = .false.
              do i = 1, nun_excl
                if (ie2 == iun_excl(i)) then
                  lexcl_ie2 = .true.
                  exit
                endif
              enddo
              if (.not.lexcl_ie2) cycle 
              mask_eps(ie2,ie1+ncg_p) = wt_excl 
            enddo ! ie2
          enddo  ! ie1

        else if(noc_excl.gt.0.and.nun_excl.eq.0) then 

          do ie1 = 1, nomx
            do i = 1, noc_excl
              if (ie1 == ioc_excl(i)) then
                mask_eps(:,ie1+ncg_p) = wt_excl 
                exit 
              endif
            enddo
          enddo

        elseif(noc_excl.eq.0.and.nun_excl.gt.0) then 

          do ie2 = numn, nbmaxpol
            do i = 1, nun_excl
              if (ie2 == iun_excl(i)) then
                mask_eps(ie2,:) = wt_excl 
                exit
              endif
            enddo
          enddo

        endif 

      elseif (iop_mask_eps.eq.2) then 

        do ie2=numn,nbmaxpol
          if(ie2<nbmin_wf.or.ie2>nbmax_wf) cycle
          do ie1=1,nomx
            if(ie1<nbmin_wf.or.ie1>nbmax_wf) cycle
            wt = weight_corr(ie1,irk,isp)*weight_corr(ie2,jrk,isp)
            mask_eps(ie2,ie1+ncg_p) = 1.0 - wt
            write(6,101) ie1,ie2,mask_eps(ie2,ie1)
          enddo
        enddo

      elseif (iop_mask_eps.eq.-1) then  !! determine the mask in terms of energy 
        do ie1=1,nomx
          e1=bande(ie1,irk,isp)
          if(e1<occ_win(1).or.e1>occ_win(2)) cycle 
          do ie2=numn,nbmaxpol
            e2=bande(ie2,jrk,isp)
            if(e2<unocc_win(1).or.e2>unocc_win(2)) cycle   
            mask_eps(ie2,ie1+ncg_p) = wt_excl 
            !write(6,102) irk,ie1,e1,jrk,ie2,e2,wt_excl 
          enddo
        enddo 
      endif 
  100 format(" weitht transition: ",i4," -> ", i4)
  101 format(" Weight the transition: ",i4," -> ", i4, " by ", f10.4)
  102 format(" Weight the transition: irk=",i4,' n=',i4,' e1=',f12.4,  &
     & " --->  irk=",i4,' m=', i4, ' e2=',f12.4," by ", f10.4)
      end subroutine 

