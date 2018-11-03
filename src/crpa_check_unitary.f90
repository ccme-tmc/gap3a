      subroutine crpa_chk_unitary(nkp,nsp) 
      use crpa, only: pln,nlmorb,nbk_wf 
      implicit none 
      integer,intent(in)::nkp,nsp

      integer:: ik,isp
      real(8):: err_ortho,s,err_max
      integer:: ib,ilm1,ilm2,nb1,nb2
      character(20):: sname='crpa_chk_unitary' 

      call linmsg(6,'-',sname) 
      err_max = 0.d0 
      do isp=1,nsp
        do ik=1,nkp 
          nb1 = nbk_wf(1,ik,isp)
          nb2 = nbk_wf(2,ik,isp)
          err_ortho  = 0.d0
          do ilm2=1,nlmorb
            do ilm1=1,nlmorb
              s = 0.d0
              do ib=nb1,nb2
                s = s+ conjg(pln(ilm1,ib,ik,isp))*pln(ilm2,ib,ik,isp)
              enddo
              if(ilm1.eq.ilm2) then
                err_ortho = err_ortho + abs(s-1.0)
              else
                err_ortho = err_ortho + abs(s)
              endif
            enddo ! ilm1
          enddo ! ilm2
          if(err_ortho.gt.1.e-6) then
            write(6,100) ik,err_ortho
          endif
          err_max=max(err_ortho,err_max) 
        enddo ! ik
      enddo ! isp
      
      if(err_max.gt.1.e-4) then
        write(6,*) "WARNING: Significant unitarity error found!"
      endif  
         
100   format("Pln unitarity error for ik=",i4,f12.6) 
      end subroutine

