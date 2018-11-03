      subroutine crpa_calcumat(iq, iomfirst, iomlast)
      use bands,     only: nspin 
      use bzinteg,   only: kwt_bz,singc1,singc2
      use freq,      only: iop_freq
      use kpoints,   only: qlist,idvq
      use constants, only: cone,czero,fourpi,pi
      use crpa,      only: umat,nlmorb,mill,ncell,rcell
      use dielmat,   only: eps,head,epsw1,epsw2
      use mixbasis,  only: matsiz
      use struk,     only: vi
      implicit none
      integer,intent(in):: iq, iomfirst, iomlast

      integer:: ierr       
      integer:: i,icell 
      integer:: is1,is2,is12         
      integer:: iom                 !! index for frequency
      integer:: ill1,ill2           !! index for (\lambda1,\lambda2)
      integer:: ilm1,ilm2,ilm3,ilm4 !! index for (a,l,m)
      integer:: lldim               !! dimension of (\lambda1,\lambda2)
      real(8):: d42,d13 
      complex(8):: phs
      real(8):: qvec(3),arg
      complex(8):: coef,coefs2,coefs1
      logical:: ldbg=.false.
      complex(8),allocatable:: tmat(:,:),w2m(:),mw1(:)
      character(20):: sname="crpa_calcumat"

      lldim = nlmorb*nlmorb
      coefs2= cone*singc2*fourpi*vi
      coefs1= cone*singc1*sqrt(fourpi*vi)

      allocate(tmat(matsiz,lldim),w2m(lldim),mw1(lldim),stat=ierr)

      do i=1,3
        qvec(i)=dble(qlist(i,iq))/dble(idvq)
      enddo

      call errmsg(ierr.ne.0,sname,"Fail to allocate tmat,em")

      do icell=1,ncell 
        arg=sum( (rcell(:,icell) - rcell(:,1))*qvec(:))
        phs=cmplx(cos(2.0d0*pi*arg),-sin(2.0d0*pi*arg),8)

        coef  = cone*kwt_bz(iq)*phs   !! weight for the summation over q

        do iom=iomfirst,iomlast
          is12 = 0 
          do is2=1,nspin 
            if(ldbg) write(6,*) "calc eps*mill"
            if(iop_freq.eq.3) then  !! imag freq
              call zhemm('l','u',matsiz,lldim,cone,eps(:,:,iom),matsiz,&
     &          mill(:,:,is2),matsiz,czero,tmat,matsiz)
            else 
              call zgemm('n','n',matsiz,lldim,matsiz,cone,eps(:,:,iom),&
     &          matsiz,mill(:,:,is2),matsiz,czero,tmat,matsiz)
            endif 
        
            if(iq.eq.1) then 
              !! prepare for the wing terms
              if(ldbg) write(6,*) "calc epsw2*mill"
              call zgemv('t',matsiz,lldim,cone,mill(:,:,is2),matsiz, &
     &            epsw2(:,iom),1,czero,w2m,1)
            endif 

            do is1=1,is2

              is12 = is12 + 1

              if(ldbg) write(6,*) "calc mill*(eps*mill)"

              call zgemm('c','n',lldim,lldim,matsiz,coef,mill(:,:,is1),&
     &         matsiz,tmat,matsiz,cone,umat(:,:,is12,iom,icell),lldim) 

              if (iq .ne. 1) cycle   

              !! iq=1 singularity

              if(ldbg) write(6,*) "calc mill*epsw1"

              call zgemv('c',matsiz,lldim,cone,mill(:,:,is1),matsiz, &
     &            epsw1(:,iom),1,czero,mw1,1)

              if(ldbg) write(6,*) "add singular contributions"

              ill2=0
              do ilm2=1,nlmorb
                do ilm4=1,nlmorb 
                  ill2 = ill2 + 1
                  d42 = f_delta(ilm4,ilm2) 

                  ill1 = 0
                  do ilm3=1,nlmorb
                    do ilm1=1,nlmorb
                      ill1 = ill1 + 1
                      d13 = f_delta(ilm1,ilm3) 
                      umat(ill1,ill2,is12,iom,icell) = &
     &                 umat(ill1,ill2,is12,iom,icell)&
     &                 + coefs2*head(iom)*d13*d42                        &
     &                 + coefs1*(d13*w2m(ill2) + mw1(ill1)*d42) 
                    enddo 
                  enddo
                enddo ! ilm4
              enddo ! ilm2

            enddo ! is1
          enddo  ! is2
        enddo ! iomeg
      enddo ! icell 
      deallocate(tmat,w2m,mw1)

      contains 
        real(8)  function f_delta(i1,i2)
        integer:: i1,i2
        if(i1.eq.i2) then
          f_delta = 1.0
        else
          f_delta = 0.0
        endif 
        end function 
      end subroutine
