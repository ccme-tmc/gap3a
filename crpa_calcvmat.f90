      subroutine crpa_calcvmat(iq)
!
!     This subroutine calculates the matrix elements of bare Coulomb
!     interaction as represented by the Wannier functions 
! 
      use crpa,      only: vmat,mill,nlmorb,ncell,rcell
      use bzinteg,   only: kwt_bz,singc1,singc2
      use constants, only: cone,czero,fourpi,pi
      use kpoints,   only: qlist, idvq
      use mixbasis,  only: matsiz
      use struk,     only: vi
      use bands,     only: nspin

      implicit none
      integer,intent(in):: iq

      integer:: i
      integer:: ill1,ill2           !! index for (\lambda1,\lambda2)
      integer:: is1,is2,is12,icell
      integer:: spin 
      integer:: ilm1,ilm2,ilm3,ilm4 !! index for (a,l,m)
      real(8):: d42,d13 
      real(8):: qvec(3),arg
      complex(8):: phs 
      complex(8):: coef,coefs2
      complex(8),external:: zdotc
      character(20):: sname="crpa_calcvmat"
      logical:: ldbg=.false.

      do i=1,3
        qvec(i)=dble(qlist(i,iq))/dble(idvq)
      enddo
      if(ldbg) then 
        write(6,100) qvec
      endif 
 100  format("q=",3f8.4)

      coefs2= cone*singc2*fourpi*vi
      if(iq.ne.1) coefs2 = 0.d0 
    
      do icell = 1, ncell
        arg=sum((rcell(:,icell)-rcell(:,1))*qvec(:))
        phs=cmplx(cos(2.0d0*pi*arg),-sin(2.0d0*pi*arg),8)
        if(ldbg) then 
          write(6,*) "arg,phs=",arg,phs
        endif 

        coef  = cone*kwt_bz(iq)*phs    !! weight for the summation over q including the phase factor 

        is12 = 0 
        do is2=1,nspin
          do is1=1,is2
            is12 = is12 + 1
 
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
                    vmat(ill1,ill2,is12,icell)=vmat(ill1,ill2,is12,icell)     &
     &               + coef*zdotc(matsiz,mill(:,ill1,is1),1,mill(:,ill2,is2),1) & 
                     + coefs2*d13*d42
                  enddo 
                enddo
              enddo ! ilm4
            enddo ! ilm2
          enddo  ! is1
        enddo ! is2
      enddo 

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
