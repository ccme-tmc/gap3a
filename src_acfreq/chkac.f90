      subroutine chkac(en,sc,omega,nomeg,npar,apar,scfit,dscfit)
      implicit none 
      integer(4),intent(in):: npar,nomeg,ne
      real(8),intent(in) ::en(ne),omega(nomeg)
      complex(8),intent(in):: sc(nomeg),apar(npar)
      complex(8),intent(out)::scfit(nomeg),dscfit(nomeg) 
      
      integer(4)::iom

      do iom=1,nomeg
        call ratfc(0.d0,omega(iom),apar,scfit(iom),dscfit(iom),npar)
