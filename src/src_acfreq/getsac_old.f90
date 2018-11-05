      subroutine getsac(iop,nomeg,npar,en,ein,omega,apar,sc,dsc) 

      implicit none 
      integer(4),intent(in)::iop,nomeg,npar
      real(8),intent(in)::en
      real(8),intent(in)::omega(nomeg)
      complex(8),intent(in)::ein,apar(npar) 
      complex(8),intent(out)::sc,dsc

      complex(8)::comega(npar)
      integer(4)::ip,iw,iwpa(npar)

      if(iop.eq.1) then 
        call acrgn(npar,ein,apar,sc,dsc)
      elseif(iop.eq.0) then   
        call setwpa(nomeg,npar,omega,iwpa)
        do ip=1,npar
          iw=iwpa(ip)
          comega(ip)=cmplx(0.d0,dsign(omega(iw),en))
        enddo 
        call acpatrd(npar,ein,comega,apar,sc,dsc)
      elseif(iop.eq.2) then  
         call setwpa(nomeg,npar,omega,iwpa)
         do ip=1,npar
           iw=iwpa(ip)
           comega(ip)=cmplx(-en,dsign(omega(iw),en))
         enddo
         call acpatrd(npar,ein-en,comega,apar,sc,dsc)
      endif 

      end subroutine 
