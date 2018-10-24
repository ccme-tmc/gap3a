!BOP
! 
! !ROUTINE: calcloabc
!
! !INTERFACE:
      subroutine calcloabc(iat,l,ilo,lapw,isp)

! !DESCRIPTION:
!
!   calcloabc calculates the cofficients a,b,c of the local orbitals (lo)
!
! !USES:
      use lapwlo, only : abcelo,umt,umtlo              
      use struk,  only : rmt
      use task,   only : fid_outmb
      
! !INPUT PARAMETERS:

      implicit none

      integer(4), intent(in) :: iat
      integer(4), intent(in) :: l
      integer(4), intent(in) :: ilo
      logical,    intent(in) :: lapw
      integer(4), intent(in) :: isp

! !DEFINED PARAMETERS:
      real(8), parameter :: cutoff = 200.0d+0

! !LOCAL VARIABLES:  
      real(8) :: xac
      real(8) :: xbc
      real(8) :: alonorm
      real(8) :: plo, dplo, pe12lo,pi12lo, e,p,pe,dp,dpe,pei,alo,blo,clo
!
! !INTRINSIC ROUTINES:
!
      intrinsic          sqrt
!EOP
!BOC
      

      p=umt(2,l,iat,isp) 
      pe=umt(3,l,iat,isp)
      dp=umt(4,l,iat,isp)
      dpe=umt(5,l,iat,isp)
      pei=umt(6,l,iat,isp) 
   
      if(ilo.eq.0) then  ! APW+lo
        alonorm=sqrt(1.d0+(p/pe)**2 * pei)   
        alo = 1.d0 /alonorm 
        blo = -p/pe/alonorm
        clo = 0.d0

      else
        plo    = umtlo(1,ilo,l,iat,isp) 
        dplo   = umtlo(2,ilo,l,iat,isp)
        pi12lo = umtlo(3,ilo,l,iat,isp) 
        pe12lo = umtlo(4,ilo,l,iat,isp) 

        if(lapw) then ! LAPW+LO
          xac = plo*dpe -  dplo*pe            
          xac = xac*rmt(iat)*rmt(iat)
          xbc = plo*dp - dplo*p
          xbc = -xbc*rmt(iat)*rmt(iat)
          clo = xac*(xac+2.0*pi12lo)+xbc*(xbc*pei+2.0*pe12lo)+1.0
          clo = 1.0/sqrt(clo)
          clo = min(clo,cutoff)
          alo = clo*xac
          blo = clo*xbc

        else ! APW+lo + LO
          xbc=-p/plo
          xac=sqrt(1.0+xbc**2+2*xbc*pi12lo)
          alo = 1.d0/xac
          blo = 0.d0
          clo = xbc/xac
        endif

      endif

      abcelo(1,ilo,l,iat,isp) = alo
      abcelo(2,ilo,l,iat,isp) = blo
      abcelo(3,ilo,l,iat,isp) = clo

      write(fid_outmb,6000) l,ilo,lapw,abcelo(1:3,ilo,l,iat,isp) 
      return
!
 6000 format ('lo coefficient: l,ilo,lapw,a,b,c  ',i2,2x,i2,2x,l1,5x, &
     &        3f12.5)
      end subroutine calcloabc
!EOC
