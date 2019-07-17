!BOP
!
! !ROUTINE: fermi
!
! !INTERFACE:
      subroutine fermi(nik,nb,eb,nel,nsp,efer,eg,info)
!
! !DESCRIPTION:
!  This subroutine calculated the Fermi energy with the tetrahedron method
!
! !USES:
      use bzint, only: nirtet,tndi,wirtet,tvol 
      implicit none     
      
! !INPUT PARAMETERS:

      integer, intent(in) :: nik                 ! Number of irreducible k-points
      integer, intent(in) :: nb                  ! Maximum number of bands
      integer, intent(in) :: nsp                 ! =1 / 2 for spin unpolarized/polarized caculations 
      real(8), intent(in) :: eb(nb,nik,nsp)      ! Band energies
      real(8), intent(in) :: nel                 ! number of electrons
      character(len=*),intent(in):: info 
      
! !OUTPUT PARAMETERS:      
      
      real(8), intent(out)   :: efer       ! the fermi energy
      real(8), intent(out)   :: eg 

! !REVISION HISTORY:
!
! Created 10th. March 2004 by RGA
! Modified by Jiang
!
! !LOCAL VARIABLES:

      integer:: ik, ib, it,isp,vbm(2),cbm(2)
      integer:: nvm

      real(8):: evbm,ecbm,eint,ocmin,ocmax,ocint,df 
      
      real(8), external :: idos
      real(8), external :: dostet
      logical:: nfer

      logical::lprt=.false.
      integer::nitmax=100000
      real(8)::eps=1.d-8
      real(8), parameter :: HeV = 27.2113961
!EOP
!BOC

      if(lprt) write(6,*) " ---- Calculate Fermi  energy ----" 

      !! nvm is the number of bands for an insulating system 
      !! since for a system with gap, the procedure to determine the
      !! band gap can be unstable, just try first whether it is an
      !! insulating system, but such a simplistic way to determine the Fermi energy
      !! is valid only for spin un-polarized cases 
      if(nsp.eq.1) then 
        nvm  = nint(nel/2.d0) 
        evbm = maxval(eb(nvm,:,1))
        ecbm = minval(eb(nvm+1,:,1))
        eint= (evbm+ecbm)/2.d0
        ocint= idos(nsp,nik,nb,eb,nirtet,tndi,wirtet,tvol,eint)
        if(ecbm.ge.evbm.and.abs(ocint-nel).lt.eps) then 
          efer = eint
          eg = ecbm - evbm
          return 
        endif 
      endif 
     
! find the minimal and maximal band energy  
      evbm=minval(eb)
      ecbm=maxval(eb,mask=eb.lt.1.0e3)
     
      ocmin=idos(nsp,nik,nb,eb,nirtet,tndi,wirtet,tvol,evbm)
      ocmax=idos(nsp,nik,nb,eb,nirtet,tndi,wirtet,tvol,ecbm)

      if(ocmax.le.nel) then 
        write(6,'(a)')  'ERROR in fermi: not enough bands'
        write(6,'(a,f10.4,2f10.2)') ' emax,ocmax,nel= ',ecbm,ocmax,nel 
        stop
      endif 

      if(lprt)  write(6,1) '#it',"evbm","ecbm","eint","ocmin","ocmax",  &
     &        "ocint"
!
! Use bisection method to determine solver the equation N( Efermi ) = Nel
!
      nfer=.true.      
      do it=1,nitmax
        eint=evbm+(ecbm-evbm)/2
        
        ocint=idos(nsp,nik,nb,eb,nirtet,tndi,wirtet,tvol,eint)

        if(lprt)  write(6,2) it,evbm,ecbm,eint,ocmin,ocmax,ocint

        if(abs(ocint-nel).lt.eps) exit  
       
        if(ocint.gt.nel)then
          ecbm=eint
          ocmax=ocint
        else
          evbm=eint
          ocmin=ocint
        endif
      enddo

      if(it.ge.nitmax) then 
        write(6,*)  "ERROR in fermi: fail to converge"
        stop
      endif 

      df=dostet(nsp,nik,nb,eb,nirtet,tndi,wirtet,tvol,eint)
!
! For insulator (including semiconductor, set fermi energy as the middle of gap
!
      
      if(df.lt.1.0d-4) then
        if(nsp.eq.1) then 
          evbm=maxval(eb(:,:,1),mask=eb(:,:,1).lt.eint)
          ecbm=minval(eb(:,:,1),mask=eb(:,:,1).gt.eint)
        else 
          evbm=max( maxval(eb(:,:,1),mask=eb(:,:,1).lt.eint), &
     &              maxval(eb(:,:,2),mask=eb(:,:,2).lt.eint))
          ecbm=min( minval(eb(:,:,1),mask=eb(:,:,1).gt.eint), &
     &              minval(eb(:,:,2),mask=eb(:,:,2).gt.eint))   
        endif 

        efer=0.5d0*(evbm+ecbm)
        eg=ecbm-evbm
      else 
        efer=eint 
        eg=-df
      endif 
      if(lprt) write(6,'(a,f12.6)') "EFermi (Ha)=",efer 

      if(trim(info).ne.'none') then 
        write(6,100)   ":E_FERMI_"//trim(info)//"(eV)=  ",efer*Hev
        if(eg.gt.0.d0) then
          write(6,100) ':BandGap_'//trim(info)//'(eV)=  ',eg*Hev
        else
          write(6,100) ':DOS_at_Fermi= ',-eg
        endif
      endif 


  1   format(A5,6A12)
  2   format(I5,6f12.6)
  3   format(A,2f12.6) 
  4   format(a10,2i4,f12.4)
 100  format(a,f12.4)
      end subroutine fermi
!EOC        
          

      
      
            



      
