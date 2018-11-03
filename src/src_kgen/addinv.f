      subroutine addinv(iord,iz,is)
!
!    adds inversion symmetry if not present (do not use in SO calculations
      IMPLICIT REAL*8 (A-H,O-Z)                         
      dimension iz(3,3,48)
! 
      do 1 ind=1,iord 
      if(iz(1,1,ind).eq.-1.and.iz(2,2,ind).eq.-1.and.iz(3,3,ind) &
      .eq.-1.and.iz(1,2,ind).eq.0.and.iz(1,3,ind).eq.0.and. &
      iz(2,1,ind).eq.0.and.iz(2,3,ind).eq.0.and.iz(3,1,ind).eq.0 &
      .and.iz(3,2,ind).eq.0) return
 1    continue
      if(iord.gt.24) stop 'iord gt 24 without inversion'
!
      write(*,*) iord,' symmetry operations without inversion' 
      write(*,*) 'Do you want to add inversion: No=0, Yes=1 (0/1)'
      read(*,*) is

      if(is.eq.1) then
        write(*,*) 'inversion added (non-spinpolarized non-so calculation)'
        do 3 ind=1,iord
         do 2 i=1,3
         do 2 j=1,3
  2      iz(i,j,ind+iord)=-iz(i,j,ind)
 3      continue
        iord=iord*2
      endif
      return
      end
