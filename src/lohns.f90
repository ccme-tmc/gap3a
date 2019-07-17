!BOP
!
! !ROUTINE: lohns
!
! !INTERFACE:
      subroutine lohns(jneq,mult,i1,i2,l0,jlo)

! !DESCRIPTION:
!
!        lohns calculates indices of loop in hns
!
! !USES:

      use lapwlo,only: nlo,lomax
      use struk, only: nat
      
! !INPUT PARAMETERS:

      implicit none
      integer            i1, i2, jneq, l0
      integer            mult(nat)

! !LOCAL VARIABLES:
!
      integer            i,l
      integer            jlo, jlo1
!
!EOP
!BOC
      i2 = 1
      do i = 1, jneq - 1
        do l = 0, lomax
          do jlo1 = 1,nlo(l,i)
            i2 = i2 + (2*l+1)*mult(i)
          enddo
        enddo
      enddo 
      do l = 0, l0-1
        do jlo1 = 1,nlo(l,jneq)
          i2 = i2 + (2*l+1)*mult(jneq)
        enddo
      enddo 

      do jlo1 =1,jlo
         i1 = i2
         i2 = i2 + (2*l+1)*mult(jneq)
      enddo

      i2 = i2 - 1

      return

      end subroutine lohns
!EOC
