      subroutine calceps_dmft(iq,iom_f,iom_l)
      use dielmat,     only: eps
      use crpa,        only: mill,ploc_dmft

! !DESCRIPTION:
!
! This subroutine calculates add the dielectric matrix with a local
! correction that is obtained from the DMFT impurity solver 
!

       
