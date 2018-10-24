!     .....................................................ORD1.........
      SUBROUTINE ORD1(NMAX,IX,WORK)                                     
!     **                                                              **
!     **  ORD1 ORDERES THE ARRAY IX WITH SIZE NMAX ACCORDING TO       **
!     **  INCREASING NUMBER                                           **
!     **                                                              **
!     **  SUBROUTINES USED:                                           **
!     **    indexx (from Numerical recipes)                           **
!     **                                                              **
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)                  
      INTEGER IX(4,NMAX),WORK(NMAX)                          
      parameter (indexm=1200000)
      common /mist/ index(indexm)                                       
        if(nmax.gt.indexm) then
        write(*,*) 'nmax too large in ord1.f, set indexm at least ',nmax
        stop 'nmax, redimension ord1.f'
        end if
!
!     sort first index of ix
      do 300 i=1,nmax
 300  work(i)=ix(1,i)
      call indexx(work,index,nmax)
      do 310 i=1,nmax
 310  ix(1,i)=work(index(i))
      do 320 i=1,nmax
 320  work(i)=ix(2,i)
      do 325 i=1,nmax
 325  ix(2,i)=work(index(i))
      do 330 i=1,nmax
 330  work(i)=ix(3,i)
      do 335 i=1,nmax
 335  ix(3,i)=work(index(i))
      do 340 i=1,nmax
 340  work(i)=ix(4,i)
      do 345 i=1,nmax
 345  ix(4,i)=work(index(i))
! 
!     sort higher indices with trivial procedure
 501  ichang=0
      do 500 i=2,nmax
      if(ix(1,i).eq.ix(1,i-1)) then
          if(ix(2,i).lt.ix(2,i-1)) then
          i2=ix(2,i-1)
          i3=ix(3,i-1)
          i4=ix(4,i-1)
          ix(2,i-1)=ix(2,i)
          ix(3,i-1)=ix(3,i)
          ix(4,i-1)=ix(4,i)
          ix(2,i)=i2
          ix(3,i)=i3
          ix(4,i)=i4
          ichang=1
          else if(ix(2,i).eq.ix(2,i-1)) then
            if(ix(3,i).lt.ix(3,i-1)) then
            i3=ix(3,i-1)
            i4=ix(4,i-1)
            ix(3,i-1)=ix(3,i)
            ix(4,i-1)=ix(4,i)
            ix(3,i)=i3
            ix(4,i)=i4
            ichang=1
            else if(ix(3,i).eq.ix(3,i-1)) then
              if(ix(4,i).lt.ix(4,i-1)) then
              i4=ix(4,i-1)
              ix(4,i-1)=ix(4,i)
              ix(4,i)=i4
              ichang=1
              end if
            end if
          end if
      end if
 500  continue
      if(ichang.eq.1) goto 501
      RETURN                                                            
      END                                                               
