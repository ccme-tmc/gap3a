      subroutine sdefl(rbas,gbas,iio,nsym,iz,iord,lattic,ortho)
!
!    redefines symmetry operation from lapw-struct with
!    unitary transformation  u(-1) * S * u
!
      IMPLICIT REAL*8 (A-H,O-Z)                         
      LOGICAL           ORTHO
      character*4 lattic                
      dimension iio(3,3,48),iz(3,3,48),rbas(3,3),gbas(3,3),a(3,3),b(3,3)
      dimension gbas1(3,3)
      pi=acos(-1.d0)
      if(iord.ne.nsym) then 
!      write(66,*) 'iord (from lapw) .ne. nsym (from kgenb)',iord,nsym
      nsym=iord
      end if
!    
      do 5 i=1,3
      do 5 j=1,3          
      gbas(i,j)=gbas(i,j)/2.d0/pi
  5   gbas1(j,i)=gbas(i,j)
!            write(6,111) rbas,gbas1
! 
      do 1 ind=1,iord 
         do 2 i=1,3
         do 2 j=1,3
  2      a(i,j)=iz(i,j,ind)
         if(ortho) then
             call matmm(b,rbas,a)
             call matmm(a,b,gbas1)
         end if
         if(.not.ortho.and.lattic(1:3).eq.'CXZ') then
!            write(6,111) a
             call matmm(b,rbas,a)
             call matmm(a,b,gbas1)
         end if
         do 3 i=1,3
         do 3 j=1,3
!  3      iio(j,i,ind)=nint(a(i,j))
  3      iio(i,j,ind)=nint(a(i,j))
!      write(6,111) (a(1,j),j=1,3)
!      write(6,111) (a(2,j),j=1,3)
!      write(6,111) (a(3,j),j=1,3)
111   format(3f10.4)
1     continue
      DO 130 I=1,INT(FLOAT(NSYM)/4.+.9)                                 
      I1=4*I-3                                                          
      I2=4*I-2                                                          
      I3=4*I-1                                                          
      I4=4*I                                                            
      WRITE (66,120) I1,I2,I3,I4                                        
  120 FORMAT (T5,'SYMMETRY MATRIX NR.',I3,T30,'SYMMETRY MATRIX NR.'      &
       ,I3,T55,'SYMMETRY MATRIX NR.',I3,T80,'SYMMETRY MATRIX NR.',I3)   
      DO 130 J=1,3                                                      
      WRITE (66,140) (Iio(J,K,I1),K=1,3),(IiO(J,K,I2),K=1,3), &
       (IiO(J,K,I3),K=1,3),(IiO(J,K,I4),K=1,3)                          
  130 CONTINUE                                                          
  140 FORMAT (T5,3I5,T30,3I5,T55,3I5,T80,3I5)                           
      return
      end
