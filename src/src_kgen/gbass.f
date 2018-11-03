
!     .....................................................GBASS........
      SUBROUTINE GBASS(RBAS,GBAS)                                       
!     **                                                              **
!     **  CALCULATE RECIPROCAL LATTICE VECTORS FROM REAL SPACE        **
!     **  LATTICE VECTORS OR VICE VERSA                               **
!     **                                                              **
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION RBAS(3,3),GBAS(3,3)                                     
      PI=4.D0*DATAN(1.D0)                                               
      GBAS(1,1)=RBAS(2,2)*RBAS(3,3)-RBAS(3,2)*RBAS(2,3)                 
      GBAS(2,1)=RBAS(3,2)*RBAS(1,3)-RBAS(1,2)*RBAS(3,3)                 
      GBAS(3,1)=RBAS(1,2)*RBAS(2,3)-RBAS(2,2)*RBAS(1,3)                 
      GBAS(1,2)=RBAS(2,3)*RBAS(3,1)-RBAS(3,3)*RBAS(2,1)                 
      GBAS(2,2)=RBAS(3,3)*RBAS(1,1)-RBAS(1,3)*RBAS(3,1)                 
      GBAS(3,2)=RBAS(1,3)*RBAS(2,1)-RBAS(2,3)*RBAS(1,1)                 
      GBAS(1,3)=RBAS(2,1)*RBAS(3,2)-RBAS(3,1)*RBAS(2,2)                 
      GBAS(2,3)=RBAS(3,1)*RBAS(1,2)-RBAS(1,1)*RBAS(3,2)                 
      GBAS(3,3)=RBAS(1,1)*RBAS(2,2)-RBAS(2,1)*RBAS(1,2)                 
      DET=0.D0                                                          
      DO 100 I=1,3                                                      
      DET=DET+GBAS(I,1)*RBAS(I,1)                                       
100   CONTINUE                                                          
      DO 110 I=1,3                                                      
      DO 110 J=1,3                                                      
      GBAS(I,J)=GBAS(I,J)*2.D0*PI/DET                                   
110   CONTINUE                                                          
      RETURN                                                            
      END                                                               
