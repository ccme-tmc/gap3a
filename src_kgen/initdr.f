!     ..................................................................
      SUBROUTINE INITDR(NAME,LENG)                                      
      DOUBLE PRECISION NAME(LENG)                                       
      DO 100 I=1,LENG                                                   
      NAME(I)=0.D0                                                      
100   CONTINUE                                                          
      RETURN                                                            
      END                                                               
