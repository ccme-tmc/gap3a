!     .....................................................TETDIV...... 
      SUBROUTINE TETDIV(N,GBAS,TET0)                                    
!     **                                                             ** 
!     **  TETDIV DETERMINES THE DIVISION OF THE PARALLELEPIPEDS      ** 
!     **  IN TETRAHEDRONS ACCORDING TO THE SHORTEST DIAGONAL         ** 
!     **                                                             ** 
!     **  INPUT:                                                      **
!     **    N           NUMBER OF DIVISIONS OF REC. LATTICE VECTORS   **
!     **    GBAS        RECIPROCAL LATTICE VECTORS                    **
!     **  OUTPUT:                                                    ** 
!     **    TET0(I,J,K) COORDINATES (I) IN THE BASIS OF SUBLATTICE   ** 
!     **                VECTORS, OF THE 4 CORNERS (J) OF A           ** 
!     **                TETRAHEDRON (K) IN A SUBLATTICE UNIT CELL    ** 
!     **                                                             ** 
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)                               
      DIMENSION P(8,3),DIAG(4),GBAS(3,3)                                
      INTEGER IACHT(8),TET(4,6),N(3),TET0(3,4,6)                        
!     ------------------------------------------------------------------
!     -- SEARCH FOR THE SHORTEST DIAGONAL                             --
!     ------------------------------------------------------------------
      DO 10 I=0,1                                                       
      DO 10 J=0,1                                                       
      DO 10 K=0,1                                                       
      ISVAR=4*I+2*J+K+1                                                 
      DO 10 L=1,3                                                       
!      P(ISVAR,L)=GBAS(L,1)*DBLE(I)/DBLE(N(1))                           
!     1          +GBAS(L,2)*DBLE(J)/DBLE(N(2))                           
!     1          +GBAS(L,3)*DBLE(K)/DBLE(N(3))                           
      P(ISVAR,L)=GBAS(1,L)*DBLE(I)/DBLE(N(1))                           &
                +GBAS(2,L)*DBLE(J)/DBLE(N(2))                           &
                +GBAS(3,L)*DBLE(K)/DBLE(N(3))                          
10     CONTINUE                                                         
!*                                                                      
       DO 20 I=1,4                                                      
       DIAG(I)=0.D0                                                     
       DO 20 J=1,3                                                      
       DIAG(I)=DIAG(I)+(P(I,J)-P(9-I,J))**2                             
20     CONTINUE                                                         
!*                                                                      
       MNDG=1                                                           
       DO 30 I=2,4                                                      
       IF(DIAG(I).LT.DIAG(MNDG)) THEN                                   
         MNDG=I                                                         
       END IF                                                           
30     CONTINUE                                                         
!     ------------------------------------------------------------------
!     -- ROTATE PARALLELEPIPED                                        --
!     ------------------------------------------------------------------
      IF(MNDG.EQ.1)THEN                                                 
        DO 40 I=1,8                                                     
        IACHT(I)=I                                                      
40      CONTINUE                                                        
      ELSE IF(MNDG.EQ.2) THEN                                           
        DO 50 I=1,4                                                     
        IACHT(2*I-1)=2*I                                                
        IACHT(2*I)=2*I-1                                                
50      CONTINUE                                                        
      ELSE IF(MNDG.EQ.3) THEN                                           
        DO 60 I=0,1                                                     
        DO 60 J=1,2                                                     
        IACHT(4*I+J)=4*I+J+2                                            
        IACHT(4*I+J+2)=4*I+J                                            
60      CONTINUE                                                        
      ELSE IF(MNDG.EQ.4) THEN                                           
        DO 70 I=1,4                                                     
        IACHT(I)=I+4                                                    
        IACHT(I+4)=I                                                    
70      CONTINUE                                                        
      END IF                                                            
!      **  CREATION OF TETRAHEDRA  **                                   
!      **  (1248);(1438);(1378);(1758);(1568);(1628)                    
       DO 80 I=1,6                                                      
       TET(1,I)=IACHT(1)                                                
       TET(4,I)=IACHT(8)                                                
80     CONTINUE                                                         
       TET(2,1)=IACHT(2)                                                
       TET(3,1)=IACHT(4)                                                
       TET(2,2)=IACHT(4)                                                
       TET(3,2)=IACHT(3)                                                
       TET(2,3)=IACHT(3)                                                
       TET(3,3)=IACHT(7)                                                
       TET(2,4)=IACHT(7)                                                
       TET(3,4)=IACHT(5)                                                
       TET(2,5)=IACHT(5)                                                
       TET(3,5)=IACHT(6)                                                
       TET(2,6)=IACHT(6)                                                
       TET(3,6)=IACHT(2)                                                
!*                                                                      
      DO 90 I=1,4                                                       
      DO 90 J=1,6                                                       
      TET0(1,I,J)=(TET(I,J)-1)/4                                        
      TET0(2,I,J)=(TET(I,J)-TET0(1,I,J)*4-1)/2                          
      TET0(3,I,J)=TET(I,J)-TET0(1,I,J)*4-TET0(2,I,J)*2-1                
90    CONTINUE                                                          
      RETURN                                                            
      END                                                               
