!     ..................................................................
      SUBROUTINE BASDIV(N,NMSHP,GBAS,IARB)                              
!     **                                                              **
!     **  BASDIV DETERMINES DIVISION OF BASEVECTORS OF REC. LATT.     **
!     **  SO THAT THE NUMBER OF MESHPOINTS IS JUST BELOW NMSHP AND    **
!     **  TAKES INTO ACCOUNT THE DEPENDENCY BY POINTSYMMETRY          **
!     **  INPUT :                                                     **
!     **    GBAS        RECIPROCAL LATTICE VECTORS                    **
!     **    IARB        DEPENDENCIES FOR DIVISIONS                    **
!     **                OF RECIPROCAL LATTICE VECTORS                 **
!     **                ( IF IARB(1)=1 THEN 1ST AND 2ND LATTICE       **
!     **                  VECTORS ARE DIVIDED BY AN EQUAL NUMBER;     **
!     **                  IF IARB(2)=1 THEN SAME FOR 2ND AND 3RD;     **
!     **                  IF IARB(3)=1 THEN SAME FOR 3RD AND 1ST)     **
!                      2 and 3 interchanged ??
!     **    NMSHP       TOTAL NUMBER OF GRID POINTS                   **
!     **  OUTPUT:                                                     **
!     **    N           NUMBER OF DIVISIONS OF RECIPROCAL LATTICE     **
!     **                VECTORS FOR DEFINITION OF SUBLATTICE          **
!     **    NMSHP       NUMBER OF SUBLATTICE POINTS IN A REC. UNIT CELL*
!     **                AND ON ALL ITS FACES                          **
!     **                                                              **
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)                               
      DIMENSION GBAS(3,3),IARB(3)                                       
      DIMENSION BETR(3),RN(3),N(3)                                      
      DATA OPAR/1.D-6/                                                  
!     ==================================================================
!     == FIND UPPER LIMIT FOR THE LENGTH OF SUBLATTICE VECTORS        ==
!     ==================================================================
      DO 30 I=1,3                                                       
      BETR(I)=DSQRT(GBAS(I,1)**2+GBAS(I,2)**2+GBAS(I,3)**2)            
!c      BETR(I)=DSQRT(GBAS(1,I)**2+GBAS(2,I)**2+GBAS(3,I)**2)           
30    CONTINUE                                                          
      SVAR=(DBLE(NMSHP)/(BETR(1)*BETR(2)*BETR(3)))**(1.D0/3.D0)         
      DO 40 I=1,3                                                       
      RN(I)=BETR(I)*SVAR                                                
40    CONTINUE                                                          
!     ==================================================================
!     == FIND DIVISIONS OF LATTICE VECTORS                            ==
!     ==================================================================
      IF(IARB(1).EQ.1.AND.IARB(2).EQ.1)THEN                             
        N(1)=INT((RN(1)*RN(2)*RN(3))**(1.D0/3.D0)+OPAR)                 
        N(2)=N(1)                                                       
        N(3)=N(1)                                                       
      ELSE IF(IARB(1).EQ.1) THEN                                        
        N(1)=INT(DSQRT(RN(1)*RN(2))+OPAR)                               
        N(2)=N(1)                                                       
        N(3)=INT(RN(3))                                                 
      ELSE IF(IARB(2).EQ.1) THEN                                        
        N(1)=INT(DSQRT(RN(1)*RN(3))+OPAR)                               
        N(2)=INT(RN(2))                                                 
        N(3)=N(1)                                                       
      ELSE IF (IARB(3).EQ.1) THEN                                       
        N(1)=INT(RN(1))                                                 
        N(2)=INT(DSQRT(RN(2)*RN(3))+OPAR)                               
        N(3)=N(2)                                                       
      ELSE                                                              
        N(1)=INT(RN(1)+OPAR)                                            
        N(2)=INT(RN(2)+OPAR)                                            
        N(3)=INT(RN(3)+OPAR)                                            
      END IF                                                            
      N(1)=MAX0(1,N(1))                                                 
      N(2)=MAX0(1,N(2))                                                 
      N(3)=MAX0(1,N(3))                                                 
      write(66,*) ' length of reciprocal lattic vectors:',rn
      write(*,'(" length of reciprocal lattic vectors:",6f8.3)') BETR,rn
!     ==================================================================
!     == USE THIS TO FIX THE K-MESH PY HAND                           ==
!     ==================================================================
!     PRINT*,' K-MESH PUT IN BY HAND!!!!!!!!!!'                         
!     N(1)=12                                                           
!     N(2)=12                                                           
!     N(3)=2                                                            
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      if(NMSHP.le.0) then
        write(*,*) " Specify 3 mesh-divisions (n1,n2,n3):"
        read(*,*) n
      endif
      NMSHP=(N(1)+1)*(N(2)+1)*(N(3)+1)                                  
      RETURN                                                            
      END                                                               
