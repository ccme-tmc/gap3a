!     ..................................................................
      SUBROUTINE ZUORD(NMSHP,NUM,N,ISHIFT,GBAS,IDKP,NKP,BK,bki)         
!     **                                                              **
!     **  ZUORD CREATES A RELATION BETWEEN "IRREDUCIBLE POINTS" AND   **
!     **  THEIR PLACES ON THE FILE, COUNTS THEM, AND CALCULATES THEIR **
!     **  COORDINATES IN K-SPACE                                      **
!     **                                                              **
!     **  INPUT:                                                      **
!     **    NMSHP       NUMBER OF SUBLATTICE POINTS INSIDE AND        **
!     **                ON ALL FACES OF A REC. UNIT CELL              **
!     **    NUM(I)      MAPPING FROM A GENERAL POINT (I) TO THE       **
!     **                CORRESPONDING IRREDUCIBLE POINT (NUM)         **
!     **    N           NUMBER OF DIVISIONS OF REC. LATTICE VECTORS   **
!     **    ISHIFT                                                    **
!     **    IDKP        MAXIMUM NUMBER OF IRREDUCIBLE K-POINTS        **
!     **    GBAS        RECIPROCAL LATTICE VECTORS                    **
!     **  OUTPUT :                                                    **
!     **    NKP         ACTUAL NUMBER OF IRREDUCIBLE K-POINTS         **
!     **    NUM(I)      MAPPING FROM A GENERAL POINT (I) TO THE       **
!     **                CORRESPONDING IRREDUCIBLE POINT (NUM)         **
!     **                (SEE REMARKS)                                 **
!     **    BK          IRREDUCIBLE K-POINT IN CARTESIAN COORDINATES  **
!     **                                                              **
!     **  REMARKS :                                                   **
!     **    NUM(I) REFERS TO THE ACTUAL POSITION OF THE IRREDUCIBLE   **
!     **    K-POINT ON INPUT. ON OUTPUT IT ONLY REFERS TO THE         **
!     **    SECOND INDEX OF BK(I,IKP)                                 **
!     **                                                              **
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)                               
      common /inp/ iswitch
      DIMENSION BK(3,IDKP),NUM(NMSHP),N(3),ISHIFT(3),GBAS(3,3)          
      dimension bki(3,idkp)
      NDIM=0                                                            
      I=0                                                               
      if(iswitch.eq.1)write(66,*)' internal and cartesian k-vectors:'
      DO 10 I1=1,N(1)+1                                                 
      DO 10 I2=1,N(2)+1                                                 
      DO 10 I3=1,N(3)+1                                                 
      I=I+1                                                             
      IF(I.GT.NMSHP) THEN                                               
        PRINT*,'I.GT.NMSHP,STOP',I,I1,I2,I3                             
        STOP                                                            
      END IF                                                            
      IF(I.EQ.NUM(I))THEN                                               
        NDIM=NDIM+1                                                     
        IF(NDIM.GT.IDKP) THEN                                           
          WRITE(66,1000)                                                
1000      FORMAT(1H ,'NUMBER OF INEQUIVALENT POINTS EXCEEDS IDKP')      
          STOP  'idkp too small in zuord'                     
        END IF                                                          
        NUM(I)=NDIM                                                     
        RINDA=(DBLE(I1-1)+DBLE(ISHIFT(1))/2.D0)/DBLE(N(1))              
        RINDB=(DBLE(I2-1)+DBLE(ISHIFT(2))/2.D0)/DBLE(N(2))              
        RINDC=(DBLE(I3-1)+DBLE(ISHIFT(3))/2.D0)/DBLE(N(3))              
        BK(1,NDIM)=GBAS(1,1)*RINDA+GBAS(2,1)*RINDB+GBAS(3,1)*RINDC      
        BK(2,NDIM)=GBAS(1,2)*RINDA+GBAS(2,2)*RINDB+GBAS(3,2)*RINDC      
        BK(3,NDIM)=GBAS(1,3)*RINDA+GBAS(2,3)*RINDB+GBAS(3,3)*RINDC      
!       BK(1,NDIM)=GBAS(1,1)*RINDA+GBAS(1,2)*RINDB+GBAS(1,3)*RINDC      
!       BK(2,NDIM)=GBAS(2,1)*RINDA+GBAS(2,2)*RINDB+GBAS(2,3)*RINDC      
!       BK(3,NDIM)=GBAS(3,1)*RINDA+GBAS(3,2)*RINDB+GBAS(3,3)*RINDC      
        bki(1,ndim)=rinda
        bki(2,ndim)=rindb
        bki(3,ndim)=rindc
        if(iswitch.eq.1)  &
        write(66,100) rinda,rindb,rindc,bk(1,ndim),bk(2,ndim),bk(3,ndim)
 100     format(3f10.5,10x,3f10.5)
      ELSE                                                              
        IF(NUM(I).GT.NMSHP) PRINT*,'ERROR'                              
        NUM(I)=NUM(NUM(I))                                              
      END IF                                                            
10    CONTINUE                                                          
      NKP=NDIM                                                          
      RETURN                                                            
      END                                                               
