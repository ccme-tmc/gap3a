!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     ----                                                         -----
!     ----  BLOCK ARBMSH:                                          -----
!     ----                                                         -----
!     ----  READS SYMMETRYELEMENTS AND FINDS IRR. K-POINTS         -----
!     ----  AND TETRAHEDRA, USED IN THE BRILLOUIN ZONE INTEGRATION -----
!     ----                                                         -----
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     .....................................................ARBMSH.......
      SUBROUTINE ARBMSH(RBAS,gbas,IDKP,NKP,BK,NSYM,IIO,IARB,W,NWX,N, &
       weight,wsum,sumwgt,bki)    
!     **                                                              **
!     **  CALCULATE IRREDUCIBLE K-POINTS                              **
!     **  AND FINDS INEQUIVALENT TETRAHEDRA                           **
!     **  FOR BRILLOUIN ZONE INTEGRATION                              **
!     **                                                              **
!     **  INPUT :                                                     **
!     **    RBAS        LATTICE VECTORS                               **
!     **    IDKP        MAXIMUM NUMBER OF IRREDUCIBLE K-POINTS        **
!     **                ( DIMENSION)                                  **
!     **    NKP         UPPER LIMIT FOR TOTAL NUMBER OF K-POINTS      **
!     **                ( PROGRAM SEARCHES THE LARGEST "POSSIBLE"     **
!     **                  NUMBER BELOW THIS FOR K-POINT MESH )        **
!     **    NSYM        NUMBER OF SYMMETRY OPERATIONS                 **
!     **    IIO         SYMMETRY OPERATIONS                           **
!     **    IARB        DEPENDENCIES FOR DIVISIONS                    **
!     **                OF RECIPROCAL LATTICE VECTORS                 **
!     **                ( IF IARB(1)=1 THEN 1ST AND 2ND LATTICE       **
!     **                  VECTORS ARE DIVIDED BY AN EQUAL NUMBER;     **
!     **                  IF IARB(3)=1 THEN SAME FOR 2ND AND 3RD;     **
!     **                  IF IARB(2)=1 THEN SAME FOR 3RD AND 1ST)     **
!     **    W           (INTEGER) WORK ARRAY OF LENGTH NWX            **
!     **    NWX         LENGTH OF WORK ARRAY                          **
!     **                                                              **
!     **  OUTPUT :                                                    **
!     **    NKP         NUMBER OF IRREDUCIBLE K-POINTS                **
!     **    BK          IRREDUCIBLE K-POINTS                          **
!     **                                                              **
!     **  OUTPUT ON FILE 15 :                                         **
!     **   NDIM         NUMBER OF INEQUIVALENT K-POINTS               **
!     **   ITTFL        STORES FOR EACH DIFFERENT TETRAHEDRON FIRST   **
!     **                THE NUMBER OF APPEARENCE OF THIS TETRAHEDRON  **
!     **                AND THEN FOUR NUMBERES OF EDGEPOINTS          **
!     **   V            =1/(NUMBER OF TETRAHEDRA)                     **
!     **   NTT          NUMBER OF DIFFERENT TETRAHEDRA                **
!     **                                                              **
!     **   AUTHOR: PETER E. BLOECHL                                   **
!     **                                                              **
!     **   SUBROUTINES USED:                                          **
!     **     DEF0,GBASS,BASDIV,REDUZ,ZUORD,TETDIV,TETCNT,ORD1         **
!     **     (ORDERS - CRAY ROUTINE,OPTIONAL,NOT IN THIS PACKAGE)     **
!     **                                                              **
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)                               
      IMPLICIT INTEGER (O)                                              
      INTEGER W(NWX)                                                    
      DIMENSION RBAS(3,3),GBAS(3,3)                                     
      DIMENSION N(3),IARB(3),IIO(3,3,NSYM)                              
      DIMENSION BK(3,IDKP),weight(idkp),wsum(nkp),bki(3,idkp)
      DATA IPR/0/IWNUM/0/IWBK/1/                                        
      CALL DEF0(NWX)                                                    
      INV=0                                                             
!     ------------------------------------------------------------------
!     -- DEFINE MESH                                                  --
!     ------------------------------------------------------------------
      NMSHP=NKP                                                         
      write(66,*) '   G1        G2        G3'
      DO 281 J=1,3                                                      
  281 WRITE (66,290) (GBAS(i,j),I=1,3)                                 
      CALL GBASS(RBAS,GBAS)                                             
      write(66,*) '   G1        G2        G3'
      DO 280 J=1,3                                                      
  280 WRITE (66,290) (GBAS(i,j),I=1,3)                                 
  290 FORMAT (1H ,3F10.6)                                 
      CALL BASDIV(N,NMSHP,GBAS,IARB)                                    
!     ------------------------------------------------------------------
!     -- FIND IRREDUCIBLE K-POINTS                                    --
!     ------------------------------------------------------------------
      CALL DEFI(ONUM,NMSHP)                                             
      CALL DEFI(OISHIF,3)                                               
      CALL DEFI(OIN,3*NMSHP)                                            
      CALL REDUZ(N,NMSHP,W(OISHIF),NSYM,IIO,W(ONUM),W(OIN),idkp,weight, &
       wsum,sumwgt)             
      CALL RLSE(OIN)                                                    
      IF(W(OISHIF).NE.0.OR.W(OISHIF+1).NE.0.OR.W(OISHIF+2).NE.0)INV=0   
      CALL ZUORD(NMSHP,W(ONUM),N,W(OISHIF),GBAS,IDKP,NKP,BK,bki)        
      CALL RLSE(OISHIF)                                                 
!     ------------------------------------------------------------------
!     -- PRINTOUT OF SYMMETRY OPERATIONS                              --
!     ------------------------------------------------------------------
      IF(IPR.EQ.1) THEN                                                 
        DO 1111 ISYM=1,NSYM                                             
        WRITE(66,6003)ISYM,((IIO(I,J,ISYM),J=1,3),I=1,3)                
6003    FORMAT(1H ,'SYMMETRYMATRIX NR. : ',I5/3(1H ,3I10/))             
1111    CONTINUE                                                        
      END IF                                                            
!     ----------------------------------------------------------------- 
!     --  PRINTOUT OF MAPPING TO IRREDUCIBLE K-POINTS                -- 
!     ----------------------------------------------------------------- 
      IF(IWNUM.EQ.1.AND.N(3)+1.LE.25) THEN                              
        DO 90 I=0,N(1)                                                  
        DO 80 J=N(2),0,-1                                               
        INDEX1=I*(N(2)+1)*(N(3)+1)+J*(N(3)+1)+1                         
        INDEX2=INDEX1+N(3)                                              
        WRITE(66,1013)(W(ONUM-1+K),K=INDEX1,INDEX2)                     
1013    FORMAT(25I5)                                                    
80      CONTINUE                                                        
        WRITE(66,1014)                                                  
1014    FORMAT(1H0)                                                     
90      CONTINUE                                                        
      END IF                                                            
!     ----------------------------------------------------------------- 
!     --  PRINTOUT OF IREDUCIBLE K-POINTS                            -- 
!     ----------------------------------------------------------------- 
      WRITE(66,1012)NKP                                                 
1012  FORMAT(1H ,' NO. OF INEQUIVALENT K-POINTS ',I5)                   
      IF(IWBK.NE.0) THEN                                                
        WRITE(66,6000)                                                  
6000    FORMAT(1H ,' INEQUIVALENT BLOCH VECTORS')                       
        NDIMM=NKP/2+1                                                   
        DO 20 I=1,NDIMM                                                 
        NDI1=(I-1)*2+1                                                  
        NDI2=NDI1+1                                                     
        IF(I.EQ.NDIMM)NDI2=NKP                                          
        WRITE(66,6010)(K,(BK(J,K),J=1,3),K=NDI1,NDI2)                   
20      CONTINUE                                                        
6010    FORMAT(1H ,I4,'(',3F10.6,')',I4,' (',3F10.6,')')                
      END IF                                                            
!     ----------------------------------------------------------------- 
!     --  CHOOSE TETRAHEDRA                                          -- 
!     ----------------------------------------------------------------- 
      CALL DEFI(OTET0,72)                                               
      CALL TETDIV(N,GBAS,W(OTET0))                                      
      IF(IPR.EQ.1) THEN                                                 
        WRITE(66,6001)1,2,3                                             
        WRITE(66,6002)                                                   &
               (((W(OTET0-1+J+(K-1)*3+(I-1)*12),K=1,4),I=1,3),J=1,3)    
        WRITE(66,6001)4,5,6                                             
        WRITE(66,6002)                                                   &
               (((W(OTET0-1+J+(K-1)*3+(I-1)*12),K=1,4),I=4,6),J=1,3)    
6001    FORMAT(1H ,3(I2,'-TES NORMTETRAHEDRON  '))                      
6002    FORMAT(1H ,4I5,5X,4I5,5X,4I5)                                   
      END IF                                                            
!     ----------------------------------------------------------------- 
!     --  FIND INEQUIVALENT TETRAHEDRA                               -- 
!     ----------------------------------------------------------------- 
! CAD  Mwrit von 100 auf 101, fixes lapw2-TETRA problems
      MWRIT=101                                                         
      CALL DEFI(OIY,6*4*NMSHP)                                          
      CALL DEFI(OITTFL,5*MWRIT)                                         
      N3=N(3)                                                           
      CALL DEFI(OITET,4*N3*6)                                           
      CALL TETCNT(NMSHP,W(ONUM),W(OTET0),N,INV,NKP,MWRIT,W(OITTFL)       &
                                   ,W(OIY),N3,W(OITET),W)               
      CALL RLSE(OTET0)                                                  
      CALL RLSE(ONUM)                                                   
      RETURN                                                            
      END                                                               
