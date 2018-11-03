!     .....................................................TETCNT.......
      SUBROUTINE TETCNT(NMSHP,NUM,TET0,N,INV,NKP,MWRIT,ITTFL,IY          &
                                                      ,N3,ITET,W)       
!     **                                                              **
!     **  TETCNT CALCULATES ALL DIFFERENT TETRAHEDRA AND COUNTS THEM  **
!     **  INPUT :                                                     **
!     **    NMSHP       NUMBER OF SUBLATTICE POINTS INSIDE AND        **
!     **                ON ALL FACES OF A REC. UNIT CELL              **
!     **    NUM(I)      MAPPING FROM A GENERAL POINT (I) TO THE       **
!     **                CORRESPONDING IRREDUCIBLE POINT (NUM)         **
!     **    TET0(I,J,K) COORDINATES (I) IN THE BASIS OF SUBLATTICE   ** 
!     **                VECTORS, OF THE 4 CORNERS (J) OF A           ** 
!     **                TETRAHEDRON (K) IN A SUBLATTICE UNIT CELL    ** 
!     **    N           NUMBER OF DIVISIONS OF REC. LATTICE VECTORS   **
!     **    INV         DUMMY NUMBER (MUST BE 0)                      **
!     **    NKP         NUMBER OF IRREDUCIBLE K-POINTS                **
!     **    MWRIT       INFORMATION FOR MWRIT TETRAHEDRA ARE WRITTEN  **
!     **                AT ONE TIME.                                  **
!     **    W           INTEGER WORK ARRAY (INITIALIZED IN ARBMSH)    **
!     **                                                              **
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)                               
      IMPLICIT INTEGER (O)                                              
      DIMENSION N(3)                                                    
      INTEGER NUM(NMSHP),TET0(3,4,6)                                     &
           ,ITET(4,N3*6),IY(4,6*NMSHP)                                  
      DIMENSION ITTFL(5*MWRIT)                                          
      INTEGER W(1)  
      DATA ICHK/1/IPR/0/                                                
!      IBGGST=2*(2**30-1)+1                                             
!      IF(NKP**4.GT.IBGGST) THEN                                        
!        PRINT*,' TEST FOR BIGGEST INTEGER FAILED'                      
!        PRINT*,' CHECK WHETHER IBGGST IS THE BIGGEST INTEGER      '    
!        PRINT*,' IBGGST=',IBGGST                                       
!        PRINT*,' IF NOT, CORRECT IT IN ROUTINE TETCNT             '    
!        PRINT*,' ELSE REDUCE NUMBER OF K-POINTS                   '    
!        PRINT*,' OR CHANGE TETCNT                                 '    
!        STOP                                                           
!      END IF                                                           
      NTMAX=N(1)*N(2)*N(3)*6                                            
      NTMAX=NTMAX/(1+INV)                                               
      IPP=0                                                             
      DO 90 K1=1,N(1)                                                   
      DO 90 K2=1,N(2)                                                   
      IND=0                                                             
      DO 40 K3=1,N(3)                                                   
      IP=K3+(N(3)+1)*((K2-1)+(N(2)+1)*(K1-1))                           
                                                                        
      DO 40 I=1,6                                                       
      IND=IND+1                                                         
      DO 30 J=1,4                                                       
      IXX=TET0(1,J,I)*(N(2)+1)*(N(3)+1)                                  &
         +TET0(2,J,I)*(N(3)+1)                                           &
         +TET0(3,J,I)                                                   
      ITET(J,IND)=IP+IXX                                                
30    CONTINUE                                                          
40    CONTINUE                                                          
                                                                        
!     --  TRANSFORM THE EDGEPOINTS ONTO THE IRREDUCIBLE POINTS          
      DO 50 M=1,4                                                       
      DO 50 J=1,N(3)*6                                                  
      ITET(M,J)=NUM(ITET(M,J))                                          
50    CONTINUE                                                          
                                                                        
!     --  ORDER THE POINTS OF EACH TETR. ACC. TO INREASING NUMBER       
      DO 61 K=1,3                                                       
      DO 61 J=K+1,4                                                     
      DO 61 L=1,N(3)*6                                                  
      ISVAR1=ITET(K,L)                                                  
      ISVAR2=ITET(J,L)                                                  
      ITET(K,L)=MIN0(ISVAR1,ISVAR2)                                     
      ITET(J,L)=MAX0(ISVAR1,ISVAR2)                                     
61    CONTINUE                                                          
                                                                        
!     --  IDENTIFY THE TETRAHEDRA WITH INTEGERS                         
      DO 80 M=1,N(3)*6                                                  
      IPP=IPP+1                                                         
!      IY(IPP)=ITET(1,M)*(NKP+1)**3                                     
!     1       +ITET(2,M)*(NKP+1)**2                                     
!     2       +ITET(3,M)*(NKP+1)                                        
!     3       +ITET(4,M)                                                
      IY(1,IPP)=ITET(1,M)
      IY(2,IPP)=ITET(2,M)
      IY(3,IPP)=ITET(3,M)
      IY(4,IPP)=ITET(4,M)

80    CONTINUE                                                          
      IF(IPP.GE.NTMAX) THEN                                             
        IPP=NTMAX                                                       
        GOTO 100                                                        
      END IF                                                            
!*                                                                      
90    CONTINUE                                                          
      PRINT*,'UNNORMAL END OF LOOP.......................STOP IN TETCNT'
      STOP                                                              
                                                                        
100   CONTINUE                                                          
!     ------------------------------------------------------------------
!     --  ORDER TETRAHEDRA                                            --
!     ------------------------------------------------------------------
      NTET=IPP                                                          
      write(66,*)'tetrahedra to sort: (min. size of parameter INDEXM in ord1.f', ntet
      CALL DEFI(OWORK,NTET)                                             
      write(66,*)'owork: (min. size of parameter NWX in main.f)',owork  
      CALL ORD1(NTET,IY,W(OWORK))                                       
!
      CALL RLSE(OWORK)                                                  
!     ==  CHECK ORDERING AND CALCULATE NUMBER OF INEQUIVALENT TETRAHEDRA
      NTT=1                                                             
      DO 111 I=1,NTET-1                                                 
!     IF(IY(I+1).GT.IY(I)) THEN                                         
!       NTT=NTT+1                                                       
!     END IF                                                            
!     IF(IY(I+1).LT.IY(I)) THEN                                         
!       PRINT*,' I ',I,' IY ',IY(I)                                     
!     END IF                                                            
      NTT=NTT+1                                                         
      IF(IY(4,I+1).eq.IY(4,I)) THEN                                     
      IF(IY(3,I+1).eq.IY(3,I)) THEN                                     
      IF(IY(2,I+1).eq.IY(2,I)) THEN                                     
      IF(IY(1,I+1).eq.IY(1,I)) THEN                                     
      ntt=ntt-1
      end if
      end if
      end if
      end if
!
111   CONTINUE                                                          
      WRITE(66,1018)NTT                                                 
1018  FORMAT(1H ,'NUMBER OF DIFFERENT TETRAHEDRA :',I5)                 
!     ------------------------------------------------------------------
!     --  WRITE ON FILE                                               --
!     ------------------------------------------------------------------
      if((ntt/mwrit)*mwrit.eq.ntt) then
      NREC=NTT/MWRIT  
      else
      NREC=NTT/MWRIT+1
      end if                                                            
      V=1.D0/DBLE(6*N(1)*N(2)*N(3))                                     
      SUM=V*DBLE(NTET*(1+INV))                                          
      IF(DABS(SUM-1.D0).GT.1.D-5) THEN                                  
        PRINT*,'SUMRULE NOT FULLFILLED...................STOP IN TETCNT'
        PRINT*,' SUM ',SUM,' SHOUD BE EQUAL TO 1'                       
        STOP                                                            
      END IF                                                            
      REWIND 15                                                         
      WRITE(15,1234)NKP,NTT,V,MWRIT,NREC                                
      NTT=0                                                             
      NREC=0                                                            
      CALL INITI(ITTFL,5*MWRIT)                                         
      DO 140 I=1,NTET                                                   
      IF(I.NE.1.AND.IY(1,I).EQ.IY(1,I-1).AND.IY(2,I).EQ.IY(2,I-1) &
       .AND.IY(3,I).EQ.IY(3,I-1).AND.IY(4,I).EQ.IY(4,I-1)) THEN         
        NI=5*(NTT-1-NREC*MWRIT)+1                                       
        ITTFL(NI)=ITTFL(NI)+1+INV                                       
      ELSE                                                              
        NTT=NTT+1                                                       
                                                                        
        IF(NTT.GT.MWRIT*(NREC+1)) THEN                                  
           NREC=NREC+1                                                  
           WRITE(15,1235)(ITTFL(J),J=1,5*MWRIT)                         
           CALL INITI(ITTFL,5*MWRIT)                                    
        END IF                                                          
                                                                        
        NI=5*(NTT-1-NREC*MWRIT)                                         
        ITTFL(NI+1)=1+INV                                               
!        ITTFL(NI+2)=IY(I)/(NKP+1)**3                                   
!        ITTFL(NI+3)=(IY(I)-ITTFL(NI+2)*(NKP+1)**3)/(NKP+1)**2          
!        ITTFL(NI+4)=                                                   
!     1            (IY(I)-ITTFL(NI+2)*(NKP+1)**3-ITTFL(NI+3)*(NKP+1)**2)
!     1            /(NKP+1)                                             
!        ITTFL(NI+5)=IY(I)-ITTFL(NI+2)*(NKP+1)**3-ITTFL(NI+3)*(NKP+1)**2
!     1                                  -ITTFL(NI+4)*(NKP+1)           
        ITTFL(NI+2)=IY(1,I)
        ITTFL(NI+3)=IY(2,I)
        ITTFL(NI+4)=IY(3,I)
        ITTFL(NI+5)=IY(4,I)
      END IF                                                            
140   CONTINUE                                                          
      NREC=NREC+1                                                       
      WRITE(15,1235)ITTFL                                               
                                                                        
      IF(ICHK.EQ.0) RETURN                                              
!     ------------------------------------------------------------------
!     --  CHECK OF I/O                                                --
!     ------------------------------------------------------------------
      REWIND 15                                                         
      READ(15,1234)NKP,NTT,V,MWRIT,NREC                                 
 1234 format(2i10,e20.12,2i10) 
      IF(IPR.EQ.1) THEN                                                 
        WRITE(66,6000)NKP,NTT,V,MWRIT,NREC                              
6000    FORMAT(1H ,' NKP= ',I5,' NTT= ',I5,' V= ',F7.5,' MWRIT= ',I5     &
                  ,' NREC= ',I5)                                        
      END IF                                                            
      SUM=0.D0                                                          
      DO 200 I=1,NREC                                                   
      CALL INITI(ITTFL,5*MWRIT)                                         
      READ(15,1235)(ITTFL(J),J=1,5*MWRIT)                               
 1235 format(6i10)
      JMAX=5*MIN0(MWRIT,NTT-(I-1)*MWRIT)                                
      IF(IPR.EQ.1) THEN                                                 
        PRINT*,'REC = ',I,'--------------------------------------'      
        WRITE(66,6001)(ITTFL(J),J=1,JMAX)                               
6001    FORMAT(1H ,5I5,5X,5I5)                                          
      END IF                                                            
      DO 200 J=1,JMAX/5                                                 
      SUM=SUM+DBLE(ITTFL(5*(J-1)+1))*V                                  
200   CONTINUE                                                          
      IF(DABS(SUM-1.D0).GT.1.D-5) THEN                                  
        PRINT*,'SUMRULE NOT FULLFILLED...................STOP IN TETCNT'
        PRINT*,' SUM ',SUM,' SHOUD BE EQUAL TO 1'                       
        STOP                                                            
      END IF                                                            
      RETURN                                                            
      END                                                               
