!     ..................................................................
      SUBROUTINE REDUZ(N,NMSHP,ISHIFT,NSYM,IO,NUM,IN,idkp,weight,wsum, &
       sumwgt)                   
!     **                                                              **
!     **  REDUZ CREATES THE RELATION BETWEEN THE MESHPOINTS AND       **
!     **  THE POINTS IN THE "IRREDUCIBLE ZONE"                        **
!     **  INPUT :                                                     **
!     **    N           NUMBER OF DIVISIONS OF REC. LATTICE VECTORS   **
!     **    NMSHP       NUMBER OF SUBLATTICE POINTS INSIDE AND        **
!     **                ON ALL FACES OF A REC. UNIT CELL              **
!     **    NSYM        NUMBER OF SYMMETRY OPERATIONS                 **
!     **    IO          SYMMETRY OPERATIONS (BASIS ARE REC. LATT. VEC.)*
!     **    IN          WORK ARRAY                                    **
!     **  OUTPUT :                                                    **
!     **    NUM(I)      MAPPING FROM A GENERAL POINT (I) TO THE       **
!     **                CORRESPONDING IRREDUCIBLE POINT (NUM)         **
!     **  REMARKS :                                                   **
!     **    THE MAPPING FROM COORDINATES TO NUMBERS IS GIVEN BY :     **
!     **   (X,Y,Z)=RBAS*(I,J,K)                                       **
!     **   (I,J,K)  <->  NUM = I*(N(2)+1)*(N(3)+1)+J*(N(3)+1)+K+1     **
!     **                                                              **
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)                               
      common /inp/ iswitch
      DIMENSION N(3),IO(3,3,NSYM),ISHIFT(3),NUM(NMSHP),IN(3,NMSHP)
      dimension wsum(nmshp),weight(idkp)
!     ------------------------------------------------------------------
!     -- TEST WHETHER SHIFT OF THE SUBLATTICE BY (1/2,1/2,1/2) IS     --
!     --  COMPATIBLE WITH THE SYMMETRYGROUP                           --
!     ------------------------------------------------------------------
      ISHIFT(1)=1                                                       
      ISHIFT(2)=1                                                       
      ISHIFT(3)=1                                                       
      DO 10 I=1,8                                                       
      IN(1,I)=(I-1)/4                                                   
      IN(2,I)=(I-IN(1,I)*4-1)/2                                         
      IN(3,I)=I-IN(1,I)*4-IN(2,I)*2-1                                   
10    CONTINUE                                                          
      DO 20 I=1,NSYM                                                    
      DO 20 J=1,8                                                       
      I1=2*IN(1,J)+ISHIFT(1)                                            
      I2=2*IN(2,J)+ISHIFT(2)                                            
      I3=2*IN(3,J)+ISHIFT(3)                                            
      J1=IO(1,1,I)*I1+IO(1,2,I)*I2+IO(1,3,I)*I3                         
      J2=IO(2,1,I)*I1+IO(2,2,I)*I2+IO(2,3,I)*I3                         
      J3=IO(3,1,I)*I1+IO(3,2,I)*I2+IO(3,3,I)*I3                         
      IF(DMOD(DBLE(J1-ISHIFT(1)),2.D0).NE.0.D0.OR.                       &
         DMOD(DBLE(J2-ISHIFT(2)),2.D0).NE.0.D0.OR.                       &
         DMOD(DBLE(J3-ISHIFT(3)),2.D0).NE.0.D0) THEN                    
        ISHIFT(1)=0                                                     
        ISHIFT(2)=0                                                     
        ISHIFT(3)=0                                                     
        ISYM=I                                                          
        write(66,*) 'SUBMESH SHIFT CONFLICTS WITH POINTGROUP'           
        WRITE(66,6003)ISYM,((IO(K,L,ISYM),K=1,3),L=1,3)                 
6003    FORMAT(1H ,'SYMMETRYMATRIX NR. : ',I5/3(1H ,3I10/))             
        GOTO 30                                                         
      END IF                                                            
20    CONTINUE                                                          
30    CONTINUE                                                          
      IF(ISHIFT(1).EQ.1.OR.ISHIFT(2).EQ.1.OR.ISHIFT(3).EQ.1) THEN       
        write(*,*) ' Shift of k-mesh allowed. Do you want to shift:', &
        ' (0=no, 1=shift)'
      read(*,*) ishift(1)
      ishift(2)=ishift(1)
      ishift(3)=ishift(1)
        write(66,*) ' SUBMESH SHIFTED; SHIFT: ',ISHIFT                  
      ELSE                                                              
        write(66,*)' SUBMESH NOT SHIFTED; SHIFT: ',ISHIFT               
      END IF        
!     ==================================================================
!     ==  INITIALIZE                                                  ==
!     ==================================================================
      DO 40 I=1,NMSHP
      wsum(i)=0.0                                                       
      NUM(I)=I                                                          
      IN(1,I)=(I-1)/((N(3)+1)*(N(2)+1))                                 
      IN(2,I)=(I-IN(1,I)*(N(2)+1)*(N(3)+1)-1)/(N(3)+1)                  
      IN(3,I)=I-IN(1,I)*(N(2)+1)*(N(3)+1)-IN(2,I)*(N(3)+1)-1            
40    CONTINUE                                                          
!     ==================================================================
!     ==  REDUCTION OF REDUCIBLE K-POINTS                             ==
!     ==================================================================
      DO 50 I=1,NSYM                                                    
      DO 50 J=1,NMSHP                                                   
      I1=2*IN(1,J)+ISHIFT(1)                                            
      I2=2*IN(2,J)+ISHIFT(2)                                            
      I3=2*IN(3,J)+ISHIFT(3)                                            
      J1=MOD(IO(1,1,I)*I1+IO(1,2,I)*I2+IO(1,3,I)*I3,2*N(1))             
      J2=MOD(IO(2,1,I)*I1+IO(2,2,I)*I2+IO(2,3,I)*I3,2*N(2))             
      J3=MOD(IO(3,1,I)*I1+IO(3,2,I)*I2+IO(3,3,I)*I3,2*N(3))             
      J1=J1+(1-ISIGN(1,J1))*N(1)                                        
      J2=J2+(1-ISIGN(1,J2))*N(2)                                        
      J3=J3+(1-ISIGN(1,J3))*N(3)                                        
      J1=(J1-ISHIFT(1))/2                                               
      J2=(J2-ISHIFT(2))/2                                               
      J3=(J3-ISHIFT(3))/2                                               
      NUM(J)=MIN0(NUM(J),J1*(N(2)+1)*(N(3)+1)                            &
                                  +J2*(N(3)+1)+J3+1) 
50    CONTINUE  
!
      WRITE(66,6000)N(1)*N(2)*N(3),(N(I),I=1,3)                         
 6000 FORMAT(1H ,' NO. OF MESH POINTS IN THE BRILLOUIN ZONE =',I6/       &
      '  DIVISION OF RECIPROCAL LATTICE VECTORS (INTERVALS)=',3I5)      
       if(iswitch.eq.1)write(66,*) ' point    coordinates     relation'
       sumwgt=0. 
       do 51 j=1,nmshp
       wgt=1.
       if(mod(in(1,j),n(1)).eq.0) wgt=wgt/2.                            
       if(mod(in(2,j),n(2)).eq.0) wgt=wgt/2.                            
       if(mod(in(3,j),n(3)).eq.0) wgt=wgt/2.                            
       wsum(num(j))=wsum(num(j))+wgt
       sumwgt=sumwgt+wgt                                                
       if(iswitch.eq.1) &
       write(66,100) j,in(1,j),in(2,j),in(3,j),num(j),wgt,wsum(num(j))
100    format(i5,5x,3i4,i8,2f10.5)
  51   continue 
      wsumme=0
      nirr=0 
      write(66,150)
      do 151 j=1,nmshp
      if(num(j).eq.j) then
      nirr=nirr+1 
        if(nirr.gt.idkp) then
          write(*,*) 'IDKP in main.f too small.'
          write(*,*) ' Set it to expected # of k-points in IBZ (gt.',nirr,')'
          stop 'reduz.f: IDKP in main.f too small'
        endif
      weight(nirr)=2.*wsum(num(j))/sumwgt
      if(iswitch.eq.1) &
      write(66,152) nirr,num(j),wsum(num(j)),weight(nirr)
      wsumme=wsumme+weight(nirr)
      endif       
  151 continue   
  150 format('  weights of k-points:')
  152 format(2i8,2f12.6)  
  153 format('  sum of weights: ',f12.6)
      if(iswitch.eq.1) write(66,153) wsumme 
      RETURN                                                            
      END                                                               
