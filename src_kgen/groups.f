      SUBROUTINE GROUPS(io,nsym,ibrava,npgr)                            
!*---------------------------------------------------------------------*
!*      IN SUBR. GROUPS ,                                               
!*                                                                      
!*      READ(6,5000)NPGR,INV                                            
!* 5000 FORMAT(A4,I2)                                                   
!*                                                                      
!*  NPGRP   : NAME OF THE POINT-GROUP.                                  
!*         SHOULD BE:      SCHOENFLIS SYMBOL      CRYSTAL SYSTEM        
!*                                                                      
!*               1    <======>   C1                                     
!*              -1               CI                (TRICLINIC)          
!*                                                                      
!*               2               C2                                     
!*               M               CS                                     
!*              2/M              C2H               (MONOCLINIC)         
!*                                                                      
!*              222              D2                                     
!*              MM2              C2V                                    
!*              MMM              D2H               (ORTHOROMBIC)        
!*                                                                      
!*               4               C4                                     
!*              -4               S4                                     
!*              4/M              C4H                                    
!*              422              D4                                     
!*              4MM              C4V                                    
!*              -42M             D2D                                    
!*              4MMM             D4H               (TETRAGONAL)         
!*                                                                      
!*               3               C3                                     
!*              -3               S6                                     
!*              32               D3                                     
!*              3M               C3V                                    
!*              -3M              D3D               (RHOMBOHEDRAL)       
!*                                                                      
!*               6               C6                                     
!*              -6               C3H                                    
!*              6/M              C6H                                    
!*              622              D6                                     
!*              6MM              C6V                                    
!*              -62M             D3H                                    
!*              6MMM             D6H               (HEXAGONAL)          
!*                                                                      
!*               23              T                                      
!*               M3              TH                                     
!*              432              O                                      
!*              -43M             TD                                     
!*              M3M              OH                (CUBIC)              
!*                                                                      
!*  INV     : =1 LATTICE HAS INVERSION SYM., =0 NO INVERSION            
!*                                                                      
!*----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)                                         
!     **                                                                
!     **  CONSTRUCTION OF TRANSFORMATION MATRICES WHICH TRANSFORM     **
!     **  THE RECIPROCAL LATTICE VECTORS UNDER THE OPERATIONS OF ONE  **
!     **  OF THE POINT GROUPS.                                        **
!     **  INPUT IS THE NAME OF THE POINT GROUP SEE UNDER LABEL IN     **
!     **  TABLE 1.3 OF BRADLEY AND CRACKNELL : THE MATHEMATICAL       **
!     **  THEORY OF SYMMETRY IN SOLIDS (OXFORD)                       **
!     **  INPUT : INV = 0 IF (NO INVERSION SYMMETRY AND SPIN-ORBIT    **
!     **  COUPLING INCLUDED) ELSE INV = 1                             **
!     **  TRANSFORMATION MATRICES STORED IN IO ;                      **
!     **  G(I) = SUM(K) ( G(K) * IO(K,I) )                            **
!     **  G(I) IS THE I'TH RECIPROCAL TRANSLATION VECTOR              **
!     **  SUBROUTINES USED:                                           **
!     **    NONE                                                      **
!     ..................................................................
!      COMMON /TRAN/ IO(3,3,50), NSYM                                   
      CHARACTER*4 NPGRP,NPGR                                            
      dimension IO(3,3,48)
      DIMENSION NPGRP(32), NGELEM(32), IEQUIV(32), NBRAVA(14)           
      DIMENSION NSYMEL(4,15,24), ITRANS(9,68), IIO(9,48)                
!      EQUIVALENCE (IO(1,1,1),IIO(1,1))                                 
!     THESE ARRAYS CONTAIN:                                             
!     NPGRP  : THE NAME OF THE POINT GROUP.                             
!     NGELEM : THE ORDER OF THE GROUP - WITHOUT THE INVERSION ELEMENTS  
!     IEQUIV : THE 15 POINT GROUPS WHICH CANNOT BE OBTAINED FROM ONE    
!              OF THE OTHERS.                                           
!     NBRAVA : IN CASE OF MORE BRAVAIS LATTICES FOR  A GIVEN CRYSTAL    
!              CLASS THIS LABELS THE BRAVAIS LATTICE.                   
!     NSYMEL : LABLE FOR THE SYMMETRY OPERATIONS                        
!     ITRANS : SYMMETRY OPERATIONS.                                     
      DATA NPGRP /'1','-1','2','M','2/M','222','MM2','MMM','4','-4','4/M' &
      ,'422','4MM','-42M','4MMM','3','-3','32','3M','-3M','6','-6','6/M' &
      ,'622','6MM','-62M','6MMM','23','M3','432','-43M','M3M'/         
      DATA IEQUIV /1,1,2,3,2,4,5,4,6,7,6,6,8,7,6,9,9,9,10,9,11,12,11,11, &
       13,12,11,14,14,14,15,14/,             IWO /1/                    
!    1 13,12,11,14,14,14,15,14/, IO /450*0/, IWO /0/                    
      DATA NBRAVA /1,1,2,1,2,3,4,1,2,1,1,1,2,3/                         
      DATA NGELEM /1,1,2,2,2,4,4,4,4,4,4,8,8,8,8,3,3,6,6,6,6,6,6,12,12,  &
       12,12,12,12,24,24,24/                                             
      DATA ((ITRANS(I,J),I=1,9),J=1,12) /1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0 &
       ,0,0,1,-1,0,0,0,0,-1,0,-1,0,-1,0,0,0,1,0,0,0,-1,1,0,0,0,-1,0,0,0, &
       -1,0,1,0,1,0,0,0,0,-1,0,-1,0,-1,0,0,0,0,-1,0,0,-1,0,-1,0,-1,0,0,0 &
       ,1,0,-1,1,0,0,0,1,-1,1,0,-1,0,0,0,0,1,0,-1,0,1,-1,0,0,0,1,1,-1,0, &
       1,0,0,0,0,1/                                                     
      DATA ((ITRANS(I,J),I=1,9),J=13,24) /-1,1,0,0,1,0,0,0,-1,1,0,0,1,-1 &
       ,0,0,0,-1,1,-1,0,0,-1,0,0,0,-1,-1,0,0,-1,1,0,0,0,-1,0,0,1,0,-1,0, &
       1,0,0,-1,0,0,0,0,1,0,1,0,0,-1,1,0,-1,0,1,-1,0,-1,0,0,-1,0,1,-1,1, &
       0,0,1,-1,1,0,-1,0,0,-1,0,1,0,-1,0,0,0,0,1,1,0,-1,1,0,0,1,-1,0,0,1 &
       ,0,0,1,-1,-1,1,0/                                                
      DATA ((ITRANS(I,J),I=1,9),J=25,36) /-1,0,1,0,-1,1,0,0,1,0,1,0,0,0, &
       1,1,0,0,0,1,0,0,0,-1,-1,0,0,0,-1,0,0,0,1,-1,0,0,0,-1,0,0,0,-1,1,0 &
       ,0,1,0,0,0,0,1,0,-1,0,0,0,-1,0,1,0,1,0,0,-1,-1,-1,1,0,0,0,0,1,1,0 &
       ,0,-1,-1,-1,0,1,0,0,0,1,0,1,0,-1,-1,-1,0,1,0,-1,-1,-1,0,0,1,1,0,0 &
       ,0,0,1,-1,-1,-1/                                                 
      DATA ((ITRANS(I,J),I=1,9),J=37,43) /-1,-1,-1,0,1,0,1,0,0,0,0,-1,-1 &
       ,0,0,1,1,1,1,1,1,-1,0,0,0,-1,0,0,-1,0,1,1,1,-1,0,0,0,-1,0,0,0,-1, &
       1,1,1,-1,0,0,1,1,1,0,0,-1,1,1,1,0,-1,0,0,0,-1/                   
!     TRICLINIC : GT                                                    
      DATA (NSYMEL(1,1,I),I=1,1) /1/                                    
!     MONOCLINIC : GM                                                   
      DATA (NSYMEL(1,2,I),I=1,2) /1,2/                                  
      DATA (NSYMEL(1,3,I),I=1,2) /1,-2/                                 
!     MONOCLINIC : GMB                                                  
      DATA (NSYMEL(2,2,I),I=1,2) /1,3/                                  
      DATA (NSYMEL(2,3,I),I=1,2) /1,-3/                                 
!     ORTHOROMBIC : GO                                                  
      DATA (NSYMEL(1,4,I),I=1,4) /1,4,5,2/                              
      DATA (NSYMEL(1,5,I),I=1,4) /1,-4,-5,2/                            
!     ORTHOROMBIC : GOB                                                 
      DATA (NSYMEL(2,4,I),I=1,4) /1,6,7,2/                              
      DATA (NSYMEL(2,5,I),I=1,4) /1,-6,-7,2/                            
!     ORTHOROMBIC : GOV                                                 
      DATA (NSYMEL(3,4,I),I=1,4) /1,19,20,21/                           
      DATA (NSYMEL(3,5,I),I=1,4) /1,-19,-20,21/                         
!     ORTHOROMBIC : GOF                                                 
      DATA (NSYMEL(4,4,I),I=1,4) /1,44,45,46/                           
      DATA (NSYMEL(4,5,I),I=1,4) /1,-44,-45,46/                         
!     TETRAGONAL : GQ                                                   
      DATA (NSYMEL(1,6,I),I=1,8) /1,22,47,2,5,4,6,7/                    
      DATA (NSYMEL(1,8,I),I=1,8) /1,22,47,2,-5,-4,-6,-7/                
      DATA (NSYMEL(1,7,I),I=1,8) /1,-22,-47,2,5,4,-6,-7/                
!     TETRAGONAL : GQV                                                  
      DATA (NSYMEL(2,6,I),I=1,8) /1,23,24,21,20,19,25,7/                
      DATA (NSYMEL(2,8,I),I=1,8) /1,23,24,21,-20,-19,-25,-7/            
      DATA (NSYMEL(2,7,I),I=1,8) /1,-23,-24,21,20,19,-25,-7/            
!     TRIGONAL   : GRH                                                  
      DATA (NSYMEL(1,9,I),I=1,6) /1,26,51,3,8,7/                        
      DATA (NSYMEL(1,10,I),I=1,6) /1,26,51,-3,-8,-7/                    
!     HEXAGONAL  : GH                                                   
      DATA (NSYMEL(1,11,I),I=1,12) /1,9,12,10,11,2,13,14,7,15,16,6/     
      DATA (NSYMEL(1,13,I),I=1,12) /1,9,12,10,11,2,-13,-14,-7,-15,-16,-6 &
       /                                                                
      DATA (NSYMEL(1,12,I),I=1,12) /1,-9,-12,10,11,-2,13,14,7,15,16,6/  
!     CUBIC      : GC                                                   
      DATA (NSYMEL(1,14,I),I=1,24) /1,5,4,2,26,27,28,29,51,52,53,54,6,7, &
       17,18,8,3,30,31,22,55,56,47/                                     
!     CUBIC      : GCF                                                  
      DATA (NSYMEL(2,14,I),I=1,24) /1,45,44,46,26,32,33,34,51,35,36,37,  &
       50,7,42,43,8,3,38,39,49,40,41,48/                                 
!     CUBIC      : GCV                                                  
      DATA (NSYMEL(3,14,I),I=1,24) /1,20,19,21,26,60,61,62,51,57,58,59,  &
       65,66,23,63,64,24,25,7,67,68,8,3/                                 
!
!     write(*,*) 'pointgroup'
!     READ(*,5000)NPGR
      INV=1
!      write(*,*) ' Inversion-symmetry?   1 ... yes    0 ... no'
!      read(*,*) inv  
!      write(*,*) ' Ibrava (1-14, see bravais), 0 nimmt standard:'
!      read(*,*) ibrav1
!      if(ibrav1.ne.0) ibrava=ibrav1
 5000 FORMAT(A4,I2)                                                     

      DO 9 I1=1,48                                                      
      DO 9 I2=1,3                                                       
      DO 9 I3=1,3                                                       
 9    IO(I3,I2,I1)=0                                                    
      DO 10 J=1,3                                                       
      DO 10 I=1,24                                                      
      IS=1                                                              
      IF (I.GT.12) IS=-1                                                
   10 NSYMEL(J,15,I)=IS*NSYMEL(J,14,I)                                  
      DO 20 I=19,43                                                     
      ITRANS(1,I+25)=ITRANS(1,I)                                        
      ITRANS(2,I+25)=ITRANS(4,I)                                        
      ITRANS(3,I+25)=ITRANS(7,I)                                        
      ITRANS(4,I+25)=ITRANS(2,I)                                        
      ITRANS(5,I+25)=ITRANS(5,I)                                        
      ITRANS(6,I+25)=ITRANS(8,I)                                        
      ITRANS(7,I+25)=ITRANS(3,I)                                        
      ITRANS(8,I+25)=ITRANS(6,I)                                        
   20 ITRANS(9,I+25)=ITRANS(9,I)                                        
!      READ (5,30) NPGR,INV                                             
   30 FORMAT (A4,I2)                                                    
!!      WRITE (66,40) NPGR,INV                                            
   40 FORMAT (1H ,' NAME OF POINT GROUP =',A4,' INVERSION =',I2)        
      DO 50 I=1,32                                                      
      IF (NPGR.EQ.NPGRP(I)) GO TO 70                                    
   50 CONTINUE                                                          
!!      WRITE (66,60)                                                     
   60 FORMAT (1H ,' ERROR : POINT GROUP NAME  NOT IN THE LIST, CHECK : BRADLEY AND CRACKNELL ')                                         
!cc
      i=1
!cc      STOP                                                           
   70 IBRAVA=NBRAVA(IBRAVA)                                             
      NSYM=NGELEM(I)                                                    
      I1=IEQUIV(I)                                                      
      IF (I.EQ.2.OR.I.EQ.5.OR.I.EQ.8.OR.I.EQ.11.OR.I.EQ.15.OR.I.EQ.17 &
      .OR.I.EQ.20.OR.I.EQ.23.OR.I.EQ.27.OR.I.EQ.29.OR.I.EQ.32) INV=INV+1  
      IF (INV.EQ.2) INV=1                                               
      DO 100 J=1,NSYM                                                   
      JSIGN=ISIGN(1,NSYMEL(IBRAVA,I1,J))                                
      DO 80 K=1,9                                                       
   80 IIO(K,J)=JSIGN*ITRANS(K,IABS(NSYMEL(IBRAVA,I1,J)))                
      IF (INV.EQ.0) GO TO 100                                           
      DO 90 K=1,9                                                       
   90 IIO(K,J+NSYM)=-IIO(K,J)                                           
  100 CONTINUE                                                          
      IF (INV.EQ.1) NSYM=NSYM+NSYM                                      
      WRITE (66,110) NSYM                                               
  110 FORMAT (1H ,' NUMBER OF POINT GROUP OPERATIONS =',I5)             
      do 131 i=1,nsym 
      i3=0
      do 131 i1=1,3
      do 131 i2=1,3
      i3=i3+1
 131  io(i2,i1,i)=iio(i3,i)
      IF (IWO.EQ.1) THEN                                                
      DO 130 I=1,INT(FLOAT(NSYM)/4.+.9)                                 
      I1=4*I-3                                                          
      I2=4*I-2                                                          
      I3=4*I-1                                                          
      I4=4*I                                                            
      WRITE (66,120) I1,I2,I3,I4                                        
  120 FORMAT (T5,'SYMMETRY MATRIX NR.',I3,T30,'SYMMETRY MATRIX NR.'      &
       ,I3,T55,'SYMMETRY MATRIX NR.',I3,T80,'SYMMETRY MATRIX NR.',I3)   
      DO 130 J=1,3                                                      
      WRITE(66,140) (IO(J,K,I1),K=1,3),(IO(J,K,I2),K=1,3),(IO(J,K,I3),K= &
       1,3),(IO(J,K,I4),K=1,3)                                          
  130 CONTINUE                                                          
  140 FORMAT (T5,3I5,T30,3I5,T55,3I5,T80,3I5)                           
      ENDIF                                                             
      RETURN                                                            
      END                                                               
