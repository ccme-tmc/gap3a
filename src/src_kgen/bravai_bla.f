      SUBROUTINE BRAVAI(latti,AX,BX,CX,rbas,gbas,afact,IARB,ibrava, &
       gamma,ortho)     
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      LOGICAL           ORTHO
!     **                                                              **
!     **  CONSTRUCTION OF TRANSLATION VECTORS : RBAS (AS COLUMNS)     **
!     **  AND RECIPROCAL LATTICE VECTORS      : GBAS (AS ROWS   )     **
!     **   V: VOLUME OF THE BRILLOUIN ZONE                            **
!     **  INPUT IS THE NAME OF THE BRAVAIS LATTICE ACCORDING TO       **
!     **  TABLE 3.3 OF BRADLEY AND CRACKNELL : THE MATHEMATICAL       **
!     **  THEORY OF SYMMETRY IN SOLIDS (OXFORD)                       **
!     **                                                              **
!     **  SUBROUTINES USED:                                           **
!     **    NONE                                                      **
!     ..................................................................
      DIMENSION     GBAS(3,3), IARB(3)                                  
      CHARACTER*4 LATTIC,LATTI                                          
      DIMENSION RBAS(3,3), EPS(3,3,3), LATTIC(14)
      DATA LATTIC /'GT','GM','GMB','GO','GOB','GOV','GOF','GQ','GQV','GR &
      H','GH','GC','GCF','GCV'/                                         
!     DATA GBAS, EPS /36*0.0/, IARB /3*1/                               
      DATA  EPS /27*0.0/
      afact=1.0
      DO 10 I=1,3                                                       
      IARB(I)=1                                                         
      DO 10 J=1,3                                                       
      GBAS(J,I)=0.0                                                     
   10 RBAS(J,I)=0.0E0                                                   
      DET=0.0E0                                                         
!     READ (5,30) LATTI                                                 
      DO 20 I=1,14                                                      
      IF (LATTI.EQ.LATTIC(I)) GO TO 50                                  
   20 CONTINUE                                                          
!     WRITE (66,40) LATTI                                               
!  30 FORMAT (A4)                                                       
!  40 FORMAT (1H ,' ERROR: THE NAME OF THE BRAVAIS LATTICE',A4,'IS NOT I
!    1N THE TABLE'/' CHECK TABLE 3.3 IN BRADLEY AND CRACKNELL  SEE COMME
!    2NT CARDS IN SUB. BRAVAIS')                                        
   50 CONTINUE                                                          
!     READ (5,60) AX,AY,AZ,BX,BY,BZ,CX,CY,CZ 
      ay=0.
      az=0.
      by=0.
      bz=0.
      cy=0.
      cz=0.
   60 FORMAT (3F15.10)                                                  
!     IBRAVA=I 
      IF(LATTI(1:1).EQ.'H') GOTO 170                                    
      IF(LATTI(1:1).EQ.'F'.AND.AX.EQ.BX) GOTO 190                       
      IF(LATTI(1:1).EQ.'F') GOTO 130                                    
!      IF(LATTI(1:1).EQ.'A'.AND.latti(2:2).eq.'O') GOTO 112          
!      IF(LATTI(1:1).EQ.'B'.AND.latti(2:2).eq.'O') GOTO 111          
!      IF(LATTI(1:1).EQ.'B'.AND.latti(2:2).eq.'M') GOTO  90          
      IF(LATTI(1:1).EQ.'B'.AND.AX.EQ.BX.AND.AX.EQ.CX) GOTO 200          
      IF(LATTI(1:1).EQ.'B'.AND.AX.EQ.BX) GOTO 150                       
      IF(LATTI(1:1).EQ.'B') GOTO 120                                    
      if(latti(1:1).eq.'P'.and.abs(gamma-1.570796d0).gt.0.0001) goto 80
      if(latti(1:1).eq.'C'.and.abs(gamma-1.570796d0).gt.0.0001) goto 90
      IF(LATTI(1:1).EQ.'S'.AND.AX.EQ.BX.AND.AX.EQ.CX) GOTO 180          
      IF(LATTI(1:1).EQ.'P'.AND.AX.EQ.BX.AND.AX.EQ.CX) GOTO 180          
      IF(LATTI(1:1).EQ.'S'.AND.AX.EQ.BX) GOTO 140                       
      IF(LATTI(1:1).EQ.'P'.AND.AX.EQ.BX) GOTO 140                       
      IF(LATTI(1:1).EQ.'S') GOTO 100                                    
      IF(LATTI(1:1).EQ.'P') GOTO 100                                    
      IF(LATTI(1:1).EQ.'C') GOTO 110                                    
      IF(LATTI(1:3).EQ.'MXZ') GOTO  90                                  
      IF(LATTI(1:3).EQ.'M  ') GOTO  80                                  
      IF(LATTI(1:1).EQ.'R') GOTO  160                                  
      GO TO (70,80,90,100,110,120,130,140,150,160,170,180,190,200), I   
!     TRICLINIC : GT                                                    
   70 RBAS(1,1)=AX                                                      
      RBAS(2,1)=AY                                                      
      RBAS(3,1)=AZ                                                      
      RBAS(1,2)=BX                                                      
      RBAS(2,2)=BY                                                      
      RBAS(3,2)=BZ                                                      
      RBAS(1,3)=CX                                                      
      RBAS(2,3)=CY                                                      
      RBAS(3,3)=CZ                                                      
      IARB(1)=0                                                         
      IARB(2)=0                                                         
      IARB(3)=0                                                         
      ibrava=1
      ORTHO=.FALSE.                                             
      GO TO 210                                                         
!     MONOCLINIC : GM                                                   
!  AX = A * SIN(GAMMA)                                                  
!  AY = A * COS(GAMMA)
 80   ay=ax*cos(gamma)
      ax=ax*sin(gamma)
      rbas(1,1)=ax
!ccc      rbas(1,2)=ay
      rbas(2,1)=ay
!ccc
      rbas(2,2)=bx
!      RBAS(2,1)=-BX                                                    
!      RBAS(1,2)=AX                                                     
!      RBAS(2,2)=AY                                                     
      RBAS(3,3)=CX                                                      
      IARB(1)=0                                                         
      IARB(2)=0                                                         
      IARB(3)=0                                                         
      ORTHO=.FALSE.                                             
      ibrava=2
      GO TO 210                                                         
!     MONOCLINIC : GMB                                                  
!     AX = A * SIN(GAMMA) / 2                                           
!     AY = A * COS(GAMMA) / 2                                           
!     CX = C / 2                                                        
 90   ay=ax*cos(gamma)/2.
      ax=ax*sin(gamma)/2. 
      cx=cx/2.
      rbas(1,1)=ax
!ccc      rbas(1,2)=ay
!ccc      rbas(1,3)=-cx
!ccc      rbas(2,2)=bx
!ccc      RBAS(3,1)=aX                                                      
!ccc      RBAS(3,2)=ay                                                      
      rbas(2,1)=ay
      rbas(3,1)=-cx
      rbas(2,2)=bx
      RBAS(1,3)=aX                                                      
      RBAS(2,3)=ay                                                      
!ccc
      RBAS(3,3)=CX                                                      
      IARB(1)=0                                                         
      IARB(3)=0                                                         
      ORTHO=.FALSE.                                             
      ibrava=3
      GO TO 210                                                         
!     ORTHORHOMBIC : GO
!     modified!                                                 
  100 RBAS(2,2)=BX                                                     
      RBAS(1,1)=AX                                                      
      RBAS(3,3)=CX                                                      
      IARB(1)=0                                                         
      IARB(2)=0                                                         
      IARB(3)=0                                                         
      ORTHO=.TRUE.                                             
      ibrava=4
      GO TO 210                                                         
!     ORTHORHOMBIC : GOB                                                
  110 IF(LATTI(2:3).EQ.'XZ') GOTO 111                                   
      IF(LATTI(2:3).EQ.'YZ') GOTO 112                                   
      RBAS(1,1)=AX*0.5E0                                                
      RBAS(1,2)=-BX*0.5E0                                               
      RBAS(2,1)=AX*0.5E0                                                
      RBAS(2,2)=BX*0.5E0                                                
      RBAS(3,3)=CX                                                      
      IARB(2)=0                                                         
      IARB(3)=0 
      ORTHO=.TRUE.                                             
      ibrava=5
      GO TO 210                                                         
 111  RBAS(1,1)=AX*0.5E0                                                
      RBAS(1,3)=-CX*0.5E0                                               
      RBAS(3,1)=AX*0.5E0                                                
      RBAS(3,3)=CX*0.5E0                                                
      RBAS(2,2)=BX                                                      
      IARB(1)=0                                                         
      IARB(3)=0                                                         
      ORTHO=.TRUE.                                             
      ibrava=5
      GO TO 210                                                         
 112  RBAS(2,2)=BX*0.5E0                                                
      RBAS(2,3)=-CX*0.5E0                                               
      RBAS(3,2)=BX*0.5E0                                                
      RBAS(3,3)=CX*0.5E0                                                
      RBAS(1,1)=AX                                                      
      IARB(1)=0                                                         
      IARB(2)=0                                                         
      ORTHO=.TRUE.                                             
      ibrava=5
      GO TO 210                                                         
!     ORTHORHOMBIC : GOV                                                
!     da det =0 war vorzeichen von (1,1),(1,2),(3,2),(2,3) geaendert !!
  120 AX=AX*0.5E0                                                       
      BX=BX*0.5E0                                                       
      CX=CX*0.5E0                                                       
      IARB(1)=0                                                        
      IARB(2)=0                                                        
      IARB(3)=0                                                        
      RBAS(1,1)=-AX                                                     
!cc      RBAS(2,1)=BX                                                      
!cc      RBAS(3,1)=CX                                                      
!cc      RBAS(1,2)=+AX                                                     
!cc      RBAS(2,2)=-BX                                                     
!cc      RBAS(3,2)=+CX                                                     
!cc      RBAS(1,3)=AX                                                      
!cc      RBAS(2,3)=+BX                                                     
      RBAS(1,2)=BX                                                      
      RBAS(1,3)=CX                                                      
      RBAS(2,1)=+AX                                                     
      RBAS(2,2)=-BX                                                     
      RBAS(2,3)=+CX                                                     
      RBAS(3,1)=AX                                                      
      RBAS(3,2)=+BX                                                     
!cc
      RBAS(3,3)=-CX                                                     
      afact=0.5
      ORTHO=.TRUE.                                             
      ibrava=6
      GO TO 210                                                         
!     ORTHOROMBIC : GOF                                                 
  130 AX=AX*0.5E0                                                       
      BX=BX*0.5E0                                                       
      CX=CX*0.5E0                                                       
      RBAS(1,1)=AX                                                      
!cc      RBAS(3,1)=CX                                                      
!cc      RBAS(2,2)=-BX                                                     
!cc      RBAS(3,2)=CX                                                      
!cc      RBAS(1,3)=AX                                                      
!cc      RBAS(2,3)=-BX                                                     
      RBAS(1,3)=CX                                                      
      RBAS(2,2)=-BX                                                     
      RBAS(2,3)=CX                                                      
      RBAS(3,1)=AX                                                      
      RBAS(3,2)=-BX                                                     
!cc
      afact=0.5
      ORTHO=.TRUE.                                             
      ibrava=7
      GO TO 210                                                         
!     TETRAGONAL : GQ                                                   
  140 RBAS(1,1)=AX                                                      
      RBAS(2,2)=AX                                                      
      RBAS(3,3)=CX                                                      
      IARB(2)=0                                                         
      IARB(3)=0                                                         
      ORTHO=.TRUE.                                             
      ibrava=8
      GO TO 210                                                         
!     TETRAGONAL : GQV                                                  
  150 AX=AX*0.5E0                                                       
      CX=CX*0.5E0                                                       
      RBAS(1,1)=-AX                                                     
!cc      RBAS(2,1)=AX                                                      
!cc      RBAS(3,1)=CX                                                      
!cc      RBAS(1,2)=AX                                                      
!cc      RBAS(2,2)=-AX                                                     
!cc      RBAS(3,2)=CX                                                      
!cc      RBAS(1,3)=AX                                                      
!cc      RBAS(2,3)=AX                                                      
      RBAS(1,2)=AX                                                      
      RBAS(1,3)=CX                                                      
      RBAS(2,1)=AX                                                      
      RBAS(2,2)=-AX                                                     
      RBAS(2,3)=CX                                                      
      RBAS(3,1)=AX                                                      
      RBAS(3,2)=AX                                                      
!cc
      RBAS(3,3)=-CX                                                     
      IARB(2)=0                                                         
      IARB(3)=0                                                         
      afact=0.5
      ORTHO=.TRUE.                                             
      ibrava=9
      GO TO 210                                                         
!     TRIGONAL : GRH                                                    
!  160 RBAS(2,1)=-AX                                                    
!      RBAS(3,1)=CX                                                     
!      RBAS(1,2)=AX*SQRT(.75E0)                                         
!      RBAS(2,2)=AX*0.5E0                                               
!      RBAS(3,2)=CX                                                     
!      RBAS(1,3)=-AX*SQRT(.75E0)                                        
!      RBAS(2,3)=AX*0.5E0                                               
!      RBAS(3,3)=CX                                                     
!      ibrava=10
!      GO TO 210                                                        
! new definition
  160 RBAS(1,1)=ax/2.d0/sqrt(3.d0)
!cc      RBAS(1,2)=-AX/2.d0                                                
!cc      RBAS(1,3)=CX/3.d0                                                 
!cc      RBAS(2,1)=AX/2.d0/SQRT(3.d0)                                      
!cc      RBAS(2,2)=AX*0.5E0                                                
!cc      RBAS(2,3)=CX/3.d0                                                 
!cc      RBAS(3,1)=-AX/SQRT(3.d0)                                         
!cc      RBAS(3,2)=0.d0                                                
      RBAS(2,1)=-AX/2.d0                                                
      RBAS(3,1)=CX/3.d0                                                 
      RBAS(1,2)=AX/2.d0/SQRT(3.d0)                                      
      RBAS(2,2)=AX*0.5E0                                                
      RBAS(3,2)=CX/3.d0                                                 
      RBAS(1,3)=-AX/SQRT(3.d0)                                         
      RBAS(2,3)=0.d0                                                
!cc
      RBAS(3,3)=CX/3.d0                                                 
      ORTHO=.FALSE.                                             
      ibrava=10
      GO TO 210                                                         
!  160 RBAS(1,1)=ax/2.d0/sqrt(3.d0)
!      RBAS(2,1)=-AX/2.d0                                               
!      RBAS(3,1)=CX/3.d0                                                
!      RBAS(1,2)=AX/2.d0/SQRT(3.d0)                                     
!      RBAS(2,2)=AX*0.5E0                                               
!      RBAS(3,2)=CX/3.d0                                                
!      RBAS(1,3)=-AX/SQRT(3.d0)                                         
!      RBAS(2,3)=0.d0                                                
!      RBAS(3,3)=CX/3.d0                                                
!      ibrava=10
!      GO TO 210                                                        
!     HEXAGONAL : GH                                                    
!c  170 RBAS(2,1)=-AX                                                   
!c      RBAS(1,2)=AX*SQRT(.75E0)                                        
!c      RBAS(2,2)=AX*0.5E0                                              
!c      RBAS(3,3)=CX                                                    
!c    new definition !!
!cc  170 RBAS(1,2)=-AX/2.                                                  
!cc      RBAS(1,1)=AX*SQRT(.75E0)                                          
  170 RBAS(2,1)=-AX/2.                                                  
      RBAS(1,1)=AX*SQRT(.75E0)                                          
!cc
      RBAS(2,2)=AX                                                
      RBAS(3,3)=CX                                                      
      IARB(2)=0                                                         
      IARB(3)=0                                                         
      ORTHO=.FALSE.                                             
      ibrava=11
      GO TO 210                                                         
!     CUBIC : GC                                                        
  180 RBAS(1,1)=AX                                                      
      RBAS(2,2)=AX                                                      
      RBAS(3,3)=AX                                                      
      ORTHO=.TRUE.                                             
      ibrava=12
      GO TO 210                                                         
!     CUBIC : GCF                                                       
  190 AX=AX*0.5E0                                                       
      RBAS(2,1)=AX                                                      
      RBAS(3,1)=AX                                                      
      RBAS(1,2)=AX                                                      
      RBAS(3,2)=AX                                                      
      RBAS(1,3)=AX                                                      
      RBAS(2,3)=AX                                                      
      afact=0.5
      ORTHO=.TRUE.                                             
      ibrava=13
      GO TO 210                                                         
!     CUBIC : GCV                                                       
  200 AX=AX*0.5E0                                                       
      RBAS(1,1)=-AX                                                     
      RBAS(2,1)=AX                                                      
      RBAS(3,1)=AX                                                      
      RBAS(1,2)=AX                                                      
      RBAS(2,2)=-AX                                                     
      RBAS(3,2)=AX                                                      
      RBAS(1,3)=AX                                                      
      RBAS(2,3)=AX                                                      
      RBAS(3,3)=-AX                                                     
      afact=0.5
      ORTHO=.TRUE.                                             
      ibrava=14
  210 CONTINUE                                                          
      WRITE (66,*) 'ortho= ',ortho                                      
!     WRITE (66,220) LATTI,I                                            
  220 FORMAT (1H ,'  NAME OF BRAVAIS LATTICE = ',A4,' NO. = ',I3/' TRANS &
      LATION VECTORS (THEY CAN BE IN UNITS OF ONE LATTICE CONSTANT E.G.  &
      A, IN A.U., IN CM OR WHATEVER )')                                 
      DO 230 J=1,3                                                      
  230 WRITE (66,240) J,(RBAS(I,J),I=1,3)                                
  240 FORMAT (1H ,' R',I1,' = ',3F10.6)                                 
      PI=4.*ATAN(1.)                                                    
      EPS(1,2,3)=1.E0                                                   
      EPS(2,3,1)=1.E0                                                   
      EPS(3,1,2)=1.E0                                                   
      EPS(1,3,2)=-1.E0                                                  
      EPS(3,2,1)=-1.E0                                                  
      EPS(2,1,3)=-1.E0                                                  
      DO 250 I=1,3                                                      
      DO 250 J=1,3                                                      
      DO 250 K=1,3                                                      
      DET=DET+EPS(I,J,K)*RBAS(1,I)*RBAS(2,J)*RBAS(3,K)                  
      GBAS(I,1)=GBAS(I,1)+EPS(I,J,K)*RBAS(2,J)*RBAS(3,K)                
      GBAS(I,2)=GBAS(I,2)+EPS(I,J,K)*RBAS(3,J)*RBAS(1,K)                
  250 GBAS(I,3)=GBAS(I,3)+EPS(I,J,K)*RBAS(1,J)*RBAS(2,K)                
      DO 260 I=1,3                                                      
      DO 260 J=1,3                                                      
  260 GBAS(J,I)=2*PI*GBAS(J,I)/DET                                      
!      V=(2*PI)**3/DET                                                  
!      WRITE (66,270)                                                   
  270 FORMAT (1H ,'  RECIPROCAL LATTICE VECTORS IN UNITS CONSISTENT WITH &
       R')                                                              
!      DO 280 J=1,3                                                     
!  280 WRITE (66,290) J,(GBAS(J,I),I=1,3)                               
   290 FORMAT (1H ,' G',I1,' = ',3F10.6)                                
      WRITE (66,300) IARB                                               
  300 FORMAT (1H ,' DEPENDENCE OF DIVISION OF TRANSLATION VECTORS IARB=' &
       ,3I3)                                                            
      RETURN                                                            
      END                                                               
