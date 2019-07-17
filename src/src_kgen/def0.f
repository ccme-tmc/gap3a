!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     ----  BLOCK MIXPIC                                            ----
!     ----  MIXED PICKLES                                           ----
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     ..................................................................
      SUBROUTINE DEF0(NMAX)                                             
!     **                                                              **
!     **  ALLOCATES SPACE ON A (INTEGER WORK ARRAY)                   **
!     **  USE DEF0 BEFORE USE TO GIVE AVAILABLE SPACE ON WORK ARRAY   **
!     **  USE DEFI TO ALLOCATE SPACE FOR INTEGER ARRAYS               **
!     **  USE DEFDR TO ALLOCATE SPACE FOR INTEGER ARRAYS              **
!     **  USE RLSE TO RELEASE SPACE                                   **
!     **                                                              **
!     **  INPUT:                                                      **
!     **    NMAX        LENGTH OF INTEGER WORK ARRAY                  **
!     **    LENG        NUMBER OF ELEMENTS IN THE ARRAY TO BE         **
!     **                MAPPED ONTO THE WORK ARRAY (ENTRY: DEFI,DEFDR)**
!     **    ONAME       ALL ARRAYS FROM POINTER ONAME ARE DROPPED     **
!     **                (ENTRY: RLSE)                                 **
!     **  OUTPUT :                                                    **
!     **    ONAME       POINTER OF ARRAY TO BE ALLOCATED              **
!     **                                          (ENTRY: DEFI,DEFDR) **
!     **  REMARKS :                                                   **
!     **    AN INTEGER NUMBER IS ASSUMED TO HAVE 4 BYTES              **
!     **    A DOUBLE PRECISION NUMBER IS ASSUMED TO HAVE 8 BYTES      **
!     **                                                              **
      IMPLICIT INTEGER (O)                                              
      SAVE OMAX,OMAXX                                                   
      OMAX=1                                                            
      OMAXX=NMAX                                                        
      RETURN                                                            
!     ==================================================================
      ENTRY DEFI(ONAME,LENG)                                            
      ONAME=OMAX                                                        
      OMAX=OMAX+LENG                                                    
      IF(OMAX.GT.OMAXX) GOTO 9999                                       
      RETURN                                                            
!     ==================================================================
      ENTRY DEFDR(ONAME,LENG)                                           
      ONAME=OMAX                                                        
      OMAX=OMAX+LENG*2                                                  
      IF(OMAX.GT.OMAXX) GOTO 9999                                       
      RETURN                                                            
!     ==================================================================
      ENTRY RLSE(ONAME)                                                 
      IF(ONAME.LE.0.OR.ONAME.GT.OMAX) GOTO 9999                         
      OMAX=ONAME                                                        
      RETURN                                                            
9999  CONTINUE                                                          
      PRINT*,'ERROR IN DEF0,oname,omaxx,omax,leng', &
       oname,omaxx,omax,leng                                            
      print*,'Redimension NWX in main.f (gt.',omax,')'
      STOP                                                              
      END                                                               
