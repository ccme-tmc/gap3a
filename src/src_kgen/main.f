!$hp9000_800 intrinsics on
!    symmetry operations read in form struct ,
!     modified by similarity transformation
!     since  iio for   2/m (ibrava=3) incorrect ?   
!     for b or latt bravais matrix changed
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!---  BLOCK : README ARBMSH                                         ----
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                                                                       
!  THE ARBMSH PACKAGE CONTAINES ROUTINES FOR BRILLOUIN ZONE INTEGRATIONS
!                                                                       
!  THE PACKAGE SHALL CONTAIN 5 SEPARATE FILES                           
!  1. TEST.FOR                                                          
!       TEST PROGRAM                                                    
!  2. ARBMSH.FOR                                                        
!       PROGRAM FOR CALCULATION OF IRREDUCIBLE K-POINTS AND TETRAHEDRA  
!  3. KINTR.FOR                                                         
!       PROGRAM FOR CALCULATION OF SAMPLING WEIGHTS FOR GIVEN ENERGIES  
!  4. MIXPIC.FOR                                                        
!       GENERAL ROUTINES USED IN ARBMSH.FOR AND KINTR.FOR               
!  5. SIC.DATA1                                                         
!       EXAMPLE INPUT DATA                                              
!                                                                       
! THE CODES IN ARBMSH(+MIXPIC) CREATE INFORMATION WHICH IS ONLY         
! DEPENDENT ON THE STRUCTURE OF THE CRYSTAL. THESE PROGRAMS RUN         
! BEFORE THE SELFCONSITENCY ITERATIONS. THEY PROVIDE THE IRREDUCIBLE    
! K-POINTS AND A FILE (15) FOR USE BY KINTR.                            
!                                                                       
! THE CODES IN KINTR(+MIXPIC) RUN IN EACH SELFCONSISTENCY ITERATION     
! THEY NEED AS INPUT THE ENERGY EIGENVALUES EB AND THE NUMBER OF        
! OCCUPIED STATES RNTOT. FOR SPIN POLARIZED CALCULATIONS TREAT          
! SPIN UP AND SPIN DOWN AS SEPARATE BANDS. OUTPUT ARE THE SAMPLING      
! WEIGHTS WGHT. THE INTEGRATION OF AN ARBITRARY FUNCTION A(N,K)         
! OVER THE OCCUPIED STATES IS PERFORMED BY SUMMATION                    
!         <A>=SUM OVER K AND N OF WGHT(N,K)*A(N,K)                      
! WHERE THE SUM RUNS OVER OCCUPIED AND! UNOCCUPIED STATES.              
! THE WEIGHTS CONTAINE BOTH THE GEOMETRICAL WEIGHT AND THE INFLUENCE    
! OF THE FERMI FUNCTION.                                                
!                                                                       
! THE EXAMPLE PRGRAM PERFORMES AN INTEGRATION OF                        
! MATRIXELEMENTS A(K)=COS(KX)+COS(KY)+COS(KZ)                           
! THE LATTICE IS SIMPLE CUBIC                                           
! FOR A HALF FILLED BAND WITH ENERGIES EB(K)=A(K)                       
! THE TOTAL NUMBER OF K-POINTS PER UNIT CELL IS 64                      
! THE NUMBER OF IRREDUCIBLE K-POINTS IS 4                               
! THE RESULT FOR <A> IS -0.5102                                         
! THE CONVERGED RESULT IS -0.5012                                       
!                                                                       
! REMARKS:                                                              
!                                                                       
! IN CASE OF SEMICONDUCTORS OR INSULATORS THE METHOD IS IDENTICAL       
! TO THE SPECIAL POINT SCHEME OF MONKHORST AND PACK.                    
!                                                                       
! FOR METALS IT IS IDENTICAL TO THE "TRADITIONAL" TETRAHEDRON METHOD    
! OF ANDERSEN AND JEPSEN, IF THE CORRECTION IS SWITCHED OFF.            
!                                                                       
! IF THE CORRECTION IS SWITCHED ON THE METHOD GIVES ALSO FOR METALS     
! RESULTS, WHICH ARE COMPARABLE TO RESULTS OF THE SPECIAL POINT         
! METHOD FOR INSULATORS.                                                
! IF THE CORRECTION IS KEPT ON THE ERROR DUE TO THE LINEAR              
! INTERPOLATION OF MATRIXELEMENTS A(N,K) IS CORRECTED FOR.              
!                                                                       
! THE CORRECTION FORMULA FOR LINEAR INTERPOLATION CAN BE SWITCHED ON    
! AND OFF BY THE PARAMETER "ICOR" IN THE SUBROUTINE "WEIGHTS"           
!                                                                       
! THE SYMMETRY OPERATIONS USED AS INPUT ARE TABULATED IN:               
!   C.J.BRADLEY AND A.P.CRACKNELL,                                      
!   THE MATHEMATICAL THEORY OF SYMMETRY IN SOLIDS,                      
!   OXFORD 1972                                                         
!                                                                       
!                                                                       
! SOMETIMES THE ROUTINE HAS PROBLEMS FINDING THE FERMI LEVEL. THIS      
! HAPPENS, IF THE FERMI LEVEL IS PINNED A A TETRAHEDRON WITH IDENTICAL  
! ENERGIES ON ALL 4 CORNERS, WHICH LEADS TO A DELTA PEAK IN THE         
! DENSITIES OF STATES. IN THIS CASE ONE SHOULD INCREASE THE K-MESH,     
! OR CHANGE THE K-MESH FROM AN EVEN NUMBER OF DIVISIONS FOR THE         
! RECIPROCAL LATTICE VECTORS TO AN ODD NUMBER OR VICE VERSA.            
!                                                                       
! THE FERMI LEVEL IS DETERMINED TO AN ACCURACY OF 1.E-5 IN THE          
! NUMBER OF STATES. THIS TOLERANCE CAN BE MODIFIED BY CHANGING THE      
! SETTING OF THE PARAMETER "TOLMAX" IN SUBROUTINE DOS.                  
!                                                                       
! THE PROGRAM CAN HANDLE UP TO APPROXIMATELY 130 IRREDUCIBLE K-POINTS   
! IF 4BYTE INTEGER NUMBERS ARE USED. IF AN INTEGER CONTAINS MORE THAN   
! 4BYTE CHANGE "IBGGST" IN SUBROUTINE TETCNT.                           
!                                                                       
! IF THE DIMENSION FOR THE NUMBER OF BANDS EXCEEDS THE ACTUAL NUMBER    
! OF BANDS, THE EIGENVALUES FOR THESE "NONEXISTING" BANDS MUST BE       
! GIVEN VALUES, WHICH LIE CERTAINLY ABOVE THE FERMI LEVEL.              
!                                                                       
! THE USE OF THE CRAY ROUTINE "ORDERS" IS OPTIONAL. THE ROUTINE         
! IS NOT INCLUDED IN THIS PACKAGE. THE CALLING STATEMENT MAY BE         
! COMMENTED OUT IF THE PARAMETER "ICRAY" IN SUBROUTINE "ORD1" IS 0.     
!                                                                       
! PETER BLOECHL                                                         
!                                                                       
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!----  BLOCK TEST                                                   ----
!----  TEST PROGRAM                                                 ----
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      PROGRAM MAIN                                                      
!                                                                       
      PARAMETER (NSYMX=50,IDKP=20000,NBX=2,NWX=20000000)                  
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)                               
      IMPLICIT INTEGER (O)                                              
      INTEGER W(NWX)   
      LOGICAL           ORTHO
      character*1 inout                                                 
      CHARACTER*4 LATTIC,ipgr                                           
      CHARACTER*11      STATUS,FORM                                     
      CHARACTER*80      FNAME                                           
      DIMENSION RBAS(3,3),AAA(3),NDIV(3),gbas(3,3)                      
      DIMENSION IARB(3),IIO(3,3,NSYMX),weight(idkp),wsum(nwx)  
      dimension klist(idkp,3),wei(idkp)   
      DIMENSION BK(3,IDKP)
      dimension iz(3,3,48),bki(3,idkp)
      common /inp/ iswitch
      save weight
!-----------------------------------------------------------------------
!                                                                       
      iarg=iargc()
      if(iarg.ne.1)STOP 'Exactly one commandline argument must be given'
      call getarg(1,fname)
      OPEN(1,FILE=fname,STATUS='OLD',ERR=8000)
 8003 READ(1,*,END=8001,ERR=8001) IUNIT,FNAME,STATUS,FORM,IRECL
      OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM, &
      ERR=8002)
      if(iunit.eq.20) then
        is=1
        do i=1,76
          if(fname(i:i+4).eq.'.ksym') then
            is=0
            exit
          endif
        enddo
      endif
      GOTO 8003
 8000 WRITE(*,*) ' ERROR IN OPENING kgen.DEF !!!!'
      STOP 'kgen.DEF'
 8002 WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
      WRITE(*,*) '       FILENAME: ',FNAME,'  STATUS: ',STATUS, &
      '  FORM:',FORM
      STOP 'OPEN FAILED'
 8001 CONTINUE 
! 
      PI=ACOS(-1.0)                                                     
!      write(*,*) 'LAPW  or ASW.....L/A  ?'
!      read(*,'(A)')  inout
       inout='L'
!      if(inout.eq.'A'.or.inout.eq.'a') goto 8888
!     ==================================================================
!     ==  INPUT  FROM FILE NFIL1=1  FOR LAPW CALCULATIONS             ==
!     ==================================================================
!     ==== RBAS(I,J) = REAL SPACE LATTICE VECTORS (I=X,Y,Z)             
      READ(20,1511) LATTIC,NAT,ipgr                                    
      READ(20,1512) A1,A2,A3,alpha,beta,gamma                           
      AAA(1)=A1                                                         
      AAA(2)=A2                                                         
      AAA(3)=A3 
      if(gamma.eq.0.d0) then
        gamma=pi/2.
      else
        gamma=gamma*pi/180.  
      end if                                                        
      if(beta.eq.0.d0) then
        beta=pi/2.
      else
        beta=beta*pi/180.  
      end if                                                        
      if(alpha.eq.0.d0) then
        alpha=pi/2.
      else
        alpha=alpha*pi/180.  
      end if                                                        
 1511 FORMAT(/,A4,23X,I3,/,30x,a4)                                      
 1512 FORMAT(6F10.7)                                                    
 1012 format(15x,i2)
 1101 FORMAT(3(3I2,/),I8)                                           
 1151 FORMAT(I4)                                                        
      DO 1 JATOM = 1,NAT                                                
         READ(20,*)
         read(20,1012) MULT 
            DO 125 MU=1,MULT-1                                   
               READ(20,*)          
 125        CONTINUE                                                    
         READ(20,*) 
         READ(20,*)
         READ(20,*)
         READ(20,*)
 1    CONTINUE                                                          
!     READ SYMMETRY-OPERATION FROM TAPE20=POTE                          
      READ(20,1151) IORD                                                
      DO 2 J=1,IORD                                                     
  2   READ(20,1101) ( (IZ(J1,J2,J),J1=1,3),J2=1,3 ),INUM   
!     ALL OF POTE IS READ                                               
!                                                                       
      GOTO 8889                       
8888  CONTINUE
!     ==================================================================
!     ==  INPUT FROM ASW - CNTRL FILE                                 ==
!     ==================================================================
! 
      READ(21,33) lattic,ipgr,a1,a2,a3,gamma                            
      AAA(1)=A1                                                         
      AAA(2)=A2                                                         
      AAA(3)=A3 
      gamma=gamma*pi/180.  
 33   FORMAT(26x,a3,2x,a4,/4F10.5)                                      
!
!     == NKP = NUMBER OF K-POINTS IN THE WHOLE UNIT CELL                
8889  CALL BRAVAI(lattic,A1,A2,A3,rbas,gbas,afact,IARB,ibrava, &
                  alpha,beta,gamma,ortho)          
      call groups(iio,nsym,ibrava,ipgr)
      call sdef(iio,nsym,lattic)
      CALL GBASS(RBAS,GBAS)   
      call addinv(iord,iz,is)                                          
      if(inout.ne.'A'.and.inout.ne.'a')  &
         call sdefl(rbas,gbas,iio,nsym,iz,iord,lattic,ortho)
!
      WRITE(*,*) ' NUMBER OF K-POINTS IN WHOLE CELL: (0 to specify&
     & 3 divisions of G)'                   
      READ(*,*)NKP
      nkp0=nkp                                                      
      write(*,*) ' SHORT OUTPUT: 0;    LONG OUTPUT: 1'
      read(*,*) iswitch
!       iswitch=1
!     ==================================================================
!     ==  FIND IRREDUCIBLE K-POINTS AND TETRAHEDRA                    ==
!     ==================================================================
      CALL ARBMSH(RBAS,gbas,IDKP,NKP,BK,NSYM,IIO,IARB,W,NWX,NDIV, &
       weight,wsum,sumwgt,bki)       
 1522 FORMAT(4f10.7)                                              
      IDIV=10000  
      if(iarb(1).eq.1.and.iarb(2).eq.1.and.iarb(3).eq.1)IDIV=NDIV(1)*2  
      if(iarb(1).eq.1.and.iarb(3).eq.0) IDIV=NDIV(1)*ndiv(3)*2          
      if(iarb(3).eq.1.and.iarb(1).eq.0) IDIV=NDIV(1)*ndiv(2)*2          
      if(iarb(2).eq.1.and.iarb(1).eq.0) IDIV=NDIV(2)*ndiv(3)*2          
      if(iarb(1).eq.0.and.iarb(2).eq.0.and.iarb(3).eq.0) IDIV=NDIV(1)* &
       ndiv(2)*ndiv(3)*2                                                
      mult=1
      write(66,*) 'NKP,NDIV,afact ',NKP,NDIV,afact                      
      write(*,*) nkp,' k-points generated, ndiv=',ndiv
      DO 300 IKP=1,NKP                                                  
      ak1=bk(1,ikp)
      ak2=bk(2,ikp)
      ak3=bk(3,ikp)
      ak1=ak1/2./pi*aaa(1)
      ak2=ak2/2./pi*aaa(2)
      ak3=ak3/2./pi*aaa(3)
 6223 WRITE(66,6221) AK1,AK2,AK3,ak1*idiv,ak2*idiv,ak3*idiv
 6221 format(3f12.5,10x,3f13.5)
      K1=NINT(AK1*IDIV)                                                 
      K2=NINT(AK2*IDIV)                                                 
      K3=NINT(AK3*IDIV)
!      write(6,*) weight(ikp),sumwgt
      wei(ikp)=weight(ikp)*sumwgt/2.
      if(inout.eq.'A'.or.inout.eq.'a') then 
!     write(*,*) aaa(1),aaa(2),aaa(3),idiv,k1,k2,k3
        xk=float(k1)/float(idiv)                              
        yk=float(k2)/float(idiv)/aaa(2)*aaa(1)                          
        zk=float(k3)/float(idiv)/aaa(3)*aaa(1)                          
        klist(ikp,1)=xk
        klist(ikp,2)=yk
        klist(ikp,3)=zk
        wei(ikp)=weight(ikp)
!cc     WRITE(8,1522) xk,yk,zk,weight(ikp)
      else                               
        if (.not.ortho) then
          klist(ikp,1)=NINT(bKi(1,ikp)*IDIV)
          klist(ikp,2)=NINT(bKi(2,ikp)*IDIV)
          klist(ikp,3)=NINT(bKi(3,ikp)*IDIV)
        else
          klist(ikp,1)=K1
          klist(ikp,2)=K2
          klist(ikp,3)=K3
        endif
       endif
 300    continue
!
        call divisi(idkp,nkp,idiv,klist)
!
        do 301 ikp=1,nkp
        if(ikp.eq.1) then 
          WRITE(8,1523) IKP,(klist(ikp,ir),ir=1,3),IDIV, &
                        wei(ikp),-7.,1.5,nkp0,ndiv
        else
          WRITE(8,1520) IKP, (klist(ikp,ir),ir=1,3),idiv,wei(ikp)
        endif                               
 301  CONTINUE                                                          
      write(8,1521)
 1520 FORMAT(I10,4I10,f5.1)                                              
 1521 format('END',/)
 1523 FORMAT(I10,4I10,3f5.1,4x,i6,' k, div: (',3i3,')')                  
      STOP 'KGEN ENDS'                                                  
      END                                                               

