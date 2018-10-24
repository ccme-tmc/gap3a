!BOP
!
! !ROUTINE: expand_evec
!
! !INTERFACE:
      subroutine expand_evec(ik,iop,ifrot,isp)

! !DESCRIPTION:

!Calculates the coeficients of the expansion of the eigenvectors of the DFT 
!calculation within the MT-spheres in the spherical basis as follows:
!
!\begin{equation}
!\Psi_{n\vec{k}}(\vec{r})=\sum\limits_a{\sum\limits_{l=0}^{lmax}{%
!\sum\limits_{m=-l}^{l}{\left[\mathcal{A}_{lm}^{na}(\vec{k})u_{al}(r)+%
!\mathcal{B}_{lm}^{na}(\vec{k})\dot{u}_{al}(r)+\mathcal{C}_{lm}^{na}(\vec{k})u_{al}(r,E_2)%
!\right]Y_{lm}(\hat{r})}}}
!\end{equation}
!where $a$ indicates the atoms. The corresponding expresions for the
!coefficients are:
!
!\begin{subequations}\label{almevec}
!\begin{align}
!\mathcal{A}_{lm}^{na}(\vec{k})=&\sum\limits_{\vec{G}}{Z^n_{\vec{k}+\vec{G}}%
!A^a_{lm}(\vec{k}+\vec{G})}\\
!\mathcal{B}_{lm}^{na}(\vec{k})=&\sum\limits_{\vec{G}}{Z^n_{\vec{k}+\vec{G}}%
!B^a_{lm}(\vec{k}+\vec{G})}\\
!\mathcal{C}_{lm}^{na}(\vec{k})=&Z^n_{\vec{k}+\vec{G}_{LO}}%
!C^a_{lm}(\vec{k}+\vec{G}_{LO})
!\end{align}
!\end{subequations}
!
! The respective coefficients are stored as:
!
!\begin{subequations}
!\begin{align}
!\mathcal{A}_{lm}^{na}(\vec{k})&\rightarrow\texttt{alfa(a,n,i)}\\
!\mathcal{B}_{lm}^{na}(\vec{k})&\rightarrow\texttt{beta(a,n,i)}\\
!\mathcal{C}_{lm}^{na}(\vec{k})&\rightarrow\texttt{gama(a,n,i)}
!\end{align}
!\end{subequations}
! 
! with $i=l^2+l+m+1$
!
! !USES:
     
      use eigenvec,    only: zzk,zzq,alm,blm,clm,alfa,beta,gama,&
     &                       alfp,betp,gamp 
      use bands,       only: nv, ibgw,nbgw,nbmax 
      use constants,   only: czero,cone
      use kpoints,     only: kpirind
      use lapwlo,      only: nlo_tot,lmax,nlomax,nt,lomax,nlo_at 
      use recipvec,    only: maxngk
      use struk,       only: mult,ndf,nat
      use task,        only: time_lapack,time_evec
      

! !INPUT PARAMETERS:

      implicit none
      
      integer, intent(in) :: ik    ! index of the k-point      
      integer, intent(in) :: iop   !  =1 --- alfa, 2 --- alfp
      logical, intent(in)    :: ifrot ! true if local rotations are taken  into account
      integer, intent(in) :: isp   ! index for spin 

! !LOCAL VARIABLES:

      integer :: iat   ! counts inequivalent atoms
      integer :: idf   ! Counts all atoms
      integer :: ieq   ! counts equivalent atoms
      integer :: ie1   ! counts bands
      integer :: irk
      integer :: ilo   
      integer :: iv1   ! run over G (lapw basis functions
      integer :: l1    ! orbital angular momentum quantum number
      integer :: ilm  ! l1^2+l1+m1+1
      integer :: m1    ! z component of the orbital ang. mom
      integer :: ngk   ! nuumbe of G-vector for a particular k-point 
      integer :: nlm
      integer    :: ierr
      logical :: ldbg=.false.

      real(8) :: time1,time2,tstart,tend
 
! !EXTERNAL ROUTINES: 

      external setabc_new

! !REVISION HISTORY:
! 
! Created Dic. 2003 by RGA
! Last Modification: July. 20th. 2004 by RGA
!
!EOP  
!BOC  
!
!     allocate the arrays to store the A_lm, B_lm and C_lm coefficients
!
      if(ldbg) write(6,*) "---- start expand_evec ---- "
      call cpu_time(tstart)
      irk=kpirind(ik)
      ngk=nv(irk)+nlo_tot

      allocate(alm(ngk,nt*nt,ndf),blm(ngk,nt*nt,ndf),stat=ierr)
      call errmsg(ierr.ne.0,'expand_evec','Fail to allocate alm etc.')
      alm = czero
      blm = czero

      if(nlomax.gt.0) then 
        allocate(clm(ngk,nLOmax,(lomax+1)**2,ndf), stat=ierr)
        call errmsg(ierr.ne.0,'expand_evec','Fail to allocate clm.')
        clm = czero 
      endif 

!
!       Calculate the coefficients A_lm, B_lm and C_lm
!
      if(ldbg) write(6,*) "call set_lapwcoef"
      call set_lapwcoef(ik,iop,ifrot,isp) 

      call cpu_time(time1)
      if(ldbg) write(6,*) "call zgemm"
      if(iop.eq.1)then
        alfa = czero
        beta = czero 
        nlm=(lmax+1)**2*ndf
        call zgemm('t','n',nbmax,nlm,ngk,cone,zzk,maxngk,alm,ngk,cone,&
     &             alfa,nbmax)
        call zgemm('t','n',nbmax,nlm,ngk,cone,zzk,maxngk,blm,ngk,cone,&
     &             beta,nbmax)

        if(nlomax.gt.0) then 
          gama = czero 
          nlm=(lomax+1)**2*ndf*nLOmax
          call zgemm('t','n',nbmax,nlm,ngk,cone,zzk,maxngk,clm,ngk,cone,&
     &             gama,nbmax)
        endif 

      else 
        alfp = czero
        betp = czero
        nlm=(lmax+1)**2*ndf
        call zgemm('t','n',nbmax,nlm,ngk,cone,zzq,maxngk,alm,ngk,cone,&
     &             alfp,nbmax)
        call zgemm('t','n',nbmax,nlm,ngk,cone,zzq,maxngk,blm,ngk,cone,&
     &             betp,nbmax)

        if(nlomax.gt.0) then 
          gamp = czero
          nlm=(lomax+1)**2*ndf*nLOmax
          call zgemm('t','n',nbmax,nlm,ngk,cone,zzq,maxngk,clm,ngk,cone,&
     &             gamp,nbmax)
        endif 

      endif 
      call cpu_time(time2)
      time_lapack=time_lapack+time2-time1

      if(ldbg) write(6,*) "---- end expand_evec ---- "
      deallocate(alm,blm)
      if(nlomax.gt.0) deallocate(clm) 

      call cpu_time(tend)
      time_evec=time_evec+tend-tstart
      end subroutine 
!EOC
