!-----------------------------------------------------------------------
!BOP
!
! !MODULE: barcoul --- bare coulomb interaction 
      module barcoul
      use constants, only: pi,czero
      use struk,     only: ndf,alat,pia,vi,br2,rbas
      use mixbasis,  only: mbsiz,lmbmax,mpwmix,mbsiz
      use kpoints,   only: nkp,nkdivs
      
! !PUBLIC VARIABLES:
      implicit none 
      integer :: iop_coul = -1    ! option to coulumb potential  
      !                     -1 -- standard bare Coulomb interaction 
      !                      0 -- truncated Coulomb interaction with truncation radius being rcut_coul
      !                      1 -- truncated Coulomb interaction for 1D system
      !                      2 -- truncated Coulomb interaction for 2D system
      !                      3 -- Thomas-Fermi screened Coulomb interaction 
      !                      4 -- erfc-screened Coulomb interaction 
      integer :: iop_coul_x = -1  ! option to coulumb potential for exchange  
      integer :: iop_coul_c = -1  ! option to coulumb potential for correlation 
      integer :: axis_cut_coul = 3  ! the axis related to truncation of Coulomb interaction
                               ! For 2D system, interaction is cutoff in the direction of axis_cut-th lattice vector
                               ! For 1D system, interaction is only preserved along axis_cut
      real(8) :: vec_cut_coul(3) ! the vector of the direction of axis_cut

      integer :: iop_coulvm     ! option to control how to calculate v-matrix 
                                !  0 -- the "standard" scheme 
                                !  1 -- via plane wave expansion (expt.)

      ! TODO check the usage of im_g0 with Prof. Jiang
      integer :: im_g0 = 1                ! the index for the basis function corresponding to G=0

      ! Experimental switch
      logical :: lcutoff_in_coul_barc = .False. ! truncate the Coulomb interaction
      ! this switch must be switched on when testing 1D and 2D Coulomb cutoff scheme

      real(8) :: barcevtol=-1.d-10        ! tolenrance to choose basis functions from bare Coulomb matrix eigenvectors 
      real(8) :: barcevtol2=0.001d0       ! barc eigenvectors that whose overlap with the G=0 plane wave   
                                          ! is larger than this value is kept for new basis set (used for iopmbq0 == 2 )
      real(8) :: rcut_coul=-1.0d0        ! cutoff radius of Coulomb interaction (for 0D)
      real(8) :: acut_coul=-1.0d0        ! cutoff radius of Coulomb interaction (for 1D)
      real(8) :: bcut_coul=-1.0d0        ! cutoff radius of Coulomb interaction (for 1D)
      real(8) :: zcut_coul=-1.0d0        ! cutoff radius of Coulomb interaction (for 2D)
                                          ! if negative, set the default value in terms of the cell size  
      real(8),allocatable::barcev(:)      ! Eigenvalues of bare Coulomb matrix 
      real(8),allocatable::barcevsq(:)    ! square root of eigenvalues of bare Coulomb matrix 
    
      real(8) :: barcsev, sqbarcsev 
      complex(8), allocatable :: barcvm(:,:)    ! transform matrix that diagonalized original bare Coulomb matrix 
      complex(8), allocatable :: barcs(:,:)     ! The matrix representation  of the bare coulomb potentianl in the mixed basis
      complex(8), allocatable :: sqbarcs(:,:)   ! The matrix representation of the bare coulomb potentianl in the mixed basis
      real(8), allocatable  :: jlam(:,:)        ! The matrix elements jlam
      real(8), allocatable ::  tilg(:)          ! The tildeg coefficients

      real(8),allocatable:: ev(:)               ! full set of the eigenvalues of barcoul matrix
      complex(8), allocatable:: vmat(:,:)       ! full set of eigenvectors of barcoul matrix 

!!    parameters related to structure constant (\Sigma)
      real(8) :: eta
      real(8) :: rcf                      ! cutoff for the real space summation
      real(8) :: gcf                      ! cutoff for the reciprocal space summation
      real(8) :: stctol                   ! convergence tolerance of the struct. const.
      complex(8), allocatable :: sgm(:,:) ! The sigma matrix
      real(8), allocatable :: rstr(:,:)   ! The lattice vectors used for the calculation of the structure constants

      real(8), private :: etx ! exp(-x)
      real(8), private :: sqx ! sqrt(x)

      interface gammaincc
        module procedure gammaincc_int
        module procedure gammaincc_hin
      end interface

      interface genrstr
        module procedure genrstr_no_cut
        module procedure genrstr_cutoff
      end interface

      contains
      
      subroutine init_barcoul(iq)
!
!     Initialize for bare Coulomb interaction 
!     iq=0: some q-independent initialization 
!     iq>0 -- q-depencent initialization 
      implicit none 
      integer,intent(in):: iq

      integer:: ierr
      integer:: ntm
      integer:: n,m
      real(8):: vtot
      character(20):: sname="init_barcoul"

      if(iq.eq.0) then !! q-independent initialization
        call linmsg(6,'-',"Initialize Coulomb matrix")
        ntm=lmbmax+1
        call calctildeg(2*ntm)
        n=ndf*(ndf+1)/2
        m=(4*ntm+1)*(4*ntm+1)
        allocate(sgm(m,n))

        eta=calceta()
        rcf=2.0d+0*rcutoff(stctol,eta,10)
        gcf=2.0d+0*gcutoff(stctol,eta,10)
        write(6,'(4x,a,3f12.6)') 'Parameters for structure constants:&
     &eta,rcf,gcf=',eta,rcf,gcf
        !! automatic set truncated/screened Coulomb interaction
        if(rcut_coul.lt.0.0) then 
          if(nint(rcut_coul).eq.-2) then
            rcut_coul=maxval(alat(1:3)*nkdivs(1:3))/2
          else
            vtot=(1.d0/vi)*nkp
            rcut_coul=(vtot/(pi*4.d0/3.d0))**(1.d0/3)
          endif
          write(6,'(4x,a,f12.4)') "set default rcut_coul=",rcut_coul
        endif 
      else  !! q-dependent initialization 
        allocate(vmat(mbsiz,mbsiz),ev(mbsiz),stat=ierr )
        call errmsg(ierr.ne.0,sname,"Fail to allocate vmat")
        vmat=0.d0
        ev = 0.d0 
      endif
      end subroutine 

      subroutine init_barcoul_2d(iq)
!
!     Initialize for bare Coulomb interaction 
!     iq=0: some q-independent initialization 
!     iq>0 -- q-depencent initialization 
      implicit none 
      integer,intent(in):: iq

      integer:: ierr
      integer:: ntm
      integer:: n,m
      real(8):: vtot
      character(20):: sname="init_barcoul_2d"

      if(iq.eq.0) then !! q-independent initialization
        call linmsg(6,'-',"Initialize Coulomb matrix (2D)")
        ntm=lmbmax+1
        call calctildeg(2*ntm)
        n=ndf*(ndf+1)/2
        m=(4*ntm+1)*(4*ntm+1)
        allocate(sgm(m,n))

        eta=calceta_2d()
        rcf=2.0d+0*rcutoff_2d(stctol,eta,10)
        gcf=2.0d+0*gcutoff_2d(stctol,eta,10)
        write(6,'(4x,a,3f12.6)') 'Parameters for structure constants:&
     &eta,rcf,gcf=',eta,rcf,gcf
        !! automatic set truncated/screened Coulomb interaction
        if(rcut_coul.lt.0.0) then 
          if(nint(rcut_coul).eq.-2) then
            rcut_coul=maxval(alat(1:3)*nkdivs(1:3))/2
          else
            vtot=(1.d0/vi)*nkp
            rcut_coul=(vtot/(pi*4.d0/3.d0))**(1.d0/3)
          endif
          write(6,'(4x,a,f12.4)') "set default rcut_coul=",rcut_coul
        endif 
        call set_coul_cutoff(2)
      else  !! q-dependent initialization 
        allocate(vmat(mbsiz,mbsiz),ev(mbsiz),stat=ierr )
        call errmsg(ierr.ne.0,sname,"Fail to allocate vmat")
        vmat=0.d0
        ev = 0.d0 
      endif
      end subroutine

      subroutine end_barcoul(iq)
        integer,intent(in):: iq
        if(iq.eq.0) then 
          deallocate(tilg,sgm)
        else 
          deallocate(vmat,ev)
        endif 
      endsubroutine 

      subroutine end_barcev
        deallocate(barcev,barcevsq,barcvm) 
        if(allocated(barcs)) deallocate(barcs,sqbarcs)
      end subroutine

      subroutine set_coul_cutoff(iop)
        ! Set up the truncation parameters by iop. 
        ! For meaning of each value
        implicit none
        integer, intent(in) :: iop
        integer :: iz, ia, ib
       
        vec_cut_coul = rbas(axis_cut_coul, :)
        if (iop.eq.-1) then
          ! very large cut-off length, for safety
          zcut_coul = 1.0D10
          acut_coul = 1.0D10
          bcut_coul = 1.0D10
        elseif ((iop.eq.0).or.(iop.eq.3).or.(iop.eq.4)) then
          ! spherical cut-off
          zcut_coul = 1.0D10
          acut_coul = 1.0D10
          bcut_coul = 1.0D10
        else
          ! 1D and 2D cut-off
          iz = mod(axis_cut_coul+2, 3)
          ia = mod(iz+1, 3)
          ib = mod(iz+2, 3)
          zcut_coul = alat(iz+1)/2.0D0
          acut_coul = alat(ia+1)/2.0D0
          bcut_coul = alat(ib+1)/2.0D0
          write(6,"(A25)") "Coulomb cutoff parameter:"
          write(6,12) "zcut", zcut_coul
          write(6,12) "acut", acut_coul
          write(6,12) "bcut", bcut_coul
        end if
12 format(10X, A8, F12.8)
      end subroutine set_coul_cutoff

      subroutine genrstr_cutoff(iop,rmax,rshift,rbs,nr)
! Generates the indexes of the lattice vectors to be included in the
!calculation of the structure constants, under the condition:
!
!\begin{equation}
!|\vec{R}+\vec{r}_{aa'}|\le R_{cutoff}
!\end{equation}
      implicit none
      integer, intent(in) :: iop  ! cut-off option
      real(8), intent(in) :: rmax ! Maximum radius
      real(8), intent(in) :: rshift(3) ! Shift of the origin
      real(8), intent(in) :: rbs(3,3) ! Bravais lattice basis
      
      integer(4), intent(out) :: nr  !number of vectors
      
      integer(4) :: i,ir,rdim,i1,i2,i3,gap,imax(3),iz,ia,ib
      
      real(8) :: lrmin              ! minimum length of the basis vectors.
      real(8) :: rleng
      
      real(8), dimension(3) :: r     ! vector belonging to the real space lattice
      real(8), dimension(3) :: lrbs ! length of the basis vectors
      real(8), dimension(3) :: rps
      real(8), dimension(4) :: rtmp
      logical :: done
      
      do i=1,3
        lrbs(i)=sqrt(sum(rbs(:,i)**2))
      enddo

      lrmin=minval(lrbs)
      imax=idint(rmax/lrmin)+1
      ! cutoff for 1D and 2D system
      iz = mod(axis_cut_coul+2,3)
      ia = mod(iz+1,3)
      ib = mod(iz+2,3)
      if (iop.eq.1) then
        imax(ia+1) = 0
        imax(ib+1) = 0
      elseif (iop.eq.2) then
        imax(iz+1) = 0
      end if
      rdim=(2*imax(1)+1)*(2*imax(2)+1)*(2*imax(3)+1)
     
      if(allocated(rstr)) deallocate(rstr) 
      allocate(rstr(4,rdim))
      
      ir=0
      do i1=-imax(1),imax(1)
        do i2=-imax(2),imax(2)
          do i3=-imax(3),imax(3)
            do i=1,3
              r(i) = dble(i1)*rbs(i,1)+dble(i2)*rbs(i,2)+               &
     &             dble(i3)*rbs(i,3)
            enddo
            rps=r+rshift
            rleng=sqrt(sum(rps*rps))

            if((rleng.le.rmax).and.(rleng.gt.1.0d-6))then
              ir=ir+1
              rstr(1:3,ir)=rps(1:3)
              rstr(4,ir)=rleng 
            endif
          enddo  
        enddo      
      enddo
      nr=ir

      !!sort by increasing length using shell algorithm
      gap=nr/2
      do while(gap.ge.1)
        done=.false.
        do while(.not.done)
          done=.true.
          do i=1,nr-gap
            if(rstr(4,i).gt.rstr(4,i+gap))then
              rtmp(1:4)=rstr(1:4,i)
              rstr(1:4,i)=rstr(1:4,i+gap)
              rstr(1:4,i+gap)=rtmp(1:4)
              done=.false.
            endif
          enddo
        enddo
        gap=gap/2
      enddo
      end subroutine genrstr_cutoff
      
      subroutine genrstr_no_cut(rmax,rshift,rbs,nr)
! Generates the indexes of the lattice vectors to be included in the
!calculation of the structure constants, under the condition:
!
!\begin{equation}
!|\vec{R}+\vec{r}_{aa'}|\le R_{cutoff}
!\end{equation}
      implicit none
      real(8), intent(in) :: rmax ! Maximum radius
      real(8), intent(in) :: rshift(3) ! Shift of the origin
      real(8), intent(in) :: rbs(3,3) ! Bravais lattice basis
      
      integer(4), intent(out) :: nr  !number of vectors
      
      integer(4) :: i,ir,rdim,i1,i2,i3,gap,imax
      
      real(8) :: lrmin              ! minimum length of the basis vectors.
      real(8) :: rleng
      
      real(8), dimension(3) :: r     ! vector belonging to the real space lattice
      real(8), dimension(3) :: lrbs ! length of the basis vectors
      real(8), dimension(3) :: rps
      real(8), dimension(4) :: rtmp
      logical :: done
      
      do i=1,3
        lrbs(i)=sqrt(sum(rbs(:,i)**2))
      enddo

      lrmin=minval(lrbs)
      imax=idint(rmax/lrmin)+1
      rdim=(2*imax+1)*(2*imax+1)*(2*imax+1)
     
      if(allocated(rstr)) deallocate(rstr) 
      allocate(rstr(4,rdim))
      
      ir=0
      do i1=-imax,imax
        do i2=-imax,imax
          do i3=-imax,imax
            do i=1,3
              r(i) = dble(i1)*rbs(i,1)+dble(i2)*rbs(i,2)+               &
     &             dble(i3)*rbs(i,3)
            enddo
            rps=r+rshift
            rleng=sqrt(sum(rps*rps))

            if((rleng.le.rmax).and.(rleng.gt.1.0d-6))then
              ir=ir+1
              rstr(1:3,ir)=rps(1:3)
              rstr(4,ir)=rleng 
            endif
          enddo  
        enddo      
      enddo
      nr=ir

      !!sort by increasing length using shell algorithm
      gap=nr/2
      do while(gap.ge.1)
        done=.false.
        do while(.not.done)
          done=.true.
          do i=1,nr-gap
            if(rstr(4,i).gt.rstr(4,i+gap))then
              rtmp(1:4)=rstr(1:4,i)
              rstr(1:4,i)=rstr(1:4,i+gap)
              rstr(1:4,i+gap)=rtmp(1:4)
              done=.false.
            endif
          enddo
        enddo
        gap=gap/2
      enddo
      end subroutine genrstr_no_cut

      function calceta() result(neta)
! This function calculates the optimal value of $\eta $ for the lattice
! summations needed to obtain the structure constants.
      implicit none

      integer(4) :: i
      real(8) :: mlg
      real(8) :: mlr
      real(8) :: neta
      real(8), dimension (3)   :: lgbs
      real(8), dimension (3)   :: lrbs
      intrinsic dsqrt
      intrinsic minval
      intrinsic isign
      do i=1,3
        lrbs(i)=dsqrt(rbas(i,1)*rbas(i,1)+rbas(i,2)*rbas(i,2)+&
     &          rbas(i,3)*rbas(i,3))
        lgbs(i)=dsqrt(br2(1,i)*br2(1,i)+br2(2,i)*br2(2,i)+&
     &         br2(3,i)*br2(3,i))
      enddo
      mlr=minval(lrbs)
      mlg=minval(lgbs)
      neta=dsqrt(2.0d0*mlr/mlg)
      end function calceta

      function calceta_2d() result(neta)
! This function calculates the optimal value of $\eta $ for the lattice
! summations needed to obtain the structure constants.
      implicit none

      integer(4) :: i
      real(8) :: mlg
      real(8) :: mlr
      real(8) :: neta
      real(8), dimension (3)   :: lgbs
      real(8), dimension (3)   :: lrbs
      intrinsic dsqrt
      intrinsic minval
      intrinsic isign
      do i=1,3
        lrbs(i)=dsqrt(rbas(i,1)*rbas(i,1)+rbas(i,2)*rbas(i,2)+&
     &          rbas(i,3)*rbas(i,3))
        lgbs(i)=dsqrt(br2(1,i)*br2(1,i)+br2(2,i)*br2(2,i)+&
     &         br2(3,i)*br2(3,i))
      enddo
      lrbs(axis_cut_coul) = 1.0d10
      lgbs(axis_cut_coul) = 1.0d10
      mlr=minval(lrbs)
      mlg=minval(lgbs)
      neta=dsqrt(2.0d0*mlr/mlg)
      end function calceta_2d

      function rcutoff(tol,eta,lambdamax) result(rcf)
!  Estimates the cutoff radius of the sums in real space for the
! calculation of the structure constants by the solving the equation:
!
!\begin{equation}
!\mathfrak{E}_{R,\lambda}^{\textrm{tol}}=\left\{%
!\begin{array}{ll}
!\frac{4\pi}{(\lambda -2)\Gamma(\lambda+\tfrac{1}{2})}%
!\left(\frac{\Gamma[\tfrac{\lambda}{2}+\tfrac{3}{2},\left(\tfrac{R_c}{\eta}\right)^2]}%
!{\eta^{\lambda-2}}-\frac{\Gamma[\lambda+\tfrac{1}{2},\left(\tfrac{R_c}{\eta}\right)^2]}%
!{R_c^{\lambda-2}}\right)&\lambda \neq 2\\
!\frac{4\pi}{\Gamma(\tfrac{5}{2})}\left[\tfrac{\eta}{R_c}%
!\Gamma[3,\left(\tfrac{R_c}{\eta}\right)^2]-\Gamma[\tfrac{5}{2},\left(\tfrac{R_c}{\eta}\right)^2]\right]&
!\lambda=2\\
!\end{array}
!\right.
!\end{equation}
!
! and taking the maximum value of $R_c$ obtained for $\lambda = 1...$
!\verb"lambdamax".      
!
! !USES:
      
      
      implicit none
      real(8),    intent(in) :: tol ! The tolerance for the convergence of the lattice sum
      real(8),    intent(in) :: eta 
      integer(4), intent(in) :: lambdamax

      real(8) :: rcf ! The maximum cutoff radius

      integer(4) :: l1
      integer(4) :: i
      real(8) :: rl
      real(8) :: x
      real(8) :: gmm
      real(8), allocatable :: eps(:)
      real(8) :: rnot
      real(8) :: gaml12
      real(8) :: gaml32
      real(8) :: prefac
      real(8), allocatable :: rct(:)
      real(8), external :: higam

      allocate(rct(lambdamax+1))
      allocate(eps(lambdamax+1))
      rnot=maxval(alat)
      rct = 5.0d+1
      do i=1,100
        x = 5.0d-1*dble(i)
        do l1=0,lambdamax
          if(l1.ne.2)then
            rl =5.0d-1*(dble(l1+1))
            gaml32=incgam(rl,x*x)
            rl = dble(l1)+5.0d-1
            gaml12=incgam(rl,x*x)
            gmm = higam(l1)
            prefac = 4.0d0*pi/(dble(l1-2)*gmm)
            eps(l1+1)=dabs(prefac*(gaml32-gaml12/(x**(l1-2)))/        &
     &              (eta**(l1-2)))
            if((eps(l1+1).lt.tol).and.(x.lt.rct(l1+1)))rct(l1+1)=x
          else
            gaml32=incgam(3.0d0,x*x)
            gaml12=incgam(2.5d0,x*x)
            gmm = higam(2)
            prefac = 4.0d0*pi/gmm
            eps(l1+1)= dabs(prefac*(gaml32/x-gaml12))
            if((eps(l1+1).lt.tol).and.(x.lt.rct(l1+1)))rct(l1+1)=x
          endif
        enddo ! l1
      enddo ! i
      rcf=maxval(rct)*eta
      
      deallocate(rct)
      deallocate(eps)
      end function rcutoff

      function rcutoff_2d(tol,eta,lambdamax) result(rcf)
! similar to rcutoff, but for 2D system
! 
! !USES:
      
      
      implicit none
      real(8),    intent(in) :: tol ! The tolerance for the convergence of the lattice sum
      real(8),    intent(in) :: eta 
      integer(4), intent(in) :: lambdamax

      real(8) :: rcf ! The maximum cutoff radius

      integer(4) :: l1
      integer(4) :: i, iz, ia, ib
      real(8) :: rl
      real(8) :: x               ! Rc/eta
      real(8) :: gmm
      real(8), allocatable :: eps(:)
      real(8) :: gaml1
      real(8) :: gaml2
      real(8) :: gaml3
      real(8) :: prefac
      real(8), allocatable :: rct(:)
      real(8), external :: higam

      allocate(rct(lambdamax+1))
      allocate(eps(lambdamax+1))
      iz = mod(axis_cut_coul+2,3)
      ia = mod(iz+1,3)+1
      ib = mod(iz+2,3)+1

      rct = 5.0d+1
      do i=1,100
        x = 5.0d-1*dble(i)
        do l1=0,lambdamax
          if(l1.ne.1)then
            rl = real(l1,8)*0.5d0+1.0d0
            gaml1=incgam(rl,x*x)
            rl = real(l1,8)+0.5d0
            gaml2=incgam(rl,x*x)
            gmm = higam(l1)
            prefac = 2.0d0*pi*vi*alat(axis_cut_coul)/gmm/real(1-l1,8)
            eps(l1+1)=abs(prefac/(eta**(l1-1))*(gaml1-gaml2/(x**(l1-1))))
            if((eps(l1+1).lt.tol).and.(x.lt.rct(l1+1)))rct(l1+1)=x
          else
            gaml1=incgam(2.0d0,x*x)
            gaml2=incgam(1.5d0,x*x)
            gmm = higam(1)
            prefac = 2.0d0*pi*vi*alat(axis_cut_coul)/gmm
            eps(l1+1)= dabs(prefac*(gaml1/x-gaml2))
            if((eps(l1+1).lt.tol).and.(x.lt.rct(l1+1)))rct(l1+1)=x
          endif
        enddo ! l1
      enddo ! i
      rcf=maxval(rct)*eta
      
      deallocate(rct)
      deallocate(eps)
      end function rcutoff_2d

      function gcutoff(tol,eta,lambdamax) result(rcf)
!  Estimates the cutoff radius of the sums in reciprocal space for the
!  calculation of the structure constants by the solving the equation:
!
!  \begin{equation}
!  \mathfrak{E}_{G,\lambda}^{\textrm{tol}}=\frac{8(\pi)^{\frac{5}{2}}}{\Omega%
!  \Gamma(\lambda+\tfrac{1}{2})\eta^{\lambda+1}}%
!  \Gamma\left[\tfrac{\lambda+1}{2},\left(\tfrac{\eta G_c}{2}\right)^2\right]
!  \end{equation}
!
!   and taking the maximum value of $G_c$ obtained for $\lambda = 1...$
!  \verb"lambdamax".      
      

! !INPUT PARAMETERS:

      implicit none

      real(8),    intent(in) :: tol ! The tolerance for the convergence of the lattice sum
      real(8),    intent(in) :: eta 
      integer(4), intent(in) :: lambdamax

      real(8) :: rcf ! The maximum cutoff radius

      integer(4) :: l1
      integer(4) :: i

      real(8) :: gaml32
      real(8) :: gmm
      real(8) :: prefac
      real(8) :: rl
      real(8) :: rnot
      real(8) :: x
      real(8), allocatable :: eps(:)
      real(8), allocatable :: rct(:)
      real(8), external :: higam

      allocate(rct(lambdamax+1))
      allocate(eps(lambdamax+1))
      rnot=maxval(pia)
      rct(:) = 5.0d+1
      do i=1,100
        x = 5.0d-1*dble(i)
        do l1=0,lambdamax
          rl =5.0d-1*(dble(l1+1))
          gaml32=incgam(rl,x*x)
          gmm = higam(l1)
          prefac = 8.0d0*pi*pi*dsqrt(pi)*vi/(gmm*eta**(l1+1))
          eps(l1+1)=dabs(prefac*gaml32)
          if((eps(l1+1).lt.tol).and.(x.lt.rct(l1+1)))rct(l1+1)=x
        enddo ! l1
      enddo ! i
      rcf=maxval(rct)*2.0d0/eta
      deallocate(rct)
      deallocate(eps)
      end function gcutoff

      function gcutoff_2d(tol,eta,lambdamax) result(rcf)
!  similar to gcutoff, but for 2d system

! !INPUT PARAMETERS:

      implicit none

      real(8),    intent(in) :: tol ! The tolerance for the convergence of the lattice sum
      real(8),    intent(in) :: eta 
      integer(4), intent(in) :: lambdamax

      real(8) :: rcf ! The maximum cutoff radius

      integer(4) :: l1
      integer(4) :: i
      integer(4) :: iz,ia,ib


      real(8) :: gaml1
      real(8) :: gamhalf
      real(8) :: gmm
      real(8) :: prefac
      real(8) :: rl
      real(8) :: rnot
      real(8) :: x
      real(8) :: b1b2    ! dot product of two periodic reciprocal basis
      real(8), allocatable :: eps(:)
      real(8), allocatable :: rct(:)
      real(8), allocatable :: yll(:)
      real(8), external :: higam

      iz=mod(axis_cut_coul+2,3)
      ia=mod(iz+1,3)+1
      ib=mod(iz+2,3)+1

      allocate(rct(lambdamax+1))
      allocate(eps(lambdamax+1))
      allocate(yll(lambdamax+1))
      ! set Y_{lam,lam}
      yll(1) = 1.0d0/sqrt(4.0d0*pi)
      do l1=1,lambdamax
        rl = real(l1,8)
        yll(l1+1) = yll(l1)*sqrt(1.0d0+0.5d0/rl)
      enddo
      ! set eps
      rct(:) = 5.0d+1
      do i=1,100
        x = 5.0d-1*real(i,8)
        gamhalf=incgam(0.5d0,x*x)
        do l1=0,lambdamax
          rl =5.0d-1*real(l1+2,8)
          gaml1=incgam(rl,x*x)
          gmm = higam(l1)
          prefac = 2.0d0/(gmm*eta**(l1+1))*yll(l1+1)/real(l1+1,8)
          eps(l1+1)=abs(prefac*(gaml1-x**(l1+1)*gamhalf))
          if((eps(l1+1).lt.tol).and.(x.lt.rct(l1+1)))rct(l1+1)=x
        enddo ! l1
      enddo ! i
      rcf=maxval(rct)*2.0d0/eta
      deallocate(rct)
      deallocate(eps)
      deallocate(yll)
      end function gcutoff_2d

      recursive function gammaincc_int(n,x) result(gmi)
      implicit none
      real(8) :: x ! value at which the gamma function is
      integer(4) :: n ! order of the gamma function
      real(8) :: gmi ! value of the gamma function
!
! Calculates the incomplete gamma function of integer order
!
      real(8), external :: e1xb
      if(n.eq.0)then
        gmi=e1xb(x)
      elseif(n.eq.1)then
        gmi = etx
      else
        gmi = (x**(n-1))*etx+dble(n-1)*gammaincc(n-1,x)
      endif
      end function gammaincc_int

      recursive function gammaincc_hin(in,x,n) result(gmh)
      implicit none
      integer(4) :: in ! order of the gamma function - 1/2
      real(8) :: x ! value at which the gamma function is calculated
      real(8) :: n ! order of the gamma function
      real(8) :: gmh ! value of the gamma function
!
!     Calculates the incomplete gamma function of halfinteger order
!
      real(8), external :: derfc
      !real(8) :: pi
      if(in.eq.0)then
        !pi = 4.0d0*datan(1.0d0)
        gmh=dsqrt(pi)*derfc(sqx)
      else
        gmh = (x**(in-1))*sqx*etx+(n-1.0d0)*&
     &       gammaincc(in-1,x,n-1.0d0)
      endif
      end function gammaincc_hin


      function incgam(n,x)
! The function works as an interface, it checks whether the index
! is integer or half integer and calls the corresponding procedure
! within the module
      implicit none
      real(8) :: x     ! Argument of the gamma function
      real(8) :: n     ! order of the gamma function
      real(8) :: incgam ! value of the gamma function (result)

      integer(4) :: in ! integer part of n

      etx = dexp(-1.0d+0*x)
      in=idint(n)
      if(n-dble(in).le.0.25)then
        incgam=gammaincc(in,x)
      else
        sqx=dsqrt(x)
        incgam=gammaincc(in,x,n)
      endif
      end function incgam

      subroutine calctildeg(lmax)
!     Calculates $\tilde{g}_{lm,l'm'}$ according to equation \ref{tildea},
!     \begin{equation}\label{calctilg}
!       \tilde{g}_{lm,l'm'}=\sqrt{4\pi}(-1)^{l}\sqrt{\frac{(l+l'+m+m')!(l+l'-m-m')!}%
!        {(2l+1)(2l'+1)[2(l+l')+1]!(l+m)!(l-m)!(l'+m')!(l'-m')!}}
!      \end{equation}
!
!     needed for the calculation of the structure constants, for l and l' = 0
!...  lmax and stores them in memory.

      implicit none
      integer(4), intent(in) :: lmax

      integer(4) :: l1,l2,m1,m2,tsize,i

      tsize=(lmax+1)*(lmax+2)*(lmax+3)*(3*lmax+2)/12
      allocate(tilg(tsize))
      tilg = 0.d0 

      i=0
      do l1=0,lmax
        do l2=0,l1
          do m1=-l1,l1
            do m2=0,l2
              i=i+1
              tilg(i)=tildeg(l1,l2,m1,m2)
            enddo ! m2
          enddo ! m1
        enddo ! l2
      enddo ! l1

      end subroutine calctildeg

      real(8) function gettildeg(l1,l2,m1,m2)
      integer(4), intent(in) :: l1
      integer(4), intent(in) :: l2
      integer(4), intent(in) :: m1
      integer(4), intent(in) :: m2

      integer(4) :: j1,j2,mj1,mj2
      integer(4) :: index1, index2, index3, index4
      integer(4) :: par,tind
      real(8) :: fact

      par=mod(abs(l1-l2),2)

      !! set the indexes j1, j2, mj1, mj2, and the multiplication factor
      if(l1.lt.l2)then
        j1=l2
        mj1=m2
        j2=l1
        mj2=m1

        !!set the multiplication factor to (-1)^(l1-l2)
        fact=-2.0d0*par+1.0d0
      else
        j1=l1
        mj1=m1
        j2=l2
        mj2=m2
        fact=1.0d0
      endif

      if(mj2.lt.0)then
        mj1=-mj1
        mj2=-mj2
      endif

      !! get the position of the value searched in the vector tilg

      index1=mj2+j2+1
      index2=(j2+1)*(mj1+j1)
      index3=j2*(j2+1)*j1+j2*(j2-1)/2
      index4=j1*(j1+1)*(j1+2)*(3*j1-1)/12
      tind=index1+index2+index3+index4

      gettildeg=fact*tilg(tind)

      end function gettildeg

      real(8) function tildeg(l1,l2,m1,m2)
!     This subroutine calculates the coefficients
!     $\tilde{g}_{lm,l'm'}$ according to equation \ref{tildea}. The factor
!     $(4\pi)^{\frac{3}{2}}$ is included into $\tilde{g}_{lm,l'm'}$.

      implicit none
      integer(4) :: l1,l2,m1,m2

      integer(4), dimension(4) :: jm        ! value of l1+m1
      integer(4) :: tjpo    ! value of 2(l1+l2)+1
      real(8) :: combjpm    ! value of (l1+l2+m1+m2)!/(l1+m1)!(l2+m2)!
      real(8) :: combjmm    ! value of (l1+l2-m1-m2)!/(l1-m1)!(l2-m2)!
      real(8) :: denom
      !real(8), parameter :: pi = 3.1415926536d0

! !EXTERNAL ROUTINES:
      real(8), external :: factr
      real(8), external :: combin
!EOP
!BOC
      jm(1)=l1+m1
      jm(2)=l2+m2
      jm(3)=l1-m1
      jm(4)=l2-m2
      tjpo=2*(l1+l2)+1
      combjpm=combin(jm(1)+jm(2),jm(1),jm(2))
      combjmm=combin(jm(3)+jm(4),jm(3),jm(4))
      denom=dble((2*l1+1)*(2*l2+1)*tjpo)
      tildeg = ((-1.0d0)**l1)*8.0d0*pi*sqrt(pi*combjmm*combjpm/denom)
      end function tildeg

      end module 
!EOP
