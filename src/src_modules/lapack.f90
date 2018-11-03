!******************** LAPACK_m *************************
!* Wrappers of some usefule lapack procedures 
module LAPACK
  implicit none
  integer,private::info
  private subErrMsg

contains 

  subroutine LA_HerEigen(A,D,N)
    integer::N
    complex(8)::A(N,N)
    real(8)::D(N)
    complex(8)::Work(2*N)
    real(8)::RWORK(3*N-2)
    integer::info,i,j

    call ZHEEV( 'V', 'L', N, A, N, D, WORK, 2*N, RWORK, INFO )
    if(info/=0 ) call subErrMSG("Failure in subHerEigen")
  end subroutine 

  subroutine LA_SymEigen(a,d,n)
    integer::N
    real(8)::a(N,N),d(N)
    real(8)::Work(N*6)
    integer::info,lwork
    character::JOBZ='V',UPLO='U'

    lwork=6*N
!    print *,'solve eigen problem by DSYEV'
    call DSYEV( JOBZ, UPLO, N, A, N, d, WORK, LWORK, INFO )
    if(info/=0) call subErrMsg("Fail to Diagonalize")
 
 
  end subroutine        
 
  subroutine LA_GTSV(x,d,du,dl,N)
    integer::N
    real(8)::x(N),d(N),du(N-1),dl(N-1)

    call DGTSV(N,1,dl,d,du,x,N,info)
    if(info/=0) call subErrMsg("Error in LAPACK_m->subGTSV")
    
  end subroutine 

  subroutine LA_SolLinEq(a,b,n)
    integer::n
    real(8)::a(n,n),b(n)
    integer::IPIV(N)
    real(8)::work(N*2)

    call DSYSV('u', N, 1, A, N, IPIV, B,N, WORK,N, INFO )
    if(info/=0) call subErrMsg("Error in LAPACK_m->subSolLinEq")
  end subroutine 

  subroutine LA_InvSym(A,N)
    integer::N
    real(8)::a(N,N)
    integer::IPIV(N),info,i,j
    real(8)::work(N)

    call DSYTRF( 'U', N,  A, N, IPIV, WORK, N, INFO )
    if(info/=0) call subErrMsg("Error in LAPACK_m->subInvSym")
    call DSYTRI( 'U', N,  A, N, IPIV, WORK, INFO )
    if(info/=0) call subErrMsg("Error in LAPACK_m->subInvSym")
    do i=1,N
      do j=1,i-1
        A(i,j)=A(j,i)
      enddo
    enddo
    
  end subroutine 

  subroutine subErrMsg(info)
  character(*):: info
  write(6,*) info
  stop 
  end subroutine 
    
end module
  
