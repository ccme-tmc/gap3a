!BOP
! !MODULE:  modmpi
! !DESCRIPTION:
!   MPI variables and interface functions
!   Can be compiled whithout mpi libraries installed
! 
!
!#include "/usr/local/include/mpif.h"

module  modmpi
  implicit none
#ifdef MPI
include  'mpif.h'
#endif
  integer,parameter::nproc_max=1000

  character(len=5):: procflag

  integer :: myrank=0           !! the rank in MPI_COMM_WORLD
  integer :: myrank_row=0       !! the rank in the row communicator 
  integer :: myrank_col=0       !! the rank in the col communicator 
  integer :: myrank_3rd=0       !! the rank in the 3rd communicator  
  integer :: myrank_ra3=0       !! the rank in the "row"+"3rd" communicator 
                                !! myrank_row should be understood as the rank in terms of "row", 
                                !! NOT the the rank in the "row". 

  integer :: mycomm
  integer :: mycomm_row         !! the communicator for "row" processes  
  integer :: mycomm_col         !! the communicator for "col" processes 
  integer :: mycomm_3rd         !! the communicator for "3rd" processes 
  integer :: mycomm_ra3         !! the communicator for "row"+"3rd" processes 

  integer :: nproc=1      !! total number of processes
  integer :: nproc_col=1, nproc_row=1, nproc_3rd=1, nproc_ra3=1

  integer :: iom_first,iom_last,iom_cnts(0:nproc_max),iom_dspl(0:nproc_max)
 
  integer,private::ierr
  real(8) :: wtime1,wtime2
  logical, private:: lsetgroup=.false.

#ifdef MPI
  interface mpi_sum_array
    module procedure mpi_sum5c
    module procedure mpi_sum4c
    module procedure mpi_sum3c
    module procedure mpi_sum2c
    module procedure mpi_sum1c
    module procedure mpi_sum1r
  end interface 

  interface mpi_gather_array
    module procedure mpi_gather4c
    module procedure mpi_gather3c
    module procedure mpi_gather2c
    module procedure mpi_gather1c
  end interface

  interface mpi_sum_scalar
    module procedure mpi_sum_real
    module procedure mpi_sum_cmplx
  end interface

#endif

contains

  subroutine init_mpi
    character(len=10),external::Int2Str
#ifdef MPI
    !        mpi init
    integer::namelen
    character(len=MPI_MAX_PROCESSOR_NAME)::processor_name

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    call MPI_GET_PROCESSOR_NAME(processor_name,namelen,ierr);
    print *, "Process ", myrank, " of ", nproc, " running on ",     &
     &     trim(processor_name)

    wtime1=MPI_WTIME()

#endif
#ifndef MPI
    nproc=1
    nproc_row=1
    nproc_col=1
    nproc_3rd=1
    myrank=0
    myrank_col=0
    myrank_row=0
    myrank_3rd=0
#endif

    procflag=''
    if(myrank.ne.0) then
      procflag='-p'//trim(int2str(myrank))
    endif

  end subroutine init_mpi

  subroutine end_mpi
#ifdef MPI
!    call MPI_barrier(MPI_COMM_WORLD,ierr)
!    call mpi_free_grp
    call MPI_FINALIZE(ierr)
#endif
  end subroutine end_mpi


  subroutine mpi_set_range(np,ip,num,i0,istart,iend,icnts,idisp)
!! This is a general subroutine that seperate i0..i0+num-1 equally to np processes 
!! in case mod(num, nranks) .eq.0, the residual is put to last mod(num, np) processes 
!! This subroutine can be used to set k-point or frequency point range  
    integer,intent(in)::np,ip,num,i0
    integer,intent(out)::istart,iend 
    integer,intent(out),optional::icnts(0:),idisp(0:)
    integer :: icount(nproc_max)
    integer :: i

!
! Since the calculation for Gamma point, iq=1, at most takes twice cpu time than
! other q-points, we take the folllowing strategy:
!    if nkpt is equal to an integral times nproc, then divide nkpt equally
!    otherwise, the residual k-points are equally assigned to non-root processes
!    starting from the last one
!
    if(np.le.1) then 
      istart=i0
      iend=i0+num-1 
      icount(1)=num
    else 
      icount=num/np
      do i=1,mod(num,np)
        icount(np-i+1)=icount(np-i+1)+1
      enddo
      istart=i0
      do i=1,ip
        istart=istart+icount(i)
      end do
      iend=istart+icount(ip+1)-1
    endif 

    if(present(icnts)) then 
      icnts(0:np-1)=icount(1:np)
      if(present(idisp)) then 
        idisp(0)=0
        do i=1,np-1
          idisp(i)=idisp(i-1)+icnts(i-1)
        enddo
      endif 
    endif 
    return
  end subroutine 

#ifdef MPI
  subroutine mpi_set_group(nq,nk) 
! this subroutine is to set up the MPI commonicators 
! all processes are divided into groups
    use task, only: fid_outgw
    implicit none
    integer, intent(in):: nq,nk
    integer myrow,mycol,my3rd,myra3
    character(len=10)::str
    integer :: nproc_col_old, nproc_row_old, nproc_3rd_old, nproc_ra3_old
    integer :: icolor 

    if(nproc.eq.1) then 
      nproc_row=1
      nproc_col=1
      nproc_3rd=1
      myrank_col=0
      myrank_row=0
      myrank_3rd=0
      myrank_ra3=0
      return 
    endif 

    write(fid_outgw,*) " Setup MPI group:"
    write(fid_outgw,*) " nq=",nq," nk=",nk

    if(lsetgroup) call mpi_free_group

    !! set the default nproc_col 
    if(nproc_col.eq.0 .or. mod(nproc,nproc_col).ne.0 &
   &    .or. mod(nq,nproc_col).ne.0) then 

      call max_common_div( nq, nproc, nproc_col) 

    endif 

    !! set the default nproc_row
   
    if(nproc_row.eq.0.or. mod(nproc/nproc_col,nproc_row).ne.0 &
   &  .or.mod(nk,nproc_row).ne.0 ) then 
      call max_common_div(nproc/nproc_col, nk, nproc_row)
    endif 

    !! set nproc_3rd
    nproc_3rd = nproc/(nproc_col*nproc_row) 
    nproc_ra3 = nproc_row*nproc_3rd 
    write(fid_outgw,*) "nproc_col/row/3rd/ra3(init) =",nproc_col,nproc_row,nproc_3rd,nproc_ra3

    my3rd = myrank/(nproc_col*nproc_row)
    mycol = (myrank - my3rd*nproc_col*nproc_row)/nproc_row 
    myrow = myrank - my3rd*nproc_col*nproc_row - mycol*nproc_row 

    !! create the communicators corresponding to row, col and 3rd
    
    !! mycomm_row : processes with identical mycol and my3rd but different myrow 
    !! mycomm_col : processes with identical myrow and my3rd but different mycol
    !! mycomm_3rd : processes with identical mycol and myrow but different my3rd 

    icolor = mycol*nproc + my3rd*nproc**2
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,icolor,myrow,mycomm_row,ierr)

    icolor = myrow + my3rd*nproc**2
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,icolor,mycol,mycomm_col,ierr)

    icolor = myrow + mycol*nproc
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,icolor,my3rd,mycomm_3rd,ierr)

    icolor = mycol 
    myra3 = my3rd*nproc_row + myrow 
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,icolor,myra3,mycomm_ra3,ierr)

    call MPI_COMM_RANK(mycomm_row,myrank_row,ierr) 
    call MPI_COMM_RANK(mycomm_col,myrank_col,ierr)
    call MPI_COMM_RANK(mycomm_3rd,myrank_3rd,ierr)
    call MPI_COMM_RANK(mycomm_ra3,myrank_ra3,ierr)
    write(fid_outgw,'(a40,5i5)') "myrank, mycol/row/3rd/ra3=",myrank,mycol,myrow,my3rd,myra3
    write(fid_outgw,'(a45,4i5)') "   myrank_col/row/3rd/ra3=",myrank_col,myrank_row,myrank_3rd,myrank_ra3
    write(fid_outgw,*) "mycomm_col/row/3rd=",mycomm_col,mycomm_row,mycomm_3rd 

    if(     mycol.ne.myrank_col &
   &   .or. myrow.ne.myrank_row &
   &   .or. my3rd.ne.myrank_3rd  & 
   &   .or. myra3.ne.myrank_ra3) then 
      write(fid_outgw,*) "WARNING: something may be wrong here!"
    endif  

!   reset nproc_row/row/3rd, since this may not be identical for each column 
    nproc_row_old = nproc_row
    nproc_col_old = nproc_col
    nproc_3rd_old = nproc_3rd 
    nproc_ra3_old = nproc_ra3
    call MPI_COMM_SIZE( mycomm_col, nproc_col, ierr)
    call MPI_COMM_SIZE( mycomm_row, nproc_row, ierr)
    call MPI_COMM_SIZE( mycomm_3rd, nproc_3rd, ierr)
    call MPI_COMM_SIZE( mycomm_ra3, nproc_ra3, ierr)

    !! some consistency check 
    if(nproc_row_old*nproc_col_old*nproc_3rd_old.eq.nproc &
   &   .and. (     nproc_col .ne. nproc_col_old         &
   &           .or.nproc_row .ne. nproc_row_old         &
   &           .or.nproc_3rd .ne. nproc_3rd_old ) ) then 
      write(fid_outgw,*) "ERROR: something wrong when setting MPI groups:"
      if(myrank.eq.0) then 
        write(fid_outgw,*) "  new nproc_col/row/3rd/ra3=",nproc_col,nproc_row,nproc_3rd,nproc_ra3
      endif 
      stop  
    endif 

    mycomm=MPI_COMM_WORLD
    write(fid_outgw,*) "nproc_col/row/3rd/ra3(final) =",nproc_col,nproc_row,nproc_3rd,nproc_ra3 

    lsetgroup = .true.

  endsubroutine 

  subroutine mpi_free_group
    call MPI_COMM_FREE(mycomm_row,ierr)
    call MPI_COMM_FREE(mycomm_col,ierr)
    call MPI_COMM_FREE(mycomm_3rd,ierr) 
    call MPI_COMM_FREE(mycomm_ra3,ierr) 
    lsetgroup = .false.
  end subroutine 

!
! This subroutine  is used to sum complex arrays of the form, a(1:n1,1:n2,1:n3) that is calculated at different 
! processes
!      iflag=0 --- reduce to the root process
!              --- each process 
!
  subroutine mpi_sum3c(iflag,a,n1,n2,n3,comm0)
    implicit none
    integer,intent(in)::  n1,n2,n3
    integer,intent(in):: iflag
    integer,intent(in):: comm0
    complex(8),intent(inout):: a(n1,n2,n3)

    integer :: ierr
    integer :: rank,comm,ndata

    comm=comm0
    if(comm.eq.0) then
      comm=MPI_COMM_WORLD
    endif

    call mpi_comm_rank( comm, rank, ierr)
    call MPI_Barrier(comm,ierr)

    if(iflag.eq.0) then
      if(rank.eq.0) then 
        call MPI_Reduce(MPI_IN_PLACE,a,n1*n2*n3,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm,ierr)
      else 
        call MPI_Reduce(a,0,n1*n2*n3,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm,ierr)
      endif 
    else
      call MPI_AllReduce(MPI_IN_PLACE,a,n1*n2*n3,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    endif
  end subroutine


  subroutine mpi_sum5c_new(iflag,a,n1,n2,n3,n4,n5,comm0)
    implicit none
    integer n1,n2,n3,n4,n5
    integer,intent(in):: iflag
    integer,intent(in):: comm0
    complex(8),intent(inout):: a(n1,n2,n3,n4,n5)

    integer :: i3,i4,i5 

    do i5=1,n5
      do i4=1,n4
        do i3=1,n3 
           call mpi_sum2c(iflag,a(:,:,i3,i4,i5),n1,n2,comm0)
        enddo 
      enddo 
    enddo 
  end subroutine 


  subroutine mpi_sum5c(iflag,a,n1,n2,n3,n4,n5,comm0)
    implicit none
    integer n1,n2,n3,n4,n5
    integer,intent(in):: iflag
    integer,intent(in):: comm0
    complex(8),intent(inout):: a(n1,n2,n3,n4,n5)

    integer :: ierr
    integer :: rank,comm,ndata

    comm=comm0
    if(comm.eq.0) then
      comm=MPI_COMM_WORLD
    endif
    ndata=n1*n2*n3*n4*n5

    call mpi_comm_rank( comm, rank, ierr)
    call MPI_Barrier(comm,ierr)

    if(iflag.eq.0) then
      if(rank.eq.0) then
        call MPI_Reduce(MPI_IN_PLACE,a,ndata,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm,ierr)
      else
        call MPI_Reduce(a,0,ndata,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm,ierr)
      endif
    else
      call MPI_AllReduce(MPI_IN_PLACE,a,ndata,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    endif
  end subroutine

  subroutine mpi_sum4c(iflag,a,n1,n2,n3,n4,comm0)
    implicit none
    integer n1,n2,n3,n4
    integer,intent(in):: iflag
    integer,intent(in):: comm0
    complex(8),intent(inout):: a(n1,n2,n3,n4)

    integer :: ierr
    integer :: rank,comm,ndata

    comm=comm0
    if(comm.eq.0) then
      comm=MPI_COMM_WORLD
    endif
    ndata=n1*n2*n3*n4

    call mpi_comm_rank( comm, rank, ierr)
    call MPI_Barrier(comm,ierr)

    if(iflag.eq.0) then
      if(rank.eq.0) then 
        call MPI_Reduce(MPI_IN_PLACE,a,ndata,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm,ierr)
      else 
        call MPI_Reduce(a,0,ndata,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm,ierr)
      endif 
    else
      call MPI_AllReduce(MPI_IN_PLACE,a,ndata,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    endif
  end subroutine

  subroutine mpi_sum2c(iflag,a,n1,n2,comm0)
    implicit none
    integer n1,n2
    integer,intent(in):: iflag
    integer,intent(in)::comm0
    complex(8),intent(inout):: a(n1,n2)

    integer :: ierr
    integer :: rank,comm

    comm=comm0
    if(comm.eq.0) then
      comm=MPI_COMM_WORLD
    endif

    call mpi_comm_rank( comm, rank, ierr)
    call MPI_Barrier(comm,ierr)

    if(iflag.eq.0) then 
      if(rank.eq.0) then 
        call MPI_reduce(MPI_IN_PLACE,a,n1*n2,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm,ierr)
      else 
        call MPI_reduce(a,0,n1*n2,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm,ierr)
      endif 
    else 
      call MPI_AllReduce(MPI_IN_PLACE,a,n1*n2,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    endif 

  end subroutine

  subroutine mpi_sum1c(iflag,a,n1,comm0)
    use constants, only: czero
    implicit none
    integer,intent(in):: iflag
    integer,intent(in):: n1
    integer,intent(in)::comm0
    complex(8),intent(inout):: a(n1)

    integer :: ierr
    integer :: rank,comm

    comm=comm0
    if(comm.eq.0) then
      comm=MPI_COMM_WORLD
    endif

    call MPI_COMM_RANK(comm, rank, ierr)
    call MPI_BARRIER(comm, ierr)

    if(iflag.eq.0) then
      if(rank.eq.0) then 
        call MPI_reduce(MPI_IN_PLACE,a,n1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm,ierr)
      else 
        call MPI_reduce(a,0,n1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm,ierr)
      endif 
    else
      call MPI_AllReduce(MPI_IN_PLACE,a,n1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    endif

  end subroutine

  subroutine mpi_sum1r(iflag,a,n1,comm0)
    implicit none
    integer,intent(in):: iflag
    integer,intent(in):: n1
    integer,intent(in)::comm0
    real(8),intent(inout):: a(n1)

    integer :: ierr
    integer :: rank,comm

    comm=comm0
    if(comm.eq.0) then
      comm=MPI_COMM_WORLD
    endif

    call mpi_comm_rank( comm, rank, ierr)
    call MPI_Barrier(comm,ierr)

    if(iflag.eq.0) then
      if(rank.eq.0) then
        call MPI_reduce(MPI_IN_PLACE,a,n1,MPI_REAL8,MPI_SUM,0,comm,ierr)
      else
        call MPI_reduce(a,0,n1,MPI_REAL8,MPI_SUM,0,comm,ierr)
      endif
    else
      call MPI_AllReduce(MPI_IN_PLACE,a,n1,MPI_REAL8,MPI_SUM,comm,ierr)
    endif
  end subroutine

!
! The following subroutines are used to gather arrays calculated on different 
! processes belonging the given communicator (comm). The distribution of the data 
! is only via the last dimension of the array 
!

  subroutine mpi_gather4c(iflag,as,ns,n1,n2,n3,ar,cnts,dspl,comm)
  implicit none
  integer::iflag
  integer :: n1, n2, n3,ns
  complex(8) :: as(n1,n2,n3,ns),ar(:,:,:,:)
  integer :: ierr
  integer :: comm
  integer :: cnts(0:nproc_max),dspl(0:nproc_max)

  integer::rstype


  call MPI_Type_contiguous(n1*n2*n3,MPI_DOUBLE_COMPLEX,rstype,ierr)
  call MPI_Type_commit(rstype,ierr)

  call MPI_Barrier(comm,ierr)
  if(iflag.eq.0) then
    call MPI_GatherV(as,ns,rstype,ar,cnts,dspl,rstype,0,comm,ierr)
  else
    call MPI_AllGatherV(as,ns,rstype,ar,cnts,dspl,rstype,comm,ierr)
  endif
  call MPI_Type_free(rstype,ierr)
  end subroutine

  subroutine mpi_gather3c(iflag,as,ns,n1,n2,ar,cnts,dspl,comm)
  implicit none
  integer::iflag
  integer :: n1, n2, ns
  complex(8) :: as(n1,n2,ns),ar(:,:,:)
  integer :: ierr
  integer :: comm
  integer :: cnts(0:nproc_max),dspl(0:nproc_max)

  integer::rstype

  call MPI_Type_contiguous(n1*n2,MPI_DOUBLE_COMPLEX,rstype,ierr)
  call MPI_Type_commit(rstype,ierr)

  call MPI_Barrier(comm,ierr)
  if(iflag.eq.0) then
    call MPI_GatherV(as,ns,rstype,ar,cnts,dspl,rstype,0,comm,ierr)
  else
    call MPI_AllGatherV(as,ns,rstype,ar,cnts,dspl,rstype,comm,ierr)
  endif
  call MPI_Type_free(rstype,ierr)
  end subroutine


  subroutine mpi_gather2c(iflag,as,ns,n1,ar,cnts,dspl,comm)
  implicit none 
  integer::iflag
  integer :: n1, ns
  complex(8) :: as(n1,ns),ar(:,:) 
  integer :: ierr
  integer :: comm 
  integer :: cnts(0:nproc_max),dspl(0:nproc_max)
  
  integer::rstype

  call MPI_Type_contiguous(n1,MPI_DOUBLE_COMPLEX,rstype,ierr)
  call MPI_Type_commit(rstype,ierr)

  call MPI_Barrier(comm,ierr)
  if(iflag.eq.0) then 
    call MPI_GatherV(as,ns,rstype,ar,cnts,dspl,rstype,0,comm,ierr)
  else 
    call MPI_AllGatherV(as,ns,rstype,ar,cnts,dspl,rstype,comm,ierr)
  endif 
  call MPI_Type_free(rstype,ierr)
  end subroutine 

  subroutine mpi_gather1c(iflag,as,ns,ar,cnts,dspl,comm)
  implicit none
  integer::iflag
  integer :: ns
  complex(8) :: as(ns),ar(:)  !! array to send (as) and received (ar)
  integer :: comm
  integer :: ierr
  integer :: cnts(0:nproc_max),dspl(0:nproc_max)

  integer::rstype

  call MPI_Barrier(comm,ierr)
  rstype=MPI_DOUBLE_COMPLEX
  if(iflag.eq.0) then
    call MPI_GatherV(as,ns,rstype,ar,cnts,dspl,rstype,0,comm,ierr)
  else
    call MPI_AllGatherV(as,ns,rstype,ar,cnts,dspl,rstype,comm,ierr)
  endif
  call MPI_Type_free(rstype,ierr)
  end subroutine


  subroutine mpi_sum_real(iflag,s,comm)
    implicit none
    integer,intent(in):: iflag
    integer::comm
    real(8),intent(inout):: s

    real(8)::sbuf(1)

    sbuf(1)=s
    call mpi_sum1r(iflag,sbuf,1,comm)
    s=sbuf(1)
  end subroutine

  subroutine mpi_sum_cmplx(iflag,s,comm)
    use constants, only: czero
    implicit none
    integer,intent(in):: iflag
    integer::comm
    complex(8),intent(inout):: s

    integer :: ierr
    complex(8)::sbuf(1)

    sbuf(1)=s
    call mpi_sum1c(iflag,sbuf,1,comm)
    s=sbuf(1)
  end subroutine

  subroutine max_common_div(n1,n2,d) 
    implicit none 
    integer, intent(in):: n1,n2
    integer, intent(out):: d
    integer:: i
    
    do i = min(n1,n2), 1, -1
      if (mod(n1,i).eq.0.and.mod(n2,i).eq.0) then 
        d = i
        exit 
      endif 
    enddo
  end subroutine 

#endif

end module modmpi
    

