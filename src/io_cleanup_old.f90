      subroutine io_cleanup
      use bands,      only: nspin 
      use eigenvec,   only: vfunit
      use task,       only: iop_scratch,fid_outdbg,fid_outmb,&
     &                   fid_outmom,fid_outkpt,fid_outqp,&
     &                   savdir
      use modmpi 
      
      implicit none 
      integer :: isp 

      write(6,*) "delete scratch files"

      if(iop_scratch.eq.1.or.(iop_scratch.eq.2.and.myrank.eq.0)) then
        close(vfunit,status='delete')
      endif
      close(fid_outmb)
      close(fid_outdbg)
      close(fid_outmom) 
      close(fid_outkpt)
      close(fid_outqp)

#ifdef MPI
      write(6,*) "Process #",myrank," finished!"
      wtime2=MPI_WTIME()
      write(6,11) wtime2-wtime1
11    format('Total wall time',40x,f16.2,' seconds')
#endif 
      close(6) 
      end subroutine
