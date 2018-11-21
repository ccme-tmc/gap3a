      subroutine io_cleanup
      use bands,      only: nspin 
      use eigenvec,   only: vfunit
      use task,       only: iop_scratch,fid_outdbg,fid_outmb,&
     &                   fid_outmom,fid_outkpt,fid_outqp,&
     &                   savdir, fid_outgw, fid_aniso
      use modmpi 
      
      implicit none 
      integer :: isp 

      write(fid_outgw,*) "delete scratch files"

      if(iop_scratch.eq.1.or.(iop_scratch.eq.2.and.myrank.eq.0)) then
        close(vfunit,status='delete')
      endif
      close(fid_outmb)
      close(fid_outdbg)
      close(fid_outmom) 
      close(fid_outkpt)
      close(fid_outqp)
      close(fid_aniso)

#ifdef MPI
      write(fid_outgw,*) "Process #",myrank," finished!"
      wtime2=MPI_WTIME()
      write(fid_outgw,11) wtime2-wtime1
11    format('Total wall time',40x,f16.2,' seconds')
#endif 
      close(fid_outgw) 
      end subroutine
