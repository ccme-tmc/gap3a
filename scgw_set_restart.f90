
!----------------------------------------------------------------------!
!                 set up restart mode                                  !
!----------------------------------------------------------------------!
        subroutine scgw_set_restart(io,isc0)
        use bands,  only: ibgw, nbgw, bande,nbmax,nspin,nbandsgw
        use kpoints, only:nirkp
        use modmpi 
        use selfenergy, only: qpwf_coef,sigm,sigm_f, &
     &                        vhmn, vh0mn, vxcmn, &
     &                        iop_dvh
        use task, only: casename,savdir
        implicit none 
        character,intent(in):: io
        integer:: isc0
        integer:: ib0,ib1,ierr
        character(len=120) :: fn_restart
        character(len=20):: sname="scgw_set_restart"
       
        fn_restart=trim(savdir)//trim(casename)//".restart"

        !!
        !! write data to the restart file
        !!
        if(io.eq.'w'.or.io.eq.'W') then 
          if(myrank.eq.0) then    
            open(99,file=trim(fn_restart),action='write',  &
     &        form='unformatted',iostat=ierr)
            call errmsg(ierr.ne.0,sname,"Fail to open "//trim(fn_restart))

            write(99) isc0+1,ibgw,nbgw 
            write(99) bande,vxcmn,qpwf_coef,sigm,sigm_f
            if(iop_dvh.gt.0) write(99) vhmn,vh0mn 
            close(99)
          endif !! myrank.eq.0

        else 
          !!
          !! Read data from restart file 
          !!
          if(myrank.eq.0) then 
            open(99,file=trim(fn_restart),action='read', &
     &        form='unformatted',iostat=ierr)
            call errmsg(ierr.ne.0,sname,"Fail to open "//trim(fn_restart))

            read(99) isc0,ib0,ib1
            if(ib0.ne.ibgw.or.ib1.ne.nbgw) then 
              write(6,*) " WARNING: inconsistent data found" 
              write(6,*) " -- can't use the restart mode"
              close(99)
              isc0=0
              return
            endif 

            read(99) bande,vxcmn,qpwf_coef,sigm,sigm_f
            if(iop_dvh.gt.0) read(99) vhmn,vh0mn
            close(99)
          endif 
#ifdef MPI
          call mpi_bcast(isc0,1,mpi_integer,0,mpi_comm_world,ierr)
          call mpi_bcast(bande,nbmax*nirkp*nspin,mpi_double_precision,&
     &                     0,mpi_comm_world,ierr)
          call mpi_bcast(qpwf_coef,nbandsgw*nbandsgw*nirkp*nspin, &
     &                     mpi_double_complex,0,mpi_comm_world,ierr)
#endif
        endif
      end subroutine 

