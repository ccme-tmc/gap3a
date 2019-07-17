!BOP
!
! !Routine : crpa_output
! !INTERFACE:
      subroutine crpa_output()  
!
!! DESCRIPTION :

!  This subroutine write the main results of CRPA and extract some main
!  information  

!! USES:
      use bands,     only: nspin
      use crpa,      only: nlorb_wf,nbmax_wf,nbmin_wf,nlmorb,pln,nbk_wf,&
     &                     nsp_crpa,info_orb,vmat,umat,iop_crpa,ncell, &
     &                     wf_centers,rcell 
      use constants, only: HeV,Bohr2Ans
      use freq,      only: nomeg,omega
      use kpoints,   only: nirkp,nkp,get_kindex,kpirind
      use eigenvec,  only: lsymvector 
      use struk,     only: alat,rbas,ortho
      use task,      only: casename 
      implicit none 

!! LOCAL VARIABLES:
      character(20):: sname="crpa_output"
      integer :: icell,ial,iom,isp,i
      integer :: ierr
      integer :: fid
      character(len=120) :: fn_u4mat,fn_umat
      real(8):: U, J, Ud

      ! first write out main results of CRPA calculations: v(m1,m2,m3,m4) and u(m1,m2,m3,m4,omega)

      fn_u4mat=trim(casename)//".U4mat"
      fn_umat=trim(casename)//".Umat"

!=====================================================================!
!     write the full U-matrix                                         !
!=====================================================================!

      fid = 998
      write(6,*) "Write cRPA U-matrix in eV to "//trim(fn_umat)
      open(fid,file=trim(fn_u4mat),action='write',iostat=ierr)
      call errmsg(ierr.ne.0,sname,"Fail to open "//trim(fn_umat))
      write(fid,'(a)'  ) "# Information on Wanner functions"
      write(fid,'(2i4)') nlorb_wf, nlmorb,ncell 
      do ial=1,nlorb_wf
        write(fid,'(5i4)') ial,info_orb(:,ial)
      enddo 
      do icell=1,ncell 
        write(fid,'(a)') "# matrix elements of bare interaction"
        do isp=1,nsp_crpa
          write(fid,'(4e24.16)') vmat(:,:,isp,icell) 
        enddo 

        if(iop_crpa.ge.0) then 
          write(fid,'(a)') "# matrix elements of screened interaction"
          do iom=1,nomeg 
            do isp=1,nsp_crpa 
              write(fid,'(4e24.16)') umat(:,:,isp,iom,icell)
            enddo 
          enddo
        endif 
      enddo 
      close(fid) 

!=====================================================================!
!             write the reduced U-matrix                              !
!=====================================================================!

      open(fid,file=trim(fn_umat),action='write',iostat=ierr)
      write(6,*) "Write the reduced interaction matrix:"
      write(6,*) "   U_{m1,m2}=U(m1,m2,m1,m2)"
      write(6,*) "   J_{m1,m2}=U(m1,m2,m2,m1)" 
      write(fid,*) "# Bare Coulomb Interactions"
      do icell=1,ncell 
        if(icell.eq.1) then
          do isp=1,nsp_crpa 
            call sub_reduced_m(vmat(:,:,isp,1)) 
          enddo   
        else 
          write(fid,*) "#Inter-cell interactions with the cell ",icell
          do isp=1,nsp_crpa
            call sub_reduced_intercell(vmat(:,:,isp,icell),icell)
          enddo
        endif  
      enddo 

      if(iop_crpa.ge.0) then  
        do iom=1,nomeg
          write(fid,'(/,a,f12.6)')"# W at omega (eV)=",omega(iom)*HeV
          do icell=1,ncell 
            if(icell.eq.1) then 
              do isp=1,nsp_crpa 
                call sub_reduced_m(umat(:,:,isp,iom,1))
              enddo ! isp 
            else 
              write(fid,*) "#Inter-cell interactions with the cell ",icell
              do isp=1,nsp_crpa
                call sub_reduced_intercell(umat(:,:,isp,iom,icell),icell)
              enddo ! isp 
            endif ! icell.eq.1
          enddo  ! icell 
        enddo ! iom 
      endif  
      close(fid) 
     
      contains 
        subroutine sub_aver_UJ(umat,jmat,nm) 
        integer,intent(in):: nm
        real(8),intent(in):: umat(nm,nm), jmat(nm,nm) 

        integer:: im

        U  = sum(umat)/(nm*nm)
        J  = sum(jmat)
        Ud = 0.d0
        do im=1,nm
          Ud = Ud + umat(im,im) 
          J  = J - jmat(im,im)
        enddo
        J = J/(nm*(nm-1))
        Ud = Ud/nm 
        end subroutine 

!=====================================================================!
!      extract reduced U-matrix for multi-center/orbital cases                 !
!=====================================================================!

        subroutine sub_reduced_m(u4)
        !! U_{m1,m2}=U(m1,m2,m1,m2) = U_{ (m1 m1), (m2 m2)}
        !! J_{m1,m2}=U(m1,m2,m2,m1) = U_{ (m1 m2), (m1 m2)} 
        complex(8),intent(in):: u4(nlmorb*nlmorb,nlmorb*nlmorb)
        integer:: nm1,nm2
        integer:: im1,im2
        integer:: ill1,ill2
        integer:: il1,il2
        integer:: im1_0,im2_0,im_0 
        real(8),allocatable:: u2(:,:),j2(:,:)
        real(8):: V_av, dist, dr(3) 

        ! first get the onsite U and J matrices 
        im2_0 = 0 
        do il2=1,nlorb_wf 
          nm2 = info_orb(3,il2) 
          if(il2.gt.1) im2_0 = im2_0 + info_orb(3,il2-1) 

          im1_0 = 0
          do il1=1,il2 
            nm1 = info_orb(3,il1) 
            if(il1.gt.1) im1_0 = im1_0 + info_orb(3,il1-1)

            allocate(u2(nm1,nm2),j2(nm1,nm2)) 

            do im2 = im2_0 + 1, im2_0 + nm2 
              do im1 = im1_0 + 1, im1_0 + nm1
                ill1 = (im1-1)*nlmorb + im1
                ill2 = (im2-1)*nlmorb + im2
                u2(im1-im1_0, im2-im2_0) = real(u4(ill1,ill2))*HeV

                ill1 = (im2-1)*nlmorb + im1
                ill2 = (im2-1)*nlmorb + im1
                j2(im1-im1_0, im2-im2_0) = real(u4(ill1,ill2))*HeV
              enddo
            enddo 
            if(il1.eq.il2) then  !! on-site interactions 
              call sub_aver_UJ(u2,j2,nm1) 
              write(fid,300) il1, U, Ud, J  
              write(fid,301) "# U"
              do im1=1,nm1
                write(fid,302) u2(im1,:)
              enddo
              write(fid,301) "# J"
              do im1=1,nm1
                write(fid,302) j2(im1,:)
              enddo

            else !! inter-atomic interations 

              dr(:)=wf_centers(:,il1)-wf_centers(:,il2) 
              if(ortho) then 
                dist = sqrt(sum( (alat(:)*dr(:))**2 ))
              else 
                dist = 0.d0
                do i=1,3
                  dist=dist+sum(dr(:)*rbas(:,i))**2 
                enddo
                dist=sqrt(dist) 
              endif 
              dist = dist*Bohr2Ans

              V_av = sum(u2)/(nm1*nm2) 
              write(fid,305) il1,il2, dist, V_av
              do im1=1,nm1
                write(fid,302) u2(im1,:)
              enddo
            endif 

            deallocate(u2,j2) 
          enddo
        enddo 
 101    format(a,3f10.4) 
 300    format("#On-site interactions for il=",i3,&
     &     ' U(full),Ud(diag) and J (in eV):',3f10.4) 
 301    format(a) 
 302    format(100f12.6) 
 305    format("#Inter-atomic V for il1, il2, d12(Ans) and V_av(eV):",&
     &          2i3,2f10.3) 


        end subroutine
 
        subroutine sub_reduced_intercell(u4,icell)
        !! U_{m1,m2}=U(m1,m2,m1,m2) = U_{ (m1 m1), (m2 m2)}
        !! J_{m1,m2}=U(m1,m2,m2,m1) = U_{ (m1 m2), (m1 m2)}
        complex(8),intent(in):: u4(nlmorb*nlmorb,nlmorb*nlmorb)
        integer,intent(in):: icell 

        integer:: nm1,nm2
        integer:: im1,im2
        integer:: ill1,ill2
        integer:: il1,il2
        integer:: im1_0,im2_0,im_0 
        real(8),allocatable:: u2(:,:),j2(:,:)
        real(8):: V_av,dist,dr(3)

        im2_0 = 0 
        do il2=1,nlorb_wf 
          nm2 = info_orb(3,il2) 
          if(il2.gt.1) im2_0 = im2_0 + info_orb(3,il2-1) 

          im1_0 = 0
          do il1=1,nlorb_wf
            nm1 = info_orb(3,il1) 
            if(il1.gt.1) im1_0 = im1_0 + info_orb(3,il1-1)

            allocate(u2(nm1,nm2),j2(nm1,nm2)) 

            do im2 = im2_0 + 1, im2_0 + nm2 
              do im1 = im1_0 + 1, im1_0 + nm1
                ill1 = (im1-1)*nlmorb + im1
                ill2 = (im2-1)*nlmorb + im2
                u2(im1-im1_0, im2-im2_0) = real(u4(ill1,ill2))*HeV
              enddo
            enddo 
            V_av = sum(u2)/(nm1*nm2) 

            dr(:)=wf_centers(:,il2)-wf_centers(:,il1) &
     &        +rcell(:,icell)-rcell(:,1)  

            if(ortho) then
              dist = sqrt(sum( (alat(:)*dr(:))**2 ))
            else
              dist = 0.d0
              do i=1,3
                dist=dist+sum(dr(:)*rbas(:,i))**2
              enddo
              dist=sqrt(dist)
            endif
            dist = dist*Bohr2Ans

            write(fid,301) il1,il2, dist, V_av
            do im1=1,nm1
              write(fid,303) u2(im1,:)
            enddo

            deallocate(u2,j2) 
          enddo
        enddo 
 300    format(a) 
 301    format("#Inter-atomic V for il1, il2, d12(Ans) and V_av(eV):",&
     &          2i3,2f10.3) 
 302    format("# ",a1,"_{m1,m2} (il1=",i2,", il2=",i2,")") 
 303    format(100f12.6) 

        end subroutine
 
      end subroutine crpa_output
