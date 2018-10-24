!BOP
!
! !ROUTINE: w2k_writevector
!
! !INTERFACE:
      subroutine w2k_writevector(gwtag,isym) 
      
! !DESCRIPTION:
!
! This subroutine writes QP energies and wavefunctions currently stored in the array bande and the 
!  direct access vectord file 
!  to an external file in the wien2k vector file format 
! !USES:
     
      use bands,       only: bande,nspin,nbmax
      use eigenvec,    only: lcmplx,zzk,kzz
      use kpoints,     only: nirkp,idikp,wkir,nkp,get_kvec
      use lapwlo,      only: lmax,lomax,nlomax,elapw,elo
      use recipvec,    only: ngkir
      use struk,       only: nat
      use task,        only: casename,spflag
      
! !LOCAL VARIABLES:

      implicit none
      character(len=*):: gwtag 
      integer,intent(in):: isym

      integer(4):: isp            ! index for spin
      integer(4):: iat            ! index for non-equivalent atoms 
      integer(4):: ib,itmp             ! index for energy band  
      integer(4):: ig             ! index for G-vector 
      integer(4):: ik,irk,iik,nk      ! index for k-points 
      integer(4):: l,k
      integer(4):: nwf
      integer(4):: ierr

      integer(4):: fout_v, fout_e            ! id for output file  
      real(8):: kvec(3) 
      character(80)::fn_en,fn_wf     ! file names for GW energy and vectors 

      character(10) :: kname=""
      character(3)  :: ipgr="" 
      character(20) :: sname="w2k_writevector"
      character(10) :: tag 

      
! !REVISION HISTORY:
!
! Created 21.06.2007 by Hong Jiang
!
!EOP
!BOC            

      if(trim(gwtag).eq.'') then 
        tag=''
      else
        tag="_"//gwtag
      endif 

      do isp=1,nspin 
!
!  set names of the new vector and energy file 
!
        fn_en=trim(casename)//".energy"//trim(spflag(isp))//trim(tag)
        fn_wf=trim(casename)//".vector"//trim(spflag(isp))//trim(tag)
        fout_e=300
        fout_v=301
!
!  write energy file
!
        open(unit=fout_e,file=trim(fn_en),action='write',iostat=ierr)
        call errmsg(ierr.ne.0,sname,"Fail to open "//trim(fn_en))

        open(unit=fout_v,file=trim(fn_wf),action='write',iostat=ierr,&
     &       form='unformatted')
        call errmsg(ierr.ne.0,sname,"Fail to open "//trim(fn_wf))

        write(6,*) "w2k_writevector: write new energy file"
        do iat=1,nat
          write(fout_e,'(100f9.5)') elapw(0:lmax,iat,isp)
          write(fout_e,'(100f9.5)') elo(0:lomax,1:nlomax,iat,isp)

          write(fout_v)  elapw(0:lmax,iat,isp)
          write(fout_v)  elo(0:lomax,1:nlomax,iat,isp)
        enddo 

        if(isym.eq.0) then 
          nk = nkp
        else
          nk = nirkp 
        endif 

        do iik=1,nk
          call get_kvec(isym,iik,ik,irk,kvec) 
          nwf = ngkir(irk)

          nwf=ngkir(irk)
          write(fout_e,100) kvec,kname,nwf,nbmax,dble(wkir(irk))
          write(fout_v) kvec,kname,nwf,nbmax,dble(wkir(irk)),ipgr

          write(fout_v) (kzz(1:3,ig,irk),ig=1,nwf)
          ik=idikp(irk)

          call readvector(ik,1,isp,0)

          do ib=1,nbmax

            write(fout_e,*) ib,bande(ib,irk,isp)*2.d0
            write(fout_v) ib,bande(ib,irk,isp)*2.d0 

            if(lcmplx) then
              write(fout_v) (zzk(ig,ib),ig=1,nwf)
            else
              write(fout_v) (real(zzk(ig,ib)),ig=1,nwf)  
            endif 
          enddo   
        enddo  
        close(fout_v) 
        close(fout_e) 
      enddo ! isp
      return

  100 format(3e19.12,a10,2i6,f5.1)
  101 format(a900)
  102 format(a)
  99  format("ERROR in w2k_writevector -- ", a) 

      end subroutine w2k_writevector
!EOC
