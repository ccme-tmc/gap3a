      subroutine calcplasmon(wpl,isym)
!
!     This subroutine calculates the plasmon frequency for metallic
!     systems 
!
      use bands,     only:fspin,nomax,numin,nomaxs,numins,nspin,nbmaxpol
      use bzinteg,   only: kwfer
      use constants, only: pi,hev
      use kpoints,   only: nirkp,nkp,idikp,kpirind,wkir
      use mommat,    only: mmatvv,init_mommat,end_mommat
      use struk,     only: vi
      use task,      only: fid_outmom
      implicit none 
      integer, intent(in) :: isym
      real(8), intent(out):: wpl 
      integer :: isp, ik, irk, ie,nomx,numn
      real(8):: coef,pnk2,wk
      complex(8):: pnk(3)
      integer:: fout=fid_outmom

      coef= fspin*4.0d0*pi*vi
      wpl = 0.d0 

      call init_mommat(numin,nomax,numin,nomax,nirkp,nspin)
      call calcmommat(0,numin,nomax,numin,nomax,1) 

      do isp=1,nspin 
        nomx = nomaxs(isp) 
        numn = numins(isp) 
        write(fout,'(3a4,5a12)') "ie","irk","ik","wk","pnk_x","pnk_y",&
     &                             "pnk_z","pnk^2"
        do irk=1,nirkp
          ik = idikp(irk)
          do ie = numn, nomx
            wk=kwfer(ie,irk,isp)*wkir(irk)
            if(wk.lt.1.e-6) cycle 
            pnk=mmatvv(:,ie,ie,irk,isp)
            pnk2 = sum(abs(pnk)**2)/3.d0
            wpl = wpl + coef*wk*pnk2
            write(fout,'(3i4,5f12.6)') ie,irk,ik,wk,real(pnk),pnk2
          enddo
        enddo
      enddo
      call end_mommat 
      wpl = sqrt(wpl)
      end subroutine 
