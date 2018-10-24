        subroutine scgw_herm_sxc(iop,sigm,bande,n1,n2,omg,nomg,sxc)
!
! This subroutine calculate an effective static and Hermitian Sxc matrix
! assuming that the Sxc_mn(i \omega) have been calculated 
! 
        use selfenergy, only: iop_ac, npar_ac 
        implicit none
        integer,    intent(in) :: iop  !! 
        integer,    intent(in) :: n1, n2
        integer,    intent(in) :: nomg
        real(8),    intent(in) :: bande(n1:n2) 
        real(8),    intent(in) :: omg(nomg) 
        complex(8), intent(in) :: sigm(n1:n2,n1:n2,0:nomg) 
      
        complex(8), intent(out):: sxc(n1:n2,n1:n2)

        integer:: ie1, ie2, iom
        real(8):: enk,emk
        complex(8):: ein,snn,dsig,smn,smn_m,smn_n,snm_m,snm_n
        complex(8):: apar(npar_ac),sc_ac(nomg)
        complex(8),allocatable:: tmat(:,:),sc(:,:)
        real(8) :: omg_ac(nomg)
!
!       Bare exchange
!
        do ie2=n1, n2 
          sxc(ie2,ie2)=real(sigm(ie2,ie2,0))
          do ie1 = n1, ie2-1
            sxc(ie1,ie2)=0.5d0*(sigm(ie1,ie2,0)+conjg(sigm(ie2,ie1,0)))
            sxc(ie2,ie1)=conjg(sxc(ie1,ie2))
          enddo
        enddo
        if(iop.eq.1) return
!
!       COHSEX correlation selfenergy
!
        if(iop.eq.2) then
          do ie2 = n1, n2
            sxc(ie2,ie2)=sxc(ie2,ie2)+real(sigm(ie2,ie2,1))
            do ie1=n1,ie2-1
              smn=0.5d0*(sigm(ie1,ie2,1)+conjg(sigm(ie2,ie1,1)))
              sxc(ie1,ie2)=sxc(ie1,ie2)+smn
              sxc(ie2,ie1)=sxc(ie2,ie1)+conjg(smn)
            enddo
          enddo
          return
        endif
!
!       GW correlation selfenergy in the FvSK scheme
!
        allocate(sc(n1:n2,n1:n2), tmat(n1:n2,n1:n2))

        do ie2=n1,n2
          enk = bande(ie2)
          if(enk.ge. 0.d0) then
            omg_ac = omg
            sc_ac=sigm(ie2,ie2,1:nomg)
          else
            omg_ac=-omg
            sc_ac=conjg(sigm(ie2,ie2,1:nomg))
          endif

          ein=cmplx(enk,0.d0)
          call calcacfreq(0,iop_ac,nomg,omg_ac,sc_ac,npar_ac, &
     &                apar,ein,snn,dsig)
          sc(ie2,ie2)=real(snn)

          do ie1=n1,ie2-1
            emk=bande(ie1)

            !! calculate Sc_mn at \omega = enk ==> smn_n
            if(enk.ge.0.d0) then
              omg_ac=omg
              sc_ac=sigm(ie1,ie2,1:nomg)
            else
              omg_ac= - omg
              sc_ac=conjg(sigm(ie2,ie1,1:nomg))
            endif
            ein=cmplx(enk,0.d0)
            call calcacfreq(0,iop_ac,nomg,omg_ac,sc_ac,npar_ac, &
     &                apar,ein,smn_n,dsig)

            !! calculate Sc_mn at \omega = emk ==> smn_m
            if(emk.ge.0.d0) then
              omg_ac=omg
              sc_ac=sigm(ie1,ie2,1:nomg)
            else
              omg_ac= - omg
              sc_ac=conjg(sigm(ie2,ie1,1:nomg))
            endif
            ein=cmplx(emk,0.d0)
            call calcacfreq(0,iop_ac,nomg,omg_ac,sc_ac,npar_ac, &
     &                apar,ein,smn_m,dsig)

            !! calculate Sc_nm at \omega = enk ==> snm_n
            if(enk.ge.0.d0) then
              omg_ac=omg
              sc_ac=sigm(ie2,ie1,1:nomg)
            else
              omg_ac= - omg
              sc_ac=conjg(sigm(ie1,ie2,1:nomg))
            endif
            ein=cmplx(enk,0.d0)
            call calcacfreq(0,iop_ac,nomg,omg_ac,sc_ac,npar_ac, &
     &                apar,ein,snm_n,dsig)

            !! calculate Sc_nm at \omega = emk ==> snm_m
            if(emk.ge.0.d0) then
              omg_ac=omg
              sc_ac=sigm(ie2,ie1,1:nomg)
            else
              omg_ac= - omg
              sc_ac=conjg(sigm(ie1,ie2,1:nomg))
            endif
            ein=cmplx(emk,0.d0)
            call calcacfreq(0,iop_ac,nomg,omg_ac,sc_ac,npar_ac, &
     &                apar,ein,snm_m,dsig)

            sc(ie1,ie2)=0.25*(smn_m+conjg(snm_m)+smn_n+conjg(snm_n))
            sc(ie2,ie1)=conjg(sc(ie1,ie2))
          enddo
        enddo

        sxc=sxc+sc
        deallocate(sc, tmat) 

        end subroutine

