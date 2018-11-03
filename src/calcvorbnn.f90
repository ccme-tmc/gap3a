!BOP
!
! !ROUTINE: calcvorbnn
!
! !INTERFACE:
      subroutine calcvorbnn(iat,idf,irk,isp)

! !DESCRIPTION:
!
!This subroutine calculates the matrix elements $v^{orb}_{nn}(\vec{k})$ 
!
! !USES:

      use xcpot,       only: vorbnn,natorb,iatorb,nlorb,lorb,vorb,uiorb,&
     &                       dmorb
      use bands,       only: ibgw,nbgw
      use bzinteg,     only: kiw 
      use eigenvec,    only: alfa, beta, gama
      use kpoints,     only: wkir,idikp
      use task,        only: fid_outdbg
      implicit none

      integer,intent(in)::iat,idf,irk,isp

! !LOCAL VARIABLES:
      integer(4) :: ie      ! Counter: run over the eigenvalues in the corresponding k-point.
      integer(4) :: m1      ! m-quantum number of the (L)APW eigenfunction  at k
      integer(4) :: m2      ! m-quantum number of the (L)APW eigenfunction   at k' = k - q
      integer(4) :: l1m1    ! Counter: Runs over MTS (L)APW basis  functions (l1,m1)
      integer(4) :: l2m2    ! Counter: Runs over MTS (L)APW basis functions l2,m2)
      integer(4) :: ia,il,l,itmp
      
      real(8) :: udud,uu2,udu2
      complex(8) :: af,bt,gm,caf,cbt,cgm

      complex(8) :: nimm(-3:3,-3:3,ibgw:nbgw)
      complex(8) :: trvd  
     
! !REVISION HISTORY:
! 
! Created  23.10.2007  by Hong Jiang
!
!EOP
!
!BOC
      do ia=1,natorb
        if(iatorb(ia) .ne.iat) cycle   

        write(fid_outdbg,*) 
        write(fid_outdbg,*) '#vorb for atom ',iat
        write(fid_outdbg,*) 

        do il=1,nlorb(ia) 
          l=lorb(il,ia) 
          
          write(fid_outdbg,'(a,6e16.6)') 'uiorb',uiorb(1:6,il,ia,isp) 
          udud=uiorb(3,il,ia,isp)
          uu2= uiorb(4,il,ia,isp)
          udu2=uiorb(5,il,ia,isp)

          do m2=-l,l
            do m1=-l,l 
              l1m1= l*l+l+m1+1
              l2m2= l*l+l+m2+1

              do ie=ibgw,nbgw 
                af=alfa(ie,l1m1,idf)
                caf=conjg(alfa(ie,l2m2,idf))
                bt=beta(ie,l1m1,idf)
                cbt=conjg(beta(ie,l2m2,idf))
                gm=gama(ie,1,l1m1,idf)
                cgm=conjg(gama(ie,1,l2m2,idf))
            
                nimm(m1,m2,ie)= af*caf  + bt*cbt*udud + gm*cgm  &
     &                         +(af*cgm+gm*caf)*uu2             &
     &                         +(bt*cgm+gm*cbt)*udu2 
                
                dmorb(m1,m2,il,ia) = dmorb(m1,m2,il,ia)         &
     &             + wkir(irk)*kiw(ie,irk,isp)*nimm(m1,m2,ie)

                vorbnn(ie,irk,isp) = vorbnn(ie,irk,isp)         &
     &                     +real(vorb(m1,m2,il,ia,isp)*nimm(m1,m2,ie))
              enddo ! ie
            enddo   ! m1
          enddo ! m2

          do ie=ibgw,nbgw 
            write(fid_outdbg,*) 
            write(fid_outdbg,*) '--- nimm for l=',l,'and i=',ie
            write(fid_outdbg,*) 
            do m2=-l,l
              do m1=-l,l
                write(fid_outdbg,'(2E16.8)') nimm(m1,m2,ie)
              enddo
            enddo
          enddo 
        enddo  ! il
      enddo ! ia

      end subroutine calcvorbnn

!EOC      


