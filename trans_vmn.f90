!
! This subroutine transforms the vxc matrix to the new vectors
!   vmat is the original matrix, tmat is the transform matrix 
!  flag indicates which kind of transform to make 
!    flag 
!      = 'n'/'N': 
!        vmat is Hermitian 
!        vmat^{out} _{m n} = \sum_{\mu,\nv} [ C_{\mu m} ]^* vmat^{in}_{\mu \nu} C_{\nu n}
!      = 'g'/'G': 
!         same as 'n' except that vmat is a general matrix, and 
!      = 'c': 

        subroutine trans_vmn(flag,vmat,tmat,nb)
        implicit none
        character, intent(in):: flag
        integer, intent(in):: nb 
        complex(8), intent(in):: tmat(nb,nb) 
        complex(8), intent(inout):: vmat(nb,nb)
        
        integer:: ierr
        complex(8), allocatable:: v(:,:),vc(:,:)
        complex(8):: czero,cone
        character(10)::sname="trans_vmat"

        czero=0.d0
        cone=1.d0
        allocate(v(nb,nb),vc(nb,nb) ,stat=ierr)
        call errmsg(ierr.ne.0,sname,"Fail to alloc v,vc ")
        v=vmat
        if(flag.eq.'n'.or.flag.eq.'N') then 
          call zhemm('l','l',nb,nb,cone,v,nb,tmat,nb,czero,vc,nb)
          call zgemm('c','n',nb,nb,nb,cone,tmat,nb,vc,nb,czero,vmat,nb)
        elseif(flag.eq.'c'.or.flag.eq.'C') then 
          call zhemm('r','l',nb,nb,cone,v,nb,tmat,nb,czero,vc,nb)
          call zgemm('n','c',nb,nb,nb,cone,vc,nb,tmat,nb,czero,vmat,nb)
        elseif(flag.eq.'g'.or.flag.eq.'G')then 
          call zgemm('n','n',nb,nb,nb,cone,v,nb,tmat,nb,czero,vc,nb)
          call zgemm('c','n',nb,nb,nb,cone,tmat,nb,vc,nb,czero,vmat,nb) 
        endif 

        deallocate(v,vc)

        end subroutine

