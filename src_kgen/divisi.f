        subroutine divisi(idkp,nkp,idiv,klist)
        parameter (nprim=16)
        parameter (niter=10)
	dimension klist(idkp,3),iprim(nprim)
!
        idivi=1
	iprim(1)=2
	iprim(2)=3
	iprim(3)=5
	iprim(4)=7
        iprim(5)=11      
        iprim(6)=13     
        iprim(7)=17     
        iprim(8)=19     
        iprim(9)=23     
        iprim(10)=29     
        iprim(11)=31     
        iprim(12)=37     
        iprim(13)=41     
        iprim(14)=43     
        iprim(15)=47     
        iprim(16)=53 
!    
	do 1 ip=1,nprim
        do 10 idummy=1,niter
!
	do 2 ik=1,nkp 
	do 3 ir=1,3
	itest=mod(klist(ik,ir),iprim(ip))
	if (itest.ne.0) goto 1
  3	continue
  2	continue
!
        idivi=idivi*iprim(ip)
!       write(*,*) idivi
        do 11 ik=1,nkp 
	do 12 ir=1,3
        klist(ik,ir)=klist(ik,ir)/iprim(ip)
 12     continue
 11     continue
!
 10     continue
  1     continue
        idiv=idiv/idivi
        if(idiv.eq.0) idiv=1
        return
	end


