	real*8 alphas,mu_f,mu_r
	integer nlf
	common/par_qcd/alphas,mu_f,mu_r,nlf

	real*8 rs,mw,mz,xmz,xmw,gamz,gamw
	common/par_gen/rs,mw,mz,xmw,xmz,gamz,gamw
	integer npart
	common/par_npart/npart

	real*8 pi,pi2,w2,gfermi,alpha0,alpha,sw,cw,mh,xmh,mtop
	real*8 sw2,cw2
        common/weak/sw,cw,sw2,cw2,mh,xmh,mtop,alpha,alpha0,gfermi
	common/para/pi,pi2,w2

	real*8 deltas,loglow,logup
	common/par_deltas/deltas,loglow,logup
	real*8 deltac
	common/par_deltac/deltac

	integer ntot(2),ncut(2),ncut_k(2),ncut_r(2)
	common/par_n/ntot,ncut,ncut_k,ncut_r

	real*8 conti(30),conti_tot(30)
	common/par_conti/conti,conti_tot

	real*8 VKM(6,6)
	data VKM(1,2),VKM(1,4),VKM(1,6)/0.975D0,0.222D0,0D0/
	data VKM(3,2),VKM(3,4),VKM(3,6)/-0.222D0,0.975D0,0D0/
	data VKM(5,2),VKM(5,4),VKM(5,6)/0D0,0D0,0D0/

	data VKM(2,1),VKM(4,1),VKM(6,1)/0.975D0,0.222D0,0D0/
	data VKM(2,3),VKM(4,3),VKM(6,3)/-0.222D0,0.975D0,0D0/
	data VKM(2,5),VKM(4,5),VKM(6,5)/0D0,0D0,0D0/
c
c mass and width of the produced weak gauge boson 
c
	real*8 mv,gamv
	common/par_mv/mv,gamv




















