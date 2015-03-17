c
c configuration common blocks: 
c
 	integer nmax
	parameter(nmax=10)

	real*8 m(nmax),ms(nmax),p(4,nmax),sinv(nmax,nmax)
	real*8 y(nmax),pt(nmax),et(nmax),phi(nmax),dR(nmax,nmax)
	common/par_config/m,ms,p,sinv,y,pt,et,phi,dR
	real*8 sp(4,nmax)
	common/par_config_smeared/sp
	real*8 b(4,2),x1,x2
	common/par_beam/b,x1,x2

	real*8 cgl(2)
	common/angle/cgl
c
c Energy of the photon in CMS:
c
	real*8 ephoton
	common/par_ephoton/ephoton




