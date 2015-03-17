	real*8 function kern(xv,wgt)

	implicit none
	integer i,j,k,icut,icutk(3)

	include 'config.inc'
	include 'common.inc'

	integer nndim,ihist,it,ncall2
	common/par_integration/nndim,ihist,it,ncall2

	real*8 ptmin1,ymax1
	common/par_cuts_1/ptmin1,ymax1

	real*8 ptmin2,ymax2
	common/par_cuts_2/ptmin2,ymax2

	real*8 xv(nndim),wgt
	real*8 xk(9),x12(2),x(5)

	real*8 dot_4,fac,sh,t1,t2,temp
	real*8 a2qqz,a2qqz_weak,a2qqz_all,a2qqw,a2qqw_check
	real*8 sig12(2),sig012(2),sigh12(2),sig0u(2),sig0d(2),
	1    sighu(2),sighd(2)
	real*8 sig_virt(2),sigu_virt(2),sigd_virt(2)
	real*8 kern_b,kern_w,kern_s,kern_c,kern_r
	real*8 s,qsq,scale,alphas2
c
c mass matrix
c
	real*8 mass(3)
c
c specification of initial and final state
c
	real*8 qf,ncf,vf,af,qiu,nci,viu,aiu,qid,vid,aid,xm
	real*8 qfs,vfs,afs
	real*8 mff,mffs,mii
	common/par_mass/mff,mffs,mii
c
	real*8 alz,alg,fermvac,dalpf
	real*8 vfeff,viueff,videff
	real*8 sl2eff,su2eff,sd2eff
	common/s2eff/sl2eff,su2eff,sd2eff
c
c only for checks: common block for photon mass	
	real*8 lambda
	common/photon/lambda
c soft form factors
	real*8 k_fac(2),k_faci(2),k_check(2),k_checki(2)
	real*8 fqedu(2),fqedd(2)
c
c o_fac: overall factor (coupling constant)
c sc_fac: scaling factor
c
	real*8 o_fac,sc_fac
c
c pdfx_fac: pdf factor
c
	real*8 pdf_fac(2),pdfc_fac(2),pdfu_fac(2),pdfd_fac(2),
	1    pdfcu_fac(2),pdfcd_fac(2),pdfc2_fac(2)
	real*8 pdfqgc_fac(2),pdfug_fac(2),pdfdg_fac(2)
        real*8 dxpdf1(-6:6),dxpdf2(-6:6)
        real*8 dxpdf1c(-6:6),dxpdf2c(-6:6)
c
c distribution functions:
c
	real*8 fpdf1(-5:5),fpdf2(-5:5)
	real*8 fpdf1c(-5:5),fpdf2c(-5:5)
	real*8 upv1,dnv1,usea1,dsea1,str1,chm1,bot1,glu1,phot1
	real*8 upv2,dnv2,usea2,dsea2,str2,chm2,bot2,glu2,phot2
c CTEQ6
	integer Iset
	real*8 Ctq6Pdf
c
c propagator squared
c
	real*8 den
	common/par_prop/den
c complex Z-propagator:
	complex*16 denzc
	common/par_cprop/denzc
c
c choice of representation:
c
	integer rep,qcd
	common/options/rep,qcd

	real*8 deltar
	common/drit/deltar
c
c switch for width
c
	integer wopt
	common/widthopt/wopt

	integer qqg
	common/par_parton/qqg
c
c change of variables
c
	integer trafo
	common/change/trafo
	real*8 y1,y2
c
c choice of collider (ppbar or pp):
c	
	integer ppswitch
	common/collider/ppswitch
c
c switch for non-resonant contribution
c
	integer qnonr
	common/par_nonres/qnonr
c
c switch for choice of PDFs
c	
	integer ipdf
	common/pdfswitch/ipdf

	real*8 smin,epmin
c
c common block for test
c
	integer test(10)
	common/par_test/test
c
c choice which cutoff-dependence is going to be tested:
c
	integer dsdc
	common/choice/dsdc
c
c choice of factorisation scheme:
c
	real*8 lfc
	common/scheme/lfc
c
c choice, if a collinear cutoff will be imposed (final state)
c
	integer collcut
	common/collinear/collcut
	real*8 z,split,splitq,fcoll
c
c choice of process (CC or NC):
c	
	integer switch
	common/process/switch
c
c common blocks for multiple photon radiation
c
	integer qfsrexp
	real*8 fsrscale
	common/multiple/fsrscale,qfsrexp
	real*8 llq,betf,oxp,xp,dxp,strucf,gllf,zc
	real*8 kern_exp,kern_qcd
	integer nboost
	common/fraction/zc,nboost
c
c common block for QCD switch and flag
c
	integer qqcd,flagqcd
	common/qcdswitch/qqcd
	common/qcdflag/flagqcd
c
c for tests with Madgraph
c
	real*8 ps1(0:3),ps2(0:3),ps3(0:3),ps4(0:3),ps5(0:3)
c
c CC: common block with weak form factor:
c
	complex*16 fvwp(2)
	common/fac_weak/fvwp
c
c initialization
c
	kern_b=0.d0
	kern_w=0.d0
	kern_s=0.d0
	kern_c=0.d0
	kern_r=0.d0
	kern_exp=0d0
	kern_qcd=0d0
	do i=1,30
	   conti(i)=0d0
	end do
c
c hadronic cm energy squared
c
	s=rs*rs
c
c counters:
c
	ntot(1)=ntot(1)+1
	if(Ihist.eq.1)then
	   ntot(2)=ntot(2)+1
	end if
c
c scaling the x's
c
	do i=1,nndim
	  xk(i)=xv(i)
	end do

	sc_fac=1.d0

	x12(1)=xk(1)
	x12(2)=xk(2)

c 0: no trafo, 1: with trafo

	if (trafo.eq.1.and.test(10).eq.0) then

	   y1 = datan(-mv/gamv)
	   y2 = datan((x12(1)*s-mv**2)/mv/gamv)
	   sh = mv**2+gamv*mv*dtan(y1+x12(2)*(y2-y1))
	   x12(2) = sh/s/x12(1)

c in config: x2=x12(2):

	   if (x12(2).gt.1d0) then
	      if ((x12(2)-1d0).le.1d-5) x12(2) = 1d0
	   end if

	   sc_fac=((sh-mv**2)**2+(mv*gamv)**2)*(y2-y1)/x12(1)/s/gamv/mv

	end if

	do i=3,7
	   x(i-2)=xk(i)
	end do

	if (dsdc.eq.1) then
c
c deltas value:
c
	   deltas=10d0**((logup-loglow)*xk(8)+loglow)

	else if (dsdc.eq.2) then
c
c deltac value:
c
	   deltac=10d0**((logup-loglow)*xk(8)+loglow)
	   
	end if

	if(test(7).eq.2) goto 10
	if(qfsrexp.eq.2.and.switch.eq.1) goto 53
c
c calculate born and k_fac:
c
	npart=2
c
c setup for configuration 
c
c convention for generation of final-state configuration:
c
c       1-> lepton2 (neutrino or electron,muon)
c	2-> lepton1 (electron,muon)
c
	mass(1)=mff
	mass(2)=mffs
	if(qqcd.ne.0)then
	   mass(1)=0d0
	   mass(2)=0d0
	endif
	call setup_config(10,npart,mass)
c
c generate the configuration:
c
	call config(npart,ptmin2,s,x12,x,fac)
c
c switch to basic convention:
c
c   p1-> electron,muon
c   p2-> neutrino or electron,muon
c
	do i=1,4
	   temp=p(i,1)
	   p(i,1)=p(i,2)
	   p(i,2)=temp
c	   ptwo(i,1) = p(i,1)
c	   ptwo(i,2) = p(i,2)
	end do
c       
c calculate the invariants :
c
	sinv(1,2)=2d0*dot_4(b(1,1),b(1,2))
	do j=1,npart
	   sinv(1,j+2)=2d0*dot_4(b(1,1),p(1,j))
	   sinv(2,j+2)=2d0*dot_4(b(1,2),p(1,j))
	end do	  
	do i=1,npart
	   do j=i,npart
	      sinv(i+2,j+2)=2d0*dot_4(p(1,i),p(1,j))
	   end do
	end do
	do i=1,npart+2
	   do j=i+1,npart+2
	      sinv(j,i)=sinv(i,j)
	   end do
	end do
c Mandelstam invariants for 2->2 process:
	sh=sinv(1,2)
	t1=-sinv(1,4)
	t2=-sinv(2,4)
c
c momenta for comparison with Madgraph
	ps1(0)=b(4,1)
	ps2(0)=b(4,2)
	ps3(0)=p(4,2)
	ps4(0)=p(4,1)
	do i=1,3
	   ps1(i)=b(i,1)
	   ps2(i)=b(i,2)
	   ps3(i)=p(i,2)
	   ps4(i)=p(i,1)
	end do 
c
c make cuts:
c
	call newcuts(npart,Icut)

	if(dabs(x1).gt.1d0.or.dabs(x2).gt.1d0)icut=1
	if(dabs(x1).lt.1d-5.or.dabs(x2).lt.1d-5)icut=1
	if(Icut.eq.1.or.fac.eq.0d0)then
	   ncut_k(1)=ncut_k(1)+1
	   if(Ihist.eq.1)ncut_k(2)=ncut_k(2)+1
	   kern_b=0d0
	   kern_w=0d0
	   kern_s=0d0
	   kern_c=0d0
	   goto 5
	end if
 
c only parton level result:
	if(test(10).eq.1)then
c CC: only u d initial state
	   if(switch.eq.1)then
	      VKM(2,5)=0d0
	      VKM(5,2)=0d0
	      VKM(4,5)=0d0
	      VKM(5,4)=0d0
	      pdf_fac(1)=1d0
	      pdf_fac(2)=0d0
	   elseif(switch.eq.2)then
c NC: only u ubar initial state
	      pdfu_fac(1)=1d0
	      pdfd_fac(1)=0d0
	      pdfu_fac(2)=0d0
	      pdfd_fac(2)=0d0
	   endif
	   goto 99
	endif

c scale for pdfs:

	if(test(8).lt.3) then
	   qsq=mu_f**2
	elseif(test(8).eq.3)then
c variable scale
	   qsq=2d0*dot_4(sp(1,1),sp(1,2))
	   mu_f=dsqrt(dabs(qsq))
	   mu_r=mu_f
	   alphas=alphas2(mu_r)
	endif
	scale = dsqrt(dabs(qsq))
	if (scale.lt.1.3d0) then
	   kern_b=0d0
	   kern_w=0d0
	   kern_s=0d0
	   kern_c=0d0
	   goto 5
	endif
c
c calculate the distribution functions for p and pbar, the conjugation
c for pbar is done here
c ppswitch=1: ppbar, ppswitch=2: pp
c pdf(i=2,4): u,c-distribution (=u_valence+u_sea)
c pdf(i=-2,-4): ubar,cbar-distribution (=u_sea)
c pdf(i=1,3,5): d,s,b-distribution
c
	if(ipdf.eq.1)then
	   call mrstqed(x1,scale,1,upv1,dnv1,usea1,dsea1,str1,chm1,bot1,
	1	glu1,phot1)
	   call mrstqed(x2,scale,1,upv2,dnv2,usea2,dsea2,str2,chm2,bot2,
	1	glu2,phot2)
	   dxpdf1(0)=glu1
	   dxpdf2(0)=glu2
	   dxpdf1(1)=dnv1+dsea1
	   dxpdf1(-1)=dsea1
	   dxpdf1(3)=str1
	   dxpdf1(-3)=str1
	   dxpdf1(5)=bot1
	   dxpdf1(-5)=bot1
	   dxpdf1(2)=upv1+usea1
	   dxpdf1(-2)=usea1
	   dxpdf1(4)=chm1
	   dxpdf1(-4)=chm1
	   dxpdf2(1)=dnv2+dsea2
	   dxpdf2(-1)=dsea2
	   dxpdf2(3)=str2
	   dxpdf2(-3)=str2
	   dxpdf2(5)=bot2
	   dxpdf2(-5)=bot2
	   dxpdf2(2)=upv2+usea2
	   dxpdf2(-2)=usea2
	   dxpdf2(4)=chm2
	   dxpdf2(-4)=chm2
	   dxpdf1(6)=0d0
	   dxpdf1(-6)=0d0
	   dxpdf2(6)=0d0
	   dxpdf2(-6)=0d0
	   do i=-nlf,nlf
	      fpdf1(i)=dxpdf1(i)/x1 
	      if (ppswitch.eq.1) fpdf2(i)=dxpdf2(-i)/x2
	      if (ppswitch.eq.2) fpdf2(i)=dxpdf2(i)/x2
	   enddo
	else 
c LO: Iset=4, NLO: Iset=400
	   Iset=400
c	   if(qqcd.gt.0.and.test(9).gt.0) Iset=400
	   Call SetCtq6(Iset)
	   do i=-nlf,nlf
	      dxpdf1(i) = Ctq6Pdf(i,x1,scale)
	      dxpdf2(i) = Ctq6Pdf(i,x2,scale)
	   enddo
	   dxpdf1(1) = Ctq6Pdf(2,x1,scale)
	   dxpdf2(1) = Ctq6Pdf(2,x2,scale)
	   dxpdf1(-1) = Ctq6Pdf(-2,x1,scale)
	   dxpdf2(-1) = Ctq6Pdf(-2,x2,scale)
	   dxpdf1(2) = Ctq6Pdf(1,x1,scale)
	   dxpdf2(2) = Ctq6Pdf(1,x2,scale)
	   dxpdf1(-2) = Ctq6Pdf(-1,x1,scale)
	   dxpdf2(-2) = Ctq6Pdf(-1,x2,scale)
	   do i=-nlf,nlf
	      fpdf1(i)=dxpdf1(i)
	      if (ppswitch.eq.1) fpdf2(i)=dxpdf2(-i)
	      if (ppswitch.eq.2) fpdf2(i)=dxpdf2(i)
	   enddo
	endif
c CC:
	if(switch.eq.1)then
c
c test(1) = 1:
c distribution functions factors for u db -> W+ -> l+ nu production
c test(1) = 2:
c distribution functions factors for ub d -> W- -> l- nu production
c
	   pdf_fac(1)=0d0
	   pdf_fac(2)=0d0
	   do i=1,nlf,2
	      do j=2,nlf,2
		 if (test(1).eq.1) then
		    pdf_fac(1)=pdf_fac(1) + fpdf1(j)*fpdf2(-i)*VKM(i,j)**2
		    pdf_fac(2)=pdf_fac(2) + fpdf1(-i)*fpdf2(j)*VKM(i,j)**2
		 elseif (test(1).eq.2) then
		    pdf_fac(2)=pdf_fac(2) + fpdf1(i)*fpdf2(-j)*VKM(i,j)**2
		    pdf_fac(1)=pdf_fac(1) + fpdf1(-j)*fpdf2(i)*VKM(i,j)**2
		 end if
	      end do
	   end do
	elseif(switch.eq.2)then
c
c distribution function factors for 
c u ubar -> Z -> l- l+ production
c
	   pdfu_fac(1)=0d0
	   pdfu_fac(2)=0d0
	   do i=2,nlf,2
	      pdfu_fac(1)=pdfu_fac(1) + fpdf1(i)*fpdf2(-i)
	      pdfu_fac(2)=pdfu_fac(2) + fpdf1(-i)*fpdf2(i)
	   end do
c
c distribution functions factors for 
c d dbar -> Z -> l- l+ production
c
	   pdfd_fac(1)=0d0
	   pdfd_fac(2)=0d0
	   do i=1,nlf,2
	      pdfd_fac(1)=pdfd_fac(1) + fpdf1(i)*fpdf2(-i)
	      pdfd_fac(2)=pdfd_fac(2) + fpdf1(-i)*fpdf2(i)
	   end do
	endif
 99	continue
c
c constant width:
c
	if (wopt.eq.1) then
	   den=(sh-mv**2)**2+(gamv*mv)**2
	   denzc=dcmplx(sinv(1,2)-mz**2,gamz*mz)
c
c s-dependent width:
c
	else if (wopt.eq.2) then
	   den=(sh-mv**2)**2+(gamv*sh/mv)**2
	   denzc=dcmplx(sinv(1,2)-mz**2,gamz*sinv(1,2)/mz)
	end if
c
c couplings for Zff vertices:
c up
	call prop(7,xm,qiu,viu,aiu,nci)
c down
	call prop(10,xm,qid,vid,aid,nci)
c electron/muon 
	if (test(2).eq.1) call prop(4,xm,qf,vf,af,ncf)
	if (test(2).eq.2) call prop(5,xm,qf,vf,af,ncf)

	if(switch.eq.1)then
c
c G_mu-representation:
c
	   if (rep.eq.2) then
	      o_fac=(w2*xmw*gfermi)**2/den*ncf/nci	 
c
c alpha_0 representation:
c
	   else
	      o_fac=(pi*alpha/sw2)**2/den*ncf/nci
	   endif
	elseif(switch.eq.2)then
	   o_fac=pi2*ncf/nci/4d0
	endif
c
c Born-matrixelements squared:
c
	if(switch.eq.1)then
	   call U_Db__e_nu(sig012)
	   kern_b=fac*sc_fac*o_fac*(pdf_fac(1)*sig012(1)+
	1	pdf_fac(2)*sig012(2))
c
	elseif(switch.eq.2)then
c
c EBA:
c	   alz=xmz/pi*w2*gfermi*sw2*cw2
c	   dalpf=-fermvac(sh)
c	   dalpf=-(0.059461d0+20d0/9d0*alpha0/pi*
c	1	dlog(sh/xmz))
c	   alg=alpha0/(1d0-dalpf)
c	   vfeff=(-1d0/2d0-2d0*sl2eff*qf)/2d0/sw/cw
c	   viueff=(1d0/2d0-2d0*su2eff*qiu)/2d0/sw/cw
c	   videff=(-1d0/2d0-2d0*sd2eff*qid)/2d0/sw/cw
c Born:
	   alz=alpha
	   alg=alpha
	   vfeff=vf
	   viueff=viu
	   videff=vid
	   sig0u(1)=a2qqz(sh,t1,t2,vfeff,af,qf,viueff,aiu,qiu,alz,alg)
	   sig0u(2)=a2qqz(sh,t2,t1,vfeff,af,qf,viueff,aiu,qiu,alz,alg)
	   sig0d(1)=a2qqz(sh,t1,t2,vfeff,af,qf,videff,aid,qid,alz,alg)
	   sig0d(2)=a2qqz(sh,t2,t1,vfeff,af,qf,videff,aid,qid,alz,alg)
	   kern_b=fac*sc_fac*o_fac*(pdfu_fac(1)*sig0u(1)+
	1	pdfu_fac(2)*sig0u(2)+pdfd_fac(1)*sig0d(1)+
	2	pdfd_fac(2)*sig0d(2))
	endif
	conti(1)=kern_b
c
c bypass k-factor calculation:
c
	if(test(9).eq.0) goto 5
	if(test(6).eq.0) goto 4
	if(qqcd.eq.2) goto 4
c
c electroweak 1-loop contribution to NC and CC processes:
c
	if (switch.eq.1) then
c
c weak 1-loop contribution:
c
c	conti(2) = kern_b*dreal(fvwp(1)+fvwp(2))
c inclusion of contributions beyond O(alpha):
c	conti(2) = kern_b*
cc	1    (((1d0+fvwp(1)/2d0)*(1d0+fvwp(2)/2d0))**2-1d0)
c	1   ((1d0+dreal(fvwp(1))+0.25d0*(dreal(fvwp(1))**2+
c	1    dimag(fvwp(1))**2))*(1d0+dreal(fvwp(2))+0.25d0*
c	1    (dreal(fvwp(2))**2+dimag(fvwp(2))**2))-1d0)
c
c add non-resonant contributions:
c	if (qnonr.eq.1.and.(dabs(sh-xmw).gt.1d0)) then
c	   sig_virt(1) = a2qqw(sh,t1,t2,sig012(1))
c	   sig_virt(2) = a2qqw(sh,t2,t1,sig012(2))
c	   conti(2) = conti(2)+
c	1	fac*sc_fac*o_fac*(pdf_fac(1)*sig_virt(1)+
c	1       pdf_fac(2)*sig_virt(2))-
c	1       kern_b*(dreal(fvwp(1))+dreal(fvwp(2)))
c	end if
c----begin old code:
	   if (qnonr.eq.0) then
	      conti(2) = kern_b*dreal(fvwp(1)+fvwp(2))
	   else if (qnonr.eq.1) then
c full virtual contribution, where the QED factors 
c have been subtracted(+delta_vs^interf):
c (i.e. for sh=xmw this coincides with kern_b*(fvwp(1)+fvwp(2)) !)
c
	      if (dabs(sh-xmw).le.1d0) then
		 conti(2) = kern_b*dreal(fvwp(1)+fvwp(2))
	      else
		 sig_virt(1) = a2qqw(sh,t1,t2,sig012(1))
		 sig_virt(2) = a2qqw(sh,t2,t1,sig012(2))
		 conti(2) = fac*sc_fac*o_fac*(pdf_fac(1)*sig_virt(1)+
	1	      pdf_fac(2)*sig_virt(2))
	      end if
c --- end of old code
c new code (from nutev calculation)
c full electroweak corrections (no separation in weak and QEd part):
c	      sig_virt(1)=a2qqw2_check(sh,t1,t2,vf,af,vfs,afs,
c	1	   vid,aid,viu,aiu,qid,qiu,qf,mii,mii,mffs,sig012(1))
c	      sig_virt(2)=a2qqw2_check(sh,t2,t1,vf,af,vfs,afs,
c	1	   vid,aid,viu,aiu,qid,qiu,qf,mii,mii,mffs,sig012(2))
c	      conti(2) = fac*sc_fac*o_fac*(pdf_fac(1)*sig_virt(1)+
c	1	   pdf_fac(2)*sig_virt(2))
	   end if
	elseif(switch.eq.2)then
	   sigu_virt(1)=a2qqz_weak(sh,t1,t2,vf,af,qf,
	1	viu,aiu,qiu)
	   sigu_virt(2)=a2qqz_weak(sh,t2,t1,vf,af,qf,
	1	viu,aiu,qiu)
	   sigd_virt(1)=a2qqz_weak(sh,t1,t2,vf,af,qf,
	1	vid,aid,qid)
	   sigd_virt(2)=a2qqz_weak(sh,t2,t1,vf,af,qf,
	1	vid,aid,qid)
c new: weak+QED (from Nutev calcul.)
c	   sigu_virt(1)=a2qqz_all(sh,t1,t2,vf,af,qf,
c	1	viu,aiu,qiu,sig0u(1))
c	   sigu_virt(2)=a2qqz_all(sh,t2,t1,vf,af,qf,
c	1	viu,aiu,qiu,sig0u(2))
c	   sigd_virt(1)=a2qqz_all(sh,t1,t2,vf,af,qf,
c	1	vid,aid,qid,sig0d(1))
c	   sigd_virt(2)=a2qqz_all(sh,t2,t1,vf,af,qf,
c	1	vid,aid,qid,sig0d(2))

	   conti(2) = fac*sc_fac*o_fac*(pdfu_fac(1)*sigu_virt(1)+
	1	pdfu_fac(2)*sigu_virt(2)+
	2	pdfd_fac(1)*sigd_virt(1)+pdfd_fac(2)*sigd_virt(2))
	endif
	kern_w = conti(2)

 4	continue
c
c bypass Born calculation:
c
	if (test(9).eq.1) then
	   conti(1) = 0d0
	   kern_b = 0d0
	end if

	if(test(6).eq.2) goto 5
	if(qqg.eq.2) goto 52

c soft photon radiation (+final-state collinear radiation+PDF CT)
c virt.+soft+PDf CT gluon radiation

	if (switch.eq.1) then
c CC
	   if (test(3).le.3) then

	      call softcollqed_cc(test(3),k_fac)
	   
	   else if(test(3).eq.4) then

	      k_fac(1) = 0d0
	      k_fac(2) = 0d0
	      k_check(1) = 0d0
	      k_check(2) = 0d0
	      do k=1,3
		 call softcollqed_cc(k,k_faci)
		 k_fac(1) = k_fac(1)+k_faci(1)
		 k_fac(2) = k_fac(2)+k_faci(2)

c 	          lambda=1d0
c	          call fvsoftcoll(k,k_checki)
c	          k_check(1) = k_check(1)+k_checki(1)
c	          k_check(2) = k_check(2)+k_checki(2)

	      end do
	   end if
c --begin checks--
c calculate YFS form factors (delta_interf._V+S has been subtracted):
c	call fvyfs(sh,t1,form)
c check if QED factor==YFS+soft(+coll.) BR contribution:	
c	write(6,*)'1',k_fac(1),k_check(1)+2d0*form
c calculate the complete virtual photon contribution
c on-shell singularities have been subtracted !
c	sig_check=a2qqw_check(sh,t1,sig12(1))
c check if complete == res+non-res. part:	
c	write(6,*)'2',(sig_check+k_check(1)*sig12(1))/sig12(1),
c	1    k_fac(1)+sig_virt(1)/sig12(1)
c --end checks--
*
	   conti(3)=fac*sc_fac*o_fac*(pdf_fac(1)*k_fac(1)*sig012(1)
	1	+pdf_fac(2)*k_fac(2)*sig012(2))
c virt+soft+soft PDF CT:
	   if(qqcd.eq.1)then
	      call virtsoftqcd_cc(k_fac)
	      conti(3)=conti(3)+
	1	   fac*sc_fac*o_fac*(pdf_fac(1)*k_fac(1)*sig012(1)
	1	   +pdf_fac(2)*k_fac(2)*sig012(2))
	   elseif(qqcd.eq.2)then
	      call virtsoftqcd_cc(k_fac)
	      conti(3)=fac*sc_fac*o_fac*(pdf_fac(1)*k_fac(1)*sig012(1)
	1	   +pdf_fac(2)*k_fac(2)*sig012(2))
	   endif

	elseif (switch.eq.2) then
c NC
	   call softcollqed_nc(sh,t1,qiu,mii,fqedu(1))
	   call softcollqed_nc(sh,t1,qid,mii,fqedd(1))
	   call softcollqed_nc(sh,t2,qiu,mii,fqedu(2))
	   call softcollqed_nc(sh,t2,qid,mii,fqedd(2))

	   conti(3)=fac*sc_fac*o_fac*(pdfu_fac(1)*sig0u(1)*fqedu(1)
	2	+pdfd_fac(1)*sig0d(1)*fqedd(1))
	end if
	kern_s=conti(3)

 66	continue
c
c collinear initial-state photon and gluon radiation:
c
c mapping of 1/(1-z):
	z=1d0-dexp(dlog(deltas)*xk(9)+(1d0-xk(9))*dlog(1d0-x2))
c	z=xk(9)*(1d0-deltas-x2)+x2
	if (x2/z.gt.1d0.or.x2/z.le.1d-5.or.dabs(z).gt.1d0) then
	   conti(4)=0d0
	   goto 51
	end if
c only parton level result:
	if(test(10).eq.1)then
c CC: only u d initial state
	   if(switch.eq.1)then
	      VKM(2,5)=0d0
	      VKM(5,2)=0d0
	      VKM(4,5)=0d0
	      VKM(5,4)=0d0
	      pdfc_fac(1)=1d0
	      pdfc_fac(2)=0d0
	   elseif(switch.eq.2)then
c NC: only u ubar initial state
	      pdfcu_fac(1)=1d0
	      pdfcd_fac(1)=0d0
	      pdfcu_fac(2)=0d0
	      pdfcd_fac(2)=0d0
	   endif
	   goto 77
	endif
	if(ipdf.eq.1)then
	   call mrstqed(x2/z,scale,1,upv2,dnv2,usea2,dsea2,str2,chm2,bot2,
	1	glu2,phot2)
	   dxpdf2c(0)=glu2
	   dxpdf2c(1)=dnv2+dsea2
	   dxpdf2c(-1)=dsea2
	   dxpdf2c(3)=str2
	   dxpdf2c(-3)=str2
	   dxpdf2c(5)=bot2
	   dxpdf2c(-5)=bot2
	   dxpdf2c(2)=upv2+usea2
	   dxpdf2c(-2)=usea2
	   dxpdf2c(4)=chm2
	   dxpdf2c(-4)=chm2
	   dxpdf2c(6)=0d0
	   dxpdf2c(-6)=0d0
	   do i=-nlf,nlf
	      if (ppswitch.eq.1) fpdf2c(i)=dxpdf2c(-i)/x2*z
	      if (ppswitch.eq.2) fpdf2c(i)=dxpdf2c(i)/x2*z
	   enddo
	else
c LO: Iset=4, NLO: Iset=400
	   Iset=400
c	   if(qqcd.gt.0) Iset=400
	   Call SetCtq6(Iset)
	   do i=-nlf,nlf
	      dxpdf2c(i) = Ctq6Pdf(i,x2/z,scale)
	   enddo
	   dxpdf2c(1) = Ctq6Pdf(2,x2/z,scale)
	   dxpdf2c(-1) = Ctq6Pdf(-2,x2/z,scale)
	   dxpdf2c(2) = Ctq6Pdf(1,x2/z,scale)
	   dxpdf2c(-2) = Ctq6Pdf(-1,x2/z,scale)
	   do i=-nlf,nlf
	      if (ppswitch.eq.1) fpdf2c(i)=dxpdf2c(-i)
	      if (ppswitch.eq.2) fpdf2c(i)=dxpdf2c(i)
	   enddo
	endif

	if(switch.eq.1)then
c
c distribution functions factors for
c u db -> W+-> l+ nu production or
c ub d -> W--> l- nu production   
c   
	   pdfc_fac(1)=0d0
	   pdfc_fac(2)=0d0
	   pdfc2_fac(1)=0d0
	   pdfc2_fac(2)=0d0
	   do i=1,nlf,2
	      do j=2,nlf,2
		 if (test(1).eq.1) then
		    pdfc_fac(1)=pdfc_fac(1)+fpdf1(j)*fpdf2c(-i)*
	1		 VKM(i,j)**2*1d0/9d0
		    pdfc_fac(2)=pdfc_fac(2)+fpdf1(-i)*fpdf2c(j)*
	1		 VKM(i,j)**2*4d0/9d0
c PDFs for QCD corrections
		    pdfc2_fac(1)=pdfc2_fac(1)+fpdf1(j)*fpdf2c(-i)*
	1		 VKM(i,j)**2
		    pdfc2_fac(2)=pdfc2_fac(2)+fpdf1(-i)*fpdf2c(j)*
	1		 VKM(i,j)**2
		 elseif (test(1).eq.2) then
		    pdfc_fac(2)=pdfc_fac(2)+fpdf1(i)*fpdf2c(-j)*
	1		 VKM(i,j)**2*4d0/9d0
		    pdfc_fac(1)=pdfc_fac(1)+fpdf1(-j)*fpdf2c(i)*
	1		 VKM(i,j)**2*1d0/9d0
c PDFs for QCD corrections
		    pdfc2_fac(2)=pdfc2_fac(2)+fpdf1(i)*fpdf2c(-j)*
	1		 VKM(i,j)**2
		    pdfc2_fac(1)=pdfc2_fac(1)+fpdf1(-j)*fpdf2c(i)*
	1		 VKM(i,j)**2
		 end if
	      end do
	   end do
	elseif(switch.eq.2)then
c
c distribution function factors for 
c u ubar -> Z -> l- l+ production
c
	   pdfcu_fac(1)=0d0
	   pdfcu_fac(2)=0d0
	   do i=2,nlf,2
	      pdfcu_fac(1)=pdfcu_fac(1) + fpdf1(i)*fpdf2c(-i)*4d0/90
	      pdfcu_fac(2)=pdfcu_fac(2) + fpdf1(-i)*fpdf2c(i)*4d0/90
	   end do
c
c distribution functions factors for 
c d dbar -> Z -> l- l+ production
c
	   pdfcd_fac(1)=0d0
	   pdfcd_fac(2)=0d0
	   do i=1,nlf,2
	      pdfcd_fac(1)=pdfcd_fac(1) + fpdf1(i)*fpdf2c(-i)*1d0/90
	      pdfcd_fac(2)=pdfcd_fac(2) + fpdf1(-i)*fpdf2c(i)*1d0/90
	   end do
	endif
 77	continue
c
c collinear initial-state photon radiation+coll. part of PDF CT
c DIS or minimal factorisation:
c
	fcoll = (1d0+z**2)/(1d0-z)*dlog((1d0-z)/z)-
	1    3d0/2d0/(1d0-z)+2d0*z+3d0
	
	split = alpha0/pi/2d0*((1d0+z**2)/(1d0-z)*
	1    dlog(sinv(1,2)/mu_f**2*(1d0-z)**2/z*deltac/2d0)+
	2    1d0-z-lfc*fcoll)

	if(switch.eq.1) then
c CC: 
	   conti(4)=fac*sc_fac*o_fac*(pdfc_fac(1)*sig012(1)+
	3	pdfc_fac(2)*sig012(2))*split*
c	2	(1d0-deltas-x2)/z
	2	(dlog(1d0-x2)-dlog(deltas))*(1d0-z)/z

	   if(qqcd.eq.1)then
c QCD PDF CT in MSbar scheme (lfc=0)
	      splitq = alphas/pi/2d0*4d0/3d0*((1d0+z**2)/(1d0-z)*
	1	   dlog(sinv(1,2)/mu_f**2*(1d0-z)**2/z*deltac/2d0)+
	2	   1d0-z)
	      conti(4)=conti(4)+
	1	   fac*sc_fac*o_fac*(pdfc2_fac(1)*sig012(1)+
	3	   pdfc2_fac(2)*sig012(2))*splitq*
	2	   (1d0-deltas-x2)/z
c	2	   (dlog(1d0-x2)-dlog(deltas))*(1d0-z)/z
	   elseif(qqcd.eq.2)then
c QCD PDF CT in MSbar scheme (lfc=0)
	      splitq = alphas/pi/2d0*4d0/3d0*((1d0+z**2)/(1d0-z)*
	1	   dlog(sinv(1,2)/mu_f**2*(1d0-z)**2/z*deltac/2d0)+
	2	   1d0-z)
	      conti(4)=fac*sc_fac*o_fac*(pdfc2_fac(1)*sig012(1)+
	3	   pdfc2_fac(2)*sig012(2))*splitq*
	2	   (1d0-deltas-x2)/z
c	2	   (dlog(1d0-x2)-dlog(deltas))*(1d0-z)/z
	   endif

	elseif(switch.eq.2) then
c NC:  
	   conti(4)=fac*sc_fac*o_fac*(pdfcu_fac(1)*sig0u(1)+
	1	pdfcu_fac(2)*sig0u(2)+pdfcd_fac(1)*sig0d(1)+
	2	pdfcd_fac(2)*sig0d(2))*alpha0/pi*split*
	2	(1d0-deltas-x2)/z
c	2	(dlog(1d0-x2)-dlog(deltas))*(1d0-z)/z
	endif

 51	continue

c mapping of 1/(1-z):
	z=1d0-dexp(dlog(deltas)*xk(9)+(1d0-xk(9))*dlog(1d0-x1))
c	z=xk(9)*(1d0-deltas-x1)+x1
	if (x1/z.gt.1d0.or.x1/z.le.1d-5.or.dabs(z).gt.1d0) then
	   conti(5)=0d0
	   goto 52
	end if
c only parton level result:
	if(test(10).eq.1)then
c CC: only u d initial state
	   if(switch.eq.1)then
	      VKM(2,5)=0d0
	      VKM(5,2)=0d0
	      VKM(4,5)=0d0
	      VKM(5,4)=0d0
	      pdfc_fac(1)=1d0
	      pdfc_fac(2)=0d0
	   elseif(switch.eq.2)then
c NC: only u ubar initial state
	      pdfcu_fac(1)=1d0
	      pdfcd_fac(1)=0d0
	      pdfcu_fac(2)=0d0
	      pdfcd_fac(2)=0d0
	   endif
	   goto 88
	endif
	if(ipdf.eq.1)then
	   call mrstqed(x1/z,scale,1,upv1,dnv1,usea1,dsea1,str1,chm1,bot1,
	1	glu1,phot1)
	   dxpdf1c(0)=glu1
	   dxpdf1c(1)=dnv1+dsea1
	   dxpdf1c(-1)=dsea1
	   dxpdf1c(3)=str1
	   dxpdf1c(-3)=str1
	   dxpdf1c(5)=bot1
	   dxpdf1c(-5)=bot1
	   dxpdf1c(2)=upv1+usea1
	   dxpdf1c(-2)=usea1
	   dxpdf1c(4)=chm1
	   dxpdf1c(-4)=chm1
	   dxpdf1c(6)=0d0
	   dxpdf1c(-6)=0d0
	   do i=-nlf,nlf
	      fpdf1c(i)=dxpdf1c(i)/x1*z 
	   enddo
	else
c LO: Iset=4, NLO: Iset=400
	   Iset=400
c	   if(qqcd.gt.0) Iset=400
	   Call SetCtq6(Iset)
	   do i=-nlf,nlf
	      dxpdf1c(i) = Ctq6Pdf(i,x1/z,scale)
	   enddo
	   dxpdf1c(1) = Ctq6Pdf(2,x1/z,scale)
	   dxpdf1c(-1) = Ctq6Pdf(-2,x1/z,scale)
	   dxpdf1c(2) = Ctq6Pdf(1,x1/z,scale)
	   dxpdf1c(-2) = Ctq6Pdf(-1,x1/z,scale)
	   do i=-nlf,nlf
	      fpdf1c(i)=dxpdf1c(i)
	   enddo
	endif

	if(switch.eq.1)then
c
c distribution functions factors for
c u db -> W+-> l+ nu production or
c ub d -> W--> l- nu production   
c   
	   pdfc_fac(1)=0d0
	   pdfc_fac(2)=0d0
	   pdfc2_fac(1)=0d0
	   pdfc2_fac(2)=0d0
	   do i=1,nlf,2
	      do j=2,nlf,2
		 if (test(1).eq.1) then
		    pdfc_fac(1)=pdfc_fac(1)+fpdf1c(j)*fpdf2(-i)*
	1		 VKM(i,j)**2*4d0/9d0
		    pdfc_fac(2)=pdfc_fac(2)+fpdf1c(-i)*fpdf2(j)*
	1		 VKM(i,j)**2*1d0/9d0
c PDFs for QCD corrections
		    pdfc2_fac(1)=pdfc2_fac(1)+fpdf1c(j)*fpdf2(-i)*
	1		 VKM(i,j)**2
		    pdfc2_fac(2)=pdfc2_fac(2)+fpdf1c(-i)*fpdf2(j)*
	1		 VKM(i,j)**2
		 elseif (test(1).eq.2) then
		    pdfc_fac(2)=pdfc_fac(2)+fpdf1c(i)*fpdf2(-j)*
	1		 VKM(i,j)**2*1d0/9d0
		    pdfc_fac(1)=pdfc_fac(1)+fpdf1c(-j)*fpdf2(i)*
	1		 VKM(i,j)**2*4d0/9d0
c PDFs for QCD corrections
		    pdfc2_fac(2)=pdfc2_fac(2)+fpdf1c(i)*fpdf2(-j)*
	1		 VKM(i,j)**2
		    pdfc2_fac(1)=pdfc2_fac(1)+fpdf1c(-j)*fpdf2(i)*
	1		 VKM(i,j)**2
		 end if
	      end do
	   end do
	elseif(switch.eq.2)then
c
c distribution function factors for 
c u ubar -> Z -> l- l+ production
c
	   pdfcu_fac(1)=0d0
	   pdfcu_fac(2)=0d0
	   do i=2,nlf,2
	      pdfcu_fac(1)=pdfcu_fac(1) + fpdf1c(i)*fpdf2(-i)*4d0/90
	      pdfcu_fac(2)=pdfcu_fac(2) + fpdf1c(-i)*fpdf2(i)*4d0/90
	   end do
c
c distribution functions factors for 
c d dbar -> Z -> l- l+ production
c
	   pdfcd_fac(1)=0d0
	   pdfcd_fac(2)=0d0
	   do i=1,nlf,2
	      pdfcd_fac(1)=pdfcd_fac(1) + fpdf1c(i)*fpdf2(-i)*1d0/90
	      pdfcd_fac(2)=pdfcd_fac(2) + fpdf1c(-i)*fpdf2(i)*1d0/90
	   end do
	endif
 88	continue
c
c collinear initial-state photon radiation+coll. part of PDF CT
c DIS or minimal factorisation:
c
	fcoll = (1d0+z**2)/(1d0-z)*dlog((1d0-z)/z)-
	1    3d0/2d0/(1d0-z)+2d0*z+3d0
	
	split = alpha0/pi/2d0*((1d0+z**2)/(1d0-z)*
	1    dlog(sinv(1,2)/mu_f**2*(1d0-z)**2/z*deltac/2d0)+
	2    1d0-z-lfc*fcoll)

	if(switch.eq.1) then
c CC: 
	   conti(5)=fac*sc_fac*o_fac*(pdfc_fac(1)*sig012(1)+
	3	pdfc_fac(2)*sig012(2))*split*
c	2	(1d0-deltas-x1)/z
	2	(dlog(1d0-x1)-dlog(deltas))*(1d0-z)/z

	   if(qqcd.eq.1)then
c QCD PDF CT in MSbar scheme (lfc=0)
	      splitq = alphas/pi/2d0*4d0/3d0*((1d0+z**2)/(1d0-z)*
	1	   dlog(sinv(1,2)/mu_f**2*(1d0-z)**2/z*deltac/2d0)+
	2	   1d0-z)
	      conti(5)=conti(5)+
	1	   fac*sc_fac*o_fac*(pdfc2_fac(1)*sig012(1)+
	3	   pdfc2_fac(2)*sig012(2))*splitq*
	2	   (1d0-deltas-x1)/z
c	2	   (dlog(1d0-x1)-dlog(deltas))*(1d0-z)/z
	   elseif(qqcd.eq.2)then
c QCD PDF CT in MSbar scheme (lfc=0)
	      splitq = alphas/pi/2d0*4d0/3d0*((1d0+z**2)/(1d0-z)*
	1	   dlog(sinv(1,2)/mu_f**2*(1d0-z)**2/z*deltac/2d0)+
	2	   1d0-z)
	      conti(5)=fac*sc_fac*o_fac*(pdfc2_fac(1)*sig012(1)+
	3	   pdfc2_fac(2)*sig012(2))*splitq*
	2	   (1d0-deltas-x1)/z
c	2	   (dlog(1d0-x1)-dlog(deltas))*(1d0-z)/z
	   endif

	elseif(switch.eq.2) then
c NC:  
	   conti(5)=fac*sc_fac*o_fac*(pdfcu_fac(1)*sig0u(1)+
	1	pdfcu_fac(2)*sig0u(2)+pdfcd_fac(1)*sig0d(1)+
	2	pdfcd_fac(2)*sig0d(2))*alpha0/pi*split*
	2	(1d0-deltas-x1)/z
c	2	(dlog(1d0-x1)-dlog(deltas))*(1d0-z)/z
	endif

 52	continue
c qg splitting:
	if (qqg.eq.1.or.switch.eq.2) goto 24
c
c initial state collinear contribution (including collinear part
c of PDF counter term !)
	z=xk(9)*(1d0-x2)+x2
	if (x2/z.gt.1d0.or.x2.gt.1d0) then
	   icut=1
	   if (qqg.eq.2) conti(4)=0d0
	   goto 23
	end if
	if(ipdf.eq.1)then
	   call mrstqed(x2/z,scale,1,upv2,dnv2,usea2,dsea2,str2,chm2,bot2,
	1	glu2,phot2)
	   dxpdf2c(0)=glu2
	   dxpdf2c(1)=dnv2+dsea2
	   dxpdf2c(-1)=dsea2
	   dxpdf2c(3)=str2
	   dxpdf2c(-3)=str2
	   dxpdf2c(5)=bot2
	   dxpdf2c(-5)=bot2
	   dxpdf2c(2)=upv2+usea2
	   dxpdf2c(-2)=usea2
	   dxpdf2c(4)=chm2
	   dxpdf2c(-4)=chm2
	   dxpdf2c(6)=0d0
	   dxpdf2c(-6)=0d0
	   do i=-nlf,nlf
	      if (ppswitch.eq.1) fpdf2c(i)=dxpdf2c(-i)/x2*z
	      if (ppswitch.eq.2) fpdf2c(i)=dxpdf2c(i)/x2*z
	   enddo
	else
c LO: Iset=4, NLO: Iset=400
	   Iset=400
c	   if(qqcd.gt.0) Iset=400
	   Call SetCtq6(Iset)
	   do i=-nlf,nlf
	      dxpdf2c(i) = Ctq6Pdf(i,x2/z,scale)
	   enddo
	   dxpdf2c(1) = Ctq6Pdf(2,x2/z,scale)
	   dxpdf2c(-1) = Ctq6Pdf(-2,x2/z,scale)
	   dxpdf2c(2) = Ctq6Pdf(1,x2/z,scale)
	   dxpdf2c(-2) = Ctq6Pdf(-1,x2/z,scale)
	   do i=-nlf,nlf
	      if (ppswitch.eq.1) fpdf2c(i)=dxpdf2c(-i)
	      if (ppswitch.eq.2) fpdf2c(i)=dxpdf2c(i)
	   enddo
	endif
*
	pdfqgc_fac(1)=0.d0
	pdfqgc_fac(2)=0.d0
c W+: dbar+sbar and u,c initial states
c W-: ubar+cbar and d,s initial states
c CKM squared add up to one
	if(test(1).eq.1)then
	   do i=1,nlf,2
	      pdfqgc_fac(2)=pdfqgc_fac(2) + fpdf1(-i)*fpdf2c(0)
	   end do
	   do i=2,nlf,2
	      pdfqgc_fac(1)=pdfqgc_fac(1) + fpdf1(i)*fpdf2c(0)
	   end do
	elseif(test(1).eq.2)then
	   do i=1,nlf,2
	      pdfqgc_fac(1)=pdfqgc_fac(1) + fpdf1(i)*fpdf2c(0)
	   end do
	   do i=2,nlf,2
	      pdfqgc_fac(2)=pdfqgc_fac(2) + fpdf1(-i)*fpdf2c(0)
	   end do
	endif
	
c g->qq splitting:
	split = alphas/2d0/pi*((z**2+(1d0-z)**2)/2d0*
	1    dlog(sinv(1,2)/mu_f**2*(1d0-z)**2/z*deltac/2d0)+z*(1d0-z))
         
	conti(4)=conti(4)+ 
	1    fac*sc_fac*o_fac*(pdfqgc_fac(1)*sig012(1)+
	2    pdfqgc_fac(2)*sig012(2))*split*
	3    (1d0-x2)/z

 23	continue

	z=xk(9)*(1d0-x1)+x1
	if (x1/z.gt.1d0.or.x1.gt.1d0) then
	   icut=1
	   if(qqg.eq.2)conti(5)=0d0
	   goto 24
	end if
	if(ipdf.eq.1)then
	   call mrstqed(x1/z,scale,1,upv1,dnv1,usea1,dsea1,str1,chm1,bot1,
	1	glu1,phot1)
	   dxpdf1c(0)=glu1
	   dxpdf1c(1)=dnv1+dsea1
	   dxpdf1c(-1)=dsea1
	   dxpdf1c(3)=str1
	   dxpdf1c(-3)=str1
	   dxpdf1c(5)=bot1
	   dxpdf1c(-5)=bot1
	   dxpdf1c(1)=upv1+usea1
	   dxpdf1c(-2)=usea1
	   dxpdf1c(4)=chm1
	   dxpdf1c(-4)=chm1
	   dxpdf1c(6)=0d0
	   dxpdf1c(-6)=0d0
	   do i=-nlf,nlf
	      fpdf1c(i)= dxpdf1c(i)/x1*z 
	   end do
	else
c LO: Iset=4, NLO: Iset=400
	   Iset=400
c	   if(qqcd.gt.0) Iset=400
	   Call SetCtq6(Iset)
	   do i=-nlf,nlf
	      dxpdf2c(i) = Ctq6Pdf(i,x1/z,scale)
	   enddo
	   dxpdf1c(1) = Ctq6Pdf(2,x1/z,scale)
	   dxpdf1c(-1) = Ctq6Pdf(-2,x1/z,scale)
	   dxpdf1c(2) = Ctq6Pdf(1,x1/z,scale)
	   dxpdf1c(-2) = Ctq6Pdf(-1,x1/z,scale)
	   do i=-nlf,nlf
	      fpdf1c(i)= dxpdf1c(i)
	   end do
	endif
*       
	pdfqgc_fac(1)=0.d0
	pdfqgc_fac(2)=0.d0
c W+: dbar+sbar and u,c initial states
c W-: ubar+cbar and d,s initial states
c CKM squared add up to one
	if(test(1).eq.1)then
	   do i=1,nlf,2
	      pdfqgc_fac(1)=pdfqgc_fac(1) + fpdf2(-i)*fpdf1c(0)
	   end do
	   do i=2,nlf,2
	      pdfqgc_fac(2)=pdfqgc_fac(2) + fpdf2(i)*fpdf1c(0)
	   end do
	elseif(test(1).eq.2)then
	   do i=1,nlf,2
	      pdfqgc_fac(2)=pdfqgc_fac(2) + fpdf2(i)*fpdf1c(0)
	   end do
	   do i=2,nlf,2
	      pdfqgc_fac(1)=pdfqgc_fac(1) + fpdf2(-i)*fpdf1c(0)
	   end do
	endif
c       
c minimal factorisation (MSbar scheme):
c g->qq splitting:
	split = alphas/2d0/pi*((z**2+(1d0-z)**2)/2d0*
	1    dlog(sinv(1,2)/mu_f**2*(1d0-z)**2/z*deltac/2d0)+z*(1d0-z))
	
	conti(4)=conti(4)+ 
	1    fac*sc_fac*o_fac*(pdfqgc_fac(1)*sig012(1)+
	2    pdfqgc_fac(2)*sig012(2))*split*
	1    (1d0-x1)/z

 24	continue      
	kern_c=conti(4)+conti(5)

 53	continue

	if (qfsrexp.eq.0.or.switch.eq.2.or.qqcd.eq.2) goto 5 
c
c Born convoluted with structure function
c
	npart=2
c
c setup for configuration 
c
c convention for generation of final-state configuration:
c
c       1-> lepton2 (neutrino or electron,muon)
c	2-> lepton1 (electron,muon)
c
	mass(1)=mff
	mass(2)=mffs
	call setup_config(10,npart,mass)
c
c generate the configuration:
c
	call config(npart,ptmin2,s,x12,x,fac)
c
c switch to basic convention:
c
c   1-> electron,muon
c   2-> neutrino
c
	do i=1,4
	   temp=p(i,1)
	   p(i,1)=p(i,2)
	   p(i,2)=temp
c	   ptwo(i,1) = p(i,1)
c	   ptwo(i,2) = p(i,2)
	end do
c
c calculate the invariants :
c
	sinv(1,2)=2d0*dot_4(b(1,1),b(1,2))
	do j=1,npart
	   sinv(1,j+2)=2d0*dot_4(b(1,1),p(1,j))
	   sinv(2,j+2)=2d0*dot_4(b(1,2),p(1,j))
	end do	  
	do i=1,npart
	   do j=i,npart
	      sinv(i+2,j+2)=2d0*dot_4(p(1,i),p(1,j))
	   end do
	end do
	do i=1,npart+2
	   do j=i+1,npart+2
	      sinv(j,i)=sinv(i,j)
	   end do
	end do
	sh=sinv(1,2)
	t1=-sinv(1,4)
	t2=-sinv(2,4)
c
	fsrscale=dsqrt(dabs(sh))
c
c make cuts:
c
	llq=dlog(fsrscale**2/mffs**2)
	betf=alpha0/pi*(llq-1d0)
	oxp=xk(9)**(1d0/betf)
	if (oxp.lt.1d-10) oxp=1d-10
	zc=1d0-oxp
	dxp=1d0/betf
	nboost=1
	if(zc.lt.dabs(mffs/p(4,1))) then
	   icut=1
	else
	   call newcuts(npart,Icut)
	endif

	if(Icut.eq.1)then
	   ncut_k(1)=ncut_k(1)+1
	   if(Ihist.eq.1)ncut_k(2)=ncut_k(2)+1
	   kern_exp=0d0
	   conti(6)=0d0
	   goto 5
	end if
	if(x1.lt.1d-5.or.x1.gt.1d0.or.
	1    x2.lt.1d-5.or.x2.gt.1d0)then
	   kern_exp=0d0
	   conti(6)=0d0
	   goto 5
	endif
c
c calculate the distribution functions for p and pbar, the conjugation
c for pbar is done here
c ppswitch=1: ppbar, ppswitch=2: pp
c pdf(i=2,4): u,c-distribution (=u_valence+u_sea)
c pdf(i=-2,-4): ubar,cbar-distribution (=u_sea)
c pdf(i=1,3,5): d,s,b-distribution
c
	if(test(8).lt.3) then
	   qsq=mu_f**2
	elseif(test(8).eq.3)then
c variable scale
	   qsq=2d0*dot_4(sp(1,1),sp(1,2))
	   mu_f=dsqrt(dabs(qsq))
	   mu_r=mu_f
	   alphas=alphas2(mu_r)
	endif
	scale = dsqrt(dabs(qsq))
c	if (scale.lt.1.3d0) then
c	   kern_exp=0d0
c	   conti(6)=0d0
c	   goto 5
c	endif
	   
	if(ipdf.eq.1)then
	   call mrstqed(x1,scale,1,upv1,dnv1,usea1,dsea1,str1,chm1,bot1,
	1	glu1,phot1)
	   call mrstqed(x2,scale,1,upv2,dnv2,usea2,dsea2,str2,chm2,bot2,
	1	glu2,phot2)
	   dxpdf1(0)=glu1
	   dxpdf2(0)=glu2
	   dxpdf1(1)=dnv1+dsea1
	   dxpdf1(-1)=dsea1
	   dxpdf1(3)=str1
	   dxpdf1(-3)=str1
	   dxpdf1(5)=bot1
	   dxpdf1(-5)=bot1
	   dxpdf1(2)=upv1+usea1
	   dxpdf1(-2)=usea1
	   dxpdf1(4)=chm1
	   dxpdf1(-4)=chm1
	   dxpdf2(1)=dnv2+dsea2
	   dxpdf2(-1)=dsea2
	   dxpdf2(3)=str2
	   dxpdf2(-3)=str2
	   dxpdf2(5)=bot2
	   dxpdf2(-5)=bot2
	   dxpdf2(2)=upv2+usea2
	   dxpdf2(-2)=usea2
	   dxpdf2(4)=chm2
	   dxpdf2(-4)=chm2
	   dxpdf1(6)=0d0
	   dxpdf1(-6)=0d0
	   dxpdf2(6)=0d0
	   dxpdf2(-6)=0d0
	   do i=-nlf,nlf
	      fpdf1(i)=dxpdf1(i)/x1 
	      if (ppswitch.eq.1) fpdf2(i)=dxpdf2(-i)/x2
	      if (ppswitch.eq.2) fpdf2(i)=dxpdf2(i)/x2
	   enddo
	else
c LO: Iset=4, NLO: Iset=400
	   Iset=400
c	   if(qqcd.gt.0) Iset=400
	   Call SetCtq6(Iset)
	   do i=-nlf,nlf
	      dxpdf1(i) = Ctq6Pdf(i,x1,scale)
	      dxpdf2(i) = Ctq6Pdf(i,x2,scale)
	   enddo
	   dxpdf1(1) = Ctq6Pdf(2,x1,scale)
	   dxpdf2(1) = Ctq6Pdf(2,x2,scale)
	   dxpdf1(-1) = Ctq6Pdf(-2,x1,scale)
	   dxpdf2(-1) = Ctq6Pdf(-2,x2,scale)
	   dxpdf1(2) = Ctq6Pdf(1,x1,scale)
	   dxpdf2(2) = Ctq6Pdf(1,x2,scale)
	   dxpdf1(-2) = Ctq6Pdf(-1,x1,scale)
	   dxpdf2(-2) = Ctq6Pdf(-1,x2,scale)
	   do i=-nlf,nlf
	      fpdf1(i)=dxpdf1(i)
	      if (ppswitch.eq.1) fpdf2(i)=dxpdf2(-i)
	      if (ppswitch.eq.2) fpdf2(i)=dxpdf2(i)
	   enddo
	endif
c
c test(1) = 1:
c distribution functions factors for u db -> W+ -> l+ nu production
c test(1) = 2:
c distribution functions factors for ub d -> W- -> l- nu production
c
        pdf_fac(1)=0.d0
        pdf_fac(2)=0.d0
	do i=1,nlf,2
	   do j=2,nlf,2
	      if (test(1).eq.1) then
		 pdf_fac(1)=pdf_fac(1) + fpdf1(j)*fpdf2(-i)*VKM(i,j)**2
		 pdf_fac(2)=pdf_fac(2) + fpdf1(-i)*fpdf2(j)*VKM(i,j)**2
	      elseif (test(1).eq.2) then
		 pdf_fac(2)=pdf_fac(2) + fpdf1(i)*fpdf2(-j)*VKM(i,j)**2
		 pdf_fac(1)=pdf_fac(1) + fpdf1(-j)*fpdf2(i)*VKM(i,j)**2
	      end if
	   end do
	end do
c
c constant width:
c
	if (wopt.eq.1) then
	   den=(sh-mv**2)**2+(gamv*mv)**2
c
c s-dependent width:
c
	else if (wopt.eq.2) then
	   den=(sh-mv**2)**2+(gamv*sh/mv)**2
	end if
c
c G_mu-representation:
c
	if (rep.eq.2) then
	   o_fac=(w2*xmw*gfermi)**2/den/3d0	 
c
c alpha_0 representation:
c
	else
	   o_fac=(pi*alpha/sw2)**2/3d0/den
	endif
c
c Born-matrix element:
c
   	call U_Db__e_nu(sig12) 
	gllf=strucf(zc,oxp,llq,betf,alpha0)
	conti(6)=fac*sc_fac*o_fac*(pdf_fac(1)*sig12(1)+
	1    pdf_fac(2)*sig12(2))*gllf*dxp

	if(qfsrexp.eq.1.and.zc.le.dabs(1d0-deltas))then
c
c subtraction of final-state coll. term when multiple final-state
c photon radiation is taken into account:
c
	   conti(6)=conti(6)-
	1	fac*sc_fac*o_fac*(pdf_fac(1)*sig12(1)+
	1	pdf_fac(2)*sig12(2))*alpha/pi/2d0*(llq-1d0)*
c	2	(1d0+zc**2)/(1d0-zc)
	2	(1d0+zc**2)*dxp*(1d0-zc)**(-betf)
	   
	endif
	kern_exp=conti(6)
c
c fill histograms and calculation of sub-contributions: 
c
 5	if(Ihist.eq.1)then
c	   
c recalculation of the invariants with smearing
c
	   if (test(5).ge.1) then
	      sinv(1,2)=2d0*dot_4(b(1,1),b(1,2))
	      do j=1,npart
		 sinv(1,j+2)=2d0*dot_4(b(1,1),p(1,j))
		 sinv(2,j+2)=2d0*dot_4(b(1,2),p(1,j))
	      end do	  
	      do i=1,npart
		 do j=i,npart
		    sinv(i+2,j+2)=2d0*dot_4(p(1,i),p(1,j))
		 end do
	      end do
	      do i=1,npart+2
		 do j=i+1,npart+2
		    sinv(j,i)=sinv(i,j)
		 end do
	      end do
	   end if

	   do i=1,6
	      conti_tot(i)=conti_tot(i)+wgt*conti(i)/it
	   end do

	   call graphs_b(kern_b,xv,xk,wgt)
	   call graphs(kern_b+kern_w+kern_s+kern_c+kern_exp,xv,xk,wgt)
	end if

 10	continue

	if(test(9).eq.0)goto 20
	if(test(7).eq.1)goto 20
c
c calculate the real, hard corrections:
c
	npart=3 
c
c generate the configuration
c
c convention for generation of phase space:
c
c       1-> neutrino or electron,muon
c       2-> electron,muon
c       3-> photon (or gluon)
c
	mass(1)=mff
	mass(2)=mffs
	mass(3)=0d0
	call setup_config(10,npart,mass)
c
c generate the phase space
c
	call config(npart,ptmin2,s,x12,x,fac)
c
c switch to basic convention:
c
c   p1-> electron,muon
c   p2-> neutrino or electron,muon
c
	do i=1,4
	   temp=p(i,1)
	   p(i,1)=p(i,2)
	   p(i,2)=temp
c	   pthree(i,1) = p(i,1)
c	   pthree(i,2) = p(i,2)
	end do
c
c calculate the invariants :
c
	sinv(1,2)=2d0*dot_4(b(1,1),b(1,2))
	do j=1,npart
	   sinv(1,j+2)=2d0*dot_4(b(1,1),p(1,j))
	   sinv(2,j+2)=2d0*dot_4(b(1,2),p(1,j))
	end do	  
	do i=1,npart
	   do j=i,npart
	      sinv(i+2,j+2)=2d0*dot_4(p(1,i),p(1,j))
	   end do
	end do
	do i=1,npart+2
	   do j=i+1,npart+2
	      sinv(j,i)=sinv(i,j)
	   end do
	end do
c momenta for comparison with Madgraph
	ps1(0)=b(4,1)
	ps2(0)=b(4,2)
	ps3(0)=p(4,2)
	ps4(0)=p(4,1)
	ps5(0)=p(4,3)
	do i=1,3
	   ps1(i)=b(i,1)
	   ps2(i)=b(i,2)
	   ps3(i)=p(i,2)
	   ps4(i)=p(i,1)
	   ps5(i)=p(i,3)
	end do

	if(qqcd.eq.2.or.qqg.eq.2) goto 37
c
c make cuts:
c
	flagqcd=0
	nboost=0
	call newcuts(npart,icut)
	   
	if(Icut.eq.1.or.fac.eq.0d0)then
	   ncut_r(1)=ncut_r(1)+1
	   if(Ihist.eq.1)ncut_r(2)=ncut_r(2)+1
	   kern_r = 0d0
	   conti(10) = 0d0
	   goto 37
	end if 

c only parton level result:
	if(test(10).eq.1)then
c CC: only u d initial state
	   if(switch.eq.1)then
	      VKM(2,5)=0d0
	      VKM(5,2)=0d0
	      VKM(4,5)=0d0
	      VKM(5,4)=0d0
	      pdf_fac(1)=1d0
	      pdf_fac(2)=0d0
	   elseif(switch.eq.2)then
c NC: only u ubar initial state
	      pdfu_fac(1)=1d0
	      pdfd_fac(1)=0d0
	      pdfu_fac(2)=0d0
	      pdfd_fac(2)=0d0
	   endif
	   goto 999
	endif

c scale for pdfs:

	if(test(8).lt.3) then
	   qsq=mu_f**2
	elseif(test(8).eq.3)then
c variable scale
	   qsq=2d0*dot_4(sp(1,1),sp(1,2))
	   mu_f=dsqrt(dabs(qsq))
	   mu_r=mu_f
	   alphas=alphas2(mu_r)
	endif
	scale = dsqrt(dabs(qsq))
c	if (scale.lt.1.3d0) then
c	   kern_r = 0d0
c	   conti(10) = 0d0
c	   goto 37
c	end if 
c
c calculate the distribution functions for p and pbar, the conjugation
c for pbar is done here
c ppswitch=1: ppbar, ppswitch=2: pp
c pdf(i=2,4): u,c-distribution (=u_valence+u_sea)
c pdf(i=-2,-4): ubar,cbar-distribution (=u_sea)
c pdf(i=1,3,5): d,s,b-distribution
c
	if(ipdf.eq.1)then
	   call mrstqed(x1,scale,1,upv1,dnv1,usea1,dsea1,str1,chm1,bot1,
	1	glu1,phot1)
	   call mrstqed(x2,scale,1,upv2,dnv2,usea2,dsea2,str2,chm2,bot2,
	1	glu2,phot2)
	   dxpdf1(0)=glu1
	   dxpdf2(0)=glu2
	   dxpdf1(1)=dnv1+dsea1
	   dxpdf1(-1)=dsea1
	   dxpdf1(3)=str1
	   dxpdf1(-3)=str1
	   dxpdf1(5)=bot1
	   dxpdf1(-5)=bot1
	   dxpdf1(2)=upv1+usea1
	   dxpdf1(-2)=usea1
	   dxpdf1(4)=chm1
	   dxpdf1(-4)=chm1
	   dxpdf2(1)=dnv2+dsea2
	   dxpdf2(-1)=dsea2
	   dxpdf2(3)=str2
	   dxpdf2(-3)=str2
	   dxpdf2(5)=bot2
	   dxpdf2(-5)=bot2
	   dxpdf2(2)=upv2+usea2
	   dxpdf2(-2)=usea2
	   dxpdf2(4)=chm2
	   dxpdf2(-4)=chm2
	   dxpdf1(6)=0d0
	   dxpdf1(-6)=0d0
	   dxpdf2(6)=0d0
	   dxpdf2(-6)=0d0
	   do i=-nlf,nlf
	      fpdf1(i)=dxpdf1(i)/x1 
	      if (ppswitch.eq.1) fpdf2(i)=dxpdf2(-i)/x2
	      if (ppswitch.eq.2) fpdf2(i)=dxpdf2(i)/x2
	   end do
	else
c LO: Iset=4, NLO: Iset=400
	   Iset=400
c	   if(qqcd.gt.0) Iset=400
	   Call SetCtq6(Iset)
	   do i=-nlf,nlf
	      dxpdf1(i) = Ctq6Pdf(i,x1,scale)
	      dxpdf2(i) = Ctq6Pdf(i,x2,scale)
	   enddo
	   dxpdf1(1) = Ctq6Pdf(2,x1,scale)
	   dxpdf2(1) = Ctq6Pdf(2,x2,scale)
	   dxpdf1(-1) = Ctq6Pdf(-2,x1,scale)
	   dxpdf2(-1) = Ctq6Pdf(-2,x2,scale)
	   dxpdf1(2) = Ctq6Pdf(1,x1,scale)
	   dxpdf2(2) = Ctq6Pdf(1,x2,scale)
	   dxpdf1(-2) = Ctq6Pdf(-1,x1,scale)
	   dxpdf2(-2) = Ctq6Pdf(-1,x2,scale)
	   do i=-nlf,nlf
	      fpdf1(i)=dxpdf1(i)
	      if (ppswitch.eq.1) fpdf2(i)=dxpdf2(-i)
	      if (ppswitch.eq.2) fpdf2(i)=dxpdf2(i)
	   end do
	endif
c CC:
	if(switch.eq.1)then
c
c test(1) = 1:
c distribution functions factors for u db -> W+ -> l+ nu production
c test(1) = 2:
c distribution functions factors for ub d -> W- -> l- nu production
c
	   pdf_fac(1)=0d0
	   pdf_fac(2)=0d0
	   do i=1,nlf,2
	      do j=2,nlf,2
		 if (test(1).eq.1) then
		    pdf_fac(1)=pdf_fac(1) + fpdf1(j)*fpdf2(-i)*VKM(i,j)**2
		    pdf_fac(2)=pdf_fac(2) + fpdf1(-i)*fpdf2(j)*VKM(i,j)**2
		 elseif (test(1).eq.2) then
		    pdf_fac(2)=pdf_fac(2) + fpdf1(i)*fpdf2(-j)*VKM(i,j)**2
		    pdf_fac(1)=pdf_fac(1) + fpdf1(-j)*fpdf2(i)*VKM(i,j)**2
		 end if
	      end do
	   end do
	elseif(switch.eq.2)then
c
c distribution function factors for 
c u ubar -> Z -> l- l+ production
c
	   pdfu_fac(1)=0d0
	   pdfu_fac(2)=0d0
	   do i=2,nlf,2
	      pdfu_fac(1)=pdfu_fac(1) + fpdf1(i)*fpdf2(-i)
	      pdfu_fac(2)=pdfu_fac(2) + fpdf1(-i)*fpdf2(i)
	   end do
c
c distribution functions factors for 
c d dbar -> Z -> l- l+ production
c
	   pdfd_fac(1)=0d0
	   pdfd_fac(2)=0d0
	   do i=1,nlf,2
	      pdfd_fac(1)=pdfd_fac(1) + fpdf1(i)*fpdf2(-i)
	      pdfd_fac(2)=pdfd_fac(2) + fpdf1(-i)*fpdf2(i)
	   end do
	endif
 999	continue
c
c couplings for the NC process:
c up
	call prop(7,xm,qiu,viu,aiu,nci)
c down
	call prop(10,xm,qid,vid,aid,nci)
c electron or muon
	if(test(2).eq.1) call prop(4,xm,qf,vf,af,ncf)
	if(test(2).eq.2) call prop(5,xm,qf,vf,af,ncf)
c
c cut to avoid soft singularities:
c
	epmin = deltas*dsqrt(sinv(1,2))/2.d0
c
c cut to avoid collinear singularities:
c	
	smin = deltac*dsqrt(sinv(1,2))*ephoton

	if(switch.eq.1) then
c CC: 
	   if(test(3).lt.4)then

	      if (test(3).eq.2.and.test(4).eq.2) epmin = deltas*mw/2d0

	      if(ephoton.lt.epmin)then
		 conti(10)=0d0
		 kern_r = 0d0
		 ncut_r(1)=ncut_r(1)+1
		 if(ihist.eq.1)ncut_r(2)=ncut_r(2)+1
		 goto 37
	      end if
c cut on quark - photon angle:
	      if(test(3).eq.1)then
		 if(dabs(sinv(1,5)).lt.smin.or.dabs(sinv(2,5)).lt.smin)then
		    conti(10)=0d0
		    kern_r = 0d0
		    ncut_r(1)=ncut_r(1)+1
		    if(ihist.eq.1)ncut_r(2)=ncut_r(2)+1
		    goto 37
		 end if
	      endif

	      if(test(3).eq.2.and.collcut.eq.1)then
c cut on lepton - photon angle:
		 smin = smin*(1d0-2d0*ephoton/dsqrt(sinv(1,2)))
		 if(dabs(sinv(3,5)).lt.smin)then
		    conti(10)=0d0
		    kern_r = 0d0
		    ncut_r(1)=ncut_r(1)+1
		    if(ihist.eq.1)ncut_r(2)=ncut_r(2)+1
		    goto 37
		 end if
	      endif
c
c matrix element squared for real, hard photon radiation:
c       
	      call mathard_cc(test(3),sigh12)
	      kern_r=sc_fac*fac*(pdf_fac(1)*sigh12(1)
	1	   +pdf_fac(2)*sigh12(2))
	   elseif(test(3).eq.4)then
c
c summation of all three contributions: initial+final+interference:
c
c	      if(ephoton.lt.epmin)then
c		 conti(10)=0d0
c		 kern_r = 0d0
c		 ncut_r(1)=ncut_r(1)+1
c		 if(ihist.eq.1)ncut_r(2)=ncut_r(2)+1
c		 goto 37
c	      end if
c cut on quark - photon angle:
c	      if(dabs(sinv(1,5)).lt.smin.or.dabs(sinv(2,5)).lt.smin)then
c		 conti(10)=0d0
c		 kern_r = 0d0
c		 ncut_r(1)=ncut_r(1)+1
c		 if(ihist.eq.1)ncut_r(2)=ncut_r(2)+1
c		 goto 37
c	      end if
c	      if(collcut.eq.1)then
c cut on lepton - photon angle:
c		 smin = smin*(1d0-2d0*ephoton/dsqrt(sinv(1,2)))
c		 if(dabs(sinv(3,5)).lt.smin)then
c		    conti(10)=0d0
c		    kern_r = 0d0
c		    ncut_r(1)=ncut_r(1)+1
c		    if(ihist.eq.1)ncut_r(2)=ncut_r(2)+1
c		    goto 37
c		 end if
c	      endif
c	      sigh12(1)=0d0
c	      sigh12(2)=0d0
c	      do k=1,3
c		 call mathard_cc(k,sig12)
c		 sigh12(1)=sigh12(1)+sig12(1)
c		 sigh12(2)=sigh12(2)+sig12(2)
c	      enddo
c	      kern_r=sc_fac*fac*(pdf_fac(1)*sigh12(1)
c	1	   +pdf_fac(2)*sigh12(2))
*
	      kern_r=0d0
	      do k = 1,3
		 icutk(k) = 1
c
c cuts to avoid soft singularities:
c
		 if(ephoton.lt.epmin)then
		    icutk(k) = 0
		    goto 18
		 end if
c
c  cuts to avoid collinear singularities:
c	
		 if(k.eq.1) then
		    if(dabs(sinv(1,5)).lt.smin.or.dabs(sinv(2,5)).lt.smin)then
		       icutk(k) = 0
		       goto 18
		    end if
		 end if
		 if(k.eq.2.and.collcut.eq.1) then
		    smin = smin*(1d0-2d0*ephoton/dsqrt(sinv(1,2)))
		    if(dabs(sinv(3,5)).lt.smin)then
		       icutk(k) = 0
		       goto 18
		    end if
		 end if
		 call mathard_cc(k,sigh12)
		 kern_r=kern_r+sc_fac*fac*(pdf_fac(1)*sigh12(1)
	1	      +pdf_fac(2)*sigh12(2))

 18		 continue
	      end do
	      if ((icutk(1)+icutk(2)+icutk(3)).eq.0) then
		 conti(10)=0d0
		 kern_r=0d0
		 ncut_r(1)=ncut_r(1)+1
		 if(ihist.eq.1)ncut_r(2)=ncut_r(2)+1
		 goto 37
	      endif
	   endif
	elseif(switch.eq.2) then
c NC: 
	   if(ephoton.lt.epmin)then
	      conti(10)=0d0
	      kern_r = 0d0
	      ncut_r(1)=ncut_r(1)+1
	      if(ihist.eq.1)ncut_r(2)=ncut_r(2)+1
	      goto 37
	   end if
c cut on quark - photon angle:
	   if(dabs(sinv(1,5)).lt.smin.or.dabs(sinv(2,5)).lt.smin)then
	      conti(10)=0d0
	      kern_r = 0d0
	      ncut_r(1)=ncut_r(1)+1
	      if(ihist.eq.1)ncut_r(2)=ncut_r(2)+1
	      goto 37
	   end if

	   if(collcut.eq.1)then
c cut on lepton - photon angle:
	      smin = smin*(1d0-2d0*ephoton/dsqrt(sinv(1,2)))
	      if(dabs(sinv(3,5)).lt.smin.or.dabs(sinv(4,5)).lt.smin)then
		 conti(10)=0d0
		 kern_r = 0d0
		 ncut_r(1)=ncut_r(1)+1
		 if(ihist.eq.1)ncut_r(2)=ncut_r(2)+1
		 goto 37
	      end if
	   endif

	   do i=1,2
	      sighu(i)=0d0
	      sighd(i)=0d0
	   enddo
	   do k=1,3
	      call mathard_nc(k,viu,aiu,qiu,vf,af,qf,sig12)
	      sighu(1)=sighu(1)+sig12(1)
	      sighu(2)=sighu(2)+sig12(2)
	      call mathard_nc(k,vid,aid,qid,vf,af,qf,sig12)
	      sighd(1)=sighd(1)+sig12(1)
	      sighd(2)=sighd(2)+sig12(2)
	   enddo

	   kern_r=fac*sc_fac*(pdfu_fac(1)*sighu(1)+
	1	pdfu_fac(2)*sighu(2)+pdfd_fac(1)*sighd(1)+
	2	pdfd_fac(2)*sighd(2))

	endif
	conti(10)=kern_r

 37	continue
	if(qqcd.eq.0) goto 15
c
c QCD contribution
c
c make cuts:
c
	flagqcd=1
	nboost=0
	call newcuts(npart,icut)
	   
	if(Icut.eq.1.or.fac.eq.0d0)then
	   ncut_r(1)=ncut_r(1)+1
	   if(Ihist.eq.1)ncut_r(2)=ncut_r(2)+1
	   kern_qcd = 0d0
	   conti(11)=0d0
	   goto 15
	end if 

c only parton level result:
	if(test(10).eq.1)then
c CC: only u d initial state
	   if(switch.eq.1)then
	      VKM(2,5)=0d0
	      VKM(5,2)=0d0
	      VKM(4,5)=0d0
	      VKM(5,4)=0d0
	      pdf_fac(1)=1d0
	      pdf_fac(2)=0d0
	   elseif(switch.eq.2)then
c NC: only u ubar initial state
	      pdfu_fac(1)=1d0
	      pdfd_fac(1)=0d0
	      pdfu_fac(2)=0d0
	      pdfd_fac(2)=0d0
	   endif
	   goto 888
	endif

c scale for pdfs:
	if(test(8).lt.3) then
	   qsq=mu_f**2
	elseif(test(8).eq.3)then
c variable scale
	   qsq=2d0*dot_4(sp(1,1),sp(1,2))
	   mu_f=dsqrt(dabs(qsq))
	   mu_r=mu_f
	   alphas=alphas2(mu_r)
	endif
	scale = dsqrt(dabs(qsq))
c
c calculate the distribution functions for p and pbar, the conjugation
c for pbar is done here
c ppswitch=1: ppbar, ppswitch=2: pp
c pdf(i=2,4): u,c-distribution (=u_valence+u_sea)
c pdf(i=-2,-4): ubar,cbar-distribution (=u_sea)
c pdf(i=1,3,5): d,s,b-distribution
c
	if(ipdf.eq.1)then
	   call mrstqed(x1,scale,1,upv1,dnv1,usea1,dsea1,str1,chm1,bot1,
	1	glu1,phot1)
	   call mrstqed(x2,scale,1,upv2,dnv2,usea2,dsea2,str2,chm2,bot2,
	1	glu2,phot2)
	   dxpdf1(0)=glu1
	   dxpdf2(0)=glu2
	   dxpdf1(1)=dnv1+dsea1
	   dxpdf1(-1)=dsea1
	   dxpdf1(3)=str1
	   dxpdf1(-3)=str1
	   dxpdf1(5)=bot1
	   dxpdf1(-5)=bot1
	   dxpdf1(2)=upv1+usea1
	   dxpdf1(-2)=usea1
	   dxpdf1(4)=chm1
	   dxpdf1(-4)=chm1
	   dxpdf2(1)=dnv2+dsea2
	   dxpdf2(-1)=dsea2
	   dxpdf2(3)=str2
	   dxpdf2(-3)=str2
	   dxpdf2(5)=bot2
	   dxpdf2(-5)=bot2
	   dxpdf2(2)=upv2+usea2
	   dxpdf2(-2)=usea2
	   dxpdf2(4)=chm2
	   dxpdf2(-4)=chm2
	   dxpdf1(6)=0d0
	   dxpdf1(-6)=0d0
	   dxpdf2(6)=0d0
	   dxpdf2(-6)=0d0
	   do i=-nlf,nlf
	      fpdf1(i)=dxpdf1(i)/x1 
	      if (ppswitch.eq.1) fpdf2(i)=dxpdf2(-i)/x2
	      if (ppswitch.eq.2) fpdf2(i)=dxpdf2(i)/x2
	   end do
	else
c LO: Iset=4, NLO: Iset=400
	   Iset=400
	   Call SetCtq6(Iset)
	   do i=-nlf,nlf
	      dxpdf1(i) = Ctq6Pdf(i,x1,scale)
	      dxpdf2(i) = Ctq6Pdf(i,x2,scale)
	   enddo
	   dxpdf1(1) = Ctq6Pdf(2,x1,scale)
	   dxpdf2(1) = Ctq6Pdf(2,x2,scale)
	   dxpdf1(-1) = Ctq6Pdf(-2,x1,scale)
	   dxpdf2(-1) = Ctq6Pdf(-2,x2,scale)
	   dxpdf1(2) = Ctq6Pdf(1,x1,scale)
	   dxpdf2(2) = Ctq6Pdf(1,x2,scale)
	   dxpdf1(-2) = Ctq6Pdf(-1,x1,scale)
	   dxpdf2(-2) = Ctq6Pdf(-1,x2,scale)
	   do i=-nlf,nlf
	      fpdf1(i)=dxpdf1(i)
	      if (ppswitch.eq.1) fpdf2(i)=dxpdf2(-i)
	      if (ppswitch.eq.2) fpdf2(i)=dxpdf2(i)
	   end do
	endif
c CC:
	if(switch.eq.1)then
c
c test(1) = 1:
c distribution functions factors for u db -> W+ -> l+ nu production
c test(1) = 2:
c distribution functions factors for ub d -> W- -> l- nu production
c
	   pdf_fac(1)=0d0
	   pdf_fac(2)=0d0
	   do i=1,nlf,2
	      do j=2,nlf,2
		 if (test(1).eq.1) then
		    pdf_fac(1)=pdf_fac(1) + fpdf1(j)*fpdf2(-i)*VKM(i,j)**2
		    pdf_fac(2)=pdf_fac(2) + fpdf1(-i)*fpdf2(j)*VKM(i,j)**2
		 elseif (test(1).eq.2) then
		    pdf_fac(2)=pdf_fac(2) + fpdf1(i)*fpdf2(-j)*VKM(i,j)**2
		    pdf_fac(1)=pdf_fac(1) + fpdf1(-j)*fpdf2(i)*VKM(i,j)**2
		 end if
	      end do
	   end do
c PDFs for qg initial state
c CKM squared sum up to one
	   pdfug_fac(1)=0.d0
	   pdfug_fac(2)=0.d0
	   pdfdg_fac(1)=0.d0
	   pdfdg_fac(2)=0.d0
	   if(test(1).eq.1)then
	      do i=1,nlf,2
		 pdfdg_fac(1)=pdfdg_fac(1) + fpdf1(-i)*fpdf2(0)
		 pdfdg_fac(2)=pdfdg_fac(2) + fpdf1(0)*fpdf2(-i)
	      end do
	      do i=2,nlf,2
		 pdfug_fac(1)=pdfug_fac(1) + fpdf1(i)*fpdf2(0)
		 pdfug_fac(2)=pdfug_fac(2) + fpdf1(0)*fpdf2(i)
	      end do
	   elseif(test(1).eq.2)then
	      do i=1,nlf,2
		 pdfdg_fac(1)=pdfdg_fac(1) + fpdf1(i)*fpdf2(0)
		 pdfdg_fac(2)=pdfdg_fac(2) + fpdf1(0)*fpdf2(i)
	      end do
	      do i=2,nlf,2
		 pdfug_fac(1)=pdfug_fac(1) + fpdf1(-i)*fpdf2(0)
		 pdfug_fac(2)=pdfug_fac(2) + fpdf1(0)*fpdf2(-i)
	      end do
	   endif
	elseif(switch.eq.2)then
c
c distribution function factors for 
c u ubar -> Z -> l- l+ production
c
	   pdfu_fac(1)=0d0
	   pdfu_fac(2)=0d0
	   do i=2,nlf,2
	      pdfu_fac(1)=pdfu_fac(1) + fpdf1(i)*fpdf2(-i)
	      pdfu_fac(2)=pdfu_fac(2) + fpdf1(-i)*fpdf2(i)
	   end do
c
c distribution functions factors for 
c d dbar -> Z -> l- l+ production
c
	   pdfd_fac(1)=0d0
	   pdfd_fac(2)=0d0
	   do i=1,nlf,2
	      pdfd_fac(1)=pdfd_fac(1) + fpdf1(i)*fpdf2(-i)
	      pdfd_fac(2)=pdfd_fac(2) + fpdf1(-i)*fpdf2(i)
	   end do
	endif
 888	continue
c
c couplings for the NC process:
C up
	call prop(7,xm,qiu,viu,aiu,nci)
c down
	call prop(10,xm,qid,vid,aid,nci)
c electron or muon
	if(test(2).eq.1) call prop(4,xm,qf,vf,af,ncf)
	if(test(2).eq.2) call prop(5,xm,qf,vf,af,ncf)
c
c cut to avoid soft singularities:
c
	epmin = deltas*dsqrt(sinv(1,2))/2.d0
c
c cut to avoid collinear singularities:
c	
	smin = deltac*dsqrt(sinv(1,2))*ephoton

	if(ephoton.lt.epmin.and.qqg.ne.2)then
	   conti(11)=0d0
	   kern_qcd = 0d0
	   ncut_r(1)=ncut_r(1)+1
	   if(ihist.eq.1)ncut_r(2)=ncut_r(2)+1
	   goto 15
	end if
c cut on quark - gluon angle:
	if(dabs(sinv(1,5)).lt.smin.or.dabs(sinv(2,5)).lt.smin)then
	   conti(11)=0d0
	   kern_qcd = 0d0
	   ncut_r(1)=ncut_r(1)+1
	   if(ihist.eq.1)ncut_r(2)=ncut_r(2)+1
	   goto 15
	end if

	if(switch.eq.1) then
c CC: 
c matrix element squared for real, hard gluon radiation:
c
	   if(qqg.eq.1.or.qqg.eq.3)then
	      call mathard_qcd(sigh12)
	      kern_qcd=sc_fac*fac*(pdf_fac(1)*sigh12(1)+
	1	   pdf_fac(2)*sigh12(2))
	   endif
	   if(qqg.eq.2.or.qqg.eq.3)then
	      do i=1,2
		 sighu(i)=0d0
		 sighd(i)=0d0
	      enddo
	      kern_qcd=kern_qcd+
	1	   fac*sc_fac*(pdfug_fac(1)*sighu(1)+
	2	   pdfug_fac(2)*sighu(2)+pdfdg_fac(1)*sighd(1)+
	2	   pdfdg_fac(2)*sighd(2))
	   endif
	elseif(switch.eq.2) then
c NC: 
	   kern_qcd=0d0
c	   kern_qcd=fac*sc_fac*(pdfu_fac(1)*sighu(1)+
c	1	pdfu_fac(2)*sighu(2)+pdfd_fac(1)*sighd(1)+
c	2	pdfd_fac(2)*sighd(2))

	endif

	conti(11)=kern_qcd
c
c calculate the integrated cross section and 
c fill histograms 
c
 15	if(Ihist.eq.1)then
c	   
c recalculation of the invariants with smearing
c
	   if (test(5).ge.1) then

	      sinv(1,2)=2d0*dot_4(b(1,1),b(1,2))

	      do j=1,npart
		 sinv(1,j+2)=2d0*dot_4(b(1,1),p(1,j))
		 sinv(2,j+2)=2d0*dot_4(b(1,2),p(1,j))
	      end do	  
	      do i=1,npart
		 do j=i,npart
		    sinv(i+2,j+2)=2d0*dot_4(p(1,i),p(1,j))
		 end do
	      end do
	      do i=1,npart+2
		 do j=i+1,npart+2
		    sinv(j,i)=sinv(i,j)
		 end do
	      end do

	   end if

	   do i=10,11
	      conti_tot(i)=conti_tot(i)+wgt*conti(i)/It
	   end do

	   call graphs(kern_r+kern_qcd,xv,xk,wgt)
	end if
c
c sum the contributions
c
 20	continue
	if(ihist.ne.1)then 
	  kern=dabs(kern_b)+dabs(kern_w)+dabs(kern_s)+dabs(kern_c)
	1	+dabs(kern_exp)+dabs(kern_r)+dabs(kern_qcd)
	else
	   kern=kern_b+kern_w+kern_s+kern_c+kern_r+kern_exp+kern_qcd
	end if	  

	if(kern.eq.0.d0)ncut(1)=ncut(1)+1
	if(kern.eq.0.d0.and.Ihist.eq.1)ncut(2)=ncut(2)+1
	return
	end
c
c definition of the scalar product of four vectors:
c
	real*8 function dot_4(p1,p2)
	implicit none
	real*8 p1(4),p2(4)
	dot_4=-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)+p1(4)*p2(4)
	return
	end
c-------------------------------------------------------------
c Born-matrix element squared for Z,gamma production:
c (all fermions are considered to be massless !)
c-------------------------------------------------------------
	real*8 function a2qqz(xs,xt,xu,xvf,xaf,xqf,xvi,xai,xqi,
	1    alz,alg)
	implicit none
	real*8 sigz,sigg,siggz
	real*8 alaq,blbq,t1,t2,xs,xt,xu
	real*8 xvf,xaf,xqf,xvi,xai,xqi
	real*8 den,alz,alg
	real*8 rs,mw,mz,xmz,xmw,gamz,gamw
	common/par_prop/den
	common/par_gen/rs,mw,mz,xmw,xmz,gamz,gamw
	t1=xt
	t2=xu
	alaq=(xaf**2+xvf**2)*(xai**2+xvi**2)
	blbq=-4d0*xaf*xvf*xai*xvi

	sigz=8.d0*(alaq*(t1**2+t2**2)+blbq*(t1**2-t2**2))
	sigg=8.d0*(t1**2+t2**2)
	siggz=8.d0*(2d0*xvf*xvi*(t1**2+t2**2)-
	1    2d0*xaf*xai*(t1**2-t2**2))
	a2qqz=16d0*(alz**2*sigz/den+
	1    alg**2*(xqf*xqi)**2*sigg/xs**2+
	2    alz*alg*xqf*xqi*(xs-xmz)/den/xs*siggz)
	return
	end
c
c Born-matrix element for W production:
c
c       I (Stephane) followed Phys. Rev. D. 43, 2892, Baer and Reno; and
c       Nucl. Phys. B185, 274., Aurenche and Lindfors.
c       There is a k^2 missing in the second paper.  I checked
c       the first one by comparing the matrix elements with the 
c       appropriate crossing.
c (all fermions are considered to be massless !)
c
	subroutine U_Db__e_nu(sig)
	implicit none
	real*8 sig(2)
	real*8 alaq,blbq,t1,t2
	include 'config.inc'
	alaq=1.d0/4.d0
	blbq=1.d0/4.d0
	t1=-sinv(1,3)
	t2=-sinv(2,3)
	sig(1)=8.d0*(alaq*(t1**2+t2**2)+blbq*(t1**2-t2**2))
	t1=-sinv(2,3)
	t2=-sinv(1,3)
	sig(2)=8.d0*(alaq*(t1**2+t2**2)+blbq*(t1**2-t2**2))
	return
	end

c
c structure function for multiple photon radiation:
c
	double precision function strucf(xx,dxx,l,bet,xalp)
	implicit none
	real*8 xx,l,bet,xalp,pi
	real*8 ge,gamma,dgammf,gamf,lxx,dxx
	complex*16 spence,ieps
	ieps=dcmplx(0d0,1d-20)
	pi=4d0*datan(1d0)
	ge=0.5772115664901532d0
c	gamma=dgammf(1d0+bet)
	gamma=gamf(1d0+bet)
	dxx=1d0-xx
	if (xx.lt.0d0.or.dxx.lt.0d0.or.xx.lt.1d-5) then
	   write(6,*)'stop in strucf',xx,dxx
	   stop
	endif
	if (dabs(dxx).lt.1d-8) then
	   lxx=-1d0-dxx/2.d0-dxx**2/3.d0-dxx**3/4.d0
c       lxx=-1d0
	else
	   lxx=dlog(xx)/(dxx)
	end if
c	strucf=dexp(bet*(-ge+3d0/4d0))/gamma*bet
c	1    +(-bet/2d0*(1d0+xx)
c	2    -bet**2/8d0*((1d0+3d0*xx*xx)*lxx
c	3    +4d0*(1d0+xx)*dlog(dxx)+5d0+xx))
c	1    *(dxx)**(1d0-bet)
c
	strucf=dexp(bet*(-ge+3d0/4d0))/gamma*bet
	1    +(-bet/2d0*(1d0+xx)
	2    -bet**2/8d0*((1d0+3d0*xx*xx)*lxx
	3    +4d0*(1d0+xx)*dlog(dxx)+5d0+xx)
	4    -bet**3/48d0*((1d0+xx)*(6d0*dreal(spence(xx+ieps))
	5    +12d0*dlog(dxx)**2-3d0*pi**2)
	6    +3d0/2d0*(1d0+8d0*xx+3d0*xx**2)*lxx+6d0*(xx+5d0)*dlog(dxx)
	7    +12d0*(1d0+xx**2)*dlog(dxx)*lxx
	8    -(1d0+7d0*xx**2)/2d0*lxx*dlog(xx)
	9    +(-24d0*xx-15d0*xx**2+39d0)/dxx/4d0))*(dxx)**(1d0-bet)

c	strucf=dexp(bet*(-ge+3d0/4d0))/gamma*bet*dxx**(bet-1d0)
c	1    +(-bet/2d0*(1d0+xx)
c	2    -bet**2/8d0*((1d0+3d0*xx*xx)*lxx
c	3    +4d0*(1d0+xx)*dlog(dxx)+5d0+xx)
c	4    -bet**3/48d0*((1d0+xx)*(6d0*dreal(spence(xx+ieps))
c	5    +12d0*dlog(dxx)**2-3d0*pi**2)
c	6    +3d0/2d0*(1d0+8d0*xx+3d0*xx**2)*lxx+6d0*(xx+5d0)*dlog(dxx)
c	7    +12d0*(1d0+xx**2)*dlog(dxx)*lxx
c	8    -(1d0+7d0*xx**2)/2d0*lxx*dlog(xx)
c	9    +15d0*xx+39d0))

c test: YFS exponentiation:
c       strucf=dexp(-bet*ge)/gamma*bet
c       $     *dexp(xalp/pi*(l/2d0-1d0+pi**2/3d0))
c       $     *(1d0+bet/2d0+(xalp/pi*l)**2/2d0
c       $     -(1d0-xx**2)/2d0+xalp/pi*l*(-(1d0+3d0*xx**2)/4d0*dlog(xx)
c       $     -1d0+xx))
	return
	end
c**********************************************************************
c gamma function
c code written by M.Roth
c**********************************************************************
	real*8 function gamf(x)
	implicit none
c local variable
	integer max
	parameter(max=10)
	real*8 x,omx,a(max),prod,pi
	integer i1
	pi=4d0*datan(1d0)
c a(1)=gammae
	data a/0.5772156649015328606d0,0d0,
     *       0.4006856343865314285d0,0d0,
     *       0.2073855510286739853d0,0d0,
     *       0.1440498967688461181d0,0d0,
     *       0.1113342658695646905d0,0d0/ 
c good convergence for small values of 1-x
	omx=x-1d0
	gamf=dexp(0.5d0*dlog(omx*pi/dsin(pi*omx)))
	do i1=1,max,2
	   prod=dexp(-a(i1)*omx**i1)
	   gamf=gamf*prod
	   if(dabs(prod-1d0).lt.1d-10)return
	enddo
	write(*,'(a)')' gamma: converges badly'
	stop
	return
	end
