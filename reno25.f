	subroutine reno(itmx,ncall,ndim,fxn,avgi,sd,chi2a,ierr,n_proc,my_id)
c
c   subroutine performs ndim-dimensional monte carlo integ'n
c	- by g.p. lepage    sept 1976/(rev)aug 1979
c	- algorithm described in j comp phys 27,192(1978)
c
c modifications for parallelization based on R. Kreckel's PVEGAS algorithm
c PVEGAS algorithm described in arXiv:physics/9710028v1 [physics.comp-ph]
c
	implicit real*8 (a-h,o-z)
	integer*4 ncall,npg,nd,ng
c
c new variables for parallel version
c
	integer ierr,n_proc,my_id
c
c
c
c status field for RNG
c
	integer*4 gfsr_m(1279)
c
	integer*4 gfsr_q
c
c current time will be used as the intial seed for RNG
c
	real*8 time
c
c flag to show if RNG has been initialized
c
	integer*4 gfsr_initialized
c
c counter for RNG
c
	integer*4 rdum
c
c norm for RNG
c
	real*8 gfsr_norm
	real*8 rand_num
c
c
c
c	common/bveg1/ncall,itmx,nprn,ndev,xl(25),xu(25),acc
	common/bveg1/nprn,ndev,xl(25),xu(25),acc
	common/bveg2/it,ndo,si,swgt,schi,xi(50,25)
	common/bveg3/alph,ndmx,mds
	common/bveg4/calls,ti,tsi
	dimension d(50,25),di(50,25),xin(50),r(50),dx(25),ia(25),
	1    kg(25),dt(25),x(25)
	data nprn/0/,acc/-1.d0/,
	1    xl/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
	1    0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
	1    0.d0,0.d0,0.d0,0.d0,0.d0/,
	1    xu/1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,
	1    1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,
	1    1.d0,1.d0,1.d0,1.d0,1.d0/,
	1    alph/1.5d0/,ndmx/50/,mds/1/,ndev/6/,
	1    ndo/1/,xi/1250*1.d0/,it/0/,si,swgt,schi/3*0.d0/
	real rand(25)
c   factor for error w/ 90 % confidence level
	data factor/1.65d0/
c
c
	data gfsr_initialized/0/
c
c
c       common block for interaction with integration, 
c
	real*8 stoti(50),sdeli(50)
	common/par_integration_reno/stoti,sdeli
c
c
c
	time=abs(MPI_WTIME())
c
c
	if (gfsr_initialized .EQ.0)then
	call   gfsr_init(time,gfsr_m,gfsr_q,gfsr_initialized,rdum,my_id)
	endif
	
c
	ndo=1
	do 6901 j=1,ndim
	   xi(1,j)=1.d0
 6901	continue
c
	entry reno1(itmx,ncall,ndim,fxn,avgi,sd,chi2a,ierr,n_proc,my_id)
c       - initializes cummulative variables, but not grid
c
c
	it=0
	si=0.d0
	swgt=0.d0
	schi=0.d0
c
	entry reno2(itmx,ncall,ndim,fxn,avgi,sd,chi2a,ierr,n_proc,my_id)
c	   - no initialization
c
c
	nd=ndmx
	ng=1
	if(mds.ne.0) then
	   ng=(ncall/2.d0)**(1.d0/ndim)
	   mds=1
	   if((2*ng-ndmx).ge.0) then
	      mds=-1
	      npg=ng/ndmx+1
	      nd=ng/npg
	      ng=npg*nd
	   end if
	end if
	k=ng**ndim
	npg=ncall/k
	if(npg.lt.2) npg=2
	calls=npg*k
	dxg=1.d0/ng
	dv2g=(calls*dxg**ndim)**2/dble(npg)/dble(npg)/dble(npg-1)
	xnd=nd
	ndm=nd-1
	dxg=dxg*xnd
	xjac=1.d0/calls
	do  6902 j=1,ndim
	   dx(j)=xu(j)-xl(j)
	   xjac=xjac*dx(j)
 6902	continue
c
c   rebin, preserving bin density
c
	if(nd.ne.ndo) then
	   rc=ndo/xnd
	   do 6903 j=1,ndim
	      xo=0.d0
	      xn=xi(1,j)
	      k=1
	      dr=1.d0
	      do 6904 i =1,ndm
 6905		 if  (rc.gt.dr) then
		    k=k+1
		    dr=dr+1.d0
		    xo=xn
		    xn=xi(k,j)
		    go to 6905
		 end if
		 dr=dr-rc
		 xin(i)=xn-(xn-xo)*dr
 6904	      continue
	      do  6906 i=1,ndm
		 xi(i,j)=xin(i)
 6906	      continue
	      xi(nd,j)=1.d0
 6903	   continue
	   ndo=nd
	end if
c
c	if(nprn.ge.0) write(ndev,200) ndim,calls,it,itmx,acc,nprn,
c	1		  alph,mds,nd,(xl(j),xu(j),j=1,ndim)
c
	entry reno3(itmx,ncall,ndim,fxn,avgi,sd,chi2a,ierr,n_proc,my_id)
c       - main integration loop
c
c
9	it=it+1
c
c	modification of alph
c
	alphp=alph*(1.d0-real(it)/real(itmx))
c	print*,' alph=',alphp
	ti=0.d0
	tsi=0.d0
	do  6911 j=1,ndim
	   kg(j)=1
	   do 6912 i=1,nd
	      d(i,j)=0.d0
	      di(i,j)=0.d0
 6912	   continue
 6911	continue
c
11	fb=0.d0
	f2b=fb
	do 6913 k=1,npg
	   wgt=xjac
	   do 6914 j=1,ndim
	      call gfsr_rand(gfsr_m,gfsr_q,rand_num)
	      xn=(kg(j)-rand_num)*dxg+1.d0
	      ia(j)=xn
	      if(ia(j).eq.1) then
		 xo=xi(ia(j),j)
		 rc=(xn-ia(j))*xo
	      else
		 xo=xi(ia(j),j)-xi(ia(j)-1,j)
		 rc=xi(ia(j)-1,j)+(xn-ia(j))*xo
	      end if
	      x(j)=xl(j)+rc*dx(j)
	      wgt=wgt*xo*xnd
 6914	   continue
c
	   f=wgt
	   f=f*fxn(x,wgt)
	   f2=f*f
	   fb=fb+f
	   f2b=f2b+f2
	   do 6915 j=1,ndim
	      di(ia(j),j)=di(ia(j),j)+f
	      if(mds.ge.0) d(ia(j),j)=d(ia(j),j)+f2
 6915	   continue
 6913	continue
c       
	f2b=dsqrt(f2b*npg)
	f2b=(f2b-fb)*(f2b+fb)
	ti=ti+fb
	tsi=tsi+f2b
	if(mds.lt.0) then
	   do 6916 j=1,ndim
	      d(ia(j),j)=d(ia(j),j)+f2b
 6916	   continue
	end if
	k=ndim
19	kg(k)=mod(kg(k),ng)+1
	if(kg(k).ne.1) go to 11
	k=k-1
	if(k.gt.0) go to 19
c
c   compute final results for this iteration
c
	tsi=tsi*dv2g
	ti2=ti*ti
	wgt=1.d0/tsi
	si=si+ti*wgt
	swgt=swgt+wgt
	schi=schi+ti2*wgt
	avgi=si/swgt
	chi2a=(schi-si*avgi)/(it-.9999d0)
	sd=sqrt(1.d0/swgt)
	fnlerr=factor*sd
c
c
	
	if(nprn.ge.0) then
	   tsi=sqrt(tsi)
	   err=factor*tsi
	   if(my_id=0) then
c           write(ndev,201) it,ti,err,avgi,fnlerr,chi2a
	   write(ndev,201) it,ti,tsi,avgi,sd,chi2a
	   endif
	   stoti(it)=ti
	   sdeli(it)=tsi
	   if(nprn.ne.0) then
	      do 6917 j=1,ndim
c       write(ndev,202) j,(xi(i,j),di(i,j),
c	1	   i=1+nprn/2,nd)
 6917	      continue
	   end if
	end if
	endif
	
	if (itmx.eq.1) return
	
c       
c   refine grid
c
	do 6918 j=1,ndim
	   xo=d(1,j)
	   xn=d(2,j)
	   d(1,j)=(xo+xn)/2.d0
	   dt(j)=d(1,j)
	   do 6919 i=2,ndm
	      d(i,j)=xo+xn
	      xo=xn
	      xn=d(i+1,j)
	      d(i,j)=(d(i,j)+xn)/3.d0
	      dt(j)=dt(j)+d(i,j)
 6919	   continue
	   d(nd,j)=(xn+xo)/2.d0
	   dt(j)=dt(j)+d(nd,j)
 6918	continue
c       
	do  6920 j=1,ndim
	   rc=0.d0
	   do 6921 i=1,nd
	      if(d(i,j).le.0.d0) then
		 r(i)=0.d0
	      else
		 xoln=log(dt(j)) - log(d(i,j))
		 if (xoln.gt.70.d0) then
		    r(i)=(1.d0/xoln)**alphp
		 else
		    xo=exp(xoln)
		    r(i)=((xo-1.d0)/xo/xoln)**alphp
		 end if
		 rc=rc+r(i)
	      end if
 6921	   continue
	   rc=rc/xnd
	   xo=0.d0
	   xn=xi(1,j)
	   k=1
	   dr=r(1)
	   do 6922 i =1,ndm
 6923	      if (rc.gt.dr) then
		 k=k+1
		 dr=dr+r(k)
		 xo=xn
		 xn=xi(k,j)
		 go to 6923
	      end if
	      dr=dr-rc
	      xin(i)=xn-(xn-xo)*dr/r(k)
 6922	   continue
	   do 6924 i=1,ndm
	      xi(i,j)=xin(i)
 6924	   continue  
	   xi(nd,j)=1.d0
 6920	continue
c
	if(it.lt.itmx.and.acc*abs(avgi).lt.sd) go to 9
	return
c
c
	if(my_id==0) then
200	format(/' input parameters for reno:  ndim=',i3,'  ncall=',f8.0
	1    /28x,'  it=',i5,'  itmx=',i5/28x,'  acc=',g9.3
	1    /28x,'  nprn=',i3,'  alph=',f5.2/28x,'  mds=',i3,' nd=',i4
	1    /28x,'  (xl,xu)=',t40,'( ',g12.6,' , ',g12.6,' )')
c       201	format(///' integration by reno'//' iteration no.',i3,
c	1	':   integral =',g14.8/21x,'error    =',g10.4/
c	1	' accumulated results:',24x,'final value =',g14.8/
c	1	48x,'error    =',g10.4/48x,'chi**2 per it''n =',g10.4)
 201	format(i3,' : ',g12.6,' +- ', g10.4,' tot:', g12.6,' +- ',
	1    g8.2, ' chi2:', g10.4)
 202	format(/' data for axis ',i2/'    x       delta i    ',
	1    4('  x      delta i       '),'  x      delta i'
	1    /(1x,f7.6,g11.4,f11.6,g11.4,f11.6,g11.4,
	1    f11.6,g11.4,f11.6,g11.4,f11.6,g11.4))
	endif
	end
c
c       parallel vegas handles RNG differently, so new subroutines have been added, but the old RNG remains because it is used by other subprograms in the build
c
	subroutine randa(n,rand)
c
c   subroutine generates uniformly distributed random no's x(i),i=1,n
	real rand(n)
c
	call ranmar(rand,n)
c
	return
	end
c
c the new RNG does not call an external function
c
	subroutine gfsr_rand(w,q,rand_num)
c
	real*8 rand_num
	integer*4 w(*)
	integer*4 q
	integer*4 j
	real*8 gfsr_norm
c
	q=q+1
	gfsr_norm=-1/dble(2**(8*4-1))
	if (q .GE. 1279) q=0;
	j=q+418
	if (j .GE. 1279) j=j-1279
	w(q)=XOr(w(q),w(j))
	rand_num=abs(dble(w(q))*dble(gfsr_norm))
c	print *,'Value of rand:'
c	print *,rand_num
	return 
	end
c
c
	subroutine kw_rand(rdum)
c
	integer*4 rdum
c
	rdum=(1812433253*rdum+314159265)
	return
	end
c
	subroutine gfsr_init(gfsr_seed,w,gfsr_q,gfsr_initialized,rdum,my_id)
c    
	integer*4 rdum
	integer my_id
	integer*4 w(*)
	integer*4 gfsr_q
	integer*4 gfsr_initialized
	real*8 gfsr_seed
	integer*4 i
c
	rdum=int(gfsr_seed);
	do 30 i = 1, 1279*(my_id+1)
	   call kw_rand(rdum)
 30	continue
c
	do 31 i = 1, 1279
	   call kw_rand(rdum)
	   w(i)=rdum
 31	continue
	gfsr_q=-1
	gfsr_initialized=1
	return
	end
	
