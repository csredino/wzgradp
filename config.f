	subroutine setup_config(iout,npart,mm)
c
c	set up for the configuration generator
c	     iout: unit for any output
c	     npart: number of particles
c	     mm(i): mass of the ith particle
c	calculates the ms(i) and put them in the common block
c	along with the masses
c	
	implicit none
	integer iout,npart
	real*8 mm(npart)
	integer i
	include 'config.inc'
	do i=1,npart
	  m(i)=mm(i)
	end do
	ms(npart)=m(npart)
	do i=npart-1,1,-1
	  ms(i)=ms(i+1)+m(i)
	end do
	return
	end	

	subroutine config(npart,ptmin,s,x12,x,fac)
c
c	generator of phase space configuration in the laboratory for 
c	pp or ppbar colliders
c	10/06/94
c
c	npart=number of particles
c	ptmin=minimum pt of the first particle
c	s=square of the center of mass energy
c	x12(2): random numbers for x1 and x2
c	x(3*npart-4) random numbers 
c       (3 for npart-2, and 2 for the npart-1th)
c	fac: overall weight of the event (in nbarns)
c
	implicit none
	integer npart
	real*8 ptmin,s,x12(2),x(3,npart-1),fac,tau,taumin,etaq
	real*8 dot_4
c
c common blocks:
c
	include 'config.inc'
c
c conversion factor, the answer is in nbarns:
c
	real*8 conv
	parameter(conv=389379.66d0)

	real*8 shmin,xmin,sh,rsh,wps
	integer i
	
	real*8 mcut
	common/masscut/mcut

	integer test(10)
	common/par_test/test
c
c energy and momentum of center of mass in lab for boost:
c
	real*8 ecm,pcm

	x1 = x12(1)
	x2 = x12(2)
	
	fac = 1d0

c enhance threshold region:
	taumin=(m(1)+m(2))**2/s
c	taumin=mcut**2/s
	if (npart.eq.3) taumin=(m(1)+m(2))**2/s
c	if (test(10).eq.0) then
cc         etaqmin=-1d1
cc         etaqmax=dlog10(1d0/taumin-1d0)
cc         etaq=(etaqmax-etaqmin)*x1+etaqmin
cc         tau=taumin*(1d0+10d0**etaq)
cc         x1=(1d0-tau)*x2+tau
cc         x2=tau/x1
cc jacobian factor:
cc         fac=(etaqmax-etaqmin)*taumin*10d0**etaq*dlog(10d0)
cc     $        *(1d0-tau)/x1
c	   tau=(1d0-taumin)*x1+taumin
c	   etaq=-dlog(tau)*x2
c	   x1=dexp(-etaq)
c	   x2=tau/x1
c jacobian factor:
c	   fac=-(1d0-taumin)*dlog(tau)
c	end if

	if (test(10).eq.0) then
	   sh=x1*x2*s
	else if (test(10).eq.1) then
	   sh=200d0**2
	end if
	rsh=dsqrt(sh)
	if (rsh.le.(m(1)+m(2))) then
	   fac=0d0
c	   write(6,*)rsh,m(1)+m(2)
	   return
	end if

c calculation of the four momenta:

c	if(npart.eq.2) call congen_new(npart,rsh,x,wps)
	if(npart.eq.2) call congen_n(npart,rsh,x,wps)
	if(npart.eq.3) call congen_n(npart,rsh,x,wps)
c
c factor:
c the answer is in pbarn                  
	fac=fac*wps*conv/2.d0/sh*1d3
c
c beam momenta
c
	b(1,1)=0.d0
	b(2,1)=0.d0
	b(3,1)=rsh/2.d0
	b(4,1)=rsh/2.d0
	
	b(1,2)=0.d0
	b(2,2)=0.d0
	b(3,2)=-rsh/2.d0
	b(4,2)=rsh/2.d0
c
c Energy of the photon in cm frame needs to be transmitted to kern:
c
	ephoton=p(4,3)
c cos theta(l) in parton cms
c tests:
c        cgl(2)=(b(1,1)*p(1,1)+b(2,1)*p(2,1)+b(3,1)*p(3,1))/
c     $       dsqrt(p(1,1)**2+p(2,1)**2+p(3,1)**2)/
c     $       dsqrt(b(1,1)**2+b(2,1)**2+b(3,1)**2)
c        cgl(1)=(b(1,1)*p(1,2)+b(2,1)*p(2,2)+b(3,1)*p(3,2))/
c     $       dsqrt(p(1,2)**2+p(2,2)**2+p(3,2)**2)/
c     $       dsqrt(b(1,1)**2+b(2,1)**2+b(3,1)**2)
c	cgl(2)=p(3,1)/dsqrt(p(1,1)**2+p(2,1)**2+p(3,1)**2)
c	cgl(1)=p(3,2)/dsqrt(p(1,2)**2+p(2,2)**2+p(3,2)**2)
	cgl(1)=(b(4,1)*p(4,1)-dot_4(b(1,1),p(1,1)))/
     $       dsqrt(p(1,1)**2+p(2,1)**2+p(3,1)**2)/
     $       dsqrt(b(1,1)**2+b(2,1)**2+b(3,1)**2)
        cgl(2)=(b(4,1)*p(4,2)-dot_4(b(1,1),p(1,2)))/
     $       dsqrt(p(1,2)**2+p(2,2)**2+p(3,2)**2)/
     $       dsqrt(b(1,1)**2+b(2,1)**2+b(3,1)**2)
c
c boost to laboratory frame
c
	if (test(10).eq.0) then
	   ecm=(x1+x2)*dsqrt(s)/2.d0
	   pcm=(x1-x2)*dsqrt(s)/2.d0
	   do i=1,2
	      call boost_zz(ecm,pcm,b(1,i))
	   end do	  
	   do i=1,npart
	      call boost_zz(ecm,pcm,p(1,i))
	   end do
	end if
	return
	end

	subroutine boost_zz(ecm,pcm,p)
c
c this subroutine performs a boost along the z direction
c       
	implicit none
	real*8 ecm,pcm,p(4)
	real*8 m,gamma,bgamma,E,Pp

	m=dsqrt(ecm**2-pcm**2)

	gamma = ecm/m		  
	bgamma = pcm/m
	E=p(4)
	Pp=p(3)
	p(4)= gamma*E+bgamma*Pp
	p(3)= gamma*Pp+bgamma*E 

	return
	end

	subroutine invboost_zz(ecm,pcm,p)
c
c this subroutine performs a boost along the z direction
c       
	implicit none
	real*8 ecm,pcm,p(4)
	real*8 m,gamma,bgamma,E,Pp

	m=dsqrt(ecm**2-pcm**2)

	gamma = ecm/m		  
	bgamma = pcm/m
	E=p(4)
	Pp=p(3)
	p(3)= (gamma*Pp-bgamma*E)/(gamma**2-bgamma**2)
	p(4)= (gamma*E-bgamma*Pp)/(gamma**2-bgamma**2)

	return
	end

	subroutine congen_new(npart,rsh,x,wps)
c	
c	inputs:
c		npart: number of particles (10 maximum)
c		rsh: center of mass energy.
c		x(3*npart-4): random numbers.
c		m(npart): masses of the particles.
c	outputs:
c		p(4,npart): momenta of the particles (4 is the energy).
c		wps: weight of the phase space
c
c	D.Wackeroth 07/06/99
c
	implicit none
	integer i,npart
	real*8 rsh,x(3,npart-1),wps
	real*8 alam,pi,p0,pb
	real*8 ct,st,ph

	include 'config.inc'

	pi=4d0*datan(1d0)
	wps=1d0

c generation of invariants and angles:

	ct=2d0*x(1,1)-1d0
	st=dsqrt(dabs((1d0-ct)*(1d0+ct)))
	ph=2d0*pi*x(1,2)

c decay a+b->f+fbar in the CMS vec(a+b)=vec(0) 
c with vec(a) is pos. z-axis:
	p0=(rsh**2+m(1)**2-m(2)**2)/2d0/rsh
	pb=alam(rsh,m(1),m(2))/2d0/rsh

	p(1,1)=pb*st*dcos(ph)
	p(2,1)=pb*st*dsin(ph)
	p(3,1)=pb*ct
	p(4,1)=p0

	p0=(rsh**2+m(2)**2-m(1)**2)/2d0/rsh
	p(1,2)=-pb*st*dcos(ph)
	p(2,2)=-pb*st*dsin(ph)
	p(3,2)=-pb*ct
	p(4,2)=p0

	wps=wps*alam(rsh,m(1),m(2))/8d0/rsh**2*4d0*pi
c
c normalization of the phase space weight:
c
	wps=wps/(2.d0*pi)**2
	return
	end
********************************************************
	subroutine boost_z(ecm,pcm,p4,p3)
********************************************************
c this subroutine performs a boost along the z direction
********************************************************
	implicit none
	real*8 ecm,pcm,p4,p3
	real*8 gamma,bgamma,E,Pp
	gamma=ecm
	bgamma=pcm
	E=p4
	Pp=p3
	p4=gamma*E+bgamma*Pp
	p3=gamma*Pp+bgamma*E 
	return
	end
********************************************************
	subroutine rotation_1(vz,st,ct,p1,p3)
********************************************************
c this subroutine performs one rotation
********************************************************
	implicit none
	real*8 st,ct,p1,p3
	real*8 p1s,p3s,vz
	p1s=p1
	p3s=p3
	p1=vz*(ct*p1s+st*p3s)
	p3=vz*(-st*p1s+ct*p3s)
	return
	end
********************************************************
	subroutine rotation_2(vz,st,ct,sp,cp,p1,p2,p3)
********************************************************
c this subroutine performs two rotations
********************************************************
	implicit none
	real*8 st,ct,sp,cp,p1,p2,p3
	real*8 p1s,p2s,p3s,vz
	p1s=p1
	p2s=p2
	p3s=p3
	p1=vz*(cp*(ct*p1s+st*p3s)-sp*p2s)
	p2=vz*(sp*(ct*p1s+st*p3s)+cp*p2s)
	p3=vz*(-st*p1s+ct*p3s)
	return
	end
********************************************************
	real*8 function alam(m,m1,m2)
	implicit none
	real*8 m,m1,m2
	real*8 s,s1,s2
	real*8 aux
	s=m**2
	s1=m1**2
	s2=m2**2
	aux=s**2+s1**2+s2**2-2*s*s1-2*s*s2-2*s1*s2
	if(aux.lt.0.d0)aux=0.d0	
	alam=dsqrt(aux)
	return
	end	

	subroutine congen_n(npart,rsh,x,wps)
c	
c	inputs:
c		npart: number of particles (10 maximum)
c		rsh: center of mass energy.
c		x(3*npart-4): random numbers.
c		m(npart): masses of the particles.
c	outputs:
c		p(4,npart): momenta of the particles (4 is the energy).
c		wps: weight of the phase space
c
c	S.Keller 07/06/94
c
	implicit none
	
	integer npart
	real*8 rsh,x(3,npart-1),wps
c
c common blocks 
c
	include 'config.inc'

	real*8 mt,pm(3),wpsi,xnpart_1(3)
	integer i
	real*8 pi
	parameter(pi=3.1415926535987932)
                                                     
	mt=rsh                       
	wps=1.d0
	pm(1)=0.d0
	pm(2)=0.d0
	pm(3)=0.d0

	if(npart.eq.2) goto 10
        do i=1,npart-2
	   call congen_2(mt,pm,x(1,i),ms(i),m(i),p(1,i),wpsi)
	   wps=wps*wpsi
	end do
c
c last two vectors:
c	
10	wps=wps/(mt-ms(npart-1))
	xnpart_1(1)=x(1,npart-1)
	xnpart_1(2)=x(2,npart-1)
	xnpart_1(3)=0.d0
	call congen_2(mt,pm,xnpart_1,ms(npart-1),
	1    m(npart-1),p(1,npart-1),wpsi)
	wps=wps*wpsi
	p(1,npart)=pm(1)
	p(2,npart)=pm(2)
	p(3,npart)=pm(3)
	p(4,npart)=dsqrt(pm(1)**2+pm(2)**2+pm(3)**2+m(npart)**2)
c
c normalization of the phase space weight:
c
	wps=wps/rsh/(4.d0*pi**2)**(npart-2)/4.d0/pi

	return
	end

	subroutine congen_2(m,pm,x,ms1,m1,p1,wps)
	implicit none
	real*8 m,pm(3),x(3),ms1,m1,p1(4),wps

	real*8 ph,ct,st,mt,p,p1r(4),p2r(4),pm1(3)
	real*8 e,gam,beta,ppar,pparlab
	real*8 alam
	integer i
	real*8 pi
 	parameter(pi=3.1415926535987932)
	ph=2.d0*pi*x(1)
	ct=1.d0-2.d0*x(2)
	st=dsqrt(1.d0-ct**2)
	mt=ms1-m1+x(3)*(m-ms1)
	p=alam(m,m1,mt)/2.d0/m
	wps=p*(m-ms1)

        p1r(1)=p*st*cos(ph)
	p1r(2)=p*st*sin(ph)
	p1r(3)=p*ct
	p1r(4)=dsqrt(p**2+m1**2)
	p2r(1)=-p1r(1)
	p2r(2)=-p1r(2)
	p2r(3)=-p1r(3)
	p2r(4)=dsqrt(p**2+mt**2)
C
C	Boost
C
	p=sqrt(pm(1)**2+pm(2)**2+pm(3)**2)
	if(p.eq.0)then
	   do i=1,3
	      p1(i)=p1r(i)
	      pm(i)=p2r(i)
	   end do
	   p1(4)=p1r(4)
	else
	   pm1(1)=pm(1)/p
	   pm1(2)=pm(2)/p
	   pm1(3)=pm(3)/p
	   e=dsqrt(p**2+m**2)
	   gam=e/m
	   beta=p/e
	   
	   ppar=pm1(1)*p1r(1)+pm1(2)*p1r(2)+pm1(3)*p1r(3)
	   pparlab=gam*(ppar+beta*p1r(4))
	   
	   p1(1)=p1r(1)+(pparlab-ppar)*pm1(1)
	   p1(2)=p1r(2)+(pparlab-ppar)*pm1(2)
	   p1(3)=p1r(3)+(pparlab-ppar)*pm1(3)
	   p1(4)=gam*(p1r(4)+beta*ppar)
	   
	   ppar=pm1(1)*p2r(1)+pm1(2)*p2r(2)+pm1(3)*p2r(3)
	   pparlab=gam*(ppar+beta*p2r(4))
	   pm(1)=p2r(1)+(pparlab-ppar)*pm1(1)
	   pm(2)=p2r(2)+(pparlab-ppar)*pm1(2)
	   pm(3)=p2r(3)+(pparlab-ppar)*pm1(3)
	   
	end if
	
	m=mt

	return
	end

