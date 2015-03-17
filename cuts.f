 	subroutine setup_cuts(iout,iflag)
	implicit none
	integer iout,iflag

	integer switch
	common/process/switch

	real*8 mcut
	common/masscut/mcut

	real*8 etmin1,ymax1
	common/par_cuts_1/etmin1,ymax1
	real*8 etmin2,ymax2
	common/par_cuts_2/etmin2,ymax2
	real*8 etminP,ymaxP
	common/par_cuts_J/etminP,ymaxP
	data etmin1,ymax1/10.d0,100.d0/
	data etmin2,ymax2/10.d0,100.d0/
	data etminP,ymaxP/10.d0,100.d0/

	if(iflag.eq.0)then
	   write(iout,*) ' Use of default cuts'
	   print*,' Use of default cuts'
	else
	  print*,' Cuts for first lepton:'
	  print*,' values of Etmin, and ymax?'
	  read*,etmin1,ymax1
	  print*,' '
	  print*,' Cuts for the second lepton:'
	  print*,' values of Etmin, ymax?'
	  read*,etmin2,ymax2
	  print*,' '
	  print*,' Cuts for the photon:'
	  print*,' values of Etmin, ymax?'
	  read*,etminP,ymaxP
	end if

	write(iout,*)' Cuts for lepton1:'
	write(iout,1)' Etmin',Etmin1,' ymax',Ymax1
	write(iout,*)' Cuts for lepton2:'
	write(iout,1)' Etmin',Etmin2,' ymax',Ymax2
	write(iout,*)' Cuts for the hard photon:'
	write(iout,*)' none applied at the moment'
	write(iout,1)' Etmin',EtminP,' ymax',YmaxP
 1	format(4(A10,f10.2,' '))

	write(iout,*)' Transverse (CC) or Invariant mass cut (NC): '
	if(switch.eq.1)
	1    print*,' Value of MT cut (GeV):'
	if(switch.eq.2)
	1    print*,' Value of Minv cut (GeV):'
	read*,mcut
	if(switch.eq.1)write(iout,*)' MTmin= ',mcut
	if(switch.eq.2)write(iout,*)' Minvmin= ',mcut

	print*,' Setup of the cuts done'
	
	return
	end

	subroutine newcuts(npart,icut)
	implicit none
	integer npart,icut
	integer i,j

	real*8 etmin1,ymax1
	common/par_cuts_1/etmin1,ymax1
	real*8 etmin2,ymax2
	common/par_cuts_2/etmin2,ymax2
	real*8 etminP,ymaxP
        common/par_cuts_J/etminP,ymaxP

	integer switch
	common/process/switch

	integer qqcd,flagqcd
	common/qcdswitch/qqcd
	common/qcdflag/flagqcd

	real*8 mcut
	common/masscut/mcut

	integer test(10)
	common/par_test/test

        real*8 dot_4,pi,dphi
	integer nboost
	real*8 z,st,ct,ecm,pcm,ps(4,1)
	common/fraction/z,nboost

        real*8 p1(4),sp1(4),p2(4),sp2(4),p3(4),sp3(4)
	integer flag,flagr

        character*1 lep
	real*8 scalet
	common/scalet/scalet

	integer ppswitch
	common/collider/ppswitch

        include 'config.inc'

	pi=4d0*datan(1d0)
c
c a priori the event is fine
c
	Icut=0

	do i=1,npart
	  call par_momentum(p(1,i),y(i),pt(i),et(i),phi(i),icut)
	end do

	do i=1,3
	   do j=1,4
	      sp(j,i)=p(j,i)
	   enddo
	enddo

c mFSR: 
	if(nboost.eq.1)then
	   ct=p(3,1)/dsqrt(p(1,1)**2+p(2,1)**2+p(3,1)**2)
	   st=dsqrt(1d0-ct**2)
c	   do i=1,4
c	      p(i,1)=z*p(i,1)
c	   enddo
	   call rotation_1(1d0,st,ct,p(1,1),p(3,1))
c	   do i=1,4
c	      write(6,*)'2',p(i,1)
c	   enddo
c	   ecm=z*dsqrt(sinv(1,2))/2d0
c	   pcm=dsqrt(ecm**2-m(2)**2)
c	   call boost_zz(ecm,pcm,p(1,1))
	   ps(4,1) = z*p(4,1)
	   ps(3,1)=dsqrt(dabs(ps(4,1)**2-m(2)**2))
	   p(4,1) = ps(4,1)
	   p(3,1) = ps(3,1)
c	   do i=1,4
c	      write(6,*)'3',p(i,1)
c	   enddo
	   call rotation_1(1d0,-st,ct,p(1,1),p(3,1))
c	   do i=1,4
c	      write(6,*)'4',p(i,1)
c	   enddo
	endif

c recombination and smearing (TeV4LHC workshop setup)
c
c smearing of the electron/muon, photon and neutrino momentum:
c
c sp(i,3) is the smeared photon momentum (local quantity)
c
	do i=1,4
	   p1(i) = p(i,1)
	   p2(i) = p(i,2) 
	   if (npart.eq.3) p3(i) = p(i,3)
	end do

	if (test(2).eq.1) then
	   lep = 'e'
	elseif (test(2).eq.2) then 
	   lep = 'm'
	endif

	if (test(5).ge.1) then
	   scalet=dsqrt(p1(1)**2+p1(2)**2)
	   if(ppswitch.eq.1) then
	      call d0upgrsmear(p1,sp1,lep)
	      if(switch.eq.1) call metsmear(p2,sp2)
	      if(switch.eq.2) call d0upgrsmear(p2,sp2,lep)
	      if(npart.eq.3.and.flagqcd.eq.0) 
	1	   call d0upgrsmear(p3,sp3,'e')
	   else
	      call atlsmear(p1,sp1,lep)
	      if(switch.eq.1) call atlsmear(p2,sp2,'n')
	      if(switch.eq.2) call atlsmear(p1,sp1,lep)
	      if (npart.eq.3.and.flagqcd.eq.0) 
	1	   call atlsmear(p3,sp3,'e')
	   endif
           do i=1,4
              p(i,1) = sp1(i)
              p(i,2) = sp2(i)
	      if (npart.eq.3.and.flagqcd.eq.0) sp(i,3) = sp3(i)
           end do
        end if

c        do i=1,npart
        do i=1,2
          call par_momentum(p(1,i),y(i),pt(i),et(i),phi(i),icut)
        end do
	if (npart.eq.3) then 
	   if (test(5).ge.1.and.flagqcd.eq.0) then
	      call par_momentum(sp(1,3),y(3),pt(3),et(3),phi(3),icut)
	   else
	      call par_momentum(p(1,3),y(3),pt(3),et(3),phi(3),icut)
	   endif
	endif

	if(test(5).lt.2.or.switch.eq.2.or.flagqcd.ne.0) goto 15

c with recombination:
        flag = 0
        flagr = 0

        if (npart.eq.3) then

           flag = 1
           dphi = dabs(phi(1)-phi(3))
           if (dphi.gt.pi) dphi = 2d0*pi-dphi
           dR(1,3) = dsqrt(dphi**2+(y(1)-y(3))**2)
c
c isolation and recombination cut:
c
c -- electrons
	   if (test(2).eq.1) then
c --Tevatron and LHC
	      if (dR(1,3).lt.0.1d0) then
		 do i=1,4
		    p(i,1) = p(i,1)+sp(i,3)
		    sp(i,3) = 0d0
		    p(i,3) = 0d0
		 end do
		 flagr = 1
		 goto 10
	      end if

	      if ((dR(1,3).ge.0.1d0).and.(dR(1,3).le.0.4d0)) then
		 
		 if (sp(4,3).gt.0.1d0*p(4,1)) then
		    icut = 1
		    return
		 end if
		 
	      endif
c
c -- muons
	   elseif (test(2).eq.2) then
c
c --Tevatron and LHC
c
	      if ((dR(1,3).lt.0.1d0).and.(sp(4,3).gt.2d0)) then
		 icut=1
		 return
	      endif
	      if ((dR(1,3).ge.0.1d0).and.(dR(1,3).le.0.4d0)) then
		 if (sp(4,3).gt.0.1d0*p(4,1)) then
		    icut = 1
		    return
		 end if
	      endif
	      
	   end if
 10        continue

           if (flagr.eq.1) then
              npart = 2
              dR(1,3) = 0d0
              y(3) = -100d0
              pt(3) = -100d0
              et(3) = -100d0
              phi(3) = 0d0
           end if

c           do i=1,npart
           do i=1,2
              call par_momentum(p(1,i),y(i),pt(i),et(i),phi(i),icut)
           end do

        end if

 15	continue

c cuts
	if(et(1).lt.etmin1.or.et(2).lt.etmin2)then
	  icut=1
	  return
	end if
	if(dabs(y(1)).gt.ymax1.or.dabs(y(2)).gt.ymax2)then
	  icut=1
	  return
	end if
	if(switch.eq.1)then
c transverse mass cut:	
	   if(dsqrt(dabs(2d0*pt(1)*pt(2)*(1d0-dcos(phi(2)-phi(1)))))
	1	.lt.mcut)then
	      icut=1
	      return
	   endif
	elseif(switch.eq.2)then
c invariant lepton-pair mass cut:	
	   if(dsqrt(2d0*dot_4(p(1,1),p(1,2))).lt.mcut) then
	      icut=1
	      return
	   endif
	endif

	return
	end

	subroutine par_momentum(p,y,pt,et,phi,icut)
c	
c       input: p(4)
c	output: y pseudo-rapidity
c	        pt transverse momentum
c		et transverse energy 
c		phi azimutal angle
c		icut control flag
c
	implicit none
	real*8 p(4),y,pt,et,phi
	integer icut
	real*8 pm,theta
	real*8 temp,pi

	pi=4d0*datan(1d0)
	pm=dsqrt(p(1)**2+p(2)**2+p(3)**2)
	if(p(1).eq.0.d0.and.p(2).eq.0.d0.or.dabs(1d0-p(3)/pm).lt.1d-8)
	1    then
c       
c       y is infinite
c       
	   pt=0.d0
	   et=0.d0
	   phi=0.d0
	   icut=1
	   return
	end if
c   	y=0.5*log( (p(4)+p(3)) / (p(4)-p(3)) )
c	pm=sqrt(p(1)**2+p(2)**2+p(3)**2)
	theta=dacos(p(3)/pm)
	temp=dtan(theta/2.d0)
	if(temp.le.0.d0)then
	   print*, 'Precision problem in par_momentum'
	   print*,temp,p
	   icut=1
	   return
	end if
	y=-dlog(temp)
	pt=dsqrt(p(1)**2+p(2)**2)
	et=p(4)*dabs(dsin(theta))
	phi=datan2(p(2),p(1))
c
c convert phi to an angle between 0 and 2*pi
c
        if(phi.lt.0.d0)phi=phi+2d0*pi
	
 10	continue
	return
	end




