        subroutine setup_graphs(iout,ndim)
        integer i,iout,ndim
        integer hmax
        parameter(hmax=55)
        integer parts(hmax)
        real*8 mini(hmax),maxi(hmax),stepi(hmax)
        common/hloc/mini,maxi,stepi,parts

	integer switch
	common/process/switch

        do i=1,2
c pt(l)
           if(switch.eq.1)then
              maxi(9+i)=55d0
              mini(9+i)=25d0
              parts(9+i)=120
           else
              maxi(9+i)=65d0
              mini(9+i)=25d0
              parts(9+i)=160
           endif
c eta(l)
           maxi(11+i)=3d0
           mini(11+i)=-3d0
           parts(11+i)=60
c cosine of lepton scattering angle
           maxi(13+i)=1d0
           mini(13+i)=-1d0
           parts(13+i)=100
        enddo
c pt(W,Z)
        maxi(16)=100d0
        mini(16)=0d0
        parts(16)=100
c pt(W,Z)
        maxi(17)=25d0
        mini(17)=0d0
        parts(17)=100
c y(W,Z)
        maxi(18)=3d0
        mini(18)=-3d0
        parts(18)=60
        if(switch.eq.1)then
c MT(l,nu)
           maxi(19)=100d0
           mini(19)=50d0
           parts(19)=100
        else
c Minv(l,l)
           maxi(19)=200d0
           mini(19)=50d0
           parts(19)=150
        endif
c X_M
        maxi(41)=1.2d0
        mini(41)=0.6d0
        parts(41)=100
c X_pt
        maxi(42)=1.2d0
        mini(42)=0.6d0
        parts(42)=100
c photon energy
        maxi(31)=50d0
        mini(31)=0d0
        parts(31)=100
c photon pt
        maxi(32)=100d0
        mini(32)=0d0
        parts(32)=200
c R_lgamma
        maxi(33)=1d0
        mini(33)=0d0
        parts(33)=50
c log10(R_lgamma)
        maxi(34)=0d0
        mini(34)=-6d0
        parts(34)=60
c y photon
        maxi(35)=1d0
        mini(35)=0d0
        parts(35)=50
c log10(y photon)
        maxi(36)=0d0
        mini(36)=-4d0
        parts(36)=40
c
c LO
c
        do i=1,2
c pt(l)
           if(switch.eq.1)then
              maxi(19+i)=55d0
              mini(19+i)=25d0
              parts(19+i)=120
           else
              maxi(19+i)=65d0
              mini(19+i)=25d0
              parts(19+i)=160
           endif
c eta(l)
           maxi(21+i)=3d0
           mini(21+i)=-3d0
           parts(21+i)=60
c cosine of lepton scattering angle
           maxi(23+i)=1d0
           mini(23+i)=-1d0
           parts(23+i)=100
        enddo
c pt(W,Z)
        maxi(26)=100d0
        mini(26)=0d0
        parts(26)=100
c pt(W,Z)
        maxi(27)=25d0
        mini(27)=0d0
        parts(27)=100
c y(W,Z)
        maxi(28)=3d0
        mini(28)=-3d0
        parts(28)=60
        if(switch.eq.1)then
c MT(l,nu)
           maxi(29)=100d0
           mini(29)=50d0
           parts(29)=100
        else
c Minv(l,l)
           maxi(29)=200d0
           mini(29)=50d0
           parts(29)=150
        endif
c X_M
        maxi(51)=1.2d0
        mini(51)=0.6d0
        parts(51)=100
c X_pt
        maxi(52)=1.2d0
        mini(52)=0.6d0
        parts(52)=100

        do i=10,19
           stepi(i)=(maxi(i)-mini(i))/parts(i)
           stepi(i+10)=(maxi(i+10)-mini(i+10))/parts(i+10)
        end do

        do i=31,36
           stepi(i)=(maxi(i)-mini(i))/parts(i)
        end do

        do i=41,42
           stepi(i)=(maxi(i)-mini(i))/parts(i)
           stepi(i+10)=(maxi(i+10)-mini(i+10))/parts(i+10)
        end do

        return
        end

        subroutine graphs_b(dsig,xv,xk,wgt)
        implicit none
        real*8 resx,resy,dot_4
        integer hmax
        parameter(hmax=55)
        integer parts(hmax)
        real*8 mini(hmax),maxi(hmax),stepi(hmax),wi(hmax)
        common/hloc/mini,maxi,stepi,parts

        integer i,ndim,ihist,it,ncall2,icut
        common/par_integration/ndim,ihist,it,ncall2
        real*8 dsig,xv(ndim),xk(ndim),wgt
        real*8 pv(4),yv,ptv,etv,phiv
        real*8 ecm,pcm,beta

        include 'config.inc'
        include 'common.inc'

	integer switch
	common/process/switch

        do i=20,29
           wi(i)=wgt*dsig/(it*stepi(i))
        enddo
        do i=51,52
           wi(i)=wgt*dsig/(it*stepi(i))
        enddo
c pt(l) and eta(l)
        call dfill(20,pt(1),0.d0,wi(20))
        call dfill(21,pt(2),0.d0,wi(21))
        call dfill(22,y(1),0.d0,wi(22))
        call dfill(23,y(2),0.d0,wi(23))
c tests:
c boost to CMS frame
c
c        ecm=(x1+x2)*rs/2.d0
c        pcm=(x1-x2)*rs/2.d0
c        do i=1,2
c           call invboost_zz(ecm,pcm,b(1,i))
c        end do	  
c        do i=1,2
c           call invboost_zz(ecm,pcm,p(1,i))
c        end do
c
c cos theta(l) in LAB
c        cgl(1)=(b(1,1)*p(1,1)+b(2,1)*p(2,1)+b(3,1)*p(3,1))/
c     $       dsqrt(p(1,1)**2+p(2,1)**2+p(3,1)**2)/
c     $       dsqrt(b(1,1)**2+b(2,1)**2+b(3,1)**2)
c        cgl(2)=(b(1,1)*p(1,2)+b(2,1)*p(2,2)+b(3,1)*p(3,2))/
c     $       dsqrt(p(1,2)**2+p(2,2)**2+p(3,2)**2)/
c     $       dsqrt(b(1,1)**2+b(2,1)**2+b(3,1)**2)
c cos theta(l) in parton CMS:
c        beta=(x1-x2)/(x1+x2)
c        cgl(1)=(cgl(1)-beta)/(1d0-beta*cgl(1))
c        cgl(2)=(cgl(2)-beta)/(1d0-beta*cgl(2))
c        cgl(1)=1d0-4d0/(x1*x2*rs**2)*b(4,1)*p(4,1)*(1d0-cgl(1))
c        cgl(2)=1d0-4d0/(x1*x2*rs**2)*b(4,1)*p(4,2)*(1d0-cgl(2))
c
        call dfill(24,cgl(1),0.d0,wi(24))
        call dfill(25,cgl(2),0.d0,wi(25))
c pt(W,Z)
	do i=1,4
	   pv(i)=p(i,1)+p(i,2)
	end do

	call par_momentum(pv,yv,ptv,etv,phiv,icut)

        call dfill(26,ptv,0.d0,wi(26))
        call dfill(27,ptv,0.d0,wi(27))
c
c       pseudo-rapidity of the V is infinite
c       replace y by the rapidity
c
        yv=0.5d0*dlog( (pv(4)+pv(3)) / (pv(4)-pv(3)) )
        call dfill(28,yv,0.d0,wi(28))
c MT or Minv
        if(switch.eq.1)then
           resx=dsqrt(dabs(2d0*pt(1)*pt(2)*(1d0-dcos(phi(2)-phi(1)))))
        else
           resx=dsqrt(dabs(2d0*dot_4(p(1,1),p(1,2))))
        endif
        call dfill(29,resx,0.d0,wi(29))
c X_M
        resx=dsqrt(dabs(2d0*pt(1)*pt(2)*(1d0-dcos(phi(2)-phi(1)))))/mv
        call dfill(51,resx,0.d0,wi(51))
c X_ptl
        resx=pt(1)/mv
        call dfill(52,resx,0.d0,wi(52))

        return
        end

        subroutine graphs(dsig,xv,xk,wgt)
        implicit none
        real*8 resx,resy,dot_4
        integer hmax
        parameter(hmax=55)
        integer parts(hmax)
        real*8 mini(hmax),maxi(hmax),stepi(hmax),wi(hmax)
        common/hloc/mini,maxi,stepi,parts

        integer i,ndim,ihist,it,ncall2,icut
        common/par_integration/ndim,ihist,it,ncall2
        real*8 dsig,xv(ndim),xk(ndim),wgt,dphi
        real*8 pv(4),yv,ptv,etv,phiv

        include 'config.inc'
        include 'common.inc'

	integer switch
	common/process/switch

        do i=10,19
           wi(i)=wgt*dsig/(it*stepi(i))
        enddo
        do i=31,36
           wi(i)=wgt*dsig/(it*stepi(i))
        enddo
        do i=41,42
           wi(i)=wgt*dsig/(it*stepi(i))
        enddo
        wi(37)=wgt*dsig/(it*stepi(36)*stepi(34))

c pt(l) and eta(l)
        call dfill(10,pt(1),0.d0,wi(10))
        call dfill(11,pt(2),0.d0,wi(11))
        call dfill(12,y(1),0.d0,wi(12))
        call dfill(13,y(2),0.d0,wi(13))
        call dfill(14,cgl(1),0.d0,wi(14))
        call dfill(15,cgl(2),0.d0,wi(15))
c pt(W,Z)
	do i=1,4
	   pv(i)=p(i,1)+p(i,2)
	end do

	call par_momentum(pv,yv,ptv,etv,phiv,icut)

        call dfill(16,ptv,0.d0,wi(16))
        call dfill(17,ptv,0.d0,wi(17))
c
c       pseudo-rapidity of the V is infinite
c       replace y by the rapidity
c
        yv=0.5d0*dlog( (pv(4)+pv(3)) / (pv(4)-pv(3)) )
        call dfill(18,yv,0.d0,wi(18))
c MT or Minv
        if(switch.eq.1)then
           resx=dsqrt(dabs(2d0*pt(1)*pt(2)*(1d0-dcos(phi(2)-phi(1)))))
        else
           resx=dsqrt(dabs(2d0*dot_4(p(1,1),p(1,2))))
        endif
        call dfill(19,resx,0.d0,wi(19))
c X_M
        resx=dsqrt(dabs(2d0*pt(1)*pt(2)*(1d0-dcos(phi(2)-phi(1)))))/mv
        call dfill(41,resx,0.d0,wi(41))
c X_ptl
        resx=pt(1)/mv
        call dfill(42,resx,0.d0,wi(42))

c photon distributions:

        call dfill(31,p(4,3),0d0,wi(31))
        call dfill(32,pt(3),0d0,wi(32))

        resx=(p(4,3)/(p(4,3)+p(4,1)))**(1d0/3d0)
        call dfill(35,resx,0d0,wi(35))

        dphi = dabs(phi(1)-phi(3))
        if (dphi.gt.pi) dphi = 2d0*pi-dphi
        resx = dsqrt(dsqrt(dphi**2+(y(1)-y(3))**2))
        call dfill(33,resx,0d0,wi(33))
c
c double differential distribution
c
        resx=dlog10(p(4,3)/(p(4,3)+p(4,1)))
        call dfill(36,resx,0d0,wi(36))
        resy=dlog10(dsqrt(dphi**2+(y(1)-y(3))**2))
        call dfill(34,resy,0d0,wi(34))
        call dfill2(37,36,34,resx,resy,0d0,wi(37))

        return
        end

        subroutine dfill(i,x,y,wgt)
        implicit none
        integer i
        real*8 x,y,wgt
        integer hmax,bmax
        parameter(hmax=55)
        parameter(bmax=250)
        integer ibin,parts(hmax),eventsh(hmax),eventsb(hmax,bmax)
        real*8 mini(hmax),maxi(hmax)
        real*8 averageh(hmax,bmax),sigmah(hmax,bmax),stepi(hmax)
        common/histo/eventsb,averageh,sigmah,eventsh
        common/hloc/mini,maxi,stepi,parts

        if (x.ge.mini(i).and.x.le.maxi(i)) then
           ibin=int(dabs((x-mini(i)))/stepi(i))+1
           eventsh(i)=eventsh(i)+1
           eventsb(i,ibin)=eventsb(i,ibin)+1
           averageh(i,ibin)=wgt+averageh(i,ibin)
           sigmah(i,ibin)=wgt**2+sigmah(i,ibin)
        end if
        return
        end

        subroutine dfill2(k,i,j,x,z,y,wgt)
        implicit none
        integer i,j,k
        real*8 x,y,z,wgt
        integer hmax,bmax
        parameter(hmax=55)
        parameter(bmax=250)
        integer ibinx,ibinz,parts(hmax),eventsh2(hmax),
     $       eventsb2(hmax,bmax,bmax)
        real*8 mini(hmax),maxi(hmax)
        real*8 averageh2(hmax,bmax,bmax),
     $       sigmah2(hmax,bmax,bmax),stepi(hmax)
        common/histo2/eventsb2,averageh2,sigmah2,eventsh2
        common/hloc/mini,maxi,stepi,parts
        if (x.ge.mini(i).and.x.le.maxi(i)) then
           ibinx=int(dabs((x-mini(i)))/stepi(i))+1
           if (z.ge.mini(j).and.z.le.maxi(j)) then
              ibinz=int(dabs((z-mini(j)))/stepi(j))+1
              eventsh2(k)=eventsh2(k)+1
              eventsb2(k,ibinx,ibinz)=eventsb2(k,ibinx,ibinz)+1
              averageh2(k,ibinx,ibinz)=wgt+averageh2(k,ibinx,ibinz)
              sigmah2(k,ibinx,ibinz)=wgt**2+sigmah2(k,ibinx,ibinz)
           endif
        end if
        return
        end

