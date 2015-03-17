c---------------------------------------------------------------------
c     checked version of "integrallib.f"
c---------------------------------------------------------------------
c====================================================================
c     checked with hollik's routines ; status 16.2.93 ; 
c====================================================================
c--------------------------------------------------------------------
      subroutine fint(s,m1,m2,f)
c--------------------------------------------------------------------
c     the f-function as specified in : boehm et.al., 
c     fortschr.phys.34 ; and in : hollik, fortschr.phys. 38 
c--------------------------------------------------------------------
      implicit real*8(a-z)
      parameter (eps=1d-6)
      complex*16 f,sc,ieps
      common/cc/ieps
      sc=s+ieps
      pi = 4d0*datan(1d0)
      x1 = m1**2
      x2 = m2**2
      p = s-(dabs(m1)+dabs(m2))**2
      m = s-(dabs(m1)-dabs(m2))**2
      imagf = 0d0
      rest = 0d0
      if (abs(s).lt.eps) then
         realf = 0d0
         goto 999
      endif
c
      if (x1.lt.eps) then
         if (x2.lt.eps) then
            realf = 0d0
         else
            if (s.gt.x2+eps) then
               realf = 1d0+(1d0-x2/s)*cdlog(1d0/(sc/x2-1d0))
               imagf = pi*(1d0-x2/s)
            else if (s.lt.x2-eps) then
               realf = 1d0+(1d0-x2/s)*dlog(1d0/(1d0-s/x2))
            else
               realf = 1d0
            endif 
         endif
         goto 999
      else if (x2.lt.eps) then
         if (x1.lt.eps) then
            realf = 0d0
         else
            if (s.gt.x1+eps) then
               realf = 1d0+(1d0-x1/s)*cdlog(1d0/(sc/x1-1d0))
               imagf = pi*(1d0-x1/s)
            else if (s.lt.x1-eps) then
               realf = 1d0+(1d0-x1/s)*dlog(1d0/(1d0-s/x1))
            else
               realf = 1d0
            endif 
         endif
         goto 999
      endif
c
      if ((abs(s).lt.x1/1d3).or.(abs(s).lt.x2/1d3)) then
         if (dabs(dabs(m1)-dabs(m2)).lt.eps) then
            realf = s/6d0/x1
         else
            realf = s/(x1-x2)**2*((x1+x2)/2d0
     $           -x1*x2/(x1-x2)*dlog(x1/x2))
         endif
         goto 999
      endif
c
      if (dabs(dabs(m1)-dabs(m2)).lt.eps) then
         rest = 2d0
      else
         rest = 1d0-((x1-x2)/s-(x1+x2)/(x1-x2))*dlog(x1/x2)/2d0
      endif
c
      if (m.lt.0d0) then
         realf = dsqrt(p*m)*dlog((dsqrt(-p)+dsqrt(-m))**2
     $        /(4d0*dsqrt(x1*x2)))/s
      else if (p.lt.0d0) then
         realf = -2d0*dsqrt(-p*m)*datan(dsqrt(-m/p))/s
      else
         realf = -dsqrt(p*m)*dlog((dsqrt(p)+dsqrt(m))**2
     $        /(4d0*dsqrt(x1*x2)))/s
         imagf = dsqrt(p*m)/s*pi
      endif
 999  continue
      f = dcmplx(rest+realf,imagf)
      end
c---------------------------------------------------------------------
c---------------------------------------------------------------------
      subroutine afunc(m,a0)
      implicit real*8(a-z)
      parameter (eps=1d-4)
      common/renorm/mue,mue2,mue4
      muedr=mue
      if (m**2.lt.eps) then 
         a0 = 0d0
      else
         a0 = m**2*(1d0-dlog(m**2/muedr**2))
      endif
      end
c---------------------------------------------------------------------
      subroutine bfunc(s,m1,m2,b0,b1,b20)
c---------------------------------------------------------------------
c     b0 and b1 are the b0"quer" and b1"quer"-quantities defined 
c     in : hollik, fortschr.phys. 38 (but with log-terms incl.)
c     b20 is the g_{mu}.{nu} coefficient of the 
c     tensor-2-point-integral (log-terms incl, too)
c---------------------------------------------------------------------
      implicit real*8(a-z)
      parameter (eps=1d-4)
      complex*16 b0,b1,b20,ff,sc,ieps
      common/cc/ieps
      common/renorm/mue,mue2,mue4
      muedr=mue
      sc=s+ieps
      pi = 4d0*datan(1d0)
      x1 = m1**2
      x2 = m2**2
      if (dabs(s).lt.eps) then
         call bnull(m1,m2,realb0)
         call beins(m1,m2,realb1)
         call bzwnull(m1,m2,realb20)
         b0 = dcmplx(realb0,0d0)
         b1 = dcmplx(realb1,0d0)
         b20 = dcmplx(realb20,0d0)
         goto 999
      endif
      if ((x1.lt.eps).and.(x2.lt.eps)) goto 47 
      call fint(s,m1,m2,ff)
 47   continue
      if (x1.lt.eps) then
         if (x2.lt.eps) then
            b0 = 2d0-cdlog(sc/muedr**2)+dcmplx(0d0,1d0)*pi
            b1 = -b0/2d0
         else
            b0 = 1d0+ff-dlog(x2/muedr**2)
            b1 = -(1d0+(1d0-x2/s)*ff-dlog(x2/muedr**2))/2d0
         endif
         goto 666
      endif
      if (x2.lt.eps) then
         if (x1.lt.eps) then
            b0 = 2d0-cdlog(sc/muedr**2)+dcmplx(0d0,1d0)*pi
            b1 = -b0/2d0
         else
            b0 = 1d0+ff-dlog(x1/muedr**2)
            b1 = -(1d0+(1d0+x1/s)*ff-dlog(x1/muedr**2))/2d0
         endif
         goto 666
      endif
      if (dabs(x1-x2).lt.1d3*eps) then
         b0 = ff
         b1 = -ff/2d0
      else
         b0 = 1d0-(x1+x2)/(x1-x2)*dlog(x1/x2)/2d0+ff
         b1 = -1d0/2d0+x1/(x1-x2)*dlog(x1/x2)/2d0-(s+x1-x2)/2d0/s*ff
      endif
      b0 = b0 - dlog(x1*x2/muedr**4)/2d0
      b1 = b1 + dlog(x2/muedr**2)/2d0 
c     
 666  continue
      call afunc(m2,a0)
      b20 = a0/6d0+x1/3d0*b0
     $     +(s+x1-x2)/6d0*b1+(x1+x2-s/3d0)/6d0
 999  continue
      end
c---------------------------------------------------------------------
      subroutine bderiv(x,m1,m2,db0,db1,db20)
c---------------------------------------------------------------------
c     real parts of d(b0(s,m1,m2))/ds , d(b1(s,m1,m2))/ds
c     and d(b20(s,m1,m2))/ds
c     (s -> x in the code)
c---------------------------------------------------------------------
      implicit real*8(a-z)
      complex*16 cf,b0,b1,b20
      parameter (eps=1d-1)
*
      pi = 4d0*datan(1d0)
      xm1 = m1**2
      xm2 = m2**2
*
      if (x.lt.eps) then
         if (dabs(xm1-xm2).lt.eps) then
            db0 = 1d0/6d0/xm2
            db1 = -1d0/12d0/xm2
         else if (xm1.lt.eps) then
            db0 = 1d0/2d0/xm2
            db1 = -1d0/6d0/xm2
         else if (xm2.lt.eps) then
            db0 = 1d0/2d0/xm1
            db1 = -1d0/6d0/xm1
         else
            db0 = (xm1**2-xm2**2-2d0*xm1*xm2*dlog(xm1/xm2))
     $           /2d0/(xm1-xm2)**3
            fss0 = (xm1**3-xm2**3+9d0*xm1*xm2*(xm1-xm2)
     $           -6d0*xm1*xm2*dlog(xm1/xm2)*(xm1+xm2))
     $           /3d0/(xm1-xm2)**5
            db1 = -db0/2d0-(xm1-xm2)/4d0*fss0
         endif
         goto 11
      endif
 13   continue
      if ((x.lt.xm1/1d3).or.(x.lt.xm2/1d3)) then
         if (dabs(xm1-xm2).lt.eps) then
            db0 = 1d0/6d0/xm1
         else
            db0 = ((xm1+xm2)/2d0-xm1*xm2/(xm1-xm2)*dlog(xm1/xm2))
     $           /(xm1-xm2)**2
         endif
         db1 = -db0/2d0
         goto 11
      endif
c
      if ((xm1.lt.eps).and.(xm2.lt.eps)) then
         db0 = -1d0/x
         db1 = -db0/2d0
         goto 11
      endif
c
      if (xm2.lt.eps) then
         if (x.gt.xm1+eps) then
            deriv = -(1d0+xm1/x*dlog(x/xm1-1d0))/x
         else if (x.lt.xm1-eps) then
            deriv = -(1d0+xm1/x*dlog(1d0-x/xm1))/x
         else
            deriv = 0d0
         endif 
         goto 10
      endif
c
      if (xm1.lt.eps) then
         if (x.gt.xm2+eps) then
            deriv = -(1d0+xm2/x*dlog(x/xm2-1d0))/x
         else if (x.lt.xm2-eps) then
            deriv = -(1d0+xm2/x*dlog(1d0-x/xm2))/x
         else
            deriv = 0d0
         endif 
         goto 10
      endif
c
      sm = xm1+xm2
      dm = xm2-xm1
      sm12 = (dabs(m1)+dabs(m2))**2
      dm12 = (dabs(m1)-dabs(m2))**2
      lm = dlog(xm2/xm1)/2d0
*
      if (x.lt.dm12) then
         s = dsqrt(sm12-x)
         d = dsqrt(dm12-x)
         fact = dlog((d+s)**2/(4d0*dabs(m1)*dabs(m2)))
         deriv = (dm*lm/x-((d**2+s**2+2d0*d**2*s**2/x)
     $        /(2d0*s*d))*fact-1d0)/x
      else if (x.lt.sm12) then
         if (dabs(xm1-xm2).lt.eps) then
            if (x.lt.eps) then
               deriv = 1d0/6d0/xm1
            else
               sx = 4d0*xm1/x
               deriv = (sx/dsqrt(sx-1d0)
     $              *datan(1d0/dsqrt(sx-1d0))-1d0)/x
            endif
         else
*
c Achtung: fuer x=sm12 geht deriv -> infty! 
*
            s = dsqrt(sm12/x-1d0)
            d = dsqrt(1d0-dm12/x)
            fact = datan(d/s)
            deriv = (dm*lm/x+((d**2-s**2+2d0*d**2*s**2)
     $           /(s*d))*fact-1d0)/x
         endif
      else
         s = dsqrt(x-sm12)
         d = dsqrt(x-dm12)
         fact = dlog((d+s)**2/(4d0*dabs(m1)*dabs(m2)))
         deriv = (dm*lm/x-((d**2+s**2-2d0*d**2*s**2/x)
     $        /(2d0*s*d))*fact-1d0)/x
      endif
*
 10   continue
c--------------------------------------------------------------------
      db0 = deriv
      call fint(x,m1,m2,cf)
      f = dreal(cf)
      db1 = -(x+xm1-xm2)/2d0/x*deriv+(xm1-xm2)/2d0/x**2*f
c     
 11   continue
c
      call bfunc(x,m1,m2,b0,b1,b20)
      db20 = xm1/3d0*db0+dreal(b1)/6d0+(x+xm1-xm2)/6d0*db1
     $     -1d0/18d0
c
 99   continue
      end
c---------------------------------------------------------------------
      subroutine cfuncnew(m,s,m1,m2,m3,c0,c1p,c1m,c20,c2p,c2m)
c---------------------------------------------------------------------
c     contains the scalar coefficient functions of the 
c     tensor three point integral
c---------------------------------------------------------------------
      implicit real*8(a-z)
      complex*16 c0,c1p,c1m,c20,c2p,c2m,C0_
     $     ,b0x32,b0x31,b0s12,b1x32,b1x31,adummy,bdummy
      x = m**2
      x1 = m1**2
      x2 = m2**2
      x3 = m3**2
c
      call bfunc(s,m1,m2,b0s12,adummy,bdummy)
      call bfunc(x,m3,m2,b0x32,b1x32,adummy)
      call bfunc(x,m3,m1,b0x31,b1x31,adummy)
c new C0 function added on Aug.21 2004
c      call c0int(m,s,m1,m2,m3,c0)
      c0=C0_(m**2,s,m**2,m3,m1,m2,0)
c
      c1m = ((b0x32+b0x31)/2d0-b0s12
     $     -(x+x3-(x1+x2)/2d0)*c0)/(4d0*x-s)
      c1p = ((b0x32-b0x31)/2d0+(x1-x2)/2d0*c0)/s
      c20 = b0s12/4d0+x3/2d0*c0-(x1-x2)/4d0*c1p
     $     -(x1+x2-2d0*(x+x3))/4d0*c1m+1d0/4d0
      c2m = (-c20+(b1x32+b1x31)/4d0+b0s12/2d0
     $     -(x+x3-(x1+x2)/2d0)*c1m)/(4d0*x-s)
      c2p = (-c20-(b1x32+b1x31)/4d0+(x1-x2)/2d0*c1p)/s
      end
c---------------------------------------------------------------------
      subroutine c0int(mf,s,m1,m2,m3,c0)
c**************************************************************
c                                                             *
c  the scalar vertex integral with equal external masses mf   *
c                                                             *
c-------------------------------------------------------------------
c     s = momentum transfer; m1,m2,m3  are the internal masses
c
      implicit real*8 (a-y)
      complex*16 z1,z2,z11,z12,z21,z22,cl1,cl2,cl3,cspen,spence,
     &     int,c0
      xmf=mf*mf
      if (mf.lt.1d-1) then
         mfstrich = 1d-1
         xmf = mfstrich**2
      endif
c     xm's : are fermion and boson masses squared
      xm1=m1*m1
      xm2=m2*m2
      xm3=m3*m3
c     the t'hooft-veltman parameters
      a=1.d0
      b=xmf/s
      c=-1.d0
      d=xm1-xm2-s
      e=xm3-xm1-xmf+s
      f=xm2/s
      d=d/s
      e=e/s
c     discriminante for alpha-equation
      disc=c*c-4.d0*a*b
      if (disc .lt. 0.d0) goto 500
      al=(-c-dsqrt(disc))/2.d0/b
      nenner=c+2.d0*al*b
c..the first integral.............................................
      y0=-(d+e*al+2.d0*a+c*al)/nenner
      y01=y0-1.d0
      d1=(c+e)**2-4.d0*b*(a+d+f)
      x1=-(c+e)/2.d0/b
      if (d1.gt.0.d0) goto 10
c.......complex zeroes of logarithms
      sq1=dsqrt(-d1)
      x2=sq1/2.d0/b
      z1=dcmplx(x1,x2)
      z2=dcmplx(x1,-x2)
      z11=y0/(y0-z1)
      z12=y01/(y0-z1)
      z21=y0/(y0-z2)
      z22=y01/(y0-z2)
      cl1=spence(z11)-spence(z12)+spence(z21)-spence(z22)
      goto 15
10    continue
c........real zeroes
      sq1=dsqrt(d1)
      x2=sq1/2.d0/b
      y1=x1+x2
      y2=x1-x2
      sig1= y0/dabs(y0)
      sig2= y01/dabs(y01)
      y11=y0/(y0-y1)
      y12=y01/(y0-y1)
      y21=y0/(y0-y2)
      y22=y01/(y0-y2)
      cl1=cspen(y11,sig1)-cspen(y12,sig2)+cspen(y21,-sig1)
     &   -cspen(y22,-sig2)
15    continue
c..the second integral............................................
      y0=-(d+e*al)/nenner/(1.d0-al)
      y01=y0-1.d0
      d2=(e+d)**2-4.d0*f*(a+b+c)
      x1=-(e+d)/2.d0/(a+b+c)
      if(d2.gt.0.d0) goto 20
c.......complex zeroes of logarithms
      sq2=dsqrt(-d2)
      x2=sq2/2.d0/(a+b+c)
      z1=dcmplx(x1,x2)
      z2=dcmplx(x1,-x2)
      z11=y0/(y0-z1)
      z12=y01/(y0-z1)
      z21=y0/(y0-z2)
      z22=y01/(y0-z2)
      cl2=spence(z11)-spence(z12)+spence(z21)-spence(z22)
      goto 25
20    continue
c........real zeroes
      x2=dsqrt(d2)/2.d0/(a+b+c)
      y1=x1+x2
      y2=x1-x2
      y11=y0/(y0-y1)
      y12=y01/(y0-y1)
      y21=y0/(y0-y2)
      y22=y01/(y0-y2)
      sig1= y0/dabs(y0)
      sig2= y01/dabs(y01)
      cl2=cspen(y11,sig1)-cspen(y12,sig2)+cspen(y21,-sig1)
     &   -cspen(y22,-sig2)
25    continue
c..the third integral............................................
      y0=(d+e*al)/nenner/al
      y01=y0-1.d0
      d3=d*d-4.d0*a*f
      x1=-d/2.d0/a
      if (d3.gt.0.d0) goto 30
c........complex zeroes of logarithms
      sq3=dsqrt(-d3)
      x2=sq3/2.d0/a
      z1=dcmplx(x1,x2)
      z2=dcmplx(x1,-x2)
      z11=y0/(y0-z1)
      z12=y01/(y0-z1)
      z21=y0/(y0-z2)
      z22=y01/(y0-z2)
      cl3=spence(z11)-spence(z12)+spence(z21)-spence(z22)
      goto 35
30    continue
c........real zeroes
      x2=dsqrt(d3)/2.d0/a
      y1=x1+x2
      y2=x1-x2
 31   format(1h ,3e12.4)
      y11=y0 /(y0-y1)
      y12=y01/(y0-y1)
      y21=y0/(y0-y2)
      y22=y01/(y0-y2)
      sig1= y0/dabs(y0)
      sig2= y01/dabs(y01)
      cl3=cspen(y11,sig1)-cspen(y12,sig2)+cspen(y21,-sig1)
     &   -cspen(y22,-sig2)
35    continue
c..summation of the 3 integrals ....................................
      int = -cl1+cl2-cl3
      c0 = int/nenner/s
      goto 501
500   continue
c..error message for complex alpha................................
      write(6,21)
21    format(1h ,'  i cannot handle a complex alpha')
501   return
      end
c--------------------------------------------------------------------
c     the x"null" subroutines calculate the b0-, b1-, c0-,
c     d0- and d20-integrals at zero external momenta  
c--------------------------------------------------------------------
      subroutine bnull(a,b,b0)
      implicit real*8(a-z)
      parameter (eps=1d-8)      
      common/renorm/mue,mue2,mue4
      muedr=mue
      xa = a**2
      xb = b**2
      x = xa+xb
      if (x.lt.eps) then
c        write(6,*)'all args zero is not allowed for b0-function !'
         b0 = 0d0
         goto 2
      endif
      if (xa*xb.eq.0d0) then
         zw = 1d0-dlog(x/muedr**2)
      else
         zw = -dlog(xa*xb/muedr**4)/2d0
         if (dabs(dabs(a)-dabs(b)).gt.eps) then
            zw = zw+1d0-(xa+xb)/(xa-xb)*dlog(xa/xb)/2d0
         endif
      endif
      b0 = zw
 2    continue
      end
c-----------------------------------------------------------------     
      subroutine beins(a,b,b1)
      implicit real*8(a-z)
      parameter (eps=1d-8)
      common/renorm/mue,mue2,mue4
      muedr=mue
      xa = a**2
      xb = b**2
      x = xa+xb
      if (x.lt.eps) then
c        write(6,*)'all args zero is not allowed for b1-function !'
         b1 = 0d0
         goto 2
      endif      
      if (xa.eq.0d0) then
         zw = -(1d0/2d0-dlog(xb/muedr**2))/2d0
      else if (xb.eq.0d0) then
         zw = -(3d0/2d0-dlog(xa/muedr**2))/2d0
      else
         zw = dlog(xb/muedr**2)/2d0
         if (dabs(dabs(a)-dabs(b)).gt.eps) then
            zw = zw -(xa+xb)/(xa-xb)/4d0-1d0/2d0
     $           +xa/(xa-xb)*dlog(xa/xb)/2d0
     $           +xa*xb/(xa-xb)**2*dlog(xa/xb)/2d0
         endif
      endif
      b1 = zw
 2    continue
      end
c-----------------------------------------------------------------     
      subroutine bzwnull(a,b,b20)
      implicit real*8(a-z)
      common/renorm/mue,mue2,mue4
      muedr=mue
      xa = a**2
      xb = b**2
      call bnull(a,b,b0)
      call beins(a,b,b1)
      zw = xa*b0/3d0+(xa-xb)*b1/6d0+(xa+xb)/6d0
      if (xb.ne.0d0) then 
         zw = zw+xb*(1d0-dlog(xb/muedr**2))/6d0
      endif
      b20 = zw
      end
c-----------------------------------------------------------------
      subroutine cnull(a,b,c,c0)
      implicit real*8(a-z)
      parameter (eps=1d-8)
      xa = a**2
      xb = b**2
      xc = c**2
      if (dabs(dabs(a)-dabs(b)).lt.eps) then
         if (dabs(dabs(a)-dabs(c)).lt.eps) then
            zw = -1d0/2d0/xa
         else 
            zw = (-1d0+xc/(xa-xc)*dlog(xa/xc))/(xa-xc)
         endif
      else 
         call bnull(a,c,bac)
         call bnull(b,c,bbc)
         zw = (bac-bbc)/(xa-xb)
      endif
      c0 = zw
      end
c-----------------------------------------------------------------     
      subroutine dnull(a,b,c,d,d0)
      implicit real*8(a-z)
      integer i,j
      real*8 m(4)
      parameter (eps=1d-8)
      m(1) = a
      m(2) = b
      m(3) = c
      m(4) = d
      do 30 i = 1,4
         do 40 j = i+1,4
            if (dabs(dabs(m(i))-dabs(m(j))).lt.eps) then
               m(i) = m(i)+eps
               m(j) = m(j)-eps
            endif
 40      continue
 30   continue
      call cnull(a,b,c,c0abc)
      call cnull(a,b,d,c0abd)
      zw = (c0abc-c0abd)/(c**2-d**2)
      d0 = zw
      end
c-----------------------------------------------------------------
      subroutine dzwnull(a,b,c,d,d20)
      implicit real*8(a-z)
      call cnull(b,c,d,c0bcd)
      call dnull(a,b,c,d,d0abcd)
      zw = c0bcd+a**2*d0abcd
      d20 = zw
      end
c-----------------------------------------------------------------     
************************************************************************
        function C0_(q01,q12,q20,m0,m1,m2,ext)                            
************************************************************************
*	scalar 3-point function
*-----------------------------------------------------------------------
*	General result from A.Denner, Fortschr. Phys. 41 (1993) 307
* d0=q^2-m0^2, d1=(q+p1)^2-m1^2, d2=(q+p2)^2-m2^2
* q01=p1^2, q20=p2^2, q12=(p1-p2)^2
*-----------------------------------------------------------------------
************************************************************************
        implicit real*8 (a-z)                                         
	integer i,j,ext
	complex*16 c0_,cspens,etass,ieps,alpha,alp(0:2),sc,xs
	complex*16 y0(0:2),y(0:2,-1:1),x(0:2,-1:1)
	real*8 thp(0:2),thy(0:2)
	common/photon/lambda

	kappa2(a,b,c) = a**2+b**2+c**2-2d0*a*b-2d0*a*c-2d0*b*c
        lambda2=lambda*lambda
        pi   = 4d0*datan(1d0)                                               
	eps  = 1d-17
	ieps = dcmplx(0d0,eps)
	m02 = m0**2
	m12 = m1**2
	m22 = m2**2
	p01 = q01
	p12 = q12
	p20 = q20
C*** Check for IR divergence
c    Permutate until m0=0
10	continue
	if ((m02*m12*m22.eq.0d0).and.(m02.ne.0d0)) then
	  m   = m02
	  m02 = m12
	  m12 = m22
	  m22 = m
	  p   = p01
	  p01 = p12
	  p12 = p20
	  p20 = p
	else
	  goto 20
	endif
	goto 10  
20	continue
	if ((p01.eq.m12).and.(p20.eq.m22).and.(m02.eq.0d0)) goto 500

C****** Regular C0 function
c    Permutate until p01=0
11	continue
	if ((p01*p12*p20.eq.0d0).and.(p01.ne.0d0)) then
	  m   = m02
	  m02 = m12
	  m12 = m22
	  m22 = m
	  p   = p01
	  p01 = p12
	  p12 = p20
	  p20 = p
	else
	  goto 21
	endif
	goto 11  
21	continue
	if (p01.eq.0d0) goto 600

C*** Regular C0 function with p01,p12,p20 =/= 0
	alpha  = sqrt( abs(kappa2(p01,p12,p20)) )
	alp(0) = sqrt( kappa2(p12,m12,m22)*(1d0+ieps*sign(1d0,p12)) )
	alp(1) = sqrt( kappa2(p20,m22,m02)*(1d0+ieps*sign(1d0,p20)) )
	alp(2) = sqrt( kappa2(p01,m02,m12)*(1d0+ieps*sign(1d0,p01)) )

	do 99 i=0,2
	  if (alp(i).eq.dcmplx(0d0,0d0)) alp(i) = ieps*abs(alpha)
99	continue

	y0(0)  = ( p12*(p12-p01-p20+2d0*m02-m12-m22)
     &	  - (p20-p01)*(m12-m22)+alpha*(p12-m12+m22) )/2d0/alpha/p12
	y0(1)  = ( p20*(p20-p12-p01+2d0*m12-m22-m02)
     &	  - (p01-p12)*(m22-m02)+alpha*(p20-m22+m02) )/2d0/alpha/p20
	y0(2)  = ( p01*(p01-p20-p12+2d0*m22-m02-m12)
     &	  - (p12-p20)*(m02-m12)+alpha*(p01-m02+m12) )/2d0/alpha/p01

	do 100 j=-1,1,2
	  x(0,j) = (p12-m12+m22+j*alp(0))/2d0/p12
	  x(1,j) = (p20-m22+m02+j*alp(1))/2d0/p20
	  x(2,j) = (p01-m02+m12+j*alp(2))/2d0/p01
	do 100 i=0,2
	  y(i,j) = y0(i)-x(i,j)
100	continue

	do 200 i=0,2
	  thp(i) = 0d0
	  thy(i) = 0d0
	  if (dimag(y(i,+1)*y(i,-1)).le.0d0) thy(i) = 1d0
200	continue
	if (p12.le.0d0) thp(0) = 1d0
	if (p20.le.0d0) thp(1) = 1d0
	if (p01.le.0d0) thp(2) = 1d0

	c0_ = 0d0
	do 300 i=0,2
	do 400 j=-1,1,2
	  c0_ = c0_ + cspens((y0(i)-1d0)/y(i,j)) - cspens(y0(i)/y(i,j))
     &	         + etass(1d0-x(i,j),1d0/y(i,j))*log((y0(i)-1d0)/y(i,j))
     &	            - etass(   -x(i,j),1d0/y(i,j))*log(y0(i)/y(i,j))
400	continue
	  c0_ = c0_ - log((y0(i)-1d0)/y0(i))*(
     &		        etass(-x(i,+1),-x(i,-1))-etass(y(i,+1),y(i,-1))
     &		      - dcmplx(0d0,2d0*pi)*thp(i)*thy(i) )
300	continue
	c0_ = c0_/alpha

	return
600	continue

	if (m02*m12.eq.0d0) goto 700
C*** Regular C0 function with p01=0 and m0,m1=/=0
	alpha   = 1d0+(p12-p20)/(m02-m12-ieps*(m02+m12))
	alp(0)  = sqrt( kappa2(p20,m02,m22)+ieps*p20*(p20-m02-m22) )
	alp(1)  = sqrt( kappa2(p12,m12,m22)+ieps*p12*(p12-m12-m22) )
	x(0,+1) = (p20-m02-m22+alp(0))/2d0/m02
	x(0,-1) = (p20-m02-m22-alp(0))/2d0/m02
	x(1,+1) = (p12-m12-m22+alp(1))/2d0/m12
	x(1,-1) = (p12-m12-m22-alp(1))/2d0/m12
	c0_ = 0d0
	do 610 i=-1,1,2
	do 620 j=0,1
	c0_ = c0_ + (1d0-2d0*j)*( 
     &		cspens(1d0+x(j,i)/alpha) - cspens(1d0+x(j,i)) 
     &		+ etass(1d0/alpha,-x(j,i))*log(1d0+x(j,i)/alpha) )
620	continue
610	continue
	c0_ = c0_ + log(alpha)*log(m02/m12)
	c0_ = c0_/(p20-p12) 

	return
700	continue

C*** Regular C0 function with p01=0 and m0=0
	if (m12.eq.0d0) then
	  m12 = m02
	  m02 = 0d0
	  p   = p12
	  p12 = p20
	  p20 = p
	endif
	alpha   = 1d0+(p12-p20)/(-m12-ieps*m12)
	alp(1)  = sqrt( kappa2(p12,m12,m22)+ieps*p12*(p12-m12-m22) )
	x(0,-1) = m22/(p20-m22+ieps*m12)
	x(1,+1) = (p12-m12-m22+alp(1))/2d0/m12
	x(1,-1) = (p12-m12-m22-alp(1))/2d0/m12
	c0_ = 0d0
	do 710 i=-1,1,2
	c0_ = c0_ - cspens(1d0+x(1,i)/alpha) + cspens(1d0+x(1,i)) 
     &		  - etass(1d0/alpha,-x(1,i))*log(1d0+x(1,i)/alpha) 
710	continue
	c0_ = c0_ + cspens(1d0+x(0,-1)/alpha) - cspens(1d0+x(0,-1)) 
     &		  + etass(1d0/alpha,-x(0,-1))*log(1d0+x(0,-1)/alpha) 
	c0_ = c0_ + log(alpha)*log((m22-ieps*m12-p20)/m12)
     &		  - log(alpha)**2/2d0
	c0_ = c0_/(p20-p12) 

	return
500	continue

C*** IR-divergent C0 function
        write(6,*)'warning: C0 is IR divergent'
c        stop
	sc  = p12+abs(p12)*ieps
	mm1 = sqrt(m12)
	mm4 = sqrt(m22)
	xs = -4d0*mm1*mm4/(sc-(mm1-mm4)**2) /
     &	     ( sqrt(1d0-4d0*mm1*mm4/( sc-(mm1-mm4)**2))+1d0 )**2 
	c0_ = xs/mm1/mm4/(1d0-xs**2)*(
     &	    log(xs)*( -log(xs)/2d0+2d0*log(1d0-xs**2)
     &		    -log(lambda2/mm1/mm4) )
     &	  - pi**2/6d0+cspens(xs**2)+log(mm1/mm4)**2/2d0
     &	  + cspens(1d0-xs*mm1/mm4) + cspens(1d0-xs*mm4/mm1) )
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        FUNCTION CSPENS(Z)                                              
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SPENCE-FUNKTION KOMPLEX, FREI NACH HOLLIK                     C
C---------------------------------------------------------------------C
C       20.07.83    LAST CHANGED 10.05.89        ANSGAR DENNER        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        COMPLEX*16 CSPENS,W,SUM,Z,U                                     
        REAL*8 RZ,AZ,A1                                                
        REAL*8 B(9)/                                                   
     1   0.1666666666666666666666666667D0,                             
     2  -0.0333333333333333333333333333D0,                             
     3   0.0238095238095238095238095238D0,                             
     4  -0.0333333333333333333333333333D0,                             
     5   0.0757575757575757575757575758D0,                             
     6  -0.2531135531135531135531135531D0,                             
     7   1.1666666666666666666666666667D0,                             
     8  -7.09215686274509804D0         ,                               
     9  54.97117794486215539D0         /                               
C     BEACHTE:                 B(N)=B2N                                
C     B(1)=1./6.                                                       
C     B(2)=-1./30.                                                     
C     B(3)=1./42.                                                      
C     B(4)=-1./30.                                                     
C     B(5)=5./66.                                                      
C     B(6)=-691./2730.                                                 
C     B(7)=7./6.                                                       
C     B(8)=-3617./510.                                                 
C     B(9)=43867./798.                                                 
C     B(10)=-174611./330.                                              
C     B(11)=854513./138.                                               
C     PI=3.1415926535897932384                                         
C     PI*PI/6.=1.6449..., PI*PI/3=3.28986...                           
C                                                                      
      Z =Z*DCMPLX(1D0)                                                 
      RZ=DREAL(Z)                                                      
      AZ=CDABS(Z)                                                      
      A1=CDABS(1D0-Z)                                                  
C     IF((SNGL(RZ) .EQ. 0.0) .AND. (SNGL(DIMAG(Z)) .EQ. 0.0)) THEN     
C ---> CHANGED  10.5.89                                                
      IF(AZ .LT. 1D-20) THEN                                           
        CSPENS=-CDLOG(1D0-Z)                                            
        RETURN                                                         
      END IF                                                           
c      IF((SNGL(RZ) .EQ. 1.0) .AND. (SNGL(DIMAG(Z)) .EQ. 0.0)) THEN     
c ---> changed 5.7.94
       IF( (ABS(RZ-1D0).LT.1D-18) .AND. (ABS(DIMAG(Z)).LT.1D-18) ) THEN     
        CSPENS=1.64493406684822643D0                                    
        RETURN                                                         
      END IF                                                           
      IF(RZ.GT.5D-1) GOTO 20                                           
      IF(AZ.GT.1D0) GOTO 10                                            
      W=-CDLOG(1D0-Z)                                                  
      SUM=W-0.25D0*W*W                                                 
      U=W                                                              
      IF(CDABS(U).LT.1D-10) GOTO 2                                     
      DO 1 K=1,9                                                       
      U=U*W*W/DFLOAT(2*K*(2*K+1))                                      
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 2                            
      SUM=SUM+U*B(K)                                                   
 1    CONTINUE                                                         
 2    CSPENS=SUM                                                        
      RETURN                                                           
10    W=-CDLOG(1D0-1D0/Z)                                              
      SUM=W-0.25D0*W*W                                                 
      U=W                                                              
      IF(CDABS(U).LT.1D-10) GOTO 12                                    
                                                                       
      DO 11 K=1,9                                                      
      U=U*W*W/DFLOAT(2*K*(2*K+1))                                      
      IF(CDABS(B(K)*U/SUM).LT.1D-20) GOTO 12                           
      SUM=SUM+U*B(K)                                                   
11    CONTINUE                                                         
12    CSPENS=-SUM-1.64493406684822643D0-.5D0*CDLOG(-Z)**2               
      RETURN                                                           
20    IF(A1.GT.1D0) GOTO 30                                            
      W=-CDLOG(Z)                                                      
      SUM=W-0.25D0*W*W                                                 
      U=W                                                              
      IF(CDABS(U).LT.1D-10) GOTO 22                                    
      DO 21 K=1,9                                                      
      U=U*W*W/DFLOAT(2*K*(2*K+1))                                      
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 22                           
      SUM=SUM+U*B(K)                                                   
21    CONTINUE                                                         
22    CSPENS=-SUM+1.64493406684822643D0-CDLOG(Z)*CDLOG(1D0-Z)           
      RETURN                                                           
30    W=CDLOG(1D0-1D0/Z)                                               
      SUM=W-0.25D0*W*W                                                 
      U=W                                                              
      IF(CDABS(U).LT.1D-10) GOTO 32                                    
      DO 31 K=1,9                                                      
      U=U*W*W/DFLOAT(2*K*(2*K+1))                                      
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 32                           
      SUM=SUM+U*B(K)                                                   
31    CONTINUE                                                         
32    CSPENS=SUM+3.28986813369645287D0                                  
     *               +.5D0*CDLOG(Z-1D0)**2-CDLOG(Z)*CDLOG(1D0-Z)       
50    CONTINUE                                                         
      END                                                            
***********************************************************************
        FUNCTION ETASS(C1,C2)                                            
***********************************************************************
*       COMPLEX ETA-FUNKTION                                           
*---------------------------------------------------------------------*
*       8.06.90    ANSGAR DENNER                                       
***********************************************************************
        IMPLICIT   LOGICAL(A-Z)                                        
        COMPLEX*16 ETASS,C1,C2                                           
        REAL*8     PI,IM1,IM2,IM12                                     
                                                                       
        PI     = 4D0*DATAN(1D0)                                        
        IM1    = DIMAG(C1)                                             
        IM2    = DIMAG(C2)                                             
        IM12   = DIMAG(C1*C2)                                          

	if (((IM1.eq.0d0).and.(DREAL(C1).lt.0d0)).or.
     &	    ((IM2.eq.0d0).and.(DREAL(C2).lt.0d0)).or.
     &	    ((IM12.eq.0d0).and.(DREAL(C1*C2).lt.0d0))) then
	  write(*,*) 'etass function on cut !!!'
	  write(*,*) 'C1    = ',C1
	  write(*,*) 'C2    = ',C2
	  write(*,*) 'C1*C2 = ',C1*C2
	  stop
	endif
                                                                       
        IF(IM1.LT.0D0.AND.IM2.LT.0D0.AND.IM12.GT.0D0) THEN             
            ETASS = DCMPLX(0D0,2D0*PI)                                   
        ELSE IF (IM1.GT.0D0.AND.IM2.GT.0D0.AND.IM12.LT.0D0) THEN       
            ETASS = DCMPLX(0D0,-2D0*PI)                                  
        ELSE                                                           
            ETASS = DCMPLX(0D0)                                          
        END IF                                                         
        END                                                            

***********************************************************************
        FUNCTION ETAS(Y,R,RS)                                            
***********************************************************************
*       MODIFIED ETA-FUNKTION                                           
*---------------------------------------------------------------------*
***********************************************************************
        IMPLICIT   LOGICAL(A-Z)                                        
        COMPLEX*16 ETA,ETAS,Y,R,RS
        REAL*8     PI,IMY,IMRS
                                                                       
        PI     = 4D0*DATAN(1D0)                                        

	IF( DIMAG(R).NE.0D0 ) THEN
	    ETAS = ETA(Y,R)
	ELSE	    
	    IF( DREAL(R).GT.0D0 ) THEN
		ETAS = DCMPLX(0D0,0D0)
	    ELSE
	 	IMY  = DIMAG(Y)
		IMRS = DIMAG(RS)
		ETAS = 2D0*DCMPLX(0D0,PI)*(
     *			(1D0+SIGN(1D0,-IMY))*(1D0+SIGN(1D0,-IMRS))-
     *			(1D0+SIGN(1D0, IMY))*(1D0+SIGN(1D0, IMRS))
     *					  )/4D0
	    ENDIF
	ENDIF
        END                                                            

***********************************************************************
        FUNCTION SQE(A,B,C)                                            
***********************************************************************
*       SOLUTION OF QUADRATIC EQUATION				      *
*---------------------------------------------------------------------*
***********************************************************************
        IMPLICIT REAL*8 (A-Z)                                        
        COMPLEX*16 A,B,C,SQE,X1,X2

	X1=(-B+SQRT(B**2-4D0*A*C))/2D0/A
	X2=(-B-SQRT(B**2-4D0*A*C))/2D0/A

	IF (ABS(X1).GT.ABS(X2)) THEN
	   SQE=X1
	ELSE
	   SQE=X2
	ENDIF

        END                                                            

************************************************************************
        FUNCTION D0_(P1,P2,P3,P4,P12,P23,M1,M2,M3,M4,ext) 
************************************************************************
*  SCALAR 4-POINT FUNCTION                                             *
*  P1,P2,P3,P4 = SQUARED EXTERNAL MOMENTA			       *
*  P12 = (p1+p2)**2,  P23 = (p2+p3)**2				       *
*----------------------------------------------------------------------*
*	General result from                                            *
*        A.Denner, U.Nierste and R.Scharf, Nucl. Phys. B367 (1991) 637 *
*	IR-divergent case from                                         *
*        W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349         *
*----------------------------------------------------------------------*
************************************************************************
        IMPLICIT REAL*8 (A-Z) 
	REAL*8 M(4),P(4,4),K(4,4)                                     
	COMPLEX*16 A1,A2,A3,A4,SWAP
	COMPLEX*16 SS(4), XX(2), X(2,4),RS(4,4)
	COMPLEX*16 S0(4),XX0(2),X0(2,4), R(4,4),G(2)
        COMPLEX*16 D0_,D0_ext,CSPENS,ETA,SQE,ETAS
	COMPLEX*16 AA,BB,CC,DD,IEPS,H,HH,L1,L2,L3,L4
	COMPLEX*16 SC,TC,XS,X2,X3,q2c,q3c,Y
	INTEGER I,J,ext

	common/photon/lambda

	if (ext.ne.0) then
	  D0_ = D0_ext(P1,P2,P3,P4,P12,P23,M1,M2,M3,M4,ext)
	  return
	endif

        lambda2=lambda*lambda
        PI = 4D0*DATAN(1D0)                                               

        MM1=M1    
        MM2=M2    
        MM3=M3    
        MM4=M4    
        M12=M1*M1 
        M22=M2*M2 
        M32=M3*M3 
        M42=M4*M4 
        Q1=P1 
        Q2=P2   
        Q3=P3   
	Q4=P4
        Q12=P12   
        Q23=P23
C	IS AT LEAST ONE MASS ZERO ???
	IF (MM1*MM2*MM3*MM4.NE.0D0) GOTO 130

C--->	****** IR-divergent CASE ******
	IF ( ((Q1.EQ.M12).AND.(Q2.EQ.M32).AND.(M22.EQ.0D0)).OR.
     *	     ((Q2.EQ.M22).AND.(Q3.EQ.M42).AND.(M32.EQ.0D0)).OR.
     *	     ((Q3.EQ.M32).AND.(Q4.EQ.M12).AND.(M42.EQ.0D0)).OR.
     *	     ((Q4.EQ.M42).AND.(Q1.EQ.M22).AND.(M12.EQ.0D0)) ) goto 50

C--->	****** REGULAR CASE with at least one mass zero ******

C	PERMUTATE UNTIL MM3=0D0
	GOTO 20
10	CONTINUE
	MM0=MM1
	MM1=MM2
	MM2=MM3
	MM3=MM4
	MM4=MM0
	M02=M12
	M12=M22
	M22=M32
	M32=M42
	M42=M02
	Q00=Q12
	Q12=Q23
	Q23=Q00
	Q0=Q1
	Q1=Q2
	Q2=Q3
	Q3=Q4
	Q4=Q0
20	IF (MM3.NE.0D0) GOTO 10
C	ONLY MM3 IS ZERO
	IF (MM1*MM2*MM4.NE.0D0) GOTO 30
C	ONLY MM3 AND MM4 ARE ZERO ==> 3->2, 4->3...   
	IF ((MM1*MM2.NE.0D0).AND.(MM4.EQ.0D0)) GOTO 10
C	ONLY MM2 AND MM3 ARE ZERO
	IF ((MM1*MM4.NE.0D0).AND.(MM2.EQ.0D0)) GOTO 40
C	ONLY MM1 AND MM3 ARE ZERO ==> m1 <-> m2
	IF ((MM2*MM4.NE.0D0).AND.(MM1.EQ.0D0)) then
	  mm0 = mm1
	  mm1 = mm2
	  mm2 = mm0
	  q0  = q2
	  q2  = q12
	  q12 = q0
	  q0  = q4
	  q4  = q23
	  q23 = q0
	  goto 40
	endif
C 	check whether all masses are zero
	if ((mm1.eq.0d0).and.(mm2.eq.0d0).and.(mm4.eq.0d0)) then
	  WRITE(*,*) 'D0 case with all mi=0 not implemented !'
	  stop
	endif
C 	permutate until mm1 is non-zero
	if (mm1.eq.0d0) goto 10
	goto 70

C	****** NO MASS EQUAL TO ZERO ******
130	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)

	IF( ABS((MM1**2+MM3**2-Q12)/MM1/MM3).LT.2D0 ) THEN
C	R13 WOULD BE NOT REAL. -> PERMUTATION! -> R(2,4) IS NOT REAL.
	   M(1)=MM2
	   M(2)=MM3
	   M(3)=MM4
	   M(4)=MM1
	   P(1,2)=Q2
	   P(1,3)=Q23
	   P(1,4)=Q1
	   P(2,3)=Q3
	   P(2,4)=Q12
	   P(3,4)=Q4
	ELSE
C	R(1,3) IS REAL.
	   M(1)=MM1
	   M(2)=MM2
	   M(3)=MM3
	   M(4)=MM4
	   P(1,2)=Q1
	   P(1,3)=Q12
	   P(1,4)=Q4
	   P(2,3)=Q2
	   P(2,4)=Q23
	   P(3,4)=Q3
	ENDIF

	DO 11 J=2,4
	DO 11 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
11	CONTINUE

	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	S0(1)=R(1,2)
	S0(2)=R(2,3)
	S0(3)=R(3,4)
	S0(4)=R(1,4)
	AA=K(3,4)/R(2,4)+R(1,3)*K(1,2)-K(1,4)*R(1,3)/R(2,4)-K(2,3)
	BB=(R(2,4)-1D0/R(2,4))*(R(1,3)-1D0/R(1,3))
     *		+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)/R(1,3)+R(2,4)*K(3,4)-K(1,4)*R(2,4)/R(1,3)-K(2,3)
	DD=K(2,3)-R(1,3)*K(1,2)-R(2,4)*K(3,4)+R(1,3)*R(2,4)*K(1,4)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	XX0(1)=SQE(AA,BB,CC)
	XX0(2)=CC/AA/XX0(1)

c	XX(1)=XX0(1)-IEPS*DD/(XX0(1)-XX0(2))
c	XX(2)=XX0(2)+IEPS*DD/(XX0(1)-XX0(2))

c	IF (ABS(DREAL(XX0(1)-XX(2))).LT.ABS(DREAL(XX0(1)-XX(1)))) THEN
	IF (ABS(XX0(1)-XX(2)).LT.ABS(XX0(1)-XX(1))) THEN
	  SWAP  =XX0(1)
	  XX0(1)=XX0(2)
	  XX0(2)=SWAP
	ENDIF

	DO 12 I=1,2
	G(I)  =SIGN( 1D0,DREAL(AA*(XX(I)-XX(3-I))) )
	 X(I,1)= XX(I)/R(2,4)
	X0(I,1)=XX0(I)/R(2,4)
	 X(I,2)= XX(I)/R(2,4)*R(1,3)
	X0(I,2)=XX0(I)/R(2,4)*R(1,3)
	 X(I,3)= XX(I)*R(1,3)
	X0(I,3)=XX0(I)*R(1,3)
	 X(I,4)= XX(I)
	X0(I,4)=XX0(I)
12	CONTINUE

	D0_ = DCMPLX(0D0,0D0)
	DO 13 I=1,2
	DO 13 J=1,4
	A1 = 1D0+X0(I,J)*S0(J) + ABS(1D0+X0(I,J)*S0(J))*IEPS*
     *				  SIGN(1D0,DIMAG(X(I,J)*SS(J)))
	A2 = 1D0+X0(I,J)/S0(J) + ABS(1D0+X0(I,J)/S0(J))*IEPS*
     *				  SIGN(1D0,DIMAG(X(I,J)/SS(J)))
	D0_ = D0_ + (-1D0)**(I+J)*(
     *		CSPENS(A1)+ETA(-X(I,J),SS(J))*LOG(A1)
     *	       +CSPENS(A2)+ETA(-X(I,J),1D0/SS(J))*LOG(A2)     )
13	CONTINUE

	IF( DIMAG(R(1,3)).EQ.0D0 ) THEN
	DO 14 I=1,2
	   A1 = (K(1,3)-2D0*R(1,3))/XX0(I)
     *		      -R(1,3)*K(1,4)+K(3,4)
     	   A2 = ((K(2,4)-2D0*R(2,4))*R(1,3)*XX0(I)
     *		      -R(2,4)*K(3,4)+K(2,3))/DD
	   A3 = (K(1,3)-2D0*R(1,3))*R(2,4)/XX0(I)
     *		      -R(1,3)*K(1,2)+K(2,3)
	   A4 = ((K(2,4)-2D0*R(2,4))*XX0(I)
     *		      -R(2,4)*K(1,4)+K(1,2))/DD
	   L1 = LOG( A1-ABS(A1)*IEPS )
     	   L2 = LOG( A2+ABS(A2)*IEPS*G(I)*SIGN(1D0,DREAL(R(1,3))
     *				        	  *DIMAG(RS(2,4))) ) 
	   L3 = LOG( A3-ABS(A3)*IEPS )
	   L4 = LOG( A4+ABS(A4)*IEPS*G(I)*SIGN(1D0,DIMAG(RS(2,4))) ) 

	   D0_ = D0_ + (3D0-2D0*I)*(
     *		 ETAS(-XX(I),R(1,3),RS(1,3))
     *		   *( LOG(R(1,3)*XX(I)) + L1 + L2 )
     *		+ETAS(-XX(I),1D0/R(2,4),1D0/RS(2,4))
     *		   *( LOG(XX(I)/R(2,4)) + L3 + L4 )
     *		-( ETAS(-XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4))
     *		  +ETA(RS(1,3),1D0/RS(2,4)) )
     *		   *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 + L2 )
     *	  	+ETA(RS(1,3),1D0/RS(2,4))
     *		   *ETAS(-XX(I),-R(1,3)/R(2,4),-RS(1,3)/RS(2,4))   )
14	CONTINUE
	ELSE
	DO 15 I=1,2
	   L1 = LOG( R(2,4)/XX0(I)+XX0(I)/R(2,4)+K(1,2)
     *		     -XX0(I)/R(2,4)*EPS*BB*G(I) )
	   L2 = LOG( R(1,3)*XX0(I)+1D0/XX0(I)/R(1,3)+K(3,4)
     *		     -XX0(I)*R(1,3)*EPS*BB*G(I) )
	   L3 = LOG( R(1,3)/R(2,4)*XX0(I)+R(2,4)/XX0(I)/R(1,3)+K(2,3)
     *		     -XX0(I)*R(1,3)/R(2,4)*EPS*BB*G(I) )
	   D0_ = D0_ + (3D0-2D0*I)*(
     *		+ETA(-XX(I),1D0/R(2,4))
     *		   *( LOG(XX(I)/R(2,4)) + L1 )
     *		+ETA(-XX(I),R(1,3))
     *		   *( LOG(R(1,3)*XX(I)) + L2 )
     *		-( ETA(-XX(I),R(1,3)/R(2,4))
     *		  +ETA(R(1,3),1D0/R(2,4)) )
     *		   *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 )
     *	  	+ETA(R(1,3),1D0/R(2,4))
     *		   *ETA(-XX(I),-R(1,3)/R(2,4))
     *		   *(1D0-G(I)*SIGN(1D0,DREAL(BB)))	    )
15	CONTINUE
	ENDIF

	D0_ = D0_/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	RETURN

C	****** ONLY MM3 IS ZERO ******
30	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)
	M(1)=MM1
	M(2)=MM2
	M(3)=10D0
	M(4)=MM4
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	DO 1 J=2,4
	DO 1 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
1	CONTINUE
	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	AA=K(3,4)/R(2,4)-K(2,3)
	BB=K(1,3)*(1D0/R(2,4)-R(2,4))+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(1,3)*K(1,4)*R(2,4)+R(2,4)*K(3,4)-K(2,3)
	DD=K(2,3)-R(2,4)*K(3,4)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	DO 2 I=1,2
	X(I,1)=XX(I)/R(2,4)
	X(I,2)=XX(I)/R(2,4)*R(1,3)
	X(I,3)=XX(I)*R(1,3)
	X(I,4)=XX(I)
2	CONTINUE
	D0_ = DCMPLX(0D0,0D0)
	DO 3 I=1,2
	D0_ = D0_ + (2D0*I-3D0)*(
     *		CSPENS(1D0+SS(4)*X(I,4))
     *	       -CSPENS(1D0+SS(1)*X(I,1))
     *	       +CSPENS(1D0+X(I,4)/SS(4))
     *	       -CSPENS(1D0+X(I,1)/SS(1))
     *	       +ETA(-X(I,4),SS(4))*LOG(1D0+SS(4)*X(I,4))
     *	       -ETA(-X(I,1),SS(1))*LOG(1D0+SS(1)*X(I,1))
     *	       +ETA(-X(I,4),1D0/SS(4))*LOG(1D0+X(I,4)/SS(4))
     *	       -ETA(-X(I,1),1D0/SS(1))*LOG(1D0+X(I,1)/SS(1))
     *	       -CSPENS(1D0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       +CSPENS(1D0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))
     *	       -ETA(-X(I,4),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       +ETA(-X(I,1),(K(2,3)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))   )
	IF (DIMAG(R(2,4)).NE.0D0) THEN
	   H=ETA(-1D0/XX(I),R(2,4))
	ELSE
	   H=DCMPLX(0D0,0D0)
	   IF (DREAL(R(2,4)).LT.0D0) THEN
	      HH=-1D0/XX(I)
	      IM1=DIMAG(HH)
	      IM2=DIMAG(RS(2,4))
	      IF ((IM1.GT.0D0).AND.(IM2.GT.0D0)) THEN
	         H=-DCMPLX(0D0,2D0*PI)
	      ENDIF
	      IF ((IM1.LT.0D0).AND.(IM2.LT.0D0)) THEN
	         H=+DCMPLX(0D0,2D0*PI)
	      ENDIF
	   ENDIF
	ENDIF
	D0_ = D0_ + (2D0*I-3D0)*
     *	          H*( LOG( (K(1,2)-R(2,4)*K(1,4)
     *			  +XX(I)*(1D0/R(2,4)-R(2,4)))/DD )
     *		     +LOG(K(1,3)-IEPS) )
3	CONTINUE
	D0_ = D0_/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	RETURN

C	****** ONLY MM2 AND MM3 ARE ZERO ******
40	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)

	M(1)=MM1
	M(2)=10D0
	M(3)=10D0
	M(4)=MM4
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	DO 4 J=2,4
	DO 4 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	IF (I.EQ.2) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.2) K(I,J)=K(I,J)-M(J)/M(I)
	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
4	CONTINUE
	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	AA=K(2,4)*K(3,4)-K(2,3)
	BB=K(1,3)*K(2,4)+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(2,3)
	DD=K(2,3)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	DO 5 I=1,2
	X(I,1)=XX(I)/R(2,4)
	X(I,2)=XX(I)/R(2,4)*R(1,3)
	X(I,3)=XX(I)*R(1,3)
	X(I,4)=XX(I)
5	CONTINUE
	D0_ = DCMPLX(0D0,0D0)
	DO 6 I=1,2
	D0_ = D0_ + (2D0*I-3D0)*(
     *		CSPENS(1D0+SS(4)*X(I,4))
     *	       +CSPENS(1D0+X(I,4)/SS(4))
     *	       +ETA(-X(I,4),SS(4))*LOG(1D0+SS(4)*X(I,4))
     *	       +ETA(-X(I,4),1D0/SS(4))*LOG(1D0+X(I,4)/SS(4))
     *	       -CSPENS(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -CSPENS(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	       -ETA(-XX(I),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -ETA(-XX(I),(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS)) 
     *	       +LOG(-XX(I))*( LOG(K(1,2)-IEPS)
     *			     +LOG(K(1,3)-IEPS)-LOG(K(2,3)-IEPS) ) )
6	CONTINUE
	D0_ = D0_/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))

	return

C	****** ONLY MM1 IS NON-ZERO ******
70	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)

	M(1)=MM1
	M(2)=10D0
	M(3)=10D0
	M(4)=10D0
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	k(1,2) = (m(1)**2-p(1,2))/m(1)/m(2)
	k(1,3) = (m(1)**2-p(1,3))/m(1)/m(3)
	k(1,4) = (m(1)**2-p(1,4))/m(1)/m(4)
	k(2,3) = -p(2,3)/m(2)/m(3)
	k(2,4) = -p(2,4)/m(2)/m(4)
	k(3,4) = -p(3,4)/m(3)/m(4)
	DO 74 J=2,4
	DO 74 I=1,J-1
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
74	CONTINUE
	AA=K(2,4)*K(3,4)
	BB=K(1,3)*K(2,4)+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(2,3)
	DD=K(2,3)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	D0_ = DCMPLX(0D0,0D0)
	DO 76 I=1,2
	D0_ = D0_ + (2D0*I-3D0)*(
     *	        CSPENS(1D0+(K(1,4)-IEPS)*xx(i))
     *	       +eta(-xx(i),K(1,4)-IEPS)*log(1D0+(K(1,4)-IEPS)*xx(i))
     *	       -CSPENS(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -CSPENS(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	       -ETA(-XX(I),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -ETA(-XX(I),(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS)) 
     *	       +LOG(-XX(I))*( LOG(K(1,2)-IEPS)
     *			     +LOG(K(1,3)-IEPS)-LOG(K(2,3)-IEPS) ) )
76	CONTINUE
	D0_ = D0_/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	return

C	****** general IR-divergent D0-function ******
50	CONTINUE
        write(6,*)'warning: D0 is IR divergent'
c        stop
C	PERMUTATE UNTIL MM1 IS THE PHOTON
	GOTO 52
51	CONTINUE
	MM0=MM1
	MM1=MM2
	MM2=MM3
	MM3=MM4
	MM4=MM0
	M02=M12
	M12=M22
	M22=M32
	M32=M42
	M42=M02
	Q00=Q12
	Q12=Q23
	Q23=Q00
	Q0=Q1
	Q1=Q2
	Q2=Q3
	Q3=Q4
	Q4=Q0
52	IF ((mm1.ne.0d0).or.(q1.ne.m22).or.(q4.ne.m42)) GOTO 51

	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)
	sc  = q23+abs(q23)*ieps
	tc  = q12+abs(q12)*ieps
	q2c = q2+abs(q2)*ieps
	q3c = q3+abs(q3)*ieps
	xs = -4d0*mm2*mm4/(sc-(mm2-mm4)**2) /
     &	     ( sqrt(1d0-4d0*mm2*mm4/( sc-(mm2-mm4)**2))+1d0 )**2 

	if (mm3.eq.0d0) goto 60

C	*** general case ***
	if (q2.ne.(mm2-mm3)**2) then
	  x2 = -4d0*mm2*mm3/(q2c-(mm2-mm3)**2) /
     &	     ( sqrt(1d0-4d0*mm2*mm3/( q2c-(mm2-mm3)**2))+1d0 )**2 
	else
	  x2 = 1d0
	endif
	if (q3.ne.(mm4-mm3)**2) then
	  x3 = -4d0*mm4*mm3/(q3c-(mm4-mm3)**2) /
     &	     ( sqrt(1d0-4d0*mm4*mm3/( q3c-(mm4-mm3)**2))+1d0 )**2 
	else
	  x3 = 1d0
	endif

	d0_ = xs/mm2/mm4/(q12-m32)/(1d0-xs**2)*(
     &	 2d0*cdlog(xs)*(cdlog(1d0-xs**2)-cdlog(lambda*mm3/(m32-tc)))
     &	+pi**2/2d0+cspens(xs**2)+cdlog(x2)**2+cdlog(x3)**2
     &	-cspens(xs*x2*x3)-(cdlog(xs)+cdlog(x2)+cdlog(x3))
     $       *cdlog(1d0-xs*x2*x3)
     &	-cspens(xs*x2/x3)-(cdlog(xs)+cdlog(x2)-cdlog(x3))
     $       *cdlog(1d0-xs*x2/x3)
     &	-cspens(xs/x2*x3)-(cdlog(xs)-cdlog(x2)+cdlog(x3))
     $       *cdlog(1d0-xs/x2*x3)
     &	-cspens(xs/x2/x3)-(cdlog(xs)-cdlog(x2)-cdlog(x3))
     $       *cdlog(1d0-xs/x2/x3) )
	return

60	continue
C	*** special case: mass mm3 opposite to photon is 0 ***
	if ((q2.eq.m22).or.(q3.eq.m42)) goto 61
	Y = mm2/mm4*(q3c-m42)/(q2c-m22)
	d0_ = xs/mm2/mm4/q12/(1d0-xs**2)*(
     &	 log(xs)*( -log(xs)/2d0+2d0*log(1d0-xs**2)-log(lambda2/mm2/mm4)
     &		   -log((q2-m22)/tc)-log((q3-m42)/tc) )
     &	+pi**2/6d0+cspens(xs**2)+log(y)**2/2d0
     &	-cspens(xs*y)-(log(xs)+log(y))*log(1d0-xs*y)
     &	-cspens(xs/y)-(log(xs)-log(y))*log(1d0-xs/y) )
	return

61	continue
C	*** special case: doubly IR-divergent D0 ***
	if ((q2.eq.m22).and.(q3.eq.m42)) then
	d0_ = -xs/mm2/mm4/q12/(1d0-xs**2)*2d0*log(xs)*log(-lambda2/tc)
	else
	  write(*,*) 'Special case of IR-divergent D0 not implemented!'
          stop
	endif

	END
************************************************************************
        function D0_ext(P1,P2,P3,P4,P12,P23,M1,M2,M3,M4,ext)
************************************************************************
*       scalar 4-point function
*-----------------------------------------------------------------------
*	special cases specified by "ext" 
*	default = general result (ext=0)
*-----------------------------------------------------------------------
************************************************************************
        implicit real*8 (a-z)
        complex*16 D0_,D0_ext,ieps,z1,z2,cspens
        integer ext

	common/photon/lambda
        common/param/pi,el,alpha,alpha0,alphaz,alphas,GF,cw,cw2,sw,sw2,
     &               mw,mw2,gw,mz,mz2,gz,mh,mh2,
     &               mxe,mxe2,ml(3),ml2(3),mqp(3),mqp2(3),mqm(3),mqm2(3)

	ieps = dcmplx(0d0,1d-20)
        lambda2=lambda*lambda
        
	if (ext.eq.401) then
c***    --------------------
c***	D0(me^2,kp^2,me^2,km^2,u,t,me,0,MW,0)  with me -> 0
	  if ((p1.eq.mxe2).and.(p3.eq.mxe2).and.(m1.eq.mxe).and.
     &	      (m2.eq.0d0).and.(m3.eq.MW).and.(m4.eq.0d0)) then
    	    if ((p2.eq.MW2).and.(p4.eq.MW2)) then
	      u = p12
	      t = p23
	      goto 401
	    endif
	  endif
c***	D0(me^2,kp^2,me^2,km^2,u,t,MW,0,me,0)  with me -> 0
	  if ((p1.eq.mxe2).and.(p3.eq.mxe2).and.(m1.eq.MW).and.
     &	      (m2.eq.0d0).and.(m3.eq.mxe).and.(m4.eq.0d0)) then
    	    if ((p2.eq.MW2).and.(p4.eq.MW2)) then
	      u = p12
	      t = p23
	      goto 401
	    endif
	  endif
	  write(*,*) 'Inconsistent call of D0_ext(...,401) !'
	  stop
401	  continue
	  D0_ext = (-log(ml(1)/mw)**2
     &		    +log((mw2-u)/mw/ml(1)-ieps)*log(mw2/lambda2)
     &	            +log((mw2-u)/mw2-ieps)**2
     &	            +2d0*log((mw2-u)/ml(1)/mw-ieps)*log(t/mw2+ieps)
     &		    +cspens(1d0-(u-mw2)/mw2-ieps) )/t/(u-mw2)
	  return

	elseif (ext.eq.402) then
c***    --------------------
c***	D0(kp^2,me^2,me^2,km^2,t,s,0,me,0,me)  with me -> 0
	  if ((p2.eq.mxe2).and.(p3.eq.mxe2).and.(m1.eq.0d0).and.
     &	      (m2.eq.mxe).and.(m3.eq.0d0).and.(m4.eq.mxe)) then
    	    if ((p1.eq.MW2).and.(p4.eq.MW2)) then
	      t = p12
	      s = p23
	      goto 402
	    endif
	  endif
	  write(*,*) 'Inconsistent call of D0_ext(...,402) !'
	  stop
402	  continue
	  D0_ext = ( -log(-s/ml2(1)-ieps)**2/2d0
     &		     +log(-s/ml2(1)-ieps)*( 
     &		        +log(-s/lambda2-ieps)
     &			+2d0*log(t/mw2+ieps) )
     &		    -pi**2/6d0 )/s/t
	  return

	elseif (ext.eq.403) then
c***    --------------------
c***	D0(kp^2,me^2,me^2,km^2,t,s,MW,0,me,0)  with me -> 0
	  if ((p2.eq.mxe2).and.(p3.eq.mxe2).and.(m1.eq.MW).and.
     &	      (m2.eq.0d0).and.(m3.eq.mxe).and.(m4.eq.0d0)) then
    	    if ((p1.eq.MW2).and.(p4.eq.MW2)) then
	      t = p12
	      s = p23
	      goto 403
	    endif
	  endif
	  write(*,*) 'Inconsistent call of D0_ext(...,402) !'
	  stop
403	  continue
	  D0_ext = 2d0*log((mw2-t)/ml(1)/mw-ieps)
     &		   *log(-s/lambda2-ieps)/s/(t-mw2)
	  return

	elseif (ext.eq.404) then
c***    --------------------
c***	D0(kp^2,me^2,me^2,km^2,t,s,MW,0,me,MZ)  with me -> 0
	  if ((p2.eq.mxe2).and.(p3.eq.mxe2).and.(m1.eq.MW).and.
     &	      (m2.eq.0d0).and.(m3.eq.mxe).and.(m4.eq.MZ)) then
    	    if ((p1.eq.MW2).and.(p4.eq.MW2).and.(p12.lt.0d0)) then
	      t = p12
	      s = p23
	      goto 404
	    endif
c***	D0(kp^2,me^2,me^2,km^2,t,s,MW,MZ,me,0)  with me -> 0
	  elseif ((p2.eq.mxe2).and.(p3.eq.mxe2).and.(m1.eq.MW).and.
     &	          (m2.eq.MZ).and.(m3.eq.mxe).and.(m4.eq.0d0)) then
    	    if ((p1.eq.MW2).and.(p4.eq.MW2).and.(p12.lt.0d0)) then
	      t = p12
	      s = p23
	      goto 404
	    endif
	  endif
	  write(*,*) 'Inconsistent call of D0_ext(...,404) !'
	  stop
404	  continue
	  z1     = dcmplx(1d0,+sqrt(4d0*mw2/mz2-1d0))/2d0
	  z2     = dcmplx(1d0,-sqrt(4d0*mw2/mz2-1d0))/2d0
	  D0_ext = (-log((mw2-t)/ml(1)/mw-ieps)**2
     &		    +log((mw2-t)/ml(1)/mw-ieps)
     &		     *( log(lambda2/ml2(1))-2d0*log(1d0-s/mz2-ieps) )
     &		    -cspens(1d0-(mw2-t)/mw2*z1)+pi**2/6d0
     &		    -cspens(1d0-(mw2-t)/mw2*z2)	    )/(t-mw2)/(mz2-s)
	  return

	endif

	D0_ext = D0_(P1,P2,P3,P4,P12,P23,M1,M2,M3,M4,0)

	end
