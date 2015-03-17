      subroutine dmn(q1,q2,q3,q4,q5,q6,m1,m2,m3,m4,d1,d2,d3,ch1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       vierpunktfkt. mit bis zu 4 integrationsimpulsen in zaehler
c       realteil
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       05.05.87 sa
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8(a-y)
      implicit complex*16(z)
      real*8 c1123(0:2),c2123(0:3),c1134(0:2),c2134(0:3),d2(0:6)
      real*8 c1234(0:2),c2234(0:3),c1124(0:2),c2124(0:3),d1(0:3)
      real*8 d3(0:3,0:3),ch1(0:2)
      complex*16 d0gen,d0reg
*
      q7 = q1+q2+q3+q6-q4-q5
      call chutmn(q2,q3,q7,m2,m3,m4,c1234,c2234)
      call chutmn(q1,q2,q4,m1,m2,m3,c1123,c2123)
      call chutmn(q1,q7,q6,m1,m2,m4,c1124,c2124)
      call chutmn(q4,q3,q6,m1,m3,m4,c1134,c2134)
      ch1(0) = c1234(0)
      ch1(1) = -c1234(0)-c1234(1)-c1234(2)
      ch1(2) = c1234(2)
      m12 = m1*m1
      m22 = m2*m2
      m32 = m3*m3
      m42 = m4*m4
      p12 = (q1+q4-q2)/2d0
      p13 = (q1+q6-q7)/2d0
      p23 = (q4+q6-q3)/2d0
      det = q1*q4*q6+2d0*p12*p13*p23-q1*p23*p23-q4*p13*p13-q6*p12*p12
      mat11 = q4*q6-p23*p23
      mat12 = p13*p23-q6*p12
      mat13 = p12*p23-q4*p13
      mat21 = mat12
      mat22 = q1*q6-p13*p13
      mat23 = p12*p13-q1*p23
      mat31 = mat13
      mat32 = mat23
      mat33 = q1*q4-p12*p12
      cf1 = q1+m12-m22
      cf2 = q4+m12-m32
      cf3 = q6+m12-m42
c      if ((m22.lt.q1).and.(m32.lt.q1).and.(m42.lt.q1)) goto 50
c         d1(0) = ggttf(q7,q4,q1,m12,m22)
c      if (m12.ge.4d0*m22) then
c         d1(0) = dreal(d0gen(q1,q2,q3,q4,q5,q6,m1,m2,m3,m4))
c      end if
c      goto 100
c--->   mass singular 3 fermion box
c 50   continue
c      if (dabs(m1-dsqrt(q1)).lt.1d-2) then
c         m1 = m1+1d-2
c         m12 = m1**2
c      end if
c      d1(0) = ggttb(q7,q4,q1,m12,m22)
c      goto 100
 100  d0new = dreal(d0reg(q1,q2,q3,q6,q7,q4,m1,m2,m3,m4))
      d1(0) = d0new
      s11 = (c1134(0)-c1234(0)-cf1*d1(0))/2d0
      s12 = (c1124(0)-c1234(0)-cf2*d1(0))/2d0
      s13 = (c1123(0)-c1234(0)-cf3*d1(0))/2d0
      d1(1) = (mat11*s11+mat12*s12+mat13*s13)/det
      d1(2) = (mat21*s11+mat22*s12+mat23*s13)/det
      d1(3) = (mat31*s11+mat32*s12+mat33*s13)/det
      ccc = c1234(1)+c1234(2)+c1234(0)
      s20 = c1234(0)+m12*d1(0)
      s211 = (ccc-cf1*d1(1))/2d0
      s212 = (c1134(1)-c1234(1)-cf1*d1(2))/2d0
      s213 = (c1134(2)-c1234(2)-cf1*d1(3))/2d0
      s221 = (c1124(1)+ccc-cf2*d1(1))/2d0
      s222 = -(c1234(1)+cf2*d1(2))/2d0
      s223 = (c1124(2)-c1234(2)-cf2*d1(3))/2d0
      s231 = (c1123(1)+ccc-cf3*d1(1))/2d0
      s232 = (c1123(2)-c1234(1)-cf3*d1(2))/2d0
      s233 = -(c1234(2)+cf3*d1(3))/2d0
      d2(0) = s20-s211-s222-s233
      d2(1) = (mat11*(s211-d2(0))+mat12*s221+mat13*s231)/det
      d2(2) = (mat21*s212+mat22*(s222-d2(0))+mat23*s232)/det
      d2(3) = (mat31*s213+mat32*s223+mat33*(s233-d2(0)))/det
c---  >   d2(1,2)
      d2(4) = (mat11*s212+mat12*(s222-d2(0))+mat13*s232)/det
c---  >   d2(1,3)
      d2(5) = (mat11*s213+mat12*s223+mat13*(s233-d2(0)))/det
c---  >   d2(2,3)=d2(2,1) with m2<-->m4 for square boxes
      d2(6) = (mat21*s213+mat22*s223+mat23*(s233-d2(0)))/det
      s310 = (c2134(0)-c2234(0)-cf1*d2(0))/2d0
      s320 = (c2124(0)-c2234(0)-cf2*d2(0))/2d0
      s330 = (c2123(0)-c2234(0)-cf3*d2(0))/2d0
      ccc = ccc+(c2234(1)+c2234(2)-c1234(0))/2d0+c2234(3)
      s311 = -ccc-cf1*d2(1)/2d0
      s321 = (c2124(1)-cf2*d2(1))/2d0-ccc
      s331 = (c2123(1)-cf3*d2(1))/2d0-ccc
      s312 = (c2134(1)-c2234(1)-cf1*d2(2))/2d0
      s322 = -(c2234(1)+cf2*d2(2))/2d0
      s332 = (c2123(2)-c2234(1)-cf3*d2(2))/2d0
      s313 = (c2134(2)-c2234(2)-cf1*d2(3))/2d0
      s323 = (c2124(2)-c2234(2)-cf2*d2(3))/2d0
      s333 = -(c2234(2)+cf3*d2(3))/2d0
      s3113 = (c2234(3)+c2234(2)+c1234(2)-cf1*d2(5))/2d0
      s3213 = (c2124(3)+c2234(3)+c2234(2)+c1234(2)-cf2*d2(5))/2d0
      s3313 = (c2234(3)+c2234(2)+c1234(2)-cf3*d2(5))/2d0
      d3(0,1) = (mat11*s310+mat12*s320+mat13*s330)/det
      d3(0,2) = (mat21*s310+mat22*s320+mat23*s330)/det
      d3(0,3) = (mat31*s310+mat32*s320+mat33*s330)/det
      d3(1,1) = (mat11*(s311-2d0*d3(0,1))+mat12*s321+mat13*s331)/det
      d3(3,3) = (mat31*s313+mat32*s323+mat33*(s333-2d0*d3(0,3)))/det
      d3(1,2) = (mat21*(s311-2d0*d3(0,1))+mat22*s321+mat23*s331)/det
      d3(1,3) = (mat31*(s311-2d0*d3(0,1))+mat32*s321+mat33*s331)/det
      d3(2,1) = (mat11*s312+mat12*(s322-2d0*d3(0,2))+mat13*s332)/det
      d3(2,2) = (mat21*s312+mat22*(s322-2d0*d3(0,2))+mat23*s332)/det
      d3(2,3) = (mat31*s312+mat32*(s322-2d0*d3(0,2))+mat33*s332)/det
      d3(3,1) = (mat11*s313+mat12*s323+mat13*(s333-2d0*d3(0,3)))/det
      d3(3,2) = (mat21*s313+mat22*s323+mat23*(s333-2d0*d3(0,3)))/det
      d3(0,0) = (mat21*(s3113-d3(0,3))+mat22*s3213+mat23*
     1     (s3313-d3(0,1)))/det
c      write(6,*)q1,q2,q3,q6,q7,q4,q7,m1,m2,m3,m4
c      write(6,*)d0new,d2(0)
c      write(6,*)d1
c      write(6,*)d2
c      write(6,*)d3
      return
      end
***********************************************************************
      complex*16 function fs(q2,rm,rn)
***********************************************************************
*       skalares einschleifenintegral, doppeltlang, 'regulaerer anteil'
*       f(q2,rm,rn)=b0(q2,rm,rn)-b0(0d0,rm,rn)  'subtrahiertes f'
*       q2=quadrat des die schleife durchlaufenden impulses
*       rm,rn: massen der teilchen auf beiden armen
*       d o u b l e p r e c i s i o n
*-----------------------------------------------------------------------
*       19.10.83
************************************************************************
      real*8 m,n,pi,a,s,t,b,c,d,q2,rm,rn,u,v,w,m2,n2
      data pi/3.1415926535897932384626438d0/
      m=rm
      n=rn
      m2 = m**2
      n2 = n**2
      if (dabs(m) .eq. dabs(n) ) goto 30
      if (n2 .eq. 0d0) goto 310
      if (m2 .eq. 0d0) goto 300
c---- >  allgemeiner fall
      if (q2 .ne. 0d0) goto 520
      b=0d0
      a=0d0
      goto 560
 520  u=m*m+n*n
      v=m*m-n*n
      w=m*n
      if (dabs(q2/v).le.1d-4) then
         b = u/v/v/2d0+2d0*w*w/v/v/v*dlog(n2/m2)/2d0
         b = b*q2
         a = 0d0
         goto 570
      end if
      s=dabs(m)+dabs(n)
      t=dabs(m)-dabs(n)
      c=dsqrt(dabs(s*s-q2))
      d=dsqrt(dabs(t*t-q2))
      b=1d0+(v/q2-u/v)*dlog(n2/m2)/2d0
      if (2d0*w .le. dabs(q2-u)) goto 550
      b=b-2d0*c*d/q2*datan(d/c)
      a=0d0
      goto 560
 550  a=c*d/q2
      b=b-dsign(1d0,q2-u)*a*dlog((c+d)*(c+d)/(4d0*w))
      a=pi*a
      if (q2 .ge. u) goto 560
      a=0d0
 560  continue
 570  fs=dcmplx(b,a)
      return
c---- >  gleiche massen
 30   if (q2 .ne. 0d0) goto 40
      b=0d0
      a=0d0
      goto 560
 40   u=4d0*m*m
      v=dsqrt(dabs(1d0-u/q2))
      if ((q2 .ge. 0d0) .and. (q2 .lt. u)) goto 50
      b=2d0-v*dlog((v+1d0)*(v+1d0)/u*dabs(q2))
      a=pi*v
      if (q2 .ge. u) goto 560
      a=0d0
      goto 560
 50   b=2d0-2d0*v*datan(1d0/v)
      a=0d0
      goto 560
c---- >  eine masse null
 300  m=n
 310  if (q2 .ne. m*m) goto 320
      a=0d0
      b=1d0
      goto 560
 320  b=1d0
      if(q2 .eq. 0d0) b=0d0
      a=b*(1d0-m*m/(q2+(1d0-b)))
      b=b-a*dlog(dabs(1d0-q2/m/m))
      a=pi*a
      if (q2 .gt. m*m) goto 560
      a=0d0
      goto 560
      end
***********************************************************************
      complex*16 function b0(q2,rm,rn)
***********************************************************************
*       skalares einschleifenintegral b0 minus
*       divergenter anteil (only (delta-log mue**2) subtracted)
*       rm,rn: massen der teilchen auf beiden armen
*       d o u b l e p r e c i s i o n
*-----------------------------------------------------------------------
*       04.02.87 sa
************************************************************************
      implicit real*8(a,c,g,q,m,s,r)
      complex*16 fs
*      
      b0 = dcmplx(0d0)
      r12 = rm*rm
      r22 = rn*rn
      if (dabs(rm).eq.dabs(rn)) goto 110
      if ((r22.eq.0d0).or.(r12.eq.0d0)) goto 120
      b0 = dcmplx(-dlog(r12)/2d0-dlog(r22)/2d0+1d0-(r12+r22)/(r12-r22)
     1     *dlog(r12/r22)/2d0,0d0)
      goto 200
c---  >   beide massen gleich
 110  if (r12.eq.0d0) goto 300
      b0 = dcmplx(-dlog(r12))
      goto 200
c---  >   eine masse null
 120  b0 = 1d0-dcmplx(dlog(r12+r22))
 200  b0 = b0+fs(q2,rm,rn)
      return
c---  >   beide massen gleich null
 300  b0 = 2d0-cdlog(dcmplx(q2,1d-20))
      return
      end
***********************************************************************
      complex*16 function b1(q2,rm,rn)
***********************************************************************
*     skalares einschleifenintegral b1 minus
*     divergenter anteil (only -1/2(delta-log(mue**2)) subtracted)
*     rm,rn: massen der teilchen auf beiden armen
*     d o u b l e p r e c i s i o n
*-----------------------------------------------------------------------
*       04.02.87 sa
************************************************************************
      implicit real*8(a,c,g,q,m,s,r)
      complex*16 fs,b0,xb0,xb1
*
      r12= rm*rm
      r22= rn*rn
      b1 = dcmplx(0d0)
      if (dabs(rm).eq.dabs(rn)) goto 200
      if (r12.eq.0d0) goto 120
      if (r22.eq.0d0) goto 130
      if (q2.eq.0d0) goto 140
      b1 = (r22-r12)*b0(q2,rm,rn)+r12*(1d0-dlog(r12))
     1     -r22*(1d0-dlog(r22))
      b1 = b1/2d0/q2
      goto 200
 120  b1 = r22*fs(q2,rm,rn)/2d0/q2
      goto 200
 130  b1 = -r12*fs(q2,rm,rn)/2d0/q2
 200  b1 = b1-b0(q2,rm,rn)/2d0
      goto 300
c Achtung: original code has been changed here!
 140  call bquer(q2,rm,rn,xb0,xb1)
      b1 = xb1+dlog(r22)/2d0
 300  continue
      return
      end
***********************************************************************
      complex*16 function b20(q2,rm,rn)
***********************************************************************
*       skalares einschleifenintegral b0 minus
*       divergenter anteil
*       rm,rn: massen der teilchen auf beiden armen
*       d o u b l e p r e c i s i o n
*-----------------------------------------------------------------------
*       04.02.87 sa
************************************************************************
      implicit real*8(r,m,q)
      complex*16 b0,b1,b21
*      
      r12= rm*rm
      r22= rn*rn
      b20 = (q2+r12-r22)*b1(q2,rm,rn)+2d0*r12*b0(q2,rm,rn)
     1     +r12+2d0*r22-q2/3d0
      if (r22.ne.0d0) b20 = b20-r22*dlog(r22)
      b20 = b20/6d0
      return
      entry b21(q2,rm,rn)
      r12= rm*rm
      if ((dabs(rm).eq.dabs(rn)).and.(q2.eq.0d0)) goto 100
      r22= rn*rn
      b21 = 2d0*(r12-r22+q2)*b1(q2,rm,rn)+r12*b0(q2,rm,rn)
     1     +(r12-r22-q2/3d0)/2d0
      if (r22.ne.0d0) b21 = b21+r22*dlog(r22)
      b21 = -b21/3d0/q2
      return
 100  b21 = -dlog(r12)/3d0
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function c0(cq1,cq2,cq3,mm1,mm2,mm3)
***********************************************************************
*       skalare dreipunktfunktion, doppeltlang,
*       q1,q2,q3=quadrate der impulse
*       m1,m2,m3: massenquadrate der teilchen auf inneren linien
*       d o u b l e p r e c i s i o n
*-----------------------------------------------------------------------
*       20.01.87 sa
************************************************************************
      implicit real*8(m,q,c)
      implicit complex*16(z)
      complex*16 y(3,3),spence,cywur
      real*8 cj(3)
*
      j = 3
      cj(1) = -1d0
      cj(2) = 1d0
      cj(3) = -1d0
      m1 = mm1*mm1
      m2 = mm2*mm2
      m3 = mm3*mm3
      q1 = cq1
      q2 = cq2
      q3 = cq3
      if (q2 .eq. 0d0) then
         q2 = cq1
         q1 = cq2
         m3 = mm1*mm1
         m1 = mm3*mm3
      end if
      if (q1 .eq. 0d0) then
         j = 2
      end if
      c1 = 1d0+(q1-q3)/q2
      if (c1*c1-4d0*q1/q2 .lt. 0d0) then
         write(6,*) ' neg. determinante !'
      end if
      cwurz = dsqrt(c1*c1-4d0*q1/q2)
      calphs = (c1-cwurz)/2d0
      c0 = 0d0
      ca = q1
      cb = q2
      cc = q3-q2-q1
      cd = m2-m1-q1
      ce = q1-q3+m3-m2
      cf = m1
      cnum = -(cd+ce*calphs)
      cdenom = cc+2d0*calphs*cb
      y(1,3) = dcmplx((cnum-2d0*ca-cc*calphs)/cdenom,0d0)
      cy0 = -cc-ce
      cywur = cdsqrt(dcmplx(cy0*cy0-4d0*cb*(ca+cd+cf),1d-20))
      y(1,1) = (cy0+cywur)/2d0/cb
      y(1,2) = (cy0-cywur)/2d0/cb
      cy0 = -ce-cd
      cywur = cdsqrt(dcmplx(cy0*cy0-4d0*cf*(ca+cb+cc),1d-20))
      y(2,3) = dcmplx(cnum/cdenom/(1d0-calphs),0d0)
      y(2,1) = (cy0+cywur)/2d0/(ca+cb+cc)
      y(2,2) = (cy0-cywur)/2d0/(ca+cb+cc)
      cywur = cdsqrt(dcmplx(cd*cd-4d0*ca*cf,1d-20))
      if (j .eq. 3) then
         y(3,3) = dcmplx(-cnum/calphs/cdenom,0d0)
         y(3,1) = (-cd+cywur)/2d0/ca
         y(3,2) = (-cd-cywur)/2d0/ca
      end if
      do 100 i=1,j
         c0 = c0+cj(i)*dreal(spence(y(i,3)/(y(i,3)-y(i,1)))
     1        -spence((y(i,3)-1d0)/(y(i,3)-y(i,1)))
     2        +spence(y(i,3)/(y(i,3)-y(i,2)))
     3        -spence((y(i,3)-1d0)/(y(i,3)-y(i,2))))
 100  continue
      c0 = c0/cdenom
      return
      end
***********************************************************************
      subroutine cmue(q1,q2,q3,m1,m2,m3,c1)
***********************************************************************
*     nicht-skalare dreipunktfunktion, doppeltlang,
*     q1,q2,q3=quadrate der impulse
*     m1,m2,m3: massen der teilchen auf inneren linien
*     nicht ir divergent
*     d o u b l e p r e c i s i o n
*-----------------------------------------------------------------------
*       20.01.87 sa
************************************************************************
      implicit real*8(a,c,g,q,m,s,d)
      implicit complex*16(z)
      complex*16 b0,spence,cscal,C0_
      real*8 c1(0:2)
*
      m12 = m1*m1
      m22 = m2*m2
      m32 = m3*m3
c--->   two gluons on shell, 3 equal internal fermions
      if ((q1.eq.0d0).and.(q2.eq.0d0).and.(q3.gt.0d0)) goto 230
      goto 240
 230  if (4d0*m12.ge.q3) then
         awur = dsqrt(4d0*m12/q3-1d0)
         z1 = dcmplx(1d0,awur)/2d0
         z2 = dcmplx(1d0,-awur)/2d0
      else
         awur = dsqrt(1d0-4d0*m12/q3)
         al1 = 1d0+awur
         al2 = 1d0-awur
         z1 = dcmplx(al1,al2*1d-10)/2d0
         z2 = dcmplx(al2,-al2*1d-10)/2d0
      end if
      c1(0) = -spence(1d0/z1)-spence(1d0/z2)
      c1(0) = c1(0)/q3
c     wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
c     write(6,*)'Check this point in ''boxlib.f'' at line 408'
c      write(6,*)'                                       ---'
c     wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
c new C0 function added on Aug.21 2004
c      c1(0) = dreal(cscal(q3,0d0,m1,m3,m2))
c      write(6,*)'1',q3,m1,m2,m3,c1(0)
      c1(0)=dreal(C0_(0d0,q3,0d0,m2,m1,m3,0))
c      write(6,*)'2',c1(0)
      goto 100
 240  qq2 = q2
      qq1 = q1
      qq3 = q3
      mm1 = m1
      mm2 = m2
      mm3 = m3
      if ((dabs(q1).gt.0d0).and.(q3.eq.0d0)) then
         qq3 = q1
         qq1 = q3
         mm2 = m3
         mm3 = m2
      end if
 300  c1(0) = c0(qq1,qq2,qq3,mm1,mm2,mm3)
      goto 100
 100  matr = (q3-q2-q1)/2d0
      det = q1*q2-matr*matr
      cb3 = dreal(b0(q3,m1,m3))
      cb2 = dreal(b0(q2,m2,m3))
      cb1 = dreal(b0(q1,m1,m2))
      cmp1 = (cb3-cb2-(m12-m22+q1)*c1(0))/2d0
      cmp2 = (cb1-cb3-(m22-m32+q3-q1)*c1(0))/2d0
      c1(1) = (q2*cmp1-matr*cmp2)/det
      c1(2) = (q1*cmp2-matr*cmp1)/det
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cmn(q1,q2,q3,m1,m2,m3,c1,c2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     dreipunktfkt. mit 2 integrationsimpulsen in zaehler
c     realteil
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     04.02.87 sa
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8(c,q,m)
      complex*16 b0,b1
      real*8 det,c1(0:2),c2(0:3)
*      
      m12 = m1*m1
      m22 = m2*m2
      m32 = m3*m3
      matr = (q3-q2-q1)/2d0
      det = q1*q2-matr*matr
      cf1 = m12-m22+q1
      cf2 = m22-m32+q3-q1
      cb0 = dreal(b0(q2,m2,m3))
      call cmue(q1,q2,q3,m1,m2,m3,c1)
      c2(0) = m12/2d0*c1(0)+(cb0+cf1*c1(1)+cf2*c1(2)+1d0)/4d0
      cb1 = dreal(b1(q3,m1,m3))
      cb2 = dreal(b1(q1,m1,m2))
      cmp1 = cb1+cb0-cf1*c1(1)-2d0*c2(0)
      cmp2 = cb2-cb1-cf2*c1(1)
      c2(1) = q2*cmp1-matr*cmp2
      c2(1) = c2(1)/2d0/det
      c2(3) =-matr*cmp1+q1*cmp2
      c2(3) = c2(3)/2d0/det
      cb3 = b1(q2,m2,m3)
      cmp1 = cb1-cb3-cf1*c1(2)
      cmp2 = -cb1-cf2*c1(2)-2d0*c2(0)
      c2(2) = -matr*cmp1+q1*cmp2
      c2(2) = c2(2)/2d0/det
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine chutmn(q1,q2,q3,m1,m2,m3,ch1,ch2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     dreipunktfkt. mit 2 integrationsimpulsen in zaehler
c     realteil
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     06.05.87 sa
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8(c,q,m)
      real*8 c1(0:2),c2(0:3),ch1(0:2),ch2(0:3)
*
      call cmn(q1,q2,q3,m1,m2,m3,c1,c2)
      ch1(0) = c1(0)
      ch1(1) = c1(1)-c1(2)
      ch1(2) = c1(2)
      ch2(0) = c2(0)
      ch2(1) = c2(1)-2d0*c2(3)+c2(2)
      ch2(2) = c2(2)
      ch2(3) = c2(3)-c2(2)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      complex*16 function d0gen(q1,q2,q3,q4,q5,q6,m1,m2,m3,m4)
***********************************************************************
*       skalare vierpunktfunktion, doppeltlang,
*       za-zf: parameters a-f from 't hooft veltman
*       d o u b l e p r e c i s i o n
*-----------------------------------------------------------------------
*       07.08.90 sa
************************************************************************
      implicit real*8(m,q,p)
      implicit complex*16(c,z)
      complex*16 cm2(4),cmm2(4),cl(4,4),cq(4,4),ca(4),ca2(4)
      real*8 p(4,4)
*
      cm2(1) = dcmplx(m1*m1,-1d-20)
      cm2(2) = dcmplx(m2*m2,-1d-20)
      cm2(3) = dcmplx(m3*m3,-1d-20)
      cm2(4) = dcmplx(m4*m4,-1d-20)
      p(1,2) = q1
      p(1,3) = q4
      p(1,4) = q6
      p(2,3) = q2
      p(2,4) = q1+q2+q3+q6-q4-q5
      p(3,4) = q3
      do 110 i1=1,4
         do 100 i2=i1+1,4
            cl(i1,i2) = p(i1,i2)-cm2(i1)-cm2(i2)
 100     continue
 110  continue
      ca(2) = dcmplx(1d0/m1/m1)
      ca2(2) = ca(2)*ca(2)
      ca(1) = cl(1,2)-cdsqrt(cl(1,2)*cl(1,2)-4d0*cm2(1)*cm2(2))
      ca(1) = -ca(1)/2d0/cm2(1)*ca(2)
      ca2(1) = ca(1)*ca(1)
      ca(3) = (-cm2(2)*ca2(2)+cm2(1)*ca2(1))
     1     /(cl(2,3)*ca(2)-cl(1,3)*ca(1))
      ca2(3) = ca(3)*ca(3)
      ca(4) = (-cm2(2)*ca2(2)+cm2(1)*ca2(1))
     1     /(cl(2,4)*ca(2)-cl(1,4)*ca(1))
      ca2(4) = ca(4)*ca(4)
      do 210 i1=1,4
         cmm2(i1) = -cm2(i1)*ca2(i1)
         do 200 i2=i1+1,4
            cq(i1,i2) = cl(i1,i2)*ca(i1)*ca(i2)+cm2(i1)*ca2(i1)
     1           +cm2(i2)*ca2(i2)
 200     continue
 210  continue
      caa = -cq(3,4)
      cb = -cq(2,3)
      cg = -cq(1,2)
      cc = -cq(2,4)+cq(2,3)+cq(3,4)
      ch = -cq(1,4)-cq(2,3)+cq(1,3)+cq(2,4)
      cj = -cq(1,3)+cq(1,2)+cq(2,3)
      cd = cmm2(3)-cmm2(4)+cq(3,4)
      ce = cmm2(2)-cmm2(3)+cq(2,4)-cq(3,4)
      ck = cmm2(1)-cmm2(2)+cq(1,4)-cq(2,4)
      cf = cmm2(4)
      d0gen = c0gen(caa,cb,cc,cd,ce,cf)
     1     -c0gen(caa,cb,cc,cd,(ce+ck),cf)
      d0gen = -ca(1)*ca(2)*ca(3)*ca(4)*d0gen/ck
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      complex*16 function c0gen(za,zb,zc,zd,ze,zf)
***********************************************************************
*       skalare dreipunktfunktion, doppeltlang,
*       za-zf: parameters a-f from 't hooft veltman
*       d o u b l e p r e c i s i o n
*-----------------------------------------------------------------------
*       07.08.90 sa
************************************************************************
      implicit real*8(m,q)
      implicit complex*16 (c,z)
      complex*16 y(3,3),c1(3),c2(3),yd1,yd2
      complex*16 spence
      real*8 cj(3)
      j = 3
      cj(1) = -1d0
      cj(2) = 1d0
      cj(3) = -1d0
*      
      cwurz = cdsqrt(zc*zc-4d0*za*zb)
      calphs = (-zc-cwurz)/2d0/zb
      c0 = dcmplx(0d0)
      ca = za
      cb = zb
      cc = zc
      cd = zd
      ce = ze
      cf = zf
      cnum = -(cd+ce*calphs)
      cdenom = cc+2d0*calphs*cb
      y(1,3) = (cnum-2d0*ca-cc*calphs)/cdenom
      cy0 = -cc-ce
      cywur = cdsqrt(cy0*cy0-4d0*cb*(ca+cd+cf))
      y(1,1) = (cy0+cywur)/2d0/cb
      y(1,2) = (cy0-cywur)/2d0/cb
      c1(1) = cb
      c2(1) = ca+cd+cf
      cy0 = -ce-cd
      cywur = cdsqrt(cy0*cy0-4d0*cf*(ca+cb+cc))
      y(2,3) = cnum/cdenom/(1d0-calphs)
      y(2,1) = (cy0+cywur)/2d0/(ca+cb+cc)
      y(2,2) = (cy0-cywur)/2d0/(ca+cb+cc)
      c1(2) = ca+cb+cc
      c2(1) = cf
      cywur = cdsqrt(cd*cd-4d0*ca*cf)
      y(3,3) = -cnum/calphs/cdenom
      y(3,1) = (-cd+cywur)/2d0/ca
      y(3,2) = (-cd-cywur)/2d0/ca
      c1(3) = ca
      c2(3) = cf
      do 100 i=1,3
         yd1 = y(i,3)-y(i,1)
         yd2 = y(i,3)-y(i,2)
         c0 = c0+cj(i)*(spence(y(i,3)/yd1)
     1        -spence((y(i,3)-1d0)/yd1)
     2        +spence(y(i,3)/yd2)
     3        -spence((y(i,3)-1d0)/yd2)
     5        )
 100  continue
      c0gen = c0/cdenom
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function ggttb(s,t,mex2,mib2,mif2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       mass singular scalar box integral gg-->tt
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       21.06.90 sa
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8(a-b,d-y)
      implicit complex*16(c,z)
      complex*16 spence
*
      pi = 4d0*datan(1d0)
      cs = dcmplx(s,s*1d-10)
      cmh2 = dcmplx(mib2,-mib2*1d-10)
      tmt = t-mex2
      cmht = cmh2-t
      cmhmt = cmh2-mex2
      mht = mib2-t
      mhmt = mib2-mex2
      c1 = -tmt/cmht
      c2 = -s*cmh2/cmhmt/cmhmt
      zlogi = dcmplx(dlog(mif2/s),pi)
      zlogh = cdlog(cmhmt)
      zlogs = dcmplx(dlog(mib2/s),pi)
      log1 = dlog(mht)
      log2 = dlog(mib2)
      log3 = dlog((s*mib2+mhmt*mhmt)/mhmt/mhmt)
      cd0 = 4d0*spence(c1)+spence(c2)
     1     -zlogi*(2d0*zlogh-2d0*log1+zlogi/2d0)
     2     +log3*(2d0*log2-2d0*zlogh-zlogs)
      ggttb = dreal(cd0)/s/mht
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function ggttf(s,t,mex2,mib2,mif2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       scalar box integral gg-->tt
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       11.07.90 sa
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8(a-b,d-y)
      implicit complex*16(c,z)
      complex*16 spence
*
      eps1 = 1d-25
      eps = mib2*eps1
      cs = dcmplx(s,eps1*s)
      ct = dcmplx(t,-eps1*t)
      cmh2 = dcmplx(mib2,-eps)
      ctm = ct-mex2
      tm = dreal(ctm)
      tmth = t+mex2-mib2
      tmth2 = tmth*tmth
      st = s*t+tm*tm
      wur1 = dsqrt(tmth2-4d0*mex2/s*st)
      w12 = wur1*s
      y1 = (tmth+wur1)/2d0*s/st
      y2 = (tmth-wur1)/2d0*s/st
      cy1 = dcmplx(y1,eps1)
      cy2 = dcmplx(y2,-eps1)
      y11 = (y1-1d0)/y1
      y21 = (y2-1d0)/y2
      wur3 = dsqrt(tmth2-4d0*mex2*t)
      cy3 = dcmplx(tmth+wur3,eps)/2d0/t
      cy4 = dcmplx(tmth-wur3,-eps)/2d0/t
      cy31 = (cy3-1d0)/cy3
      cy41 = (cy4-1d0)/cy4
      cy13 = y1-cy3
      cy14 = y1-cy4
      cy23 = y2-cy3
      cy24 = y2-cy4
      mth = 2d0*mex2-mib2
      if (4d0*mex2.le.mib2) then
         wur5 = mib2*dsqrt(1d0-4d0*mex2/mib2)
         cy5 = dcmplx(mth+wur5,eps)/2d0/mex2
         cy6 = dcmplx(mth-wur5,-eps)/2d0/mex2
      else
         wur5 = mib2*dsqrt(4d0*mex2/mib2-1d0)
         cy5 = dcmplx(mth,wur5)/2d0/mex2
         cy6 = dcmplx(mth,-wur5)/2d0/mex2
      end if
      cy51 = (cy5-1d0)/cy5
      cy61 = (cy6-1d0)/cy6
      cy15 = y1-cy5
      cy16 = y1-cy6
      cy25 = y2-cy5
      cy26 = y2-cy6
      cy7 = s/(s+ctm)
      cy17 = y1-cy7
      cy27 = y2-cy7
      cd0 =
     1     -spence((cy1-1d0)/cy1)
     2     -spence(y1/cy17)+spence((y1-1d0)/cy17)
     3     -spence(y1/cy13)+spence((y1-1d0)/cy13)
     4     -spence(y1/cy14)+spence((y1-1d0)/cy14)
     5     +spence(y1/cy15)-spence((y1-1d0)/cy15)
     6     +spence(y1/cy16)-spence((y1-1d0)/cy16)
     7     +dlog((y1-1d0)/y1)*( dlog(tm*mex2/t/(s+tm))+cdlog(cy1)
     8     -cdlog(cy17)-cdlog(cy13*cy14)+cdlog(cy15*cy16) )
      cd0 = cd0-(
     1     -spence((cy2-1d0)/cy2)
     2     -spence(y2/cy27)+spence((y2-1d0)/cy27)
     3     -spence(y2/cy23)+spence((y2-1d0)/cy23)
     4     -spence(y2/cy24)+spence((y2-1d0)/cy24)
     5     +spence(y2/cy25)-spence((y2-1d0)/cy25)
     6     +spence(y2/cy26)-spence((y2-1d0)/cy26)
     7     +dlog((y2-1d0)/y2)*( dlog(tm*mex2/t/(s+tm))+cdlog(cy2)
     8     -cdlog(cy27)-cdlog(cy23*cy24)+cdlog(cy25*cy26) )
     9     )
*     
      beta = dsqrt(1d0-4d0*mex2/s)
      al1 = (1d0+beta)/2d0
      al2 = (1d0-beta)/2d0
      sat = al1*s+tm
      rz1 = (-beta*tm+mib2+wur1)/2d0/sat
      rz2 = (-beta*tm+mib2-wur1)/2d0/sat
      z1 = dcmplx(rz1,eps1)
      z2 = dcmplx(rz2,-eps1)
      z3 = cmplx(al1,eps1)
      z4 = cmplx(al2,-eps1)
      z13 = rz1-z3
      z14 = rz1-z4
      z23 = rz2-z3
      z24 = rz2-z4
      z5 = -cy5+1d0
      z6 = -cy6+1d0
      z15 = rz1-al2*z5
      z16 = rz1-al2*z6
      z25 = rz2-al2*z5
      z26 = rz2-al2*z6
      cd0 = cd0+(
     1     -spence((rz1-al2)/z15)+spence(rz1/z15)
     2     -spence((rz1-al2)/z16)+spence(rz1/z16)
     3     -spence(z1/(z1-al2))
     7     -cdlog(z1/(z1-al2))*( dlog(-mex2/tm/al2)-cdlog(al2-z1)
     8     +cdlog(z15*z16) )
     7     )
      cd0 = cd0-(
     1     -spence((rz2-al2)/z25)+spence(rz2/z25)
     2     -spence((rz2-al2)/z26)+spence(rz2/z26)
     3     -spence(z2/(z2-al2))
     7     -cdlog(z2/(z2-al2))*( dlog(-mex2/tm/al2)-cdlog(al2-z2)
     8     +cdlog(z25*z26) )
     7     )
      z15 = rz1+al1*z5
      z16 = rz1+al1*z6
      z25 = rz2+al1*z5
      z26 = rz2+al1*z6
      z17 = rz1+al1*(1d0-cy7)
      z27 = rz2+al1*(1d0-cy7)
      zst = dcmplx(s+tm,eps)
      cd0 = cd0+(
     1     +spence((rz1+al1)/z15)-spence(rz1/z15)
     2     +spence((rz1+al1)/z16)-spence(rz1/z16)
     3     -spence((rz1+al1)/z17)+spence(rz1/z17)
     7     +cdlog(z1/(z1+al1))*( cdlog(mex2/zst/al1)
     8     -cdlog(-z17)+cdlog(z15*z16) )
     7     )
      cd0 = cd0-(
     1     +spence((rz2+al1)/z25)-spence(rz2/z25)
     2     +spence((rz2+al1)/z26)-spence(rz2/z26)
     3     -spence((rz2+al1)/z27)+spence(rz2/z27)
     7     +cdlog(z2/(z2+al1))*( cdlog(mex2/zst/al1)
     8     -cdlog(-z27)+cdlog(z25*z26) )
     7     )
      z15 = rz1+z3
      z25 = rz2+z3
      z16 = z15-1d0
      z26 = z25-1d0
*     
      cd0 = cd0+(
     1     +spence(z15/z16)
     2     +spence(z16/(z16+al1))-spence(z15/(z16+al1))
     3     +spence(z16/z1)-spence(z15/z1)
     7     -cdlog(z15/z16)*( cdlog(-z16)-cdlog(-z1)
     8     -cdlog(-z16-al1) )
     7     )
      cd0 = cd0-(
     1     +spence(z25/z26)
     2     +spence(z26/(z26+al1))-spence(z25/(z26+al1))
     3     +spence(z26/z2)-spence(z25/z2)
     7     -cdlog(z25/z26)*( cdlog(-z26)-cdlog(-z2)
     8     -cdlog(-z26-al1) )
     7     )
*     
 999  ggttf = dreal(cd0)/w12
      return
      end
