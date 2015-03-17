c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine ddfunc(xpi1,xpi2,xpi3,pi12,pi13,pi23,
     $     m0,mi1,mi2,mi3,d1,d2,d3,c1,c2)
*---------------------------------------------------------------------
c entwicklungskoeff. fuer das vektorielle und tensorielle vierpkt.-integral
c d0=k^2-m0^2, d1=(k+pi1)^2-m1^2, d2=(k+pi2)^2-m2^2, d3=(k+pi3)^2-m3^2
c xpi1=pi1^2,xpi2=pi2^2,xpi3=pi3^2,pi12=pi1.pi2,pi13=pi1.pi3,pi23=pi2.pi3
c Up to rank 3 !
c Only IR finite D-functions ! 
*---------------------------------------------------------------------
      implicit none
      integer i
      real*8 xpi1,xpi2,xpi3,pi12,pi13,pi23,q12,q13,q23,
     $     mi1,mi2,mi3,m0,xmi1,xmi2,xmi3,xm0,cf1,cf2,cf3,det
      real*8 mat11,mat12,mat13,mat21,mat22,mat23,mat31,mat32,mat33
      complex*16 d1(0:3),d2(0:6),d3(0:3,0:3),c1(0:2,4),c2(0:3,4)
      complex*16 D0_
      complex*16 c1234(0:2),c2234(0:3),c1134(0:2),c2134(0:3),
     $     c1124(0:2),c2124(0:3),c1123(0:2),c2123(0:3),ccc
      complex*16 s11,s12,s13,s20,s211,s212,s213,s221,s222,s223,s231,
     $     s232,s233,s311,s312,s313,s321,s322,s323,s331,s332,s333,
     $     s3113,s3213,s3313,s310,s320,s330
      xm0 = m0*m0
      xmi1 = mi1*mi1
      xmi2 = mi2*mi2
      xmi3 = mi3*mi3
      q12=xpi1+xpi2-2d0*pi12
      q13=xpi1+xpi3-2d0*pi13
      q23=xpi2+xpi3-2d0*pi23
      cf1=xpi1-xmi1+xm0
      cf2=xpi2-xmi2+xm0
      cf3=xpi3-xmi3+xm0
c Gram determinant
      det=xpi1*xpi2*xpi3-xpi1*pi23**2-xpi2*pi13**2-xpi3*pi12**2+
     $     2d0*pi12*pi13*pi23
c      if(dabs(det).lt.1d-1)then
c         write(6,*)'warning: potential instability in ddfunc'
c         write(6,*)xpi1,xpi2,xpi3,pi23,pi13,pi12
c         stop
c      endif
c matrix elements of the inverse Gram matrix x det
      mat11=xpi2*xpi3-pi23**2
      mat12=pi13*pi23-pi12*xpi3
      mat13=pi12*pi23-pi13*xpi2
      mat21=mat12
      mat22=xpi1*xpi3-pi13**2
      mat23=pi12*pi13-pi23*xpi1
      mat31=mat13
      mat32=mat23
      mat33=xpi1*xpi2-pi12**2
*
      d1(0)=D0_(xpi1,q12,q23,xpi3,xpi2,q13,
     $     m0,mi1,mi2,mi3,0)
      call ccfunc2(q12,q13,pi23-pi12-pi13+xpi1,mi1,mi2,mi3,
     $     c1234(0),c1234(1),c1234(2),
     $     c2234(1),c2234(2),c2234(0),c2234(3))
      call ccfunc2(xpi2,xpi3,pi23,m0,mi2,mi3,
     $     c1134(0),c1134(1),c1134(2),
     $     c2134(1),c2134(2),c2134(0),c2134(3))
      call ccfunc2(xpi1,xpi3,pi13,m0,mi1,mi3,
     $     c1124(0),c1124(1),c1124(2),
     $     c2124(1),c2124(2),c2124(0),c2124(3))
      call ccfunc2(xpi1,xpi2,pi12,m0,mi1,mi2,
     $     c1123(0),c1123(1),c1123(2),
     $     c2123(1),c2123(2),c2123(0),c2123(3))
      do i=0,2
         c1(i,1)=c1234(i)
         c1(i,2)=c1134(i)
         c1(i,3)=c1124(i)
         c1(i,4)=c1123(i)
      end do
      do i=0,3
         c2(i,1)=c2234(i)
         c2(i,2)=c2134(i)
         c2(i,3)=c2124(i)
         c2(i,4)=c2123(i)
      end do
*---------------------------------------------------------------------
*--entwicklungskoeff. des vektorintegrals:----------------------------
      s11 = (c1134(0)-c1234(0)-cf1*d1(0))/2d0
      s12 = (c1124(0)-c1234(0)-cf2*d1(0))/2d0
      s13 = (c1123(0)-c1234(0)-cf3*d1(0))/2d0
      d1(1) = (mat11*s11+mat12*s12+mat13*s13)/det
      d1(2) = (mat21*s11+mat22*s12+mat23*s13)/det
      d1(3) = (mat31*s11+mat32*s12+mat33*s13)/det
*---------------------------------------------------------------------
*--entwicklungskoeff. des tensorintegrals (Rank 2):-------------------
      ccc = c1234(1)+c1234(2)+c1234(0)
      s20 = c1234(0)+xm0*d1(0)
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
*---------------------------------------------------------------------
*--entwicklungskoeff. des tensorintegrals (Rank 3):-------------------
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
 99   continue
      return
      end
c---------------------------------------------------------------------
      subroutine ccfunc2(xpi1,xpi2,pi12,m0,mi1,mi2,
     $     c0,c11,c12,c21,c22,c20,c23)
*---------------------------------------------------------------------
c entwicklungskoeff. fuer das vektorielle und tensorielle dreipkt.-integral
c d0=k^2-m0^2, d1=(k+pi1)^2-m1^2, d2=(k+pi2)^2-m2^2
c xpi1=pi1^2,xpi2=pi2^2,pi12=pi1.pi2
*---------------------------------------------------------------------
      implicit none
      real*8 xpi1,xpi2,pi12,mi1,mi2,m0,xmi1,xmi2,xm0,f1,f2,q12,det
      complex*16 c0,c11,c12,c21,c22,c20,c23,
     &     b01,b02,b03,b11,b12,b13,b20,C0_
      xm0 = m0*m0
      xmi1 = mi1*mi1
      xmi2 = mi2*mi2
      q12=xpi1+xpi2-2d0*pi12
      f1=xpi1-xmi1+xm0
      f2=xpi2-xmi2+xm0
      det=xpi1*xpi2-pi12**2
c      if(dabs(det).lt.1d-1)then
c         write(6,*)'warning: potential instability in ccfunc2'
c         stop
c      endif
      if(dabs(q12).ne.0d0.and.dabs(q12).lt.1d-6)q12=0d0
      if(dabs(xpi1).ne.0d0.and.dabs(xpi1).lt.1d-6)xpi1=0d0
      if(dabs(xpi2).ne.0d0.and.dabs(xpi2).lt.1d-6)xpi2=0d0
      c0=C0_(xpi1,q12,xpi2,m0,mi1,mi2,0)
c
c UV divergent B0 and B1 functions
c UV divergences have been subtracted: N_eps=2/eps-gamma_E+ln(4pi)
c
      call bfunc(q12,mi1,mi2,b01,b11,b20)
      call bfunc(xpi2,m0,mi2,b02,b12,b20)
      call bfunc(xpi1,m0,mi1,b03,b13,b20)
*---------------------------------------------------------------------
*--entwicklungskoeff. des vektorintegrals:----------------------------
      c11=1d0/2d0/det*(pi12*(b01-b03)+xpi2*(b02-b01)
     $     +(-xpi2*f1+pi12*f2)*c0)

      c12=1d0/2d0/det*(pi12*(b01-b02)+xpi1*(b03-b01)
     $     +(-xpi1*f2+pi12*f1)*c0)
*---------------------------------------------------------------------
*------------entwicklungskoeff. des tensorintegrals--------------
      c20 = (b01+1d0)/4d0+xm0/2d0*c0+(f1*c11+f2*c12)/4d0

      c21 = 1d0/2d0/det*(pi12*(-b01-b11-b13+f2*c11)+
     $     xpi2*(b01+b11-2d0*c20-f1*c11))

      c22 = 1d0/2d0/det*(pi12*(b11-b12+f1*c12)+
     $     xpi1*(-b11-2d0*c20-f2*c12))

      c23 = 1d0/2d0/det*(pi12*(b11+2d0*c20+f2*c12)+
     $     xpi2*(-b11+b12-f1*c12))
      return
      end

