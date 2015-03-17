************************************************************************
*                                                                      *
*    routines for the calculation of d0 for arbitrary real masses      *
*    no ir singularities are allowed                                   *
*    q13 < 0 required                                                  *
*                                                                      *
*    ansgar denner    18.4.96                                          *
************************************************************************
*   formulae are from:                                                 *
*     a. denner, u. nierste, r. scharf, nucl. phys. b367 (1991) 637    *
*   conventions:                                                       *
*     denominators: (q**2 - mm1**2)*[(q+k2)**2 - mm2**2]*              *
*                   [(q+k3)**2 - mm3**2]*[(q+k4)**2 - mm4**2]          *
*     invariants:   q12=k2**2, q23=(k2-k3)**2, q34=(k3-k4)**2,         *
*                   q14=k4**2, q24=(k2-k4)**2, q13=k3**2               *
************************************************************************
*                                                                      *
*   warning: d0m0,  d0m00, d0m000,d00000 go wrong for b*b = 4*a*c !!!  *
*                                                                      *
************************************************************************
* functions:                                                           *
* d0reg, d0m0,  d0m00, d0m000,d00000                                   *
* cdln,  cspenc,cspenh,cspen_new, cspcon,cspcoe, eta,ettile,etae       *
************************************************************************
      function d0reg(q12,q23,q34,q14,q24,q13,m1,m2,m3,m4)
************************************************************************
*  scalar 4-point function  for  q13 < 0 (one r real positive)         *
*  regular case                                                        *
*  imaginary part sometimes wrong                                      *
*  qij = (pi-pj)**2   and   pi is momentum of propagator with mass mi  *
*----------------------------------------------------------------------*
*  07.01.94 ansgar denner       last changed 23.10.95 ansgar denner    *
************************************************************************
      implicit   none
      real*8     q12,q23,q34,q14,q13,q24,m1,m2,m3,m4
      real*8     k12,k13,k14,k23,k24,k34
      real*8     m12,m22,m32,m42
      real*8     ir12,ir14,ir23,ir24,ir34
      real*8     ix(2,4),is(4)
      complex*16 r12,r13,r14,r23,r24,r34
      complex*16 a,b,c
      complex*16 x(2,4),s(4)
      complex*16 d0reg,cspcoe,etae,cdln,d0m0
      integer    k,j
 
      if (m1.eq.0d0) then
        d0reg = d0m0(q23,q12,q14,q34,q24,q13,m3,m2,m4)
      else if (m2.eq.0d0) then
        d0reg = d0m0(q14,q12,q23,q34,q13,q24,m4,m1,m3)
      else if (m3.eq.0d0) then
        d0reg = d0m0(q12,q23,q34,q14,q24,q13,m1,m2,m4)
      else if (m4.eq.0d0) then
        d0reg = d0m0(q23,q34,q14,q12,q13,q24,m2,m3,m1)
      else
      m12 = m1*m1
      m22 = m2*m2
      m32 = m3*m3
      m42 = m4*m4
      k12 = (m12+m22-q12)/m1/m2
      k13 = (m12+m32-q13)/m1/m3
      k14 = (m12+m42-q14)/m1/m4
      k23 = (m22+m32-q23)/m2/m3
      k24 = (m22+m42-q24)/m2/m4
      k34 = (m32+m42-q34)/m3/m4
      if (k13.lt.2d0) then
       write(*,*) ' d0reg: case not implemented'
       write(*,*) ' k13 = ',k13
       write(*,*) ' q13 = ',q13
       write(*,*) ' m1  = ',m1
       write(*,*) ' m3  = ',m3
      end if
      r12 = k12/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k12**2,0d0)))
      r13 = k13/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k13**2,0d0)))
      r13 = 1d0/r13
      r14 = k14/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k14**2,0d0)))
      r23 = k23/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k23**2,0d0)))
      r24 = k24/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k24**2,0d0)))
      r24 = 1d0/r24
      r34 = k34/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k34**2,0d0)))
      a   =  k34/r24-k23 + (k12-k14/r24)*r13
      b   =  (1d0/r13-r13)*(1d0/r24-r24)+k12*k34-k14*k23
      c   =  k34*r24-k23 + (k12-k14*r24)/r13
      x(1,4) = (-b+sqrt(b*b-4d0*a*c))/2d0/a
      x(2,4) = (-b-sqrt(b*b-4d0*a*c))/2d0/a
      if(abs(x(1,4)).gt.abs(x(2,4))) then
        x(2,4) = c/a/x(1,4)
      else
        x(1,4) = c/a/x(2,4)
      end if

      if(k12.lt.-2d0) then
        ir12 = sign(1d1,1d0-abs(r12))
      else
        ir12 = 0d0
      end if
      if(k14.lt.-2d0) then
        ir14 = sign(1d1,1d0-abs(r14))
      else
        ir14 = 0d0
      end if
      if(k23.lt.-2d0) then
        ir23 = sign(1d1,1d0-abs(r23))
      else
        ir23 = 0d0
      end if
      if(k24.lt.-2d0) then
        ir24 = sign(1d1,1d0-abs(r24))
      else
        ir24 = 0d0
      end if
      if(k34.lt.-2d0) then
        ir34 = sign(1d1,1d0-abs(r34))
      else
        ir34 = 0d0
      end if
      if(dreal(x(1,4)).gt.0d0) then
         ix(1,4) = 1d0
      else
         ix(1,4) = 0d0
      end if
      if(dreal(x(2,4)).gt.0d0) then
         ix(2,4) = -1d0
      else
         ix(2,4) = 0d0
      end if

      x(1,1) = x(1,4)/r24
      x(2,1) = x(2,4)/r24
      x(1,2) = x(1,4)/r24*r13
      x(2,2) = x(2,4)/r24*r13
      x(1,3) = x(1,4)*r13
      x(2,3) = x(2,4)*r13
      s(1)  = r12
      s(2)  = r23
      s(3)  = r34
      s(4)  = r14

      is(1)  = ir12
      is(2)  = ir23
      is(3)  = ir34
      is(4)  = ir14
c --> changed 23.10.95
      if(dreal(x(1,1)).gt.0d0) then
         ix(1,1) = ix(1,4) + ir24
      else
c        ix(1,1) = 0d0
         ix(1,1) = -ix(1,4) - ir24
      end if
      if(dreal(x(2,1)).gt.0d0) then
         ix(2,1) = ix(2,4) + ir24
      else
c        ix(2,1) = 0d0
         ix(2,1) = -ix(2,4) - ir24
      end if
c --> changed 23.10.95
      ix(1,3) = ix(1,4)
      ix(2,3) = ix(2,4)
      ix(1,2) = ix(1,1) 
      ix(2,2) = ix(2,1) 
 
      d0reg = dcmplx(0d0,0d0)
      do 20 k=1,2
         do 10 j=1,4
            d0reg = d0reg + (-1d0)**(j+k) * (
     &           cspcoe(-x(k,j),s(j),-ix(k,j),is(j))
     &           + cspcoe(-x(k,j),1d0/s(j),-ix(k,j),-is(j)) )
 10      continue
         d0reg = d0reg - (-1d0)**k* 
     &        etae(-x(k,4),1d0/r24,-ix(k,4),-ir24,-ix(k,1))*
     &        cdln((1d0+k14*x(k,4)+x(k,4)**2)/(1d0+k34*x(k,3)
     $        +x(k,3)**2),
     &        -dreal(1d0+k34*x(k,3)+x(k,3)**2))
 20   continue
      d0reg = d0reg/m1/m2/m3/m4/sqrt(b*b-4d0*a*c)
      end if
      end
************************************************************************
      function d0m0(q12,q23,q34,q14,q24,q13,m1,m2,m4)
************************************************************************
*  scalar 4-point function  for m3 = 0                                 *
*  regular case                                                        *
*----------------------------------------------------------------------*
*  29.03.92 ansgar denner       last changed 08.07.94 ansgar denner    *
************************************************************************
      implicit   none
      real*8     q12,q23,q34,q14,q13,q24,m1,m2,m4
      real*8     k12,k13,k14,k23,k24,k34
      real*8     m3,m12,m22,m32,m42
      real*8     ir12,ir14,ir24
      real*8     ix1(2),ix4(2)
      complex*16 r12,r13,r14,r24
      complex*16 a,b,c,d,det
      complex*16 x1(2),x4(2)
      complex*16 d0m0,cspcoe,ettile,d0m00,cdln
      integer    i
 
      if (m1.eq.0d0) then
        d0m0 = d0m00(q12,q13,q34,q24,q14,q23,m2,m4)
      else if (m2.eq.0d0) then
        d0m0 = d0m00(q12,q23,q34,q14,q24,q13,m1,m4)
      else if (m4.eq.0d0) then
        d0m0 = d0m00(q14,q34,q23,q12,q24,q13,m1,m2)
      else
      m3  = m1
      m12 = m1*m1
      m22 = m2*m2
      m32 = m3*m3
      m42 = m4*m4
      k12 = (m12+m22-q12)/m1/m2
      k13 = (m12    -q13)/m1/m3
      k14 = (m12+m42-q14)/m1/m4
      k23 = (m22    -q23)/m2/m3
      k24 = (m22+m42-q24)/m2/m4
      k34 = (    m42-q34)/m3/m4
      r12 = k12/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k12**2)))
      r13 = k13/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k13**2)))
      r14 = k14/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k14**2)))
      r24 = k24/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k24**2)))
      a   =  k34/r24-k23
      b   =  k13*(1d0/r24-r24)+k12*k34-k14*k23
      c   =  k13*(k12-r24*k14)+r24*k34-k23
      d   = -k34*r24+k23
      det =  k12*k12*k34*k34 + k14*k14*k23*k23 + k24*k24*k13*k13
     &     - 2d0*(k12*k23*k34*k14 + k12*k24*k34*k13 + k14*k24*k23*k13)
     &     + 4d0*(k12*k23*k13 + k23*k34*k24 + k14*k34*k13)
     &     - 4d0*(k23*k23 + k34*k34 + k13*k13)
      x4(1) = (-b+sqrt(det))/2d0/a 
      x4(2) = (-b-sqrt(det))/2d0/a 
      if(abs(x4(1)).gt.abs(x4(2))) then
        x4(2) = c/a/x4(1)
      else
        x4(1) = c/a/x4(2)
      end if
      x1(1) = x4(1)/r24
      x1(2) = x4(2)/r24

      if(k12.lt.-2d0) then
        ir12 = sign(1d1,1d0-abs(r12))
      else
        ir12 = 0d0
      end if
      if(k14.lt.-2d0) then
        ir14 = sign(1d1,1d0-abs(r14))
      else
        ir14 = 0d0
      end if
      if(k24.lt.-2d0) then
        ir24 = sign(1d1,1d0-abs(r24))
      else
        ir24 = 0d0
      end if
      ix4(1) = -sign(1d0,dreal(d)) 
      ix4(2) =  sign(1d0,dreal(d))
c  choice of impart avoids  i*pi*log(1-(1+eps)) terms:
c     ix4(2) = -sign(1d0,dreal(d))
c  but yields wrong imaginary part!  16.03.95
      ix1(1) =  sign(1d0,ix4(1)*dreal(r24))
      ix1(2) =  sign(1d0,ix4(2)*dreal(r24))
 
      d0m0 = dcmplx(0d0)
      do 10 i=1,2
      d0m0 = d0m0 + (2*i-3) * (
     &       cspcoe(-x4(i),r14,-ix4(i),ir14)
     &     + cspcoe(-x4(i),1d0/r14,-ix4(i),-ir14)
     &     - cspcoe(-x1(i),r12,-ix1(i),ir12) 
     &     - cspcoe(-x1(i),1d0/r12,-ix1(i),-ir12)
     &     - cspcoe(-x4(i),dcmplx(k34/k13),-ix4(i),-k13) 
     &     + cspcoe(-x1(i),dcmplx(k23/k13),-ix1(i),-k13)
     &     - ettile(-x4(i),1d0/r24,-ix4(i),-ir24) *
     &       (cdln((k12-r24*k14-(r24-1d0/r24)*x4(i))/d,
     &            dreal(-(r24-1d0/r24)*ix4(i))/d)
     &       +cdln(dcmplx(k13),-1d0)) ) 
10    continue
      d0m0 = d0m0/m1/m2/m3/m4/a/(x4(1)-x4(2))
      end if
      end
************************************************************************
      function d0m00(q12,q23,q34,q14,q24,q13,m1,m4)
************************************************************************
*  scalar 4-point function  for m2 = m3 = 0                            *
*  regular case                                                        *
*----------------------------------------------------------------------*
*  10.04.92 ansgar denner       last changed 13.07.95 ansgar denner    *
************************************************************************
      implicit   none
      real*8     q12,q23,q34,q14,q13,q24,m1,m2,m4
      real*8     k12,k13,k14,k23,k24,k34
      real*8     m3,m12,m22,m32,m42
      real*8     eps
      complex*16 r12,r13,r14,r23,r24,r34
      complex*16 k12e,k13e,r14e,k23e,k24e,k34e
      complex*16 a,b,c,d,cd
      complex*16 x4(2)
      complex*16 d0m00,cspcon,d0m000
      integer    i
c     common /eps/    eps
 
      eps = 1d-20
      if (m1.eq.0d0) then
        d0m00  = d0m000(q24,q34,q13,q12,q14,q23,m4)
      else if (m4.eq.0d0) then
        d0m00  = d0m000(q12,q23,q34,q14,q24,q13,m1)
      else
c     eps = 1d-15
      m3  = m1
      m2  = m1
      m12 = m1*m1
      m22 = m2*m2
      m32 = m3*m3
      m42 = m4*m4
      k12 = (m12    -q12)/m1/m2
      k13 = (m12    -q13)/m1/m3
      k14 = (m12+m42-q14)/m1/m4
      k23 = (       -q23)/m2/m3
      k24 = (    m42-q24)/m2/m4
      k34 = (    m42-q34)/m3/m4
c --> changed 13.07.94
      r12 = k12/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k12**2)))
      r13 = k13/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k13**2)))
      r14 = k14/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k14**2)))
      r23 = k23/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k23**2)))
      r24 = k24/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k24**2)))
      r34 = k34/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k34**2)))
c --> changed 13.07.94
      a   =  k34*k24-k23
      b   =  k13*k24+k12*k34-k14*k23
      c   =  k13*k12-k23
      d   =  k23
      cd  = c+dcmplx(0d0,eps)*d
      x4(1) = (-b+sqrt(b*b-4d0*a*cd))/2d0/a
      x4(2) = (-b-sqrt(b*b-4d0*a*cd))/2d0/a
      if(abs(x4(1)).gt.abs(x4(2))) then
        x4(2) = cd/a/x4(1)
      else
        x4(1) = cd/a/x4(2)
      end if
      k12e  = k12*dcmplx(1d0,-sign(eps,k12))
      k13e  = k13*dcmplx(1d0,-sign(eps,k13))
      k23e  = k23*dcmplx(1d0,-sign(eps,k23))
      k24e  = k24*dcmplx(1d0,-sign(eps,k24))
      k34e  = k34*dcmplx(1d0,-sign(eps,k34))
      r14e  = r14*dcmplx(1d0,sign(eps,dreal(1d0/r14-r14)))
 
      d0m00 = dcmplx(0d0)
      do 10 i=1,2
      d0m00 = d0m00 + (2*i-3) * (
     &       cspcon(-x4(i),r14e) + cspcon(-x4(i),1d0/r14e)
     &     - cspcon(-x4(i),k34e/k13e) - cspcon(-x4(i),k24e/k12e)
     &     + log(-x4(i))*(log(k12e)+log(k13e)-log(k23e))  )
10    continue
      d0m00 = d0m00/m1/m2/m3/m4/a/(x4(1)-x4(2))
      end if
      end
************************************************************************
      function d0m000(q12,q23,q34,q14,q24,q13,m1)
************************************************************************
*  scalar 4-point function  for m2 = m3 = m4 = 0                       *
*  regular case                                                        *
*----------------------------------------------------------------------*
*  23.04.92 ansgar denner       last changed 17.03.95 ansgar denner    *
************************************************************************
      implicit   none
      real*8     q12,q23,q34,q14,q13,q24,m1
      real*8     k12,k13,k14,k23,k24,k34
      real*8     m2,m3,m4,m12,m22,m32,m42
      real*8     eps
      complex*16 k12e,k13e,k14e,k23e,k24e,k34e
      complex*16 a,b,c,d,cd
      complex*16 x4(2)
      complex*16 d0m000,d00000,cspcon
      integer    i
c     common /eps/    eps
 
      eps = 1d-20
      if (m1.eq.0d0) then
        d0m000 = d00000(q12,q23,q34,q14,q24,q13)
      else
c     eps = 1d-15
      m4  = m1
      m3  = m1
      m2  = m1
      m12 = m1*m1
      m22 = m2*m2
      m32 = m3*m3
      m42 = m4*m4
      k12 = (m12    -q12)/m1/m2
      k13 = (m12    -q13)/m1/m3
      k14 = (m12    -q14)/m1/m4
      k23 = (       -q23)/m2/m3
      k24 = (       -q24)/m2/m4
      k34 = (       -q34)/m3/m4
      a   =  k34*k24
      b   =  k13*k24+k12*k34-k14*k23
      c   =  k13*k12-k23
      d   =  k23
      cd  = c+dcmplx(0d0,eps)*d
      x4(1) = (-b+sqrt(b*b-4d0*a*cd))/2d0/a
      x4(2) = (-b-sqrt(b*b-4d0*a*cd))/2d0/a
      if(abs(x4(1)).gt.abs(x4(2))) then
        x4(2) = cd/a/x4(1)
      else
        x4(1) = cd/a/x4(2)
      end if
      k12e  = k12*dcmplx(1d0,-sign(eps,k12))
      k13e  = k13*dcmplx(1d0,-sign(eps,k13))
      k23e  = k23*dcmplx(1d0,-sign(eps,k23))
      k24e  = k24*dcmplx(1d0,-sign(eps,k24))
      k34e  = k34*dcmplx(1d0,-sign(eps,k34))
      k14e  = k14*dcmplx(1d0,-sign(eps,k14))
 
      d0m000 = dcmplx(0d0)
      do 10 i=1,2
      d0m000 = d0m000 + (2*i-3) * (
     &       cspcon(-x4(i),k14e)
     &     - cspcon(-x4(i),k34e/k13e) - cspcon(-x4(i),k24e/k12e)
     &     + log(-x4(i))*(log(k12e)+log(k13e)-log(k23e))  )
10    continue
      d0m000 = d0m000/m1/m2/m3/m4/a/(x4(1)-x4(2))
      end if
      end
************************************************************************
      function d00000(q12,q23,q34,q14,q24,q13)
************************************************************************
*  scalar 4-point function  for m1 = m2 = m3 = m4 = 0                  *
*  regular case                                                        *
*----------------------------------------------------------------------*
*  20.01.95 ansgar denner       last changed 16.05.93 ansgar denner    *
************************************************************************
      implicit   none
      real*8     q12,q23,q34,q14,q13,q24
      real*8     k12,k13,k14,k23,k24,k34
      real*8     m2
      real*8     eps
      complex*16 k12e,k13e,k14e,k23e,k24e,k34e
      complex*16 a,b,c,d,cd
      complex*16 x4(2)
      complex*16 d00000,cspcon
      integer    i
c     common /eps/    eps
 
      eps = 1d-20
      m2  = abs(q24)
      k12 = (       -q12)/m2
      k13 = (       -q13)/m2
      k14 = (       -q14)/m2
      k23 = (       -q23)/m2
      k24 = (       -q24)/m2
      k34 = (       -q34)/m2
      a   =  k34*k24
      b   =  k13*k24+k12*k34-k14*k23
      c   =  k13*k12
      d   =  k23
      cd  = c+dcmplx(0d0,eps)*d
      x4(1) = (-b+sqrt(b*b-4d0*a*cd))/2d0/a
      x4(2) = (-b-sqrt(b*b-4d0*a*cd))/2d0/a
      if(abs(x4(1)).gt.abs(x4(2))) then
        x4(2) = cd/a/x4(1)
      else
        x4(1) = cd/a/x4(2)
      end if
      k12e  = k12*dcmplx(1d0,-sign(eps,k12))
      k13e  = k13*dcmplx(1d0,-sign(eps,k13))
      k23e  = k23*dcmplx(1d0,-sign(eps,k23))
      k24e  = k24*dcmplx(1d0,-sign(eps,k24))
      k34e  = k34*dcmplx(1d0,-sign(eps,k34))
      k14e  = k14*dcmplx(1d0,-sign(eps,k14))
 
      d00000 = dcmplx(0d0)
      do 10 i=1,2
      d00000 = d00000 + (2*i-3) * (
     &     - log(-x4(i))**2/2d0
     &     - cspcon(-x4(i),k34e/k13e) - cspcon(-x4(i),k24e/k12e)
     &     + log(-x4(i))*(log(k12e)+log(k13e)-log(k14e)-log(k23e))  )
10    continue
      d00000 = d00000/m2/m2/a/(x4(1)-x4(2))
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function spen(x)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       spence function
c       realteil
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       02.10.89 ansgar denner
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        implicit logical(a-z)
        complex*16 cspen_new,cx
        real*8     spen,x
 
        cx=dcmplx(x)
        spen=dreal(cspen_new(cx))
        end
************************************************************************
      function d0m0ne(q12,q23,q34,q14,q24,q13,m1,m2,m4)
************************************************************************
*  scalar 4-point function  for m3 = 0                                 *
*  regular case                                                        *
*----------------------------------------------------------------------*
*  29.03.92 ansgar denner
************************************************************************
      implicit logical (a-z)
      real*8     q12,q23,q34,q14,q13,q24,m1,m2,m4
      real*8     k12,k13,k14,k23,k24,k34
      real*8     m3,m12,m22,m32,m42
      real*8     eps
      complex*16 r12,r13,r14,r23,r24,r34
      complex*16 r12e,k13e,r14e,k23e,r24e,k34e
      complex*16 a,b,c,d,cd
      complex*16 x1(2),x4(2)
      complex*16 d0m0ne,cspcon,ettile
      integer    i
 
      eps = 1d-15
      m3  = m1
      m12 = m1*m1
      m22 = m2*m2
      m32 = m3*m3
      m42 = m4*m4
      k12 = (m12+m22-q12)/m1/m2
      k13 = (m12    -q13)/m1/m3
      k14 = (m12+m42-q14)/m1/m4
      k23 = (m22    -q23)/m2/m3
      k24 = (m22+m42-q24)/m2/m4
      k34 = (    m42-q34)/m3/m4
      r12 = k12/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k12**2)))
      r13 = k13/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k13**2)))
      r14 = k14/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k14**2)))
      r23 = k23/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k23**2)))
      r24 = k24/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k24**2)))
      r34 = k34/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k34**2)))
      a   =  k34/r24-k23
      b   =  k13*(1d0/r24-r24)+k12*k34-k14*k23
      c   =  k13*(k12-r24*k14)+r24*k34-k23
      d   = -k34*r24+k23
      cd  = c+dcmplx(0d0,eps)*d
      x4(1) = (-b+sqrt(b*b-4d0*a*cd))/2d0/a
      x4(2) = (-b-sqrt(b*b-4d0*a*cd))/2d0/a
      if(abs(x4(1)).gt.abs(x4(2))) then
        x4(2) = cd/a/x4(1)
      else
        x4(1) = cd/a/x4(2)
      end if
      x1(1) = x4(1)/r24
      x1(2) = x4(2)/r24
      k13e  = k13*dcmplx(1d0,-sign(eps,k13))
      k23e  = k23*dcmplx(1d0,-sign(eps,k23))
      k34e  = k34*dcmplx(1d0,-sign(eps,k34))
      r12e  = r12*dcmplx(1d0,sign(eps,dreal(1d0/r12-r12)))
      r14e  = r14*dcmplx(1d0,sign(eps,dreal(1d0/r14-r14)))
      r24e  = r24*dcmplx(1d0,sign(eps,dreal(1d0/r24-r24)))
 
      d0m0ne = dcmplx(0d0)
      do 10 i=1,2
      d0m0ne = d0m0ne + (2*i-3) * (
     &       cspcon(-x4(i),r14e) + cspcon(-x4(i),1d0/r14e)
     &     - cspcon(-x1(i),r12e) - cspcon(-x1(i),1d0/r12e)
     &     - cspcon(-x4(i),k34e/k13e) + cspcon(-x1(i),k23e/k13e)
     &     - ettile(-x4(i),1d0/r24,1d0/r24e) *
     &       (log((k12-r24*k14-(r24-1d0/r24)*x4(i))/d)+log(k13e))  )
10    continue
      d0m0ne = d0m0ne/m1/m2/m3/m4/a/(x4(1)-x4(2))
      end
************************************************************************
      function cdln(cz,eps)
************************************************************************
*     complex logarithm of cz + i*eps                                  *
*----------------------------------------------------------------------*
*     09.01.90 ansgar denner                                           *
************************************************************************
      implicit   none
      real*8     eps
      complex*16 cdln,cz,cdlog
      real*8     pi

      pi = 4d0*datan(1d0)
      if(dimag(cz).eq.0d0.and.dreal(cz).le.0.d0)then
        if (eps.eq.0d0) then
          write(80,*) 'cdln:  argument on cut '
          write(80,*) 'cdln:  eps = 0'
          write(80,*) 'cdln:  cz  = ',cz
        end if
        cdln=cdlog(-cz)+dcmplx(0d0,pi)*dsign(1d0,eps)
      else
        cdln=cdlog(cz)
      end if
      end
************************************************************************
      function cspenc(cz,eps)
************************************************************************
*       complex spence function  of cz + i*eps                         *
*       calculated by mapping on the area where there is a quickly     *
*       convergent series                                              *
*----------------------------------------------------------------------*
*     08.01.90 ansgar denner                                           *
************************************************************************
      implicit   none
      real*8     eps
      complex*16 cspenc,cz
      real*8     pi,pi26
      real*8     az,rz,az1
      complex*16 cz1,cspenh,cdln

      pi = 4d0*datan(1d0)
      pi26 = pi*pi/6d0
      cz1     = 1d0-cz
      az1     = cdabs(cz1)
      az      = cdabs(cz)
      rz      = dreal(cz)

      if (az1.lt.1d-15) then
         cspenc = pi26
      else if (rz.lt.0.5d0) then
         if (az.lt.1d0) then
            cspenc = cspenh(cz,eps)
         else
            cspenc = -pi26 - .5d0*cdln(-cz,-eps)**2
     &           - cspenh(1d0/cz,-eps)
         end if
      else
         if (az1.lt.1d0) then
            cspenc =  pi26 - cdln(cz,eps)*cdln(cz1,-eps)
     &           - cspenh(cz1,-eps)
         else
            cspenc = 2d0*pi26 + .5d0*cdln(-cz1,-eps)**2
     &           - cdln(cz,eps)*cdln(cz1,-eps)
     &           + cspenh(1d0/cz1,eps)
         end if
      end if
      end
************************************************************************
      function cspenh(cz,eps)
************************************************************************
*       complex spence function of cz + i*eps                          *
*       in convergence region                                          *
*       calculation of bernoulli series                                *
*----------------------------------------------------------------------*
*     09.01.90 ansgar denner                                           *
************************************************************************
      implicit   none
      complex*16 cspenh,cdln,cz,x,x2
      real*8     eps
c     real*8     b(11),eps
      real*8 b(11)/
     1   0.1666666666666666666666666667d0,
     2  -0.0333333333333333333333333333d0,
     3   0.0238095238095238095238095238d0,
     4  -0.0333333333333333333333333333d0,
     5   0.0757575757575757575757575758d0,
     6  -0.2531135531135531135531135531d0,
     7   1.1666666666666666666666666667d0,
     8  -7.0921568627450980392156862745d0,
     9  54.97117794486215538847117794486d0,
     +  -529.124242424242424242424242424242d0,
     1  6192.123188405797101449275362318d0  /
c     beachte:                 b(n)=b2n
c     b(1)=1./6.
c     b(2)=-1./30.
c     b(3)=1./42.
c     b(4)=-1./30.
c     b(5)=5./66.
c     b(6)=-691./2730.
c     b(7)=7./6.
c     b(8)=-3617./510.
c     b(9)=43867./798.
c     b(10)=-174611./330.
c     b(11)=854513./138.
c     pi=3.1415926535897932384
c     pi*pi/6.=1.6449..., pi*pi/3=3.28986...
c
      integer    j
      real*8     factor
      complex*16 power,term,csp

      b(11)  =    854513d0/ 138d0
      b(10)  =  - 174611d0/ 330d0
      b(9)   =     43867d0/ 798d0
      b(8)   =  -   3617d0/ 510d0
      b(7)   =         7d0/   6d0
      b(6)   =  -    691d0/2730d0
      b(5)   =         5d0/  66d0
      b(4)   =  -      1d0/  30d0
      b(3)   =         1d0/  42d0
      b(2)   =  -      1d0/  30d0
      b(1)   =         1d0/   6d0
      x      =  -cdln(1d0-cz,-eps)
c     write(80,*)  'cspenh'
      x2     =  x*x
      power  =  x
      factor =  1d0
      cspenh =  x - x2/4d0
      do 10 j=2,22,2
         factor = factor / j / (j+1)
         power  = power * x2
         term   = b(j/2) * factor * power
         csp    = cspenh + term
         if (csp.eq.cspenh) return
         cspenh = csp
10    continue
      if (cdabs(term/csp).gt.1d-15) then
        write(80,*) 'cspenh converges badly  ',cz,x
        write(80,*) 'cspenh converges badly  ',csp-term,cspenh,term
      end if 
      end
************************************************************************
      function cspen_new(cz)
************************************************************************
*       complex spence function                                        *
*----------------------------------------------------------------------*
*     08.07.94 ansgar denner                                           *
************************************************************************
      implicit   none
      complex*16 cspen_new,cspenc,cz

      if((dimag(cz).eq.0d0).and.(dreal(cz-1d0).gt.5d-13))then
        write(80,*) 'cspen_new:  argument on cut '
        write(80,*) 'cspen_new:  cz  = ',cz
      end if
      cspen_new = cspenc(cz,0d0)
      end
************************************************************************
      function cspcon(z1,z2)
************************************************************************
*  complex spence function plus continuation terms                     *
*----------------------------------------------------------------------*
*  29.03.92 ansgar denner                                              *
************************************************************************
      implicit   none
      complex*16 cspcon,cspen_new,eta,z1,z2
      real*8     pi,pi26
 
      pi = 4d0*datan(1d0)
      pi26 = pi*pi/6d0
      if(dreal(z1*z2).gt.0d0) then
        cspcon = cspen_new(1d0-z1*z2) + eta(z1,z2)*log(1d0-z1*z2)
      else
        cspcon = pi26-cspen_new(z1*z2)
     &           - (log(z1)+log(z2))*log(1d0-z1*z2)
      end if
      end
************************************************************************
      function cspcoe(z1,z2,i1,i2)
************************************************************************

*  complex spence function plus continuation terms                     *
*  i@ is assumed to dominate i1                                        *
*----------------------------------------------------------------------*
*  08.07.94 ansgar denner      last changed  17.03.95                  *
************************************************************************
      implicit   none
      complex*16 cspcoe,cspenc,etae,cdln,z1,z2,z12,etaa
      real*8     i1,i2,i12
      real*8     pi,pi26

      pi = 4d0*datan(1d0)
      pi26 = pi*pi/6d0
      z12 = z1*z2
      i12 = i2*sign(1d0,dreal(z1))
      if(dreal(z12).gt.0.5d0) then
        cspcoe = cspenc(1d0-z12,0d0) 
        etaa   = etae(z1,z2,i1,i2,i12)
        if(etaa.ne.0d0) then
          cspcoe = cspcoe + etaa*cdln(1d0-z12,-i12)
        end if
      else if(abs(z12).lt.1d-4) then
         cspcoe = pi26-cspenc(z12,0d0)
     &        + (cdln(z1,i1)+cdln(z2,i2))
     &        *z12*(1d0+z12/2d0+z12*z12/3d0+z12*z12*z12/4d0)
      else
         cspcoe = pi26-cspenc(z12,0d0)
     &        - (cdln(z1,i1)+cdln(z2,i2))
     &        *cdln(1d0-z12,-0d0)
      end if

      end
************************************************************************
      function eta(c1,c2)
************************************************************************
*     complex eta-function
*----------------------------------------------------------------------*
*     8.06.90    ansgar denner       last changed   11.07.94
************************************************************************
      implicit     none
      complex*16 eta,c1,c2
      real*8     im1,im2,im12,re1,re2
      real*8     pi

      pi = 4d0*datan(1d0)
      re1    = dreal(c1)
      re2    = dreal(c2)
      im1    = dimag(c1)
      im2    = dimag(c2)
      im12   = dimag(c1*c2)
 
      if(im1.lt.0d0.and.im2.lt.0d0.and.im12.gt.0d0) then
          eta = dcmplx(0d0,2d0*pi)
      else if (im1.gt.0d0.and.im2.gt.0d0.and.im12.lt.0d0) then
          eta = dcmplx(0d0,-2d0*pi)
      else
          eta = dcmplx(0d0,0d0)
          if(.not.(im2.eq.0d0.and.re2.gt.0d0.or.
     &             im1.eq.0d0.and.re1.gt.0d0).and.
     &       (im1.eq.0.and.re1.lt.0d0 .or.
     &        im2.eq.0.and.re2.lt.0d0 .or.
     &        im12.eq.0.and.dreal(c1*c2).lt.0d0)) then
             write(80,*) ' eta not defined '
             write(80,*) ' eta:  c1  = ',c1
             write(80,*) ' eta:  c2  = ',c2
             write(80,*) ' eta:  c12 = ',c1*c2
          end if
      end if
      end
************************************************************************
      function etae(c1,c2,i1,i2,i12)
************************************************************************
*     complex eta-function
*----------------------------------------------------------------------*
*     8.06.90    ansgar denner       last changed   11.07.94
************************************************************************
      implicit     none
      complex*16 etae,c1,c2
      real*8     i1,i2,i12
      real*8     im1,im2,im12
      real*8     pi

      pi = 4d0*datan(1d0)
      im1    = dimag(c1)
      im2    = dimag(c2)
      im12   = dimag(c1*c2)
      if(im1 .eq.0d0) im1  = i1
      if(im2 .eq.0d0) im2  = i2
      if(im12.eq.0d0) im12 = i12
 
      if(im1.lt.0d0.and.im2.lt.0d0.and.im12.gt.0d0) then
         etae = dcmplx(0d0,2d0*pi)
      else if (im1.gt.0d0.and.im2.gt.0d0.and.im12.lt.0d0) then
         etae = dcmplx(0d0,-2d0*pi)
      else
         etae = dcmplx(0d0,0d0)
         if(im1.eq.0.and.dreal(c1).lt.0d0 .or.
     &        im2.eq.0.and.dreal(c2).lt.0d0 .or.
     &        im12.eq.0.and.dreal(c1*c2).lt.0d0) then
            write(80,*) ' eta not defined '
            write(80,*) ' eta:  c1  = ',c1
            write(80,*) ' eta:  c2  = ',c2
            write(80,*) ' eta:  c12 = ',c1*c2
         end if
      end if
      end
************************************************************************
      function ettile(c1,c2,i1,i2)
************************************************************************
*  complex eta-tilde   for 16-spence-d0                                *
*----------------------------------------------------------------------*
*  29.03.92 ansgar denner        last changed 14.07.94                 *
************************************************************************
      implicit     none
      complex*16 ettile,etae,c1,c2
      real*8     i1,i2
      real*8     im1,im2,re2
      real*8     pi

      pi = 4d0*datan(1d0)
      im1    = dimag(c1)
      if(im1.eq.0d0) im1 = i1
      im2    = dimag(c2)
      re2    = dreal(c2)
      if(im2.ne.0d0) then
         ettile = etae(c1,c2,i1,0d0,0d0)
      else if (re2.gt.0d0) then
         ettile = dcmplx(0d0,0d0)
      else if (im1.gt.0d0.and.i2.gt.0d0) then
         ettile = dcmplx(0d0,-2d0*pi)
      else if (im1.lt.0d0.and.i2.lt.0d0) then
         ettile = dcmplx(0d0, 2d0*pi)
      else
         ettile = dcmplx(0d0,0d0)
         if(im1.eq.0.and.dreal(c1).lt.0d0 .or.
     &        im2.eq.0.and.i2.eq.0d0.and.dreal(c1*c2).lt.0d0) then
            write(80,*) ' ettile not defined '
            write(80,*) ' ettile:  c1  = ',c1
            write(80,*) ' ettile:  c2  = ',c2
            write(80,*) ' ettile:  i1  = ',i1
            write(80,*) ' ettile:  i2  = ',i2
         end if
      end if
      end
