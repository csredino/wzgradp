c
c********************* atlsmear **************************************
c
      subroutine atlsmear(p,sp,lep)
c*********************************************************************
c
c       this subroutne smears the 3-momentum component of the
c       4-vector p.
c
c       this subroutne is suitable for lhc (atlas)
c
c       input:
c      p   : 4-vector to be smeared (real*8)
c      lep : paricle type, 'e' for electron or photon, 'h' for hadron
c    common /iseed1/ iseed1 : integer seed for random number generator
c
c       output:
c       sp  : smeared 4-vector
c
c*********************************************************************
c
      implicit none
      character*1 lep
      integer i,j,iseed1
      real*8 scalet
      common/ scalet/ scalet
      common/iseed1/ iseed1
      real*8 p(4),sp(4)
      real*8 diff,sum,eta,arg,et,del,ret,dpx,dpy,dp
c
      diff = p(4)-p(3)
      sum = p(4)+p(3)
      if (dabs(diff).lt.1.d-5.or.dabs(sum).lt.1.d-5) then
         eta = 20.d0
      else
         arg = sum/diff
         if (arg.le.0.d0) then
            eta = 20.d0
         elseif (arg.gt.0.d0) then
            eta = dabs(0.5d0*dlog((p(4)+p(3))/(p(4)-p(3))))
         endif
      endif
      et = dsqrt(p(1)**2+p(2)**2)
      if (lep.eq.'e') then
         if (eta.lt.2.5d0) then
            del = dsqrt(0.095d0**2/p(4) + 0.005d0**2 ) * p(4)
         elseif (eta.ge.2.5d0) then
            del = 0.d0
         endif
      endif
      if (lep.eq.'m') then
         if (eta.lt.2.d0) then
            del = dsqrt((5.d-4*et)**2+0.012d0**2)*et*
     $           (1.d0+eta**10/7000.d0)
         elseif (eta.ge.2.d0) then
            del = 0.d0
         endif
      endif
      if (lep.eq.'n') then
         ret = scalet
         del = 0.45d0*dsqrt(ret)
         call rgauss(del,dpx)
         call rgauss(del,dpy)
         sp(1)=p(1)+dpx
         sp(2)=p(2)+dpy
         sp(3)=p(3)
         sp(4)=dsqrt(sp(1)**2+sp(2)**2+sp(3)**2)
         goto 15
      endif
      if (lep.eq.'h') then
         if (eta.lt.3.d0) then
            del = dsqrt(0.5d0**2/p(4) + 0.03d0**2 ) * p(4)
         elseif (eta.ge.3.d0.and.eta.lt.4.5d0) then
            del = dsqrt(1.d0**2/p(4) + 0.07d0**2 ) * p(4)
         elseif (eta.ge.4.5d0) then
            del = 0.d0
         endif
      endif
 9    call rgauss(del,dp)
      if (p(4).gt.0.d0.and.dp.gt.-p(4) ) then
         sp(4) = p(4) + dp
         do 10 i=1,3
            sp(i)=p(i)*(1.d0+dp/p(4))
 10      continue
      else
         sp(4) = 1d-15 * p(4)
         do 11 j=1,3
            sp(j)=1d-15*p(j)
 11      continue
      endif
 15   return
      end
c
c***************************************************************************
c
      subroutine rgauss(del,dp)
c-------yields gaussian distribution of random numbers
      implicit none
      integer iseed1
      common/iseed1/ iseed1
      data iseed1 /23147/
      real*4 x(2)
      real*8 del,dp,gauss

      if (del.le.0.d0) then
        dp = 0.d0
        return
      endif 

 5    call randa(2,x)
      
      dp=3d0*del*(-1.d0+2d0*x(1))

      gauss=dexp(-dp**2/2.d0/del**2)
      if (1d0*x(2)-gauss) 10,10,5
 10   return
      end
c
c*********************** nsmear  ******************************************
c
      subroutine nsmear(p,sp,lep)
c
c***************************************************************************
c
c       this subroutne does not smear
c
c***************************************************************************
*
c
      implicit none
      integer i
      character*1 lep
      real*8 p(4),sp(4)
      do 1 i=1,4
    1 sp(i) = p(i)
      return
      end

      subroutine metsmear(p,sp)
c
c***************************************************************************
c
c       this subroutine smears the pt-momentum component of the 4-vector p. of
c       a neutrino
c     
c       this subroutine is suitable for tevatron (d0) inspired)
c
c       input:
c       p   : 4-vector to be smeared (real*8)
c       sumet: sum of the et in the event
c
c       output:
c       sp  : smeared 4-vector
c
c***************************************************************************
c
      implicit none
      integer i
      real*8 p(4),sp(4),pi
      real*8 et,del,dp
      real*4 x
      parameter(pi=3.1415926535987932d0)
c
      et=sqrt(p(1)**2+p(2)**2)
      del=5.d0

      call rgauss(del,dp)
c
c     if et=0, then smear in random direction
c
      if(et.eq.0d0)then
 5       call randa(1,x)
         sp(1)=dp*dcos(x*pi)
         sp(2)=dp*dsin(x*pi)
      else
         do i=1,2
            sp(i)=p(i)*(1.d0+dp/et)
         end do
      end if
      sp(3)=p(3)
      sp(4)=dsqrt(sp(3)**2+sp(2)**2+sp(1)**2)
      return
      end
C
C*********************** D0UPGRSMEAR  *******************************************
C
      subroutine d0upgrsmear(p,sp,lep)
c
c***************************************************************************
c
c       this subroutine smears the 3-momentum component of the 4-vector p.
c
c       this subroutne is suitable for tevatron (upgraded d0)
c
c       input:
c       p   : 4-vector to be smeared (real*8)
c       lep : paricle type, 'e' for electron or photon, 'h' for hadron,
c             'n' for missing pt
c       common /iseed1/ iseed1 : integer seed for random number generator
c
c       output:
c       sp  : smeared 4-vector
c
c     updated 05/11/95
c
c****************************************************************************
c
      implicit none
      integer i,j,iseed1
      character*1 lep
      real*8 scalet
      common/ scalet/ scalet
      common/iseed1/ iseed1
      real*8 p(4),sp(4)
      real*8 ratio,eta,c,s,en,et,a,b,ret,del,dpx,dpy,dp
c
      if (p(4).ne.p(3)) then
        ratio = (p(4)+p(3))/(p(4)-p(3))
        if (ratio.gt.0.) then
          eta = dabs(0.5d0*dlog((p(4)+p(3))/(p(4)-p(3))))
        else
          eta = 10.d0
        endif
      else
        eta = 10.d0
      endif
      et = dsqrt(p(1)**2+p(2)**2)
c  electrons
      if (lep.eq.'e') then
        c = 0.003d0
        if (eta.lt.1.3d0) then
          s = 0.140d0
          en = 0.140d0
        elseif (eta.ge.1.3d0.and.eta.lt.4.d0) then
          s = 0.157d0
          en = 0.290d0
        endif
        if (eta.lt.4.d0) then
          del = dsqrt( (en/p(4))**2 + s**2/p(4) + c**2 )*p(4)
        elseif (eta.ge.4.d0) then
          del = 0.d0
        endif
        goto 9
      endif
c  muons
      if (lep.eq.'m') then
        a = 0.015d0
        b = 0.0016d0
        if (eta.lt.1.6d0) then
          del = dsqrt( (b*et)**2 + a**2 )*et
        else
          del = 0.d0
        endif
        goto 9
      endif
c jets
      if (lep.eq.'h') then
        del = dsqrt((2.16d0/et)**2 + 0.74d0**2/et + 0.01d0**2)*et
        goto 9
      endif
c missing et
      if (lep.eq.'n') then
        ret = scalet
        del = 1.08d0 + 0.019d0 * ret
        call rgauss(del,dpx)
        call rgauss(del,dpy)
        sp(1)=p(1)+dpx
        sp(2)=p(2)+dpy
        sp(3)=p(3)
        sp(4)=dsqrt(sp(1)**2+sp(2)**2+sp(3)**2)
        goto 15
      endif
    9 call rgauss(del,dp)
      if (p(4).gt.0.d0.and.dp.gt.-p(4) ) then
        sp(4) = p(4) + dp
        do 10 i=1,3
   10   sp(i)=p(i)*(1.d0+dp/p(4))
      else
        sp(4) = 1d-6 * p(4)
        do 11 j=1,3
   11   sp(j)=1d-6*p(j)
      endif
   15 return
      end
