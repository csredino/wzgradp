      REAL*8 FUNCTION ALPHAS2(SCALE)
C
C **********************************************************************
C *                                                                    *
C *   Program to calculate alfa strong to second order                 *
C *   as a function of lambda of 5 flavours                            *
C *   and the desired number of flavours (NFL)                         *
C *   for the selected set of structure functions which fixes lambda.  *
C *   NFL is limited to max. of 6 flavours.                            *
C *   The value of alfa is matched at the thresholds q = mq.           *
C *   When invoked with NFL < 0, it chooses NFL                        *
C *   as the number of flavours with mass less then q.                 *
C *   (mc = 1.5 GeV, mb = 4.75 GeV, mt = 100 GeV)                      *
C *                                                                    *
C *   Input:   SCALE   = QCD scale in GeV                              *
C *                                                                    *
C *   Output:  ALPHAS2 = alpha strong to second order,                 *
C *                      if not LO = 1                                 *
C *                                                                    *
C *   The variables NPTYPE,NGROUP,NSET,NFL,LO,TMAS,QCDL4,QCDL5         *
C *   should be provided by the user via a call to the                 *
C *   subroutine PDFSET at the initialization phase.                   *
C *                                                                    *
C *   Author:   H. Plothow-Besch (CERN-PPE)                            *
C *                                                                    *
C **********************************************************************
C
      implicit none
      real*8 NF,KNF
      real*8 PI,XMC,XMB,XMT,ZEROD,PONED,ONED,TWOD
      real*8 B6,BP6,B5,BP5,B4,BP4,B3,BP3,XLC,XLB,XLT,XLLC,XLLB,XLLT
      real*8 C65,C45,C35,Q,SCALE,XLQ,XLLQ,ALF
      real*8 QCDL4,QCDL5,TMAS
      common/W50512/QCDL4,QCDL5
      integer NFL,LO
      common/W50511/NFL,LO,TMAS

      DATA XMC/1.5D0/,XMB/4.75D0/,XMT/100.D0/
      DATA ZEROD/0.D0/,PONED/0.001D0/,ONED/1.D0/,TWOD/2.D0/
C
C    initialization of the parameters for the coupling constant
C
C    QCDL4  : QCD scale, which is matched to the scale used in the
C             parton distribution function chosen above for 4 flavours.
C
C    QCDL5  : QCD scale, which is matched to the scale used in the
C              parton distribution function chosen above.
C
C be aware, this lambda is given for 5 flavours (mb = 4.75 GeV/c**2) !!
C
C
C   set the default values
C
      ALPHAS2 = ZEROD

      PI=4.0D0*DATAN(ONED)
      B6  = (33.D0-2.D0*6.D0)/PI/12.D0
      BP6 = (153.D0 - 19.D0*6.D0) / PI / TWOD / (33.D0 - 2.D0*6.D0)
      B5  = (33.D0-2.D0*5.D0)/PI/12.D0
      BP5 = (153.D0 - 19.D0*5.D0) / PI / TWOD / (33.D0 - 2.D0*5.D0)
      B4  = (33.D0-2.D0*4.D0)/PI/12.D0
      BP4 = (153.D0 - 19.D0*4.D0) / PI / TWOD / (33.D0 - 2.D0*4.D0)
      B3  = (33.D0-2.D0*3.D0)/PI/12.D0
      BP3 = (153.D0 - 19.D0*3.D0) / PI / TWOD / (33.D0 - 2.D0*3.D0)
      XLC = TWOD * DLOG( XMC/QCDL5)
      XLB = TWOD * DLOG( XMB/QCDL5)
      XLT = TWOD * DLOG( XMT/QCDL5 * TMAS/XMT)
      XLLC = DLOG( XLC)
      XLLB = DLOG( XLB)
      XLLT = DLOG( XLT)
      C65  =  ONED/( ONED/(B5 * XLT) - XLLT*BP5/(B5 * XLT)**2 )
     +     - ONED/( ONED/(B6 * XLT) - XLLT*BP6/(B6 * XLT)**2 )
      C45  =  ONED/( ONED/(B5 * XLB) - XLLB*BP5/(B5 * XLB)**2 )
     +     - ONED/( ONED/(B4 * XLB) - XLLB*BP4/(B4 * XLB)**2 )
      C35  =  ONED/( ONED/(B4 * XLC) - XLLC*BP4/(B4 * XLC)**2 )
     +     - ONED/( ONED/(B3 * XLC) - XLLC*BP3/(B3 * XLC)**2 ) + C45
*     
      Q   = SCALE
      XLQ = TWOD *  DLOG( Q/QCDL5 )
      XLLQ =  DLOG( XLQ )
      KNF = NFL*1d0
      NF = KNF
      IF  ( NF .LT. ZEROD) THEN
         IF      ( Q .GT. XMT * TMAS/XMT) THEN
            NF = 6.D0
         ELSEIF  ( Q .GT. XMB ) THEN
            NF = 5.D0
         ELSEIF  ( Q .GT. XMC ) THEN
            NF = 4.D0
         ELSE
            NF = 3.D0
         ENDIF
      ENDIF
      IF(NF .GT. 6.D0) NF = 6.D0
      IF      ( NF .EQ. 6.D0 ) THEN
         ALF = ONED/(ONED/(ONED/(B6*XLQ)- BP6/(B6*XLQ)**2*XLLQ) + C65)
         IF (LO.EQ.1) ALF = ONED/B6/XLQ
      ELSEIF  ( NF .EQ. 5.D0 ) THEN
         ALF = ONED/(B5 * XLQ) -  BP5/(B5 * XLQ)**2 * XLLQ
         IF (LO.EQ.1) ALF = ONED/B5/XLQ
      ELSEIF  ( NF .EQ. 4.D0 ) THEN
         ALF = ONED/(ONED/(ONED/(B4*XLQ)- BP4/(B4*XLQ)**2*XLLQ) + C45)
         IF (LO.EQ.1) ALF = ONED/B4/XLQ
      ELSEIF  ( NF .EQ. 3.D0 ) THEN
         ALF = ONED/(ONED/(ONED/(B3*XLQ)- BP3/(B3*XLQ)**2*XLLQ) + C35)
         IF (LO.EQ.1) ALF = ONED/B3/XLQ
      ELSE
         write(6,*)'error in alphas2',NFL
         STOP
      ENDIF
      ALPHAS2 = ALF
      RETURN
      END
