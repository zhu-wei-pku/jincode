      SUBROUTINE REVERSE_T1
     +     (ainp, nparmin, nparmout, b, q, ra, dec, qlat, qlong,
     +     t0par, t1par, t2par, tbinary, te_0, te)
c     Subroutine to find te given t1
      implicit none
      integer nparmout, nparmin
      real*8 ainp(nparmin)
      real*8 aoutp(nparmout)
      real*8 t2, u0, t1, rhos, piex, piey, theta, db, dep, b, q !input
      real*8 te, t0 !output

      real*8 t0par, t1par, t2par, tbinary, ra, dec, qlat, qlong

      real*8 te_lo, te_hi, te_try, te_0
      real*8 t1_lo, t1_hi, t1_try !trial t2s

      integer i, n
      real*8 tol
      

      tol = 10.d0**(-6) !required precision
      n = 11 !number of loops

c     Initialize variables
      t2    = ainp(1)
      u0    = ainp(2)
      t1    = ainp(3)
      rhos  = ainp(4)
      piex  = ainp(5)
      piey  = ainp(6)
      theta = ainp(7)
      db    = ainp(8)
      dep   = ainp(9)

      if(te_0.lt.1d0)te_0=35.0

      aoutp(1) = t2
      aoutp(2) = u0

      aoutp(4) = rhos
      aoutp(5) = piex 
      aoutp(6) = piey
      aoutp(7) = theta
      aoutp(8) = db
      aoutp(9) = dep
      aoutp(10)= b
      aoutp(11)= q

      te_lo = te_0 - 0.1d0
      te_hi = te_0 + 0.1d0

 100  continue
      aoutp(3) = te_lo
      call caustic1
     +     (aoutp, nparmin, nparmout, ra, dec, qlat, qlong, 
     +     t0par, t1par, t2par, tbinary, t1_lo)
c      write(6,*) 't1',t1, t1_lo

      if (t1_lo.lt.t1) then
         te_lo = te_lo-0.1d0
         te_hi = te_hi-0.1d0
         go to 100
      endif


 200  continue
      aoutp(3) = te_hi
      call caustic1
     +     (aoutp, nparmin, nparmout, ra, dec, qlat, qlong, 
     +     t0par, t1par, t2par, tbinary, t1_hi)

      if(t1_hi.gt.t1) then
         te_hi = te_hi+0.1d0
         goto 200
      endif

c      n = int(dlog((te_hi-te_lo)/tol)/dlog(2d0))
c     if (n.lt.15) n=15

      do i = 1, n 
         te_try = (te_lo+te_hi)/2d0

         aoutp(3) = te_try
         call caustic1
     +        (aoutp, nparmin, nparmout, ra, dec, qlat, qlong, 
     +        t0par, t1par, t2par, tbinary, t1_try)

         if (t1_try.lt.t1) then
            te_hi = te_try
         else
            te_lo  = te_try
         endif
         

      enddo

      te = (te_lo+te_hi)/2d0
      END


      SUBROUTINE CAUSTIC1
     +     (aoutp, nparmin, nparmout, ra, dec, qlat, qlong, 
     +     t0par, t1par, t2par, tbinary, t1)
c     Given u0, calculate the time of the first caustic crossing
      implicit none
      integer nparmout, nparmin
      real*8 aoutp(nparmout)
      real*8 ainp(nparmin)

      real*8 t2, u0, te, rhos, piex, piey, theta, db, dep, b_0, q !input
      real*8 b, ep, ep_0, tbinary, dcm
      real*8 t1 !Output
      real*8 t0 !internal

      real*8 t1_lo, t1_hi !Trial t1s.
      integer n_lo, n_hi !number of images for trial t1s.
      real*8 t1_try
      integer n_try
      integer i,n

      real*8 offset_x, offset_y 
c     offset of center of mass from b/2
      real*8 tau, xs, ys, xcm, ycm

      real*8 magnew
      integer max_sol
      parameter (max_sol=5)
      real*8 am(max_sol), xic(max_sol), yic(max_sol)!other getbinp parameters
      real*8 qn, qe, qnp, qep, taup, betap
      real*8 dtau, dbeta

      real*8 ra, dec, qlat, qlong, t0par, t1par, t2par
      real*8 tol

c-----Fixed paramters
      tol = 10.d0**(-6) !required precision
      n = 11 !number of loops

      ep_0=0d0 !reference binary rotation
c      tbinary = 4716.7d0 !reference time (where the map is calculated)
c---------------------

      t2    = aoutp(1)
      u0    = aoutp(2)
      te    = aoutp(3)
      rhos  = aoutp(4)
      piex  = aoutp(5)
      piey  = aoutp(6)
      theta = aoutp(7)
      db    = aoutp(8)
      dep   = aoutp(9)
      b_0   = aoutp(10)
      q     = aoutp(11)

      ainp(1) = t2
      ainp(2) = u0
      ainp(3) = te
      ainp(4) = rhos
      ainp(5) = piex
      ainp(6) = piey
      ainp(7) = theta
      ainp(8) = db
      ainp(9) = dep

c      write(6,*) 'caustic 1: t2, te: ', t2, te

c     Find t0
      call reverse_t2
     +     (ainp, nparmin, nparmout, b_0, q, ra, dec, qlat, qlong,
     +      t2par, t0par, tbinary, t0)
      
      t1 = 0d0

      t1_lo = -0.1d0
      t1_hi = 0.2d0

c     Position of the center of mass
      offset_y = 0.d0

c     Check that the boundaries are good.
 100  continue
         tau = ((t1_lo-t0)+(t1par-t0par))/te

         if (t1_lo+t1par.lt.4660.)then
            b = b_0+db*(4660.-tbinary)/365.25
            ep = ep_0 + dep*(4660.-tbinary)/365.25
         else if (t1_lo+t1par.gt.4800.)then
            b = b_0+db*(4800.-tbinary)/365.25
            ep = ep_0 + dep*(4800.-tbinary)/365.25
         else
            b = b_0+db*(t1_lo+t1par-tbinary)/365.25
            ep = ep_0 + dep*(t1_lo+t1par-tbinary)/365.25
         endif

         offset_x = b*(-0.5d0+1d0/(1d0+q))

         call geta(qn, qe, t1_lo+t1par, ra, dec, t0par)
         call gett(qnp, qep, t1_lo+t1par, qlat, qlong,
     +        ra, dec)

         qn = qn + qnp
         qe = qe + qep
         dtau = piex*qn + piey*qe
         dbeta = -piex*qe + piey*qn
         dbeta = -dbeta
         taup  = tau + dtau
         betap = u0 + dbeta

         xcm= taup*dcos(theta)+betap*dsin(theta)
         ycm=-taup*dsin(theta)+betap*dcos(theta)

         xs = xcm*dcos(ep)-ycm*dsin(ep)-offset_x
         ys = xcm*dsin(ep)+ycm*dcos(ep)

         call getbinp(magnew,am,xic,yic,xs,ys,q,b,n_lo)

         if (n_lo.ne.3) then
c     Shift the interval so that t1_lo is correct
            t1_lo = t1_lo-0.1d0
            t1_hi = t1_hi-0.1d0
            go to 100
         endif

 200     continue
         tau = ((t1_hi-t0)+(t1par-t0par))/te

         if (t1_hi+t1par.lt.4660.)then
            b = b_0+db*(4660.-tbinary)/365.25
            ep = ep_0 + dep*(4660.-tbinary)/365.25
         else if (t1_hi+t1par.gt.4800.)then
            b = b_0+db*(4800.-tbinary)/365.25
            ep = ep_0 + dep*(4800.-tbinary)/365.25
         else
            b = b_0+db*(t1_hi+t1par-tbinary)/365.25
            ep = ep_0 + dep*(t1_hi+t1par-tbinary)/365.25
         endif
         
         offset_x = b*(-0.5d0+1d0/(1d0+q))

         call geta(qn, qe, t1_hi+t1par, ra, dec, t0par)
         call gett(qnp, qep, t1_hi+t1par, qlat, qlong,
     +        ra, dec)

         qn = qn + qnp
         qe = qe + qep
         dtau = piex*qn + piey*qe
         dbeta = -piex*qe + piey*qn
         dbeta = -dbeta
         taup  = tau + dtau
         betap = u0 + dbeta

         xcm= taup*dcos(theta)+betap*dsin(theta)
         ycm=-taup*dsin(theta)+betap*dcos(theta)

         xs = xcm*dcos(ep)-ycm*dsin(ep)-offset_x
         ys = xcm*dsin(ep)+ycm*dcos(ep)

         call getbinp(magnew,am,xic,yic,xs,ys,q,b,n_hi)

         if (n_hi.ne.5) then 
c     Expand the interval so that t1_hi is correct
            t1_hi = t1_hi+0.1d0
            go to 200
         endif

c      n = int(dlog((t1_hi-t1_lo)/tol)/dlog(2d0))
c      if (n.lt.15) n=15

c     Use a binary search to find t1 given t0.
      do i = 1, n
         t1_try = (t1_lo+t1_hi)/2.d0
         tau = ((t1_try-t0)+(t1par-t0par))/te

         if (t1_try+t1par.lt.4660.)then
            b = b_0+db*(4660.-tbinary)/365.25
            ep = ep_0 + dep*(4660.-tbinary)/365.25
         else if (t1_try+t1par.gt.4800.)then
            b = b_0+db*(4800.-tbinary)/365.25
            ep = ep_0 + dep*(4800.-tbinary)/365.25
         else
            b = b_0+db*(t1_try+t1par-tbinary)/365.25
            ep = ep_0 + dep*(t1_try+t1par-tbinary)/365.25
         endif

         offset_x = b*(-0.5d0+1d0/(1d0+q))

         call geta(qn, qe, t1_try+t1par, ra, dec, t0par)
         call gett(qnp, qep, t1_try+t1par, qlat, qlong,
     +        ra, dec)

         qn = qn + qnp
         qe = qe + qep
         dtau = piex*qn + piey*qe
         dbeta = -piex*qe + piey*qn
         dbeta = -dbeta
         taup  = tau + dtau
         betap = u0 + dbeta

         xcm= taup*dcos(theta)+betap*dsin(theta)
         ycm=-taup*dsin(theta)+betap*dcos(theta)

         xs = xcm*dcos(ep)-ycm*dsin(ep)-offset_x
         ys = xcm*dsin(ep)+ycm*dcos(ep)

         call getbinp(magnew,am,xic,yic,xs,ys,q,b,n_try)

         if (n_try.eq.5) then 
            t1_hi = t1_try
         else
            t1_lo = t1_try
         endif
      enddo

      t1 = (t1_lo+t1_hi)/2.d0
      END



c-------------------------------

      SUBROUTINE REVERSE_T2
     +     (ainp, nparmin, nparmout, b, q, ra, dec, qlat, qlong, 
     +     t2par, t0par, tbinary, t0)
c     Subroutine to find t0 given t2
      implicit none
      integer nparmout, nparmin
      real*8 ainp(nparmin)
      real*8 aoutp(nparmout)
      real*8 t2, u0, te, rhos, piex, piey, theta, db, dep, b, q !input
      real*8 t0 !output

      real*8 t0par, t2par, tbinary, ra, dec, qlat, qlong

      real*8 t0_lo, t0_hi, t0_try
      real*8 t2_lo, t2_hi, t2_try !trial t2s

      integer i, n
      real*8 tol

      tol = 10.d0**(-6) !required precision
      n = 11 !number of loops
c     Initialize variables
      t2    = ainp(1)
      u0    = ainp(2)
      te    = ainp(3)
      rhos  = ainp(4)
      piex  = ainp(5)
      piey  = ainp(6)
      theta = ainp(7)
      db    = ainp(8)
      dep   = ainp(9)

      aoutp(2) = u0
      aoutp(3) = te
      aoutp(4) = rhos
      aoutp(5) = piex 
      aoutp(6) = piey
      aoutp(7) = theta
      aoutp(8) = db
      aoutp(9) = dep
      aoutp(10)= b
      aoutp(11)= q

      t0_lo = -0.2d0
      t0_hi = +0.2d0

      t2_lo = 0d0
      t2_hi = 0d0

 100  continue
      aoutp(1) = t0_lo
      call caustic2
     +     (aoutp, nparmout, ra, dec, qlat, qlong, t2par, t0par,
     +     tbinary, t2_lo)
      
      if (t2.lt.t2_lo) then
         t0_lo = t0_lo-0.1d0
         t0_hi = t0_hi-0.1d0
         go to 100
      endif

 200  continue
      aoutp(1) = t0_hi
      call caustic2
     +     (aoutp, nparmout, ra, dec, qlat, qlong, t2par, t0par,
     +     tbinary, t2_hi)

      if (t2.gt.t2_hi) then
         t0_hi = t0_hi+0.1d0
         go to 200
      endif

c      n = int(dlog((t0_hi-t0_lo)/tol)/dlog(2d0))
c      if (n.lt.15) n=15

c     Binary search for t0
      do i = 1, n 
         t0_try = (t0_hi+t0_lo)/2d0
         aoutp(1) = t0_try
         call caustic2
     +        (aoutp, nparmout, ra, dec, qlat, qlong, t2par, t0par,
     +        tbinary, t2_try)
         
         if (t2_try.gt.t2) then
            t0_hi = t0_try
         else
            t0_lo = t0_try
         endif

      enddo
      
      t0 = (t0_hi+t0_lo)/2d0
      END

      SUBROUTINE CAUSTIC2 
     +     (aoutp, nparmout, ra, dec, qlat, qlong, t2par, t0par,tbinary,
     +     t2)
c     subroutine to find the time of the 2nd caustic crossing given t0.

c     Input: aoutp, nparmout
c     Output: t2
      implicit none
      integer nparmout
      real*8 aoutp(nparmout)

      real*8 t0, u0, te, rhos, piex, piey, theta, db, dep, b_0, q !input
      real*8 b, ep, ep_0, tbinary, dcm
      real*8 t2 !Output

      real*8 t2_lo, t2_hi !Trial t2s.
      integer n_lo, n_hi !number of images for trial t2s.
      real*8 t2_try
      integer n_try
      integer i,n

      real*8 offset_x, offset_y 
c     offset of center of mass from b/2
      real*8 tau, xs, ys, xcm, ycm

      real*8 magnew
      integer max_sol
      parameter (max_sol=5)
      real*8 am(max_sol), xic(max_sol), yic(max_sol)!other getbinp parameters
      real*8 qn, qe, qnp, qep, taup, betap
      real*8 dtau, dbeta

      real*8 ra, dec, qlat, qlong, t0par, t2par
      real*8 tol

c-----Fixed parameters
      tol = 10.d0**(-6) !required precision
      n = 11 !number of loops

      ep_0=0d0 !reference binary rotation
c      tbinary = 4716.7d0 !reference time (where the map is calculated)
c---------------------

      t0    = aoutp(1)
      u0    = aoutp(2)
      te    = aoutp(3)
      rhos  = aoutp(4)
      piex  = aoutp(5)
      piey  = aoutp(6)
      theta = aoutp(7)
      db    = aoutp(8)
      dep   = aoutp(9)
      b_0   = aoutp(10)
      q     = aoutp(11)

      t2 = 0d0

      t2_lo = -0.3d0
      t2_hi = 0d0

c     Position of the center of mass
      offset_y = 0.d0

c     Check that the boundaries are good.
 100  continue
         tau = ((t2_lo-t0)+(t2par-t0par))/te

         if (t2_lo+t2par.lt.4660.)then
            b = b_0+db*(4660.-tbinary)/365.25
            ep = ep_0 + dep*(4660.-tbinary)/365.25
         else if (t2_lo+t2par.gt.4800.)then
            b = b_0+db*(4800.-tbinary)/365.25
            ep = ep_0 + dep*(4800.-tbinary)/365.25
         else
            b = b_0+db*(t2_lo+t2par-tbinary)/365.25
            ep = ep_0 + dep*(t2_lo+t2par-tbinary)/365.25
         endif

         offset_x = b*(-0.5d0+1d0/(1d0+q))

         call geta(qn, qe, t2_lo+t2par, ra, dec, t0par)
         call gett(qnp, qep, t2_lo+t2par, qlat, qlong,
     +        ra, dec)

         qn = qn + qnp
         qe = qe + qep
         dtau = piex*qn + piey*qe
         dbeta = -piex*qe + piey*qn
         dbeta = -dbeta
         taup  = tau + dtau
         betap = u0 + dbeta

         xcm= taup*dcos(theta)+betap*dsin(theta)
         ycm=-taup*dsin(theta)+betap*dcos(theta)

         xs = xcm*dcos(ep)-ycm*dsin(ep)-offset_x
         ys = xcm*dsin(ep)+ycm*dcos(ep)

         call getbinp(magnew,am,xic,yic,xs,ys,q,b,n_lo)

         if (n_lo.ne.5) then
            t2_lo = t2_lo-0.1d0
            t2_hi = t2_hi-0.1d0
            go to 100
         endif

 200     continue
         tau = ((t2_hi-t0)+(t2par-t0par))/te

         if (t2_hi+t2par.lt.4660.)then
            b = b_0+db*(4660.-tbinary)/365.25
            ep = ep_0 + dep*(4660.-tbinary)/365.25
         else if (t2_hi+t2par.gt.4800.)then
            b = b_0+db*(4800.-tbinary)/365.25
            ep = ep_0 + dep*(4800.-tbinary)/365.25
         else
            b = b_0+db*(t2_hi+t2par-tbinary)/365.25
            ep = ep_0 + dep*(t2_hi+t2par-tbinary)/365.25
         endif

         offset_x = b*(-0.5d0+1d0/(1d0+q))

         call geta(qn, qe, t2_hi+t2par, ra, dec, t0par)
         call gett(qnp, qep, t2_hi+t2par, qlat, qlong,
     +        ra, dec)

         qn = qn + qnp
         qe = qe + qep
         dtau = piex*qn + piey*qe
         dbeta = -piex*qe + piey*qn
         dbeta = -dbeta
         taup  = tau + dtau
         betap = u0 + dbeta

         xcm= taup*dcos(theta)+betap*dsin(theta)
         ycm=-taup*dsin(theta)+betap*dcos(theta)

         xs = xcm*dcos(ep)-ycm*dsin(ep)-offset_x
         ys = xcm*dsin(ep)+ycm*dcos(ep)

         call getbinp(magnew,am,xic,yic,xs,ys,q,b,n_hi)

         if (n_hi.ne.3) then 
            t2_hi = t2_hi+0.1d0
            go to 200
         endif

c      n = int(dlog((t2_hi-t2_lo)/tol)/dlog(2d0))
c      if (n.lt.15) n=15

c     Use a binary search to find t2 given t0.
      do i = 1, n
         t2_try = (t2_lo+t2_hi)/2.d0
         
         tau = ((t2_try-t0)+(t2par-t0par))/te
         
         if (t2_try+t2par.lt.4660.)then
            b = b_0+db*(4660.-tbinary)/365.25
            ep = ep_0 + dep*(4660.-tbinary)/365.25
         else if (t2_try+t2par.gt.4800.)then
            b = b_0+db*(4800.-tbinary)/365.25
            ep = ep_0 + dep*(4800.-tbinary)/365.25
         else
            b = b_0+db*(t2_try+t2par-tbinary)/365.25
            ep = ep_0 + dep*(t2_try+t2par-tbinary)/365.25
         endif

         offset_x = b*(-0.5d0+1d0/(1d0+q))

         call geta(qn, qe, t2_try+t2par, ra, dec, t0par)
         call gett(qnp, qep, t2_try+t2par, qlat, qlong,
     +        ra, dec)

         qn = qn + qnp
         qe = qe + qep
         dtau = piex*qn + piey*qe
         dbeta = -piex*qe + piey*qn
         dbeta = -dbeta
         taup  = tau + dtau
         betap = u0 + dbeta

         xcm= taup*dcos(theta)+betap*dsin(theta)
         ycm=-taup*dsin(theta)+betap*dcos(theta)

         xs = xcm*dcos(ep)-ycm*dsin(ep)-offset_x
         ys = xcm*dsin(ep)+ycm*dcos(ep)

         call getbinp(magnew,am,xic,yic,xs,ys,q,b,n_try)

         if (n_try.eq.3) then 
            t2_hi = t2_try
         else
            t2_lo = t2_try
         endif

      enddo

      t2 = (t2_lo+t2_hi)/2.d0

      END
