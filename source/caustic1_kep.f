      SUBROUTINE REVERSE_T1
     +     (ainp, nparmin, nparmout, b, q, ra, dec, qlat, qlong,
     +     t0par, t1par, t2par, tbinary, thetastar, te_0, te, pass)
c     Subroutine to find te given t1
      implicit none
      integer nparmout, nparmin
      real*8 ainp(nparmin)
      real*8 aoutp(nparmout)
      real*8 t2, u0, t1, rhos, piex, piey, theta, db, dep, b, q !input
      real*8 ri3, vi3
      real*8 te, t0 !output

      real*8 t0par, t1par, t2par, tbinary, ra, dec, qlat, qlong
      real*8 pis, thetastar

      real*8 te_lo, te_hi, te_try, te_0
      real*8 t1_lo, t1_hi, t1_try !trial t2s

      integer i, n
      real*8 tol
      logical pass
      
      pass = .true.

c      write(6,*) 'r1'
      
      tol = 10.d0**(-6) !required precision
      n = 12 !number of loops

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

      ri3   = ainp(12)
      vi3   = ainp(13)
      pis   = 1d0/ainp(14)

      if(te_0.lt.1d0)te_0=33.0

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
      aoutp(12)= ri3
      aoutp(13)= vi3
      aoutp(14)= ainp(14)

      te_lo = te_0 - 0.2d0 
      te_hi = te_0 + 0.2d0

c$$$      if(te_lo.lt.1.0d0) then
c$$$         te_lo=1d0
c$$$         te_hi=11d0
c$$$      endif

 100  continue
c      write(6,*) 'te', te_lo, te_hi
      if(pass) then

      aoutp(3) = te_lo
      call caustic1
     +     (aoutp, nparmin, nparmout, ra, dec, qlat, qlong, 
     +     t0par, t1par, t2par, tbinary, thetastar, t1_lo,pass)

      if (t1_lo.lt.t1) then
         te_lo = te_lo-0.2d0
         te_hi = te_hi-0.2d0
         if(te_lo.lt.1.0d0) then
            pass = .false.
         endif
         go to 100
      endif

      endif

 200  continue
      if(pass) then 
c      write(6,*) 'te', te_lo, te_hi
      aoutp(3) = te_hi
      call caustic1
     +     (aoutp, nparmin, nparmout, ra, dec, qlat, qlong, 
     +     t0par, t1par, t2par, tbinary, thetastar, t1_hi,pass)

      if(t1_hi.gt.t1) then
         te_hi = te_hi+0.2d0
         if(te_hi.gt.500.) then
            pass = .false.
         endif
         goto 200
      endif
      endif

c      n = int(dlog((te_hi-te_lo)/tol)/dlog(2d0))
c     if (n.lt.15) n=15

      if(pass) then
      do i = 1, n 
         if (pass) then
         te_try = (te_lo+te_hi)/2d0

         aoutp(3) = te_try
         call caustic1
     +        (aoutp, nparmin, nparmout, ra, dec, qlat, qlong, 
     +        t0par, t1par, t2par, tbinary, thetastar, t1_try,pass)

         if (t1_try.lt.t1) then
            te_hi = te_try
         else
            te_lo  = te_try
         endif
         
         endif
      enddo
      endif

      te = (te_lo+te_hi)/2d0

c      write(6,*) 'r1 done'

      END


      SUBROUTINE CAUSTIC1
     +     (aoutp, nparmin, nparmout, ra, dec, qlat, qlong, 
     +     t0par, t1par, t2par, tbinary, thetastar, t1,pass)
c     Given u0, calculate the time of the first caustic crossing
      implicit none
      real*8 PI 
      parameter (PI=3.141592653589793238462643d0)

      integer nparmout, nparmin
      real*8 aoutp(nparmout)
      real*8 ainp(nparmin)

      real*8 t2, u0, te, rhos, piex, piey, theta, db, dep, b_0, q !input
      real*8 b, ep, ep_0, tbinary, dcm
      real*8 pis, thetastar, psi
      real*8 t1 !Output
      real*8 t0 !internal

      real*8 ri(3), vi(3)
      real*8 rot(3,3), a, omega, ecc, tperi
      real*8 qm, dl, thetae
      real*8 kappa
      parameter(kappa=8.14d0)
      real*8 p1, p2 !x1 and x2 for keppos (name change b/c x2 is in use)


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
      logical pass

      pass = .true.

c      write(6,*) 'c1'
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
      ri(1) = b_0
      ri(2) = 0.
      ri(3) = aoutp(12)
      vi(1) = db/(365.25)
      vi(2) = -dep*b_0/(365.25)
      vi(3) = aoutp(13)
      pis   = 1d0/aoutp(14)

      thetae = thetastar/rhos
      qm = thetae/(kappa*sqrt(piex**2+piey**2))
      dl = 1d0/(thetae*sqrt(piex**2+piey**2)+pis)

      call kepler(rot,a,omega,ecc,tperi,qm,dl,psi,ri,vi,piex,piey,
     +     thetae,pis,tbinary)

      ainp(1) = t2
      ainp(2) = u0
      ainp(3) = te
      ainp(4) = rhos
      ainp(5) = piex
      ainp(6) = piey
      ainp(7) = theta
      ainp(8) = db
      ainp(9) = dep

      ainp(12)= ri(3)
      ainp(13)= vi(3)
      ainp(14)= aoutp(14)

c      write(6,*) 'caustic 1: t2, te: ', t2, te

c     Find t0
      call reverse_t2
     +     (ainp, nparmin, nparmout, b_0, q, ra, dec, qlat, qlong,
     +      t2par, t0par, tbinary, thetastar, t0,pass)

      t1 = 0d0

      t1_lo = -0.1d0
      t1_hi = 0.2d0

c     Position of the center of mass
      offset_y = 0.d0

c     Check that the boundaries are good.
 100  continue

      if (pass) then 
         tau = ((t1_lo-t0)+(t1par-t0par))/te

c            if (t1_lo+t1par.lt.4660.) then
c               call keppos(p1,p2,a,omega,ecc,tperi,rot,4660.d0)
c            else if (t1_lo+t1par.gt.4800.) then
c               call keppos(p1,p2,a,omega,ecc,tperi,rot,4800.d0)
c            else
               call keppos(p1,p2,a,omega,ecc,tperi,rot,t1_lo+t1par)
c            endif

            b = sqrt(p1**2+p2**2)
            ep = atan(p2/p1)
            if (p1.lt.0d0) ep=ep+PI !Check quadrant

c$$$         if (t1_lo+t1par.lt.4660.)then
c$$$            b = b_0+db*(4660.-tbinary)/365.25
c$$$            ep = ep_0 + dep*(4660.-tbinary)/365.25
c$$$         else if (t1_lo+t1par.gt.4800.)then
c$$$            b = b_0+db*(4800.-tbinary)/365.25
c$$$            ep = ep_0 + dep*(4800.-tbinary)/365.25
c$$$         else
c$$$            b = b_0+db*(t1_lo+t1par-tbinary)/365.25
c$$$            ep = ep_0 + dep*(t1_lo+t1par-tbinary)/365.25
c$$$         endif

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
            if(t1_lo.lt.-50.) then
               pass = .false.
            endif
            go to 100
         endif

      endif

 200  continue
      if (pass) then 

         tau = ((t1_hi-t0)+(t1par-t0par))/te

c            if (t1_hi+t1par.lt.4660.) then
c               call keppos(p1,p2,a,omega,ecc,tperi,rot,4660.d0)
c            else if (t1_hi+t1par.gt.4800.) then
c               call keppos(p1,p2,a,omega,ecc,tperi,rot,4800.d0)
c            else
               call keppos(p1,p2,a,omega,ecc,tperi,rot,t1_hi+t1par)
c            endif

            b = sqrt(p1**2+p2**2)
            ep = atan(p2/p1)
            if (p1.lt.0d0) ep=ep+PI !Check quadrant

c$$$         if (t1_hi+t1par.lt.4660.)then
c$$$            b = b_0+db*(4660.-tbinary)/365.25
c$$$            ep = ep_0 + dep*(4660.-tbinary)/365.25
c$$$         else if (t1_hi+t1par.gt.4800.)then
c$$$            b = b_0+db*(4800.-tbinary)/365.25
c$$$            ep = ep_0 + dep*(4800.-tbinary)/365.25
c$$$         else
c$$$            b = b_0+db*(t1_hi+t1par-tbinary)/365.25
c$$$            ep = ep_0 + dep*(t1_hi+t1par-tbinary)/365.25
c$$$         endif
         
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
            if(t1_hi.gt.50) then
               pass = .false.
            endif
            go to 200
         endif
      endif

c      n = int(dlog((t1_hi-t1_lo)/tol)/dlog(2d0))
c      if (n.lt.15) n=15

c     Use a binary search to find t1 given t0.
      if(pass) then
      do i = 1, n
         t1_try = (t1_lo+t1_hi)/2.d0
         
         tau = ((t1_try-t0)+(t1par-t0par))/te

c            if (t1_try+t1par.lt.4660.) then
c               call keppos(p1,p2,a,omega,ecc,tperi,rot,4660.d0)
c            else if (t1_try+t1par.gt.4800.) then
c cc              call keppos(p1,p2,a,omega,ecc,tperi,rot,4800.d0)
c            else
               call keppos(p1,p2,a,omega,ecc,tperi,rot,t1_try+t1par)
c            endif

            b = sqrt(p1**2+p2**2)
            ep = atan(p2/p1)
            if (p1.lt.0d0) ep=ep+PI !Check quadrant

c$$$         if (t1_try+t1par.lt.4660.)then
c$$$            b = b_0+db*(4660.-tbinary)/365.25
c$$$            ep = ep_0 + dep*(4660.-tbinary)/365.25
c$$$         else if (t1_try+t1par.gt.4800.)then
c$$$            b = b_0+db*(4800.-tbinary)/365.25
c$$$            ep = ep_0 + dep*(4800.-tbinary)/365.25
c$$$         else
c$$$            b = b_0+db*(t1_try+t1par-tbinary)/365.25
c$$$            ep = ep_0 + dep*(t1_try+t1par-tbinary)/365.25
c$$$         endif

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

      endif

      t1 = (t1_lo+t1_hi)/2.d0

c      write(6,*) 'c1 done'
      END



c-------------------------------

      SUBROUTINE REVERSE_T2
     +     (ainp, nparmin, nparmout, b, q, ra, dec, qlat, qlong, 
     +     t2par, t0par, tbinary, thetastar, t0,pass)
c     Subroutine to find t0 given t2
      implicit none
      integer nparmout, nparmin
      real*8 ainp(nparmin)
      real*8 aoutp(nparmout)
      real*8 t2, u0, te, rhos, piex, piey, theta, db, dep, b, q !input
      real*8 ri3, vi3
      real*8 t0 !output

      real*8 t0par, t2par, tbinary, ra, dec, qlat, qlong
      real*8 pis, thetastar

      real*8 t0_lo, t0_hi, t0_try
      real*8 t2_lo, t2_hi, t2_try !trial t2s

      integer i, n
      real*8 tol
      logical pass

      pass = .true.
c      write(6,*) 'r2'

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

      ri3 = ainp(12)
      vi3 = ainp(13)
      pis = 1d0/ainp(14)


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
      aoutp(12)= ri3
      aoutp(13)= vi3
      aoutp(14)= ainp(14)

      t0_lo = -0.2d0
      t0_hi = +0.2d0

 100  continue
c      write(6,*) 't0', t0_lo, t0_hi
      if (pass) then 
      aoutp(1) = t0_lo
      call caustic2
     +     (aoutp, nparmout, ra, dec, qlat, qlong, t2par, t0par,
     +     tbinary, thetastar, t2_lo,pass)
      
      if (t2.lt.t2_lo) then
         t0_lo = t0_lo-0.1d0
         t0_hi = t0_hi-0.1d0

         if(t0_lo.lt.-50.)then
            pass = .false.
         endif

         go to 100
      endif

      endif

 200  continue
      if (pass) then
c      write(6,*) 't0', t0_lo, t0_hi
      aoutp(1) = t0_hi
      call caustic2
     +     (aoutp, nparmout, ra, dec, qlat, qlong, t2par, t0par,
     +     tbinary, thetastar, t2_hi,pass)

      if (t2.gt.t2_hi) then
         t0_hi = t0_hi+0.1d0
         if(t0_hi.gt.50.)then
            pass = .false.
         endif
         go to 200
      endif

      endif
c      n = int(dlog((t0_hi-t0_lo)/tol)/dlog(2d0))
c      if (n.lt.15) n=15

c     Binary search for t0
      if (pass) then
      do i = 1, n 
         if(pass) then
            t0_try = (t0_hi+t0_lo)/2d0
            aoutp(1) = t0_try
            call caustic2
     +           (aoutp, nparmout, ra, dec, qlat, qlong, t2par, t0par,
     +           tbinary, thetastar, t2_try,pass)
            
            if (t2_try.gt.t2) then
               t0_hi = t0_try
            else
               t0_lo = t0_try
            endif
         endif
      enddo
      
      endif

      t0 = (t0_hi+t0_lo)/2d0

c      write(6,*) 'r2 done'
      END

      SUBROUTINE CAUSTIC2 
     +     (aoutp, nparmout, ra, dec, qlat, qlong, t2par, t0par,tbinary,
     +      thetastar, t2,pass)
c     subroutine to find the time of the 2nd caustic crossing given t0.

c     Input: aoutp, nparmout
c     Output: t2
      implicit none
      integer nparmout
      real*8 PI 
      parameter (PI=3.141592653589793238462643d0)

      real*8 aoutp(nparmout)

      real*8 t0, u0, te, rhos, piex, piey, theta, db, dep, b_0, q !input
      real*8 b, ep, ep_0, tbinary, dcm
      real*8 pis, thetastar, psi
      real*8 t2 !Output

      real*8 ri(3), vi(3)
      real*8 rot(3,3), a, omega, ecc, tperi
      real*8 qm, dl, thetae
      real*8 kappa
      parameter(kappa=8.14d0)
      real*8 p1, p2 !x1 and x2 for keppos (name change b/c x2 is in use)

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

      real*8 b_test, ep_test
      logical pass

      pass = .true.

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
      ri(1) = b_0
      ri(2) = 0.d0
      ri(3) = aoutp(12)
      vi(1) = db/(365.25d0)
      vi(2) = -dep*b_0/(365.25d0)
      vi(3) = aoutp(13)
      pis   = 1d0/aoutp(14)

      thetae = thetastar/rhos
      qm = thetae/(kappa*sqrt(piex**2+piey**2))
      dl = 1d0/(thetae*sqrt(piex**2+piey**2)+pis)

      call kepler(rot,a,omega,ecc,tperi,qm,dl,psi,ri,vi,piex,piey,
     +     thetae,pis,tbinary)

c      call keppos(p1,p2,a,omega,ecc,tperi,rot,tbinary)
c$$$      
c$$$         b = sqrt(p1**2+p2**2)
c$$$         ep = atan(p2/p1)
c$$$         if (p1.lt.0d0) ep=ep+PI !Check quadrant
c      if(abs(b_0-b).gt.0.00001) then
c         write(6,*) aoutp
c         write(6,*) b_0, b, ep
c         write(6,*) vi
c         write(6,*) p1, p2
c      endif

      t2 = 0d0

      t2_lo = -0.3d0
      t2_hi = 0d0

c     Position of the center of mass
      offset_y = 0.d0


c     Check that the boundaries are good.
 100  continue
c      write(6,*) 't2', t2_lo, t2_hi
      if (pass) then 
         tau = ((t2_lo-t0)+(t2par-t0par))/te

c         if (t2_lo+t2par.lt.4660.) then
c            call keppos(p1,p2,a,omega,ecc,tperi,rot,4660.d0)
c         else if (t2_lo+t2par.gt.4800.) then
c            call keppos(p1,p2,a,omega,ecc,tperi,rot,4800.d0)
c         else
            call keppos(p1,p2,a,omega,ecc,tperi,rot,t2_lo+t2par)
c         endif

         b = sqrt(p1**2+p2**2)
         ep = atan(p2/p1)
         if (p1.lt.0d0) ep=ep+PI !Check quadrant

c$$$         if (t2_lo+t2par.lt.4660.)then
c$$$            b_test = b_0+db*(4660.-tbinary)/365.25
c$$$            ep_test = ep_0 + dep*(4660.-tbinary)/365.25
c$$$         else if (t2_lo+t2par.gt.4800.)then
c$$$            b_test = b_0+db*(4800.-tbinary)/365.25
c$$$            ep_test = ep_0 + dep*(4800.-tbinary)/365.25
c$$$         else
c$$$            b_test = b_0+db*(t2_lo+t2par-tbinary)/365.25
c$$$            ep_test = ep_0 + dep*(t2_lo+t2par-tbinary)/365.25
c$$$         endif

c$$$         write(6,1000) b_0, t2_lo, b-b_test, ep-ep_test, b, ep, b_test, 
c$$$     +        ep_test
c$$$c         write(6,1000) t2_lo+t2par-tbinary, b_0,b, vi
 1000    format(f8.4, 2f9.5,4f9.5)


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
            if(t2_lo.lt.-50.)then
               pass = .false.
            endif
            go to 100
         endif

      endif


 200     continue

      if (pass) then
c         write(6,*) 't2', t2_lo, t2_hi
         tau = ((t2_hi-t0)+(t2par-t0par))/te

c         if (t2_hi+t2par.lt.4660.) then
c            call keppos(p1,p2,a,omega,ecc,tperi,rot,4660.d0)
c         else if (t2_hi+t2par.gt.4800.) then
c            call keppos(p1,p2,a,omega,ecc,tperi,rot,4800.d0)
c         else
            call keppos(p1,p2,a,omega,ecc,tperi,rot,t2_hi+t2par)
c         endif

         b = sqrt(p1**2+p2**2)
         ep = atan(p2/p1)
         if (p1.lt.0d0) ep=ep+PI !Check quadrant

c$$$         if (t2_hi+t2par.lt.4660.)then
c$$$            b_test = b_0+db*(4660.-tbinary)/365.25
c$$$            ep_test = ep_0 + dep*(4660.-tbinary)/365.25
c$$$         else if (t2_hi+t2par.gt.4800.)then
c$$$            b_test = b_0+db*(4800.-tbinary)/365.25
c$$$            ep_test = ep_0 + dep*(4800.-tbinary)/365.25
c$$$         else
c$$$            b_test = b_0+db*(t2_hi+t2par-tbinary)/365.25
c$$$            ep_test = ep_0 + dep*(t2_hi+t2par-tbinary)/365.25
c$$$         endif
c$$$
c$$$         write(6,1000) t2_hi, b-b_test, ep-ep_test, b, ep, b_test, 
c$$$     +        ep_test
c$$$ 1000    format(f8.4, 2f9.5,4f9.5)

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
            if(t2_hi.gt.50.)then
               pass = .false.
            endif
            go to 200
         endif

      endif

c      n = int(dlog((t2_hi-t2_lo)/tol)/dlog(2d0))
c      if (n.lt.15) n=15

c     Use a binary search to find t2 given t0.
c$$$      if((t2_lo.lt.-10.0).or.(t2_hi.gt.10.0)) then
c$$$         write(6,*) 'Potato!', t2_lo, t2_hi
c$$$      endif

      if(pass) then
      do i = 1, n
         t2_try = (t2_lo+t2_hi)/2.d0
         
         tau = ((t2_try-t0)+(t2par-t0par))/te
         
c         if (t2_try+t2par.lt.4660.) then
c            call keppos(p1,p2,a,omega,ecc,tperi,rot,4660.d0)
c         else if (t2_try+t2par.gt.4800.) then
c            call keppos(p1,p2,a,omega,ecc,tperi,rot,4800.d0)
c         else
            call keppos(p1,p2,a,omega,ecc,tperi,rot,t2_try+t2par)
c         endif
         
         b = sqrt(p1**2+p2**2)
         ep = atan(p2/p1)
         if (p1.lt.0d0) ep=ep+PI !Check quadrant

c$$$         if (t2_try+t2par.lt.4660.)then
c$$$            b = b_0+db*(4660.-tbinary)/365.25
c$$$            ep = ep_0 + dep*(4660.-tbinary)/365.25
c$$$         else if (t2_try+t2par.gt.4800.)then
c$$$            b = b_0+db*(4800.-tbinary)/365.25
c$$$            ep = ep_0 + dep*(4800.-tbinary)/365.25
c$$$         else
c$$$            b = b_0+db*(t2_try+t2par-tbinary)/365.25
c$$$            ep = ep_0 + dep*(t2_try+t2par-tbinary)/365.25
c$$$         endif

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

      endif

      t2 = (t2_lo+t2_hi)/2.d0


      END
