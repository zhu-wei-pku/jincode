      program map_mcmc
      implicit none
      real*8 chi2tot
      integer nparmout
      parameter (nparmout = 11) !number of model parameters
      real*8 aoutp(nparmout)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer nobs, nob, ndata
      parameter (ndata=10001)
      parameter (nobs =1)
      real*8    date(ndata,nobs),flux(ndata,nobs),err(ndata,nobs),
     + qmag(ndata,nobs), tol(ndata,nobs), mgamma(nobs)
      integer ndat(nobs)
      logical bad(ndata,nobs)
      real*8 fs(nobs), fb(nobs)
      real*8 chi2(nobs)
      integer kmax(nobs)
      real*8 chi2max(nobs)
      logical iflag
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 tmin, tmax, width, radius12, radius22
      real*8 mrhomin,mgridsize,s_r_max, thickness
c      real*8 ma, da, mb, db, phi
      real*8 b, q, alpha
      real*8 b_0, ep_0, ep, dcm, db, dep, tbinary
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical magfail, magfail_loop
      real*8 offset_x, offset_y
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer nmodel
      parameter (nmodel = 1)
      integer mmodel
      parameter(mmodel = 6775) !number of points in model lc.
      real*8 date_model(mmodel,nmodel)
      real*8 mag_model(mmodel,nmodel), gamma_model(nmodel)
      logical bad_model(mmodel,nmodel)
      integer ndat_model(nmodel)

      integer traj_port_model(nmodel), traj_port(nobs)

      real*8 cchi2

c     !Added for Jan zroots
      integer MAX_SOL
      parameter (MAX_SOL=5)
      complex*16 zr(max_sol) 


      integer k, i, itemp, knttot
      real*8 fl, res, e, flpre

      real*8 PI
      parameter (PI = 3.1415926535d0)
      real*8 fsogle, fbogle

      real*8 tdraw_min, tdraw_max

c     Han model from April 29, 2013
c      data aoutp/6406.35,0.0029,24.91,0.00430,
c     +     0.0, 0.0, 3.38,
c     +     0.0, 0.0, 0.205,  1.603/!,  0.905,  0.000044d0,  -11.344/

c     Best-fit from Arjuna.
c      data aoutp/6406.36214931,  0.0021764795,  24.9106218,
c     +     0.004279057045,  0.000000000,  0.000000000,  3.400965000,
c     +     0.000000,  0.000000,  0.204963,  1.602995 /

c     Input model:
c     t0, u0, te
c     rho, pien, piee, alpha (rad)
c     dsdt, depdt, s, q
      data aoutp/ 23.38,  0.008516,  15.68,
     +     0.004,  0.000000000,  0.000000000,  3.425261,
     +     0.000000,  0.000000,  0.853404,  7.1564e-4/

      real*8 qlat(nobs), qlong(nobs)
      real*8 qlat_model(nmodel), qlong_model(nmodel)
      real*8 ra, dec, t0par, te_0


      logical binflag(ndata,nobs)

      real*8 start, finish
      call cpu_time(start)


      t0par = 6410.d0 !reference time for parallax
      tbinary = t0par

c      aoutp(1) = aoutp(1)+t0par

      write(6,306) (aoutp(itemp),itemp=1,nparmout)

c      call getlatlon(qlat, qlong, nobs)


      traj_port_model(1) = 78

      do nob = 1,nobs
         traj_port(nob) = 78 + nob
c     define the file port for trajectory files
      enddo

c     Added for Jan zroots
      do i = 1, max_sol
         zr(i) = 0d0
      enddo

      tdraw_min = 0
      tdraw_max = 47.041667
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Initializing...
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do nob=1,nobs
         ndat(nob)=0
         do i=1,ndata
            date(i,nob)=0d0
            flux(i,nob)=0d0
            err(i,nob)=0d0
            qmag(i,nob)=0d0
         enddo
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call getdat(date,flux,err,bad,ndat,ndata,nobs,tol,mgamma, 
     +     qlat, qlong, ra, dec)

c      write(6,*) 'gamma', mgamma

      qlat_model(1)  = qlat(1)
      qlong_model(1) = qlong(1)

c      write(6,306) (aoutp(itemp),itemp=1,nparmout)
      
      b = aoutp(10)
      q = aoutp(11)

      write(6,*) 'Getting Magnification...'

c      write(6,306) (aoutp(itemp),itemp=1,nparmout)
 306  format(f12.7, 2f10.6, 1f11.7, 2f10.5, 1f11.7, 2f10.6,2f11.7)
      call cpu_time(start)

      call getmags_binary
     + (qmag,aoutp,nparmout,date,bad,ndata,nobs,ndat,
     +   tmin,tmax,
     +   mgamma,
     +   magfail,
     +   traj_port, 
     +   qlat, qlong, ra, dec, t0par,tbinary,
     +   binflag,zr)
      open(22,file='test-lc')
      do nob=1,nobs
        do k=1,ndat(nob)
            write(22,*) date(k,nob),qmag(k,nob)
        enddo
      enddo

      call cpu_time(finish)
      print '("Time to getmag = ",f10.3," seconds.")',finish-start

      ep = 0.0
      offset_x=b*(-0.5d0+1d0/(1d0+q))
      offset_y = 0.d0
      call gencaustics(b, q, ep, offset_x, offset_y)
c      return
      call getchi2(chi2,fs,fb,qmag,flux,err,bad,ndata,nobs,ndat,
     *      kmax,chi2max,iflag)

      knttot = 0
      do nob=1,nobs
         do k = 1,ndat(nob)
            write(61,302) date(k,nob), qmag(k,nob)
         enddo
         knttot=knttot+ndat(nob)
      enddo
 302  format(f15.8, f10.6)
      close(61)

      write(6,*) '# of data points',  knttot

      chi2tot = 0

      do 720 nob=1,nobs
         chi2tot = chi2tot + chi2(nob)
720    continue
 
      fs(1) = 0.912006
      fb(1) = 1.0-fs(1)
      write(6,*) '--------------------------'
      write(6,*) 'Chi^2 = '
      write(6,310) chi2tot
      write(6,312) (fs(nob),fb(nob),nob=1,nobs)
      write(6,311) (chi2(nob),nob=1,nobs)
      write(6,313) (ndat(nob), nob=1,nobs)
      write(6,*) '--------------------------'
 310  format(f12.4)
 312  format(18f10.6)
 311  format(9f10.2)
 313  format(9i10)



      if (nobs.lt.9) then
         fsogle = fs(1)
         fbogle = fb(1)
      else
         fsogle = fs(8)
         fbogle = fb(8)
      endif
      fsogle = 0.912006
      fbogle = 1-fsogle
c$$$      write(6,*) 'fs, fb ogle ', fsogle, fbogle

      do i=1,nobs
         do k=1,ndat(i)
c            if(.not.(bad(k,i))) then
            if(0.ne.1) then
               fl=flux(k,i)
               res = -2.5d0*dlog10(fl/(fs(i)*qmag(k,i)+fb(i)))
               cchi2 = ((fl - fs(i)*qmag(k,i) - fb(i))/err(k,i))**2
c               if(i.eq.2) then
c                  write(18,*) date(k,i), flux(k,i), fl
c               endif

               e=err(k,i)
c     use this code for making correct lightcurves
               fl=(fl-fb(i))/fs(i)*fsogle+fbogle
c               write(i+36,102) date(k,i), fl, e, res, cchi2, qmag(k,i)
               e=e/fs(i)*fsogle
               e=e/(fl*dlog(10d0)/2.5d0)

c     use this code for getting the correct chi2s.
c               e = e/fl*2.5d0/dlog(10d0)
c               fl=(fl-fb(i))/fs(i)*fsogle+fbogle
               
               fl = 18.d0 - 2.5*log10(fl)
               flpre = fsogle*qmag(k,i)+fbogle
               flpre = 18.-2.5*log10(flpre)
               res = fl-flpre
               
               write(i+36,102) date(k,i), fl, e, res, cchi2,res**2/e**2, 
     +               qmag(k,i)
c            else
c               write(6,*)k, i
            endif
         enddo
      enddo
 102     format(f15.7, 10f10.4)


      write(6,*) 'Generating Model...'
      do i = 1,mmodel

         date_model(i,1) = (tdraw_max - tdraw_min)/(mmodel - 1.d0)*(i-1)
     +        +tdraw_min     
c      date_model(i,1) = tdraw_min + 0.002d0*(i-1)

         bad_model(i,1)  = .false.
      enddo

      ndat_model(1) = mmodel
c      gamma_model(1) = 0.35d0
c      gamma_model(1) = 0.d0
      gamma_model(1) = mgamma(1) !Set model LD to OGLE

      call getmags_binary
     + (mag_model,aoutp,nparmout,date_model,bad_model,
     +  mmodel,nmodel,ndat_model,
     +   tmin,tmax,
     +   gamma_model,
     +   magfail,
     +   traj_port_model,
     +   qlat_model, qlong_model, ra, dec, t0par,tbinary,
     +   binflag,zr)
      write(6,*) '--------------------------'

      write(6,*) 'traj port', traj_port_model

      do i=1,mmodel
         if (binflag(i,1)) then
            fl= mag_model(i,1)*fsogle + fbogle
            fl= 18.d0 - 2.5*log10(fl)
            write(35,101) date_model(i,1), fl, mag_model(i,1)
 101        format(f15.7, f10.5, f10.5)
         endif
      enddo
      write(6,*) '--------------------------'

c$$$      if(b.lt.1)then
c$$$         offset_x=b*(-0.5d0+1d0/(1d0+q))
c$$$      else
c$$$         offset_x=b/2d0 -q/(1d0+q)/b
c$$$      endif
c$$$
c$$$      offset_y = 0.d0

c      tbinary = 4716.85
      b_0 = aoutp(10)
      ep_0 = 0.d0
      db = aoutp(8)
      dep = aoutp(9)

c     First caustic crossing
      b = b_0+db*(6386.2-tbinary)/365.25
      ep = ep_0 + dep*(6386.2-tbinary)/365.25

      offset_x=b*(-0.5d0+1d0/(1d0+q))
      offset_y = 0.d0

      write(6,*) '--------------------------'

      write(6,*) b, q

      WRITE(6,*)'#',offset_x, offset_y

      call gencaustics(b, q, ep, offset_x, offset_y)

1001  continue

                                ! put code to test here
      call cpu_time(finish)
      print '("Time = ",f10.3," seconds.")',finish-start

      end
      
      include './source/getdat.f'
      Include './source/mcmc_sub.f'
      Include './source/jinsubs.f'
      Include './source/getbinp.f'
      Include './source/parsubs.f'
      Include './source/dcaustic_orb.f'
      Include './source/roots10.f'
     
      subroutine getmags_binary
     +           (qmag,aoutp,nparm,date,bad,ndata,nobs,ndat,
     +           tmin,tmax,
     +           mgamma,
     +           magfail,
     +           traj_port,
     +           qlat, qlong, ra, dec, t0par, tbinary, binflag,zr)

      implicit none
      integer*4 nobs,nob,k,ndata,nparm
      real*8 mapgridsize
c     loopgridsize
      real*8 xorigin,yorigin
      integer*4 errorflag
      integer error
      real*8 PI 
      parameter (PI=3.141592653589793238462643d0)
      real*8 rhos,magnew
      real*8 magbps, magbfs
      real*8 mgamma(nobs),gamma,xs,ys,fstol
      real*8 theta,u0,t,m1,m2
      real*8 t0,te,tau
      real*8 qmag(ndata,nobs),date(ndata,nobs),aoutp(nparm)
      logical bad(ndata,nobs)
      integer*4 ndat(nobs)
      real*8 tmin,tmax
      integer MAX_SOL
      parameter (MAX_SOL=5)
      complex*16 zr(max_sol) !Added for Jan zroots
      real*8 xic(MAX_SOL), yic(MAX_SOL), am(MAX_SOL)
      integer*4 nsol
      integer*4 ntab 
      real*8 b0,b1,db0,db1

      integer i

      real*8 offset_x, offset_y
c      real*8 da, db, ma, mb, phi
      real*8 b, q
      real*8 db, dep, b_0, ep_0, ep, tbinary, xcm, ycm

      logical magfail
      integer traj_port(nobs)

      real*8 qlat(nobs), qlong(nobs)
      real*8 qn, qe, piex, piey, qnp, qep, taup, betap
      real*8 dtau, dbeta
      real*8 ra, dec, t0par
      double precision npts ! number of points in the jin loop

      real*8 magtemp

      logical binflag(ndata,nobs)

c     Added for Jan zroots
      do i = 1, max_sol
         zr(i) = 0d0
      enddo

c-----Fixed paramters
      ep_0=0d0 !reference binary rotation
c      tbinary = 4716.85d0 !reference time (where the map is calculated)
c---------------------

      magfail = .false.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
c$$$      b     = aoutp(8)
c$$$      q     = aoutp(9)


      offset_y = 0.d0

      do 100 nob=1,nobs
         gamma = mgamma(nob)
c         write(6,*) 'gamma', gamma
c         if (nob.eq.1) pause
         do 50 k=1,ndat(nob)
c            write(6,*) nob, k, date(k,nob)
            binflag(k,nob)=.true.
c            if(bad(k,nob))go to 50
            t = date(k,nob)
            tau=(t-t0)/te

c$$$            if (t.lt.4660.)then
c$$$               b = b_0+db*(4660.-tbinary)/365.25
c$$$               ep = ep_0 + dep*(4660.-tbinary)/365.25
c$$$            else if (t.gt.4800.)then
c$$$               b = b_0+db*(4800.-tbinary)/365.25
c$$$               ep = ep_0 + dep*(4800.-tbinary)/365.25
c$$$            else
               b = b_0+db*(t-tbinary)/365.25
               ep = ep_0 + dep*(t-tbinary)/365.25
c            endif

            offset_x = b*(-0.5d0+1d0/(1d0+q))

            call geta(qn, qe, t, ra, dec, t0par) !orbital parallax
            call gett(qnp, qep, t, qlat(nob), qlong(nob),
     +           ra, dec) !terrestrial parallax

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

c      write(traj_port(nob),200) t, xs+offset_x, ys+offset_y, nob, rhos
      write(traj_port(nob),200) t, xcm, ycm, nob, rhos
200   format(f12.6, 2f13.7, i6, f13.7)

c     First interval is not in regular code. Only here for plotting purposes.
c       if(t.gt.6386.2d0.and.t.lt.6388.d0) then
c      if ((t.ge.6405.d0).and.(t.le.6407.6))then
cc     Use Jin's code
c
c         if((t.ge.6406.71.and.t.le.6406.78).or.
c     +        (t.ge.6406.87.and.t.le.6406.96))then
c            npts=10000.
c         else if(t.gt.6406.5.and.t.lt.6406.99) then
c            npts=5000.
c         else if(t.ge.6406.99.and.t.lt.6407.03) then
c            npts=500.
c         else
c            npts=0.
c         endif
cc     Use Jin's code
c                  call getmag_jin(magbps, magbfs, xs, ys, b, q, rhos, 
c     +                 gamma, npts, error,zr)
c                  if(error.ne. 0)write(6,*) 'jin error: ', t, error
c                  magnew = magbfs
c      else
c333           continue
c              binflag(k,nob)=.false.

c             if (t.ge.6404.5 .and. t.le.6407.9) then
                call taylor_2(magnew,xs,ys,q,b,gamma,rhos,zr,4)
c             else if  (t.gt.6402. .and. t.lt.6411.) then
c                call taylor_2(magnew,xs,ys,q,b,gamma,rhos,zr,2)
c             else
c                call getbinp(magnew,am,xic,yic,xs,ys,q,b,nsol,zr)
c             endif
c      endif
      qmag(k,nob) = magnew
 50      continue
 100  continue
      return
      end

      subroutine getchi2(chi2,fs,fb,qmag,flux,err,bad,ndata,nobs,ndat,
     *      kmax,chi2max,flag)
      implicit real*8 (a-h,o-z)
      dimension chi2(nobs),fs(nobs),fb(nobs),cov(nobs,2)
      dimension qmag(ndata,nobs)
      dimension flux(ndata,nobs),err(ndata,nobs),ndat(nobs)
      dimension kmax(nobs),chi2max(nobs)
      logical bad(ndata,nobs)
      parameter (ntop=2)
      real*8 a(ntop),b(ntop,ntop),c(ntop,ntop),d(ntop)
     *    ,dum(ntop,ntop),f(ntop)

      real*8 col, qi
      logical flag


      do 110 nob=1,nobs
         do 100 iloop=1,2

            do 10 i=1,ntop
               d(i) = 0
               do 10 j=1,ntop
                  b(i,j) = 0
 10         continue

            chi2(nob) = 0
            chi2max(nob) = 0

            do 85 k =1,ndat(nob)
               if(bad(k,nob))go to 85

               y = flux(k,nob)
               sig2 = err(k,nob)**2
               f(1) = 1d0
               f(2) = qmag(k,nob)
 11            format(i4,4f12.5)
               ypre = 0d0

               do 80 i = 1,ntop
                  ypre = ypre + a(i)*f(i)
                  d(i) = d(i) + y*f(i)/sig2
                  do 80 j = 1,ntop
                     b(i,j) = b(i,j) + f(i)*f(j)/sig2
 80            continue

               dif = y - ypre
               chi2add = dif**2/sig2

               if(chi2add.gt.chi2max(nob))then
                  chi2max(nob) = chi2add
                  kmax(nob) = k
               endif

               chi2(nob) = chi2(nob) + chi2add
 85         continue

            if(iloop.eq.2)go to 100
            call inv2(b,c,dum,ntop,ntop,flag)

            if(flag) then
               do k=1,nobs
                  chi2(k)=999999
               enddo
               return
            endif

            do 90 i=1,ntop
               a(i) = 0d0
               do 90 j=1,ntop
                  a(i) = a(i) + c(i,j)*d(j)
 90         continue

 100     continue

c         write(6,*) c

         fb(nob) = a(1)
         fs(nob) = a(2)

         cov(nob,1) = c(1,1)
         cov(nob,2) = c(2,2)
        
 110  continue

c     3 = CTIO I, 4 = CTIO V
c$$$      nI=5
c$$$      nV=6
c$$$      col = -2.5*dlog10(fs(nV)/fs(nI))
c$$$      ec = sqrt(cov(nV,2)/fs(nV)**2+cov(nI,2)/fs(nI)**2)
c$$$      qi = 18 - 2.5*dlog10(fs(nI))
c$$$      eq = sqrt(cov(nI,2))/fs(nI)
c$$$      write(6,88) 'Source: ', col, ec, qi, eq
c$$$
c$$$      col = -2.5*dlog10(fb(nV)/fb(nI))
c$$$      ec = sqrt(cov(nV,1)/fb(nV)**2+cov(nI,1)/fb(nI)**2)
c$$$      qi = 18 - 2.5*dlog10(fb(nI))
c$$$      eq = sqrt(cov(nI,1))/fb(nI)
c$$$      write(6,88) 'Blend: ', col, ec, qi, eq

 88   format(a8, 2(f7.3,f7.4))
      

      return
      end

        subroutine inv2(a,ainv,adum,n,nmat,flag)
        real*8 a(nmat,nmat),adum(nmat,nmat),ainv(nmat,nmat)
        real*8 hdum(1000),hinv(1000),Q
        data tol/0.00000000000001/
        logical flag
        flag=.false.
        if(n.gt.1000)then
                write(97,*)n,' bigger than 1000'
		call flush(97)
                flag = .true.
		return
        endif

        do  10 i=1,n
           do 10 j=1,n 
              adum(i,j) = a(i,j)
              ainv(i,j) = 0d0
   10   continue

        do 15 i=1,n
           ainv(i,i)=1
 15     continue

        do 100 j=n,1,-1
           if(abs(adum(j,j)).gt.tol) go to 40

           do 20 j1=j-1,1,-1
              if(abs(adum(j,j1)).gt.tol) go to 21
 20        continue

cc      write(6,*)j, 'degenerate matrix'
        flag=.true.
        return
   21   continue

        do 25 i=1,n
           hdum(i) = adum(i,j)
           hinv(i) = ainv(i,j)
   25   continue

        do 30 i=1,n
           adum(i,j) = adum(i,j1)
           ainv(i,j) = ainv(i,j1)
   30   continue

        do 35 i=1,n
           adum(i,j1) = hdum(i)
           ainv(i,j1) = hinv(i)
   35   continue

   40   continue

        q = 1/adum(j,j)
        do 45 i=1,j
           adum(i,j) = adum(i,j)*q
           hdum(i) = -adum(i,j)
 45     continue

        do 50 i=1,n
           ainv(i,j) = ainv(i,j)*q
           hinv(i) = -ainv(i,j)
   50   continue

        do 70 j1 = 1,j-1
           q = adum(j,j1)
           if(q.eq.0.0)go to 69
           do 55 i=1,j
              adum(i,j1) = adum(i,j1) + q*hdum(i)
 55        continue

           do 60 i=1,n
              ainv(i,j1) = ainv(i,j1) + q*hinv(i)
 60        continue

 69        continue
 70     continue

 100  continue

      do 200 j=1,n-1
         do 110 i=1,n
            hinv(i) = ainv(i,j)
 110     continue

         do 130 j1=j+1,n
            q = -adum(j,j1)
            if(q.eq.0.0)go to 129
            do 120 i=1,n
               ainv(i,j1) = ainv(i,j1) + q*hinv(i)
 120        continue

 129        continue

 130     continue
 200  continue
      
      do 220 i=1,n
         do 220 j=1,n
            adum(i,j)=0
            do 220 k=1,n
               adum(i,j) = adum(i,j) + a(i,k)*ainv(k,j)
  220   continue

        return 
        end 


c



       subroutine trimstring(string,ndim,iend)
       implicit none
       integer iend, ndim
       character string(ndim)
       integer i
       iend = 1
       do i = ndim,1,-1
       if(string(i).ne.' ') then
       iend = i
       goto 10
       endif
       enddo
10     continue
       return
       end 

      subroutine getqfac(qfac,b)
      implicit real*8 (a-h,o-z)
      bt = b+1/b
      eps = bt/2 - 1
      cosphi = 3/4.*bt*(1-sqrt(1 - 32/9./bt**2))
      sinphi = sqrt(1-cosphi**2)
      qfac = 4*sinphi**3/(bt - 2*cosphi)**2
c      write(6,8)eps,cosphi,1-3*eps,qfac,sqrt(27/32./eps) !+18*eps**2
 8    format(6f10.6)
      return
      end


      subroutine taylor_2(qmag, x, y, q, d, gamma,
     + rho,zr,ntaylor)
      implicit real*8 (a-h,o-z)
      real*8 q, d
      integer nsol
      integer MAX_SOL
      parameter(MAX_SOL = 5)
      real*8 xic(MAX_SOL), yic(MAX_SOL), am(MAX_SOL)
ccc      call
ccc     + getimt(xic, yic, am, 
ccc     + nsol, da, db, phi, ma, mb, x, y, qmag)
      call getbinp(qmag,am,xic,yic,x,y,q,d,nsol,zr)
      if(ntaylor.eq.0)return
      qmag0 = qmag
      call getquad(qmag1,x,y,q,d,rho,zr,0)
      qmag1 = qmag1 - qmag0
      if(ntaylor.eq.2)then
         qmag = qmag0 + qmag1/2.d0*(1.d0 - gamma/5.d0)
c         write(6,*) 'taylor 2', qmag, qmag0+qmag1/2d0, gamma
         return
      endif
      if(ntaylor.eq.4)then
         call getquad(qmag1p,x,y,q,d,rho,zr,1)
        qmag1p = qmag1p - qmag0
         call getquad(qmaghalf,x,y,q,d,rho/2,zr,0)
         qmaghalf = qmaghalf - qmag0
         a = (16*qmaghalf - qmag1)/3
         b = (qmag1 + qmag1p)/2 - a
         qmag = qmag0 + a/2.d0*(1.d0 - gamma/5.d0)
     *        + b/3.d0*(1.d0 - gamma*11.d0/35.d0)
c         write(6,*) 'taylor 4', qmag, qmag0+a/2d0+b/3.d0, gamma
         return
      endif
      write(6,*)ntaylor,' bad ntaylor'
      stop
      end
c
c
      subroutine getquad(qmag1,x0,y0,q,b,rho,zr,nrot)
      implicit real*8 (a-h,o-z)
      dimension npos(5)
      data npos/1,0,-1,0,1/
      real*8 q, b 
      integer nsol
      integer MAX_SOL
      parameter(MAX_SOL = 5)
      real*8 xic(MAX_SOL), yic(MAX_SOL), am(MAX_SOL)
      qmag1 = 0
      if(nrot.eq.0)then
         do 10 i=1,4 
            x = x0 + npos(i)*rho
            y = y0 + npos(i+1)*rho
c            call getpt(qmag,x,y,da,db,ma,mb,phi)
ccc       call
ccc     + getimt(xic, yic, am,
ccc     + nsol, da, db, phi, ma, mb, x, y, qmag)
       call getbinp(qmag,am,xic,yic,x,y,q,b,nsol,zr)

            qmag1 = qmag1 + qmag
 10      continue
         qmag1 = qmag1/4d0
         return
      endif
      if(nrot.eq.1)then
         fac = sqrt(0.5d0)
         do 20 i=-1,1,2
            do 20 j=-1,1,2
               x = x0 + i*rho*fac
               y = y0 + j*rho*fac
c               call getpt(qmag,x,y,da,db,ma,mb,phi)
cc       call
cc     + getimt(xic, yic, am,
cc     + nsol, da, db, phi, ma, mb, x, y, qmag)
       call getbinp(qmag,am,xic,yic,x,y,q,b,nsol,zr)
               qmag1 = qmag1 + qmag
 20      continue
         qmag1 = qmag1/4
         return
      endif
      write(6,*)nrot,' bad nrot'
      stop
      end
