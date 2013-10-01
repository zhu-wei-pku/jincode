      SUBROUTINE GETMAG_JIN(magbps0, magbfs, xs, ys, b, q, rhos,
     + gamma_lld, npts, errorflag,zr)
c     Packaging for Jin's code. To be expanded to include limb-darkening.
      implicit none
      double precision magbps, magbfs
      double precision magbps0, magbfs0
      double precision magbps1, magbfs1
      double precision magbps2, magbfs2
      double precision b, q
      double precision xs, ys, rhos, gamma_lld
      double precision xs_n !convert to Jin's geometry
      double precision npts

      integer MAX_SOL
      parameter(MAX_SOL = 5)
      complex*16 zr(MAX_SOL)

      integer n
      parameter (n=1) !number of rings used to approximate the limb-darkening
      double precision r(n), s(n), A(n)
c     outer radius of each ring, SB weight, magnification contribution

      integer i

      double precision PI
      parameter (PI=3.141592653589793238462643d0)

      integer errorflag

      xs_n = -xs

      magbfs = 0d0
c      magbps = 0d0
      call getampfb(magbfs0, magbps0, xs_n, ys, b, q, rhos,
     +     rhos, npts, errorflag,zr)     !errorflag = 0 means no error.

      call rings(r, rhos, n)

      call surfb(s,r,n,gamma_lld)

      do i=1, n
         if (i.eq.1) then
            call getampfb(magbfs2, magbps1, xs_n, ys, b, q, r(i)*rhos,
     +           rhos,npts,errorflag,zr)     !errorflag = 0 means no error.
            A(i) = s(i)*magbfs2*r(i)**2d0
         else
            call getampfb(magbfs2, magbps1, xs_n, ys, b, q, r(i)*rhos,
     +           rhos,npts,errorflag,zr)     !errorflag = 0 means no error.
c$$$            call getampfb(magbfs2, magbps2, xs_n, ys, b, q, r(i-1)*rhos,
c$$$     +           errorflag)     !errorflag = 0 means no error.
            A(i) =s(i)*(magbfs2*r(i)**2d0-magbfs1*r(i-1)**2d0)
         endif
         magbfs1 = magbfs2
         magbfs = magbfs+A(i)
      enddo

c      write(6,*) magbps0, magbfs0, magbfs, magbfs0-magbfs

      END

      SUBROUTINE RINGS(r, rhos, n)
c     set the outer fractional radius of each ring
      implicit none
      integer n
      double precision rhos, r(n)
      integer i


      do i=1,n
         r(i) = (1d0-(dble(n-i)/dble(n))**2d0)
      enddo

c      write(6,*) 'r ', r

      END

      SUBROUTINE SURFB(s,r,n,gamma)
c     Determine the weight of each ring based on it's surface brightness.
      implicit none
      integer n
      double precision s(n), r(n), gamma, tots, I0

      integer i

      double precision PI
      parameter (PI=3.141592653589793238462643d0)

      tots = 0d0

c      I0 = 1d0 - gamma*(1d0 - (3.d0/2.d0))

      do i=1, n
c     Constant SB, no limb-darkening. Stupid trial function. Must be replaced.
c         s(i) = 1d0

         if (i.eq.1) then
            s(i) = ( 1.0d0-gamma*( 1d0 + ((1d0-r(i)**2d0)**(3d0/2d0)
     +        -(1d0)**(3d0/2d0)) / (r(i)**2d0) ) )
c            s(i) = 1d0
            s(i) = s(i)
            tots = tots+s(i)*(r(i)**2d0)
         else
            s(i) = ( 1.0d0-gamma*( 1d0 + ((1d0-r(i)**2d0)**(3d0/2d0)
     +           -(1d0-r(i-1)**2d0)**(3d0/2d0)) 
     +           / (r(i)**2d0-r(i-1)**2d0) ) )
c            s(i) = 1d0
            s(i) = s(i)
            tots = tots+s(i)*(r(i)**2d0-r(i-1)**2d0)
         endif
c         s(i) = s(i)/I0
      enddo
c         write(6,*) 's ', gamma, tots!, I0
c         write(6,*) s
c         write(6,*) r
      END

      SUBROUTINE getampfb(magfs,magpt,xsi0,eta0,b,q,rhostar,rho0,
     + npts,err,zr)
      ! Binary lens + finite source effect
CU    uses getim3,cmpexp,area1,finecross,(getmagf)
	! magfs: magnification (finite source)
	! magpt: magnification (point source)
	! (xsi0,eta0): source position 
	! b: binary separation
	! q: mass ratio
	! rhostar: source size/r_E
      double precision PI
      parameter (PI=3.141592653589793238462643d0)
      parameter (nstep = 50000)
      complex*16 zcen,cexp,zetak
      complex*16 z(nstep),zbd(nstep)
      complex*16 zeta(nstep,5)
      complex*16 crsimg(5*nstep)
      double precision magfs,magpt,xsi0,xsi1,eta0,rhostar,rho0,b,q
      double precision mag,theta,xsi,eta,arso,difm,prt,ar2,dif
      double precision xim(5),yim(5),amp(5),thetaadd,am(5)
      complex*16 zr(5)
      integer poscrs(10),posend(10,2)
      integer sols(nstep)
      integer parity(nstep,5)
      integer nsolc,nsol0,nsol,nss
      integer ncrs,posst,posed,ptk,npatch,mstep,mcrs
      integer i,j,k,km,ip,kp,im,in,ia,ib,ipp,iii
      logical check
      logical used(10,5)
c      logical err
      integer err
      double precision npts

c      err = .false.
      err=0.

      arso = PI*rhostar**2
      magfs = 0.d0
      zcen = cmplx(xsi0,eta0)
      call getim3(xim,yim,amp,nsolc,b,q,xsi0,eta0,magpt,zr)
      iii = 0

      do i=1,5
         am(i) = 0.d0
      enddo

      xsi = xsi0 + rhostar
      eta = eta0
      call getim3(xim,yim,amp,nsol0,b,q,xsi,eta,mag,zr)
      ncrs = 0
      check = .false.
      theta = 0.d0
      mstep = 0
      do 1 i=1,nstep
         if (npts.eq.0)then
            thetaadd = 2d0*PI/(1.5d2*rhostar/rho0)
         else
            thetaadd = 2.d0*PI/npts
         endif
         theta = theta + thetaadd
         if(theta.ge.2.d0*PI)then
            check = .true.
            theta = 0.d0
         endif
         mstep = mstep + 1
         call cmpexp(cexp,theta)
         z(mstep) = zcen + rhostar*cexp
         xsi = dreal(z(mstep))
         eta = dimag(z(mstep))
         call getim3(xim,yim,amp,nsol,b,q,xsi,eta,mag,zr)
         k = 0
         do 2 j=1,nsol
            if(amp(j).ne.0.)then
               k = k + 1
               zeta(mstep,k) = cmplx(xim(j),yim(j))
               if(amp(j).gt.0)then
                  parity(mstep,k) = +1
               else
                  parity(mstep,k) = -1
               endif
            endif
 2       continue
         if(k.ne.nsol)then
c            write(6,*)'getim3 fails discovered by getmagf'
c            stop
         endif
         sols(mstep) = nsol
         if(nsol.ne.nsol0)then
            ncrs = ncrs + 1
            call finecross(theta,thetaadd,mstep,z,zeta,parity,
     *      sols,b,q,zcen,rhostar,nsol0,mcrs,nstep,zr)
            poscrs(ncrs) = mcrs
         endif
         nsol0 = sols(mstep)
         if(check)go to 333
 1    continue
c      write(6,*)'not enough count'
      call getmagf(magfs,magpt,xsi0,eta0,rhostar,b,q,zr)
      return
 333  continue
      if(mstep.gt.nstep)then
c         write(6,*)'too much count'
         call getmagf(magfs,magpt,xsi0,eta0,rhostar,b,q,zr)
         return
      endif
      check = .false.
      
      if(abs(1.*ncrs/2.-aint(ncrs/2.)).gt.0.3
     *                 .or.ncrs.gt.10)then
c         write(6,*)' Wrong crossing points',ncrs
c         stop
      endif

      call area1(ar2,z,mstep,nstep)

      if(ncrs.eq.0)then
         posst = 1
      else
         do 8 i=1,10
            do 9 j=1,5
               used(i,j) = .true.
 9          continue
 8       continue
         if(sols(poscrs(1)).eq.3)then
            posst = poscrs(1)
         else
            posst = poscrs(2)
         endif
         posend(1,1) = posst
         if(posst.eq.1)then
            posend(ncrs,2) = mstep
         else
            posend(ncrs,2) = posst - 1
         endif
         do 24 j=1,3
            used(1,j) = .false.
 24       continue
      endif

      in = 1
      do 3 i=1,mstep-1
         ip = posst + i
         im = ip - 1
         if(ip.gt.mstep)ip=ip-mstep
         if(im.gt.mstep)im=im-mstep
         if(ncrs.ne.0)then
         do 25 j=1,ncrs
            if(poscrs(j).eq.ip)then
               posend(in,2) = im
               in = in + 1
               if(in.gt.ncrs)then
c                write(6,*)'crossing points counting error in getmagff'
c                stop
               endif
               posend(in,1) = ip
               do 10 k=1,sols(ip)
                  used(in,k) = .false.
 10            continue
            endif
 25      continue
         endif
         do 4 j = 1,min(sols(im),sols(ip))
            difm = 1.d2
            km = 0
            do 5 k = j,sols(ip)
               if(parity(ip,k)*parity(im,j).ne.1)go to 5
               dif = abs(zeta(ip,k) - zeta(im,j))
               if(dif.lt.difm)then
                  difm = dif
                  km = k
               endif
 5          continue
            if(km.eq.0)go to 4
            zetak = zeta(ip,km)
            ptk = parity(ip,km)
            zeta(ip,km) = zeta(ip,j)
            parity(ip,km) = parity(ip,j)
            zeta(ip,j) = zetak
            parity(ip,j) = ptk
 4       continue
 3    continue

      if(ncrs.eq.0)then

         do 21 i=1,sols(posst)
            difm = abs(zeta(mstep,i) - zeta(posst,i))
            do 22 j=1,sols(posst)
               if(j.eq.i)go to 22
               dif = abs(zeta(mstep,i) - zeta(posst,j))
               if(dif.lt.difm)then
c                  write(6,*)'wrong image topology A'
c                  magfs = magpt*1.d2
                  magfs = magpt
c                  err = .true.
                  err=1.
                  return
               endif
 22         continue
 21      continue

         do 6 i=1,sols(posst)
            do 7 j=1,mstep
               zbd(j) = zeta(j,i)
 7          continue
            call area1(prt,zbd,mstep,nstep)
            iii = iii + 1
            am(iii) = prt*parity(posst,i)/ar2
            magfs = magfs + prt*parity(posst,i)
 6       continue

         magfs = magfs/ar2

         return

      endif

 150  continue

      do 1000 ia = 1,ncrs
      do 1000 ib = 1,5 
      if(used(ia,ib))go to 1000

      npatch = 0
      ptk = parity(posend(ia,1),ib)
      nss = sols(posend(ia,1))
      j = (3 - ptk)/2
      posst = posend(ia,j)
      posed = posend(ia,3-j)
      ip = ia
      kp = ib
      i = posst

 32   continue
      npatch = npatch + 1
 40   continue
      if(i.gt.mstep)then
         i = i - mstep
         go to 40
      endif
      if(i.le.0)then
         i = i + mstep
         go to 40
      endif
      crsimg(npatch) = zeta(i,kp)
      if(i.ne.posed)then
         i = i + ptk
         go to 32
      endif

 36   continue
      difm = 1.d0
      ipp = ip + ptk
      if(ipp.gt.ncrs)ipp=ipp-ncrs
      if(ipp.le.0)ipp=ipp+ncrs
      do 34 k=1,5
         if(used(ipp,k))go to 34
         if(ptk*parity(posend(ipp,1),k).lt.0)go to 34
         j = (3 - ptk)/2
         posst = posend(ipp,j)
         dif = abs(crsimg(npatch)-zeta(posst,k))
         if(dif.lt.difm)then
            difm = dif
            im = ipp
            km = k
         endif
 34   continue
      if(nss.eq.5)then
         do 33 k=1,5
            if(k.eq.kp)go to 33
            if(used(ip,k))go to 33
            if(ptk*parity(posend(ip,1),k).gt.0)go to 33
            j = (3 + ptk)/2
            posst = posend(ip,j)
            dif = abs(crsimg(npatch)-zeta(posst,k))
            if(dif.lt.difm)then
               difm = dif
               im = ip
               km = k
            endif
 33      continue
      endif
      
      if(difm.ge.1.d0)then
c         write(6,*)'cannot find the connected patch in getmagff'
c          magfs = 1.d2*magpt
          magfs = magpt
c          err = .true.
          err = 2.
          return 
c         stop
      endif

      used(im,km) = .true.
      if(im.eq.ia.and.km.eq.ib)then
         if(npatch.gt.3)then
            call area1(prt,crsimg,npatch,5*nstep)
            iii = iii + 1
            am(iii) = prt/ar2
            magfs = magfs + prt
         endif
         go to 1000
      endif

      if(im.eq.ip)then
         ptk = -ptk
      endif
      nss = sols(posend(im,1))
      j = (3 - ptk)/2
      posst = posend(im,j)
      posed = posend(im,3-j)
      i = posst
      ip = im
      kp = km
      go to 40

1000  continue

      magfs = magfs/ar2

      return

      end
*
***   funecross
*
      SUBROUTINE finecross(theta,thetaadd,mstep,z,zeta,parity,
     *           sols,b,q,zcen,rhostar,nsol0,mcrs,nstep,zr)
      parameter(nbi=5)
      complex*16 zcen,cexp
      complex*16 z(nstep),zeta(nstep,5)
      complex*16 z2(nbi),zeta2(nbi,5)
      integer parity(nstep,5),sols(nstep)
      integer parity2(nbi,5),sols2(nbi)
      double precision theta,thetaadd,b,q,xsi,eta
      double precision rhostar,mag,thetau,thetad,mthe
      double precision xi(5),yi(5),am(5),thetain(nbi)
      complex*16 zr(5)
      integer mstep,nsol0,mcrs,nsol,im
      logical nused(nbi)

      z(mstep+nbi) = z(mstep)
      sols(mstep+nbi) = sols(mstep)
      do i=1,sols(mstep+nbi)
         zeta(mstep+nbi,i) = zeta(mstep,i)
         parity(mstep+nbi,i) = parity(mstep,i)
      enddo

      thetau = theta
      thetad = theta - thetaadd

      do 1 i=1,nbi
      thetain(i) = (thetau + thetad)/2.d0
      nused(i) = .true.
      call cmpexp(cexp,thetain(i))
      z2(i) = zcen + rhostar*cexp
      xsi = dreal(z2(i))
      eta = dimag(z2(i))

      call getim3(xi,yi,am,nsol,b,q,xsi,eta,mag,zr)

      k = 0
      do 2 j=1,nsol
         if(am(j).ne.0)then
            k = k + 1
            zeta2(i,k) = cmplx(xi(j),yi(j))
            if(am(j).gt.0)then
               parity2(i,k) = +1
            else
               parity2(i,k) = -1
            endif
         endif
 2    continue
      if(k.ne.nsol)then
c         write(6,*)'getim3 fails discovered by finecross',k,nsol
c         stop
      endif
      sols2(i) = nsol
      if(nsol.ne.nsol0)then
         thetau = thetain(i)
      else
         thetad = thetain(i)
      endif
 1    continue

      do 4 j=0,nbi-1
         mthe = 10.
         do 3 i=1,nbi
            if(nused(i))then
               if(thetain(i).le.mthe)then
                  mthe = thetain(i)
                  im = i
               endif
            endif
 3       continue
         nused(im) = .false.
         z(mstep+j) = z2(im)
         sols(mstep+j) = sols2(im)
         do 5 i=1,sols(mstep+j)
            zeta(mstep+j,i) = zeta2(im,i)
            parity(mstep+j,i) = parity2(im,i)
 5       continue
 4    continue

      do i=0,nbi
         if(nsol0.ne.sols(mstep+i))then
            mcrs = mstep + i
            mstep = mstep + nbi
            return
         endif
      enddo

c      stop
      end
*
***    area1
*
      subroutine area1(ar,bound,num,nmt)
CU    uses conj_jin
      double precision tol
      parameter (tol=1.d-9)
      complex*16 bound(nmt)
      complex*16 local(num)
      complex*16 ztop,zbot,z,zb,delz,delzb,avrz,avrzb
      complex*16 integ
      complex*16 sixi
      double precision ar,ai
      integer i,j,k

      sixi = cmplx(0.d0,6.d0)
      integ = cmplx(0.d0,0.d0)

      do 1 i=1,num
         local(i) = bound(i)
 1    continue

      do 2 i=1,num
         j = i + 1
         if(j.gt.num)j=j-num
         k = i - 1
         if(k.lt.1)k=k+num
         if(i.gt.num)then
         z = local(i-num)
         else
         z = local(i)
         endif
         call conj_jin(zb,z)
         zbot = local(k)
         ztop = local(j)
         delz = (ztop - zbot) / 2.d0
         call conj_jin(delzb,delz)
         avrz = (ztop + zbot) / 2.d0
         call conj_jin(avrzb,avrz)
         integ = integ+(2.d0*zb+avrzb)*delz+2.d0*(avrz-z)*delzb
 2    continue

      integ = integ / sixi
      ar = dreal(integ)
      ai = dimag(integ)

      if(abs(ai).gt.tol)then
c         write(6,*)'cannot find image area'
c         stop
      endif

      return
      end
*
***   cmpexp
*
      subroutine cmpexp(cexp,theta)
      complex*16 cexp
      double precision theta
      cexp = cmplx(dcos(theta),dsin(theta))
      return
      end
*
*
      subroutine conj_jin(zb,z)
      complex*16 z,zb
      double precision a,b
      a = dreal(z)
      b = dimag(z)
      zb = cmplx(a,-b)
      return
      end
*
*** lenseq_jin
*
      subroutine lenseq_jin(zetap,z,z1,z2,m1,m2)
      double precision m1,m2
      complex*16 zetap,z,z1,z2,z1b,z2b,zb
      call conj_jin(z1b,z1)
      call conj_jin(z2b,z2)
      call conj_jin(zb,z)
      zetap = z+(m1/(z1b-zb))+(m2/(z2b-zb))
      return
      end
*
*
      subroutine solchk(derr,z,c,m)
      complex*16 z,c(m+1),err
      double precision derr
      err = c(m+1)
      do n=m,1,-1
         err = err*z + c(n)  
      enddo
      derr = abs(err)
      return
      end
*
***  getim0
*
      subroutine getim0(xi,yi,am,nsol,b,q,xsi,eta,magtot,z)
CU    uses zroots_jin,laguer_jin,lenseq_jin,dmag_jin,conj_jin
      double precision tol
      double precision b,q
      double precision xsi,eta
      double precision xi(5),yi(5),am(5)
      integer nsol
      double precision magtot
      double precision x1,m1,m2,mtot,dm,cm,mtotsq,dmsq,md,cd
      complex*16 z1,z2,z1sq,z1cu,z1fo,z1b,z2b
      complex*16 zeta,zetab,zetabsq
      complex*16 c(6),z(5),zp(5)
      logical polish
      logical polish_only, first_n_minus_2_roots_order_changed
      integer m
      integer loop,loop1,loop2,mloop
      complex*16 dumz,dumzb,dumzp,err,zetap,dzeta,dzetab
      double precision zero,detJ,derr
      integer soln(5)

c
c-----geometry configuration (dreal coordinate)
      x1 = b/2.d0
      m1 = 1.d0/(1.d0+q)
      m2 = q/(1.d0+q)
      dm = (m2-m1)/2.d0
      mtot = (m1+m2)/2.d0
      cm = -x1*(dm/mtot)
      zero = 0.d0
c     
c-----geometry configuration (complex plane)
      z1 = cmplx(x1,zero)
      z2 = -z1
      z1b = z1
      z2b = z2
      z1sq = z1*z1
      z1cu = z1sq*z1
      z1fo = z1cu*z1
c
      zeta = cmplx(xsi,eta)
      zetab =cmplx(xsi,-eta)
      zetabsq = zetab*zetab
c
      mtotsq = mtot*mtot
      dmsq = dm*dm
c
c-----calculate coefficients and solutions
      c(6) = z1sq - zetabsq
      c(5) = zeta*zetabsq-zeta*z1sq-2.d0*mtot*zetab-2.d0*dm*z1
      c(4) = 4.d0*mtot*zeta*zetab + 4.d0*dm*zetab*z1
     *     + 2.d0*zetabsq*z1sq - 2.d0*z1fo
      c(3) = 4.d0*mtotsq*zeta+4.d0*mtot*dm*z1-4.d0*dm*zeta*zetab*z1
     *     + 4.d0*dm*z1cu + 2.d0*zeta*z1fo - 2.d0*zeta*zetabsq*z1sq
      c(2) = z1fo*z1sq - z1fo*zetabsq - 4.d0*mtot*zeta*zetab*z1sq
     *     - 4.d0*dm*zetab*z1cu-4.d0*mtotsq*z1sq-4.d0*dmsq*z1sq
     *     - 8.d0*mtot*dm*zeta*z1
      c(1) = z1sq*(4.d0*dmsq*zeta + 4.d0*mtot*dm*z1 - zeta*z1fo
     *     + 4.d0*dm*zeta*zetab*z1 + 2.d0*mtot*zetab*z1sq
     *     + zeta*zetabsq*z1sq - 2.d0*dm*z1cu )
c
      m = 5
      polish = .true.
      if (z(1).eq.0d0) then
         polish_only = .false.
      else
         polish_only=.true.
      endif
c      call zroots_jin(c,m,z,polish)
      call cmplx_roots_n
     * (z, first_n_minus_2_roots_order_changed, c, polish_only)
c-----test to see if roots are correct
      nsol = 0
      do 5 loop=1,5
         dumz = z(loop)
         dumzp = zetab + m1/(dumz-z1) + m2/(dumz-z2)
         call conj_jin(dumzb,dumzp)
         zp(loop) = dumzb
 5    continue
      do 10 loop1=1,5
         md = 1.d1
         do 11 loop2=1,5
            cd = abs(z(loop1)-zp(loop2))
            if(cd.ge.md)go to 11
            md = cd
            mloop = loop2
 11      continue
         if(mloop.ne.loop1)then
            soln(loop1) = 0
         else
            nsol = nsol + 1
            soln(loop1) = 1
         endif
 10   continue
c
      if((nsol.ne.3).and.(nsol.ne.5))then
         write(6,*)'Wrong nsol', nsol
         do 15 loop=1,5
            dumz = z(loop)
            call solchk(tol,dumz,c,m)
c            zetap = lenseq_jin(dumz,z1,z2,m1,m2)
            call lenseq_jin(zetap,dumz,z1,z2,m1,m2)
            err = zeta - zetap
            derr = abs(err)
c            write(6,*)zeta,zetap
c            write(6,*)z(loop),zp(loop)
c            write(6,*)loop,derr,tol
 15      continue
c         write(6,*)'enter new tol'
c         read(5,*)tol
      endif
c
c-----find magnifications
      do 20 loop=1,5
         if(soln(loop).eq.1)then
            dumz = z(loop)
            xi(loop) = dreal(dumz)
            yi(loop) = dimag(dumz)
c            dumzb = conj_jin(dumz)
            call conj_jin(dumzb,dumz)
            dzeta = (mtot-dm)/(z1b-dumzb)**2
     *            + (mtot+dm)/(-z1b-dumzb)**2
c            dzetab = conj_jin(dzeta)
            call conj_jin(dzetab,dzeta)
            detJ = 1.d0 - dreal(dzeta*dzetab)
            am(loop) = 1.d0/detJ
         else
            am(loop) = 0.d0
         endif
 20   continue
c
c-----calculate total magnification
      magtot = 0.d0
      do 25 loop=1,5
         magtot = magtot + abs(am(loop))
 25   continue
c
c-----return
      return

      end
*
***  getim2
*
      subroutine getim2(xi,yi,am,nsol,b,q,xsi,eta,magtot,zr)
CU    uses getim,(or getim0)
      double precision b,q
      double precision xsi,eta
      double precision xi(5),yi(5),am(5)
      complex*16 zr(5)
      double precision magtot
      integer nsol
      integer knt
      knt = 0
 1    continue
      call getim0(xi,yi,am,nsol,b,q,xsi,eta,magtot,zr)
c      call getim(xi,yi,am,nsol,b,q,xsi,eta,magtot)
      if(nsol.lt.3)then
c         write(6,*)'bad nsol',nsol
c         write(6,*)'stop by getim2'
c         stop
      endif
      if(nsol.ne.4)return
      xsi = xsi + 1.d-8
      knt = knt + 1
      if(knt.gt.100)then
c         write(6,*)'stop by getim2'
c         stop
      endif
c      write(6,*)knt
      go to 1
      end
*
***   getim3
*
      subroutine getim3(xii,yii,amm,nsol,b,q,xsi,eta,magtot,zr)
CU    uses getim2
      double precision b,q
      double precision xsi,eta
      double precision xi(5),yi(5),am(5)
      double precision xii(5),yii(5),amm(5)
      double precision magtot,amin,amax
      complex*16 zr(5)
c      double precision xik,yik,amk
      logical used(5)
      integer nsol
      integer i,j,k

      call getim2(xi,yi,am,nsol,b,q,xsi,eta,magtot,zr)

c      write(6,*) zr

      amax = 0.
      do 1 i=1,5
c         write(6,*) abs(am(i)), am(i), amax
         if(abs(am(i)).gt.amax)amax = abs(am(i))
         used(i) = .false.
 1    continue

      do 2 i=1,5
         amin = 10.d0*amax
c         write(6,*) amin
         do 3 j=1,5
            if(used(j))go to 3
            if(abs(am(j)).lt.amin)then
               amin = abs(am(j))
               k = j
            endif
 3       continue
c         xik = xi(i)
c         yik = yi(i)
c         amk = am(i)
c         write(6,*) i, j, k
         xii(i) = xi(k)
         yii(i) = yi(k)
         amm(i) = am(k)
         used(k) = .true.
c         xi(k) = xik
c         yi(k) = yik
c         am(k) = amk
 2    continue

      if(nsol.eq.5)return

C      if(nsol.eq.5)then
C         do 7 i=5,3,-1
C            if(amm(i).gt.0)then
C               xik = xii(5)
C               yik = yii(5)
C               amk = amm(5)
C               xii(5) = xii(i)
C               yii(5) = yii(i)
C               amm(5) = amm(i)
C               xii(i) = xik
C               yii(i) = yik
C               amm(i) = amk
C               return
C            endif
C 7       continue
C         write(6,*)'error by getim2 found in getim3'
C         stop
C      endif

      do 4 i=1,2
         if(amm(i).ne.0.)then
c            write(6,*)'error detected by getim3'
c            write(6,*)xii(i),yii(i),amm(i)
c            stop
         endif
 4    continue

      do 5 i=1,3
         xii(i) = xii(i+2)
         yii(i) = yii(i+2)
         amm(i) = amm(i+2)
 5    continue

      do 6 i=4,5
         xii(i) = 0.d0
         yii(i) = 0.d0
         amm(i) = 0.d0
 6    continue

      return

      end
*
*** getmagf
*
      subroutine getmagf(magfs,magpt,xsi0,eta0,rhostar,b,q,zr)
CU    uses getim3,cmpexp,area1
      double precision PI
      parameter (PI=3.141592653589793238462643d0)
      parameter (nstep = 7500)
      complex*16 zcen,cexp,zetak
      complex*16 z(nstep),zbd(nstep)
      complex*16 zeta(nstep,5)
      complex*16 crsimg(5*nstep)
      double precision magfs,magpt,xsi0,eta0,rhostar,b,q
      double precision mag,theta,xsi,eta,arso,difm,prt,ar2,dif
      double precision xim(5),yim(5),amp(5)
      complex*16 zr(5)
      integer poscrs(10),posend(10,2)
      integer sols(nstep)
      integer parity(nstep,5)
      integer nsolc,nsol0,nsol,nss
      integer ncrs,posst,posed,ptk,npatch
      integer i,j,k,km,ip,kp,im,in,ia,ib,ipp
      logical check
      logical used(10,5)

      arso = PI*rhostar**2
      magfs = 0.d0
      zcen = cmplx(xsi0,eta0)
      call getim3(xim,yim,amp,nsolc,b,q,xsi0,eta0,magpt,zr)

      xsi = xsi0 + rhostar
      eta = eta0
      call getim3(xim,yim,amp,nsol0,b,q,xsi,eta,mag,zr)
      ncrs = 0
      check = .false.
      do 1 i=1,nstep
c         theta = 2.d0*PI*dreal(i)/real(nstep)
         theta = 2.d0*PI*dble(i)/real(nstep)
         call cmpexp(cexp,theta)
         z(i) = zcen + rhostar*cexp
         xsi = dreal(z(i))
         eta = dimag(z(i))
         call getim3(xim,yim,amp,nsol,b,q,xsi,eta,mag,zr)
         k = 0
         do 2 j=1,nsol
            if(amp(j).eq.0.d0)go to 2
            k = k + 1
            zeta(i,k) = cmplx(xim(j),yim(j))
            if(amp(j).gt.0.d0)then
               parity(i,k) = +1
            else
               parity(i,k) = -1
            endif
 2       continue
         if(k.ne.nsol)then
c            write(6,*)'getim3 fails discovered by getmagf'
c            stop
         endif
         sols(i) = nsol
         if(nsol.ne.nsol0)then
            ncrs = ncrs + 1
            poscrs(ncrs) = i
         endif
         nsol0 = sols(i)
 1    continue
      
c      if(abs(1.d0*dreal(ncrs)/2.d0-(real(ncrs)/2.d0))
      if(abs(1.d0*dble(ncrs)/2.d0-(real(ncrs)/2.d0))
     c         .gt.0.3d0.or.ncrs.gt.10)then
c         write(6,*)' Wrong crossing points',ncrs
c         stop
      endif

      call area1(ar2,z,nstep,nstep)

      if(ncrs.eq.0)then
         posst = 1
      else
         do 8 i=1,10
            do 9 j=1,5
               used(i,j) = .true.
 9          continue
 8       continue
         if(sols(poscrs(1)).eq.3)then
            posst = poscrs(1)
         else
            posst = poscrs(2)
         endif
         posend(1,1) = posst
         if(posst.eq.1)then
            posend(ncrs,2) = nstep
         else
            posend(ncrs,2) = posst - 1
         endif
         do 24 j=1,3
            used(1,j) = .false.
 24       continue
      endif

      in = 1
      do 3 i=1,nstep-1
         ip = posst + i
         im = ip - 1
         if(ip.gt.nstep)ip=ip-nstep
         if(im.gt.nstep)im=im-nstep
         if(ncrs.ne.0)then
         do 25 j=1,ncrs
            if(poscrs(j).eq.ip)then
               posend(in,2) = im
               in = in + 1
               if(in.gt.ncrs)then
c                write(6,*)'crossing points counting error in getmagf'
c                stop
               endif
               posend(in,1) = ip
               do 10 k=1,sols(ip)
                  used(in,k) = .false.
 10            continue
            endif
 25      continue
         endif
         do 4 j = 1,min(sols(im),sols(ip))
            difm = 100.d0
            km = 0
            do 5 k = j,sols(ip)
               if(parity(ip,k)*parity(im,j).ne.1)go to 5
               dif = abs(zeta(ip,k) - zeta(im,j))
               if(dif.lt.difm)then
                  difm = dif
                  km = k
               endif
 5          continue
            if(km.eq.0)go to 4
            zetak = zeta(ip,km)
            ptk = parity(ip,km)
            zeta(ip,km) = zeta(ip,j)
            parity(ip,km) = parity(ip,j)
            zeta(ip,j) = zetak
            parity(ip,j) = ptk
 4       continue
 3    continue

      if(ncrs.eq.0)then

         do 21 i=1,sols(posst)
            difm = abs(zeta(nstep,i) - zeta(posst,i))
            do 22 j=1,sols(posst)
               if(j.eq.i)go to 22
               dif = abs(zeta(nstep,i) - zeta(posst,j))
               if(dif.lt.difm)then
c                  write(6,*)'wrong image topology A'
                  magfs = magpt*1.d2
                  return
               endif
 22         continue
 21      continue

         do 6 i=1,sols(posst)
            do 7 j=1,nstep
               zbd(j) = zeta(j,i)
 7          continue
            call area1(prt,zbd,nstep,nstep)
            magfs = magfs + prt*parity(posst,i)
 6       continue

         magfs = magfs/ar2

         return

      endif

 150  continue

      do 1000 ia = 1,ncrs
      do 1000 ib = 1,5 
      if(used(ia,ib))go to 1000

      npatch = 0
      ptk = parity(posend(ia,1),ib)
      nss = sols(posend(ia,1))
      j = (3 - ptk)/2
      posst = posend(ia,j)
      posed = posend(ia,3-j)
      ip = ia
      kp = ib
      i = posst

 32   continue
      npatch = npatch + 1
 40   continue
      if(i.gt.nstep)then
         i=i-nstep
         go to 40
      endif
      if(i.le.0)then
         i=i+nstep
         go to 40
      endif
      crsimg(npatch) = zeta(i,kp)
      if(i.ne.posed)then
         i = i + ptk
         go to 32
      endif

 36   continue
      difm = 1.
      ipp = ip + ptk
      if(ipp.gt.ncrs)ipp=ipp-ncrs
      if(ipp.le.0)ipp=ipp+ncrs
      do 34 k=1,5
         if(used(ipp,k))go to 34
         if(ptk*parity(posend(ipp,1),k).lt.0)go to 34
         j = (3 - ptk)/2
         posst = posend(ipp,j)
         dif = abs(crsimg(npatch)-zeta(posst,k))
         if(dif.lt.difm)then
            difm = dif
            im = ipp
            km = k
         endif
 34   continue
      if(nss.eq.5)then
         do 33 k=1,5
            if(k.eq.kp)go to 33
            if(used(ip,k))go to 33
            if(ptk*parity(posend(ip,1),k).gt.0)go to 33
            j = (3 + ptk)/2
            posst = posend(ip,j)
            dif = abs(crsimg(npatch)-zeta(posst,k))
            if(dif.lt.difm)then
               difm = dif
               im = ip
               km = k
            endif
 33      continue
      endif
      
      if(difm.ge.1.)then
c         write(6,*)'cannot find the connected patch in getmagf'
          magfs = 1.d2*magpt
          return 
c         stop
      endif

      used(im,km) = .true.
      if(im.eq.ia.and.km.eq.ib)then
         if(npatch.gt.3)then
            call area1(prt,crsimg,npatch,3*nstep)
            magfs = magfs + prt
         endif
         go to 1000
      endif

      if(im.eq.ip)then
         ptk = -ptk
      endif
      nss = sols(posend(im,1))
      j = (3 - ptk)/2
      posst = posend(im,j)
      posed = posend(im,3-j)
      i = posst
      ip = im
      kp = km
      go to 40

1000  continue

      magfs = magfs/ar2

      return

      end
*
*** zroots_jin
*
CCCCCC solve the m-th order complex polynomial
CCCCCC a(m+1)*z^m + ... + a(1) = 0
CCCCCC by Laguerre method
      subroutine zroots_jin(a,m,roots,polish)
CU    uses laguer_jin
      integer MAXM
      double precision EPS
      parameter (MAXM=15)
      parameter (EPS=1.d-14)
      integer m
      complex*16 a(m+1),roots(m)
      complex*16 ad(MAXM)
      complex*16 x,b,c
      integer i,j,jj,its
      logical polish

      do 10 j=1,m+1
         ad(j) = a(j)
 10   continue 

      zero=0.d0

      do 12 j=m,1,-1
         x = cmplx(0.,0.)
         call laguer_jin(ad,j,x,its)
         if(abs(dimag(x)).le.2.*EPS**2*abs(dreal(x)))then
            x = cmplx(dreal(x),zero)
         endif
         roots(j) = x
         b = ad(j+1)
         do 11 jj=j,1,-1
            c = ad(jj)
            ad(jj) = b
            b = x*b + c
 11      continue
 12   continue

      do 13 j=1,m+1
         ad(j) = a(j)
 13   continue 

      if(polish)then
         do 14 j=1,m
            call laguer_jin(ad,m,roots(j),its)
 14      continue
      endif

      do 17 j=2,m
         x = roots(j)
         do 15 i=j-1,1,-1
            if(dreal(roots(i)).le.dreal(x))go to 16
            roots(i+1) = roots(i)
 15      continue
         i = 0
 16      roots(i+1) = x
 17   continue

      return

      end
*
**** laguer_jin
*
      subroutine laguer_jin(ad,m,x,its)
      integer MAXM,MAXIT,MR,MT
      double precision EPSS
      parameter (MAXM=15,MR=8,MT=10000,MAXIT=MT*MR)
      parameter (EPSS=1.d-14)
      integer m,its
      complex*16 a(m+1),x
      integer i,iter,j
      double precision abx,abp,abm,err,frac(MR),az,bz
      complex*16 ad(MAXM),dx,x1,b,d,f,g,h,gp,gm,cz
      save frac
      data frac /.5d0,.25d0,.75d0,.13d0,.38d0,.62d0,.88d0,1.d0/

      do 10 i=1,m+1
         a(i) = ad(i)
 10   continue

      do 12 iter=1,MAXIT
         its = iter
         b = a(m+1)
         err = abs(b)
         d = cmplx(0.,0.)
         f = cmplx(0.,0.)
         abx = abs(x)
         do 11 j=m,1,-1
            f = x*f + d
            d = x*d + b
            b = x*b + a(j)
            err = abs(b) + abx*err
 11      continue
         err=EPSS*err
         if(abs(b).le.err)then
            return
         else
            g = d/b
            h = g*g - f/b - f/b
            gp = g + sqrt(m*m*h+g*g-m*h-m*g*g)
            gm = g - sqrt(m*m*h+g*g-m*h-m*g*g)
            abp = abs(gp)
            abm = abs(gm)
            if(abp.lt.abm)gp=gm
            if(max(abp,abm).gt.0.)then
               dx = m/gp
            else
               az=dlog(1.+abx)
c               bz=dreal(iter)
               bz=real(iter)
               cz=cmplx(az,bz)
               dx = exp(cz)
            endif
         endif
         x1 = x - dx
         if(x.eq.x1)return
         if(mod(iter,MT).ne.0)then
            x = x1
         else
            x = x - dx*frac(iter/MT)
         endif
 12   continue
c      write(6,*) 'too many iterations in SUBROUTINE laguer_jin'
      return
      end






        FUNCTION SIND(x)
        implicit    none
        double precision  sind,x,rad,pi
        parameter         (pi=3.141592654)
        rad = (x/180.0)*pi
        sind = sin(rad)
        return
        end


        FUNCTION COSD(x)
        implicit    none
        double precision  cosd,x,rad,pi
        parameter         (pi=3.141592654)
        rad = (x/180.0)*pi
        cosd = cos(rad)
        return
        end


        function sgn(x)
        implicit          none
        double precision  x,sgn
        sgn = x/abs(x)
        return
        end



        ! =============================================================
        SUBROUTINE getampfs(ampfs,amppt,ux0,uy0,rho)
c       magnification of a single lens event with a finite source
c       by using the Stoke theorem
        implicit          none
        integer           i,n
        parameter         (n=100)
        double precision  theta,dtheta,thetamin,thetamax
        double precision  sum,factor,pi,rho,area0
        double precision  ux0,uy0,u,ux,uy,u0,amppt
        double precision  xplus(n),yplus(n),xminus(n),yminus(n)
        double precision  ampplus,ampminus,ampfs,areaplus,areaminus

        pi = 3.141592654

        u0 = sqrt(ux0**2+uy0**2)
        amppt = (u0**2+2.)/(u0*sqrt(u0**2+4.))

        area0 = pi*rho**2

        thetamin = 0.
        thetamax = 2.0*pi
c        dtheta = (thetamax-thetamin)/dreal(n-1)
        dtheta = (thetamax-thetamin)/dble(n-1)

        theta = thetamin - dtheta
        do 21 i = 1,n
           theta = theta + dtheta
c          --- along the source contour
           ux = ux0 + rho*cos(theta)
           uy = uy0 + rho*sin(theta)
           u = sqrt(ux**2+uy**2)
c          --- image contours
           xplus(i) = 0.5*(u+sqrt(u**2+4.0))*(ux/u)
           yplus(i) = 0.5*(u+sqrt(u**2+4.0))*(uy/u)
           xminus(i) = 0.5*(u-sqrt(u**2+4.0))*(ux/u)
           yminus(i) = 0.5*(u-sqrt(u**2+4.0))*(uy/u)
21      continue
        call STOKE(xplus,yplus,n,areaplus)
        call STOKE(xminus,yminus,n,areaminus)
        ampplus = areaplus/area0
        ampminus = areaminus/area0
        ampfs = ampplus - ampminus
        return
        end



        SUBROUTINE STOKE(x,y,n,area)
        implicit          none
        integer           i,n
        double precision  x(n),y(n),sum,factor,area
        sum = 0.0
        do 22 i = 1,n-1
           factor = x(i)*(y(i+1)-y(i)) - y(i)*(x(i+1)-x(i))
           sum = sum + factor
22      continue
        area = sum/2.
        return
        end

