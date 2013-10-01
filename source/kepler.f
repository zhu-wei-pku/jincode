       subroutine kepler(rot,a,omega,ecc,tperi,qm,dl,psi,ri,vi,pien,
     *     piee,thetae,pis,tspec)
       implicit real*8 (a-h,o-z)
       real*8 vi(3),ri(3),rot(3,3),el(3),apse(3),bpse(3)
       real*8 v(3),r(3),acr(3)
       real*8 kappa
       data pi/3.1415926535d0/
       data kappa/8.14/
       data year/365.25d0/
c     input units r=einstein radii, v= einstein radii/day
c     calculational units: r=AU, v=v_earth
c     final units r=Einstein radii, v= einstein radii/day
c     distances in kpc, angles in mas
c     tspec is time when velocities vi and positions ri are measured
c     get mass and distance
       pie = sqrt(pien**2 + piee**2)
       qm = thetae/(kappa*pie)
       pirel = pie*thetae
       pil = pis + pirel
       dl = 1/pil
c     convert from input units to calculation units
       do 10 i=1,3
          r(i) = ri(i)*dl*thetae
          v(i) = vi(i)*dl*thetae*year/(2*pi)
  10   continue
c     get energy and angular momentum
       call dot(r2,r,r)
       call dot(v2,v,v)
       energy = v2/2 - qm/sqrt(r2)
       call cross(el,r,v)
       call dot(el2,el,el)
c       write(6,11)r,v,el
 11    format(3f10.4)
c     get semi-major axis, period, eccentrcity
       a = -qm/(2*energy)
       p = sqrt(a**3/qm)
       ecc2 = 1 - el2/(a*qm)
       if(ecc2.lt.0.000001)then
          ecc = 0
       else
          ecc = sqrt(ecc2)
       endif
c     get apse and bpse vectors,
       call cross(apse,v,el)
c       write(6,*)apse
       do 20 i=1,3
          if(ecc.eq.0)then
             apse(i) = r(i)
          else
             apse(i) = apse(i) - qm*r(i)/sqrt(r2)
          endif
  20   continue
c       write(6,*)apse,' apse'
       call cross(bpse,el,apse)
c     normalize orientation vectors
       call dot(apse2,apse,apse)
       call dot(bpse2,bpse,bpse)
       do 30 i = 1,3
          apse(i) = apse(i)/sqrt(apse2)
          bpse(i) = bpse(i)/sqrt(bpse2)
          el(i) = el(i)/sqrt(el2)
  30   continue
c     get psi
       if(ecc.eq.0)then
          theta = 0
          phi = 0
       else
          call cross(acr,apse,r)
          call dot(accr,acr,el)
          call dot(adr,apse,r)
          sinth = accr/sqrt(r2)
          costh = adr/sqrt(r2)
          cospsi = (costh + ecc)/(1 + ecc*costh)
          psi = acos(cospsi)
          if(sinth.lt.0)psi = - psi
       endif
       omega = 2*pi/(p*year)
       a = a/(dl*thetae)
       tperi = tspec - (psi - ecc*sin(psi))/omega
       do 40 i=1,3
          rot(i,1) = apse(i)
          rot(i,2) = bpse(i)
          rot(i,3) = el(i)
 40    continue
       return
       end
c
      subroutine keppos(x1,x2,a,omega,ecc,tperi,rot,t)
      implicit real*8 (a-h,o-z)
      real*8 rot(3,3),rin(3),rout(3),vin(3),vout(3)

      ksw = 0
      if(1.gt.2)ksw = 1

      phi = omega*(t-tperi)
      call getpsi(psi,phi,ecc)
      rin(1) = a*(cos(psi) - ecc)
      rin(2) = a*sin(psi)*sqrt(1-ecc**2)
      rin(3) = 0

      do 10 i=1,3
         rout(i) = 0
         do 9 j=1,3
            rout(i) = rout(i) + rot(i,j)*rin(j)
 9       continue
 10   continue

      x1 = rout(1)
      x2 = rout(2)

      if(ksw.eq.0)return

      dpsidt = omega/(1- ecc*cos(psi))
      vin(1) = -a*sin(psi)*dpsidt
      vin(2) = a*cos(psi)*sqrt(1-ecc**2)*dpsidt
      vin(3) = 0
 
      do 18 i=1,3
         vout(i) = 0
         do 17 j=1,3
            vout(i) = vout(i) + rot(i,j)*vin(j)
 17      continue
 18   continue

      write(6,*)vout
      return
      end
