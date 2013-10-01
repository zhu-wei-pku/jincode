c$$$	subroutine getpsi(psi,phi,ecc)
c$$$	implicit real*8 (a-h,o-z)
c$$$	pi = 3.14159265
c$$$	psi= phi
c$$$	do 10 i=1,4
c$$$	   fun = psi - ecc*sin(psi)
c$$$	   dif = phi - fun
c$$$	   der = 1 - ecc*cos(psi)
c$$$	   psi = psi + dif/der
c$$$ 10	continue
c$$$	return
c$$$	end
	subroutine getpsi(psi,phi,ecc)
 	implicit real*8 (a-h,o-z)
 	pi = 3.14159265d0

        if(ecc.gt.0.8)then
           difmin = 1e8

           do 5 npsi = 0,100
              psi = npsi*2*pi/100 + phi
              fun = psi - ecc*sin(psi)
              dif = abs(phi - fun)
              if(dif.lt.difmin)then
                 difmin = dif
                 psimin = psi
              endif
 5         continue

           psi = psimin
        else
           psi= phi
        endif
	psi=phi+ecc*sin(psi) !JS

 	do 10 i=1,10
 	   fun = psi - ecc*sin(psi)
 	   dif = phi - fun
 	   der = 1 - ecc*cos(psi)
 	   psi = psi + dif/der
  10	continue

 	return
 	end
c
c
	subroutine geta(qn,qe,hjd,alpha,delta,t0)
	implicit real*8 (a-h,o-z)
	real*8 sun(3),xpos(3),ypos(3),rad(3),north(3),east(3)
	real*8 spring(3),summer(3)
	data spring/1.,0.,0./
	data summer/0.,0.9174,0.3971/
	data pi/3.14159265/
	ecc = 0.0167
c	vernal = 2719.0
	vernal = 2719.55
c	vernal = 2719.0 + 1000.
	offset = 75
c	offset = 0
	peri   = vernal - offset
	phi = (1 - offset/365.25)*2*pi
c	phi = (1 - 91/365.25)*2*pi
c	phi = 0
	call getpsi(psi,phi,ecc)
c	write(6,*)psi,phi
	costh = (cos(psi) - ecc)/(1-ecc*cos(psi))
	sinth = -sqrt(1-costh**2)
	do 3 i = 1,3
c	   xpos(i) = spring(i)*cos(psi) + summer(i)*sin(psi)
c	   ypos(i) =-spring(i)*sin(psi) + summer(i)*cos(psi)
	   xpos(i) = spring(i)*costh + summer(i)*sinth
	   ypos(i) =-spring(i)*sinth + summer(i)*costh
 3	continue
c	write(6,4)xpos
c	write(6,4)ypos
c	read(5,*)xyz
 4	format(3f10.4)
	north(1) = 0
	north(2) = 0
	north(3) = 1
	radian = 180/3.14159265
c	alpha = (18 +  05./60 + 40.00/3600.)*15
c	delta = -(32 + 56./60 + 08.6/3600.)
c	t0 = 2771
c	t0 = t0-1
	rad(1) = cos(alpha/radian)*cos(delta/radian)
	rad(2) = sin(alpha/radian)*cos(delta/radian)
	rad(3) = sin(delta/radian)
	call cross(east,north,rad)
	call dot(e2,east,east)
	do 5 i=1,3
	   east(i) = east(i)/sqrt(e2)
 5	continue
	call cross(north,rad,east)
 6	format(3f7.3)
c	theta = (t0+1 - peri)/365.25*360.
	phi   = (t0+1 - peri)/365.25*2*pi
	call getpsi(psi,phi,ecc)
	qn2 = 0
	qe2 = 0
	do 10 i=1,3
c	   sun(i) = xpos(i)*cos(theta/radian)+ypos(i)*sin(theta/radian)
	   sun(i) = xpos(i)*(cos(psi)-ecc) +
     *              ypos(i)*sin(psi)*sqrt(1-ecc**2)
	   qn2 = qn2 + sun(i)*north(i)
	   qe2 = qe2 + sun(i)*east(i)
 10	continue
c	theta = (t0-1 - peri)/365.25*360.
	phi   = (t0-1 - peri)/365.25*2*pi
	call getpsi(psi,phi,ecc)
	qn1 = 0
	qe1 = 0
	do 20 i=1,3
c	   sun(i) = xpos(i)*cos(theta/radian)+ypos(i)*sin(theta/radian)
	   sun(i) = xpos(i)*(cos(psi)-ecc) +
     *              ypos(i)*sin(psi)*sqrt(1-ecc**2)
	   qn1 = qn1 + sun(i)*north(i)
	   qe1 = qe1 + sun(i)*east(i)
 20	continue
c	theta = (t0 - peri)/365.25*360.
	phi   = (t0 - peri)/365.25*2*pi
	call getpsi(psi,phi,ecc)
	qn0 = 0
	qe0 = 0
	do 30 i=1,3
c	   sun(i) = xpos(i)*cos(theta/radian)+ypos(i)*sin(theta/radian)
	   sun(i) = xpos(i)*(cos(psi)-ecc) +
     *              ypos(i)*sin(psi)*sqrt(1-ecc**2)
	   qn0 = qn0 + sun(i)*north(i)
	   qe0 = qe0 + sun(i)*east(i)
 30	continue
	vn0 = (qn2-qn1)/2
	ve0 = (qe2-qe1)/2
	factor = 365.25*4.74
c	write(6,*)qn0,qe0,vn0*factor,ve0*factor
c	read(5,*)xyz
	t = hjd
c	theta = (t - peri)/365.25*360.
	phi   = (t - peri)/365.25*2*pi
	call getpsi(psi,phi,ecc)
	qn = -qn0 - vn0*(t-t0)
	qe = -qe0 - ve0*(t-t0)
	do 40 i=1,3
c	   sun(i) = xpos(i)*cos(theta/radian)+ypos(i)*sin(theta/radian)
	   sun(i) = xpos(i)*(cos(psi)-ecc) +
     *              ypos(i)*sin(psi)*sqrt(1-ecc**2)
	   qn = qn + sun(i)*north(i)
	   qe = qe + sun(i)*east(i)
 40	continue
c	write(31,11)njd,qn,qe
 11	format(i6,2f9.5)
 100	continue
	return
	end
c       
	subroutine cross(c,a,b)
	implicit real*8 (a-h,o-z)
	dimension a(3),b(3),c(3)
	c(1) = a(2)*b(3) - b(2)*a(3)
	c(2) = a(3)*b(1) - b(3)*a(1)
	c(3) = a(1)*b(2) - b(1)*a(2)
	return
	end
c
      subroutine dot(c,a,b)
      implicit real*8 (a-h,o-z)
      dimension a(3),b(3)
      c = a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      return
      end

C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.

	subroutine gett(qn,qe,hjd,qlat,qlong,alpha,delta)
	implicit real*8 (a-h,o-z)
	real*8 sun(3),xpos(3),ypos(3),rad(3),north(3),east(3)
	real*8 spring(3),summer(3)
	data pi/3.14159265d0/
	radian = 180/pi
        rearth = 20000/pi/1.5e8 !/100   (remove "/100" to implement)
	vernal = 2719.0
	north(1) = 0
	north(2) = 0
	north(3) = 1
c	alpha = (18 +  05./60 + 41.05/3600.)*15
c	delta = -(28 + 45./60 + 36.2/3600.)
c	t0 = 4233.66
	rad(1) = cos(alpha/radian)*cos(delta/radian)
	rad(2) = sin(alpha/radian)*cos(delta/radian)
	rad(3) = sin(delta/radian)
	call cross(east,north,rad)
	call dot(e2,east,east)
	do 5 i=1,3
	   east(i) = east(i)/sqrt(e2)
 5	continue
	call cross(north,rad,east)
c
	qn = 0
	qe = 0
	qr = 0
        phase = (hjd-vernal)*366.25/365.25 + qlong/360
c        write(6,*)phase
        sun(1) = -cos(phase*2*pi)*cos(qlat/radian)
        sun(2) = -sin(phase*2*pi)*cos(qlat/radian)
        sun(3) = -sin(qlat/radian)
	do 30 i=1,3
	   qn = qn + sun(i)*north(i)*rearth
	   qe = qe + sun(i)*east(i)*rearth
	   qr = qr + sun(i)*rad(i)*rearth
 30	continue
	return
	end
c
c
	subroutine getlatlon(qlat,qlong,nobs)
	implicit real*8 (a-h,o-z)
	dimension qlat(nobs),qlong(nobs)
	do 10 nob=1,nobs
c ---   ogle
	if(nob.eq.0)then    
	   qlong(nob) = -70.702
	   qlat(nob) = -29.0083
	endif
c ---   moa 170 27.9East 43 59.2south  (MOA)
	if(nob.eq.1)then
	   qlong(nob) = +(170 + 27.9/60.)
	   qlat(nob) =  -(43 + 59.2/60.)
	endif
c ---   ctio
	if(nob.eq.2)then
	   qlong(nob) =  -70.815
	   qlat(nob)  =  -30.165
	endif
c ---   TASMANIA
	if(nob.eq.0)then
	   qlong(nob) =  +(147 +26/60. + 21/3600.)
	   qlat(nob)  =  -(42 +48/60. +18/3600.)
	endif
c ---   Mt Lemmon
	if(nob.eq.0)then
	   qlong(nob) =  -110.76
	   qlat(nob)  =  +32.43
	endif
cc ---   farmcove 174:53:37E 36:53:37S  ???
c ---   lt
	if(nob.eq.0)then
c	   qlong(nob) = +(174 + 53/60. + 37/3600.)
c	   qlat(nob) =   -(36 + 53/60. + 43/3600.)
	   qlong(nob) = -(17  + 52/60. + 45/3600.)
	   qlat(nob) =  +(28 + 45/60. + 44.8/3600.)
	endif
c ---   bronberg
	if(nob.eq.3)then
	   qlong(nob) = +(28 + 26/60. + 44/3600.)
	   qlat(nob)  = -(25 + 54/60. + 48/3600.)
	endif
 10	continue
	return
	end
