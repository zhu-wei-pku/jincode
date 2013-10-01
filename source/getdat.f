c     getdat.f
c     Written by Subo Dong/Andy Gould
c     Modified and commented by Jennifer Yee
c     Incorporates code written by Andy Gould for sfit000.f
c
c     This subroutine read the data from disk into vectors and passes it back 
c     to the main program.
c
c     Passed variables:
c     date   - OUT: array of dates for each datapoint
c     flux   - OUT: array of fluxes for each datapoint
c     err    - OUT: errors for each datapoint
c     see    - OUT: seeing for each datapoint
c     bad    - OUT: bad datapoints (flux=0.)
c     ndat   - OUT: number of datapoints/observatory
c     ndata  - IN: maximum number of data points per observatory
c     nobs   - IN: number of observatories
c     tol    - OUT:?
c     kgamma - OUT: filter of each dataset; determines limb-darkening 
c                   coefficients. Unfiltered light --> R
c     qlat, qlong - OUT: latitude and longitude of each observatory
c     alpha, delta - OUT: RA/DEC
c
c     Things you need to directly tinker with:
c     1) character*n direc - 'n'=number of characters in the directory name
c     2) Set the RA and DEC of the event
c     3) direc = Name of the directory where the data are kept
c     4) if using the raw MOA data files, set nobmoa = the dataset number of 
c        the MOA data
c     5) set the names of the data files
c     6) set the error scaling for each data set
c     7) set data points to exclude from the fit
c

      subroutine getdat(date,flux,err,bad,ndat,ndata,nobs,tol,
     +           mgamma,qlat, qlong, alpha, delta)
c      implicit real*8 (a-h,o-z)
      implicit none
      character*21 filename
      character*80 line
      character*3  rootname
      character*1  bandname
      character*1 char
      integer nob, i, ndata, nobs, nform,k
      character*12 nnn
      integer nobmoa
      double precision dat, qi, e, se
      double precision qlat1, qlong1
      logical bad(ndata,nobs)
      logical getlat
      double precision date(ndata,nobs),flux(ndata,nobs),
     +                 err(ndata,nobs),see(ndata, nobs)
      double precision tol(ndata,nobs)
      integer ndat(nobs)
      double precision kgamma(nobs), qlat(nobs), qlong(nobs),
     +     mgamma(nobs)
      double precision alpha, delta, offmoa, hjd
      character*20 imname, str_dat, str_qi, str_e
      character*4 seeing
      real*8 ld(3), gammai, gammav, gammah
      data ld/0.0, 0.0,  0.0/ !
C     #1 TINKER HERE
      character*200 direc 
      integer length
      real*8 date_test,mag_test,err_test,seeing_test,Baseline

      Baseline = 18.707612

c     parameter related to the positions of the observatories
      getlat = .false.
c      if(getlat)call getlatlon(qlat,qlong,nobs)

c     #2 TINKER HERE

c     RA and DEC of the event
      alpha =  (17d0 + 52d0/60d0 + 7.46d0/3600d0)*15d0
      delta = -(29d0 + 50d0/60d0 + 45.72d0/3600d0)

C     #3 TINKER HERE
      direc = '../jincode/data/'
      length=Index(direc, 'data/')
      length=length+4

C     #4 TINKER HERE
      nobmoa = 9
      offmoa = 0.d0 !fixed offset for the MOA data

c     Read in data
      do 10 nob=1,nobs
c         pause
C     #5 TINKER HERE
c     For each observatory, add a line such as:
c        if(nob.eq.1)filename='CT13_OB130341I.bin'
         if(nob.eq.1)filename='curve-4-5451.dat'
         if(nob.eq.2)filename='CT13_OB130341V.pysis'
         if(nob.eq.3)filename='AO_OB130341R.bin'
         if(nob.eq.4)filename='FCO_OB130341U.bin'
         if(nob.eq.5)filename='Tur_OB130341R.bin'
         if(nob.eq.6)filename='PEST_OB130341U.bin'
         if(nob.eq.7)filename='IAC80_OB130341I.bin'
         if(nob.eq.8)filename='phot.bin'
         if(nob.eq.9)filename='moa.phot'
         if(nob.eq.10)filename='Pos_OB130341U.bin'
c         if(nob.eq.4)filename='CT13_OB11V.pho'

         write(6,*) direc(1:length)//filename
         open(1,file=direc(1:length)//filename,status='old')
c         read(1,*) date_test,mag_test,err_test,seeing_test
c         write(*,*) date_test,mag_test,err_test,seeing_test

c     Set the file formatting and record the filter
         nform = 1                  ! OGLE-uFUN
         if(nob.eq.nobmoa)nform = 2 ! MOA
         call getrootband(qlat1,qlong1,bandname,rootname,filename)
         if(.not.getlat)then
            qlat(nob) = qlat1
            qlong(nob)= qlong1
         endif

         if(rootname.eq.'raw')nform = 3 ! PLANET
         if(rootname.eq.'pys')nform = 1 ! PLANET
         kgamma(nob)                    = 1 ! I band
         if(bandname.eq.'V')kgamma(nob) = 2 ! V band
         if(bandname.eq.'H')kgamma(nob) = 3 ! H band
         if(bandname.eq.'R')kgamma(nob) = 4 ! R band
         if(bandname.eq.'U')kgamma(nob) = 4 ! unfiltered = R 
         nform = 0  ! simulated data
         kgamma(nob) = 1

c     Read in the data
         i = 0
 3       continue

         read(1,1,end=4)line
c         write(*,*)line
 1       format(a80)
         read(line,2)char
 2       format(a1)
         if(char.eq.'#')go to 3 !Skip commented lines

c     Compensates for different file formats
         if(nob.eq.nobmoa)then
            read(line,*)str_dat,str_qi,str_e, imname, seeing
            if(seeing.ne.'nan ') then
               read(seeing,*)se
            else
               go to 3
            endif
            if(str_qi.ne.'nan') then
               read(str_dat,*) dat
               read(str_qi,*) qi
               read(str_e,*) e
            else
               go to 3
            endif
            see(i+1,nob)= se - 3d0
         else
            if(nform.eq.0)read(line,*)dat,qi,e
            if(nform.eq.1)read(line,*)dat,qi,e
            if(nform.eq.3)read(line,*)qi,e,dat,se
            if(nform.eq.4)read(line,*)nnn,qi,e,dat,se
            se = 3.0
            see(i+1,nob)=se-1.5d0
         endif


c         write(6,*) filename, dat, qi, e
c     Set errors
C     #6 TINKER HERE
c     For each observatory, add a line such as:
c         if(nob.eq.1)e = e*1.0
c     '1.0' may be replaced with the normalized error.
         if(nob.eq.1)e = e*1.d0
         if(nob.eq.2)e = e*1.d0
         if(nob.eq.3)e = e*1.d0
         if(nob.eq.4)e = e*1.d0

         if(dat.gt.10000d0)dat = dat - 2450000.d0 !abbreviate dates
         !eliminate extraneous data before the start of the current season
c         if(dat.lt.6294.) go to 3 

c     Record the data
         i=i+1
         bad(i,nob) = .false.

         if(nob.eq.nobmoa)then
            flux(i,nob) = qi/6500.d0 + offmoa
            err(i,nob) = e/6500.d0
            call hjdcor(hjd,dat,alpha,delta)
            dat = hjd
         else
            flux(i,nob) = 10.d0**(0.4d0*(Baseline-qi)) !convert mag to flux
            err(i,nob) = e*flux(i,nob)*dlog(10.d0)/2.5d0 !convert errors to flux
         endif
         date(i,nob) = dat

         tol(i,nob) = dlog(10.d0)*e/2.5d0*0.03d0 !convert errors to ?

         if(qi.eq.0.d0)bad(i,nob) = .true.

         if(date(i,nob).lt.0.)bad(i,nob)=.true.
c         if(date(i,nob).lt.5200.)bad(i,nob)=.true.!Restrict extraneous data pts

c     Find null data points
c$$$         if(date(i,nob).lt.1000.)write(6,*) nob,i,' date=',date(i,nob)
c$$$         if(flux(i,nob).lt.0.)write(6,*) nob,i,' flux=', flux(i,nob)
c$$$         if(err(i,nob).lt.0)write(6,*)nob,i,' err=', err(i,nob)
         go to 3

 4       continue
         ndat(nob) = i !set number of data points
c         write(6,*)nob, ndat(nob)
         close(1)
 10   continue

c     #7 TINKER HERE
c     Write bad data statements below, such as:
c     bad(1037,1) = .true. ! 3492.0041

c     Set linear Limb-darkening Parameter:
      gammai=ld(1)
      gammav=ld(2)
      gammah=ld(3)
	do k=1, nobs 
	   if(kgamma(k).eq.1)then
	      mgamma(k) = gammai
	   endif
	   if(kgamma(k).eq.2)then
	      mgamma(k) = gammav
	   endif
	   if(kgamma(k).eq.3)then
	      mgamma(k) = gammah
	   endif
	   if(kgamma(k).eq.4)then
	      mgamma(k) = (gammai+gammav)/2
	   endif
	enddo
      
      return
      end

c-----------------------Supplementary Subroutintes-------------------
c     HJDCOR
c     subroutine written by Andy Gould to correct the HJDs of the MOA data.
	subroutine hjdcor(hjd,date,alpha,delta)
	implicit double precision (a-h,o-z)
	double precision sun(3),xpos(3),ypos(3),rad(3)
	double precision spring(3),summer(3)
	data spring/1.d0,0.d0,0.d0/
	data summer/0.d0,0.9174d0,0.3971d0/
	data pi/3.14159265d0/
	radian = 180d0/pi
	ecc = 0.0167d0
	vernal = 2719.55d0
	offset = 75d0
	peri   = vernal - offset
	phi = (1 - offset/365.25d0)*2d0*pi
	call getpsi(psi,phi,ecc)
	costh = (cos(psi) - ecc)/(1d0-ecc*cos(psi))
	sinth = -sqrt(1-costh**2d0)
	do 3 i = 1,3
	   xpos(i) = spring(i)*costh + summer(i)*sinth
	   ypos(i) =-spring(i)*sinth + summer(i)*costh
 3	continue
	rad(1) = cos(alpha/radian)*cos(delta/radian)
	rad(2) = sin(alpha/radian)*cos(delta/radian)
	rad(3) = sin(delta/radian)
	phi   = (date - peri)/365.25d0*2d0*pi
	call getpsi(psi,phi,ecc)
 7      qr0 = 0
	do 30 i=1,3
	   sun(i) = xpos(i)*(cos(psi)-ecc) +
     *              ypos(i)*sin(psi)*sqrt(1-ecc**2)
	   qr0 = qr0 - sun(i)*rad(i)*8.31d0/1440.d0
 30	continue
	hjd = date + qr0
 100	continue
	return
	end

c     GETROOTBAND
c     Subroutine written by Andy Gould to get the filter and formatting from
c     the data file name. It also returns the location of the observatory.
c
c     Returns:
c     qlat and qlong - the latitude and longitude of the observatory on the 
c                      Earth.
c     bandname - the last letter before the '.' and is the name of the filter.
c     rootname - the file ending after the '.' (first 3 letters) used to 
c                determine the data formatting
      subroutine getrootband(qlat,qlong,bandname,rootname,filename)
      implicit double precision (a-h,o-z)
      character*1 filechar(21)
      character*21 filename
      character*3  rootname
      character*1  bandname
      character*8 word
      logical found
      read(filename,1)filechar
 1    format(21a1)
      do 10 i=20,2,-1
         if(filechar(i).eq.'.')go to 12
 10   continue
      write(6,*)'bad filename ',filename
      stop
 12   continue
      bandname = filechar(i-1)
      write(rootname,13)filechar(i+1),filechar(i+2),filechar(i+3)
 13   format(3a1)
c     For simulated light curves
      word = 'curve'
      length = 5
      call vetword(found,word,length,filechar,0)
      if(found) then
         qlong = 0.0
         qlat  = 0.0
         bandname = 'I'
         return
      endif
c     SURVEY TEAMS
c ---   ogle
      word = 'phot.dat'
      length = 8
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -70.702d0
         qlat  = -29.0083d0
         bandname = 'I'
         return
      endif
      word = 'phot.bin'
      length = 8
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -70.702d0
         qlat  = -29.0083d0
         bandname = 'I'
         return
      endif
      word = 'phot2.dat'
      length = 8
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -70.702d0
         qlat  = -29.0083d0
         bandname = 'I'
         return
      endif
c ---   moa
      word = 'moa'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(170d0 + 27.9d0/60.d0)
         qlat  = -(43d0 + 59.2d0/60.d0)
         bandname = 'I'
         return
      endif
      if(rootname.ne.'pho'.and.rootname.ne.'pys'.and.rootname.ne.'bin')
     +     go to 20
c     uFUN normal
c ---   ctio
      word = 'CT13'
      length = 4
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -70.815d0
         qlat  = -30.165d0
         return
      endif
c ---   Mt Lemmon
      word = 'LOAO'
      length = 4
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -110.76d0
         qlat  = +32.43d0
         return
      endif
c ---  Farmcove
      word = 'FCO'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(174d0 + 53d0/60.d0 + 37d0/3600.d0)
         qlat  =  -(36d0 + 53d0/60.d0 + 43d0/3600.d0)
         return
      endif
c ---  Auckland
      word = 'AO'
      length = 2
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(174d0 + 46d0/60.d0 + 37d0/3600.d0)
         qlat  =  -(36d0 + 54d0/60.d0 + 22d0/3600.d0)
         return
      endif
c ---   CAO
      word = 'CAO'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -(68d0 + 10d0/60.d0 + 49d0/3600.d0)
         qlat  =  -(22d0 +57d0/60.d0 + 10d0/3600.d0)
         return
      endif
c ---   Bronberg
      word = 'Bron'
      length = 4
      call vetword(found,word,length,filechar,0)
 8    continue
      if(found)then
         qlong = +(28d0 + 26d0/60.d0 + 44d0/3600.d0)
         qlat  = -(25d0 + 54d0/60.d0 + 48d0/3600.d0)
         return
      endif
c ---   Hunters Hill
      word = 'HHO'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(149d0 + 06d0/60.d0 + 36d0/3600.d0)
         qlat  = -(35d0 + 09d0/60.d0 + 45d0/3600.d0)
         return
      endif
c --- Klein Karoo
      word = 'KKO'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(21d0 + 40d0/60.d0 + 00d0/3600.d0)
         qlat  = -(33d0 + 32d0/60.d0 + 00d0/3600.d0)
         return
      endif
c ---   Kumeu
      word = 'Kumeu'
      length = 5
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(174d0 + 31d0/60.d0 + 29d0/3600.d0)
         qlat  = -(36d0 + 48d0/60.d0 + 23d0/3600.d0)
         return
      endif
c ---   IAC80
      word = 'IAC80'
      length = 5
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong  = -(16d0 + 30d0/60.d0 + 33.44d0/3600.d0)
         qlat = (28d0 + 17d0/60.d0 + 53.6d0/3600.d0)
         return
      endif
c ---   Hereford
      word = 'HAO'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -(110d0 + 14d0/60.d0 + 15d0/3600.d0)
         qlat  = +(31d0 + 27d0/60.d0 + 08d0/3600.d0)
         return
      endif
c ---   Molehill
      word = 'MAO'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(174d0 + 43d0/60.d0 + 43d0/3600.d0)
         qlat  = -(36d0 + 47d0/60.d0 + 55d0/3600.d0)
         return
      endif
c ---   MDM
      word = 'MDM'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -(111d0 + 36d0/60.d0 + 56d0/3600.d0)
         qlat  = +(31d0 + 57d0/60.d0 + 05d0/3600.d0)
         return
      endif
c ---   OPD - Brazil
      word = 'OPD'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -(45d0 + 34d0/60.d0 + 57d0/3600.d0)
         qlat  = -(22d0 + 32d0/60.d0 + 04d0/3600.d0)
         return
      endif
c ---   Palomar
      word = 'PAL'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -(116d0 + 51d0/60.d0 + 36d0/3600.d0)
         qlat  = +(33d0 + 21d0/60.d0 + 26d0/3600.d0)
         return
      endif
c ---   PEST
      word = 'PEST'
      length = 4
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(115d0 + 47d0/60.d0 + 53d0/3600.d0)
         qlat  = -(31d0 + 59d0/60.d0 + 34d0/3600.d0)
         return
      endif
c ---   Perth
      word = 'Perth'
      length = 5
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(115d0 + 49d0/60.d0 + 00d0/3600.d0)
         qlat  = -(31d0 + 58d0/60.d0 + 00d0/3600.d0)
         return
      endif
c ---   Possum
      word = 'Pos'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(177d0 + 53d0/60.d0 + 29d0/3600.d0)
         qlat  = -(38d0 + 37d0/60.d0 + 26d0/3600.d0)
         return
      endif
c ---   Southern Stars
      word = 'SSO'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -(149d0 + 35d0/60.d0 + 15d0/3600.d0)
         qlat  = -(17d0 + 33d0/60.d0 + 04d0/3600.d0)
         return
      endif
c ---   Turitea
      word = 'Tur'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(175d0 + 39d0/60.d0 + 14d0/3600.d0)
         qlat  = -(40d0 + 24d0/60.d0 + 43d0/3600.d0)
         return
      endif
c ---   Vintage Lane
      word = 'VLO'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(173d0 + 50d0/60.d0 + 21d0/3600.d0)
         qlat  = -(41d0 + 29d0/60.d0 + 30d0/3600.d0)
         return
      endif
c ---   Wise
      word = 'Wise'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(34d0 + 45d0/60.d0 + 44d0/3600.d0)
         qlat  = +(30d0 + 35d0/60.d0 + 50d0/3600.d0)
         return
      endif
      word = 'WC18'
      length = 4
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(34d0 + 45d0/60.d0 + 44d0/3600.d0)
         qlat  = +(30d0 + 35d0/60.d0 + 50d0/3600.d0)
         return
      endif
      go to 120
 20   continue
      if (rootname.eq.'dat') go to 30
c ---  Bronberg (alternate)      
      word = 'bron'
      length = 4
      call vetword(found,word,length,filechar,1)
      if(found)then
         bandname='U'
         go to 8
      endif
      if(rootname.ne.'raw'.and.rootname.ne.'pys')go to 40
c     PLANET
c ---   Tasmania
      word = 'U'
      length = 1
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(147d0 + 26d0/60.d0 + 21d0/3600.d0)
         qlat  = -(42d0 + 48d0/60.d0 + 18d0/3600.d0)
         return
      endif
c ---   SAAO
      word = 'A'
      length = 1
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(16d0 + 12d0/60.d0 + 00d0/3600.d0)
         qlat  = -(23d0 + 18d0/60.d0 + 00d0/3600.d0)
         return
      endif
c ---   Perth
      word = 'W'
      length = 1
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(116d0 + 8d0/60.d0 + 00d0/3600.d0)
         qlat  = -(32d0 + 0d0/60.d0 + 00d0/3600.d0)
         return
      endif
c ---   Danish
      word = 'Z'
      length = 1
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(289d0 + 15d0/60.d0 + 55.457d0/3600.d0)
         qlat  = -(29d0 + 15d0/60.d0 + 15.433d0/3600.d0)
         return
      endif

 30   continue
c     ROBONET
c ---   FTN
      word = 'H'
      length = 1
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -(156d0 + 15d0/60.d0 + 21d0/3600.d0)
         qlat  = +(20d0 + 42d0/60.d0 + 27d0/3600.d0)
         return
      endif
c ---   Liverpool (cararies
      word = 'L'
      length = 1
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -(17d0 + 53d0/60.d0 + 45d0/3600.d0)
         qlat  = +(28d0 + 45d0/60.d0 + 45d0/3600.d0)
         return
      endif
      word = 'J'
      length = 1
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -(17d0 + 53d0/60.d0 + 45d0/3600.d0)
         qlat  = +(28d0 + 45d0/60.d0 + 45d0/3600.d0)
         return
      endif
c     FTS -Robonet
      word = 'I'
      length = 1
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(149d0 + 03d0/60.d0 + 58d0/3600.d0)
         qlat  = -(31d0 + 16d0/60.d0 + 37d0/3600.d0)
         return
      endif
 40   continue
      go to 120
 120  continue
      write(6,121)filechar
      write(*,*) filechar
 121  format(21a1,1x,'lat-long not found')
      read(5,*)xyz
      qlat = 0d0
      qlong= 0d0
      return
      end
c
      subroutine vetword(found,word,length,filechar,iloop)
      character*8 word
      character*1 wordchar(8),filechar(21)
      logical found
      found = .false.
      read(word,1)wordchar
 1    format(8a1)
      do 10 i=0,(21-length)*iloop
         do 5 j=1,length
            if(wordchar(j).ne.filechar(i+j))go to 10
 5       continue
         found = .true.
         return
 10   continue
      return
      end
