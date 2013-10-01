      subroutine probtran(ainp, aintrial, adj, nparmin, idum)
      implicit none
      integer nparmin
      real*8 ainp(nparmin), adj(nparmin), aintrial(nparmin)
      integer idum, i
      real RAN2, gasdev
c      external RAN2, gasdev


      do i = 1, nparmin
      aintrial(i) = ainp(i) + gasdev(idum)*adj(i)
      enddo

      end

      subroutine probtran_cov(ainp, aintrial, nparmin, idum,
     + e_vec, e_val)
      implicit none
      integer nparmin
      real*8 ainp(nparmin), adj(nparmin), aintrial(nparmin)
      integer idum, i, j
      real RAN2, gasdev
      real*8 cov(nparmin, nparmin)
      real*8 e_vec(nparmin, nparmin), e_val(nparmin)
      real*8 random(nparmin)
c      external RAN2, gasdev

      do i = 1, nparmin
      random(i) = gasdev(idum)
      enddo

      do i = 1, nparmin
      aintrial(i) = ainp(i)
         do j=1,nparmin
      aintrial(i) = aintrial(i) + random(j)*e_vec(i,j)*e_val(j)
         enddo
      enddo

      end

      FUNCTION RAN2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END


      function get_min(matrix, nchain, last_element, nelement)
      implicit none
      integer nchain
      real*8 get_min, matrix(nchain)
      real*8 chi2min
      integer i, last_element, nelement, first_element

      first_element = last_element - nelement + 1

      chi2min = matrix(first_element)
      do i = first_element,last_element
      if(matrix(i).lt.chi2min) chi2min = matrix(i)
      enddo
      get_min = chi2min

      return
      end

      subroutine  get_cov(ainp_large,chi2_large,
     + eigen_value, eigen_vector, nchain_large, ndim, chi2min
     + ,last_element, nchain)
      implicit none
      integer ndim, nchain
      integer nchain_large
      real*8 ainp_large(ndim,nchain_large), chi2_large(nchain_large)

      real*8 ainp(ndim,nchain), chi2(nchain), chi2tot(nchain)
      real*8 eigen_value(ndim), eigen_vector(ndim, ndim)

      real*8 ex(ndim), pro(ndim,ndim), cov(ndim,ndim)
      real*8 dum(ndim,ndim), sigma(ndim), b(ndim,ndim)
      real*8 prob_sum, chi2min
      integer j1, i1
      logical fix(ndim)
      integer nfix
      integer nrot

      integer last_element, first_element, i_large
      integer i, j, k

      first_element = last_element - nchain + 1

      do i_large=first_element, last_element
      i = i_large - first_element + 1
      chi2(i)    = chi2_large(i_large)
         do j = 1, ndim
      ainp(j, i) = ainp_large(j, i_large)
         enddo
      enddo

306   format(f13.6,f10.6,f10.5,f8.3,f9.6,f9.4, f9.6,
     + 2f9.5, f14.2)

      do i=1,ndim
      eigen_value(i) = 0.d0
      ex(i)    = 0.d0
      do j=1,ndim
      cov(i,j)  = 0.d0
      pro(i,j)   = 0.d0
      eigen_vector(i,j) = 0.d0
      enddo
      enddo

      prob_sum = 0.d0

      do i = 1,nchain
      chi2tot(i)= -(chi2(i) - chi2min)*0.5d0
      enddo

      do i = 1,nchain
      if ( chi2tot(i).ge.-8d0 ) then
ccc      prob_sum = prob_sum + exp(chi2tot(i))
      prob_sum = prob_sum + 1.d0 
      endif
      enddo

      do j = 1,ndim
                do i=1,nchain
      if ( chi2tot(i).ge.-8d0 ) then
ccc          ex(j) = ex(j) + ainp(j,i)*exp(chi2tot(i))/prob_sum
          ex(j) = ex(j) + ainp(j,i)/prob_sum
      endif
                enddo
      enddo

      do j = 1,ndim
             do k = 1,ndim
                do i=1,nchain
      if ( chi2tot(i).ge.-8d0 ) then
ccc          pro(j,k) = pro(j,k) +
ccc     + ainp(j,i)*ainp(k,i)*exp(chi2tot(i))/prob_sum
          pro(j,k) = pro(j,k) +
     + ainp(j,i)*ainp(k,i)/prob_sum
      endif
                enddo
             enddo
          enddo

      do j = 1,ndim
             do k = 1,ndim
          cov(j,k) = pro(j,k) - ex(j)*ex(k)
c          write(6,*) j,k,pro(j,k),ex(j),ex(k),cov(j,k)
             enddo
         enddo

      call jacobi(cov,ndim,ndim,eigen_value,
     + eigen_vector,nrot)

      do j=1,ndim
      if(eigen_value(j).ge.0) then
      eigen_value(j) = sqrt(eigen_value(j))
      else
      eigen_value(j) = 0.0000001d0
      endif
      enddo

      return
      end


      SUBROUTINE jacobi(a,n,np,d,v,nrot)
      implicit none
      INTEGER n,np,nrot,NMAX
      REAL*8 a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=500)
      INTEGER i,ip,iq,j
      REAL*8 c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.
11      continue
        v(ip,ip)=1.
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))
14        continue
15      continue
        if(sm.eq.0.)return
        if(i.lt.4)then
          tresh=0.2*sm/n**2
        else
          tresh=0.
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+
     *g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5*h/a(ip,iq)
                t=1./(abs(theta)+sqrt(1.+theta**2))
                if(theta.lt.0.)t=-t
              endif
              c=1./sqrt(1+t**2)
              s=t*c
              tau=s/(1.+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
23      continue
24    continue
      pause 'too many iterations in jacobi'
      return
      END
      subroutine
     + checkparm(ainp_trial, ainp_min, ainp_max, pass, nparmin)
      integer nparmin, i
      real*8 ainp_trial(nparmin), ainp_min(nparmin), ainp_max(nparmin)
      logical pass
      real*8 te
      pass = .true. 
      do i = 1, nparmin
      if(ainp_trial(i).le.ainp_min(i).or.
     + ainp_trial(i).ge.ainp_max(i)) then
      pass = .false.
      return
      endif
      enddo

      te = ainp_trial(3)/ainp_trial(2)
      if(te.le.50.d0.or.te.ge.300.d0) then
      pass = .false.
      endif

      return
      end

      subroutine
     + checkparm_old(ainp_trial, ainp_min, ainp_max, pass, nparmin)
      integer nparmin, i
      real*8 ainp_trial(nparmin), ainp_min(nparmin), ainp_max(nparmin)
      logical pass
      pass = .true. 
      do i = 1, nparmin
      if(ainp_trial(i).le.ainp_min(i).or.
     + ainp_trial(i).ge.ainp_max(i)) then
      pass = .false.
      return
      endif
      enddo

      return
      end


      FUNCTION gasdev(idum)
      INTEGER idum
      REAL gasdev
CU    USES ran1
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran2
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.*ran2(idum)-1.
        v2=2.*ran2(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END
