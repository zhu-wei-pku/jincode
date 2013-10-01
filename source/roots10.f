!  Copyright 2012 Jan Skowron & Andrew Gould
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!       http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.
!
!-------------------------------------------------------------------!
!
! The authors also release this file under the GNU Lesser General 
! Public License version 2 or any later version, as well as under
! a "customary scientific license", which implies 
! that if this code was important in the scientific process or
! for the results of your scientific work, we ask for the 
! appropriate citation of the Paper (Skowron & Gould 2012).
!
!-------------------------------------------------------------------!
!
!    No    Subroutine
!
!     1   cmplx_roots_gen               - general polynomial solver, works for random degree, not as fast or robust as cmplx_roots_5
!     2   cmplx_roots_5                 - complex roots finding algorithm taylored for 5th order polynomial (with failsafes for polishing)
!     3   sort_5_points_by_separation   - sorting of an array of 5 points, 1st most isolated, 4th and 5th - closest
!     4   sort_5_points_by_separation_i - sorting same as above, returns array of indicies rather than sorted array
!     5   find_2_closest_from_5         - finds closest pair of 5 points
!     6   cmplx_laguerre                - Laguerre's method with simplified Adams' stopping criterion 
!     7   cmplx_newton_spec             - Newton's method with stopping criterion calculated every 10 steps
!     8   cmplx_laguerre2newton         - three regime method: Laguerre's, Second-order General method and Newton's
!     9   solve_quadratic_eq            - quadratic equation solver
!    10   solve_cubic_eq                - cubic equation solver based on Lagrange's method
!    11   divide_poly_1                 - division of the polynomial by (x-p)
!
! fortran 77 code
!
! Paper:  Skowron & Gould 2012  
!         "General Complex Polynomial Root Solver and Its Further Optimization for Binary Microlenses"
!
! for a full text see:
!     http://www.astrouw.edu.pl/~jskowron/cmplx_roots_sg/
!     or http://arxiv.org/find/astro-ph
!     or http://www.adsabs.harvard.edu/abstract_service.html
! see also file NOTICE and LICENSE
!
!-------------------------------------------------------------------!
! _1_                    CMPLX_ROOTS_GEN                            !
!-------------------------------------------------------------------!
      subroutine cmplx_roots_gen(roots, poly, degree, 
     *            polish_roots_after, use_roots_as_starting_points)
      ! This subroutine finds roots of a complex polynomial. 
      ! It is general, however less fast or robust than cmplx_roots_5
      ! which contains failsafe checks in the polishing stage, but is
      ! designed only for 5th order polynomials.
      ! It uses a new dynamic root finding algorithm (see the Paper).
      !
      ! It can use Laguerre's method (subroutine cmplx_laguerre)
      ! or Laguerre->SG->Newton method (subroutine 
      ! cmplx_laguerre2newton - this is default choice) to find 
      ! roots. It divides polynomial one by one by found roots. At the 
      ! end it finds last root from Viete's formula for quadratic 
      ! equation. Finally, it polishes all found roots using a full
      ! polynomial and Newton's or Laguerre's method (default is
      ! Laguerre's - subroutine cmplx_laguerre). 
      ! You can change default choices by commenting out and uncommenting
      ! certain lines in the code below.
      !
      ! Note:
      ! - we solve for the last root with Viete's formula rather 
      !   than doing full Laguerre step (which is time consuming
      !   and unnecessary)
      ! - we do not introduce any preference to real roots
      ! - in Laguerre implementation we omit unneccesarry calculation of
      !   absolute values of denominator
      ! - we do not sort roots. If you need to sort 
      !   roots - we have provided sorting subroutine called:
      !   sort_n_points_by_separation, which sorts points from most 
      !   isolated to most close. Algorithm in this routine can be 
      !   easily used for number of points different than 10.
      !
      implicit none
      ! roots  - array which will hold all roots that had been found. 
      !          If the flag 'use_roots_as_starting_points' is set to 
      !          .true., then instead of point (0,0) we use value from
      !          this array as starting point for cmplx_laguerre
      ! poly -   is an array of polynomial cooefs, length = degree+1, 
      !          poly(1) is a constant term:
      !               1              2             3
      !          poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + ...
      ! degree - degree of the polynomial and size of 'roots' array
      ! polish_roots_after - after all roots have been found by dividing
      !          original polynomial by each root found,
      !          you can opt in to polish all roots using full  
      !          polynomial
      ! use_roots_as_starting_points - usually we start Laguerre's 
      !          method from point (0,0), but you can decide to use the 
      !          values of 'roots' array as starting point for each new
      !          root that is searched for. This is useful if you have
      !          very rough idea where some of the roots can be. 
      !
      integer degree ! intent(in)
      complex*16 poly(degree+1) ! intent(in), coeffs of the polynomial
      complex*16 roots(degree) ! intent(inout)
      logical polish_roots_after, use_roots_as_starting_points ! intent(in)

      complex*16 poly2(degree+1)
      complex*16 zero
      parameter(zero=(0d0,0d0))
      integer i, n, iter
      logical success
      complex*16 coef, prev


      do i=1,degree+1
        poly2(i)=poly(i)
      enddo
      
      ! initialize starting points
      if(.not.use_roots_as_starting_points) then
        do i=1,degree
          roots(i)=zero
        enddo
      endif

      ! skip small degree polynomials from doing Laguerre's method
      if(degree.le.1)then
        if(degree.eq.1) roots(1)=-poly(2)/poly(1)
        return
      endif


      do n=degree, 3, -1

        ! find root with Laguerre's method
        !call cmplx_laguerre(poly2, n, roots(n), iter, success) 
        ! or
        ! find root with (Laguerre's method -> SG method -> Newton's method)
        call cmplx_laguerre2newton(poly2, n, roots(n), 
     *                                     iter, success, 2)
        if(.not.success) then
          roots(n)=zero
          call cmplx_laguerre(poly2, n, roots(n), iter, success)
        endif

        ! divide the polynomial by this root
        coef=poly2(n+1)
        do i=n,1,-1
          prev=poly2(i)
          poly2(i)=coef
          coef=prev+roots(n)*coef
        enddo
        ! coef now holds a remainder - should be close to 0

      enddo

      ! find all but last root with Laguerre's method
      !call cmplx_laguerre(poly2, 2, roots(2), iter, success)
      ! or
      call cmplx_laguerre2newton(poly2, 2, roots(2), 
     *                                   iter, success, 2)
      if(.not.success) then
        call solve_quadratic_eq(roots(2),roots(1),poly2)
      else
        ! calculate last root from Viete's formula
        roots(1)=-(roots(2)+poly2(2)/poly2(3))
      endif



      if(polish_roots_after)then
        do n=1, degree ! ! polish roots one-by-one with a  full polynomial 
          call cmplx_laguerre(poly, degree, roots(n), iter, success) 
          ! or
          !call cmplx_newton_spec(poly, degree, roots(n), iter, success)
        enddo
      endif  

     
      return
      end


!-------------------------------------------------------------------!
! _2_                     CMPLX_ROOTS_n                             !
!-------------------------------------------------------------------!
      subroutine cmplx_roots_n
     * (roots, first_n_minus_2_roots_order_changed, poly, polish_only)
      implicit none
      ! Subroutine finds or polishes roots of a complex polynomial 
      ! (degree=10)
      ! This routine is especially tailored for solving binary lens 
      ! equation in form of 10th order polynomial. 
      !
      ! Use of this routine, in comparission to 'cmplx_roots_gen' can yield
      ! consideribly faster code, because it makes polishing of the roots 
      ! (that come in as a guess from previous solutions) secure by
      ! implementing additional checks on the result of polishing. 
      ! If those checks are not satisfied then routine reverts to the 
      ! robust algorithm. These checks are designed to work for 10th order 
      ! polynomial originated from binary lens equation.
      !
      ! Usage:
      !
      ! polish_only == false - I do not know the roots, routine should  
      !                find them from scratch. At the end it
      !                sorts roots from the most distant to closest.
      !                Two last roots are the closest (in no particular
      !                order).
      ! polish_only = true - I do know the roots pretty well, for example
      !                I have changed the coefficiens of the polynomial 
      !                only a bit, so the two closest roots are 
      !                most likely still the closest ones.
      !                If the output flag 'first_n_minus_2_roots_order_changed'
      !                is returned as 'false', then first 3 returned roots
      !                are in the same order as initialy given to the 
      !                routine. The last two roots are the closest ones, 
      !                but in no specific order (!).
      !                If 'first_n_minus_2_roots_order_changed' is 'true' then
      !                it means that all roots had been resorted.
      !                Two last roots are the closest ones. First is most 
      !                isolated one.
      !
      !
      ! If you do not know the position of the roots just use flag
      ! polish_only=.false. In this case routine will find the roots by
      ! itself.
      
      ! Returns all five roots in the 'roots' array.
      !
      ! poly  - is an array of polynomial cooefs, length = degree+1 
      !       poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + poly(4) x^3 + ...
      ! roots - roots of the polynomial ('out' and optionally 'in')
      !
      !
      ! Jan Skowron 2011
      !
      integer degree
c      parameter(degree=10) ! specifying this speeds up optimization a little
      parameter(degree=5)
      complex*16 roots(degree)                     ! intent(inout) 
      logical first_n_minus_2_roots_order_changed          ! intent(out) 
      complex*16 poly(degree+1)                    ! intent(in)
      logical polish_only                          ! intent(in)
      !------------------------------------------------------------------
      !
      complex*16  remainder, roots_robust(degree)
      real*8 d2min
      integer iter, m, root_n_minus_1, root_n,kk,go_to_robust,i,i2,loops
      complex*16 poly2(degree+1)
      !integer sorted_ind(degree)
      complex*16 zero
      parameter(zero=(0d0,0d0))
      logical succ


      !---------------------------------------
      go_to_robust=0
      if(.not.polish_only) then
        do kk=1,degree
          roots(kk)=zero
        enddo
        go_to_robust=1
      endif
      first_n_minus_2_roots_order_changed=.false.

      do loops=1,3
 111    continue

        ! ROBUST
        ! (we do not know the roots)
        if(go_to_robust>0) then 

          if(go_to_robust>2)then  ! something is wrong
            do kk=1,degree
              roots(kk)=roots_robust(kk)    ! return not-polished roots, because polishing creates errors
            enddo
            return
          endif

          do kk=1,degree+1
            poly2(kk)=poly(kk) ! copy coeffs
          enddo
          do m=degree,4,-1 ! find the roots one-by-one (until 3 are left to be found)
            call cmplx_laguerre2newton(poly2, m, roots(m), 
     *                                         iter, succ, 2)
            if(.not.succ)then
              roots(m)=zero
              call cmplx_laguerre(poly2, m, roots(m), iter, succ)
            endif
            ! divide polynomial by this root
            call divide_poly_1(poly2, remainder, roots(m), poly2, m)
          enddo
          ! find last 3 roots with cubic euqation solver (Lagrange's method)
          call solve_cubic_eq(roots(1),roots(2),roots(3),poly2)
          ! all roots found
          
          ! sort roots - first will be most isolated, last two will be the closest
          call sort_n_points_by_separation(roots) 
          ! copy roots in case something will go wrong during polishing
          do kk=1,degree
            roots_robust(kk)=roots(kk)
          enddo

          ! set flag, that roots have been resorted
          first_n_minus_2_roots_order_changed=.true.
        endif  ! go_to_robust>0

        ! POLISH 
        ! (we know the roots approximately, and we guess that last two are closest)
        !---------------------
          do kk=1,degree+1
            poly2(kk)=poly(kk) ! copy coeffs
          enddo

          do m=1,degree-2
          !do m=1,degree                      ! POWN - polish only with Newton (option)

            ! polish roots with full polynomial
            call cmplx_newton_spec
     *           (poly2, degree, roots(m), iter, succ)

            if(.not.succ)then
              ! go back to robust
              go_to_robust=go_to_robust+1
              do kk=1,degree  ! restart
                roots(kk)=zero
              enddo
              goto 111 ! go back
            endif

          enddo ! m=1,degree-2

          ! comment out division and quadratic if you (POWN) polish with Newton only
          do m=1,degree-2
            call divide_poly_1(poly2, remainder, roots(m), 
     *                        poly2, degree-m+1)
          enddo
          ! last two roots are found with quadratic equation solver 
          ! (this is faster and more robust, although little less accurate)
          call solve_quadratic_eq(roots(degree-1), 
     *                            roots(degree  ),poly2)
          ! all roots found and polished

          ! TEST ORDER
          ! test closest roots if they are the same pair as given to polish
          call find_2_closest_from_n(root_n_minus_1,root_n,d2min, roots)

          ! check if the closest roots are not too close, this could happen
          ! when using polishing with Newton only, when two roots erroneously 
          ! colapsed to the same root. This check is not needed for polishing
          ! 3 roots by Newton and using quadratic  for the remaining two.
          ! If the real roots are so close indeed (very low probability), this will just 
          ! take more time and the unpolished result be returned at the end
          !if(d2min<1d-18) then             ! POWN - polish only with Newton 
          !  go_to_robust=go_to_robust+1    ! POWN - polish only with Newton
          !else                             ! POWN - polish only with Newton

          if((root_n_minus_1.lt.degree-1).or.(root_n.lt.degree-1)) then
            ! after polishing some of the n-2 far roots become one of the 2 closest ones
            ! go back to robust
            if(go_to_robust.gt.0)then
              ! if came from robust 
              ! copy two most isolated roots as starting points for new robust
              do i=1,degree-3
                roots(degree-i+1)=roots_robust(i)
              enddo
            else
              ! came from users initial guess
              ! copy first 2 roots (except the closest ones)
              i2=degree
              do i=1,degree
                if((i.ne.root_n_minus_1).and.(i.ne.root_n))then
                  roots(i2)=roots(i)
                  i2=i2-1
                endif
                if(i2.le.3) exit ! do not copy those that will be done by cubic in robust
              enddo
            endif
            go_to_robust=go_to_robust+1
          else
            ! root_n_minus_1 and root_n comes from the initial last pair
            ! most common case
            return
          endif
          !endif                            ! POWN - polish only with Newton
        !---------------------
      enddo ! loops

      return
      end

!-------------------------------------------------------------------!
! _3_                 SORT_n_POINTS_BY_SEPARATION                   !
!-------------------------------------------------------------------!
      subroutine sort_n_points_by_separation(points)
      ! Sort array of n points 
      ! Most isolated point will become the first point in the array
      ! The closest points will be the last two points in the array
      ! 
      ! Algorithm works well for all dimensions. We put n=10 as 
      ! a hardcoded value just for optimization purposes.
      implicit none
      integer n
      parameter(n=5) !  works for different n as well, but is faster for n as constant (optimization)
!          could make n a subroutine input variable
      complex*16 points(n) ! intent(inout)

      integer sorted_points(n)
      complex*16 savepoints(n)
      integer i

      call sort_n_points_by_separation_i(sorted_points, points) 
      do i=1,n
        savepoints(i)=points(i)
      enddo
      do i=1,n
        points(i)=savepoints(sorted_points(i))
      enddo
      return
      end

!-------------------------------------------------------------------!
! _4_             SORT_n_POINTS_BY_SEPARATION_I                     !
!-------------------------------------------------------------------!
      subroutine sort_n_points_by_separation_i(sorted_points, points)
      ! Return index array that sorts array of n points 
      ! Index of the most isolated point will appear on the first place 
      ! of the output array.
      ! The indices of the closest 2 points will be at the last two 
      ! places in the 'sorted_points' array
      !
      ! Algorithm works well for all dimensions. We put n=10 as 
      ! a hardcoded value just for optimization purposes.
      implicit none
      integer n
      parameter(n=5) !  works for different n as well, but is faster for n as constant (optimization)
!         could make the above "n" a subroutine variable!
      integer sorted_points(n)      ! intent(out)
      complex*16 points(n)          ! intent(in)
       
      real*8 dmin, d1, d2, d
      real*8 distances2(n,n) 
      integer ki, kj, ind2, put
      real*8 neigh1st(n), neigh2nd(n)
      complex*16 p

      do kj=1, n
        distances2(kj,kj)=1d100
      enddo
      dmin=1d100

      do kj=1, n               
        do ki=1, kj-1    
          p=points(ki)-points(kj)
          d=conjg(p)*p
          distances2(ki,kj)=d
          distances2(kj,ki)=d
        enddo
      enddo

      ! find neighbours  
      do kj=1, n
        neigh1st(kj)=1d100
        neigh2nd(kj)=1d100
      enddo
      do kj=1, n
        do ki=1, n
          d=distances2(kj,ki)
          if(d.lt.neigh2nd(kj))then
            if(d.lt.neigh1st(kj))then
              neigh2nd(kj)=neigh1st(kj)
              neigh1st(kj)=d
            else
              neigh2nd(kj)=d
            endif
          endif
        enddo
      enddo    

      ! initialize sorted_points
      do ki=1,n
        sorted_points(ki)=ki
      enddo   
     
      ! sort the rest 1..n-2
      do kj=2,n
        d1=neigh1st(kj)
        d2=neigh2nd(kj)
        put=1
        do ki=kj-1,1,-1
          ind2=sorted_points(ki)
          d=neigh1st(ind2)
          if(d.ge.d1) then
            if(d.eq.d1)then
              if(neigh2nd(ind2)>d2)then
                put=ki+1
                exit
              endif
            else
              put=ki+1
              exit
            endif
          endif
          sorted_points(ki+1)=sorted_points(ki)
        enddo
        sorted_points(put)=kj
      enddo
        
      return
      end

!-------------------------------------------------------------------!
! _5_                     FIND_2_CLOSEST_FROM_n                     !
!-------------------------------------------------------------------!
      subroutine find_2_closest_from_n(i1,i2,d2min, points) 
      ! Returns indices of the two closest points out of array of 10
      implicit none
      integer n
      parameter(n=5) ! will work for other n too, but it is faster with n as constant
!   could make n an input variable
      integer i1,i2        ! intent(out) 
      complex*16 points(n) ! intent(in) 
      real*8 d2min         ! intent(out)
       
      real*8 d2min1, d2
      integer i,j
      complex*16 p

      d2min1=1d100
      do j=1,n
        do i=1,j-1
          p=points(i)-points(j)
          d2=conjg(p)*p
          if(d2.le.d2min1)then
            i1=i
            i2=j
            d2min1=d2
          endif
        enddo
      enddo
      d2min=d2min1
      end


!-------------------------------------------------------------------!
! _6_                     CMPLX_LAGUERRE                            !
!-------------------------------------------------------------------!
      subroutine cmplx_laguerre
     *           (poly, degree, root, iter, success)
      implicit none
      ! Subroutine finds one root of a complex polynomial using 
      ! Laguerre's method. In every loop it calculates simplified
      ! Adams' stopping criterion for the value of the polynomial.
      !
      ! Uses 'root' value as a starting point (!!!!!)
      ! Remember to initialize 'root' to some initial guess or to 
      ! point (0,0) if you have no prior knowledge.
      !
      ! poly - is an array of polynomial cooefs
      !        length = degree+1, poly(1) is constant 
      !               1              2             3
      !          poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + ...
      ! degree - a degree of the polynomial
      ! root - input: guess for the value of a root
      !        output: a root of the polynomial
      ! iter - number of iterations performed (the number of polynomial
      !        evaluations and stopping criterion evaluation)
      ! success - is false if routine reaches maximum number of iterations
      !
      ! For a summary of the method go to: 
      ! http://en.wikipedia.org/wiki/Laguerre's_method
      !
      integer MAX_ITERS, FRAC_JUMP_EVERY, FRAC_JUMP_LEN
      parameter(MAX_ITERS=200)  ! Laguerre is used as a failsafe 
      ! constants needed to break cycles in the scheme
      parameter(FRAC_JUMP_EVERY=10,FRAC_JUMP_LEN=10)
      real*8 FRAC_JUMPS(FRAC_JUMP_LEN)
      real*8 pi
      parameter(pi = 3.141592653589793d0)
      real*8 faq ! jump length
      real*8 FRAC_ERR
      parameter(FRAC_ERR = 2.0d-15)  ! fractional error for kind=8 (see. Adams 1967 Eqs 9 and 10)

      ! subroutine parameters
      integer degree
      complex*16 poly(degree+1) ! intent(in)
      integer iter              ! intent(out)
      complex*16 root           ! intent(inout)
      logical success           ! intent(out)


      complex*16 p         ! value of polynomial
      complex*16 dp        ! value of 1st derivative
      complex*16 d2p_half  ! value of 2nd derivative 
      integer i, k
      logical  good_to_go
      complex*16 denom, denom_sqrt, dx, newroot 
      real*8 ek, absroot, abs2p
      complex*16 fac_netwon, fac_extra, F_half, c_one_nth
      real*8 one_nth, n_1_nth, two_n_div_n_1
      complex*16 zero, c_one
      parameter(zero=(0d0,0d0),c_one=(1d0,0d0))
      real*8 stopping_crit2

      !---------------------------------------
      data FRAC_JUMPS /0.64109297d0, 0.91577881d0, 0.25921289d0,
     *                 0.50487203d0, 0.08177045d0, 0.13653241d0,
     *                 0.306162d0  , 0.37794326d0, 0.04618805d0,
     *                 0.75132137d0/     ! some random numbers


 
      iter=0
      success=.true.

      ! next if-endif block is an EXTREME failsafe, not usually needed, and thus turned off in this version.
      if(.false.)then ! change false-->true if you would like to use caution about having first coefficient == 0
        if(degree.lt.0) then
          write(*,*) 'Error: cmplx_laguerre: degree<0'
          return
        endif
        if(poly(degree+1).eq.zero) then
          if(degree.eq.0) return
          call cmplx_laguerre_wrapper
     *         (poly, degree-1, root, iter, success)
          return
        endif
        if(degree.le.1)then
          if(degree.eq.0) then  ! we know from previous check than poly(1) not equal zero
            success=.false.
            write(*,*) 'Warning: cmplx_laguerre: '//
     *                 'degree=0 and poly(1)/=0, no roots'
            return
          else
            root=-poly(1)/poly(2)
            return
          endif
        endif
      endif
      !  end EXTREME failsafe        

      good_to_go=.false.
      one_nth=1d0/degree
      n_1_nth=(degree-1d0)*one_nth
      two_n_div_n_1=2d0/n_1_nth
      c_one_nth=dcmplx(one_nth,0d0)


      do i=1,MAX_ITERS
        ! prepare stoping criterion
        ek=abs(poly(degree+1))
        absroot=abs(root)
        ! calculate value of polynomial and its first two derivatives
        p  =poly(degree+1)
        dp =zero
        d2p_half=zero
        do k=degree,1,-1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
          d2p_half=dp + d2p_half*root
          dp =p + dp*root
          p  =poly(k)+p*root    ! b_k
          ! Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
          ! Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
          ! ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
          ! Eq 8.
          ek=absroot*ek+abs(p)
        enddo
        iter=iter+1
        
        abs2p=conjg(p)*p
        if(abs2p.eq.0d0) return
        stopping_crit2=(FRAC_ERR*ek)**2
        if(abs2p.lt.stopping_crit2) then ! (simplified a little Eq. 10 of Adams 1967) 
          ! do additional iteration if we are less than 10x from stopping criterion
          if(abs2p.lt.0.01d0*stopping_crit2) then
            return ! return immediately, because we are at very good place
          else
            good_to_go=.true. ! do one iteration more
          endif
        else
          good_to_go=.false. ! reset if we are outside the zone of the root
        endif
      
        faq=1d0
      
        fac_netwon=p/dp
        fac_extra=d2p_half/dp
        F_half=fac_netwon*fac_extra
        
        denom_sqrt=sqrt(c_one-two_n_div_n_1*F_half)

        !G=dp/p  ! gradient of ln(p)
        !G2=G*G
        !H=G2-2d0*d2p_half/p  ! second derivative of ln(p)
        !denom_sqrt=sqrt( (degree-1)*(degree*H-G2) )

        ! NEXT LINE PROBABLY CAN BE COMMENTED OUT 
        if(dble(denom_sqrt).ge.0d0)then ! use realpart or dble in g77 or ifort
          ! real part of a square root is positive for probably all compilers. You can 
          ! test this on your compiler and if so, you can omit this check
          denom=c_one_nth+n_1_nth*denom_sqrt
        else
          denom=c_one_nth-n_1_nth*denom_sqrt
        endif
        if(denom.eq.zero)then !test if demoninators are > 0.0 not to divide by zero
          dx=(absroot+1d0)*exp(dcmplx(0d0,
     *       FRAC_JUMPS(mod(i,FRAC_JUMP_LEN)+1)*2*pi)) ! make some random jump
        else
          dx=fac_netwon/denom
          !dx=degree/denom
        endif
        
        newroot=root-dx
        if(newroot.eq.root) return ! nothing changes -> return
        if(good_to_go)then         ! this was jump already after stopping criterion was met
          root=newroot
          return
        endif

        if(mod(i,FRAC_JUMP_EVERY).eq.0) then ! decide whether to do a jump of modified length (to break cycles)
          faq=FRAC_JUMPS(mod(i/FRAC_JUMP_EVERY-1,FRAC_JUMP_LEN)+1)
          !write(*,'(3i5,f11.6)') i, i/FRAC_JUMP_EVERY, mod(i/FRAC_JUMP_EVERY-1,FRAC_JUMP_LEN)+1, faq
          newroot=root-faq*dx ! do jump of some semi-random length (0<faq<1)
        endif
        root=newroot
      enddo
      success=.false.
      ! too many iterations here  
      end


!-------------------------------------------------------------------!
! _7_                     CMPLX_NEWTON_SPEC                         !
!-------------------------------------------------------------------!
      subroutine cmplx_newton_spec(poly, degree, root, iter, success)
      implicit none
      ! Subroutine finds one root of a complex polynomial using 
      ! Newton's method. It calculates simplified Adams' stopping 
      ! criterion for the value of the polynomial once per 10 iterations,
      ! after initial iteration. This is done to speed up calculations
      ! when polishing roots that are known preety well, and stopping
      ! criterion does not significantly change in their neighborhood.
      !
      ! Uses 'root' value as a starting point (!!!!!)
      ! Remember to initialize 'root' to some initial guess.
      ! Do not initilize 'root' to point (0,0) if the polynomial 
      ! coefficients are strictly real, because it will make going 
      ! to imaginary roots impossible.
      !
      ! poly - is an array of polynomial cooefs
      !        length = degree+1, poly(1) is constant 
      !               1              2             3
      !          poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + ...
      ! degree - a degree of the polynomial
      ! root - input: guess for the value of a root
      !        output: a root of the polynomial
      ! iter - number of iterations performed (the number of polynomial
      !        evaluations)
      ! success - is false if routine reaches maximum number of iterations
      !
      ! For a summary of the method go to: 
      ! http://en.wikipedia.org/wiki/Newton's_method
      !
      integer MAX_ITERS, FRAC_JUMP_EVERY, FRAC_JUMP_LEN
      parameter(MAX_ITERS=50)
      ! constants needed to break cycles in the scheme
      parameter(FRAC_JUMP_EVERY=10,FRAC_JUMP_LEN=10)
      real*8 FRAC_JUMPS(FRAC_JUMP_LEN)
      real*8 pi
      parameter(pi = 3.141592653589793d0)
      real*8 faq ! jump length
      real*8 FRAC_ERR
      parameter(FRAC_ERR = 2.0d-15)  ! fractional error for kind=8 (see. Adams 1967 Eqs 9 and 10)

      ! subroutine parameters
      integer degree
      complex*16 poly(degree+1) ! intent(in)
      integer iter              ! intent(out)
      complex*16 root           ! intent(inout)
      logical success           ! intent(out)


      complex*16 p    ! value of polynomial
      complex*16 dp   ! value of 1st derivative 
      integer i, k
      logical  good_to_go
      complex*16 dx, newroot
      real*8 ek, absroot, abs2p
      complex*16 zero
      parameter(zero=(0d0,0d0))
      real*8 stopping_crit2

      !---------------------------------------
      data FRAC_JUMPS /0.64109297d0, 0.91577881d0, 0.25921289d0,  
     *                 0.50487203d0, 0.08177045d0, 0.13653241d0,  
     *                 0.306162d0  , 0.37794326d0, 0.04618805d0,  
     *                 0.75132137d0/     ! some random numbers


      iter=0
      success=.true.

      ! next if-endif block is an EXTREME failsafe, not usually needed, and thus turned off in this version.
      if(.false.)then ! change false-->true if you would like to use caution about having first coefficient == 0
        if(degree.lt.0) then
          write(*,*) 'Error: cmplx_newton_spec: degree<0'
          return
        endif
        if(poly(degree+1).eq.zero) then
          if(degree.eq.0) return
          call cmplx_newton_spec_wrapper
     *         (poly, degree-1, root, iter, success)
          return
        endif
        if(degree.le.1)then
          if(degree.eq.0) then  ! we know from previous check than poly(1) not equal zero
            success=.false.
            write(*,*) 'Warning: cmplx_newton_spec: degree=0 and'//
     *                 ' poly(1)/=0, no roots'
            return
          else
            root=-poly(1)/poly(2)
            return
          endif
        endif
      endif
      ! end EXTREME failsafe
 
      good_to_go=.false.

      do i=1,MAX_ITERS
        faq=1d0

        ! prepare stoping criterion
        ! calculate value of polynomial and its first two derivatives
        p  =poly(degree+1)
        dp =zero

        if(mod(i,10).eq.1) then ! calculate stopping criterion every tenth iteration
          ek=abs(poly(degree+1))
          absroot=abs(root)
          do k=degree,1,-1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
            dp =p + dp*root
            p  =poly(k)+p*root    ! b_k
            ! Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
            ! Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
            ! ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
            ! Eq 8.
            ek=absroot*ek+abs(p)
          enddo
          stopping_crit2=(FRAC_ERR*ek)**2
        else                  ! calculate just the value and derivative
          do k=degree,1,-1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
            dp =p + dp*root
            p  =poly(k)+p*root    ! b_k
          enddo
        endif
        iter=iter+1

        
        abs2p=conjg(p)*p
        if(abs2p.eq.0d0) return

        if(abs2p.lt.stopping_crit2) then ! (simplified a little Eq. 10 of Adams 1967)
          if(dp.eq.zero) return ! if we have problem with zero, but we are close to the root, just accept
          ! do additional iteration if we are less than 10x from stopping criterion
          if(abs2p.lt.0.01d0*stopping_crit2) then
            return ! return immediately, because we are at very good place
          else
            good_to_go=.true. ! do one iteration more
          endif
        else
          good_to_go=.false. ! reset if we are outside the zone of the root
        endif


        if(dp.eq.zero)then
          ! problem with zero
          dx=(abs(root)+1d0)*
     *       exp(dcmplx(0d0,FRAC_JUMPS(mod(i,FRAC_JUMP_LEN)+1)*2*pi)) ! make some random jump
        else
          dx=p/dp  ! Newton method, see http://en.wikipedia.org/wiki/Newton's_method
        endif

      

        newroot=root-dx
        if(newroot.eq.root) return ! nothing changes -> return
        if(good_to_go)then         ! this was jump already after stopping criterion was met
          root=newroot
          return
        endif

        if(mod(i,FRAC_JUMP_EVERY).eq.0) then ! decide whether to do a jump of modified length (to break cycles)
          faq=FRAC_JUMPS(mod(i/FRAC_JUMP_EVERY-1,FRAC_JUMP_LEN)+1)
          newroot=root-faq*dx ! do jump of some semi-random length (0<faq<1)
        endif
        root=newroot
      enddo
      success=.false.
      return
      ! too many iterations here  
      end

!-------------------------------------------------------------------!
! _8_                     CMPLX_LAGUERRE2NEWTON                     !
!-------------------------------------------------------------------!
      subroutine cmplx_laguerre2newton
     *           (poly, degree, root, iter, success, starting_mode)
      implicit none
      ! Subroutine finds one root of a complex polynomial using 
      ! Laguerre's method, Second-order General method and Newton's
      ! method - depending on the value of function F, which is a 
      ! combination of second derivative, first derivative and
      ! value of polynomial [F=-(p"*p)/(p'p')].
      ! 
      ! Subroutine has 3 modes of operation. It starts with mode=2
      ! which is the Laguerre's method, and continues until F 
      ! becames F<0.50, at which point, it switches to mode=1, 
      ! i.e., SG method (see paper). While in the first two
      ! modes, routine calculates stopping criterion once per every
      ! iteration. Switch to the last mode, Newton's method, (mode=0) 
      ! happens when becomes F<0.05. In this mode, routine calculates
      ! stopping criterion only once, at the beginning, under an 
      ! assumption that we are already very close to the root.
      ! If there are more than 10 iterations in Newton's mode, 
      ! it means that in fact we were far from the root, and
      ! routine goes back to Laguerre's method (mode=2).
      !
      ! Uses 'root' value as a starting point (!!!!!)
      ! Remember to initialize 'root' to some initial guess or to 
      ! point (0,0) if you have no prior knowledge.
      !
      ! poly - is an array of polynomial cooefs
      !        length = degree+1, poly(1) is constant 
      !               1              2             3
      !          poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + ...
      ! degree - a degree of the polynomial
      ! root - input: guess for the value of a root
      !        output: a root of the polynomial
      ! iter - number of iterations performed (the number of polynomial
      !        evaluations and stopping criterion evaluation)
      ! success - is false if routine reaches maximum number of iterations
      ! starting_mode - this should be by default = 2. However if you  
      !                 choose to start with SG method put 1 instead. 
      !                 Zero will cause the routine to 
      !                 start with Newton for first 10 iterations, and 
      !                 then go back to mode 2.
      !                 
      !
      ! For a summary of the method see the paper: Skowron & Gould (2012)
      !

      integer MAX_ITERS, FRAC_JUMP_EVERY, FRAC_JUMP_LEN
      parameter(MAX_ITERS=50)
      ! constants needed to break cycles in the scheme
      parameter(FRAC_JUMP_EVERY=10,FRAC_JUMP_LEN=10)
      real*8 FRAC_JUMPS(FRAC_JUMP_LEN)
      real*8 pi
      parameter(pi = 3.141592653589793d0)
      real*8 faq ! jump length
      real*8 FRAC_ERR
      parameter(FRAC_ERR = 2.0d-15)  ! fractional error for kind=8 (see. Adams 1967 Eqs 9 and 10)

      ! subroutine parameters
      integer degree
      complex*16 poly(degree+1) ! intent(in)
      integer iter              ! intent(out)
      complex*16 root           ! intent(inout)
      logical success           ! intent(out)
      integer starting_mode     ! intent(in)


      complex*16 p         ! value of polynomial
      complex*16 dp        ! value of 1st derivative
      complex*16 d2p_half  ! value of 2nd derivative 
      integer i, j, k
      logical  good_to_go
      complex*16 denom, denom_sqrt, dx, newroot
      real*8 ek, absroot, abs2p, abs2_F_half
      complex*16 fac_netwon, fac_extra, F_half, c_one_nth
      real*8 one_nth, n_1_nth, two_n_div_n_1
      integer mode
      complex*16 zero, c_one
      parameter(zero=(0d0,0d0),c_one=(1d0,0d0))
      real*8 stopping_crit2

      !---------------------------------------
      data FRAC_JUMPS /0.64109297d0, 0.91577881d0, 0.25921289d0,
     *                 0.50487203d0, 0.08177045d0, 0.13653241d0,
     *                 0.306162d0  , 0.37794326d0, 0.04618805d0,
     *                 0.75132137d0/     ! some random numbers

     
      iter=0
      success=.true.

      ! next if-endif block is an EXTREME failsafe, not usually needed, and thus turned off in this version.
      if(.false.)then ! change false-->true if you would like to use caution about having first coefficient == 0
        if(degree.lt.0) then
          write(*,*) 'Error: cmplx_laguerre2newton: degree<0'
          return
        endif
        if(poly(degree+1).eq.zero) then
          if(degree.eq.0) return 
          call cmplx_laguerre2newton_wrapper
     *          (poly, degree-1, root, iter, success, starting_mode)
          return
        endif
        if(degree.le.1)then
          if(degree.eq.0) then  ! we know from previous check than poly(1) not equal zero
            success=.false.
            write(*,*) 'Warning: cmplx_laguerre2newton: '//
     *                 'degree=0 and poly(1)/=0, no roots'
            return
          else
            root=-poly(1)/poly(2)
            return
          endif  
        endif
      endif
      ! end EXTREME failsafe

      j=1
      good_to_go=.false.

      mode=starting_mode  ! mode=2 full laguerre, mode=1 SG, mode=0 newton

      do ! infinite loop, just to be able to come back from newton, if more than 10 iteration there

      !------------------------------------------------------------- mode 2
      if(mode.ge.2) then  ! LAGUERRE'S METHOD
        one_nth=1d0/degree
        n_1_nth=(degree-1d0)*one_nth
        two_n_div_n_1=2d0/n_1_nth
        c_one_nth=dcmplx(one_nth,0d0)
      
        do i=1,MAX_ITERS  !
          faq=1d0

          ! prepare stoping criterion
          ek=abs(poly(degree+1))
          absroot=abs(root)
          ! calculate value of polynomial and its first two derivatives
          p  =poly(degree+1)
          dp =zero
          d2p_half=zero
          do k=degree,1,-1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
            d2p_half=dp + d2p_half*root
            dp =p + dp*root
            p  =poly(k)+p*root    ! b_k
            ! Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
            ! Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
            ! ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
            ! Eq 8.
            ek=absroot*ek+abs(p)
          enddo
          abs2p=conjg(p)*p !abs(p)
          iter=iter+1
          if(abs2p.eq.0d0) return
          
          stopping_crit2=(FRAC_ERR*ek)**2
          if(abs2p.lt.stopping_crit2) then ! (simplified a little Eq. 10 of Adams 1967) 
            ! do additional iteration if we are less than 10x from stopping criterion
            if(abs2p.lt.0.01d0*stopping_crit2) then ! ten times better than stopping criterion
              return ! return immediately, because we are at very good place
            else
              good_to_go=.true. ! do one iteration more
            endif
          else
            good_to_go=.false. ! reset if we are outside the zone of the root
          endif
        
        
          fac_netwon=p/dp
          fac_extra=d2p_half/dp
          F_half=fac_netwon*fac_extra

          abs2_F_half=conjg(F_half)*F_half
          if(abs2_F_half.le.0.0625d0)then     ! F<0.50, F/2<0.25
            ! go to SG method
            if(abs2_F_half.le.0.000625d0)then ! F<0.05, F/2<0.025
              mode=0 ! go to Newton's
            else
              mode=1 ! go to SG
            endif
          endif
          
          
          denom_sqrt=sqrt(c_one-two_n_div_n_1*F_half)

          ! NEXT LINE PROBABLY CAN BE COMMENTED OUT
          if(dble(denom_sqrt).ge.0d0)then  ! use realpart or dble for g77 or ifort
            ! real part of a square root is positive for probably all compilers. You can 
            ! test this on your compiler and if so, you can omit this check
            denom=c_one_nth+n_1_nth*denom_sqrt
          else
            denom=c_one_nth-n_1_nth*denom_sqrt
          endif
          if(denom.eq.zero)then !test if demoninators are > 0.0 not to divide by zero
            dx=(abs(root)+1d0)*exp( 
     *          dcmplx(0d0,FRAC_JUMPS(mod(i,FRAC_JUMP_LEN)+1)*2*pi)) ! make some random jump
          else
            dx=fac_netwon/denom
          endif
     
          newroot=root-dx
          if(newroot.eq.root) return ! nothing changes -> return
          if(good_to_go)then         ! this was jump already after stopping criterion was met
            root=newroot
            return
          endif

          if(mode.ne.2) then 
            root=newroot
            j=i+1    ! remember iteration number
            exit     ! go to Newton's or SG 
          endif

          if(mod(i,FRAC_JUMP_EVERY).eq.0) then ! decide whether to do a jump of modified length (to break cycles)
            faq=FRAC_JUMPS(mod(i/FRAC_JUMP_EVERY-1,FRAC_JUMP_LEN)+1)
            newroot=root-faq*dx ! do jump of some semi-random length (0<faq<1)
          endif
          root=newroot
        enddo ! do mode 2

        if(i.ge.MAX_ITERS) then
          success=.false.
          return
        endif

      endif ! if mode 2


      !------------------------------------------------------------- mode 1
      if(mode.eq.1) then  ! SECOND-ORDER GENERAL METHOD (SG)

        do i=j,MAX_ITERS  !
          faq=1d0

          ! calculate value of polynomial and its first two derivatives
          p  =poly(degree+1)
          dp =zero
          d2p_half=zero
          if(mod(i-j,10).eq.0)then
            ! prepare stoping criterion
            ek=abs(poly(degree+1))
            absroot=abs(root)
            do k=degree,1,-1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
              d2p_half=dp + d2p_half*root
              dp =p + dp*root
              p  =poly(k)+p*root    ! b_k
              ! Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
              ! Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
              ! ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
              ! Eq 8.  
              ek=absroot*ek+abs(p)
            enddo
            stopping_crit2=(FRAC_ERR*ek)**2
          else
            do k=degree,1,-1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
              d2p_half=dp + d2p_half*root
              dp =p + dp*root
              p  =poly(k)+p*root    ! b_k
            enddo
          endif


          abs2p=conjg(p)*p !abs(p)**2
          iter=iter+1
          if(abs2p.eq.0d0) return


          if(abs2p.lt.stopping_crit2) then ! (simplified a little Eq. 10 of Adams 1967) 
            if(dp.eq.zero) return
            ! do additional iteration if we are less than 10x from stopping criterion
            if(abs2p.le.0.01d0*stopping_crit2) then ! ten times better than stopping criterion
              return ! return immediately, because we are at very good place
            else
              good_to_go=.true. ! do one iteration more
            endif
          else
            good_to_go=.false. ! reset if we are outside the zone of the root
          endif
          
          if(dp.eq.zero)then !test if demoninators are > 0.0 not to divide by zero
            dx=(abs(root)+1d0)*exp(
     *          dcmplx(0d0,FRAC_JUMPS(mod(i,FRAC_JUMP_LEN)+1)*2*pi)) ! make some random jump
          else
            fac_netwon=p/dp
            fac_extra=d2p_half/dp
            F_half=fac_netwon*fac_extra

            abs2_F_half=conjg(F_half)*F_half
            if(abs2_F_half.le.0.000625d0)then ! F<0.05, F/2<0.025
              mode=0 ! set Newton's, go there after jump
            endif
          
            dx=fac_netwon*(c_one+F_half)  ! SG
          endif
          
          newroot=root-dx
          if(newroot.eq.root) return ! nothing changes -> return
          if(good_to_go)then         ! this was jump already after stopping criterion was met
            root=newroot
            return
          endif

          if(mode.ne.1) then 
            root=newroot
            j=i+1    ! remember iteration number
            exit     ! go to Newton's
          endif

          if(mod(i,FRAC_JUMP_EVERY).eq.0) then ! decide whether to do a jump of modified length (to break cycles)
            faq=FRAC_JUMPS(mod(i/FRAC_JUMP_EVERY-1,FRAC_JUMP_LEN)+1)
            newroot=root-faq*dx ! do jump of some semi-random length (0<faq<1)
          endif
          root=newroot



        enddo ! do mode 1

        if(i.ge.MAX_ITERS) then
          success=.false.
          return
        endif
     
      endif ! if mode 1


      !------------------------------------------------------------- mode 0
      if(mode.eq.0) then  ! newton

        do i=j,j+10  ! do only 10 iterations the most, then go back to full Laguerre's
          faq=1d0

          
          ! calculate value of polynomial and its first two derivatives
          p  =poly(degree+1)
          dp =zero
          if(i.eq.j)then ! calculate stopping crit only once at the begining
            ! prepare stoping criterion
            ek=abs(poly(degree+1))
            absroot=abs(root)
            do k=degree,1,-1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
              dp =p + dp*root
              p  =poly(k)+p*root    ! b_k
              ! Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
              ! Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
              ! ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
              ! Eq 8.  
              ek=absroot*ek+abs(p)
            enddo
            stopping_crit2=(FRAC_ERR*ek)**2
          else
            do k=degree,1,-1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
              dp =p + dp*root
              p  =poly(k)+p*root    ! b_k
            enddo
          endif
          abs2p=conjg(p)*p !abs(p)**2
          iter=iter+1
          if(abs2p.eq.0d0) return


          if(abs2p.lt.stopping_crit2) then ! (simplified a little Eq. 10 of Adams 1967)
            if(dp.eq.zero) return 
            ! do additional iteration if we are less than 10x from stopping criterion
            if(abs2p.lt.0.01d0*stopping_crit2) then ! ten times better than stopping criterion
              return ! return immediately, because we are at very good place
            else
              good_to_go=.true. ! do one iteration more
            endif
          else
            good_to_go=.false. ! reset if we are outside the zone of the root
          endif
        
          if(dp.eq.zero)then ! test if demoninators are > 0.0 not to divide by zero
            dx=(abs(root)+1d0)*exp(
     *          dcmplx(0d0,FRAC_JUMPS(mod(i,FRAC_JUMP_LEN)+1)*2*pi)) ! make some random jump
          else
            dx=p/dp
          endif
          
          newroot=root-dx
          if(newroot.eq.root) return ! nothing changes -> return
          if(good_to_go)then         
            root=newroot
            return
          endif

          ! this loop is done only 10 times. So skip this check
          !if(mod(i,FRAC_JUMP_EVERY).eq.0) then ! decide whether to do a jump of modified length (to break cycles)
          !  faq=FRAC_JUMPS(mod(i/FRAC_JUMP_EVERY-1,FRAC_JUMP_LEN)+1)
          !  newroot=root-faq*dx ! do jump of some semi-random length (0<faq<1)
          !endif
          root=newroot

        enddo ! do mode 0 10 times

        if(iter.ge.MAX_ITERS) then
          ! too many iterations here
          success=.false.
          return
        endif
        mode=2 ! go back to Laguerre's. This happens when we were unable to converge in 10 iterations with Newton's

      endif ! if mode 0

      enddo ! end of infinite loop

      !------------------------------------------------------------- 
      success=.false.
      end


!-------------------------------------------------------------------!
! _9_                     SOLVE_QUADRATIC_EQ                        !
!-------------------------------------------------------------------!
      subroutine solve_quadratic_eq(x0,x1,poly)
      ! Quadratic equation solver for complex polynomial (degree=2)
      implicit none
      complex*16 x0,x1    ! intent(out)
      complex*16 poly(*)  ! intent(in)  ! coeffs of the polynomial
      ! poly - is an array of polynomial cooefs, length = degree+1, poly(1) is constant 
      !             1              2             3
      !        poly(1) x^0 + poly(2) x^1 + poly(3) x^2
      complex*16 a, b, c, b2, delta

      complex*16  val, x
      integer  i
      complex*16 zero
      parameter(zero=(0d0,0d0))


      a=poly(3)
      b=poly(2)
      c=poly(1)
      ! quadratic equation: a z^2 + b z + c = 0  

      b2=b*b
      delta=sqrt(b2-4d0*(a*c))
      if( dble(conjg(b)*delta).ge.0d0 )then  ! scallar product to decide the sign yielding bigger magnitude
        ! use dble() or realpart() instead of real() not to loose precision in ifort or g77
        x0=-0.5d0*(b+delta)
      else
        x0=-0.5d0*(b-delta)
      endif
      if(x0.eq.zero)then
        x1=zero
      else ! Viete's formula
        x1=c/x0
        x0=x0/a
      endif


      if(.false.)then  ! print the results

        x=x0
        val=poly(3)
        do i=2,1,-1
          val=val*x+poly(i)
        enddo
        write(*,'(2f19.15,a3,2f19.15)') x,' ->',val

        x=x1
        val=poly(3)
        do i=2,1,-1
          val=val*x+poly(i)
        enddo
        write(*,'(2f19.15,a3,2f19.15)') x,' ->',val

      endif

      end

!-------------------------------------------------------------------!
! _10_                    SOLVE_CUBIC_EQ                            !
!-------------------------------------------------------------------!
      subroutine solve_cubic_eq(x0,x1,x2,poly)
      ! Cubic equation solver for complex polynomial (degree=3)
      ! http://en.wikipedia.org/wiki/Cubic_function   Lagrange's method
      implicit none
      complex*16 x0, x1, x2   ! intent(out)
      complex*16 poly(*)      ! intent(in)  ! coeffs of the polynomial
      ! poly - is an array of polynomial cooefs, length = degree+1, poly(1) is constant 
      !             1              2             3             4
      !        poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + poly(4) x^3
      complex*16 zeta, zeta2, zero
      real*8 third
      parameter(zero=(0d0,0d0))
      parameter(zeta =(-0.5d0, 0.8660254037844386d0))  ! sqrt3(1)
      parameter(zeta2=(-0.5d0,-0.8660254037844386d0))  ! sqrt3(1)**2
      parameter(third=0.3333333333333333d0)            ! 1/3
      complex*16 s0, s1, s2
      complex*16 E1 ! x0+x1+x2
      complex*16 E2 ! x0x1+x1x2+x2x0
      complex*16 E3 ! x0x1x2
      complex*16 A, B, a_1, E12
      complex*16 delta, A2

      complex*16 val, x
      integer  i

      a_1=poly(4)**(-1)
      E1=-poly(3)*a_1
      E2=poly(2)*a_1
      E3=-poly(1)*a_1
      
      s0=E1
      E12=E1*E1
      A=2d0*E1*E12-9d0*E1*E2+27d0*E3  ! =  s1^3 + s2^3
      B=E12-3d0*E2                    !  = s1 s2
      ! quadratic equation: z^2-Az+B^3=0  where roots are equal to s1^3 and s2^3
      A2=A*A
      delta=sqrt(A2-4d0*(B*B*B))
      if( dble(conjg(A)*delta).ge.0d0 )then  ! scallar product to decide the sign yielding bigger magnitude
        ! use dble() or realpart() instead of real() not to loose precision in ifort or g77
        s1=(0.5d0*(A+delta))**third
      else
        s1=(0.5d0*(A-delta))**third
      endif
      if(s1.eq.zero)then
        s2=zero
      else
        s2=B/s1
      endif


      x0=third*(s0+s1+s2)
      x1=third*(s0+s1*zeta2+s2*zeta )
      x2=third*(s0+s1*zeta +s2*zeta2)


      if(.false.)then  ! print the results

        x=x0
        val=poly(4)
        do i=3,1,-1
          val=val*x+poly(i)
        enddo
        write(*,'(2f19.15,a3,2f19.15)') x,' ->',val

        x=x1
        val=poly(4)
        do i=3,1,-1
          val=val*x+poly(i)
        enddo
        write(*,'(2f19.15,a3,2f19.15)') x,' ->',val

        x=x2
        val=poly(4)
        do i=3,1,-1
          val=val*x+poly(i)
        enddo
        write(*,'(2f19.15,a3,2f19.15)') x,' ->',val


      endif
      end


!-------------------------------------------------------------------!
! _11_                    DIVIDE_POLY_1                             !
!-------------------------------------------------------------------!
      subroutine divide_poly_1(polyout, remainder, p, polyin, degree)
      ! Subroutine will divide complex polynomial 'polyin' by (x-p)
      ! results will be returned in polynomial 'polyout' of degree-1
      ! The remainder of the division will be returned in 'remainder'
      !
      ! You can provide same array as 'polyin' and 'polyout' - this
      ! routine will work fine, though it will not set to zero the 
      ! unused, highest coefficient in the output array. You just have
      ! remember the proper degree of a polynomial.
      implicit none
      integer degree ! intent(in) 
      complex*16 polyout(degree)  ! intent(out)
      complex*16 remainder         ! intent(out)
      complex*16 p                ! intent(in) 
      complex*16 polyin(degree+1) ! intent(in)  ! coeffs of the polynomial
      ! poly - is an array of polynomial cooefs, length = degree+1, poly(1) is constant 
      !             1              2             3
      !        poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + ...
      integer i
      complex*16 coef, prev

      coef=polyin(degree+1)
      do i=1,degree
        polyout(i)=polyin(i)
      enddo
      do i=degree,1,-1
        prev=polyout(i)
        polyout(i)=coef
        coef=prev+p*coef
      enddo
      remainder=coef
      return
      end

!-------------------------------------------------------------------!
! _12_                    EVAL_POLY                                 !
!-------------------------------------------------------------------!
      complex*16 function eval_poly(x, poly, degree, errk)
      ! Evaluation of the complex polynomial 'poly' of a given degree 
      ! at the point 'x'. This routine calculates also the simplified
      ! Adams' (1967) stopping criterion. ('errk' should be multiplied 
      ! by 2d-15 for double precisioni,real*8, arithmetic)
      implicit none
      complex*16 x              ! intent(in)
      integer degree            ! intent(in)
      real*8 errk               ! intent(out)
      complex*16 poly(degree+1) ! intent(in)  ! coeffs of the polynomial
      ! poly - is an array of polynomial cooefs, length = degree+1, poly(1) is constant 
      !             1              2             3
      !        poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + ...
      integer i
      complex*16 val
      real*8 absx

      ! prepare stoping criterion
      errk=abs(poly(degree+1))
      val=poly(degree+1)
      absx=abs(x)
      do i=degree,1,-1  ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
        val=val*x+poly(i)
        ! Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
        ! Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
        ! ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
        ! Eq 8.
        errk=errk*absx+abs(val)
      enddo
      eval_poly=val

      ! if(abs(val)<2d-15*errk) return  ! (simplified a little Eq. 10 of Adams 1967)

      return
      end

!-------------------------------------------------------------------!
! _13_                    MULTIPLY_POLY_1                           !
!-------------------------------------------------------------------!
      subroutine multiply_poly_1(polyout, p, polyin, degree)
      ! Subroutine will multiply polynomial 'polyin' by (x-p)
      ! results will be returned in polynomial 'polyout' of degree+1
      !
      ! You can provide same array as 'polyin' and 'polyout' - this
      ! routine will work fine.
      implicit none
      integer degree               ! intent(in) ! OLD degree, new will be +1
      complex*16 polyout(degree+2) ! intent(out) 
      complex*16 p                 ! intent(in) 
      complex*16 polyin(degree+1)  ! intent(in) ! coeffs of the polynomial
      ! poly - is an array of polynomial cooefs, length = degree+1, poly(1) is constant 
      !             1              2             3
      !        poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + ...
      integer i

      do i=1,degree+1         ! copy
        polyout(i)=polyin(i)
      enddo

      polyout(degree+2)=polyout(degree+1)
      do i=degree+1,2,-1
        polyout(i)=polyout(i-1)-polyout(i)*p
      enddo
      polyout(1)=-polyout(1)*p
      return
      end


!-------------------------------------------------------------------!
! _14_                    CREATE_POLY_FROM_ROOTS                    !
!-------------------------------------------------------------------!
      subroutine create_poly_from_roots(poly, degree, a, roots)
      ! Routine will build polynomial from a set of points given in 
      ! the array 'roots'. These points will be zeros of the resulting 
      ! polynomial.
      !
      ! poly - is an array of polynomial cooefs, length = degree+1, poly(1) is constant 
      !             1              2             3
      !        poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + ...
      ! degree - is and integer denoting size of the 'roots' array
      ! a - gives the leading coefficient of the resutling polynomial
      ! roots - input array of points, size=degree
      !
      ! This subroutine works, but it is not optimal - it will work  
      ! up to a polynomial of degree~50, if you would like to have 
      ! more robust routine, up to a degree ~ 2000, you should 
      ! split your factors into a binary tree, and then multiply 
      ! leafs level by level from down up with subroutine like: 
      ! multiply_poly for arbitraty polynomial multiplications
      ! not multiply_poly_1.
      !
      implicit none
      integer degree            ! intent(in) 
      complex*16 poly(degree+1) ! intent(out)
      complex*16 a              ! intent(in) 
      complex*16 roots(degree)  ! intent(in)
      !
      integer i
      complex*16 zero
      parameter(zero=(0d0,0d0))
      ! 
      do i=1,degree+1
        poly(i)=zero
      enddo
      poly(1)=a  ! leading coeff of the polynomial
      do i=1,degree
        call multiply_poly_1(poly, roots(i), poly, i-1)
      enddo
      end


!-------------------------------------------------------------------!
! _15_                    CMPLX_LAGUERRE_WRAPPER                    !
!-------------------------------------------------------------------!
      subroutine cmplx_laguerre_wrapper
     *          (poly, degree, root, iter, success)
      implicit none
      ! subroutine parameters
      integer degree
      complex*16 poly(degree+1) ! intent(in)
      integer iter              ! intent(out)
      complex*16 root           ! intent(inout)
      logical success           ! intent(out)

      call cmplx_laguerre(poly, degree, root, iter, success)
      end


!-------------------------------------------------------------------!
! _16_                    CMPLX_NEWTON_SPEC_WRAPPER                 !
!-------------------------------------------------------------------!
      subroutine cmplx_newton__spec_wrapper 
     *           (poly, degree, root, iter, success)
      implicit none
      integer degree
      complex*16 poly(degree+1) ! intent(in)
      integer iter              ! intent(out)
      complex*16 root           ! intent(inout)
      logical success           ! intent(out)

      call cmplx_newton_spec(poly, degree, root, iter, success)
      end


!-------------------------------------------------------------------!
! _17_             CMPLX_LAGUERRE2NEWTON_WRAPPER                 !
!-------------------------------------------------------------------!
      subroutine cmplx_laguerre2newton_wrapper
     *           (poly, degree, root, iter, success, starting_mode)
      implicit none
      ! subroutine parameters
      integer degree
      complex*16 poly(degree+1) ! intent(in)
      integer iter              ! intent(out)
      complex*16 root           ! intent(inout)
      logical success           ! intent(out)
      integer starting_mode     ! intent(in)

      call cmplx_laguerre2newton
     *           (poly, degree, root, iter, success, starting_mode)


      end

!-------------------------------------------------------------------!
