      subroutine evalfun(fun, n, x, f, info)
      implicit none
      character(len=*), intent(in) :: fun
      integer(kind=4), intent(in) :: n
      real(kind=8), intent(in) :: x(n)
      real(kind=8), intent(out) :: f
      integer(kind=4), intent(out) :: info

!      integer(kind=4), external :: mod
      character(len=len(fun)) :: ful
      real(kind=8) :: y(0 : n+1)
      integer(kind=4) :: i, j, k, hk, m, j1, j2, j3, j4, j5, ml, mu
      real(kind=8) :: tmp, tmp1, tmp2, tmp3
      real(kind=8) :: ychebq(n+1, n), dp1, dp2, pd, r, s
      real(kind=8) :: alpha, beta, gamm, delta, p1, p2, p3, p4
      real(kind=8) :: f1, f2, f3, f4, f5
      real(kind=8) :: c, betai, betaj
      real(kind=8) :: zeta, kappa, p, pi, scaling
      real(kind=8) :: ysrf(ceiling(sqrt(dble(n))), 
     +ceiling(sqrt(dble(n))))
      real(kind=8) :: ytmp(n), par(n), h, t(n), MST(n, n)
      real(kind=8) :: MaS(2*n, n), MaC(2*n, n), sinx(n), cosx(n), so(n),
     + co(n), d(2*n), theta
      real(kind=8), parameter :: INFINITY = 1.0D308

      pi = acos(-1.0D0)

      y(0) = 0.0D0
      y(n+1) = 0.0D0
      do i = 1, n
          !y(i) = x(perm(i))
          y(i) = x(i)
      enddo
      f = INFINITY
      info = 0 ! info = 0 means successful evaluation

      ful = fun
      do i = 1, len(trim(fun)) 
          k = iachar(fun(i:i)) 
          if ( k >= iachar('A') .and. k <= iachar('Z') ) then 
              k = k + iachar('a') - iachar('A') 
              ful(i:i) = achar(k) 
          endif 
      enddo       

      if (trim(ful) .eq. 'arglina') then
C See J. More, 1981, Testing Unconstrained Optimization Software,
C Problem 32, also CUTEr.
           m = 2*n
           tmp = 0.0D0
           do i = 1, n
               tmp = tmp+y(i)
           enddo
           f = 0.0D0
           do i = 1, n
               f = f+(y(i)-2.0D0/dble(m)*tmp-1.0D0)**2
           enddo
           do i = n+1, m
               f = f+(-2.0D0/dble(m)*tmp-1.0D0)**2
           enddo
      elseif (trim(ful) .eq. 'arglina4') then
CCCC !!!!!! DIFFERENT FROM ARGLINA. Quartic instead of  Quadratic. !!!!
C See J. More, 1981, Testing Unconstrained Optimization Software,
C Problem 32, also CUTEr.
          m = 2*n
          tmp = 0.0D0
          do i = 1, n
              tmp = tmp+y(i)
          enddo
          f = 0.0D0
          do i = 1, n
              f = f+(y(i)-2.0D0/dble(m)*tmp-1.0D0)**4
          enddo
          do i = n+1, m
              f = f+(-2.0D0/dble(m)*tmp-1.0D0)**4
          enddo
      elseif (trim(ful) .eq. 'arglinb') then
C See J. More, 1981, Testing Unconstrained Optimization Software,
C Problem 33, also CUTEr.
          m = 2*n
          tmp = 0.0D0
          do i = 1, n
              tmp = tmp+dble(i)*y(i)
          enddo
          f = 0.0D0
          do i = 1,m
              f = f + (dble(i)*tmp-1.0D0)**2
          enddo
      elseif (trim(ful) .eq. 'arglinc') then
C See J. More, 1981, Testing Unconstrained Optimization Software,
C Problem 34, also CUTEr.
          m = 2*n
          tmp = 0.0D0
          do i = 2, n-1
              tmp = tmp+dble(i)*y(i)
          enddo
          f = 1.0D0
          do i = 2, m-1
              f = f + (dble(i-1)*tmp-1.0D0)**2
          enddo
          f = f + 1.0D0
      elseif (trim(ful) .eq. 'argtrig') then
C See J. More, 1981, Testing Unconstrained Optimization Software,
C Problem 26, also CUTEr.
          f = 0.0D0
          do i = 1, n
              tmp = dble(n)+dble(i)*(1.0D0-cos(y(i)))-sin(y(i))
              do j = 1, n
                  tmp = tmp - cos(y(j))
              enddo
              f = f + tmp*tmp
          enddo
      elseif (trim(ful) .eq. 'arwhead') then
          f = 0.0D0
          do i = 1,n-1
              f = f + (y(i)*y(i)+y(n)*y(n))**2 - 4.0D0*y(i) +3.0D0
          enddo
      elseif (trim(ful) .eq. 'bdqrtic') then
C See Powell, 2003, On trust region methods for unconstrained
C minimization without derivatives, Problem 2. Also see CUTEr.
C NOTICE: the definitions in the above references are different,
C the difference being -4Y(I)+3 V.S. (-4Y(I)+3)^2. We call the former
C one BDQRTICP (see below). The original definition
C is in Conn, Gould Lescrenoer, Toint, 1994, Performance of a
C Multifrontal scheme for partially seperable optimization, Problem
C 61. But the paper is not available for now.
          f = 0.0D0
          do i = 1, n-4
              f = f + (y(i)**2+2.0D0*y(i+1)**2+3.0D0*y(i+2)**2
     +              + 4.0D0*y(i+3)**2+5.0D0*y(n)**2)**2
     +              + (- 4.0D0*y(i) + 3.0D0)**2
          enddo
      elseif (trim(ful) .eq. 'bdqrticp') then
C See Powell, 2003, On trust region methods for unconstrained
C minimization without derivatives, Problem 2. Also see CUTEr.
C NOTICE: the definitions in the above references are different,
C the difference being -4Y(I)+3 V.S. (-4Y(I)+3)^2. Here is Powell's
C version.
          f = 0.0D0
          do i = 1,n-4
              f = f + (y(i)**2+2.0D0*y(i+1)**2+3.0D0*y(i+2)**2
     +              + 4.0D0*y(i+3)**2+5.0D0*y(n)**2)**2
     +              - 4.0D0*y(i) + 3.0D0
          enddo
      elseif (trim(ful) .eq. 'bdvalue') then
C See J. More, 1981, Testing Unconstrained Optimization Software,
C Problem 28, also CUTEr.
          f = 0.0D0
          h = 1.0D0/dble(n+1)
          do i = 1, n
          f = f+(2.0D0*y(i)-y(i-1)-y(i+1)
     +         +0.5D0*h*h*(y(i)+dble(i)/dble(n+1)+1.0D0)**3)**2
          enddo
      elseif (trim(ful) .eq. 'brownal') then
C See J. More, 1981, Testing Unconstrained Optimization Software,
C Problem 27, also CUTEr.
          f = 0.0D0
          tmp1 = 0.0D0
          tmp2 = 1.0D0
          do i = 1, n
              tmp1 = tmp1 + y(i)
              tmp2 = tmp2*y(i)
          enddo
          f = (tmp2 - 1.0D0)**2
          do i = 1, n-1
              f = f + (y(i) + tmp1 - dble(n+1))**2
          enddo
      elseif (trim(ful) .eq. 'broydn3d') then
C See CUTEr.
          f = 0.0D0
          do i = 1, n
               f = f+((3.0D0-2.0D0*y(i))*y(i)-y(i-1)-2*y(i+1)+1.0D0)**2
          enddo
      elseif (trim(ful) .eq. 'broydn7d') then
C See Ph. L. Toint "Some Numerical Results Using a Sparse Matrix
C Updating Formula in Unconstrained Optimization", Problem 3.4, also
C CUTEr.
          f = 0.0D0
          do i = 1, n
          f =f+dabs(y(i-1)-y(i)*(3.0D0-0.5D0*y(i))+2.0D0*y(i+1)-1.0D0)
     +        **(7.0D0/3.0D0)
          enddo
          do i = 1, n/2
              f = f + dabs(y(i)+y(i+n/2))**(7.0D0/3.0D0)
          enddo
      elseif (trim(ful) .eq. 'brybnd') then
C See J. More, 1981, Testing Unconstrained Optimization Software,
C Problem 31, also CUTEr.
          ml = 5
          mu = 1
          f = 0.0D0
          do i = 1, n
              tmp = 0.0D0
              do j = max(1, i-ml), min(n, i+mu)
                  if (j .ne. i) then
                      tmp = tmp+y(j)*(1.0D0+y(j))
                  endif
              enddo
              f = f + (y(i)*(2.0D0+5.0D00*y(i)*y(i))+1.0D0-tmp)**2
          enddo
      elseif (trim(ful) .eq. 'chainwoo') then
C See CUTEr.
          f = 1.0D0
          do i = 1, n/2-1
              j = 2*i
              f = f + 1.0D2*(y(j)-y(j-1)**2)**2+(1.0D0-y(j-1))**2
     +              + 9.0D1*(y(j+2)-y(j+1)**2)**2 + (1.0D0-y(j+1))**2
     +              + 1.0D1*(y(j)+y(j+2)-2.0D0)**2 +
     +                1.0D-1*(y(j)-y(j+2))**2
          enddo
      elseif (trim(ful) .eq. 'chebquad') then
C See J. More, 1981, Testing Unconstrained Optimization Software,
C Problem 35, also CUTEr.
          ychebq(1, 1:n) = 1.0D0
          ychebq(2, 1:n) = 2.0D0*y(1:n) - 1.0D0
          do i=2, n
              ychebq(i+1, 1:n) = 2.0D0*ychebq(2, 1:n)*ychebq(i, 1:n)
     +                           - ychebq(i-1, 1:n)
          enddo
          f = 0.0D0
          do i = 1, n+1
              tmp = sum(ychebq(i, 1:n))/dble(n)
              if (mod(i, 2) .eq. 1) then
                  tmp=tmp+1.0D0/dble(i*i-2*i)
              endif
              f = f + tmp*tmp
          enddo
      elseif (trim(ful) .eq. 'chpowellb') then
C Chained version of Powell Badly Scaled Function.
C See J. More, 1981, Testing Unconstrained Optimization Software,
C Problem 3.
          f = 0.0D0
          do i = 1, n-1
              f = f + (1.0D4*y(i)*y(i+1)-1.0D0)**2 +
     +            (exp(-y(i))+exp(-y(i+1)) - 1.0001D0)**2
          enddo
      elseif (trim(ful) .eq. 'chpowells') then
C Chained version of Powell Singular Function.
C See J. More, 1981, Testing Unconstrained Optimization Software,
C Problem 13, and Ying-Jie Li, Dong-Hui Li, Truncated regularized Newton
C Method for Convex Minimization, Problem 6, and Luksan, Vicek,
C Sparse and partially separable test problems for unconstrained and
C equality constrained optimization.
          f = 0.0D0
          do j = 1, (n-2)/2
              i = 2*j
              f = f + (y(i-1)+10.0D0*y(i))**2 +
     +            5.0D0*(y(i+1)-y(i+2))**2 +
     +            (y(i)-2.0D0*y(i+1))**4 +
     +            10.0D0*(y(i-1)-y(i+2))**4
          enddo
      elseif (trim(ful) .eq. 'chnrosnb') then
C See CUTEr. MODIFIED by Luksan, Matonoha, Vlcek, "Modified CUTE Probelems for
C Unconstrained Optimzation." It is another version of chained
C ROSENBROCK. Compare with CHROSEN and ROSENBROCK.
          alpha = 16.0D0*(1.5D0+sin(dble(i)))**2
          f = 0.0D0
          do i = 2, n
              f = f + alpha*(y(i-1)-y(i)**2)**2+(y(i)-1)**2
          enddo
      elseif (trim(ful) .eq. 'chrosen') then
          f = 0.0D0
          do i = 1,n-1
              f = f + (4.0D0)*(y(i)-y(i+1)*y(i+1))*
     +            (y(i)-y(i+1)*y(i+1)) +
     +            (1.0D0-y(i+1))*(1.0D0-y(i+1))
          enddo
      elseif (trim(ful) .eq. 'cosine') then
C See CUTEr.
          f = 0.0D0
          do i = 1, n-1
              f = f + cos(y(i)**2 - 0.5D0*y(i+1))
          enddo
      elseif (trim(ful) .eq. 'cragglvy') then
C See CUTEr.
          f = 0.0D0
          do j = 1, (n-2)/2
              i = 2*j
              f = f + (exp(y(i-1))-y(i))**4
     +              + 1.0D2*(y(i)-y(i+1))**6
     +              + (tan(y(i+1)-y(i+2))+y(i+1)-y(i+2))**4
     +              + y(i-1)**8 + (y(i+2)-1.0D0)**2
          enddo
      elseif (trim(ful) .eq. 'cube') then
C See CUTEr.
          f = (y(1)-1.0D0)**2
          do i = 2, n
              f = f + 100.0D0 * (y(i)-y(i-1)**3)**2
          enddo
      elseif (trim(ful) .eq. 'curly10' .or. trim(ful) .eq. 'curly20'
     +  .or.  trim(ful) .eq. 'curly30') then
C See CUTEr.
          if (trim(ful) .eq. 'curly10') then
              k = 10
          elseif (trim(ful) .eq. 'curly20') then
              k = 20
          else
              k = 30
          endif
          f = 0.0D0
          do i = 1, n
              tmp = 0.0D0
              do j = i, min(i+k,n)
                  tmp = tmp + y(j)
              enddo
           f = f + tmp*(tmp*(tmp**2-2.0D1)-1.0D-1)
          enddo
      elseif (trim(ful) .eq. 'dixmaane' .or. trim(ful) .eq. 'dixmaanf'
     +   .or. trim(ful) .eq. 'dixmaang' .or. trim(ful) .eq. 'dixmaanh'
     +   .or. trim(ful) .eq. 'dixmaani' .or. trim(ful) .eq. 'dixmaanj'
     +   .or. trim(ful) .eq. 'dixmaank' .or. trim(ful) .eq. 'dixmaanl'
     +   .or. trim(ful) .eq. 'dixmaanm' .or. trim(ful) .eq. 'dixmaann'
     +   .or. trim(ful) .eq. 'dixmaano' .or. trim(ful) .eq. 'dixmaanp')
     +   then
          m = n/3
          if (trim(ful) .eq. 'dixmaane' .or. trim(ful) .eq. 'dixmaani'
     +        .or. trim(ful) .eq. 'dixmaanm') then
              alpha = 1.0D0
              beta = 0.0D0
              gamm = 0.125D0
              delta = 0.125D0
          elseif (trim(ful) .eq. 'dixmaanf' .or. trim(ful).eq.'dixmaanj'
     +        .or. trim(ful) .eq. 'dixmaann') then
              alpha = 1.0D0
              beta = 0.0625D0
              gamm = 0.0625D0
              delta = 0.0625D0
          elseif (trim(ful) .eq. 'dixmaang' .or. trim(ful).eq.'dixmaank'
     +        .or. trim(ful) .eq. 'dixmaano') then
              alpha = 1.0D0
              beta = 0.125D0
              gamm = 0.125D0
              delta = 0.125D0
          elseif (trim(ful) .eq. 'dixmaanh' .or. trim(ful).eq.'dixmaanl'
     +        .or. trim(ful) .eq. 'dixmaanp') then
              alpha = 1.0D0
              beta = 0.26D0
              gamm = 0.26D0
              delta = 0.26D0
          endif
          if (trim(ful) .eq. 'dixmaane' .or. trim(ful) .eq. 'dixmaanf'
     +    .or. trim(ful) .eq. 'dixmaang' .or. trim(ful) .eq. 'dixmaanh')
     +    then
              p1 = 1
              p2 = 0
              p3 = 0
              p4 = 1
          elseif (trim(ful) .eq. 'dixmaani' .or. trim(ful).eq.'dixmaanj'
     +    .or. trim(ful) .eq. 'dixmaank' .or. trim(ful) .eq. 'dixmaanl')
     +    then
              p1 = 2
              p2 = 0
              p3 = 0
              p4 = 2
          elseif (trim(ful) .eq. 'dixmaanm' .or. trim(ful).eq.'dixmaann'
     +    .or. trim(ful) .eq. 'dixmaano' .or. trim(ful) .eq. 'dixmaanp')
     +    then
              p1 = 2
              p2 = 1
              p3 = 1
              p4 = 2
          endif
          f1 = 0.0D0
          do i = 1, n
              f1 = f1 + (dble(i)/dble(n))**p1*y(i)**2
          enddo
          f2 = 0.0D0
          do i = 1, n-1
              f2 = f2 + (dble(i)/dble(n))**p2*
     +                y(i)**2*(y(i+1)+y(i+1)**2)**2
          enddo
          f3 = 0.0D0
          do i = 1, 2*m
              f3 = f3 + (dble(i)/dble(n))**p3*y(i)**2*y(i+m)**4
          enddo
          f4 = 0.0D0
          do i = 1, m
              f4 = f4 + (dble(i)/dble(n))**p4*y(i)*y(i+2*m)
          enddo
          f = 1.0D0
          f = f + alpha*f1 + beta*f3 + gamm*f3 + delta*f4
      elseif (trim(ful) .eq. 'dqrtic') then
C See CUTEr.
          f = 0.0D0
          do i = 1, n
              f = f + (y(i)-dble(i))**4
          enddo
      elseif (trim(ful) .eq. 'edensch') then
C See CUTEr.
          f = 16.0D0
          do i = 1, n-1
              f = f + (y(i)-2.0D0)**4 + (y(i)*y(i+1)-2.0D0*y(i+1))**2
     +              + (y(i+1)+1.0D0)**2
          enddo
      elseif (trim(ful) .eq. 'eg2') then
C See CUTEr.
          f = 0.0D0
          do i = 1, n-1
              f = f + sin(y(1)+y(i)**2-1.0D0)
          enddo
          f = f + 0.5D0*sin(y(n)**2)
      elseif (trim(ful) .eq. 'engval1') then
C See CUTEr. Compare with ARWHEAD.
          f = 0.0D0
          do i = 1, n-1
              f = f + (y(i)**2+y(i+1)**2)**2 - 4.0D0*y(i) + 3.0D0
          enddo
      elseif (trim(ful) .eq. 'errinros') then
C See CUTEr. MODIFIED by Luksan, Matonoha, Vlcek, "Modified CUTE Probelems for
C Unconstrained Optimzation." It is another version of chained
C ROSENBROCK. Compare with CHNROSNB, CHROSEN, and ROSENBROCK.
          alpha = 16.0D0*(1.5D0+sin(dble(i)))**2
          f = 0.0D0
          do i = 2, n
              f = f + (y(i-1)-alpha*y(i)**2)**2+(y(i)-1)**2
          enddo
      elseif (trim(ful) .eq. 'expsum') then
C See http://www.opt.uni-duesseldorf.de/~jarre/dot/f_cx.m
          alpha = 4.0D0
          beta = 2.0D0*sqrt(alpha)
          f = 0.0D0
          do i = 1, n
              f = f + dble(i**2)*(exp(dble(i)*sum(y(1:i)))
     +          + alpha*exp(-dble(i)*sum(y(1:i))) - beta)
          enddo
      elseif (trim(ful) .eq. 'extrosnb') then
C See CUTEr and Luksan, Matonoha, Vlcek, "Modified CUTE Probelems for
C Unconstrained Optimzation." It is another version of chained
C ROSENBROCK. Compare with CHNROSNB, CHROSEN, ERRINROS, and ROSENBROCK.
          f = (y(1)-1.0D0)**2
          do i = 2, n
              f = f + 1.0D2*(y(i)-y(i-1)**2)**2
          enddo
      elseif (trim(ful) .eq. 'exttet') then
C See N. Andrei, 2008, An Unconstrained Optimization Test Functions
C Collection, Extended Three Term Exponentials.
          f = 0.0D0
          do j = 1, n/2
              i = 2*j
              f = f + exp(y(i-1)+3.0D0*y(i)-0.1D0)
     +              + exp(y(i-1)-3.0D0*y(i)-0.1D0)
     +              + exp(-y(i-1)-0.1D0)
              enddo
      elseif (trim(ful) .eq. 'firose') then
C Five-Diagonal ROSENBROCK. The Jacobian of the corresponding nonlinear
C equations is five-diagonal. See Example 5.3 of "The Secant/Finite
C Difference Algorithms for Solving Sparse Nonlinear Systems of Equations"
C by Guangye Li (SIAM Journal on Numerical Analysis, 1988).
          f = 0.0D0
          if (n .ge. 4) then
          f = (4.0D0*(y(1)-y(2)**2) + y(2)-y(3)**2)**2
     +        + (8.0D0*y(2)*(y(2)**2-y(1))- 2.0D0*(1.0D0-y(2))
     +        + 4.0D0*(y(2)-y(3)**2) + y(3)-y(4)**2)**2
          do i = 3, n-2
              f = f + (8.0D0*y(i)*(y(i)**2-y(i-1))-2.0D0*(1.0D0-y(i))
     +            + 4.0D0*(y(i)-y(i+1)**2)+y(i-1)**2-y(i-2)
     +            + y(i+1)-y(i+2)**2)**2
          enddo
          f = f + (8.0D0*y(n-1)*(y(n-1)**2-y(n-2))
     +          - 2.0D0*(1.0D0-y(n-1))+4.0D0*(y(n-1)-y(n)**2)
     +          + y(n-2)**2-y(n-3))**2
          f = f + (8.0D0*y(n)*(y(n)**2-y(n-1))-2.0D0*(1.0D0-y(n))
     +          + y(n-1)**2-y(n-2))**2
          endif
      elseif (trim(ful) .eq. 'fletcbv2') then
C See CUTEr.
          h = 1.0D0/dble(n+1)
          tmp1 = y(1)**2
          do i = 1, n-1
              tmp1 = tmp1 + (y(i)-y(i+1))**2
          enddo
          tmp1 = tmp1 + y(n)**2
          tmp2 = 0.0D0
          do i = 1, n
              tmp2 = tmp2 + 2.0D0*y(i)+cos(y(i))
          enddo
          f = 0.5D0*tmp1 - h**2*tmp2 -y(n)
      elseif (trim(ful) .eq. 'fletcbv3') then
C See CUTEr. MODIFIED by Luksan, Matonoha, Vlcek, "Modified CUTE Probelems for
C Unconstrained Optimzation."
          p = 1.0D-8
          h = 1.0D0/dble(n+1)
          kappa = 1.0D0
          tmp1 = y(1)**2
          do i = 1, n-1
              tmp1 = tmp1 + (y(i)-y(i+1))**2
          enddo
          tmp1 = tmp1 + y(n)**2
          tmp2 = 0.0D0
          do i = 1,n
              tmp2 = tmp2 +
     +             1.0D2*sin(y(i)*1.0D-2)*(1.0D0+2.0D0/h**2) +
     +             1.0D0/h**2*kappa*cos(y(i))
          enddo
          f = 0.5D0*p*tmp1 - p*tmp2
      elseif (trim(ful) .eq. 'fletchcr') then
C See CUTEr.
          f = 0.0D0
          do i = 1, n-1
              f = f + (y(i+1)-y(i)+1.0D0-y(i)**2)**2
          enddo
          f = 1.0D2*f
      elseif (trim(ful) .eq. 'fminsrf2') then
C See CUTEr
          f = 0.0D0
          ysrf = 0.0D0
          k = floor(sqrt(dble(n)))
          hk = floor(dble(k)/2.0D0)
          if (k .le. 1) then
              info = 2
          endif
          do i = 1, k
              do j = 1, k
                ysrf(i, j) = x((i-1)*k+j)
              enddo
          enddo
          do i = 1, k-1
              do j = 1, k-1
                  f = f+sqrt(1.0D0+0.5D0*dble(k-1)**2*
     +                (((ysrf(i,j)-ysrf(i+1,j+1))**2
     +                +(ysrf(i+1,j)-ysrf(i,j+1))**2)))
              enddo
          enddo
          f = f/dble(k-1)**2 + ysrf(hk, hk)**2/dble(n)
      elseif (trim(ful) .eq. 'freuroth') then
          tmp1 = 0.0D0
          tmp2 = 0.0D0
          do i = 1, n-1
              tmp1 = tmp1 + ((5.0D0-y(i+1))*y(i+1)**2 + y(i) -
     +           2.0D0*y(i+1)-13.0D0)**2
              tmp2 = tmp2 +
     +           ((1.0D0+y(i+1))*y(i+1)**2+y(i)-14.0D0*y(i+1)-29.0D0)**2
          enddo
          f = tmp1 + tmp2
      elseif (trim(ful) .eq. 'genbrown') then
C Generalized Brown Function 1. See Ying-Jie Li, Dong-Hui Li, Truncated
C regularized Newton Method for Convex Minimization, Problem 6,
C and Luksan, Vicek, Sparse and partially separable test problems for
C unconstrained and equality constrained optimization.
          f = 0.0D0
          do i = 2, n
              f = f + (y(i-1)-3.0D0)**2 + (y(i-1)-y(i))**2
     +            + exp(20.0D0*(y(i-1)-y(i)))
          enddo
      elseif (trim(ful) .eq. 'genhumps') then
C See CUTEr.
          zeta = 2.0D0
          f = 0.0D0
          do i = 1, n-1
              f = f + sin(zeta*y(i))**2*sin(zeta*y(i+1))**2
     +              + 0.05D0*(y(i)**2+y(i+1)**2)
          enddo
      elseif (trim(ful) .eq. 'genrose') then
C See CUTEr.
C Compare with CHROSEN and ROSENBROCK.
C The starting point is (I/(N+1)).
          f = 1.0D0
          do i = 2, n
              f = f + 1.0D2*(y(i)-y(i-1)**2)**2+(y(i)-1.0D0)**2
          enddo
      elseif (trim(ful) .eq. 'indef') then
C See CUTEr. MODIFIED by Luksan, Matonoha, Vlcek, "Modified CUTE Probelems for
C Unconstrained Optimzation."
          tmp1 = 0.0D0
          do i = 1, n
              tmp1 = tmp1 + 1.0D2*sin(1.0D-2*y(i))
          enddo
          tmp2 = 0.0D0
          do i = 2, n-1
               tmp2 = tmp2 + cos(2.0D0*y(i)-y(n)-y(1))
          enddo
          f = tmp1 + 0.5D0*tmp2
      elseif (trim(ful) .eq. 'integreq') then
C See J. More, 1981, Testing Unconstrained Optimization Software,
C Problem 29 and also CUTEr.
          h = 1.0D0/dble(n+1)
          do i =1, n
              t(i) = dble(i)/dble(n+1)
          enddo
          f = 0.0D0
          do i = 1, n
              tmp1 = 0.0D0
              do j = 1, i
                  tmp1 = tmp1 + t(j)*(y(j) + t(j) + 1.0D0)**3
              enddo
              tmp2 = 0.0D0
              do j = i+1, n
              tmp2 = tmp2 + (1.0D0-t(j))*(y(j)+t(j)+1.0D0)**3
              enddo
              f = f +
     +             (y(i) + 0.5D0*h*((1-t(i))*tmp1+t(i)*tmp2))**2
          enddo
      elseif (trim(ful) .eq. 'liarwhd') then
          f = 0.0D0
          do i =1 ,n
              f = f + 4.0D0*(y(i)**2-y(1))**2+(y(i)-1.0D0)**2
          enddo
      elseif (trim(ful) .eq. 'lilifun3') then
C See Ying-Jie Li, Dong-Hui Li, Truncated
C regularized Newton Method for Convex Minimization, Problem 3.
          f = 0.0D0
          do i = 2, n
              f = f + exp(y(i)-y(i-1))**2 + (y(i)-y(i-1))**2
     +              + 2.0D0*y(i)**4+4.0D0*y(i-1)**4
          enddo
      elseif (trim(ful) .eq. 'lilifun4') then
C See Ying-Jie Li, Dong-Hui Li, Truncated
C regularized Newton Method for Convex Minimization, Problem 4.
          f = 0.0D0
          do i = 2, n
              f = f + 0.5D0*(y(i)-y(i-1))**2 + sin(y(i)-y(i-1))
     +              + 2.0D0*(2.0D0*y(i)+3.0D0*y(i-1)-15.0D0)**4
          enddo
      elseif (trim(ful) .eq. 'morebv' .or. trim(ful).eq. 'morebvl') then
C See CUTEr. MOREBVL uses the start point suggested by Luksan, Matonoha,
C Vlcek, "Modified CUTE Probelems for
C Unconstrained Optimzation."
          h = 1.0D0/dble(n+1)
          f = 0.0D0
          do i = 1, n
              f = f + (2.0D0*y(i)-y(i-1)-y(i+1)
     +              + 0.5D0*h*h*(y(i)+dble(i)*h+1.0D0)**3)**2
          enddo
      elseif (trim(ful) .eq. 'ncb20') then
C See CUTEr. Corrected by Luksan, Matonoha, Vlcek, "Modified CUTE Probelems for
C Unconstrained Optimzation."
          f = 0.0D0
          if (n .ge. 20) then
              tmp1 = 0.0D0
              do i = 1, n-30
                  tmp = 0.0D0
                  do j = 1, 20
                      tmp = tmp + y(i+j-1)/(1+y(i+j-1)**2)
                  enddo
                  tmp1 = tmp1 + 10.0D0/dble(i)*tmp**2
                  tmp = 0.0D0
                  do j = 1, 20
                      tmp = tmp + y(i+j-1)
                  enddo
                  tmp1 = tmp1 - 0.2D0*tmp
              enddo
              tmp2 = 0.0D0
              do i = 1, n-10
                  tmp2 = tmp2 + y(i)**4+2.0D0
              enddo
              tmp3 = 0.0D0
              do i = 1, 10
                  tmp3 = tmp3 +
     +                   y(i)*y(i+10)*y(i+n-10) +2.0D0*y(i+n-10)**2
              enddo
              f = 2.0D0 + tmp1 + tmp2 + 1.0D-4*tmp3
          endif
      elseif (trim(ful) .eq. 'ncb20b') then
C See CUTEr. Corrected by Luksan, Matonoha, Vlcek, "Modified CUTE Probelems for
C Unconstrained Optimzation."
          f = 0.0D0
          if (n .ge. 20) then
              tmp1 = 0.0D0
              do i = 1, n-19
                  tmp = 0.0D0
                  do j = 1, 20
                      tmp = tmp + y(i+j-1)/(1+y(i+j-1)**2)
                  enddo
                  tmp1 = tmp1 + 10.0D0/dble(i)*tmp**2
                  tmp = 0.0D0
                  do j = 1, 20
                      tmp = tmp + y(i+j-1)
                  enddo
                  tmp1 = tmp1 - 0.2D0*tmp
              enddo
              tmp2 = 0.0D0
              do i = 1, n
                  tmp2 = tmp2 + (1.0D2*y(i)**4+2.0D0)
              enddo
              f = tmp1 +tmp2
          endif
      elseif(trim(ful) .eq. 'noncvxu2') then
C See CUTEr.
          f = 0.0D0
          do i = 1, n
              j = mod(3*i-2,n)
              k = mod(7*i-3,n)
              tmp = y(i) + y(j+1) + y(k+1)
              f = f + tmp**2 + 4.0D0*cos(tmp)
          enddo
      elseif(trim(ful) .eq. 'noncvxun') then
C See CUTEr.
          f = 0.0D0
          do i = 1, n
              j = mod(2*i-1,n)
              k = mod(3*i-1,n)
              tmp = y(i) + y(j+1) + y(k+1)
              f = f + tmp**2 + 4.0D0*cos(tmp)
          enddo
      elseif (trim(ful) .eq. 'nondia') then
C See CUTEr.
          f = (y(1)-1.0D0)**2
          do i = 2, n
              f = f + (1.0D2*y(1)-y(i-1)**2)**2
          enddo
      elseif (trim(ful) .eq. 'nondquar') then
          f = (y(1)-y(2))**2 + (y(n-1)-y(n))**2
          do i = 1, n-2
              f = f + (y(i) + y(i+1) + y(n))**4
          enddo
      elseif (trim(ful) .eq. 'penalty1') then
          dp1 = 0.0D0
          dp2 = 0.0D0
          do i = 1, n
              dp1 = dp1 + (y(i)-1.0D0)**2
              dp2 = dp2 + y(i)**2
          enddo
          f=1.0D-5*dp1 +(0.25D0-dp2)**2
      elseif (trim(ful) .eq. 'penalty2') then
          f = 0.0D0
          tmp = 0.0D0
          do i = 2,n
              f = f +
     +            (exp(y(i-1)*0.1D0)+exp(y(i)*0.1D0)-
     +            exp(dble(i-1)*0.1D0)
     +            -exp(dble(i)*0.1D0))**2 + (exp(y(i)*0.1D0)
     +            -exp(-0.1D0))**2
          enddo
          do i =1,n
              tmp = tmp + dble(n-i+1)*y(i)*y(i)
          enddo
          f = f + (1.0D0-tmp)*(1.0D0-tmp) + (y(1)-0.2D0)*(y(1)-0.2D0)
      elseif (trim(ful) .eq. 'penalty3' .or. trim(ful) .eq. 'penalty3p')
     +       then
          f = 0.0D0
          r = 0.0D0
          s = 0.0D0
          pd = 0.0D0
          do i = 1, n
              pd = pd + (y(i)**2-dble(n))
          enddo
          if (trim(ful) .eq. 'penalty3') then
              alpha = 1.0D0
          else
              alpha = 1.0D-3  ! Suggested by Powell, "The NEWUOA Software
C            for unconstrained optimization without derivatives", 2004,
C            section 8.
          endif
          do i = 1, n-2
              r = r + (y(i) + 2.0D0* y(i+1) + 1.0D1*y(i+2) - 1.0D0)**2
              s = s + (2.0D0*y(i) + y(i+1) - 3.0D0)**2
          enddo
          do i = 1, n/2
              f = f + (y(i) - 1.0D0)**2
          enddo
          f = f + pd**2 + alpha*(1.0D0+r*exp(y(n)) + s*exp(y(n-1)) +r*s)
      elseif (trim(ful) .eq. 'powellsg') then
          f = 0.0D0
          do i = 1, n/4
             j = 4*(i-1) + 1
             f = f + (y(j)+1.0D1*y(j+1))**2 +5.0D0*(y(j+2)-y(j+3))**2
     +             + (y(j+1)-2.0D0*y(j+2))**4+1.0D1*(y(j)-y(j+3))**4
          enddo
      elseif (trim(ful) .eq. 'power') then
c see cuter.
          f = 0.0D0
          do i = 1, n
              f = f + (dble(i)*y(i))**2
          enddo
      elseif (trim(ful) .eq. 'rosenbrock') then
C In claasical Rosenbrock function, alpha = 100.0D0
C When alpha = 4.0D0, this function is essentially the same as Chrosen,
C except the order of the variables.
          alpha = 100.0D0
          f = 0.0D0
          do i = 1, n-1
              f = f + (1.0D0-y(i))*(1.0D0-y(i)) +
     +            alpha*(y(i+1)-y(i)*y(i))*(y(i+1) - y(i)*y(i))
          enddo
      elseif (trim(ful) .eq. 'sbrybnd' .or. trim(ful) .eq. 'sbrybndl')
     +       then
C See CUTEr. And see Luksan, Matonoha, Vlcek, "Modified CUTE Probelems for
C Unconstrained Optimzation." The scaling parameter of Luksan is 6
C instead of 12. There is a typo in Luksan: the square is missed, which
C makes the function unbounded from below.
          ml = 5
          mu = 1
          if (trim(ful) .eq. 'sbrybnd') then
              scaling = 12.0D0
          else
              scaling = 6.0D0
          endif
          f = 0.0D0
          if (n .ge. 2) then
              do i = 1, n
                  ytmp(i) =exp(scaling*(dble(i-1)/dble(n-1)))*y(i)
              enddo
              do i = 1, n
                  tmp = 0.0D0
                  do j = max(1, i-ml), min(n,i+mu)
                      if (j .ne. i) then
                          tmp = tmp+ytmp(j)*(1.0D0+ytmp(j))
                      endif
                  enddo
                  f = f + (ytmp(i)*(2.0D0+5.0D00*ytmp(i)*ytmp(i))
     +                 + 1.0D0-tmp)**2
              enddo
          endif
      elseif (trim(ful) .eq. 'schmvett') then
C See CUTEr.
          f = 0.0D0
          do i = 1, n -2
              f = f - 1.0D0/(1.0D0+(y(i)-y(i+1))**2)
     +              - sin(0.5D0*(pi*y(i+1)+y(i+2)))
     +              - exp(-((y(i)+y(i+2))/y(i+1)-2.0D0)**2)
          enddo
      elseif (trim(ful) .eq. 'scosine' .or.trim(ful).eq.'scosinel') then
C See CUTEr. And see Luksan, Matonoha, Vlcek, "Modified CUTE Probelems for
C Unconstrained Optimzation." The scaling parameter of Luksan is 6
C instead of 12.
          if (trim(ful) .eq. 'scosine') then
              scaling = 12.0D0
          else
              scaling = 6.0D0
          endif
          f = 0.0D0
          if (n .ge. 2) then
              do i = 1, n
                  ytmp(i) =exp(scaling*(dble(i-1)/dble(n-1)))*y(i)
              enddo
              do i = 1, n-1
                  f = f +cos(ytmp(i)**2-0.5D0*ytmp(i+1))
              enddo
          endif
      elseif (trim(ful) .eq. 'serose') then
C Seven-Diagonal ROSENBROCK. The Jacobian of the corresponding nonlinear
C equations is 7-diagonal. See Example 5.4 of "The Secant/Finite
C Difference Algorithms for Solving Sparse Nonlinear Systems of Equations"
C by Guangye Li (SIAM Journal on Numerical Analysis, 1988).
C The formulation in the paper seems not correct. The cases with j=3 and
C j=n-2 should be defined seperately as well.
          f = 0.0D0
          if (n .ge. 6) then
          f = (4.0D0*(y(1)-y(2)**2) + y(2)-y(3)**2 + y(3)-y(4)**2)**2
          f = f + (8.0D0*y(2)*(y(2)**2-y(1))-2.D0*(1.0D0-y(2))+
     +          4.0D0*(y(2)-y(3)**2)+y(3)-y(4)**2+y(4)-y(5)**2)**2
          f = f + (8.0D0*y(3)*(y(3)**2-y(2))-2.D0*(1.0D0-y(3))+
     +   4.0D0*(y(3)-y(4)**2)+y(2)**2-y(1)+y(4)-y(5)**2+y(5)-y(6)**2)**2
          do i = 4, n-3
              f = f + (8.0D0*y(i)*(y(i)**2-y(i-1))-2.0D0*(1.0D0-y(i))
     +              + 4.0D0*(y(i)-y(i+1)**2)
     +              + y(i-1)**2-y(i-2)+y(i+1)-y(i+2)**2
     +              + y(i-2)**2-y(i-3)+y(i+2)-y(i+3)**2)**2
          enddo
          f = f + (8.0D0*y(n-2)*(y(n-2)**2-y(n-3))-2.0D0*(1.0D0-y(n-2))+
     +     4.0D0*(y(n-2)-y(n-1)**2)+y(n-3)**2-y(n-4)+y(n-1)-y(n)**2+
     +     y(n-4)**2-y(n-5))**2
          f = f + (8.0D0*y(n-1)*(y(n-1)-y(n-2))-2.0D0*(1.0D0-y(n-1))+
     +     4.0D0*(y(n-1)-y(n)**2)+y(n-2)**2-y(n-3)+y(n-3)**2-y(n-4))**2
          f = f + (8.0D0*y(n)*(y(n)**2-y(n-1))-2.0D0*(1.0D0-y(n))
     +          + y(n-1)**2-y(n-2)+y(n-2)**2-y(n-3))**2
          endif
      elseif(trim(ful) .eq. 'sinquad') then
C See CUTEr.
          f = (y(1)-1.0D0)**4
          do i = 2, n-1
              f = f + (sin(y(i)-y(n))-y(1)**2+y(i)**2)**2
          enddo
          f = f +(y(n)**2-y(1)**2)**2
      elseif (trim(ful) .eq. 'sphrpts') then
          f = 0.0D0
          if (mod(n,2) .ne. 0) then
              info = 2
          else
              do i = 2, n/2
                  do j = 1, i-1
                      f = f + 1.0D0/((cos(y(2*i-1))*cos(y(2*i))
     +                  - cos(y(2*j-1))*cos(y(2*j)))**2
     +                  + (sin(y(2*i-1))*cos(y(2*i))
     +                  - sin(y(2*j-1))*cos(y(2*j)))**2
     +                  + (sin(y(2*i))-sin(y(2*j)))**2)
                  enddo
              enddo
          endif
      elseif (trim(ful) .eq. 'sparsine') then
          f = 0.0D0
          do i = 1, n
              j1 = mod(2*i-1,n)+1
              j2 = mod(3*i-1,n)+1
              j3 = mod(5*i-1,n)+1
              j4 = mod(7*i-1,n)+1
              j5 = mod(11*i-1,n)+1
              f = f + dble(i)*(sin(y(i))+sin(y(j1))+sin(y(j2))
     +              + sin(y(j3))+sin(y(j4))+sin(y(j5)))**2
          enddo
          f = 0.5D0*f
      elseif (trim(ful) .eq. 'sparsqur') then
          f = 0.0D0
          do i = 1, n
              j1 = mod(2*i-1,n)+1
              j2 = mod(3*i-1,n)+1
              j3 = mod(5*i-1,n)+1
              j4 = mod(7*i-1,n)+1
              j5 = mod(11*i-1,n)+1
              f = f + dble(i)*(y(i)**2+y(j1)**2+y(j2)**2
     +              + y(j3)**2+y(j4)**2+y(j5)**2)**2
          enddo
          f = 0.125D0*f
      elseif (trim(ful) .eq. 'spmsrtls') then
          do i = 1, n
              par(i) = sin(dble(i)**2)
          enddo
          m = (n+2)/3
          f = 0.0D0
          do i = 1, m
              f1 = 0.0D0
              f2 = 0.0D0
              f3 = 0.0D0
              f4 = 0.0D0
              f5 = 0.0D0
              j = 3*(i-1) + 1
              if (i .gt. 2) then
                  f1 = (y(j-4)*y(j-1)-par(j-4)*par(j-1))**2
              endif
              f3 = (y(j)**2-par(j)**2)**2
              if (i .gt. 1) then
                  f2 = (y(j-3)*y(j-1)+y(j-1)*y(j)
     +                 - par(j-3)*par(j-1)-par(j-1)*par(j))**2
                  f3 = f3 + (y(j-2)*y(j-1)-par(j-2)*par(j-1))**2
              endif
              if (i .lt. m) then
                  f3 = f3 + (y(j+2)*y(j+1)-par(j+2)*par(j+1))**2
                  f4 = (y(j+3)*y(j+1)+y(j+1)*y(j)
     +                 - par(j+3)*par(j+1)-par(j+1)*par(j))**2
              endif
              if (i .lt. m-1) then
                  f5 = (y(j+4)*y(j+1)-par(j+4)*par(j+1))**2
              endif
              f = f + f1 + f2 + f3 + f4 + f5
          enddo
      elseif (trim(ful) .eq. 'srosenbr') then
          f = 0.0D0
          do i = 1, n/2
              f = f + 1.0D2*(y(2*i)-y(2*i-1)**2)**2
     +              + (y(2*i-1)-1.0D0)**2
          enddo
      elseif (trim(ful) .eq. 'stmod') then
C Modified Stybilsky-Tang-function. 
C See http://www.opt.uni-duesseldorf.de/~jarre/dot/f_stmod.m
          f = 0.0D0
          MST = 0.0D0
          do i = 1, n
              MST(1:i, i) = dble(i);
              MST(i, i) = MST(i, i) + 1.0D0
          enddo
          y(1:n) = matmul(MST, y(1:n))
          do i = 1, n
              f = f + y(i)**4 - 1.6D1*y(i)**2 + 3.5D1*y(i) 
          enddo
          f = f + dble(n)*1.7119854821721323D2
      elseif (trim(ful) .eq. 'tointgss') then
          f = 0.0D0
          do i = 1, n-2
              f = f + (1.0D1/dble(n+2)+y(i+2)**2)*(2.0D0-
     +              exp(-(y(i)-y(i+1))**2/(1.0D-1+y(i+2)**2)))
          enddo
      elseif (trim(ful) .eq. 'tointtrig') then
C See Ph. L. Toint, "Some Numerical Results Using a Sparse Matrix
C Updating Formula in Unconstrained Optimization", Example 3.6.
          f = 0.0D0
          do i = 2, n
              betai = 1.0D0 + 1.0D-1*dble(i)
              do j = 1, i-1
                  alpha = dble(5*(1+mod(i,5)+mod(j,5)))
                  betaj = 1.0D0 + 1.0D-1*dble(j)
                  c = dble(i+j)*1.0D-1
                  f = f + alpha*sin(betai*y(i)+betaj*y(j)+c)
              enddo
          enddo
      elseif (trim(ful) .eq. 'tquartic') then
          f = (y(1)-1.0D0)**2
          do i = 1, n-1
              f = f + (y(1)**2-y(i+1)**2)**2
          enddo
      elseif (trim(ful) .eq. 'trigsabs') then
C See Powell, "The NEWUOA Software for unconstrained optimization without derivatives", 2004, section 8.
C The random numbers are replaced by "irregular" numbers generated with trigonometrics.
          do i = 1, 2*n
              do j = 1, n
                  MaS(i,j) = 1.0D2*sin(dble(i*j*n*163)**(2.5D-1))
                  MaC(i,j) = 1.0D2*cos(dble(i+j+n+42)/4.0D0)
              enddo
          enddo
          do j = 1, n
              sinx(j) = sin(y(j))
              cosx(j) = cos(y(j))
              so(j) = sin(pi*cos(dble(n+j+163)/3.0D0))
              co(j) = cos(pi*cos(dble(n+j+163)/3.0D0))
          enddo
          d = matmul(MaS, sinx-so) + matmul(MaC, cosx-co)
          f = sum(abs(d))
      elseif (trim(ful) .eq. 'trigssqs') then
C See Powell, "The NEWUOA Software for unconstrained optimization without derivatives", 2004, section 8.
C The random numbers are replaced by "irregular" numbers generated with trigonometrics.
          do i = 1, 2*n
              do j = 1, n
                  MaS(i,j) = 1.0D2*sin(dble(i*j*n*163)**(2.5D-1))
                  MaC(i,j) = 1.0D2*cos(dble(i+j+n+42)/4.0D0)
              enddo
          enddo
          do j = 1, n
              theta = (1.0D0+sin(dble(n*j*42)**(1.0D0/3.0D0)))/2.0D0
              theta = 1.0D1**(-theta)
              sinx(j) = sin(theta*y(j))
              cosx(j) = cos(theta*y(j))
              so(j) = sin(pi*cos(dble(n+j+163)/3.0D0))
              co(j) = cos(pi*cos(dble(n+j+163)/3.0D0))
          enddo
          d = matmul(MaS, sinx-so) + matmul(MaC, cosx-co)
          f = dot_product(d, d) 
      elseif (trim(ful) .eq. 'trirose1') then
C Tri-Diagonal ROSENBROCK. The Jacobian of the corresponding nonlinear
C equations is 3-diagonal. See Example 5.1 of "The Secant/Finite
C Difference Algorithms for Solving Sparse Nonlinear Systems of Equations"
C by Guangye Li (SIAM Journal on Numerical Analysis, 1988).
          f = 0.0D0
          if (n .ge. 2) then
          f = 64.0D0*(y(1)-y(2)**2)**2
          do i = 2, n-1
              f = f + (16.0D0*y(i)*(y(i)**2-y(i-1)) - 2.0D0*(1.0D0-y(i))
     +              + 8.0D0*(y(i)-y(i+1)**2))**2
          enddo
          f = f + (16.0D0*y(n)*(y(n)**2-y(n-1)) - 2.0D0*(1.0D0-y(n)))**2
          endif
      elseif (trim(ful) .eq. 'trirose2') then
C Tri-Diagonal ROSENBROCK. The Jacobian of the corresponding nonlinear
C equations is 3-diagonal. See Example 5.2 of "The Secant/Finite
C Difference Algorithms for Solving Sparse Nonlinear Systems of Equations"
C by Guangye Li (SIAM Journal on Numerical Analysis, 1988).
          f = 0.0D0
          if (n .ge. 2) then
          f = 16.0D0*(y(1)-y(2)**2)**2
          do i = 2, n-1
              f = f + (8.0D0*y(i)*(y(i)**2-y(i-1)) - 2.0D0*(1.0D0-y(i))
     +              + 4.0D0*(y(i)-y(i+1)**2))**2
          enddo
          f = f + (8.0D0*y(n)*(y(n)**2-y(n-1)) - 2.0D0*(1.0D0-y(n)))**2
          endif
      elseif (trim(ful) .eq. 'vardim') then
          f = 0.0D0
          do i = 1, n
              f = f + dble(i)*(y(i)-1.0D0)
          enddo
          f = f*f
          f = f + f*f
          do i = 1, n
              f = f + (y(i) - 1.0D0)*(y(i) - 1.0D0)
          enddo
      elseif (trim(ful) .eq. 'woods') then
          f = 0.0D0
          do i = 1, n/4
              j = 4*i
              f = f + 1.0D2*(y(j-2)-y(j-3)**2)**2+(1.0D0-y(j-3))**2
     +              + 9.0D1*(y(j)-y(j-1)**2)**2 + (1.0D0-y(j-1))**2 +
     +              1.0D1*(y(j-2)+y(j)-2.0D0)**2+1.0D-1*(y(j-2)-y(j))**2
          enddo
      else
          info = 1; ! info = 1 means unknown function name
      endif

      return
      endsubroutine evalfun
