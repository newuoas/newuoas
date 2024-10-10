      subroutine setuptest (fun, n, x, rhobeg, fopt, info)
! Set up x0, rhobeg, and fopt for the test.
! This version is not inteded to be released. It is only for test. 
!
! All rights reserved. 
!
! ZHANG Zaikun, 08/08/2016
! Department of Applied Mathematics, The Hong Kong Polytechnic University

      implicit none
      character(len=*), intent(in) :: fun 
      integer(kind=4), intent(in) :: n
      real(kind=8), intent(out) :: x(n), rhobeg, fopt
      integer(kind=4), intent(out) :: info

!      integer(kind=4), external :: mod 
      character(len=len(fun)) :: ful
      integer(kind=4) :: i, j, k 
      real(kind=8) :: ysrf(ceiling(sqrt(dble(n))),
     +ceiling(sqrt(dble(n)))), scaling, xo, pt, pi, theta

      pi = acos(-1.0D0)

      x = 1.0D308
      rhobeg = 1.0D0
      fopt = 1.0D308
      info = 0

      ful = fun
      do i = 1, len(trim(fun)) 
          k = iachar(fun(i:i)) 
          if ( k >= iachar('A') .and. k <= iachar('Z') ) then 
              k = k + iachar('a') - iachar('A') 
              ful(i:i) = achar(k) 
          endif 
      enddo       

      if (trim(ful) .eq. 'arglina') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'arglina4') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'arglinb') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'arglinc') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'argtrig') then
          do i = 1, n
              x(i) = 1.0D0/dble(n)
          enddo
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'arwhead') then
          x = 1.0D0
          rhobeg = 0.5D0
      elseif (trim(ful) .eq. 'bdqrtic') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'bdqrticp') then
          x = 1.0D0
          rhobeg = 0.5D0
      elseif (trim(ful) .eq. 'bdvalue') then
          do i = 1, n
              x(i) = dble(i*(i-n-1))/dble((n+1)*(n+1)) 
          enddo
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'brownal') then
          x = 5.0D-1
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'broydn3d') then
          x = -1.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'broydn7d') then
          x = -1.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'brybnd') then
          x = -1.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'chainwoo') then
          do i = 1, n
              if (i .eq.1 .or. i.eq.3) then 
                  x(i) = -3.0D0
              elseif (i .eq.2 .or. i.eq.4) then 
                  x(i) = -1.0D0
              else
                  x(i) = -2.0D0 
              endif
          enddo
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'chebquad') then
          do i = 1, n
              x(i)=dble(i)/dble(n+1)
          enddo
          rhobeg=0.2D0*x(1)
      elseif (trim(ful) .eq. 'chpowellb') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'chpowells') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'chnrosnb') then
          x = -1.0D0
          !rhobeg = 0.5D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'chrosen') then
          x = -1.0D0
          rhobeg = 0.5D0
      elseif (trim(ful) .eq. 'cosine') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'cragglvy') then
          x = 2.0D0
          x(1) = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'curly10' .or. trim(ful) .eq. 'curly20'
     +       .or. trim(ful) .eq. 'curly30') then
          do i = 1, n
              x(i) = (1.0D-4)/(dble(n+1)) 
          enddo
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'cube') then
          x(1) = -1.2D0
          if (n .ge. 2) x(2) = 1.0D0
          do i = 3, n
              x(i) = x(i-2)
          enddo
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'dixmaane' .or. trim(ful) .eq. 'dixmaanf' 
     +   .or. trim(ful) .eq. 'dixmaang' .or. trim(ful) .eq. 'dixmaanh'
     +   .or. trim(ful) .eq. 'dixmaani' .or. trim(ful) .eq. 'dixmaanj'
     +   .or. trim(ful) .eq. 'dixmaank' .or. trim(ful) .eq. 'dixmaanl'
     +   .or. trim(ful) .eq. 'dixmaanm' .or. trim(ful) .eq. 'dixmaann'
     +   .or. trim(ful) .eq. 'dixmaano' .or. trim(ful) .eq. 'dixmaanp') 
     +   then
          do i = 1, n
              x(i) = 2.0D0 
          enddo
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'dqrtic') then
          x = 2.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'edensch') then
          x = 0.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'eg2') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'engval1') then
          x = 2.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'errinros') then
          x = -1.0D0
          !rhobeg = 0.5D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'expsum') then
          x = 0.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'extrosnb') then
          x = -1.0D0
          !rhobeg = 0.5D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'exttet') then
          x = 1.0D-1 
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'firose') then
          x = -2.0D0
          rhobeg = 0.5D0
      elseif (trim(ful) .eq. 'fletcbv2' .or. trim(ful) .eq. 'fletcbv3') 
     +    then
          do i = 1, n
              x(i) = dble(i)/dble(n+1)
          enddo
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'fletchcr') then
          do i = 1, n
              x(i) = 0.0D0
          enddo
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'fminsrf2') then
          x = 0.0D0 
          ysrf = 0.0D0
          k = floor(sqrt(dble(n)))
          if (k .le. 1) then
              info = 2
          endif
          do i = 1, k
              ysrf(1, i) = dble(4*(i-1))/dble(k-1) + 1.0D0
              ysrf(i, 1) = dble(8*(i-1))/dble(k-1) + 9.0D0
              ysrf(k, i) = dble(4*(i-1))/dble(k-1) + 5.0D0
              ysrf(i, k) = dble(8*(i-1))/dble(k-1) + 1.0D0
          enddo
          do i = 1, k
              do j = 1, k
                  x((i-1)*k+j) = ysrf(i, j)
              enddo
          enddo
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'freuroth') then
          x = 0.0D0
          x(1) = 0.5D0
          if (n .ge. 2) then
              x(2) = -2.0D0
          endif
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'genbrown') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'genhumps') then
          x = -5.062D2
          x(1) = -5.060D2
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'genrose') then
          do i = 1, n
              x(i) = dble(i)/dble(n+1) 
          enddo
C            rhobeg = 0.5D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'indef') then
          do i = 1, n
              x(i) = dble(i)/dble(n+1) 
          enddo
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'integreq') then
          do i = 1, n
              x(i) = dble(i*(i-n-1))/dble((n+1)*(n+1)) 
          enddo
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'liarwhd') then
          x = 4.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'lilifun3') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'lilifun4') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'morebv') then
          do i = 1, n
              x(i) = dble(i*(i-n-1))/dble((n+1)*(n+1)) 
          enddo
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'morebvl') then
          x = 0.0D0 ! Suggested by  Luksan, Matonoha, Vlcek, "Modified
                    ! CUTE probelems for unconstrained optimization.
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'ncb20') then
          if (n .lt. 20) then
              x = 0.0D0
          else 
              x(1:n-10) = 0.0D0
              x(n-9:n) = 1.0D0
          endif
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'ncb20b') then
          x = 0.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'noncvxu2') then
          do i = 1,n
              x(i) = dble(i)
          enddo
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'noncvxun') then
          do i = 1,n
              x(i) = dble(i)
          enddo
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'nondia') then
          x = -1.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'nondquar') then
          do i = 1, n
              if (mod(i,2) .eq. 1) x(i) = 1.0D0    
              if (mod(i,2) .eq. 0) x(i) = -1.0D0    
          enddo
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'penalty1') then
          do i = 1, n
              x(i) = dble(i)
          enddo
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'penalty2') then
          x = 0.5D0
          rhobeg = 0.1D0
      elseif (trim(ful) .eq. 'penalty3') then
          do i = 1, n
              x(i) = dble(i)/dble(n+1)
          enddo
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'penalty3p') then
          x = 0.0D0  ! Suggested by Powell, "The NEWUOA software 
C            for unconstrained optimization without
C            derivatives", 2004, Section 8.                       
          rhobeg = 0.1D0
      elseif(trim(ful) .eq. 'powellsg') then
          do i = 1, n
              if (mod(i,4) .eq. 1) x(i) = 3.0D0
              if (mod(i,4) .eq. 2) x(i) = -1.0D0
              if (mod(i,4) .eq. 3) x(i) = 0.0D0
              if (mod(i,4) .eq. 0) x(i) = 1.0D0
          enddo
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'power') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'rosenbrock') then
          x = -1.0D0
          rhobeg = 0.5D0
      elseif (trim(ful) .eq. 'sbrybnd') then
          scaling = 1.2D1
          do i = 1, n
              x(i) = exp(scaling*dble(1-i)/dble(n-1))
          enddo
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'sbrybndl') then
          scaling = 1.2D1
          do i = 1, n
              x(i) = exp(scaling*dble(1-i)/dble(n-1))
          enddo
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'schmvett') then
          x = 3.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'scosine') then
          scaling = 1.2D1
          do i = 1, n
              x(i) = exp(scaling*dble(1-i)/dble(n-1))
          enddo
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'scosinel') then
          scaling = 1.2D1
          do i = 1, n
              x(i) = exp(scaling*dble(1-i)/dble(n-1))
          enddo
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'serose') then
          x = -2.0D0
          rhobeg = 0.5D0
      elseif (trim(ful) .eq. 'sinquad') then
          x = 0.1D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'sparsine') then
          x = 0.5D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'sparsqur') then
          x = 0.5D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'spmsrtls') then
          do i = 1, n
              x(i) = 0.2D0*sin(dble(i)**2)
          enddo
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'sphrpts') then
          if (mod(n,2)  .ne.  0) then
              info = 2
          else
              do i =1,n/2
                  x(2*i) = 0.0D0
                  x(2*i-1) = 8.0D0*acos(0.0D0)*dble(i)/dble(n)
              enddo
          endif
          rhobeg = 1.0D0/dble(n)
      elseif (trim(ful) .eq. 'srosenbr') then
          do i = 1, n
              if (mod(i,2) .eq. 1) x(i) = -1.2D0
              if (mod(i,2) .eq.0) x(i) = 1.0D0
          enddo
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'stmod') then
          x = 0.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'tointgss') then
          x = 3.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'tointtrig') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'tquartic') then
          x = 0.1D0
          rhobeg = 1.0D0
      elseif (trim(ful) .eq. 'trigsabs' .or. trim(ful) .eq. 'trigssqs') 
     +    then
          do i = 1, n
              xo = pi*cos(dble(n+i+163)/3.0D0)
              pt = pi*sin(dble(n*i*42)**(1.0D0/3.0D0))
              if (trim(ful) .eq. 'trigssqs') then
                  theta = (1.0D0+sin(dble(n*i*42)**(1.0D0/3.0D0)))/2.0D0
                  theta = 1.0D1**(-theta)
                  x(i) = (xo+pt*1.0D-1)/theta
              else
                  x(i) = xo+pt*1.0D-1
              endif
          enddo
          rhobeg = 1.0D-1
      elseif (trim(ful) .eq. 'trirose1') then
          x = -1.0D0
          rhobeg = 0.5D0
      elseif (trim(ful) .eq. 'trirose2') then
          x = 12.0D0
          rhobeg = 2.0D0
      elseif (trim(ful) .eq. 'vardim') then
          do i = 1, n
              x(i) = dble(n-i)/dble(n)
          enddo
          rhobeg = 1.0D0/dble(2*n)
      elseif (trim(ful) .eq. 'woods') then
          do i = 1, n
              if (mod(i,2) .eq. 1) x(i) = -3.0D0
              if (mod(i,2) .eq. 0) x(i) = -1.0D0
          enddo
          rhobeg = 1.0D0
      else
          info = 1
      endif 
      endsubroutine setuptest
