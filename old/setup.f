      subroutine setup (fun, n, x, rhobeg, fopt, info)
      implicit none
      character(len=*), intent(in) :: fun 
      integer(kind=4), intent(in) :: n
      real(kind=8), intent(out) :: x(n), rhobeg, fopt
      integer(kind=4), intent(out) :: info

!      integer(kind=4), external :: mod 

      integer(kind=4) :: i, j, k 
      real(kind=8) :: ysrf(ceiling(dble(n)), ceiling(dble(n))), scaling

      x = 1.0D308
      rhobeg = 1.0D0
      fopt = 1.0D308
      info = 0

      if (trim(fun) .eq. 'arglina') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'arglina4') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'arglinb') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'arglinc') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'argtrig') then
          do i = 1, n
              x(i) = 1.0D0/dfloat(n)
          enddo
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'arwhead') then
          x = 1.0D0
          rhobeg = 0.5D0
      elseif (trim(fun) .eq. 'bdqrtic') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'bdqrticp') then
          x = 1.0D0
          rhobeg = 0.5D0
      elseif (trim(fun) .eq. 'bdvalue') then
          do i = 1, n
              x(i) = dfloat(i*(i-n-1))/dfloat((n+1)*(n+1)) 
          enddo
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'brownal') then
          x = 5.0D-1
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'broydn3d') then
          x = -1.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'broydn7d') then
          x = -1.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'brybnd') then
          x = -1.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'chainwoo') then
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
      elseif (trim(fun) .eq. 'chebquad') then
          do i = 1, n
              x(i)=dfloat(i)/dfloat(n+1)
          enddo
          rhobeg=0.2D0*x(1)
      elseif (trim(fun) .eq. 'chpowellb') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'chpowells') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'chnrosnb') then
          x = -1.0D0
          !rhobeg = 0.5D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'chrosen') then
          x = -1.0D0
          rhobeg = 0.5D0
      elseif (trim(fun) .eq. 'cosine') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'cragglvy') then
          x = 2.0D0
          x(1) = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'curly10' .or. trim(fun) .eq. 'curly20'
     +       .or. trim(fun) .eq. 'curly30') then
          do i = 1, n
              x(i) = (1.0D-4)/(dfloat(n+1)) 
          enddo
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'cube') then
          x(1) = -1.2D0
          if (n .ge. 2) x(2) = 1.0D0
          do i = 3, n
              x(i) = x(i-2)
          enddo
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'dixmaane' .or. trim(fun) .eq. 'dixmaanf' 
     +   .or. trim(fun) .eq. 'dixmaang' .or. trim(fun) .eq. 'dixmaanh'
     +   .or. trim(fun) .eq. 'dixmaani' .or. trim(fun) .eq. 'dixmaanj'
     +   .or. trim(fun) .eq. 'dixmaank' .or. trim(fun) .eq. 'dixmaanl'
     +   .or. trim(fun) .eq. 'dixmaanm' .or. trim(fun) .eq. 'dixmaann'
     +   .or. trim(fun) .eq. 'dixmaano' .or. trim(fun) .eq. 'dixmaanp') 
     +   then
          do i = 1, n
              x(i) = 2.0D0 
          enddo
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'dqrtic') then
          x = 2.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'edensch') then
          x = 0.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'eg2') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'engval1') then
          x = 2.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'errinros') then
          x = -1.0D0
          !rhobeg = 0.5D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'extrosnb') then
          x = -1.0D0
          !rhobeg = 0.5D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'exttet') then
          x = 1.0D-1 
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'firose') then
          x = -2.0D0
          rhobeg = 0.5D0
      elseif (trim(fun) .eq. 'fletcbv2' .or. trim(fun) .eq. 'fletcbv3') 
     +    then
          do i = 1, n
              x(i) = dfloat(i)/dfloat(n+1)
          enddo
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'fletchcr') then
          do i = 1, n
              x(i) = 0.0D0
          enddo
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'fminsrf2') then
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
      elseif (trim(fun) .eq. 'freuroth') then
          x = 0.0D0
          x(1) = 0.5D0
          if (n .ge. 2) then
              x(2) = -2.0D0
          endif
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'genbrown') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'genhumps') then
          x = -5.062D2
          x(1) = -5.060D2
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'genrose') then
          do i = 1, n
              x(i) = dfloat(i)/dfloat(n+1) 
          enddo
C            rhobeg = 0.5D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'indef') then
          do i = 1, n
              x(i) = dfloat(i)/dfloat(n+1) 
          enddo
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'integreq') then
          do i = 1, n
              x(i) = dfloat(i*(i-n-1))/dfloat((n+1)*(n+1)) 
          enddo
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'liarwhd') then
          x = 4.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'lilifun3') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'lilifun4') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'morebv') then
          do i = 1, n
              x(i) = dfloat(i*(i-n-1))/dfloat((n+1)*(n+1)) 
          enddo
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'morebvl') then
          x = 0.0D0 ! Suggested by  Luksan, Matonoha, Vlcek, "Modified
                    ! CUTE probelems for unconstrained optimization.
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'ncb20') then
          if (n .lt. 20) then
              x = 0.0D0
          else 
              x(1:n-10) = 0.0D0
              x(n-9:n) = 1.0D0
          endif
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'ncb20b') then
          x = 0.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'noncvxu2') then
          do i = 1,n
              x(i) = dfloat(i)
          enddo
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'noncvxun') then
          do i = 1,n
              x(i) = dfloat(i)
          enddo
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'nondia') then
          x = -1.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'nondquar') then
          do i = 1, n
              if (mod(i,2) .eq. 1) x(i) = 1.0D0    
              if (mod(i,2) .eq. 0) x(i) = -1.0D0    
          enddo
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'penalty1') then
          do i = 1, n
              x(i) = dfloat(i)
          enddo
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'penalty2') then
          x = 0.5D0
          rhobeg = 0.1D0
      elseif (trim(fun) .eq. 'penalty3') then
          do i = 1, n
              x(i) = dfloat(i)/dfloat(n+1)
          enddo
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'penalty3p') then
          x = 0.0D0  ! Suggested by Powell, "The NEWUOA software 
C            for unconstrained optimization without
C            derivatives", 2004, Section 8.                       
          rhobeg = 0.1D0
      elseif(trim(fun) .eq. 'powellsg') then
          do i = 1, n
              if (mod(i,4) .eq. 1) x(i) = 3.0D0
              if (mod(i,4) .eq. 2) x(i) = -1.0D0
              if (mod(i,4) .eq. 3) x(i) = 0.0D0
              if (mod(i,4) .eq. 0) x(i) = 1.0D0
          enddo
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'power') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'rosenbrock') then
          x = -1.0D0
          rhobeg = 0.5D0
      elseif (trim(fun) .eq. 'sbrybnd') then
          scaling = 1.2D1
          do i = 1, n
              x(i) = exp(scaling*dfloat(1-i)/dfloat(n-1))
          enddo
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'sbrybndl') then
          scaling = 1.2D1
          do i = 1, n
              x(i) = exp(scaling*dfloat(1-i)/dfloat(n-1))
          enddo
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'schmvett') then
          x = 3.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'scosine') then
          scaling = 1.2D1
          do i = 1, n
              x(i) = exp(scaling*dfloat(1-i)/dfloat(n-1))
          enddo
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'scosinel') then
          scaling = 1.2D1
          do i = 1, n
              x(i) = exp(scaling*dfloat(1-i)/dfloat(n-1))
          enddo
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'serose') then
          x = -2.0D0
          rhobeg = 0.5D0
      elseif (trim(fun) .eq. 'sinquad') then
          x = 0.1D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'sparsine') then
          x = 0.5D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'sparsqur') then
          x = 0.5D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'spmsrtls') then
          do i = 1, n
              x(i) = 0.2D0*sin(dfloat(i)**2)
          enddo
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'sphrpts') then
          if (mod(n,2)  .ne.  0) then
              info = 2
          else
              do i =1,n/2
                  x(2*i) = 0.0D0
                  x(2*i-1) = 8.0D0*acos(0.0D0)*dfloat(i)/dfloat(n)
              enddo
          endif
          rhobeg = 1.0D0/dfloat(n)
      elseif (trim(fun) .eq. 'srosenbr') then
          do i = 1, n
              if (mod(i,2) .eq. 1) x(i) = -1.2D0
              if (mod(i,2) .eq.0) x(i) = 1.0D0
          enddo
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'tointgss') then
          x = 3.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'tointtrig') then
          x = 1.0D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'tquartic') then
          x = 0.1D0
          rhobeg = 1.0D0
      elseif (trim(fun) .eq. 'trirose1') then
          x = -1.0D0
          rhobeg = 0.5D0
      elseif (trim(fun) .eq. 'trirose2') then
          x = 12.0D0
          rhobeg = 2.0D0
      elseif (trim(fun) .eq. 'vardim') then
          do i = 1, n
              x(i) = dfloat(n-i)/dfloat(n)
          enddo
          rhobeg = 1.0D0/dfloat(2*n)
      elseif (trim(fun) .eq. 'woods') then
          do i = 1, n
              if (mod(i,2) .eq. 1) x(i) = -3.0D0
              if (mod(i,2) .eq. 0) x(i) = -1.0D0
          enddo
          rhobeg = 1.0D0
      else
          info = 1
      endif 
      endsubroutine setup
