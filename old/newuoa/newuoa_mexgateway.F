#include "fintrf.h"

      module global_var
      implicit none
      save
      mwPointer :: fun_ptr
      integer(kind=4) :: nf
      real(kind=8) :: fopt
      real(kind=8), allocatable :: fhist(:)
      real(kind=8), parameter :: INFINITY = 1.0D308
      end module global_var 
C======================================================================
C
C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C [xopt, fopt, nf, fhist] = newuoa(fun, x, rhobeg, rhoend, maxfun)

C     Declarations
      use global_var
      implicit none

C     mexFunction arguments:
      integer(kind=4), intent(in) :: nlhs, nrhs
      mwPointer, intent(in) :: prhs(nrhs)
      mwPointer, intent(inout) :: plhs(nlhs)

C     Function declarations:
      mwPointer :: mxCreateDoubleMatrix
      mwPointer :: mxCreateDoubleScalar
      mwPointer :: mxCreateNumericMatrix
      mwPointer :: mxGetPr
      mwSize :: mxGetM !, mxGetN
C      mwPointer :: mxCreateString
      integer(kind=4) :: mxClassIDFromClassName
      logical :: mxIsClass
C      logical :: mxIsDouble
C      logical :: mexPrintf
C      integer(kind=4) :: mxIsChar

C     Problem dimension 
      mwSize :: n

C     Workspace and variables
      real(kind=8), allocatable :: x(:), w(:)
      real(kind=8) :: rhobeg, rhoend, maxfun_r
      integer(kind=4) :: npt, iprint, maxfun

C-----------------------------------------------------------------------
C     Check for proper number of arguments.
      if (nrhs .ne. 5) then
         call mexErrMsgIdAndTxt ('MATLAB:interface:nInput',
     +                           '5 input argument required.')
      endif
      if (nlhs .gt. 4) then
         call mexErrMsgIdAndTxt ('MATLAB:interface:nOutput',
     +                           'At most 3 output arguments.')
      endif

CC     Validate inputs
CC     Check data type of input argument
C      if (mxIsDouble(prhs(1)) .eq. 0) then
C         call mexErrMsgIdAndTxt ('MATLAB:interface:NonDouble',
C     +                      'Input argument must be of type double.')
C      endif
C      if (mxIsComplex(prhs(1)) .eq. 1) then
C         call mexErrMsgIdAndTxt ('MATLAB:interface:NonReal',
C     +                           'Input argument must be real only')
C      endif

C     Get the size and pointers to input data
C     The first input must be a string.
!      if(mxIsClass(prhs(1), 'function_handle')) then
!          info = mexCallMatlab(1, plhs, 1, prhs, 'fun2str')
!      elseif(~mxIsClass(prhs(1), 'char')) then
!          call mexErrMsgIdAndTxt ('MATLAB:newuoa:InvalidFun',
!     +         'The first input must be a string or a function handle.')
!      endif
!!C     The input must be a row vector.
!!      elseif (mxGetM(prhs(1)) .ne. 1) then
!!         call mexErrMsgIdAndTxt ('MATLAB:newuoa:NonRowVector',
!!     +                          'The first input must be a row vector.')
!!      endif
!    
!C     Get the length of the input string.
!      fun_len = mxGetN(prhs(1))
!      allocate(character(len=fun_len) :: fun)
!      str_status = mxGetString(prhs(1), fun, fun_len)
!C     Check if mxGetString is successful.
!      if (str_status .ne. 0) then 
!         call mexErrMsgIdAndTxt ('MATLAB:newuoa:ReadError',
!     +                           'Error reading FunName.')
!      endif
      
      fun_ptr = prhs(1)
      n = mxGetM(prhs(2))
C      npt = (n+1)*(n+2)/2
C      npt = 2*n + 1
C Should be an input!!!!!!
      if (n .le. 5) then
          npt = (n+1)*(n+2)/2
      else
          npt = 2*n + 1
      endif
      if (allocated(x)) deallocate (x)
      allocate (x(n))
      call mxCopyPtrToReal8(mxGetPr(prhs(2)), x, n)
      call mxCopyPtrToReal8(mxGetPr(prhs(3)), rhobeg, 1)
      call mxCopyPtrToReal8(mxGetPr(prhs(4)), rhoend, 1)
      if (mxIsClass(prhs(5), 'double') .or. mxIsClass(prhs(5),'single'))
     + then 
         call mxCopyPtrToReal8(mxGetPr(prhs(5)), maxfun_r, 1)
         maxfun = int(maxfun_r)
      else
         call mxCopyPtrToInteger4(mxGetPr(prhs(5)), maxfun, 1)
      endif
      iprint = 1 
     
      if (allocated(w)) deallocate (w)
      allocate (w((npt+13)*(npt+n)+3*n*(n+3)/2))

      nf = 0
      fopt = INFINITY
      if (allocated(fhist)) deallocate (fhist)
      allocate (fhist(maxfun))
      fhist = INFINITY 

      call NEWUOA(n, npt, x, rhobeg, rhoend, iprint, maxfun, w)

C     Output
      plhs(1) = mxCreateDoubleMatrix(n, 1, 0)
      call mxCopyReal8ToPtr(x, mxGetPr(plhs(1)), n)
      plhs(2) = mxCreateDoubleScalar(fopt)
      plhs(3) = mxCreateNumericMatrix(1, 1,
     +          mxClassIDFromClassName('int32'), 0)
      call mxCopyInteger4ToPtr(nf, mxGetPr(plhs(3)), 1)
      plhs(4) = mxCreateDoubleMatrix(1, nf, 0)
      call mxCopyReal8ToPtr(fhist(1:nf), mxGetPr(plhs(4)), nf)

      !!!! FREE MEMORY !!!!!
      deallocate (x)
      deallocate (w)
      deallocate (fhist)

      return
      end subroutine mexFunction

      subroutine calfun(n, x, funval)
      use global_var
      implicit none
      integer(kind=4), intent(in) :: n
      real(kind=8), intent(in) :: x(n)
      real(kind=8), intent(out) :: funval

      integer(kind=4) :: mexCallMatlab
      mwPointer :: mxCreateDoubleMatrix
      mwPointer :: mxCreateDoubleScalar
      mwPointer :: mxGetPr

      integer(kind=4) :: info
      mwPointer :: plhs(1), prhs(2) 

      plhs(1) = mxCreateDoubleScalar(INFINITY)
      prhs(1) = fun_ptr ! fun_ptr is a global variable
      prhs(2) = mxCreateDoubleMatrix(n, 1, 0)
      call mxCopyReal8ToPtr(x, mxGetPr(prhs(2)), n)

      info = mexCallMatlab(1, plhs, 2, prhs, 'feval')
      if (info .ne. 0) then
          call mexErrMsgTxt('mex calls fun unsuccessfully.')
      endif

      call mxCopyPtrToReal8(mxGetPr(plhs(1)), funval, 1)

      !!!! FREE MEMORY !!!!!
      call mxDestroyArray(plhs(1))
      call mxDestroyArray(prhs(2))

      ! Update global variables
      nf = nf + 1
      fopt = min(funval, fopt)
      fhist(nf) = funval
      return
      end subroutine calfun
